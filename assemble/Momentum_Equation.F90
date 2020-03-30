!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

   module momentum_equation

      use spud
      use fldebug
      use global_parameters, only: FIELD_NAME_LEN
      use futils, only: int2str
      use parallel_tools
      use element_numbering, only: FAMILY_SIMPLEX
      use sparse_tools
      use linked_lists
      use parallel_fields, only: zero_non_owned
      use fields
      use profiler
      use sparse_tools_petsc
      use state_module
      use field_options
      use boundary_conditions
      use sparsity_patterns_meshes
      use sparse_matrices_fields
      use vtk_interfaces
      use dgtools, only: dg_apply_mass
      use state_fields_module
      use field_priority_lists
      use solvers
      use diagnostic_fields, only: calculate_diagnostic_variable
      use multiphase_module
      use divergence_matrix_cv
      use divergence_matrix_cg
      use coordinates
      use tidal_module
      use boundary_conditions_from_options
      use free_surface_module
      use petsc_solve_state_module
      use state_matrices_module
      use rotated_boundary_conditions
      use momentum_cg
      use slope_limiters_dg
      use momentum_dg
      use assemble_cmc
      use momentum_diagnostic_fields, only: calculate_momentum_diagnostics
      use compressible_projection
      use full_projection
      use hydrostatic_pressure
      use geostrophic_pressure
      use vertical_balance_pressure
      use foam_drainage, only: calculate_drainage_source_absor
      use oceansurfaceforcing
      use drag_module
      use pressure_dirichlet_bcs_cv
      use shallow_water_equations

      implicit none

      private
      public :: solve_momentum, momentum_equation_check_options

      ! The timestep
      real :: dt

      ! Are we going to form the Diagonal Schur complement preconditioner?
      logical :: get_diag_schur
      ! Do we need the scaled pressure mass matrix?
      logical :: get_scaled_pressure_mass_matrix
      ! Do we need an auxiliary matrix for full_projection solve?
      logical :: assemble_schur_auxiliary_matrix

      ! Do we want to use the compressible projection method?
      logical :: compressible_eos
      ! are we solving the shallow water equations (which partly follows the compressible projection path)
      logical :: shallow_water_projection
      ! Are we doing a full Schur solve?
      logical :: full_schur
      ! Are we lumping mass or assuming consistent mass?
      logical, dimension(:), allocatable :: lump_mass
      ! are we using a cv pressure 
      logical :: cv_pressure 
      ! for a CG pressure are we testing the continuity with cv
      logical :: cg_pressure_cv_test_continuity

      ! Do we need to reassemble the C^T or CMC matrices?
      logical :: reassemble_all_ct_m, reassemble_all_cmc_m

      ! Do we want to apply a theta weighting to the pressure gradient term?
      logical :: use_theta_pg

      ! Is a theta-weighting term present in the velocity-divergence?
      logical :: use_theta_divergence

      ! Are we using a discontinuous Galerkin discretisation?
      logical, dimension(:), allocatable :: dg
      ! True if advection-subcycling is performed
      logical, dimension(:), allocatable :: subcycle

      ! Apply KMK stabilisation?
      logical :: apply_kmk

      logical :: diagonal_big_m
      logical :: pressure_debugging_vtus

      ! Increased each call to momentum equation, used as index for pressure debugging vtus
      integer, save :: pdv_count = -1

      logical, dimension(:), allocatable :: sphere_absorption

      ! Are we running a multi-phase simulation?
      logical :: multiphase

      ! Do we have a prognostic free surface (currently only in 
      ! combination with a no_normal_stress free_surface)
      logical :: implicit_prognostic_fs, explicit_prognostic_fs, standard_fs

   contains

      subroutine solve_momentum(state, at_first_timestep, timestep)
         !!< Construct and solve the momentum and continuity equations
         !!< using Chorin's projection method (Chorin, 1968)

         ! An array of buckets full of fields
         ! The whole array is needed for the sake of multimaterial assembly
         type(state_type), dimension(:), intent(inout) :: state
         logical, intent(in) :: at_first_timestep
         integer, intent(in) :: timestep

         ! Counter iterating over each state
         integer :: istate 

         ! The pressure gradient matrix (extracted from state)
         type(block_csr_matrix_pointer), dimension(:), allocatable :: ct_m
         ! The pressure projection matrix (extracted from state)
         type(csr_matrix), pointer :: cmc_m

         ! logical to indicate whether ct_m and cmc_m need reassembling
         ! (used for each state within the assembly loop)
         logical :: reassemble_ct_m, reassemble_cmc_m
         ! is there a pressure in state?
         logical :: have_pressure
         ! Are we solving a Poisson pressure equation?
         logical :: poisson_p

         ! Matrix sparsity patterns for the matrices we allocate locally
         type(csr_sparsity), pointer :: u_sparsity

         !! Locally allocated matrices:
         ! Momentum LHS
         type(petsc_csr_matrix), dimension(:), allocatable, target :: big_m
         ! Matrix for split explicit advection
         type(block_csr_matrix), dimension(:), allocatable :: subcycle_m
         ! Pointer to matrix for full projection solve:
         type(petsc_csr_matrix_pointer), dimension(:), allocatable :: inner_m
         ! Pointer to preconditioner matrix for full projection solve:
         type(csr_matrix), pointer :: full_projection_preconditioner
         ! Auxiliary matrix for full_projection solve
         type(csr_sparsity), pointer :: schur_auxiliary_matrix_sparsity
         type(csr_matrix) :: schur_auxiliary_matrix
         ! Scaled pressure mass matrix - used for preconditioning full projection solve:
         type(csr_matrix), target :: scaled_pressure_mass_matrix
         type(csr_sparsity), pointer :: scaled_pressure_mass_matrix_sparsity
         ! Left hand matrix of CMC. For incompressible flow this points to ct_m as they are identical, 
         ! unless for CG pressure with CV tested continuity case when this matrix will be the 
         ! CV divergence tested matrix and ct_m the CG divergence tested matrix (right hand matrix of CMC).
         ! For compressible flow this differs to ct_m in that it will contain the variable density.
         type(block_csr_matrix_pointer), dimension(:), allocatable :: ctp_m
         ! The lumped mass matrix (may vary per component as absorption could be included)
         type(vector_field), dimension(1:size(state)) :: inverse_masslump
         ! Mass matrix
         type(petsc_csr_matrix), dimension(1:size(state)), target :: mass
         ! For DG:
         type(block_csr_matrix), dimension(1:size(state)):: inverse_mass

         ! Momentum RHS
         type(vector_field), dimension(1:size(state)):: mom_rhs
         ! Projection RHS
         type(scalar_field) :: projec_rhs
         ! RHS for continuity equation
         type(scalar_field), dimension(1:size(state)):: ct_rhs
         ! RHS for subcycling containing advection bc terms
         type(vector_field), dimension(1:size(state)):: subcycle_rhs

         ! Do we want to assemble the KMK stabilisation matrix?
         logical :: assemble_kmk

         ! Change in pressure
         type(scalar_field) :: delta_p
         ! Change in velocity
         type(vector_field) :: delta_u

         ! Dummy fields
         type(scalar_field), pointer :: dummyscalar, dummydensity, dummypressure

         ! Pressure and density
         type(scalar_field), pointer :: p, density
         ! this is the mesh on which the pressure projection is performed
         ! usually just p%mesh, but in case of a no_normal_stress free_surface
         ! (and in future hydrostatic projection) these are different
         type(mesh_type), pointer :: p_mesh
         ! Velocity and space
         type(vector_field), pointer :: u, x
         ! with a no_normal_stress free_surface we use a prognostic free surface
         type(scalar_field), pointer:: free_surface

         ! with free-surface or compressible pressure projection pressures 
         ! are at integer time levels and we apply a theta weighting to the
         ! pressure gradient term
         real :: theta_pg, theta_u
         ! With free-surface or compressible-projection the velocity divergence is
         ! calculated at time n+theta_divergence instead of at the end of the timestep
         real :: theta_divergence
         ! in this case p_theta=theta_pg*p+(1-theta_pg)*old_p
         type(scalar_field), pointer :: old_p, p_theta
         type(vector_field), pointer :: old_u
         ! all of this only applies if use_theta_pg .eqv. .true.
         ! without a free surface, or with a free surface and theta==1
         ! use_theta_pg .eqv. .false. and p_theta => p

         ! What is the equation type?
         character(len=FIELD_NAME_LEN) :: equation_type, poisson_scheme, schur_scheme, pressure_pmat

         integer :: stat

         ! The list of stiff nodes
         ! This is saved because the list is only formed when cmc is assembled, which
         ! isn't necessarily every time this subroutine is called but the list is
         ! still needed to fix the rhs (applying the fix to cmc itself wipes out the
         ! information that would be required to recompile the list)
         type(ilist), save :: stiff_nodes_list

         !! Variables for multi-phase flow model
         integer :: prognostic_count
         ! Do we have a prognostic pressure field to solve for?
         logical :: prognostic_p
         ! Prognostic pressure field's state index (if present)
         integer :: prognostic_p_istate 
         ! The 'global' CMC matrix (the sum of all individual phase CMC matrices)
         type(csr_matrix), pointer :: cmc_global
         ! An array of submaterials of the current phase in state(istate).
         type(state_type), dimension(:), pointer :: submaterials
         ! The index of the current phase (i.e. state(istate)) in the submaterials array
         integer :: submaterials_istate
         ! The full list of indices between submaterials and state
         integer, dimension(:), pointer :: submaterials_indices
         ! Do we have fluid-particle drag between phases?
         logical :: have_fp_drag

         ewrite(1,*) 'Entering solve_momentum'


         !! Get diagnostics (equations of state, etc) and assemble matrices

         ! Get some options that are independent of the states
         call get_option("/timestepping/timestep", dt)

         ! Are we running a multi-phase simulation?
         prognostic_count = option_count("/material_phase/vector_field::Velocity/prognostic")
         if(prognostic_count == 0) then
            return ! We don't have a velocity field to solve for, so exit.
         else if(prognostic_count > 1) then
            multiphase = .true.
         else
            multiphase = .false.
         end if
         ! Do we have fluid-particle drag (for multi-phase simulations)?
         have_fp_drag = have_option("/multiphase_interaction/fluid_particle_drag")

         ! Get the pressure p^{n}, and get the assembly options for the divergence and CMC matrices
         ! find the first non-aliased pressure
         do istate=1, size(state)
           p => extract_scalar_field(state(istate), "Pressure", stat)
           if (stat/=0) cycle
           if (.not. aliased(p)) exit
         end do
         have_pressure = istate<=size(state)

         if(.not. have_pressure) then
           ! Allocate a dummy scalar field in case we have no pressure
           allocate(dummypressure)
           ! pull a random velocity out of any of the states
           u => extract_vector_field(state, "Velocity")
           call allocate(dummypressure, u%mesh, "DummyPressure", field_type=FIELD_TYPE_CONSTANT)
           call zero(dummypressure)
           dummypressure%option_path = ""

           p => dummypressure
           p_mesh => u%mesh
           old_p => dummypressure
           free_surface => dummypressure

           prognostic_p=.false.
           standard_fs=.false.
           implicit_prognostic_fs=.false.
           explicit_prognostic_fs=.false.

         else

           p_mesh => p%mesh
           nullify(dummypressure)
           old_p => extract_scalar_field(state, "OldPressure", stat)
           if(stat/=0) then
             old_p => p
           end if
           prognostic_p = have_option(trim(p%option_path)//"/prognostic")
           if (prognostic_p) then
             prognostic_p_istate = istate
           end if

           u => extract_vector_field(state, "Velocity")
           standard_fs = has_standard_free_surface_bc(u)
           implicit_prognostic_fs = has_implicit_viscous_free_surface_bc(u)
           explicit_prognostic_fs = has_explicit_viscous_free_surface_bc(u)

           if(implicit_prognostic_fs.or.explicit_prognostic_fs) then
             free_surface => extract_scalar_field(state(istate), "FreeSurface")
             assert(have_option(trim(free_surface%option_path)//"/prognostic"))
           else
             free_surface => dummypressure
           end if

           if (implicit_prognostic_fs) then
             p_mesh => get_extended_pressure_mesh_for_viscous_free_surface(state(istate), &
                p%mesh, free_surface)
           end if

         end if

         !! Get some pressure options
         call get_pressure_options(p)

         ! Allocate arrays for the N states/phases
         allocate(big_m(size(state)))
         allocate(ct_m(size(state)))
         allocate(ctp_m(size(state)))
         allocate(subcycle_m(size(state)))
         allocate(inner_m(size(state)))

         nullify(cmc_global)

         ! Allocate arrays for phase-dependent options
         allocate(dg(size(state)))
         allocate(subcycle(size(state)))
         allocate(lump_mass(size(state)))
         allocate(sphere_absorption(size(state)))

         call profiler_tic("assembly_loop")
         assembly_loop: do istate = 1, size(state)

            ! Get the velocity u^{n}
            u => extract_vector_field(state(istate), "Velocity", stat)
            ! If there's no velocity then cycle
            if(stat/=0) cycle
            ! If this is an aliased velocity then cycle
            if(aliased(u)) cycle
            ! If the velocity isn't prognostic then cycle
            if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

            ! Calculate equations of state, etc.
            call profiler_tic("momentum_diagnostics")

            ! This sets up an array of the submaterials of a phase.
            ! NB: The submaterials array includes the current state itself, at index submaterials_istate.
            call get_phase_submaterials(state, istate, submaterials, submaterials_istate, submaterials_indices)
            call calculate_momentum_diagnostics(state, istate, submaterials, submaterials_istate, submaterials_indices)
            deallocate(submaterials)
            deallocate(submaterials_indices)

            call profiler_toc("momentum_diagnostics")

            ! Print out some statistics for the velocity
            ewrite_minmax(u)

            x => extract_vector_field(state(istate), "Coordinate")

            !! Get some velocity options:
            call get_velocity_options(state, istate, u)

            if(.not. have_pressure) then

               ! Don't bother solving for pressure if a pressure field doesn't exist
               nullify(ct_m(istate)%ptr)
               reassemble_ct_m = .false.
               nullify(cmc_m)
               reassemble_cmc_m = .false.

            else

               call profiler_tic(p, "assembly")
               ! Get the pressure gradient matrix (i.e. the divergence matrix)
               ! reassemble_ct_m is set to true if it does not already exist in state(i) 
               if (implicit_prognostic_fs) then
                 ct_m(istate)%ptr => get_extended_velocity_divergence_matrix(state(istate), u, free_surface, p_mesh, get_ct=reassemble_ct_m) 
               else
                 ct_m(istate)%ptr => get_velocity_divergence_matrix(state(istate), get_ct=reassemble_ct_m) 
               end if
               reassemble_ct_m = reassemble_ct_m .or. reassemble_all_ct_m
               
               ! For the CG pressure with CV tested continuity case 
               ! get the CV tested pressure gradient matrix (i.e. the divergence matrix)
               ! if required with a different unique name. Note there is no need
               ! to again decide reassemble_ct_m as ctp_m for this case is assembled when ct_m is.
               if (.not. (compressible_eos .or. shallow_water_projection) .and. cg_pressure_cv_test_continuity) then
                  if (implicit_prognostic_fs) then
                    ctp_m(istate)%ptr => get_extended_velocity_divergence_matrix(state(istate), u, free_surface, p_mesh, ct_m_name = "CVTestedExtendedVelocityDivergenceMatrix") 
                  else
                    ctp_m(istate)%ptr => get_velocity_divergence_matrix(state(istate), ct_m_name = "CVTestedVelocityDivergenceMatrix")
                  end if
               end if

               ! Get the pressure poisson matrix (i.e. the CMC/projection matrix)
               if (implicit_prognostic_fs) then
                 cmc_m => get_extended_pressure_poisson_matrix(state(istate), ct_m(istate)%ptr, p_mesh, get_cmc=reassemble_cmc_m)
               else
                 cmc_m => get_pressure_poisson_matrix(state(istate), get_cmc=reassemble_cmc_m)
               end if
               reassemble_cmc_m = reassemble_cmc_m .or. reassemble_all_cmc_m

               call profiler_toc(p, "assembly")

            end if
            ewrite_minmax(p)

            allocate(dummydensity)
            call allocate(dummydensity, x%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
            call set(dummydensity, 1.0)
            dummydensity%option_path = ""

            allocate(dummyscalar)
            call allocate(dummyscalar, x%mesh, "DummyScalar", field_type=FIELD_TYPE_CONSTANT)
            call zero(dummyscalar)
            dummyscalar%option_path = ""

            ! Depending on the equation type, extract the density or set it to some dummy field allocated above
            call get_option(trim(u%option_path)//"/prognostic/equation[0]/name", &
                           equation_type)
            select case(equation_type)
               case("LinearMomentum")
                  density=>extract_scalar_field(state(istate), "Density")
                  reassemble_cmc_m = reassemble_cmc_m .or. .not.constant_field(density)
               case("Boussinesq", "ShallowWater", "Drainage")
                  density=>dummydensity
               case default
                  ! developer error... out of sync options input and code
                  FLAbort("Unknown equation type for velocity")
            end select
            ewrite_minmax(density)

            if(full_schur) then
               ! Check to see whether pressure cmc_m preconditioning matrix is needed:
               call get_option(trim(p%option_path)//&
                           &"/prognostic/scheme/use_projection_method"//&
                           &"/full_schur_complement/preconditioner_matrix[0]/name", pressure_pmat)
             
               ! this is an utter mess, Rhodri, please clean up! 
               select case(pressure_pmat)
                  case("LumpedSchurComplement")
                     full_projection_preconditioner => cmc_m
                  case("DiagonalSchurComplement")
                     reassemble_cmc_m = .false.
                     get_diag_schur = .true.
                     full_projection_preconditioner => cmc_m
                  case("ScaledPressureMassMatrix")
                     reassemble_cmc_m = .false.
                     get_scaled_pressure_mass_matrix = .true.
                     full_projection_preconditioner => scaled_pressure_mass_matrix
                  case("NoPreconditionerMatrix")
                     reassemble_cmc_m = .false.
                     full_projection_preconditioner => cmc_m
                  case default
                     ! Developer error... out of sync options input and code
                     FLAbort("Unknown Matrix Type for Full_Projection")
               end select

               ! Decide on configuration of inner_m for full_projection solve:
               call get_option(trim(p%option_path)//&
                           &"/prognostic/scheme/use_projection_method"//&
                           &"/full_schur_complement/inner_matrix[0]/name", schur_scheme)
               select case(schur_scheme)
                  case("FullMassMatrix")
                     inner_m(istate)%ptr => mass(istate)
                  case("FullMomentumMatrix")
                     inner_m(istate)%ptr => big_m(istate)
                  case default
                     ! Developer error... out of sync options input and code
                     FLAbort("Unknown Matrix Type for Full_Projection")
               end select
            end if

            if (standard_fs .or. implicit_prognostic_fs .or. compressible_eos .or. shallow_water_projection) then
               ! this needs fixing for multiphase theta_pg could in principle be chosen
               ! per phase but then we need an array and we'd have to include theta_pg
               ! in cmc_m, i.e. solve for theta_div*dt*dp instead of theta_div*theta_pg*dt*dp
               ! theta_div can only be set once, but where and what is the default?
               if (has_boundary_condition(u, "free_surface") .and. multiphase) then
                 FLExit("Multiphase does not work with a free surface.")
               end if

               call get_option( trim(u%option_path)//'/prognostic/temporal_discretisation/theta', &
                     theta_u)
               ! With free surface or compressible-projection pressures are at integer
               ! time levels and we apply a theta-weighting to the pressure gradient term
               ! Also, obtain theta-weighting to be used in divergence term
               call get_option( trim(u%option_path)//'/prognostic/temporal_discretisation/theta_pressure_gradient', &
                     theta_pg, default=theta_u)
               use_theta_pg = (theta_pg/=1.0)
               call get_option( trim(u%option_path)//&
                     '/prognostic/temporal_discretisation/theta_divergence', &
                     theta_divergence, default=theta_u)
               use_theta_divergence = (theta_divergence/=1.0)
               ewrite(2,*) "Pressure gradient is evaluated at n+theta_pg"
               ewrite(2,*) "theta_pg: ", theta_pg
               ewrite(2,*) "Velocity divergence is evaluated at n+theta_divergence"
               ewrite(2,*) "theta_divergence: ", theta_divergence
               
               ! Note: Compressible multiphase simulations work, but only when use_theta_pg and use_theta_divergence
               ! are false. This needs improving - see comment above.
               if((compressible_eos .or. shallow_water_projection) .and. multiphase .and. (use_theta_pg .or. use_theta_divergence)) then
                  ewrite(-1,*) "Currently, for compressible multiphase flow simulations, the"
                  ewrite(-1,*) "temporal_discretisation/theta and temporal_discretisation/theta_divergence values"
                  ewrite(-1,*) "for each Velocity field must be set to 1.0."
                  FLExit("Multiphase does not work when use_theta_pg or use_theta_divergence are true.")
               end if
            else
               ! Pressures are, as usual, staggered in time with the velocities
               use_theta_pg=.false.
               use_theta_divergence=.false.
               theta_divergence=1.0
               theta_pg=1.0
            end if

            if (implicit_prognostic_fs) then
               allocate(p_theta)
               ! allocate p_theta on the extended mesh:
               call allocate(p_theta, p_mesh, "PressureAndFreeSurfaceTheta")
               p_theta%option_path=p%option_path ! Use p's solver options
               call copy_to_extended_p(p, free_surface, theta_pg, p_theta)
            else if (use_theta_pg) then
               allocate(p_theta)
               call allocate(p_theta, p_mesh, "PressureTheta")

               ! p_theta = theta*p + (1-theta)*old_p
               call set(p_theta, p, old_p, theta_pg)
               p_theta%option_path=p%option_path ! Use p's solver options
            else
               p_theta => p
               theta_pg=1.0
            end if

            call profiler_tic(u, "assembly")
            ! Allocation of big_m
            if(dg(istate)) then
               call allocate_big_m_dg(state(istate), big_m(istate), u)

               if(subcycle(istate)) then
                  u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
                  ! subcycle_m currently only contains advection, so diagonal=.true.
                  call allocate(subcycle_m(istate), u_sparsity, (/u%dim, u%dim/), &
                     diagonal=.true., name = "subcycle_m")
                  call allocate(subcycle_rhs(istate), u%dim, u%mesh, "SubCycleMomentumRHS")
               end if
            else
               ! Create a sparsity if necessary or pull it from state:
               u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
               ! and then allocate
               call allocate(big_m(istate), u_sparsity, (/u%dim, u%dim/), group_size=(/u%dim, u%dim/),&
                                       diagonal=diagonal_big_m, name="BIG_m")
            end if

            ! Initialise the big_m, ct_m and ctp_m matrices
            call zero(big_m(istate))
            if(reassemble_ct_m) then
               call zero(ct_m(istate)%ptr)         
               if (.not.(compressible_eos .or. shallow_water_projection) .and. cg_pressure_cv_test_continuity) then
                  call zero(ctp_m(istate)%ptr)
               end if
            end if

            ! Allocate the momentum RHS
            call allocate(mom_rhs(istate), u%dim, u%mesh, "MomentumRHS")
            call zero(mom_rhs(istate))
            ! Allocate the ct RHS
            call allocate(ct_rhs(istate), p_mesh, "DivergenceRHS")
            call zero(ct_rhs(istate))
            call profiler_toc(u, "assembly")

            if(has_scalar_field(state(istate), hp_name)) then
               call calculate_hydrostatic_pressure(state(istate))
            end if
            if(has_vector_field(state(istate), hpg_name)) then
               call calculate_hydrostatic_pressure_gradient(state(istate))
            end if
            if(has_scalar_field(state(istate), gp_name)) then
               call calculate_geostrophic_pressure_options(state(istate))
            end if

            if (has_vector_field(state(istate), "VelocityDrainageK1")) then
               call calculate_drainage_source_absor(state(istate))
            endif 

            ! Assemble the momentum equation
            call profiler_tic(u, "assembly")
            if(dg(istate)) then
               if(subcycle(istate)) then
                  call construct_momentum_dg(u, p, density, x, &
                     big_m(istate), mom_rhs(istate), state(istate), &
                     inverse_masslump=inverse_masslump(istate), &
                     inverse_mass=inverse_mass(istate), &
                     include_pressure_bcs=.not. cv_pressure, &
                     subcycle_m=subcycle_m(istate), subcycle_rhs=subcycle_rhs(istate))
               else
                  call construct_momentum_dg(u, p, density, x, &
                     big_m(istate), mom_rhs(istate), state(istate), &
                     inverse_masslump=inverse_masslump(istate), &
                     inverse_mass=inverse_mass(istate), &
                     include_pressure_bcs=.not. cv_pressure)
               end if
               if(has_scalar_field(state(istate), gp_name)) then
                  call subtract_geostrophic_pressure_gradient(mom_rhs(istate), state(istate))
               end if
            else
               ! This call will form the ct_rhs, which for compressible_eos
               ! or cg_pressure_cv_test_continuity is formed for a second time later below.
               call construct_momentum_cg(u, p, density, x, &
                     big_m(istate), mom_rhs(istate), ct_m(istate)%ptr, &
                     ct_rhs(istate), mass(istate), inverse_masslump(istate), &
                     state(istate), &
                     assemble_ct_matrix_here=reassemble_ct_m .and. .not. cv_pressure, &
                     include_pressure_and_continuity_bcs=.not. cv_pressure)
            end if
            
            ! If CV pressure then add in any dirichlet pressure BC integrals to the mom_rhs.
            if (cv_pressure) then
               call add_pressure_dirichlet_bcs_cv(mom_rhs(istate), u, p, state(istate))
            end if
            
            ! Add in multiphase interactions (e.g. fluid-particle drag) if necessary
            ! Note: this is done outside of construct_momentum_cg/dg to keep things
            ! neater in Momentum_CG/DG.F90, since we would need to pass around multiple phases 
            ! and their fields otherwise.
            if(multiphase .and. have_fp_drag) then
               call add_fluid_particle_drag(state, istate, u, x, big_m(istate), mom_rhs(istate))
            end if
            
            call profiler_toc(u, "assembly")

            if(has_scalar_field(state(istate), hp_name)) then
               call subtract_hydrostatic_pressure_gradient(mom_rhs(istate), state(istate))
            end if
            if(has_vector_field(state(istate), hpg_name)) then
               call subtract_hydrostatic_pressure_gradient(mom_rhs(istate), state(istate))
            end if
            if(has_scalar_field(state(istate), vbp_name)) then
               call calculate_vertical_balance_pressure(state(istate))
               call subtract_vertical_balance_pressure_gradient(mom_rhs(istate), state(istate))
            end if

            call profiler_tic(u, "assembly")
            if (has_boundary_condition(u, "wind_forcing")) then
               call wind_forcing(state(istate), mom_rhs(istate))
            end if

            if (has_boundary_condition(u, "drag")) then
               call drag_surface(big_m(istate), mom_rhs(istate), state(istate), density)
            end if

            call profiler_toc(u, "assembly")

            call profiler_tic(p, "assembly")

            ! Assemble divergence matrix C^T.
            ! At the moment cg does its own ct assembly. We might change this in the future.
            ! This call will form the ct_rhs, which for compressible_eos
            ! or cg_pressure_cv_test_continuity is formed for a second time later below.
            if(dg(istate) .and. .not. cv_pressure) then
               call assemble_divergence_matrix_cg(ct_m(istate)%ptr, state(istate), ct_rhs=ct_rhs(istate), &
                 test_mesh=p_theta%mesh, field=u, get_ct=reassemble_ct_m)
            end if

            if(cv_pressure) then
               ! This call will form the ct_rhs, which for compressible_eos
               ! is formed for a second time later below.
               call assemble_divergence_matrix_cv(ct_m(istate)%ptr, state(istate), ct_rhs=ct_rhs(istate), &
                                             test_mesh=p_theta%mesh, field=u, get_ct=reassemble_ct_m)
               if (implicit_prognostic_fs .and. reassemble_ct_m) then
                 call add_implicit_viscous_free_surface_integrals_cv(state(istate), &
                   ct_m(istate)%ptr, u, p_mesh, free_surface)
               end if
               if (explicit_prognostic_fs) then
                 call add_explicit_viscous_free_surface_integrals_cv(state(istate), &
                   ct_m(istate)%ptr, reassemble_ct_m, &
                   u, p_mesh, free_surface, mom_rhs=mom_rhs(istate))
               end if
            else
               if (implicit_prognostic_fs .and. reassemble_ct_m) then
                 call add_implicit_viscous_free_surface_integrals(state(istate), &
                   ct_m(istate)%ptr, u, p_mesh, free_surface)
               end if
               if (explicit_prognostic_fs) then
                 call add_explicit_viscous_free_surface_integrals(state(istate), &
                   mom_rhs(istate), ct_m(istate)%ptr, reassemble_ct_m, &
                   u, p_mesh, free_surface)
               end if
            end if

            call profiler_toc(p, "assembly")

            call profiler_tic(u, "assembly")
            if (have_rotated_bcs(u)) then
               ! Rotates big_m, rhs and the velocity field at strong, surface_aligned dirichlet bcs
               call rotate_momentum_equation(big_m(istate), mom_rhs(istate), u, state(istate), dg(istate))
               if (reassemble_ct_m) then
                  call rotate_ct_m(ct_m(istate)%ptr, u)
               end if
            end if
            if (sphere_absorption(istate)) then
               ! On the sphere inverse_masslump can currently only be assembled
               ! in the rotated frame. Thus we need to rotate anything that will
               ! interact with this.
               call rotate_momentum_to_sphere(big_m(istate), mom_rhs(istate), u, state(istate), dg(istate))
               if (reassemble_ct_m) then
                  call rotate_ct_m_sphere(state(istate), ct_m(istate)%ptr, u)
               end if
            end if

            call profiler_toc(u, "assembly")

            if (associated(ct_m(istate)%ptr)) then
              ewrite_minmax(ct_m(istate)%ptr)
            end if

            ! Do we want to solve for pressure?
            call profiler_tic(p, "assembly")
            
            if (prognostic_p) then
               
               ! Set up the left C matrix in CMC
               
               if(compressible_eos) then
                  ! Note: If we are running a compressible multiphase simulation then the C^T matrix for each phase becomes:
                  ! rho*div(vfrac*u) for each incompressible phase
                  ! rho*div(vfrac*u) + vfrac*u*grad(rho) for the single compressible phase.

                  allocate(ctp_m(istate)%ptr)
                  call allocate(ctp_m(istate)%ptr, ct_m(istate)%ptr%sparsity, (/1, u%dim/), name="CTP_m")
                  ! NOTE that this is not optimal in that the ct_rhs
                  ! was formed already above. The call here will overwrite those values.
                  if(cv_pressure) then
                     call assemble_compressible_divergence_matrix_cv(ctp_m(istate)%ptr, state, ct_rhs(istate))
                  else
                     call assemble_compressible_divergence_matrix_cg(ctp_m(istate)%ptr, state, istate, ct_rhs(istate))
                  end if               
               else if (shallow_water_projection) then
                 
                  assert(istate==1)
                  allocate(ctp_m(1)%ptr)
                  call allocate(ctp_m(1)%ptr, ct_m(1)%ptr%sparsity, (/1, u%dim/), name="CTP_m")
                  call assemble_swe_divergence_matrix_cg(ctp_m(1)%ptr, state(1), ct_rhs(1))

               else                  
                  ! Incompressible scenario
                  if (cg_pressure_cv_test_continuity) then
                     ! Form the CV tested divergence matrix and ct_rhs.
                     ! This will only reassemble ctp_m when ct_m 
                     ! also requires reassemble. NOTE that this is not optimal in that the ct_rhs
                     ! was formed already above. The call here will overwrite those values.
                     call assemble_divergence_matrix_cv(ctp_m(istate)%ptr, state(istate), ct_rhs=ct_rhs(istate), &
                                                        test_mesh=p_theta%mesh, field=u, get_ct=reassemble_ct_m)
                    if (implicit_prognostic_fs .and. reassemble_ct_m) then
                      call add_implicit_viscous_free_surface_integrals_cv(state(istate), &
                        ctp_m(istate)%ptr, u, p_mesh, free_surface)
                    end if
                    if (explicit_prognostic_fs) then
                      call add_explicit_viscous_free_surface_integrals_cv(state(istate), &
                        ctp_m(istate)%ptr, reassemble_ct_m, &
                        u, p_mesh, free_surface)
                    end if
                  else                  
                     ! ctp_m is identical to ct_m
                     ctp_m(istate)%ptr => ct_m(istate)%ptr
                  end if
               end if
               
               if (compressible_eos .or. shallow_water_projection .or. cg_pressure_cv_test_continuity) then
                  if (have_rotated_bcs(u)) then
                     if (dg(istate)) then
                       call zero_non_owned(u)
                     end if
                     call rotate_ct_m(ctp_m(istate)%ptr, u)
                  end if
                  if (sphere_absorption(istate)) then
                     if (dg(istate)) then
                       call zero_non_owned(u)
                     end if
                     call rotate_ct_m_sphere(state(istate), ctp_m(istate)%ptr, u)
                  end if
               end if

               ewrite_minmax(ctp_m(istate)%ptr)
               ewrite_minmax(ct_rhs(istate))

               ! Decide whether or not to form KMK stabilisation matrix:
               apply_kmk = (continuity(p_mesh) >= 0 .and. p_mesh%shape%degree == 1 &
                     & .and. p_mesh%shape%numbering%family == FAMILY_SIMPLEX &
                     & .and. continuity(u%mesh) >= 0 .and. u%mesh%shape%degree == 1 &
                     & .and. u%mesh%shape%numbering%family == FAMILY_SIMPLEX &
                     & .and. .not. have_option(trim(p%option_path) // &
                     & "/prognostic/spatial_discretisation/continuous_galerkin/remove_stabilisation_term") &
                     & .and. .not. cv_pressure)
               assemble_kmk = apply_kmk .and. &
                        ((.not. has_csr_matrix(state(istate), "PressureStabilisationMatrix")) .or. &
                        have_option(trim(p%option_path)// &
                        "/prognostic/scheme/update_discretised_equation") .or. &
                        have_option("/mesh_adaptivity/mesh_movement"))


               ! Assemble KMK stabilisation matrix if required:
               if(assemble_kmk) then
                  ewrite(2,*) "Assembling P1-P1 stabilisation"
                  call assemble_kmk_matrix(state(istate), p%mesh, x, theta_pg)
               end if

               if(full_schur) then
                  ! Decide whether we need to assemble an auxiliary matrix for full_projection solve:
                  if(apply_kmk) then
                     assemble_schur_auxiliary_matrix = .true.
                  end if
                  if (standard_fs .or. implicit_prognostic_fs) then
                     assemble_schur_auxiliary_matrix = .true.
                  end if

                  ! If schur_auxiliary_matrix is needed then assemble it:
                  if(assemble_schur_auxiliary_matrix) then
                     ! Get sparsity and assemble:
                     ewrite(2,*) "Assembling auxiliary matrix for full_projection solve"
                     if (implicit_prognostic_fs) then
                       schur_auxiliary_matrix_sparsity => get_extended_schur_auxillary_sparsity(state(istate), &
                         ct_m(istate)%ptr,  p_mesh)
                     else
                       schur_auxiliary_matrix_sparsity => get_csr_sparsity_secondorder(state(istate), p%mesh, u%mesh)
                     end if
                     call allocate(schur_auxiliary_matrix, schur_auxiliary_matrix_sparsity,&
                          name="schur_auxiliary_matrix")
                     ! Initialize matrix:
                     call zero(schur_auxiliary_matrix)
                     if(apply_kmk) then
                        ewrite(2,*) "Adding kmk stabilisation matrix to full_projection auxiliary matrix"
                        call add_kmk_matrix(state(istate), schur_auxiliary_matrix)
                     end if
                     if (standard_fs .or. implicit_prognostic_fs) then
                        ewrite(2,*) "Adding free surface to full_projection auxiliary matrix"
                        call add_free_surface_to_cmc_projection(state(istate), &
                                          schur_auxiliary_matrix, dt, theta_pg, &
                                          theta_divergence, assemble_cmc=.true., rhs=ct_rhs(istate))
                     end if
                  end if
               end if

               !! Assemble the appropriate projection matrix (CMC)
               if(reassemble_cmc_m) then
                  call zero(cmc_m)

                  if(dg(istate).and.(.not.lump_mass(istate))) then
                     call assemble_cmc_dg(cmc_m, ctp_m(istate)%ptr, ct_m(istate)%ptr, inverse_mass(istate))
                  else
                     call assemble_masslumped_cmc(cmc_m, ctp_m(istate)%ptr, inverse_masslump(istate), ct_m(istate)%ptr)

                     ! P1-P1 stabilisation
                     if (apply_kmk) then
                        ewrite(2,*) "Adding P1-P1 stabilisation matrix to cmc_m"
                        call add_kmk_matrix(state(istate), cmc_m)
                     end if
                  end if

                  if(have_option(trim(p%option_path)//"/prognostic/repair_stiff_nodes")) then
                     call repair_stiff_nodes(cmc_m, stiff_nodes_list)
                  end if

               end if ! end 'if(reassemble_cmc_m)'

               if (standard_fs .or. implicit_prognostic_fs) then
                  call add_free_surface_to_cmc_projection(state(istate), &
                           cmc_m, dt, theta_pg, theta_divergence, &
                           assemble_cmc=reassemble_cmc_m, rhs=ct_rhs(istate))
               end if
               
               if(get_diag_schur) then
                  ! Assemble diagonal Schur complement preconditioner:
                  call assemble_diagonal_schur(cmc_m, u, inner_m(istate)%ptr, ctp_m(istate)%ptr, ct_m(istate)%ptr)
                  ! P1-P1 stabilisation:
                  if (apply_kmk) then
                     ewrite(2,*) "Adding P1-P1 stabilisation to diagonal schur complement preconditioner matrix"
                     call add_kmk_matrix(state(istate), cmc_m)
                  end if
                  if (standard_fs .or. implicit_prognostic_fs) then
                     ewrite(2,*) "Adding free surface to diagonal schur complement preconditioner matrix"
                     call add_free_surface_to_cmc_projection(state(istate), &
                           cmc_m, dt, theta_pg, theta_divergence, assemble_cmc=.true.)
                  end if
               end if

               if(get_scaled_pressure_mass_matrix) then
                  ! Assemble scaled pressure mass matrix which will later be used as a
                  ! preconditioner in the full projection solve:
                  ewrite(2,*) "Assembling scaled pressure mass matrix preconditioner"
                  if (implicit_prognostic_fs) then
                    scaled_pressure_mass_matrix_sparsity => get_extended_schur_auxillary_sparsity(state(istate), &
                      ct_m(istate)%ptr, p_mesh)
                  else
                    scaled_pressure_mass_matrix_sparsity => get_csr_sparsity_firstorder(state(istate), p%mesh, p%mesh)
                  end if
                  call allocate(scaled_pressure_mass_matrix, scaled_pressure_mass_matrix_sparsity,&
                           name="scaled_pressure_mass_matrix")
                  call assemble_scaled_pressure_mass_matrix(state(istate),scaled_pressure_mass_matrix, p_mesh, dt)
                  if (implicit_prognostic_fs) then
                    call add_implicit_viscous_free_surface_scaled_mass_integrals(state(istate), scaled_pressure_mass_matrix, u, p, free_surface, dt)
                  end if
                  if (standard_fs .or. implicit_prognostic_fs) then
                     ewrite(2,*) "Adding free surface to scaled pressure mass matrix preconditioner"
                     call add_free_surface_to_cmc_projection(state(istate), &
                           scaled_pressure_mass_matrix, dt, theta_pg, theta_divergence, assemble_cmc=.true.)
                  end if
               end if

            end if ! end of prognostic pressure
            call profiler_toc(p, "assembly")
  
            if (associated(dummypressure)) then
               call deallocate(dummypressure)
               deallocate(dummypressure)
            end if
            call deallocate(dummydensity)
            deallocate(dummydensity)
            call deallocate(dummyscalar)
            deallocate(dummyscalar)

         end do assembly_loop
         call profiler_toc("assembly_loop")  ! End of Step 1 (diagnostics and matrix assembly)


         !! Obtain pressure guess p^{*} (assemble and solve a Poisson pressure equation for p^{*} if desired)


         ! Do we have a prognostic pressure field we can actually solve for?
         if(prognostic_p) then
            call profiler_tic(p, "assembly")

            u => extract_vector_field(state(prognostic_p_istate), "Velocity", stat)
            x => extract_vector_field(state(prognostic_p_istate), "Coordinate")

            ! Are we solving a Poisson pressure equation for a pressure guess p^{*}?
            call get_option(trim(p%option_path)//&
               "/prognostic/scheme/poisson_pressure_solution", poisson_scheme, &
               default="never")
            select case (poisson_scheme)
               case ("never")
                  poisson_p = .false.
               case ("only first timestep")
                  poisson_p = at_first_timestep
                  call set_option(trim(p%option_path)//&
                                    "/prognostic/scheme/poisson_pressure_solution", "never")
               case default
                  FLExit(trim(poisson_scheme)//" is not a legal poisson_pressure_solution")
            end select

            ! If desired, assemble Poisson pressure equation and get an initial guess at the pressure
            if(poisson_p) then   
               call solve_poisson_pressure(state, prognostic_p_istate, x, u, p, old_p, p_theta, theta_pg, &
                                           ct_m, ctp_m, mom_rhs, ct_rhs, inner_m, inverse_mass, &
                                           inverse_masslump, cmc_m, full_projection_preconditioner, schur_auxiliary_matrix)
            end if ! end of Poisson pressure solution


            ! Allocate RHS for pressure correction step
            call allocate(projec_rhs, p_mesh, "ProjectionRHS")
            call zero(projec_rhs)

            call profiler_toc(p, "assembly")
         end if ! end of prognostic pressure

         !! Advance velocity from u^{n} to an intermediate velocity u^{*}

         call profiler_tic("velocity_solve_loop")
         velocity_solve_loop: do istate = 1, size(state)

            ! Get the velocity u^{n}
            u => extract_vector_field(state(istate), "Velocity", stat)

            ! If there's no velocity then cycle
            if(stat/=0) cycle
            ! If this is an aliased velocity then cycle
            if(aliased(u)) cycle
            ! If the velocity isn't prognostic then cycle
            if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle
            
            if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
                 &"/continuous_galerkin").or.&
                 have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
                 &"/discontinuous_galerkin")) then

               x => extract_vector_field(state(istate), "Coordinate")

               if(use_theta_divergence) then
                  ! old_u is only used if use_theta_divergence, i.e. if theta_divergence/=1.0
                  old_u => extract_vector_field(state(istate), "OldVelocity")
                  if (old_u%aliased) then
                     ! in the case of one non-linear iteration, there is no OldVelocity,
                     ! it's just aliased to Velocity, therefore we make temp. version
                     allocate(old_u)
                     ! give it a distinct name, so we know to deallocate it
                     call allocate(old_u, u%dim, u%mesh, "TempOldVelocity")
                     call set(old_u, u)
                  end if
               end if
               
               call advance_velocity(state, istate, x, u, p_theta, big_m, ct_m, &
                    mom_rhs, subcycle_m, subcycle_rhs, inverse_mass)

               if(prognostic_p) then
                  call assemble_projection(state, istate, u, old_u, p, cmc_m, reassemble_cmc_m, cmc_global, ctp_m, &
                       ct_rhs, projec_rhs, p_theta, theta_pg, theta_divergence)
               end if

               ! Deallocate the old velocity field
               if(use_theta_divergence) then
                  if (old_u%name == "TempOldVelocity") then
                     call deallocate(old_u)
                     deallocate(old_u)
                  else if (have_rotated_bcs(u)) then
                     if (dg(istate)) then
                        call zero_non_owned(old_u)
                     end if
                     call rotate_velocity_back(old_u, state(istate))
                  end if
                  if (sphere_absorption(istate)) then
                     if (dg(istate)) then
                        call zero_non_owned(old_u)
                     end if
                     call rotate_velocity_back_sphere(old_u, state(istate))
                  end if
               end if
               
            end if ! end of prognostic velocity
            
         end do velocity_solve_loop
         call profiler_toc("velocity_solve_loop")


         !! Solve for delta_p -- the pressure correction term
         if(prognostic_p) then
            call profiler_tic(p, "assembly")
            
            ! Get the intermediate velocity u^{*} and the coordinate vector field
            u=>extract_vector_field(state(prognostic_p_istate), "Velocity", stat)
            x=>extract_vector_field(state(prognostic_p_istate), "Coordinate")
            
            if(multiphase) then
               cmc_m => cmc_global ! Use the sum over all individual phase CMC matrices
            end if
            
            call correct_pressure(state, prognostic_p_istate, x, u, p, old_p, delta_p, &
                 p_theta, free_surface, theta_pg, theta_divergence, &
                 cmc_m, ct_m, ctp_m, projec_rhs, inner_m, full_projection_preconditioner, &
                 schur_auxiliary_matrix, stiff_nodes_list)
            
            call deallocate(projec_rhs)
            call profiler_toc(p, "assembly")
         end if
         
         
         !! Correct and update velocity fields to u^{n+1} using pressure correction term delta_p
         if(prognostic_p) then
            
            call profiler_tic("velocity_correction_loop")
            velocity_correction_loop: do istate = 1, size(state)
               
               ! Get the velocity u^{*}
               u => extract_vector_field(state(istate), "Velocity", stat)
               
               ! If there's no velocity then cycle
               if(stat/=0) cycle
               ! If this is an aliased velocity then cycle
               if(aliased(u)) cycle
               ! If the velocity isn't prognostic then cycle
               if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle
               
               if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
                    &"/continuous_galerkin").or.&
                    have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
                    &"/discontinuous_galerkin")) then
                  
                  call profiler_tic(u, "assembly")
                  
                  ! Correct velocity according to new delta_p
                  if(full_schur) then
                     call correct_velocity_cg(u, inner_m(istate)%ptr, ct_m(istate)%ptr, delta_p, state(istate))
                  else if(lump_mass(istate)) then
                     call correct_masslumped_velocity(u, inverse_masslump(istate), ct_m(istate)%ptr, delta_p)
                  else if(dg(istate)) then
                     call correct_velocity_dg(u, inverse_mass(istate), ct_m(istate)%ptr, delta_p)
                  else
                     ! Something's gone wrong in the code
                     FLAbort("Don't know how to correct the velocity.")
                  end if
                  
                  if(implicit_prognostic_fs.or.explicit_prognostic_fs) then
                     call update_prognostic_free_surface(state(istate), free_surface, implicit_prognostic_fs, &
                          explicit_prognostic_fs)
                  end if
                  
                  call profiler_toc(u, "assembly")
                  
                  if(compressible_eos .or. shallow_water_projection) then
                     call deallocate(ctp_m(istate)%ptr)
                     deallocate(ctp_m(istate)%ptr)
                  end if
                  
               end if ! prognostic velocity
               
            end do velocity_correction_loop
            call profiler_toc("velocity_correction_loop")
            
            !! Deallocate some memory reserved for the pressure solve
            call deallocate(delta_p)
            
            if(assemble_schur_auxiliary_matrix) then
               ! Deallocate schur_auxiliary_matrix:
               call deallocate(schur_auxiliary_matrix)
            end if
            
            if(get_scaled_pressure_mass_matrix) then
               ! Deallocate scaled pressure mass matrix:
               call deallocate(scaled_pressure_mass_matrix)
            end if

         end if ! prognostic pressure

         !! Finalisation and memory deallocation
         call profiler_tic("finalisation_loop")
         finalisation_loop: do istate = 1, size(state)

            ! Get the velocity u^{n+1}
            u => extract_vector_field(state(istate), "Velocity", stat)

            ! If there's no velocity then cycle
            if(stat/=0) cycle
            ! If this is an aliased velocity then cycle
            if(aliased(u)) cycle
            ! If the velocity isn't prognostic then cycle
            if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

            if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
                                    &"/continuous_galerkin").or.&
               have_option(trim(u%option_path)//"/prognostic/spatial_discretisation"//&
                                    &"/discontinuous_galerkin")) then

               call finalise_state(state, istate, u, mass, inverse_mass, inverse_masslump, &
                                   big_m, mom_rhs, ct_rhs, subcycle_m, subcycle_rhs)

            end if

         end do finalisation_loop
         call profiler_toc("finalisation_loop")

         u => extract_vector_field(state, "Velocity")
         if(implicit_prognostic_fs .or. use_theta_pg) then
            call deallocate(p_theta)
            deallocate(p_theta)
         end if

         ! Deallocate arrays of matricies/fields/pointers
         deallocate(big_m)
         deallocate(ct_m)
         deallocate(ctp_m)
         deallocate(subcycle_m)
         deallocate(inner_m)

         if(multiphase .and. associated(cmc_global)) then
            call deallocate(cmc_global)
            deallocate(cmc_global)
         end if

         ! Deallocate arrays of options
         deallocate(dg)
         deallocate(subcycle)
         deallocate(lump_mass)
         deallocate(sphere_absorption)

      end subroutine solve_momentum


      subroutine get_velocity_options(state, istate, u)
         !!< Gets some velocity options from the options tree

         ! An array of buckets full of fields
         type(state_type), dimension(:), intent(inout) :: state
         integer, intent(in) :: istate

         type(vector_field), pointer :: u

         ! Local variables
         integer :: stat
         type(vector_field), pointer :: dummy_absorption
         logical :: have_viscosity, have_les, stress_form, partial_stress_form, have_coriolis
         logical :: on_sphere, have_absorption, have_vertical_stabilization

         ewrite(1,*) 'Entering get_velocity_options'


         dg(istate) = have_option(trim(u%option_path)//&
                           &"/prognostic/spatial_discretisation"//&
                           &"/discontinuous_galerkin")

         subcycle(istate) = have_option(trim(u%option_path)//&
            &"/prognostic/temporal_discretisation"//&
            &"/discontinuous_galerkin/maximum_courant_number_per_subcycle")

         ! Are we lumping the mass matrix?
         lump_mass(istate) = have_option(trim(u%option_path)//&
                           &"/prognostic/spatial_discretisation"//&
                           &"/continuous_galerkin/mass_terms"//&
                           &"/lump_mass_matrix").or.&
                     have_option(trim(u%option_path)//&
                           &"/prognostic/spatial_discretisation"//&
                           &"/discontinuous_galerkin/mass_terms"//&
                           &"/lump_mass_matrix")

         ! Here is where we try to decide how big big_m should be
         have_viscosity = have_option(trim(u%option_path)//&
            &"/prognostic/tensor_field::Viscosity")

         ! The following should include a dg option when a stress form version gets implemented
         stress_form = have_option(trim(u%option_path)//&
            &"/prognostic/spatial_discretisation/continuous_galerkin"//&
            &"/stress_terms/stress_form")

         partial_stress_form = have_option(trim(u%option_path)//&
            &"/prognostic/spatial_discretisation/continuous_galerkin"//&
            &"/stress_terms/partial_stress_form")

         have_les = have_option(trim(u%option_path)//"/prognostic/spatial_discretisation/"//&
            &"/continuous_galerkin/les_model")

         have_coriolis = have_option("/physical_parameters/coriolis")

         diagonal_big_m = .not.have_coriolis.and.(.not.((have_viscosity.or.have_les).and.(stress_form.or.partial_stress_form)))

         ! Do we want to rotate our equations to include absorption in a spherical geometry? 
         on_sphere = have_option('/geometry/spherical_earth/')
         dummy_absorption => extract_vector_field(state(istate), "VelocityAbsorption", stat)
         have_absorption = stat == 0
         have_vertical_stabilization = have_option(trim(u%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation").or. &
                                    have_option(trim(u%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy")
         sphere_absorption(istate) = on_sphere.and.(have_absorption.or.have_vertical_stabilization)

         shallow_water_projection = have_option(trim(u%option_path)//"/prognostic/equation::ShallowWater")
         ! NOTE: this relies on get_pressure_options() being called first!
         reassemble_all_cmc_m = reassemble_all_cmc_m .or. shallow_water_projection

      end subroutine get_velocity_options


      subroutine get_pressure_options(p)
         !!< Gets some pressure options from the options tree

         type(scalar_field), pointer :: p

         ewrite(1,*) 'Entering get_pressure_options'

         ! Are we using a compressible projection?
         compressible_eos = option_count("/material_phase/equation_of_state/compressible") > 0

         reassemble_all_cmc_m = have_option(trim(p%option_path)//&
                     "/prognostic/scheme/update_discretised_equation") .or. &
                     compressible_eos

         reassemble_all_ct_m = have_option(trim(p%option_path)//&
                     "/prognostic/scheme/update_discretised_equation")

         pressure_debugging_vtus = have_option(trim(p%option_path)// &
                     "/prognostic/output/debugging_vtus")
         if (pressure_debugging_vtus) then
            pdv_count = pdv_count+1
         end if

         get_diag_schur = .false.
         get_scaled_pressure_mass_matrix = .false.
         assemble_schur_auxiliary_matrix = .false.

         full_schur = have_option(trim(p%option_path)//&
                                 &"/prognostic/scheme&
                                 &/use_projection_method/full_schur_complement")

         ! Are we getting the pressure gradient matrix using control volumes?
         cv_pressure = have_option(trim(p%option_path)//&
                           "/prognostic/spatial_discretisation/control_volumes")

         ! For CG pressure are we testing the continuity with the CV dual 
         cg_pressure_cv_test_continuity = have_option(trim(p%option_path)//&
                   &"/prognostic/spatial_discretisation/continuous_galerkin&
                   &/test_continuity_with_cv_dual")

      end subroutine get_pressure_options


      subroutine solve_poisson_pressure(state, prognostic_p_istate, x, u, p, old_p, &
                                       p_theta, theta_pg, ct_m, ctp_m, &
                                       mom_rhs, ct_rhs, inner_m, inverse_mass, inverse_masslump, &
                                       cmc_m, full_projection_preconditioner, schur_auxiliary_matrix)
         !!< Solves a Poisson pressure equation for the pressure guess p^{*}

         ! An array of buckets full of fields
         type(state_type), dimension(:), intent(inout) :: state

         integer, intent(in) :: prognostic_p_istate

         type(block_csr_matrix), dimension(:), intent(in) :: inverse_mass
         type(vector_field), dimension(:), intent(in) :: inverse_masslump

         type(scalar_field), pointer :: p, p_theta, old_p
         type(vector_field), pointer :: x, u

         type(vector_field), dimension(:), intent(inout) :: mom_rhs
         type(scalar_field), dimension(:), intent(inout) :: ct_rhs

         type(csr_matrix), pointer :: cmc_m

         ! The pressure gradient matrices (extracted from state)
         type(block_csr_matrix_pointer), dimension(:), intent(inout) :: ct_m
         ! Compressible pressure gradient operator/left hand matrix of CMC
         type(block_csr_matrix_pointer), dimension(:), intent(inout) :: ctp_m

         ! Pointer to matrix for full projection solve:
         type(petsc_csr_matrix_pointer), dimension(:), intent(inout) :: inner_m
         ! Pointer to preconditioner matrix for full projection solve:
         type(csr_matrix), pointer :: full_projection_preconditioner

         type(csr_matrix), intent(in) :: schur_auxiliary_matrix

         type(vector_field), pointer :: positions

         real, intent(in) :: theta_pg

         !! Local variables
         type(scalar_field) :: poisson_rhs

         ewrite(1,*) 'Entering solve_poisson_pressure'

         call allocate(poisson_rhs, p_theta%mesh, "PoissonRHS")

         if (full_schur) then
            call assemble_poisson_rhs(poisson_rhs, ctp_m(prognostic_p_istate)%ptr, mom_rhs(prognostic_p_istate), ct_rhs(prognostic_p_istate), inner_m(prognostic_p_istate)%ptr, u, dt, theta_pg)
         else
            ! Get the RHS for the Poisson pressure equation...
            if(dg(prognostic_p_istate) .and. .not.lump_mass(prognostic_p_istate)) then
               call assemble_poisson_rhs_dg(poisson_rhs, ctp_m(prognostic_p_istate)%ptr, inverse_mass(prognostic_p_istate), mom_rhs(prognostic_p_istate), ct_rhs(prognostic_p_istate), u, dt, theta_pg)
            else
               ! Here we assume that we're using mass lumping if we're not using dg
               ! if this isn't true then this leads to inconsistent mass matrices in poisson_rhs and cmc_m
               ! but as we're only hoping to get a guesstimate of the pressure hopefully this won't be too
               ! bad.
               call assemble_masslumped_poisson_rhs(poisson_rhs, ctp_m(prognostic_p_istate)%ptr, mom_rhs(prognostic_p_istate), ct_rhs(prognostic_p_istate), inverse_masslump(prognostic_p_istate), u, dt, theta_pg)
            end if
         end if

         if (standard_fs .or. implicit_prognostic_fs) then
            call add_free_surface_to_poisson_rhs(poisson_rhs, state(prognostic_p_istate), dt, theta_pg)
         end if

         ! Apply strong dirichlet conditions
         call apply_dirichlet_conditions(cmc_m, poisson_rhs, p)

         positions => extract_vector_field(state(prognostic_p_istate), "Coordinate")
         call impose_reference_pressure_node(cmc_m, poisson_rhs, positions, trim(p%option_path))

         call profiler_toc(p, "assembly") ! Don't include Poisson solve
         if(full_schur) then
            if(assemble_schur_auxiliary_matrix) then
               call petsc_solve_full_projection(p_theta, ctp_m(prognostic_p_istate)%ptr, inner_m(prognostic_p_istate)%ptr, ct_m(prognostic_p_istate)%ptr, poisson_rhs, &
                  full_projection_preconditioner, u, state(prognostic_p_istate), u%mesh, &
                  auxiliary_matrix=schur_auxiliary_matrix)
            else
               call petsc_solve_full_projection(p_theta, ctp_m(prognostic_p_istate)%ptr, inner_m(prognostic_p_istate)%ptr, ct_m(prognostic_p_istate)%ptr, poisson_rhs, &
                  full_projection_preconditioner, u, state(prognostic_p_istate), u%mesh)
            end if
         else
            !! Go ahead and solve for the pressure guess p^{*}
            call petsc_solve(p_theta, cmc_m, poisson_rhs, state(prognostic_p_istate))
         end if
         call profiler_tic(p, "assembly")

         if (standard_fs .or. implicit_prognostic_fs) then
            ! Use this as initial pressure guess, except at the free surface
            ! where we use the prescribed initial condition
            call copy_poisson_solution_to_interior(state(prognostic_p_istate), &
               p_theta, p, old_p, u)
         end if

         if (pressure_debugging_vtus) then
            if (p_theta%mesh==p%mesh) then
              call vtk_write_fields("initial_poisson", pdv_count, x, p%mesh, &
                  sfields=(/ p_theta, p, old_p /))
            else
              ! for the prognostic fs, p_theta is on the extended mesh, so we omit it here
              call vtk_write_fields("initial_poisson", pdv_count, x, p%mesh, &
                  sfields=(/ p, old_p /))
            end if
         end if

         ewrite_minmax(p_theta)

         call deallocate(poisson_rhs)

      end subroutine solve_poisson_pressure


      subroutine advance_velocity(state, istate, x, u, p_theta, big_m, ct_m, mom_rhs, subcycle_m, subcycle_rhs, inverse_mass)
         !!< Solve momentum equation using pressure guess and advance velocity from u^{n} to u^{*}

         ! An array of buckets full of fields
         type(state_type), dimension(:), intent(inout) :: state
         integer, intent(in) :: istate

         type(vector_field), pointer :: x, u
         type(scalar_field), pointer :: p_theta

         ! Momentum LHS
         type(petsc_csr_matrix), dimension(:), target, intent(inout) :: big_m

         type(block_csr_matrix), dimension(:), intent(inout) :: inverse_mass

         ! The pressure gradient matrix (extracted from state)
         type(block_csr_matrix_pointer), dimension(:) :: ct_m

         type(vector_field), dimension(:), intent(inout) :: mom_rhs

         ! Matrix and rhs for split explicit advection
         type(block_csr_matrix), dimension(:), intent(in) :: subcycle_m
         type(vector_field), dimension(:), intent(in) :: subcycle_rhs

         !! Local variables
         ! Change in velocity
         type(vector_field) :: delta_u
         type(vector_field), pointer :: positions

         ! Fields for the subtract_out_reference_profile option under the Velocity field
         type(scalar_field), pointer :: hb_pressure
         type(scalar_field) :: combined_p
         integer :: stat

         ewrite(1,*) 'Entering advance_velocity'


         ! Allocate the momentum solution vector
         call profiler_tic(u, "assembly")
         call allocate(delta_u, u%dim, u%mesh, "DeltaU")
         delta_u%option_path = trim(u%option_path)

         ! Apply advection subcycling
         if(subcycle(istate)) then
            call subcycle_momentum_dg(u, mom_rhs(istate), subcycle_m(istate), &
              subcycle_rhs(istate), inverse_mass(istate), state(istate))
         end if

         if (associated(ct_m(istate)%ptr)) then
            ! add - ct_m^T*p to the rhs of the momentum eqn
            ! (delta_u is just used as dummy memory here)
            !
            ! despite multiplying pressure by a nonlocal operator
            ! a halo_update isn't necessary as this is just a rhs
            ! contribution
            if (have_option('/ocean_forcing/tidal_forcing') .or. &
               &have_option('/ocean_forcing/shelf')) then
            ewrite(1,*) "shelf: Entering compute_pressure_and_tidal_gradient"
               call compute_pressure_and_tidal_gradient(state(istate), delta_u, ct_m(istate)%ptr, p_theta, x)
            else if (have_option(trim(state(istate)%option_path)//'/equation_of_state/compressible/subtract_out_reference_profile')) then
               ! Splits up the Density and Pressure fields into a hydrostatic component (') and a perturbed component (''). 
               ! The hydrostatic components, denoted p' and rho', should satisfy the balance: grad(p') = rho'*g
               ! We subtract the hydrostatic component from the pressure used in the pressure gradient term of the momentum equation.
               hb_pressure => extract_scalar_field(state(istate), "HydrostaticReferencePressure", stat)
               if(stat /= 0) then
                  FLExit("When using the subtract_out_reference_profile option, please set a (prescribed) HydrostaticReferencePressure field.")
                  ewrite(-1,*) 'The HydrostaticReferencePressure field, defining the hydrostatic component of the pressure field, needs to be set.'
               end if
               call allocate(combined_p,p_theta%mesh, "PressurePerturbation")
               call set(combined_p, p_theta)
               call addto(combined_p, hb_pressure, scale=-1.0)
               call mult_T(delta_u, ct_m(istate)%ptr, combined_p)
               call deallocate(combined_p)
            else
               call mult_T(delta_u, ct_m(istate)%ptr, p_theta)
            end if

            if (dg(istate)) then
               ! We have just poluted the halo rows of delta_u. This is incorrect
               ! in the dg case due to the non-local assembly system employed.
               call zero_non_owned(delta_u)
            end if

            ewrite(2,*) 'note that delta_u = ct_m^T*p at this stage'
            ewrite_minmax(delta_u)
            call addto(mom_rhs(istate), delta_u)
         end if

         ! Impose zero guess on change in u
         call zero(delta_u)

         ! Impose any reference nodes on velocity
         positions => extract_vector_field(state(istate), "Coordinate")
         call impose_reference_velocity_node(big_m(istate), mom_rhs(istate), trim(u%option_path), positions)
         call apply_dirichlet_conditions(big_m(istate), mom_rhs(istate), u, dt)

         call profiler_toc(u, "assembly")

         !! Solve for the change in velocity
         call petsc_solve(delta_u, big_m(istate), mom_rhs(istate), state(istate))
         ewrite_minmax(delta_u)

         call profiler_tic(u, "assembly")
         ! Apply change to velocity field (Note that this gets stored in state)
         call addto(u, delta_u, dt)
         ewrite_minmax(u)

         call deallocate(delta_u)
         call profiler_toc(u, "assembly")

      end subroutine advance_velocity


      subroutine assemble_projection(state, istate, u, old_u, p, cmc_m, reassemble_cmc_m, cmc_global, ctp_m, &
                                     ct_rhs, projec_rhs, p_theta, theta_pg, theta_divergence)
         !!< Assembles the RHS for the projection solve step and, if required, the 'global' CMC matrix for multi-phase simulations.
         !!< Note that in the case of multi-phase simulations, projec_rhs contains the sum of ct_m*u over each prognostic velocity field,
         !!< and cmc_global contains the sum of the individual phase CMC matrices.

         ! An array of buckets full of fields
         type(state_type), dimension(:), intent(inout) :: state
         integer, intent(in) :: istate

         type(vector_field), pointer :: u, old_u
         type(scalar_field), pointer :: p, p_theta

         ! Compressible pressure gradient operator/left hand matrix of CMC
         type(block_csr_matrix_pointer), dimension(:), intent(inout) :: ctp_m
         ! The pressure projection matrix (extracted from state)
         type(csr_matrix), pointer :: cmc_m, cmc_global
         logical, intent(in) :: reassemble_cmc_m

         type(scalar_field), dimension(:), intent(inout) :: ct_rhs
         type(scalar_field), intent(inout) :: projec_rhs

         real, intent(in) :: theta_pg
         real, intent(in) :: theta_divergence

         ! Local variables
         type(scalar_field) :: kmk_rhs, temp_projec_rhs, compress_projec_rhs
         type(vector_field) :: delta_u
         integer :: stat

         ewrite(1,*) 'Entering assemble_projection'


         call profiler_tic(p, "assembly")
         ! Assemble the rhs
         ! If we are adding the P1-P1 stabilisation,
         ! this will have to have KMK * P added to it;
         !
         ! Despite multiplying velocity by a nonlocal operator
         ! a halo_update isn't necessary as this is just a rhs
         ! contribution
         call allocate(temp_projec_rhs, p_theta%mesh, "TempProjectionRHS")
         call zero(temp_projec_rhs)

         if (.not. use_theta_divergence) then
            ! Velocity divergence is evaluated at the end of the time step
            call mult(temp_projec_rhs, ctp_m(istate)%ptr, u)
         else
            ! Evaluate continuity at n+theta_divergence
            ! compute theta_divergence*u+(1-theta_divergence)*old_u
            call allocate(delta_u, u%dim, u%mesh, "VelocityTheta")
            if (have_rotated_bcs(u)) then
               if (dg(istate)) then
                  call zero_non_owned(old_u)
               end if
               call rotate_velocity(old_u, state(istate))
            end if
            if (sphere_absorption(istate)) then
               if (dg(istate)) then
                  call zero_non_owned(old_u)
               end if
               call rotate_velocity_sphere(old_u, state(istate))
            end if
            call set(delta_u, u, old_u, theta_divergence)
            call mult(temp_projec_rhs, ctp_m(istate)%ptr, delta_u)
            call deallocate(delta_u)
         end if

         ! Allocate the RHS
         call allocate(kmk_rhs, p_theta%mesh, "KMKRHS")
         call zero(kmk_rhs)

         if (apply_kmk) then
            call add_kmk_rhs(state(istate), kmk_rhs, p_theta, dt)
         end if
         ewrite_minmax(kmk_rhs)

         call addto(temp_projec_rhs, kmk_rhs)
         call scale(temp_projec_rhs, -1.0)
         call addto(temp_projec_rhs, ct_rhs(istate))
         ewrite_minmax(temp_projec_rhs)

         call deallocate(kmk_rhs)

         cmc_m => extract_csr_matrix(state(istate), "PressurePoissonMatrix", stat)
         
         if((compressible_eos .and. have_option(trim(state(istate)%option_path)//'/equation_of_state/compressible')) &
           .or. shallow_water_projection) then
            call allocate(compress_projec_rhs, p_theta%mesh, "CompressibleProjectionRHS")

            if (shallow_water_projection) then
               assert(istate==1)
               call assemble_shallow_water_projection(state(1), cmc_m, compress_projec_rhs, dt, &
                                                      theta_pg, theta_divergence, reassemble_cmc_m)
            else if(cv_pressure) then
               call assemble_compressible_projection_cv(state, cmc_m, compress_projec_rhs, dt, &
                                                      theta_pg, theta_divergence, reassemble_cmc_m)
            else
               call assemble_compressible_projection_cg(state, istate, cmc_m, compress_projec_rhs, dt, &
                                                      theta_pg, theta_divergence, reassemble_cmc_m)
            end if

            ewrite_minmax(compress_projec_rhs)
            ewrite_minmax(cmc_m)

            call addto(temp_projec_rhs, compress_projec_rhs)

            call deallocate(compress_projec_rhs)
         end if

         !! Add individual phase CMC matrix to 'global' CMC matrix
         if(multiphase) then
            if(.not.associated(cmc_global)) then
               ! If not yet allocated, allocate it here using the current CMC's sparsity pattern
               ! Assumes the same sparsity throughout (i.e. the same velocity and pressure mesh is used for each velocity field)
               allocate(cmc_global)
               call allocate(cmc_global, cmc_m%sparsity) 
               call zero(cmc_global)
            end if
            call addto(cmc_global, cmc_m)
         end if

         call addto(projec_rhs, temp_projec_rhs)
         ewrite_minmax(projec_rhs)

         call deallocate(temp_projec_rhs)

         call profiler_toc(p, "assembly")

      end subroutine assemble_projection

      subroutine correct_pressure(state, prognostic_p_istate, x, u, p, old_p, delta_p, &
                                 p_theta, free_surface, theta_pg, theta_divergence, &
                                 cmc_m, ct_m, ctp_m, projec_rhs, inner_m, full_projection_preconditioner, &
                                 schur_auxiliary_matrix, stiff_nodes_list)
         !!< Finds the pressure correction term delta_p needed to make the intermediate velocity field (u^{*}) divergence-free         

         ! An array of buckets full of fields
         type(state_type), dimension(:), intent(inout) :: state
         type(vector_field), pointer :: x, u
         type(scalar_field), pointer :: p, old_p, p_theta, free_surface
         type(scalar_field), intent(inout) :: delta_p

         integer, intent(in) :: prognostic_p_istate

         real, intent(inout) :: theta_pg
         real, intent(inout) :: theta_divergence

         ! The pressure projection matrix (extracted from state)
         type(csr_matrix), pointer :: cmc_m

         ! The pressure gradient matrix (extracted from state)
         type(block_csr_matrix_pointer), dimension(:), intent(inout) :: ct_m
         ! Compressible pressure gradient operator/left hand matrix of CMC
         type(block_csr_matrix_pointer), dimension(:), intent(inout) :: ctp_m

         ! Projection RHS
         type(scalar_field), intent(inout) :: projec_rhs

         ! Pointer to matrix for full projection solve:
         type(petsc_csr_matrix_pointer), dimension(:), intent(inout) :: inner_m
         ! Pointer to preconditioner matrix for full projection solve:
         type(csr_matrix), pointer :: full_projection_preconditioner

         type(csr_matrix), intent(in) :: schur_auxiliary_matrix

         type(ilist), intent(inout) :: stiff_nodes_list

         type(vector_field), pointer :: positions

         ewrite(1,*) 'Entering correct_pressure'


         ! Apply strong Dirichlet conditions
         ! we're solving for "delta_p"=theta_pg*theta_divergence*dp*dt, where dp=p_final-p_current
         ! apply_dirichlet_condition however assumes we're solving for
         ! "acceleration" dp/dt, by providing dt=1/(dt*theta_pg**2) we get what we want
         call apply_dirichlet_conditions(cmc_m, projec_rhs, p, &
                                         dt=1.0/(dt*theta_pg*theta_divergence))

         positions => extract_vector_field(state(prognostic_p_istate), "Coordinate")
         call impose_reference_pressure_node(cmc_m, projec_rhs, positions, trim(p%option_path))

         ! Allocate the change in pressure field
         call allocate(delta_p, p_theta%mesh, "DeltaP")
         delta_p%option_path = trim(p%option_path)
         call zero(delta_p)

         if(have_option(trim(p%option_path)//"/prognostic/repair_stiff_nodes")) then
            call zero_stiff_nodes(projec_rhs, stiff_nodes_list)
         end if
         call profiler_toc(p, "assembly")

         ! Solve for the change in pressure, delta_p
         if(full_schur) then
            if(assemble_schur_auxiliary_matrix) then
               call petsc_solve_full_projection(delta_p, ctp_m(prognostic_p_istate)%ptr, inner_m(prognostic_p_istate)%ptr, ct_m(prognostic_p_istate)%ptr, projec_rhs, &
                 full_projection_preconditioner, u, state(prognostic_p_istate), u%mesh, &
                 auxiliary_matrix=schur_auxiliary_matrix)
            else
               call petsc_solve_full_projection(delta_p, ctp_m(prognostic_p_istate)%ptr, inner_m(prognostic_p_istate)%ptr, ct_m(prognostic_p_istate)%ptr, projec_rhs, &
                 full_projection_preconditioner, u, state(prognostic_p_istate), u%mesh)
            end if
         else
            call petsc_solve(delta_p, cmc_m, projec_rhs, state(prognostic_p_istate))
         end if

         ewrite_minmax(delta_p)

         if (pressure_debugging_vtus) then
            ! Writes out the pressure and velocity before the correction is added in
            ! (as the corrected fields are already available in the convergence files)
            if (p%mesh==p_theta%mesh) then
              call vtk_write_fields("pressure_correction", pdv_count, x, p%mesh, &
                  sfields=(/ p, old_p, p_theta /))
              ! same thing but now on velocity mesh:
              call vtk_write_fields("velocity_before_correction", pdv_count, x, u%mesh, &
                  sfields=(/ p, old_p, p_theta /), vfields=(/ u /))
            else
              ! this is the case for the prognostic fs where p_theta is on an extended mesh
              ! we simply omit p_theta - additional output could be implemented
              call vtk_write_fields("pressure_correction", pdv_count, x, p%mesh, &
                  sfields=(/ p, old_p /))
              call vtk_write_fields("velocity_before_correction", pdv_count, x, u%mesh, &
                  sfields=(/ p, old_p /), vfields=(/ u /))
            end if
         end if

         call profiler_tic(p, "assembly")
         if (use_theta_divergence) then
            ! We've solved theta_pg*theta_divergence*dt*dp, in the velocity correction
            ! however we need theta_pg*dt*dp
            call scale(delta_p, 1.0/theta_divergence)
         end if

         if (implicit_prognostic_fs) then
           call update_pressure_and_viscous_free_surface(state(prognostic_p_istate), p, free_surface, delta_p, theta_pg)
         else
           ! Add the change in pressure to the pressure
           ! (if .not. use_theta_pg then theta_pg is 1.0)
           call addto(p, delta_p, scale=1.0/(theta_pg*dt))
         end if
         ewrite_minmax(p)

         if(compressible_eos) then
            call update_compressible_density(state)
         end if

      end subroutine correct_pressure

      subroutine finalise_state(state, istate, u, mass, inverse_mass, inverse_masslump, &
                               big_m, mom_rhs, ct_rhs, subcycle_m, subcycle_rhs)
         !!< Does some finalisation steps to the velocity field and deallocates some memory
         !!< allocated for the specified state.

         ! An array of buckets full of fields
         type(state_type), dimension(:), intent(inout) :: state

         integer, intent(in) :: istate
         type(vector_field), pointer :: u

         ! Mass matrix
         type(petsc_csr_matrix), dimension(:), target, intent(inout) :: mass
         ! For DG:
         type(block_csr_matrix), dimension(:), intent(inout) :: inverse_mass
         ! The lumped mass matrix (may vary per component as absorption could be included)
         type(vector_field), dimension(:), intent(inout) :: inverse_masslump

         ! Momentum LHS
         type(petsc_csr_matrix), dimension(:), target, intent(inout) :: big_m
         ! Momentum RHS
         type(vector_field), dimension(:), intent(inout) :: mom_rhs
         type(scalar_field), dimension(:), intent(inout) :: ct_rhs
         ! Matrix and rhs for split explicit advection
         type(block_csr_matrix), dimension(:), intent(inout) :: subcycle_m
         type(vector_field), dimension(:), intent(inout) :: subcycle_rhs

         ! Slope limiter variables:
         integer :: d
         type(scalar_field) :: u_cpt

         ewrite(1,*) 'Entering finalise_state'

         call profiler_tic(u, "assembly")
         if (have_rotated_bcs(u)) then
            if (dg(istate)) then
               call zero_non_owned(u)
            end if
            call rotate_velocity_back(u, state(istate))
         end if
         if (sphere_absorption(istate)) then
            if (dg(istate)) then
               call zero_non_owned(u)
            end if
            call rotate_velocity_back_sphere(u, state(istate))
         end if
         if (subcycle(istate)) then
            ! Filter wiggles from u
            do d = 1, mesh_dim(u)
               u_cpt = extract_scalar_field_from_vector_field(u, d)
               call limit_vb(state(istate), u_cpt)
            end do
         end if
         call profiler_toc(u, "assembly")

         if(dg(istate)) then
            if(lump_mass(istate)) then
               call deallocate(inverse_masslump(istate))
            else
               call deallocate(inverse_mass(istate))
            end if
         else
            call deallocate_cg_mass(mass(istate), inverse_masslump(istate))
         end if

         call deallocate(mom_rhs(istate))
         call deallocate(ct_rhs(istate))
         call deallocate(big_m(istate))
         if(subcycle(istate)) then
            call deallocate(subcycle_m(istate))
            call deallocate(subcycle_rhs(istate))
         end if

      end subroutine finalise_state

      subroutine momentum_equation_check_options

         integer :: i, nmat 
         character(len=FIELD_NAME_LEN) :: schur_scheme
         character(len=FIELD_NAME_LEN) :: schur_preconditioner
         character(len=FIELD_NAME_LEN) :: pressure_mesh
         character(len=FIELD_NAME_LEN) :: pressure_mesh_element_type

         ewrite(1,*) 'Checking momentum discretisation options'

         nmat = option_count("/material_phase")

         do i = 0, nmat-1

            if(have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic"//&
                                 &"/reference_node").and.&
            have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic"//&
                                 &"/solver/remove_null_space")) then
               FLExit("Can't set a pressure reference node and remove the null space.")
            end if

            if(have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic"//&
                                 &"/spatial_discretisation/discontinuous_galerkin") &
              .and. .not. have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/continuous_galerkin")) then
              FLExit("With discontinuous galerkin Pressure you need a continuous Velocity")
            end if
 
            if(have_option("/material_phase["//int2str(i)//&
                        &"]/vector_field::Velocity/prognostic/reference_node")) then
               if((.not.(have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/continuous_galerkin/mass_terms/exclude_mass_terms").and. &
                           have_option("/material_phase["//int2str(i)//&
                                 "]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/continuous_galerkin/advection_terms/exclude_advection_terms"))).and. &
                  (.not.(have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/discontinuous_galerkin/mass_terms/exclude_mass_terms").and. &
                           have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/discontinuous_galerkin/advection_scheme/none")))) then
                  ewrite(-1,*) "Error: You have set a Velocity reference node but don't appear"
                  ewrite(-1,*) "to be solving the Stokes equation."
                  ewrite(-1,*) "Setting a reference node for Velocity only makes sense if both"
                  ewrite(-1,*) "mass and advection terms are excluded from the momentum"
                  ewrite(-1,*) "equation.  Even then whether it is valid depends on your"
                  ewrite(-1,*) "boundary conditions and whether your velocity components"
                  ewrite(-1,*) "are coupled but I can't check for that."
                  FLExit("Don't set a Velocity reference_node unless solving the Stokes equation.")
               end if
            end if


            if(have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/continuous_galerkin")&
               .and.(.not.have_option("/material_phase["//int2str(i)//&
                                    &"]/vector_field::Velocity/prognostic"//&
                                    &"/spatial_discretisation/continuous_galerkin"//&
                                    &"/mass_terms/lump_mass_matrix"))) then

               if(have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic"//&
                                 &"/scheme/use_projection_method")) then
                  if(.not.have_option("/material_phase["//int2str(i)//&
                                    &"]/scalar_field::Pressure/prognostic"//&
                                    &"/scheme/use_projection_method"//&
                                    &"/full_schur_complement")) then
                     ewrite(-1,*) "Error: You're not lumping the velocity mass matrix"
                     ewrite(-1,*) "but haven't selected any schur complement options."
                     ewrite(-1,*) "Are you sure you don't want to lump the mass?"
                     ewrite(-1,*) "The consistent mass method is VERY SLOW!"
                     ewrite(-1,*) "If you really want to use it then its under"
                     ewrite(-1,*) "the projection scheme options underneath Pressure."
                     ewrite(-1,*) "Otherwise switch on mass lumping under vector_field::Velocity/"
                     ewrite(-1,*) "spatial_discretisation/continuous_galerkin/"
                     ewrite(-1,*) "mass_terms/lump_mass_matrix"
                     FLExit("Good luck!")
                  end if
               end if

            end if

            if(have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/tensor_field::Viscosity/prescribed/value"//&
                                 &"/isotropic").or. &
               have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/tensor_field::Viscosity/prescribed/value"//&
                                 &"/diagonal")) then

               if(have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/continuous_galerkin/stress_terms/stress_form").or.&
                  have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/continuous_galerkin/stress_terms/partial_stress_form")) then
                  ewrite(-1,*) "You have selected stress form viscosity but have entered an isotropic or"
                  ewrite(-1,*) "diagonal Viscosity tensor."
                  ewrite(-1,*) "Zero off diagonal entries in the Viscosity tensor do not make physical"
                  ewrite(-1,*) "sense when using stress form viscosity."
                  ewrite(-1,*) "Use tensor_form or anisotropic_symmetric Viscosity instead."
                  FLExit("Use tensor_form or anisotropic_symmetric Viscosity.")
               end if

            end if

            if(have_option("/material_phase["//int2str(i)//"]/vector_field::Velocity/prognostic/"//&
               &"spatial_discretisation/continuous_galerkin/temperature_dependent_viscosity")) then

               if(.not.have_option("/material_phase["//int2str(i)//"]/scalar_field::Temperature")) then
                  FLExit("You must have a temperature field to have a temperature dependent viscosity.")
               end if

               ewrite(-1,*) "Warning - any viscosity values set under tensor_field::Viscosity will be"
               ewrite(-1,*) "overwritten by a calculated temperature dependent viscosity. Nonetheless,"
               ewrite(-1,*) "to ensure that the viscosity tensor is simulated in the correct form, please"
               ewrite(-1,*) "select a form under tensor_field::Viscosity. Note that only partial stress and"
               ewrite(-1,*) "stress form are valid for a spatially varying viscosity field."

               if(have_option("/material_phase["//int2str(i)//&
                  &"]/vector_field::Velocity/prognostic"//&
                  &"/tensor_field::Viscosity/prescribed/value"//&
                  &"/isotropic").or.&
                  &have_option("/material_phase["//int2str(i)//&
                  &"]/vector_field::Velocity/prognostic"//&
                  &"/tensor_field::Viscosity/prescribed/value"//&
                  &"/diagonal")) then

                  ewrite(-1,*) "A spatially varying viscosity (for example a viscosity that depends"
                  ewrite(-1,*) "upon a spatiall varying temperature field) is only valid with stress"
                  ewrite(-1,*) "or partial stress form viscosity"
                  FLExit("For a spatially varying viscosity field, use stress or partial stress form viscosity")
               end if

            end if

            if(have_option("/material_phase["//int2str(i)//&
               &"]/vector_field::Velocity/prognostic"//&
               &"/tensor_field::Viscosity/prescribed/value"//&
               &"/anisotropic_symmetric").or.&
            have_option("/material_phase["//int2str(i)//&
               &"]/vector_field::Velocity/prognostic"//&
               &"/tensor_field::Viscosity/prescribed/value"//&
               &"/anisotropic_asymmetric")) then

               if(have_option("/material_phase["//int2str(i)//&
                  &"]/scalar_field::Pressure/prognostic"//&
                  &"/scheme/use_projection_method")) then

                  if(have_option("/material_phase["//int2str(i)//&
                     &"]/scalar_field::Pressure/prognostic"//&
                     &"/scheme/use_projection_method"//&
                     &"/full_schur_complement")) then

                     call get_option("/material_phase["//int2str(i)//&
                        &"]/scalar_field::Pressure/prognostic"//&
                        &"/scheme/use_projection_method"//&
                        &"/full_schur_complement/preconditioner_matrix[0]/name", schur_preconditioner)

                     select case(schur_preconditioner)
                        case("ScaledPressureMassMatrix")
                           ewrite(-1,*) "WARNING - At present, the viscosity scaling for the pressure mass matrix is"
                           ewrite(-1,*) "taken from the 1st component of the viscosity tensor. Such a scaling"
                           ewrite(-1,*) "is only valid when all components of each viscosity tensor are constant."
                     end select

                  end if

               end if

            end if

            if(have_option("/material_phase["//int2str(i)//&
                                 &"]/vector_field::Velocity/prognostic"//&
                                 &"/spatial_discretisation/discontinuous_galerkin")) then

               if(have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic"//&
                                 &"/scheme/use_projection_method")) then
                  if(have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic"//&
                                 &"/scheme/use_projection_method"//&
                                 &"/full_schur_complement")) then

                     call get_option("/material_phase["//int2str(i)//&
                                    &"]/scalar_field::Pressure/prognostic"//&
                                    &"/scheme/use_projection_method"//&
                                    &"/full_schur_complement/inner_matrix[0]/name", schur_scheme)
                     select case(schur_scheme)
                        case("FullMassMatrix")
                           FLExit("Can't do a full schur complement solve with dg velocity and a mass inner matrix.")
                     end select

                  end if
               end if

            end if

            ! Check options for case with CG pressure and
            ! testing continuity with CV dual mesh. 
            ! Will not work with compressible, free surface or 
            ! wetting and drying. Also will not work if the pressure is on a mesh that has 
            ! bubble or trace shape functions.
            if (have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic&
                                 &/spatial_discretisation/continuous_galerkin&
                                 &/test_continuity_with_cv_dual")) then
               
               ! Check that the incompressible projection is being used
               if(.not.have_option("/material_phase["//int2str(i)//&
                                 &"]/scalar_field::Pressure/prognostic&
                                 &/scheme/use_projection_method")) then

                  ewrite(-1,*) "Error: For a CG Pressure the continuity"
                  ewrite(-1,*) "can only be tested with the cv dual mesh"
                  ewrite(-1,*) "if the pressure scheme is the incompressible"
                  ewrite(-1,*) "projection method, which is given by the option"
                  ewrite(-1,*) "path material_phase/Pressure/spatial_discretisation/"
                  ewrite(-1,*) "continuous_galerkin/scheme/use_projection_method"
                  FLExit("Use incompressible projection method if wanting to test continuity with cv dual with CG pressure")                  
               end if
               
               ! Check that the wetting_and_drying model is not being used
               if(have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")) then
                  FLExit("For CG Pressure cannot test the continuity equation with CV when using the wetting and drying model")
               end if
               
               ! get the pressure mesh name
               call get_option("/material_phase["//int2str(i)//"]/scalar_field::Pressure/prognostic/mesh/name", &
                                pressure_mesh)
               
               ! check that the pressure mesh options 
               ! do NOT say bubble or trace
               call get_option("/geometry/mesh::"//trim(pressure_mesh)//"/from_mesh/mesh_shape/element_type", &
                                pressure_mesh_element_type, &
                                default = "lagranian")
               
               if (trim(pressure_mesh_element_type) == "bubble") then
                  FLExit("For CG Pressure cannot test the continuity equation with CV if the pressure mesh has element type bubble")
               end if
               
               if (trim(pressure_mesh_element_type) == "trace") then
                  FLExit("For CG Pressure cannot test the continuity equation with CV if the pressure mesh has element type trace")
               end if
               
            end if
            
            ! Check that each particle phase has a PROGNOSTIC PhaseVolumeFraction field.
            ! The fluid phase cannot have a prognostic PhaseVolumeFraction as this is not always valid.
            ! For example, since we do not include the Density in the advection-diffusion equation for the PhaseVolumeFraction,
            ! solving this equation for the compressible fluid phase would not be correct. The particle phases on the other hand
            ! are always incompressible where the density is constant.
            if((have_option("/material_phase["//int2str(i)//"]/multiphase_properties/particle_diameter") .or. have_option("/material_phase["//int2str(i)//"]/multiphase_properties/particle_dia_use_scalar_field")) .and. &
               .not.(have_option("/material_phase["//int2str(i)//"]/scalar_field::PhaseVolumeFraction/prognostic") .or. &
               have_option("/material_phase["//int2str(i)//"]/scalar_field::PhaseVolumeFraction/prescribed"))) then
               FLExit("All particle phases must have a prognostic/prescribed PhaseVolumeFraction field. The diagnostic PhaseVolumeFraction field should always be in the continuous/fluid phase.")
            end if
            
         end do

         ewrite(1,*) 'Finished checking momentum discretisation options'

      end subroutine momentum_equation_check_options

   end module momentum_equation
