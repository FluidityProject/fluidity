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
!    C.Pain@Imperial.ac.uk
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

    use fields
    use state_module
    use spud
    use fldebug
    use momentum_cg
    use divergence_matrix_cv
    use divergence_matrix_cg
    use momentum_dg
    use assemble_cmc
    use field_priority_lists
    use momentum_diagnostic_fields, only: calculate_momentum_diagnostics
    use field_options
    use compressible_projection
    use boundary_conditions
    use boundary_conditions_from_options
    use sparse_matrices_fields
    use sparse_tools
    use free_surface_module
    use solvers
    use full_projection
    use petsc_solve_state_module
    use Profiler
    use geostrophic_pressure
    use hydrostatic_pressure
    use vertical_balance_pressure
    use oceansurfaceforcing
    use drag_module
    use parallel_tools
    use linked_lists
    use sparsity_patterns_meshes
    use state_matrices_module
    use vtk_interfaces
    use rotated_boundary_conditions
    use Weak_BCs
    use reduced_model_runtime
    use state_fields_module

    implicit none

    private
    public :: momentum_loop, momentum_equation_check_options
  contains

    subroutine momentum_loop(state, at_first_timestep, timestep, POD_state)
      !!< Construct and solve the momentum and continuity equations using
      !!< a continuous galerkin discretisation.

      ! an array of buckets full of fields
      ! the whole array is needed for the sake of multimaterial assembly
      type(state_type), dimension(:), intent(inout) :: state
      logical, intent(in) :: at_first_timestep
      type(state_type), dimension(:) :: POD_state
      integer, intent(in) :: timestep

      type(vector_field), pointer :: u
      integer :: istate, stat

      ewrite(1,*) 'Entering momentum_loop'

      ! this loop is just about the limit of multiphase support in this
      ! version of the momentum solve so far!
      call profiler_tic("momentum_loop")
      state_loop: do istate = 1, size(state)

        ! get the velocity
        u=>extract_vector_field(state(istate), "Velocity", stat)
        ! if there's no velocity then cycle
        if(stat/=0) cycle  
        ! if this is an aliased velocity then cycle
        if(aliased(u)) cycle
        ! if the velocity isn't prognostic then cycle
        if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

        if(have_option(trim(u%option_path)//"/prognostic/spatial_discretisation&
                                  &/continuous_galerkin").or.&
           have_option(trim(u%option_path)//"/prognostic/spatial_discretisation&
                                  &/discontinuous_galerkin")) then

           call profiler_tic("momentum_diagnostics")
           call calculate_momentum_diagnostics(state, istate)
           call profiler_toc("momentum_diagnostics")

           call profiler_tic("momentum_solve")
           call solve_momentum(state, istate, at_first_timestep, timestep, POD_state)
           call profiler_toc("momentum_solve")

        end if

      end do state_loop
      call profiler_toc("momentum_loop")

    end subroutine momentum_loop

    subroutine solve_momentum(state, istate, at_first_timestep, timestep, POD_state)
      !!< Construct and solve the momentum and continuity equations.
      ! an array of buckets full of fields
      ! the whole array is needed for the sake of multimaterial assembly
      type(state_type), dimension(:), intent(inout) :: state
      logical, intent(in) :: at_first_timestep
      integer, intent(in) :: istate, timestep
      type(state_type), dimension(:), intent(inout) :: POD_state

      ! the pressure gradient matrix (extracted from state)
      type(block_csr_matrix), pointer :: ct_m
      ! the pressure projection matrix (extracted from state)
      type(csr_matrix), pointer :: cmc_m
      ! do we want to form the ct or cmc matrices? are we solving a poisson pressure equation?
      logical :: get_ct_m, get_cmc_m, poisson_p

      ! matrix sparsity patterns for the matrices we allocate locally
      type(csr_sparsity), pointer :: u_sparsity

      ! locally allocated matrices:
      ! momentum lhs
      type(petsc_csr_matrix), target :: big_m
      ! Pointer to matrix for full projection solve:
      type(petsc_csr_matrix), pointer :: inner_m
      ! Pointer to preconditioner matrix for full projection solve:
      type(csr_matrix), pointer :: full_projection_preconditioner
      ! Are we going to form the Diagonal Schur complement preconditioner?
      logical :: get_diag_schur
      ! Do we need the scaled pressure mass matrix?
      logical :: get_scaled_pressure_mass_matrix
      ! Do we need an auxiliary matrix for full_projection solve?
      logical :: assemble_schur_auxiliary_matrix
      type(csr_sparsity), pointer :: schur_auxiliary_matrix_sparsity
      type(csr_matrix) :: schur_auxiliary_matrix
      ! Scaled pressure mass matrix - used for preconditioning full projection solve:
      type(csr_matrix), target :: scaled_pressure_mass_matrix
      type(csr_sparsity), pointer :: scaled_pressure_mass_matrix_sparsity      
      ! compressible pressure gradient operator/left hand matrix of cmc
      type(block_csr_matrix), pointer :: ctp_m
      ! the lumped mass matrix (may vary per component as absorption could be included)
      type(vector_field) :: inverse_masslump
      ! mass matrix
      type(petsc_csr_matrix), target :: mass
      ! for DG:
      type(block_csr_matrix) :: inverse_mass

      ! momentum rhs
      type(vector_field) :: mom_rhs
      ! projection rhs
      type(scalar_field) :: ct_rhs, projec_rhs, kmk_rhs, compress_projec_rhs, poisson_rhs

      ! change in pressure
      type(scalar_field) :: delta_p
      ! change in velocity
      type(vector_field) :: delta_u

      ! a dummy pressure field
      type(scalar_field), pointer :: dummyscalar, dummydensity, dummypressure

      ! pressure and density
      type(scalar_field), pointer :: p, density
      type(mesh_type), pointer :: p_mesh
      ! velocity
      type(vector_field), pointer :: u, x

      ! do we want to use the compressible projection method?
      logical :: use_compressible_projection
      ! are we doing a full schur solve?
      logical :: full_schur
      ! are we lumping mass or assuming consistent mass?
      logical :: lump_mass
      ! pressure gradient matrix using cv or cg?
      logical :: cv_pressure, cg_pressure

      ! the timestep
      real :: dt
      
      ! with free surface pressures are at integer time levels
      ! and we apply a theta weighting to the pressure gradient term
      ! also the continuity equation is considered at time n+theta
      ! instead of at the end of the timestep      
      real :: theta_pg
      ! in this case p_theta=theta*p+(1-theta)*old_p
      type(scalar_field), pointer :: old_p, p_theta
      type(vector_field), pointer :: old_u
      ! all of this only applies if use_theta_pg .eqv. .true.
      ! without a free surface, or with a free surface and theta==1
      ! use_theta_pg .eqv. .false. and p_theta => p
      logical :: use_theta_pg
      
      ! what is the equation type?
      character(len=FIELD_NAME_LEN) :: equation_type, poisson_scheme, schur_scheme, pressure_pmat

      integer :: i, stat

      logical :: dg
      logical :: apply_kmk, assemble_kmk
      logical :: have_viscosity, stress_form, have_coriolis, diagonal
      logical :: pressure_debugging_vtus
      !! True if the momentum equation should be solved with the reduced model.
      logical :: reduced_model
      
      ! the list of stiff nodes
      ! this is saved because the list is only formed when cmc is assembled, which
      ! isn't necessarily every time this subroutine is called but the list is
      ! still needed to fix the rhs (applying the fix to cmc itself wipes out the
      ! information that would be required to recompile the list)
      type(ilist), save :: stiff_nodes_list
      
      ! increased each call to momentum equation, used as index for pressure debugging vtus
      integer,save :: pdv_count=-1

      !!for reduced model
      type(vector_field), pointer :: snapmean_velocity
      type(scalar_field), pointer :: snapmean_pressure
      integer :: d


      ewrite(1,*) 'Entering solve_momentum'
      
      ! get the velocity
      u=>extract_vector_field(state(istate), "Velocity")
      do i = 1, u%dim
        ewrite_minmax(u%val(i)%ptr(:))
      end do

      dg=have_option(trim(u%option_path)//&
                          "/prognostic/spatial_discretisation&
                          &/discontinuous_galerkin")

      x=>extract_vector_field(state(istate), "Coordinate")

      ! get some velocity options:
      ! are we lumping the mass matrix?
      lump_mass = have_option(trim(u%option_path)//&
                          "/prognostic/spatial_discretisation&
                          &/continuous_galerkin/mass_terms&
                          &/lump_mass_matrix").or.&
                  have_option(trim(u%option_path)//&
                          "/prognostic/spatial_discretisation&
                          &/discontinuous_galerkin/mass_terms&
                          &/lump_mass_matrix")
      ! here is where we try to decide how big big_m should be
      have_viscosity=have_option(trim(u%option_path)//&
          &"/prognostic/tensor_field::Viscosity")
      ! the following should include a dg option when a stress form version gets implemented
      stress_form=have_option(trim(u%option_path)//&
          &"/prognostic/spatial_discretisation/continuous_galerkin&
          &/stress_terms/stress_form") 
      have_coriolis = have_option("/physical_parameters/coriolis")
      diagonal = .not.have_coriolis.and.(.not.(have_viscosity.and.stress_form))

      reduced_model= have_option("/reduced_model/execute_reduced_model")

      ! get the pressure
      p=>extract_scalar_field(state(istate), "Pressure", stat)
      if(stat/=0) then 
        ! allocate a dummy scalar field in case we have no pressure
        allocate(dummypressure)
        call allocate(dummypressure, u%mesh, "DummyPressure", field_type=FIELD_TYPE_CONSTANT)
        call zero(dummypressure)
        dummypressure%option_path=""
        p => dummypressure
        p_mesh => u%mesh
        old_p => dummypressure
        nullify(ct_m)
        get_ct_m=.false.
        nullify(cmc_m)
        get_cmc_m=.false.
      else
        p_mesh => p%mesh
        nullify(dummypressure)
        old_p => extract_scalar_field(state, "OldPressure", stat)
        if(stat/=0) old_p => p
        call profiler_tic(p, "assembly")
        ! get the pressure gradient matrix
        ct_m => get_velocity_divergence_matrix(state, get_ct=get_ct_m)
      
        ! get the pressure poisson matrix
        cmc_m => get_pressure_poisson_matrix(state, get_cmc=get_cmc_m)
        call profiler_toc(p, "assembly")
      end if
      ewrite_minmax(p%val)
            
      allocate(dummydensity)
      call allocate(dummydensity, u%mesh, "DummyDensity", field_type=FIELD_TYPE_CONSTANT)
      call set(dummydensity, 1.0)
      dummydensity%option_path=""

      allocate(dummyscalar)
      call allocate(dummyscalar, u%mesh, "DummyScalar", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummyscalar)
      dummyscalar%option_path=""


      call get_option(trim(u%option_path)//"/prognostic/equation[0]/name", &
                      equation_type)

      select case(equation_type)
      case("LinearMomentum")
        density=>extract_scalar_field(state(istate), "Density")
        get_cmc_m = get_cmc_m .or. (.not.constant_field(density))
      case("Boussinesq")
        density=>dummydensity
      case("Drainage")
        density=>dummyscalar
      case default
        FLAbort("Unknown equation type for velocity")
      end select
      ewrite_minmax(density%val)

      ! get some pressure options:
      ! are we using a compressible projection?
      use_compressible_projection = have_option(trim(p%option_path)//&
                                    "/prognostic/scheme&
                                    &/use_compressible_projection_method")
      
      get_cmc_m = get_cmc_m .or. &
                  have_option(trim(p%option_path)//&
                  "/prognostic/scheme/update_discretised_equation") .or. &
                  use_compressible_projection
      get_ct_m = get_ct_m .or. &
                  have_option(trim(p%option_path)//&
                  "/prognostic/scheme/update_discretised_equation")
                  
      pressure_debugging_vtus = have_option(trim(p%option_path)// &
                   "/prognostic/output/debugging_vtus")
      if (pressure_debugging_vtus) then
         pdv_count = pdv_count+1
      end if

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
      
      ! If we are using the reduced model then there is no pressure projection.
      if (reduced_model) get_cmc_m=.false.

      get_diag_schur = .false.
      get_scaled_pressure_mass_matrix = .false.
      assemble_schur_auxiliary_matrix = .false.

      full_schur = have_option(trim(p%option_path)//&
          &"/prognostic/scheme&
          &/use_projection_method/full_schur_complement")
      if(full_schur) then
         ! Check to see whether pressure cmc_m preconditioning matrix is needed:
         call get_option(trim(p%option_path)//&
              "/prognostic/scheme/use_projection_method&
              &/full_schur_complement/preconditioner_matrix[0]/name", pressure_pmat)
         select case(pressure_pmat)
         case("LumpedSchurComplement")
            full_projection_preconditioner => cmc_m
         case("DiagonalSchurComplement")
            if(.not.poisson_p) get_cmc_m = .false.
            get_diag_schur = .true.
            full_projection_preconditioner => cmc_m
         case("ScaledPressureMassMatrix")
            if(.not.poisson_p) get_cmc_m = .false.
            get_scaled_pressure_mass_matrix = .true.
            full_projection_preconditioner => scaled_pressure_mass_matrix
         case("NoPreconditionerMatrix")
            if(.not.poisson_p) get_cmc_m = .false.
            full_projection_preconditioner => cmc_m
         case default
            FLAbort("Unknown Matrix Type for Full_Projection")
         end select

         ! Decide on configuration of inner_m for full_projection solve:
         call get_option(trim(p%option_path)//&
              "/prognostic/scheme/use_projection_method&
              &/full_schur_complement/inner_matrix[0]/name", schur_scheme)
         select case(schur_scheme)
         case("FullMassMatrix")
            inner_m => mass
         case("FullMomentumMatrix")
            inner_m => big_m
         case default
            FLAbort("Unknown Matrix Type for Full_Projection")
         end select
      end if

      ! are we getting the pressure gradient matrix using control volumes?
      cv_pressure = (have_option(trim(p%option_path)//&
                        "/prognostic/spatial_discretisation/control_volumes"))
      ! or using cg (we do this in every case of not having a control volume
      ! option so that prescribed pressures will work as well)
      cg_pressure = (.not.cv_pressure)
      
      call get_option("/timestepping/timestep", dt)

      if (has_boundary_condition(u, "free_surface").or.use_compressible_projection) then
         ! with free surface pressures are at integer time levels
         ! and we apply a theta weighting to the pressure gradient term
         call get_option( trim(u%option_path)//'/prognostic/temporal_discretisation/theta', &
            theta_pg, default=1.0)
         use_theta_pg= (theta_pg/=1.0)
         ewrite(2,*) "Continuity equation and pressure gradient are evaluated at n+theta_pg"
         ewrite(2,*) "theta_pg: ", theta_pg
      else
         ! pressures are as usual staggered in time with the velocities
         use_theta_pg=.false.
      end if
            
      if (use_theta_pg) then
        allocate(p_theta)
        call allocate(p_theta, p%mesh, "PressureTheta")
        
        ! p_theta = theta*p + (1-theta)*old_p
        call set(p_theta, p, old_p, theta_pg)
        p_theta%option_path=p%option_path ! use p's solver options
        old_u => extract_vector_field(state(istate), "OldVelocity")
        if (old_u%aliased) then
          ! in the case of one non-linear iteration, there is no OldVelocity,
          ! it's just aliased to Velocity, therefore we make temp. version
          allocate(old_u)
          ! give it a distinct name, so we know to deallocate it
          call allocate(old_u, u%dim, u%mesh, "TempOldVelocity")
          call set(old_u, u)
        end if
      else
        p_theta => p
        theta_pg=1.0
      end if
      
      call profiler_tic(u, "assembly")
      ! allocation of big_m
      if(dg) then
        call allocate_big_m_dg(state(istate), big_m, u)
      else
        ! create a sparsity if necessary or pull it from state:
        u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
        ! and then allocate
        call allocate(big_m, u_sparsity, (/u%dim, u%dim/), diagonal=diagonal,&
            & name="BIG_m")
      end if

      call zero(big_m)
      if(get_ct_m) call zero(ct_m)

      ! allocate the momentum rhs
      call allocate(mom_rhs, u%dim, u%mesh, "MomentumRHS")
      call zero(mom_rhs)
      ! allocate the ct rhs
      call allocate(ct_rhs, p_mesh, "DivergenceRHS")
      call zero(ct_rhs)
      call profiler_toc(u, "assembly")
      
      if(has_scalar_field(state(istate), gp_name)) then
        call calculate_geostrophic_pressure_options(state(istate))
      end if

      ! assemble the momentum equation
      call profiler_tic(u, "assembly")
      if(dg) then
         call construct_momentum_dg(u, p, density, x, &
              big_m, mom_rhs, state(istate), &
              inverse_masslump=inverse_masslump, &
              inverse_mass=inverse_mass, &
              cg_pressure=cg_pressure)
        
        if(has_scalar_field(state(istate), gp_name)) then
          call subtract_geostrophic_pressure_gradient(mom_rhs, state(istate))
        end if
      else
        call construct_momentum_cg(u, p, density, x, &
             big_m, mom_rhs, ct_m, &
             ct_rhs, mass, inverse_masslump, &
             state(istate), &
             assemble_ct_matrix=get_ct_m, &
             cg_pressure=cg_pressure)
      end if
      call profiler_toc(u, "assembly")
      
      if(has_scalar_field(state(istate), hp_name)) then
        call calculate_hydrostatic_pressure(state(istate))
        call subtract_hydrostatic_pressure_gradient(mom_rhs, state(istate))
      end if
      
      if(has_scalar_field(state(istate), vbp_name)) then
        call calculate_vertical_balance_pressure(state(istate))
        call subtract_vertical_balance_pressure_gradient(mom_rhs, state(istate))
      end if
      
      call profiler_tic(u, "assembly")
      if (has_boundary_condition(u, "wind_forcing")) then
        call wind_forcing(state(istate), mom_rhs)
      end if
      
      if (has_boundary_condition(u, "drag")) then
        call drag_surface(big_m, mom_rhs, state(istate), density)
      end if

      ! near wall treatment
      if ( has_boundary_condition(u, "near_wall_treatment") .or. &
            has_boundary_condition(u, "log_law_of_wall")) then
         call wall_functions(big_m, mom_rhs, state(istate))
      end if
      call profiler_toc(u, "assembly")

      call profiler_tic(p, "assembly")
      if(cv_pressure) then
        call assemble_divergence_matrix_cv(ct_m, state(istate), ct_rhs=ct_rhs, &
                                           test_mesh=p%mesh, field=u, get_ct=get_ct_m)
      end if

      ! At the moment cg does its own ct assembly. We might change this in
      ! the future.
      if(cg_pressure.and.dg) then
        call  assemble_divergence_matrix_cg(CT_m, state(istate), ct_rhs=ct_rhs, & 
             test_mesh=p%mesh, field=u, get_ct=get_ct_m)
      end if
      call profiler_toc(p, "assembly")

      call profiler_tic(u, "assembly")
      if (have_rotated_bcs(u)) then
        ! rotates big_m, rhs and the velocity field at strong, surface_aligned dirichlet bcs
        call rotate_momentum_equation(big_m, mom_rhs, u, state(istate))
        if (get_ct_m) call rotate_ct_m(ct_m, u)
      end if    
      
      call apply_dirichlet_conditions(big_m, mom_rhs, u, dt)
      call profiler_toc(u, "assembly")

      if (associated(ct_m)) then
        do i = 1, ct_m%blocks(2)
          ewrite_minmax(ct_m%val(1,i)%ptr(:))
        end do
      end if

      ! do we want to solve for pressure?
      call profiler_tic(p, "assembly")
      if(have_option(trim(p%option_path)//"/prognostic").and..not.reduced_model) then
        if(use_compressible_projection) then
          allocate(ctp_m)
          call allocate(ctp_m, ct_m%sparsity, (/1, u%dim/), name="CTP_m")
          if(cv_pressure) then
            call assemble_compressible_divergence_matrix_cv(ctp_m, state, ct_rhs)
          else if(cg_pressure) then
            call assemble_compressible_divergence_matrix_cg(ctp_m, state, ct_rhs)
          else
            FLAbort("Unknown pressure discretisation for compressible projection.")
          end if
          if (have_rotated_bcs(u)) then
            call rotate_ct_m(ctp_m, u)
          end if
        else
          ctp_m=>ct_m
        end if
        do i = 1, ctp_m%blocks(2)
          ewrite_minmax(ctp_m%val(1,i)%ptr(:))
        end do
        ewrite_minmax(ct_rhs%val)

        ! Allocate the RHS
        call allocate(projec_rhs, p%mesh, "ProjectionRHS")
        call zero(projec_rhs)
        call allocate(kmk_rhs, p%mesh, "KMKRHS")
        call zero(kmk_rhs)

        ! Decide whether or not to form kmk stabilization matrix:
        apply_kmk = (continuity(p_mesh) >= 0 .and. p_mesh%shape%degree == 1 .and. p_mesh%shape%numbering%family == FAMILY_SIMPLEX .and. &
              & continuity(u%mesh) >= 0 .and. u%mesh%shape%degree == 1 .and. u%mesh%shape%numbering%family == FAMILY_SIMPLEX .and. &
              & .not. have_option(trim(p%option_path) // &
              & "/prognostic/spatial_discretisation/continuous_galerkin/remove_stabilisation_term") .and. &
              & .not. cv_pressure .and.&
              & .not. reduced_model)
        assemble_kmk= apply_kmk .and. &
                  ((.not. has_csr_matrix(state(istate), "PressureStabilisationMatrix")) .or. &
                  have_option(trim(p%option_path)// &
                  "/prognostic/scheme/update_discretised_equation") .or. &
                  have_option("/mesh_adaptivity/mesh_movement"))

            
        ! Assemble kmk stabilization matrix if required:
        if(assemble_kmk) then
           ewrite(2,*) "Assembling P1-P1 stabilisation"
           call assemble_kmk_matrix(state(istate), p%mesh, x, theta_pg)
        end if

        if(full_schur) then
           ! Decide whether we need to assemble an auxiliary matrix for full_projection solve:
           if(apply_kmk) assemble_schur_auxiliary_matrix = .true.
           if (has_boundary_condition(u, "free_surface")) assemble_schur_auxiliary_matrix = .true.
           ! If schur_auxiliary_matrix is needed then assemble it:
           if(assemble_schur_auxiliary_matrix) then
              ! Get sparsity and assemble:
              ewrite(2,*) "Assembling auxiliary matrix for full_projection solve"
              schur_auxiliary_matrix_sparsity => get_csr_sparsity_secondorder(state(istate), p%mesh, u%mesh)
              call allocate(schur_auxiliary_matrix, schur_auxiliary_matrix_sparsity,&
                   name="schur_auxiliary_matrix")
              ! Initialize matrix:
              call zero(schur_auxiliary_matrix)
              if(apply_kmk) then
                 ewrite(2,*) "Adding kmk stabilisation matrix to full_projection auxiliary matrix"
                 call add_kmk_matrix(state(istate), schur_auxiliary_matrix)
              end if
              if (has_boundary_condition(u, "free_surface")) then
                 ewrite(2,*) "Adding free surface to full_projection auxiliary matrix"
                 call add_free_surface_to_cmc_projection(state(istate), &
                      schur_auxiliary_matrix, dt, theta_pg, get_cmc=.true., rhs=ct_rhs)
              end if
           end if
        end if

        ! assemble the C_{P}^{T} M^{-1} C matrix
        if(get_cmc_m) then
          call zero(cmc_m)

          if(dg.and.(.not.lump_mass)) then
            call assemble_cmc_dg(cmc_m, ctp_m, ct_m, inverse_mass)
          else
            call assemble_masslumped_cmc(cmc_m, ctp_m, inverse_masslump, ct_m)

            ! P1-P1 stabilisation
            if (apply_kmk) then
              ewrite(2,*) "Adding P1-P1 stabilisation matrix to cmc_m"
              call add_kmk_matrix(state(istate), cmc_m)
            end if
          end if
          
          if(have_option(trim(p%option_path)//"/prognostic/repair_stiff_nodes")) then
            call repair_stiff_nodes(cmc_m, stiff_nodes_list)
          end if
          
        end if
        
        if (has_boundary_condition(u, "free_surface")) then
          call add_free_surface_to_cmc_projection(state(istate), &
              cmc_m, dt, theta_pg, get_cmc=get_cmc_m, rhs=ct_rhs)
        end if
        
        ! do we want to get an initial guess at the pressure?
        if(poisson_p) then

          call allocate(poisson_rhs, p%mesh, "PoissonRHS")

          ! get the rhs for the poisson pressure equation...
          if(dg.and.(.not.lump_mass)) then
            call assemble_poisson_rhs_dg(poisson_rhs, ctp_m, inverse_mass, mom_rhs, ct_rhs, u, dt, theta_pg)
          else
            ! here we assume that we're using mass lumping if we're not using dg
            ! if this isn't true then this leads to inconsistent mass matrices in poisson_rhs and cmc_m
            ! but as we're only hoping to get a guesstimate of the pressure hopefully this won't be too
            ! bad.
            call assemble_masslumped_poisson_rhs(poisson_rhs, ctp_m, mom_rhs, ct_rhs, inverse_masslump, u, dt, theta_pg)
          end if
          
          if (has_boundary_condition(u, "free_surface")) then
            call add_free_surface_to_poisson_rhs(poisson_rhs, state(istate), dt, theta_pg)
          end if
          
          ! apply strong dirichlet conditions
          call apply_dirichlet_conditions(cmc_m, poisson_rhs, p)
          
          call impose_reference_pressure_node(cmc_m, poisson_rhs, trim(p%option_path))
          
          call profiler_toc(p, "assembly") ! don't include poisson solve
          call petsc_solve(p_theta, cmc_m, poisson_rhs, state(istate))
          call profiler_tic(p, "assembly")
            
          if (has_boundary_condition(u, "free_surface")) then
            ! use this as initial pressure guess, except at the free surface
            ! where we use the prescribed initial condition
            call copy_poisson_solution_to_interior(p_theta, p, old_p, u)
          end if
          
          if (pressure_debugging_vtus) then
             call vtk_write_fields("initial_poisson", pdv_count, x, p%mesh, &
                sfields=(/ p_theta, p, old_p /))
          end if
          
          ewrite_minmax(p_theta%val)
          
          call deallocate(poisson_rhs)

        end if ! poisson pressure solution

        if(get_diag_schur) then
           ! Assemble diagonal schur complement preconditioner:
           call assemble_diagonal_schur(cmc_m,u,inner_m,ctp_m,ct_m)
           ! P1-P1 stabilisation:
           if (apply_kmk) then
              ewrite(2,*) "Adding P1-P1 stabilisation to diagonal schur complement preconditioner matrix"
              call add_kmk_matrix(state(istate), cmc_m)
           end if
           if (has_boundary_condition(u, "free_surface")) then
              ewrite(2,*) "Adding free surface to diagonal schur complement preconditioner matrix"
              call add_free_surface_to_cmc_projection(state(istate), &
                   cmc_m, dt, theta_pg, get_cmc=.true.)
           end if
        end if

        if(get_scaled_pressure_mass_matrix) then
           ! Assemble scaled pressure mass matrix which will later be used as a 
           ! preconditioner in the full projection solve:
           ewrite(2,*) "Assembling scaled pressure mass matrix preconditioner"
           scaled_pressure_mass_matrix_sparsity => get_csr_sparsity_firstorder(state(istate), p%mesh, p%mesh)
           call allocate(scaled_pressure_mass_matrix, scaled_pressure_mass_matrix_sparsity,&
                name="scaled_pressure_mass_matrix")
           call assemble_scaled_pressure_mass_matrix(state(istate),scaled_pressure_mass_matrix)
        end if

        if (apply_kmk) then
           call add_kmk_rhs(state(istate), kmk_rhs, p_theta, dt)
        end if
        ewrite_minmax(kmk_rhs%val)

      end if ! prognostic pressure
      call profiler_toc(p, "assembly")
      
      if (.not.reduced_model) then
      
          ! allocate the momentum solution vector
          call profiler_tic(u, "assembly")
          call allocate(delta_u, u%dim, u%mesh, "DeltaU")
          delta_u%option_path = trim(u%option_path)

          if (associated(ct_m)) then
            ! add - ct_m^T*p to the rhs of the momentum eqn
            ! (delta_u is just used as dummy memory here)
            !
            ! despite multiplying pressure by a nonlocal operator
            ! a halo_update isn't necessary as this is just a rhs
            ! contribution
            call mult_T(delta_u, ct_m, p_theta)

            if (dg) then
              ! We have just poluted the halo rows of delta_u. This is incorrect
              ! in the dg case due to the non-local assembly system employed.
              call zero_non_owned(delta_u)
            end if

            ewrite(2,*) 'note that delta_u = ct_m^T*p at this stage'
            do i = 1, delta_u%dim
              ewrite_minmax(delta_u%val(i)%ptr(:))
            end do
            call addto(mom_rhs, delta_u)
          end if

          ! impose zero guess on change in u
          call zero(delta_u)

          ! impose any reference nodes on velocity
          call impose_reference_velocity_node(big_m, mom_rhs, trim(u%option_path))
       
          call profiler_toc(u, "assembly")

          ! solve for the change in velocity
          call petsc_solve(delta_u, big_m, mom_rhs, state(istate))
          do i = 1, u%dim
            ewrite_minmax(delta_u%val(i)%ptr(:))
          end do

          call profiler_tic(u, "assembly")
          ! apply change to velocity field
          call addto(u, delta_u, dt)
          do i = 1, u%dim
            ewrite_minmax(u%val(i)%ptr(:))
          end do

          call deallocate(delta_u)
          call profiler_toc(u, "assembly")

          ! do we want to solve for pressure?
          if(have_option(trim(p%option_path)//"/prognostic")) then
            call profiler_tic(p, "assembly")
            ! assemble the rhs
            ! if we are adding the P1-P1 stabilisation,
            ! this will have to have KMK * P added to it;
            !
            ! despite multiplying velocity by a nonlocal operator
            ! a halo_update isn't necessary as this is just a rhs
            ! contribution
            if (.not. use_theta_pg) then
              ! continuity is evaluated at the end of the time step
              call mult(projec_rhs, ctp_m, u)
            else
              ! evaluate continuity at n+theta
              ! compute theta*u+(1-theta)*old_u
              call allocate(delta_u, u%dim, u%mesh, "VelocityTheta")
              if (have_rotated_bcs(u)) then
                call rotate_velocity(old_u, state(istate))
              end if
              call set(delta_u, u, old_u, theta_pg)
              call mult(projec_rhs, ctp_m, delta_u)
              call deallocate(delta_u)
            end if
            call addto(projec_rhs, kmk_rhs)
            call scale(projec_rhs, -1.0)
            call addto(projec_rhs, ct_rhs)
            ewrite_minmax(projec_rhs%val(:))
            
            if(use_compressible_projection) then
              call allocate(compress_projec_rhs, p%mesh, "CompressibleProjectionRHS")
              
              if(cv_pressure) then
                call assemble_compressible_projection_cv(state, cmc_m, compress_projec_rhs, dt, theta_pg, &
                                                        get_cmc_m)
              else if(cg_pressure) then
                call assemble_compressible_projection_cg(state, cmc_m, compress_projec_rhs, dt, theta_pg, &
                                                        get_cmc_m)
              else
                FLAbort("Unknown pressure discretisation for compressible projection.")
              end if
              
              ewrite_minmax(compress_projec_rhs%val)
              ewrite_minmax(cmc_m%val)
              
              call addto(projec_rhs, compress_projec_rhs)

              call deallocate(compress_projec_rhs)
            end if

            ! apply strong dirichlet conditions
            ! we're solving for "delta_p"=theta_pg**2*dp*dt, where dp=p_final-p_current
            ! apply_dirichlet_condition however assumes we're solving for
            ! "acceleration" dp/dt, by providing dt=1/(dt*theta_pg**2) we get what we want
            call apply_dirichlet_conditions(cmc_m, projec_rhs, p, &
              dt=1.0/(dt*theta_pg**2))
            
            call impose_reference_pressure_node(cmc_m, projec_rhs, trim(p%option_path))
            
            ! allocate the change in pressure field
            call allocate(delta_p, p%mesh, "DeltaP")
            delta_p%option_path = trim(p%option_path)
            call zero(delta_p)

            if(have_option(trim(p%option_path)//"/prognostic/repair_stiff_nodes")) then
              call zero_stiff_nodes(projec_rhs, stiff_nodes_list)
            end if
            call profiler_toc(p, "assembly")
      
            ! solve for the change in pressure
            if(full_schur) then
               if(assemble_schur_auxiliary_matrix) then
                  call petsc_solve_full_projection(delta_p,ctp_m,inner_m,ct_m,projec_rhs, &
                       full_projection_preconditioner,schur_auxiliary_matrix)
               else
                  call petsc_solve_full_projection(delta_p,ctp_m,inner_m,ct_m,projec_rhs, &
                       full_projection_preconditioner)
               end if
            else
              call petsc_solve(delta_p, cmc_m, projec_rhs, state(istate))
            end if

            ewrite_minmax(delta_p%val)

            if (pressure_debugging_vtus) then
              ! writes out the pressure and velocity before the correction is added in
              ! (as the corrected fields are already available in the convergence files)
              call vtk_write_fields("pressure_correction", pdv_count, x, p%mesh, &
                  sfields=(/ delta_p, p, old_p, p_theta /))
              ! same thing but now on velocity mesh:
              call vtk_write_fields("velocity_before_correction", pdv_count, x, u%mesh, &
                  sfields=(/ delta_p, p, old_p, p_theta /), vfields=(/ u /))
            end if
            
            call profiler_tic(p, "assembly")
            if (use_theta_pg) then
              ! we've solved theta_pg**2*dt*dp, in the velocity correction
              ! however we need theta_pg*dt*dp
              call scale(delta_p, 1.0/theta_pg)
            end if        
            
            ! add the change in pressure to the pressure
            ! (if .not. use_theta_pg then theta_pg is 1.0)
            call addto(p, delta_p, scale=1.0/(theta_pg*dt))
            ewrite_minmax(p%val)
            
            if(use_compressible_projection) then
              call update_compressible_density(state)
            end if
            call profiler_toc(p, "assembly")

            call profiler_tic(u, "assembly")
            ! correct velocity according to new delta_p
            if(full_schur) then
              call correct_velocity_cg(u, inner_m, ct_m, delta_p)
            elseif(lump_mass) then
              call correct_masslumped_velocity(u, inverse_masslump, ct_m, delta_p)
            elseif(dg) then
              call correct_velocity_dg(u, inverse_mass, ct_m, delta_p)
            else
              FLAbort("Don't know how to correct the velocity.")
            end if
            call profiler_toc(u, "assembly")
            
            call deallocate(kmk_rhs)
            call deallocate(projec_rhs)
            call deallocate(delta_p)

            if(assemble_schur_auxiliary_matrix) then
               ! Deallocate schur_auxiliary_matrix:
               call deallocate(schur_auxiliary_matrix)
            end if

            if(get_scaled_pressure_mass_matrix) then
               ! Deallocate scaled pressure mass matrix:
               call deallocate(scaled_pressure_mass_matrix)
            end if

            if(use_compressible_projection) then
              call deallocate(ctp_m)
              deallocate(ctp_m)
            end if

          end if

      else ! Reduced model version
  
          ! allocate the change in pressure field
          call allocate(delta_p, p%mesh, "DeltaP")
          delta_p%option_path = trim(p%option_path)
          call zero(delta_p)
          
          call allocate(delta_u, u%dim, u%mesh, "DeltaU")
          delta_u%option_path = trim(u%option_path)
          call zero(delta_u)

          call solve_momentum_reduced(delta_u, delta_p, big_m, mom_rhs, ct_m, ct_rhs, timestep, POD_state) 

          snapmean_velocity=>extract_vector_field(POD_state(1), "SnapmeanVelocity")
          snapmean_pressure=>extract_scalar_field(POD_state(1), "SnapmeanPressure")

          if(timestep==1)then
              do d=1,snapmean_velocity%dim
                u%val(d)%ptr=snapmean_velocity%val(d)%ptr
              enddo
              p%val=snapmean_pressure%val

              call addto(u, delta_u, dt)
              call addto(p, delta_p, dt)
          else
              call addto(u, delta_u, dt)
              call addto(p, delta_p, dt)
          endif

          call deallocate(delta_p)
          call deallocate(delta_u)
          
      end if ! end of if reduced model
      
      call profiler_tic(u, "assembly")
      if (have_rotated_bcs(u)) then
        call rotate_velocity_back(u, state(istate))
      end if
      call profiler_toc(u, "assembly")
      
      if(dg)then
        if (lump_mass) then
          call deallocate(inverse_masslump)
        else
          call deallocate(inverse_mass)
        end if
      else
        call deallocate_cg_mass(mass, inverse_masslump)
      end if
      
      if (use_theta_pg) then
        call deallocate(p_theta)
        deallocate(p_theta)
        if (old_u%name=="TempOldVelocity") then
          call deallocate(old_u)
          deallocate(old_u)
        else if (have_rotated_bcs(u)) then
          call rotate_velocity_back(old_u, state(istate))
        end if
      end if

      call deallocate(mom_rhs)
      call deallocate(ct_rhs)
      if (associated(dummypressure)) then
        call deallocate(dummypressure)
        deallocate(dummypressure)
      end if
      call deallocate(dummydensity)
      deallocate(dummydensity)
      call deallocate(dummyscalar)
      deallocate(dummyscalar)
      call deallocate(big_m)
 
    end subroutine solve_momentum
      
    subroutine momentum_equation_check_options

      integer :: i, nmat
      character(len=FIELD_NAME_LEN) :: schur_scheme
      character(len=FIELD_NAME_LEN) :: schur_preconditioner
      
      ewrite(1,*) 'Checking momentum discretisation options'

      nmat = option_count("/material_phase")

      do i = 0, nmat-1

        if(have_option("/material_phase["//int2str(i)//&
                            "]/scalar_field::Pressure/prognostic&
                            &/reference_node").and.&
           have_option("/material_phase["//int2str(i)//&
                            "]/scalar_field::Pressure/prognostic&
                            &/solver/remove_null_space")) then
          FLExit("Can't set a pressure reference node and remove the null space.")
        end if

        if(have_option("/material_phase["//int2str(i)//&
                            "]/vector_field::Velocity/prognostic/reference_node")) then
          if((.not.(have_option("/material_phase["//int2str(i)//&
                              "]/vector_field::Velocity/prognostic&
                              &/spatial_discretisation/continuous_galerkin/mass_terms/exclude_mass_terms").and. &
                    have_option("/material_phase["//int2str(i)//&
                              "]/vector_field::Velocity/prognostic&
                              &/spatial_discretisation/continuous_galerkin/advection_terms/exclude_advection_terms"))).and. &
             (.not.(have_option("/material_phase["//int2str(i)//&
                              "]/vector_field::Velocity/prognostic&
                              &/spatial_discretisation/discontinuous_galerkin/mass_terms/exclude_mass_terms").and. &
                    have_option("/material_phase["//int2str(i)//&
                              "]/vector_field::Velocity/prognostic&
                              &/spatial_discretisation/discontinuous_galerkin/advection_scheme/none")))) then
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
                            "]/vector_field::Velocity/prognostic&
                            &/spatial_discretisation/continuous_galerkin")&
           .and.(.not.have_option("/material_phase["//int2str(i)//&
                            "]/vector_field::Velocity/prognostic&
                            &/spatial_discretisation/continuous_galerkin&
                            &/mass_terms/lump_mass_matrix"))) then

          if(have_option("/material_phase["//int2str(i)//&
                            "]/scalar_field::Pressure/prognostic&
                            &/scheme/use_projection_method")) then
            if(.not.have_option("/material_phase["//int2str(i)//&
                              "]/scalar_field::Pressure/prognostic&
                              &/scheme/use_projection_method&
                              &/full_schur_complement")) then
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
          else if(have_option("/material_phase["//int2str(i)//&
                            "]/scalar_field::Pressure/prognostic&
                            &/scheme/use_compressible_projection_method")) then
            ewrite(-1,*) "You must lump the velocity mass matrix with the"
            ewrite(-1,*) "compressible projection method."
            FLExit("Sorry.")
          end if

        end if

        if(have_option("/material_phase["//int2str(i)//&
                            "]/vector_field::Velocity/prognostic&
                            &/tensor_field::Viscosity/prescribed/value&
                            &/isotropic").or. &
           have_option("/material_phase["//int2str(i)//&
                            "]/vector_field::Velocity/prognostic&
                            &/tensor_field::Viscosity/prescribed/value&
                            &/diagonal")) then

          if(have_option("/material_phase["//int2str(i)//&
                            "]/vector_field::Velocity/prognostic&
                            &/spatial_discretisation/continuous_galerkin/stress_terms/stress_form")) then
            ewrite(-1,*) "You have selected stress form viscosity but have entered an isotropic or"
            ewrite(-1,*) "diagonal Viscosity tensor."
            ewrite(-1,*) "Zero off diagonal entries in the Viscosity tensor do not make physical"
            ewrite(-1,*) "sense when using stress form viscosity."
            ewrite(-1,*) "Use an tensor_form or anisotropic_symmetric Viscosity instead."
            FLExit("Use tensor_form or anisotropic_symmetric Viscosity.")
          end if

        end if

        if(have_option("/material_phase["//int2str(i)//&
             "]/vector_field::Velocity/prognostic&
             &/tensor_field::Viscosity/prescribed/value&
             &/anisotropic_symmetric").or.&
           have_option("/material_phase["//int2str(i)//&
             "]/vector_field::Velocity/prognostic&
             &/tensor_field::Viscosity/prescribed/value&
             &/anisotropic_asymmetric")) then

          if(have_option("/material_phase["//int2str(i)//&
               "]/scalar_field::Pressure/prognostic&                                                                                                                                              
               &/scheme/use_projection_method")) then

             if(have_option("/material_phase["//int2str(i)//&
                  "]/scalar_field::Pressure/prognostic&                                                                                                                                            
                  &/scheme/use_projection_method&                                                                                                                                                  
                  &/full_schur_complement")) then

                call get_option("/material_phase["//int2str(i)//&
                     "]/scalar_field::Pressure/prognostic&                                                                                                                                              
                     &/scheme/use_projection_method&                                                                                                                                                    
                     &/full_schur_complement/preconditioner_matrix[0]/name", schur_preconditioner)

                select case(schur_preconditioner)
                case("ScaledPressureMassMatrix")
                   ewrite(-1,*) "At present, the viscosity scaling for the pressure mass matrix is only"
                   ewrite(-1,*) "valid for isotropic viscosity tensors. Please use another preconditioner"
                   ewrite(-1,*) "for the Full Projection solve"
                   FLExit("Sorry!")
                end select

             end if

          end if

       end if

        if(have_option("/material_phase["//int2str(i)//&
                            "]/vector_field::Velocity/prognostic&
                            &/spatial_discretisation/discontinuous_galerkin")) then

          if(have_option("/material_phase["//int2str(i)//&
                            "]/scalar_field::Pressure/prognostic&
                            &/scheme/use_projection_method")) then
            if(have_option("/material_phase["//int2str(i)//&
                              "]/scalar_field::Pressure/prognostic&
                              &/scheme/use_projection_method&
                              &/full_schur_complement")) then
                              
              call get_option("/material_phase["//int2str(i)//&
                            "]/scalar_field::Pressure/prognostic&
                            &/scheme/use_projection_method&
                            &/full_schur_complement/inner_matrix[0]/name", schur_scheme)
              select case(schur_scheme)
              case("FullMassMatrix")
                FLExit("Can't do a full schur complement solve with dg velocity and a mass inner matrix.")
              end select
              
            end if
          end if

        end if

      end do

      ewrite(1,*) 'Finished checking momentum discretisation options'

    end subroutine momentum_equation_check_options
    
  end module momentum_equation
