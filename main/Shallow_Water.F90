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
  program shallow_water
    use advection_diffusion_dg
    use advection_diffusion_cg
    use linear_shallow_water
    use spud
    use signals
    use fields
    use state_module
    use FLDebug
    use populate_state_module
    use write_state_module
    use vtk_interfaces
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use hybridized_helmholtz
    use assemble_cmc
    use global_parameters, only: option_path_len, python_func_len, current_time, dt
    use adapt_state_prescribed_module
    use memory_diagnostics
    use reserve_state_module
    use boundary_conditions_from_options
      use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
    use iso_c_binding
    use mangle_options_tree
    use manifold_projections
    use adjoint_controls
#ifdef HAVE_ADJOINT
    use libadjoint_data_callbacks
    use shallow_water_adjoint_callbacks
    use libadjoint
    use adjoint_functional_evaluation
    use adjoint_python
    use adjoint_global_variables
    use shallow_water_adjoint_controls
    use adjoint_main_loop
    use forward_main_loop
#include "libadjoint/adj_fortran.h"
#endif
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  implicit none
#include "petsc_legacy.h"

    ! Interface blocks for the initialisation routines we need to call
    interface
      subroutine set_global_debug_level(n)
        integer, intent(in) :: n
      end subroutine set_global_debug_level

      subroutine python_init
      end subroutine python_init

      subroutine petscinitialize(s, i)
        character(len=*), intent(in) :: s
        integer, intent(out) :: i
      end subroutine petscinitialize
    end interface

    type(state_type), dimension(:), pointer :: state

    type(state_type), target :: matrices ! We collect all the cached
                                         ! matrices in this state so that
                                         ! the adjoint callbacks can use
                                         ! them
    real :: D0, g, theta, itheta
    logical :: exclude_velocity_advection, exclude_pressure_advection
    logical :: hybridized
   integer :: timestep, nonlinear_iterations
    integer :: ierr

    type(vector_field), pointer :: v_field
    !! Sparsity for matrices.
    type(csr_sparsity) :: ct_sparsity,u_sparsity,wave_sparsity

    !! Mass matrices
    type(csr_matrix) :: h_mass_mat
    type(block_csr_matrix) :: u_mass_mat
    !! Coriolis matrix
    type(block_csr_matrix) :: coriolis_mat
    !! inverse Coriolis matrix
    type(block_csr_matrix) :: inverse_coriolis_mat
    !! div matrix
    type(block_csr_matrix) :: div_mat
    !! Wave matrix
    type(csr_matrix) :: wave_mat
    !! U momentum matrix
    type(block_csr_matrix) :: big_mat
    character(len = OPTION_PATH_LEN) :: simulation_name
    logical :: adjoint
    integer, save :: dump_no=0
    real :: energy
#ifdef HAVE_ADJOINT
    ierr = adj_create_adjointer(adjointer)
    ! Register the data callbacks
    call adj_register_femtools_data_callbacks(adjointer)
    ! Register the operator callbacks
    call register_sw_operator_callbacks(adjointer)
#endif

#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init
    call read_command_line
    call mangle_options_tree_forward
    ! Establish signal handlers
    call initialise_signals()

    adjoint = have_option("/adjoint")
#ifndef HAVE_ADJOINT
    if (adjoint) then
      FLExit("Cannot run the adjoint model without having compiled fluidity --with-adjoint.")
    endif
#else
    call register_functional_callbacks()
    if (.not. adjoint) then
      ! disable the adjointer
      ierr = adj_deactivate_adjointer(adjointer) 
      call adj_chkierr(ierr)
    end if
#endif

    is_shallow_water=.true.
    timestep=0

    call populate_state(state)
    ! Read in any control variables
    call adjoint_load_controls(timestep, dt, state)
    call adjoint_register_initial_eta_condition(state)

    call insert_time_in_state(state)

    call allocate_and_insert_additional_fields(state(1))

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state)

    call get_parameters

    ! No support for multiphase or multimaterial at this stage.
    if (size(state)/=1) then
       FLExit("Multiple material_phases are not supported")
    end if

    call calculate_diagnostic_variables(state)
    call calculate_diagnostic_variables_new(state)
    call write_diagnostics(state, current_time, dt, timestep)

    ! Always output the initial conditions.
    call output_state(state)

    hybridized = .false.
    v_field => extract_vector_field(state(1),"Velocity")
    if(associated(v_field%mesh%shape%constraints)) then
       if(v_field%mesh%shape%constraints%type.ne.CONSTRAINT_NONE) hybridized =&
            & .true.
    end if

    if(hybridized) then
       !project velocity into div-conforming space
       v_field => extract_vector_field(state(1),"LocalVelocity")
       call project_to_constrained_space(state(1),v_field)
    end if

    if(hybridized) then
       call compute_energy_hybridized(state(1),energy)
    else
       call get_linear_energy(state(1),u_mass_mat,h_mass_mat,d0,g,energy)
    end if
    ewrite(2,*) 'Initial Energy:', energy

    ! Register the control variables to disk if desired
    call adjoint_write_controls(timestep, dt, state)
    timestep_loop: do
       timestep=timestep+1
       ewrite (1,*) "SW: start of timestep ", timestep, current_time

       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)

       ! evaluate prescribed fields at time = current_time+dt
       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time + (theta * dt))
       ! Read in any control variables
       call adjoint_load_controls(timestep, dt, state)

       if (has_vector_field(state(1), "VelocitySource")) then
          v_field => extract_vector_field(state(1), "VelocitySource")
          call project_cartesian_to_local(state(1), v_field)
       end if 

       call execute_timestep(state, 1, dt)

       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time + dt)

       call project_local_to_cartesian(state(1))
       call calculate_diagnostic_variables(state,&
            & exclude_nonrecalculated = .true.)
       call calculate_diagnostic_variables_new(state,&
            & exclude_nonrecalculated = .true.)
       call adjoint_register_timestep(timestep, dt, state)
       ! Save the control variables to disk if desired
       call adjoint_write_controls(timestep, dt, state)

       call advance_current_time(current_time, dt)
       if (simulation_completed(current_time, timestep)) exit timestep_loop

       if (do_write_state(current_time, timestep)) then
          call output_state(state)
       end if
  
#ifdef HAVE_ADJOINT
       call calculate_functional_values(timestep-1)
#endif
       call write_diagnostics(state,current_time, dt, timestep)

       !Update the variables
       if(hybridized) then
          call compute_energy_hybridized(state(1),energy)
       else
          call get_linear_energy(state(1),u_mass_mat,h_mass_mat,d0,g,energy)
       end if
       ewrite(2,*) 'Energy = ',energy
       
    end do timestep_loop

    if(hybridized) then
       call compute_energy_hybridized(state(1),energy)
    else
       call get_linear_energy(state(1),u_mass_mat,h_mass_mat,d0,g,energy)
    end if
    ewrite(2,*) 'Final Energy:', energy

    ! One last dump
    call output_state(state)
#ifdef HAVE_ADJOINT
    call calculate_functional_values(timestep-1)
#endif
    call write_diagnostics(state,current_time,dt,timestep)

    if(.not.hybridized) then
       call deallocate(h_mass_mat)
       call deallocate(u_mass_mat)
       call deallocate(coriolis_mat)
       call deallocate(inverse_coriolis_mat)
       call deallocate(div_mat)
       call deallocate(wave_mat)
       call deallocate(big_mat)
    end if
    call deallocate(state)
    call deallocate_transform_cache
    call deallocate_reserve_state
    call close_diagnostic_files
    call uninitialise_diagnostics
    ! Clean up registered diagnostics
    call destroy_registered_diagnostics 

    if (.not. adjoint) then
      call deallocate(matrices)
      call print_references(0)
    endif

#ifdef HAVE_ADJOINT
    if (adjoint) then
      if (have_option("/adjoint/debug/replay_forward_run")) then
        ! Let's run the forward model through libadjoint, too, for the craic
        ewrite(1,*) "Entering forward computation through libadjoint"
        call clear_options
        call read_command_line
        call mangle_options_tree_forward
        call populate_state(state)
        call allocate_and_insert_additional_fields(state(1))
        call check_diagnostic_dependencies(state)
        call compute_forward(state)
        call deallocate_transform_cache
        call deallocate_reserve_state
        call deallocate(state)
      end if

      ewrite(1,*) "Entering adjoint computation"
      call clear_options
      call read_command_line

      call mangle_options_tree_adjoint
      call populate_state(state)
      call allocate_and_insert_additional_fields(state(1), adjoint=.true.)
      call check_diagnostic_dependencies(state)
      call compute_matrix_transposes(matrices)

      dump_no = dump_no - 1

      call compute_adjoint(state, dump_no, shallow_water_adjoint_timestep_callback, c_loc(matrices))

      call deallocate_transform_cache
      call deallocate_reserve_state

      call deallocate(state)
      call deallocate(matrices)
      call print_references(0)
    else
      ewrite(1,*) "No adjoint specified, not entering adjoint computation"
    end if
#endif

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

#ifdef HAVE_ADJOINT
    ierr = adj_destroy_adjointer(adjointer)
    call adj_chkierr(ierr)
#endif

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif
  contains

    subroutine get_parameters()
      implicit none
      type(vector_field), pointer :: u, v_field, coord
      type(scalar_field), pointer :: eta
      type(vector_field) :: dummy_field

      !Get some parameters
      !gravity
      call get_option("/physical_parameters/gravity/magnitude", g)
      !theta
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/prognostic/temporal_discretisation/theta",theta)
      !itheta

      !D0
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/p&
           &rognostic/mean_layer_thickness",D0)

      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/timestep", dt)

      if(.not.hybridized) then
         call setup_wave_matrices(state(1),u_sparsity,&
              wave_sparsity,&
              ct_sparsity,&
              h_mass_mat,u_mass_mat,&
              coriolis_mat,inverse_coriolis_mat,&
              div_mat,wave_mat,big_mat, &
              dt,theta,D0,g)
      end if
      v_field => extract_vector_field(state(1), "Velocity")
      call project_cartesian_to_local(state(1), v_field)
      if (has_vector_field(state(1), "VelocitySource")) then
                v_field => extract_vector_field(state(1), "VelocitySource")
        call project_cartesian_to_local(state(1), v_field)
      end if

      call get_option("/timestepping/nonlinear_iterations"&
           &,nonlinear_iterations)
      exclude_pressure_advection = &
           have_option("/material_phase::Fluid/scalar_field::LayerThickness/pro&
           &gnostic/spatial_discretisation/continuous_galerkin/advection_terms&
           &/exclude_advection_terms")
      exclude_velocity_advection = &
           have_option("/material_phase::Fluid/vector_field::Velocity/prognost&
           &ic/spatial_discretisation/discontinuous_galerkin/advection_scheme/&
           &none")
      call get_option("/material_phase::Fluid/scalar_field::LayerThickness/pro&
           &gnostic/temporal_discretisation/relaxation",itheta)
         ! Geostrophic balanced initial condition, if required
      if(have_option("/material_phase::Fluid/vector_field::Velocity/prognostic&
           &/initial_condition::WholeMesh/balanced")) then
         if(hybridized) then
            call set_velocity_from_geostrophic_balance_hybridized(state(1))
         else
            call set_velocity_from_geostrophic_balance(state(1), &
                 div_mat, coriolis_mat, inverse_coriolis_mat)
         end if
         call adjoint_register_initial_u_condition(balanced=.true.)
      else
         call adjoint_register_initial_u_condition(balanced=.false.)
      end if

      if(.not.hybridized) then
         ! Set up the state of cached matrices
         call insert(matrices, u_mass_mat, "LocalVelocityMassMatrix")
         call insert(matrices, h_mass_mat, "LayerThicknessMassMatrix")
         call insert(matrices, coriolis_mat, "CoriolisMatrix")
         call insert(matrices, inverse_coriolis_mat,&
              & "InverseCoriolisMatrix")
         call insert(matrices, div_mat, "DivergenceMatrix")
         call insert(matrices, wave_mat, "WaveMatrix")
         call insert(matrices, big_mat, "InverseBigMatrix")
      end if
      ! Also save the velocity and pressure mesh and the dimension
      eta => extract_scalar_field(state, "LayerThickness")
      u => extract_vector_field(state, "Velocity")
      call insert(matrices, eta%mesh, "LayerThicknessMesh")
      ! Insert a dummy velocity field which will be used in the adjoint callbacks
      call allocate(dummy_field, u%dim, u%mesh, "VelocityDummy", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummy_field)
      call insert(matrices, dummy_field, "VelocityDummy")
      call deallocate(dummy_field)
      ! Insert a dummy local velocity field which will be used in the adjoint callbacks
      v_field => extract_vector_field(state(1), "LocalVelocity")
      call allocate(dummy_field, v_field%dim, v_field%mesh, "LocalVelocityDummy", field_type=FIELD_TYPE_CONSTANT)
      call zero(dummy_field)
      call insert(matrices, dummy_field, "LocalVelocityDummy")
      call deallocate(dummy_field)
      ! And don't forget Coordinate
      coord => extract_vector_field(state(1), "Coordinate")
      call insert(matrices, coord, "Coordinate")
      ! And the CartesianVelocityMassMatrix
      call assemble_cartesian_velocity_mass_matrix(state(1))
      call insert(matrices, extract_block_csr_matrix(state(1), "CartesianVelocityMassMatrix"), "CartesianVelocityMassMatrix")
      ! And the VelocitySource and LayerThicknessSource
      if (has_vector_field(state(1), "VelocitySource")) then
        call insert(matrices, extract_vector_field(state(1), "VelocitySource"), "VelocitySource")
      end if
      if (has_scalar_field(state(1), "LayerThicknessSource")) then
        call insert(matrices, extract_scalar_field(state(1), "LayerThicknessSource"), "LayerThicknessSource")
      end if  
    end subroutine get_parameters

    subroutine assemble_cartesian_velocity_mass_matrix(state)
      type(state_type), intent(inout) :: state
      type(vector_field), pointer :: positions, velocity
      type(block_csr_matrix) :: mass_matrix
      type(csr_sparsity), pointer :: sparsity
      integer :: ele

      positions => extract_vector_field(state, "Coordinate")
      velocity => extract_vector_field(state, "Velocity")
      sparsity => get_csr_sparsity_firstorder(state, velocity%mesh, velocity%mesh)
      call allocate(mass_matrix, sparsity, (/velocity%dim, velocity%dim/), name="CartesianVelocityMassMatrix", diagonal=.true., equal_diagonal_blocks=.true.)
      call zero(mass_matrix)

      do ele = 1, ele_count(positions)
        call assemble_cartesian_velocity_mass_matrix_ele(mass_matrix, velocity%mesh, positions, ele)
      end do

      ! Insert mass matrix into state
      call insert(state, mass_matrix, "CartesianVelocityMassMatrix")
      call deallocate(mass_matrix)

    end subroutine assemble_cartesian_velocity_mass_matrix

    subroutine assemble_cartesian_velocity_mass_matrix_ele(mass_matrix, mesh, positions, ele)
      type(block_csr_matrix), intent(inout) :: mass_matrix
      type(mesh_type), intent(in) :: mesh
      type(vector_field), intent(in) :: positions
      integer, intent(in) :: ele
      real, dimension(ele_loc(mesh, ele), ele_loc(mesh, ele)) :: little_mass_matrix
      real, dimension(ele_ngi(mesh, ele)) :: detwei

      call transform_to_physical(positions, ele, detwei=detwei)
      little_mass_matrix = shape_shape(ele_shape(mesh, ele), ele_shape(mesh, ele), detwei)
      call addto(mass_matrix, 1, 1, ele_nodes(mesh, ele), ele_nodes(mesh, ele), little_mass_matrix) ! we only add it once because equal_diagonal_blocks=.true.
    end subroutine assemble_cartesian_velocity_mass_matrix_ele
    
    subroutine insert_time_in_state(state)
      type(state_type), dimension(:), intent(inout) :: state

      type(scalar_field) :: aux_sfield
      type(mesh_type), pointer :: x_mesh
      real :: current_time

      ! Disgusting and vomitous hack to ensure that time is output in
      ! vtu files.
      x_mesh => extract_mesh(state, "CoordinateMesh")
      call allocate(aux_sfield, x_mesh, "Time", field_type=FIELD_TYPE_CONSTANT)
      call get_option("/timestepping/current_time", current_time)
      call set(aux_sfield, current_time)
      aux_sfield%option_path = ""
      call insert(state, aux_sfield, trim(aux_sfield%name))
      call deallocate(aux_sfield)

    end subroutine insert_time_in_state

    subroutine execute_timestep(state, istate, dt)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state
      integer, intent(in) :: istate
      real, intent(in) :: dt

      !! Layer thickness
      type(scalar_field), pointer :: D, d_src
      !! velocity
      type(vector_field), pointer :: U
      !! Source term
      type(vector_field), pointer :: source

      !!Intermediate fields
      type(scalar_field) :: d_rhs, delta_d, old_d, md_src
      type(vector_field) :: u_rhs, delta_u, advecting_u, old_u
      type(scalar_field) :: velocity_cpt, old_velocity_cpt
      integer :: dim, nit, d1
      real :: energy
      logical :: have_source

      !Pull the fields out of state(istate)
      D=>extract_scalar_field(state(istate), "LayerThickness")
      U=>extract_vector_field(state(istate), "LocalVelocity")
      have_source = .false.
      if (has_vector_field(state(istate), "LocalVelocitySource")) then
        have_source = .true.
        source => extract_vector_field(state(istate), "LocalVelocitySource")
      end if
      dim = U%dim

      call execute_timestep_setup(D,U,d_rhs,u_rhs,advecting_u, &
           old_u,old_d,delta_d,delta_u)

      call insert(state(istate),advecting_u,"NonlinearVelocity")

      call set(advecting_u,u)
      call set(old_u,u)
      call set(old_d,d)

      do nit = 1, nonlinear_iterations

         call set(u,old_u)
         call set(d,old_d)

         !velocity advection step
         if(.not.exclude_velocity_advection) then
            do d1 = 1, dim
               velocity_cpt = extract_scalar_field(u,d1)
               old_velocity_cpt = extract_scalar_field(old_u,d1)
               call insert(state(istate),velocity_cpt,"VelocityComponent")
               call insert(state(istate),old_velocity_cpt,"VelocityComponentOld")
               call solve_advection_diffusion_dg("VelocityComponent", state(istate))
            end do
         end if
         !pressure advection
         if(.not.exclude_pressure_advection) then
            call solve_field_equation_cg("LayerThickness", state, 1, dt, &
                 "NonlinearVelocity")
         end if

         if(hybridized) then
            call solve_hybridized_helmholtz(&
                 &state(istate),&
                 &U_out=U_rhs,D_out=D_rhs,&
                 &compute_cartesian=.true.,&
                 &check_continuity=.true.,output_dense=.false.)
            ewrite(1,*) 'jump in D', maxval(abs(D_rhs%val-D%val))
            ewrite(1,*) 'jump in U', maxval(abs(U_rhs%val-U%val))
            call set(D,D_rhs)
            call set(U,U_rhs)
         else
            !Wave equation step
            ! M\Delta u + \Delta t F(\theta\Delta u + u^n) + \Delta t C(\theta
            ! \Delta\eta + \eta^n) = 0
            ! M\Delta h - \Delta t HC^T(\theta\Delta u + u^n) = 0
            ! SO
            ! (M+\theta\Delta t F)\Delta u =
            !   -\theta\Delta t g C\Delta\eta + \Delta t(-Fu^n - gC\eta^n)
            ! SET r = \Delta t(-Fu^n - gC\eta^n)
            ! THEN SUBSTITUTION GIVES
            ! (M + \theta^2\Delta t^2 gH C^T(M+\theta\Delta t F)^{-1}C)\Delta\eta
            !  = \Delta t HC^T(u^n + \theta(M+\theta\Delta t F)^{-1}r)

            !Construct explicit parts of u rhs
            if (have_source) then
               call get_u_rhs(u_rhs,U,D,dt,g, &
                    coriolis_mat,div_mat,u_mass_mat, source)
            else
               call get_u_rhs(u_rhs,U,D,dt,g, &
                    coriolis_mat,div_mat)
            end if

            !Construct explicit parts of h rhs in wave equation
            call get_d_rhs(d_rhs,u_rhs,D,U,div_mat,big_mat,D0,dt,theta)

            if (has_scalar_field(state(istate), "LayerThicknessSource")) then
               d_src => extract_scalar_field(state(istate), "LayerThicknessSource")
               call allocate(md_src, d_src%mesh, "MassMatrixTimesLayerThicknessSource")
               call zero(md_src)
               call mult(md_src, h_mass_mat, d_src)
               call addto(d_rhs, md_src, scale=dt)
               call deallocate(md_src)
            endif

            !Solve wave equation for D update
            delta_d%option_path = d%option_path
            call petsc_solve(delta_d, wave_mat, d_rhs)

            !Add the new D contributions into the RHS for u
            call update_u_rhs(u_rhs,U,delta_D,div_mat,theta,dt,g)

            !Solve momentum equation for U update
            call zero(delta_u)
            call mult(delta_u, big_mat, u_rhs)

            !Check the equation was solved correctly
            if(have_option("/debug/check_solution")) then
               call check_solution(delta_u,delta_d,d,u,dt,theta,g,D0,u_mass_mat&
                    &,h_mass_mat, coriolis_mat,div_mat)
            end if
         end if

         call addto(u,delta_u)
         call addto(d,delta_d)

         call set(advecting_u,old_u)
         call scale(advecting_u,(1-itheta))
         call addto(advecting_u,u,scale=itheta)

      end do

      call deallocate(d_rhs)
      call deallocate(u_rhs)
      call deallocate(delta_d)
      call deallocate(delta_u)
      call deallocate(old_u)
      call deallocate(old_d)
      call deallocate(advecting_u)

    end subroutine execute_timestep

    subroutine execute_timestep_setup(D,U,d_rhs,u_rhs,advecting_u, &
         old_u,old_d,delta_d,delta_u)
      implicit none
      type(scalar_field), pointer :: D
      type(scalar_field), intent(inout) :: D_rhs, delta_d, old_d
      type(vector_field), intent(inout), pointer :: U
      type(vector_field), intent(inout) :: U_rhs, delta_u, advecting_u, old_u
      integer :: dim

      dim = U%dim

      !allocate RHS variables
      call allocate(d_rhs, D%mesh, "LayerThickness_RHS")
      call zero(d_rhs)
      call allocate(u_rhs, dim, U%mesh, "Velocity_RHS")
      call zero(u_rhs)
      !allocate update variables
      call allocate(delta_d, D%mesh, "LayerThickness_update")
      call zero(delta_d)
      call allocate(delta_u, dim, U%mesh, "Velocity_update")
      call zero(delta_u)
      !allocate advecting velocity
      call allocate(advecting_u,dim,U%mesh, "Advecting_velocity")
      call zero(advecting_u)
      !allocate previous timestep variables
      call allocate(old_d, D%mesh, "LayerThickness_old")
      call zero(old_d)
      call allocate(old_u, U%dim, U%mesh, "Velocity_old")
      call zero(old_u)

    end subroutine execute_timestep_setup

    subroutine allocate_and_insert_additional_fields(state, adjoint)
      !!< Allocate and insert fields not specified in schema
      type(state_type), intent(inout) :: state
      logical, intent(in), optional :: adjoint

      ! coriolis
      type(scalar_field) :: f
      ! velocity in local coordinates
      type(vector_field) :: U_local

      type(vector_field), pointer :: X, U
      character(len=PYTHON_FUNC_LEN) :: coriolis
      integer :: stat

      X=>extract_vector_field(state, "Coordinate")
      if (present_and_true(adjoint)) then
        U=>extract_vector_field(state, "AdjointVelocity")
      else
        U=>extract_vector_field(state, "Velocity")
      end if

      call allocate(f, X%mesh, "Coriolis")
      call get_option("/physical_parameters/coriolis", coriolis, stat)
      if(stat==0) then
         call set_from_python_function(f, coriolis, X, time=0.0)
      else
         call zero(f)
      end if
      call insert(state, f, "Coriolis")
      call deallocate(f)

      if (present_and_true(adjoint)) then
        call allocate(U_local, mesh_dim(U), U%mesh, "AdjointLocalVelocity")
        call zero(U_local)
        call insert(state, U_local, "AdjointLocalVelocity")
        call deallocate(U_local)
      else
        call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocity")
        call zero(U_local)
        call insert(state, U_local, "LocalVelocity")
        call deallocate(U_local)
      endif

      if (has_vector_field(state, "VelocitySource")) then
        call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocitySource")
        call zero(U_local)
        call insert(state, U_local, "LocalVelocitySource")
        call deallocate(U_local)
      end if

    end subroutine allocate_and_insert_additional_fields

    subroutine project_velocity(state)
      !!< Project the quadratic continuous VelocityInitialCondition into
      !!< the velocity space.
      !!<
      !!< This is a hack to attempt to get third order convergence for wave
      !!< problems.
      type(state_type), intent(inout) :: state

      type(vector_field), pointer :: X, U, U_in

      integer :: ele

      call print_state(state, 0)

      X=>extract_vector_field(state, "Coordinate")
      U=>extract_vector_field(state, "Velocity")
      U_in=>extract_vector_field(state, "VelocityInitialCondition")

      do ele=1, element_count(U)

         call project_velocity_ele(ele, X, U, U_in)

      end do
    end subroutine project_velocity

    subroutine project_velocity_ele(ele, X, U, U_in)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: X
      type(vector_field), intent(inout) :: U
      type(vector_field), intent(in) :: U_in

      real, dimension(ele_ngi(U,ele)) :: detwei
      real, dimension(ele_loc(U,ele), ele_loc(U,ele)) :: mass

      real, dimension(U%dim, ele_loc(U,ele)) :: rhs

      type(element_type), pointer :: U_shape
      integer :: d

      call transform_to_physical(X, ele, detwei=detwei)

      U_shape=>ele_shape(U,ele)

      mass=shape_shape(U_shape, U_shape, detwei)

      call invert(mass)

      rhs=shape_vector_rhs(U_shape, ele_val_at_quad(U_in,ele), detwei)

      do d=1,U%dim
         call set(U, d, ele_nodes(U,ele), matmul(mass,rhs(d,:)))
      end do

    end subroutine project_velocity_ele

    subroutine advance_current_time(current_time, dt)
      implicit none
      real, intent(inout) :: current_time
      real, intent(in) :: dt

      ! Adaptive timestepping could go here.

      current_time=current_time + dt
      call set_option("/timestepping/current_time", current_time)

    end subroutine advance_current_time

    subroutine output_state(state, adjoint)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state
      logical, intent(in), optional :: adjoint

      ! project the local velocity to cartesian coordinates
      call project_local_to_cartesian(state(1), adjoint)

      ! Now we're ready to call write_state

      call write_state(dump_no, state, adjoint)
    end subroutine output_state

    subroutine set_velocity_from_geostrophic_balance(state,div_mat,&
         coriolis_mat, inverse_coriolis_mat)
      implicit none
      type(state_type), intent(inout) :: state
      type(block_csr_matrix), intent(in) :: div_mat, coriolis_mat, &
           inverse_coriolis_mat
      !
      type(scalar_field), pointer :: D
      type(vector_field), pointer :: U,X,source
      type(vector_field) :: u_tmp
      logical :: have_source

      ewrite(1,*) '    subroutine set_velocity_from_geostrophic_balance()'

      !Pull the fields out of state
      D=>extract_scalar_field(state, "LayerThickness")
      U=>extract_vector_field(state, "LocalVelocity")
      X=>extract_vector_field(state, "Coordinate")
      have_source = .false.
      if (has_vector_field(state, "VelocitySource")) then
        have_source = .true.
        source => extract_vector_field(state, "LocalVelocitySource")
      end if

      ewrite(2,*) 'inside', sum(D%val)

      call allocate(u_tmp,mesh_dim(U),U%mesh,name="tempmem")

      call mult_T(u_tmp,div_mat,D)
      call scale(u_tmp, -g)

      call mult(U, inverse_coriolis_mat, u_tmp)

      if (have_source) then
        call get_u_rhs(u_tmp,U,D,dt,g, &
           coriolis_mat,div_mat,u_mass_mat,source)
     else
        call get_u_rhs(u_tmp,U,D,dt,g, &
           coriolis_mat,div_mat)
     end if

      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'
      ewrite(3,*) 'SW: infty norm of u_tmp(1)=', maxval(abs(u_tmp%val(1,:)))
      ewrite(3,*) 'SW: infty norm of u_tmp(2)=', maxval(abs(u_tmp%val(2,:)))
      ewrite(3,*) 'SW: TESTING BALANCED INITIAL CONDITION'

      call deallocate(u_tmp)

      ewrite(1,*) 'END subroutine set_velocity_from_geostrophic_balance()'

    end subroutine set_velocity_from_geostrophic_balance

    subroutine read_command_line()
      implicit none
      ! Read the input filename.

      character(len=1024) :: argument
      integer :: status, argn, level

      call set_global_debug_level(0)

      argn=1
      do

         call get_command_argument(argn, value=argument, status=status)
         argn=argn+1

         if (status/=0) then
            call usage
            stop
         end if

         if (argument=="-v") then
            call get_command_argument(argn, value=argument, status=status)
            argn=argn+1

            if (status/=0) then
               call usage
               stop
            end if

            read(argument, "(i1)", err=666) level
            call set_global_debug_level(level)

            ! Go back to pick up the command line.
            cycle
         end if

         exit
      end do

      call load_options(argument)
      if(.not. have_option("/simulation_name")) goto 666

      return

666   call usage
      stop

    end subroutine read_command_line

    subroutine usage
      implicit none

      write (0,*) "usage: shallow_water [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

    subroutine compute_matrix_transposes(matrices)
      type(state_type), intent(inout) :: matrices
      type(csr_matrix), pointer :: wave_mat
      type(csr_matrix) :: wave_mat_T
      type(block_csr_matrix), pointer :: inv_coriolis_mat
      type(block_csr_matrix) :: inv_coriolis_mat_T

      wave_mat => extract_csr_matrix(matrices, "WaveMatrix")
      wave_mat_T = transpose(wave_mat, symmetric_sparsity=.true.)
      call insert(matrices, wave_mat_T, "WaveMatrixTranspose")
      call deallocate(wave_mat_T)

      inv_coriolis_mat => extract_block_csr_matrix(matrices, "InverseCoriolisMatrix")
      inv_coriolis_mat_T = transpose(inv_coriolis_mat, symmetric_sparsity=.false.)
      call insert(matrices, inv_coriolis_mat_T, "InverseCoriolisMatrixTranspose")
      call deallocate(inv_coriolis_mat_T)
    end subroutine compute_matrix_transposes

    subroutine adjoint_register_initial_eta_condition(states)
      type(state_type), dimension(:), intent(in) :: states
#ifdef HAVE_ADJOINT
      ! Register the initial condition for eta_0.
      type(adj_block) :: I
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: eta0
      real :: start_time
      real :: dt
      integer :: nfunctionals, j
      character(OPTION_PATH_LEN) :: buf, functional_name, mesh_name
      type(adj_variable), dimension(:), allocatable :: vars
      integer :: nmeshes
      type(mesh_type), pointer :: mesh
      type(vector_field), pointer :: u
      type(scalar_field), pointer :: eta
      type(adj_vector) :: mesh_vec
      type(adj_storage_data) :: storage

      ierr = adj_create_block("LayerThicknessMassMatrix", block=I, context=c_loc(matrices))
      call adj_chkierr(ierr)

      ierr = adj_create_variable("Fluid::LayerThickness", timestep=0, iteration=0, auxiliary=.false., variable=eta0)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(variable=eta0, blocks=(/I/), targets=(/eta0/), equation=equation)
      call adj_chkierr(ierr)

      ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)

      ierr = adj_destroy_equation(equation)
      ierr = adj_destroy_block(I)

      call get_option("/timestepping/current_time", start_time)
      call get_option("/timestepping/timestep", dt)

      ! We also may as well set the option paths for each variable now, since we know them here
      ! (in fluidity this would be done as each equation is registered)
      u => extract_vector_field(states(1), "Velocity")
      eta => extract_scalar_field(states(1), "LayerThickness")
      ierr = adj_dict_set(adj_path_lookup, "Fluid::Velocity", trim(u%option_path))
      ierr = adj_dict_set(adj_path_lookup, "Fluid::LayerThickness", trim(eta%option_path))

      ierr = adj_dict_set(adj_solver_path_lookup, "Fluid::LocalVelocityDelta", trim(u%option_path))
      ierr = adj_dict_set(adj_solver_path_lookup, "Fluid::LocalVelocity", trim(u%option_path))
      ierr = adj_dict_set(adj_solver_path_lookup, "Fluid::LayerThicknessDelta", trim(eta%option_path))
      ierr = adj_dict_set(adj_solver_path_lookup, "Fluid::LayerThickness", trim(eta%option_path))

      ! We may as well set the times for this timestep now
      ierr = adj_timestep_set_times(adjointer, timestep=0, start=start_time-dt, end=start_time)
      ! And we also may as well set the functional dependencies now
      nfunctionals = option_count("/adjoint/functional")
      do j=0,nfunctionals-1
        call get_option("/adjoint/functional[" // int2str(j) // "]/name", functional_name)
        call adj_record_anything_necessary(adjointer, python_timestep=1, timestep_to_record=0, functional=trim(functional_name), states=states)
      end do
#endif
    end subroutine adjoint_register_initial_eta_condition

    subroutine adjoint_register_initial_u_condition(balanced)
      logical, intent(in) :: balanced
#ifdef HAVE_ADJOINT
      type(adj_block) :: I, L, gC, E, P
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: local_u0, cartesian_u0, eta0
      logical :: check_transposes

      ! balanced is whether we derived the u initial condition from eta0 and geostrophic balance.
      ! if balanced is false, we just read in an initial condition as usual

      check_transposes = have_option("/adjoint/debug/check_action_transposes")

      if (.not. balanced) then ! we just read in u from file
        ierr = adj_create_block("CartesianVelocityMassMatrix", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(I, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::Velocity", timestep=0, iteration=0, auxiliary=.false., variable=cartesian_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=cartesian_u0, blocks=(/I/), targets=(/cartesian_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
        call adj_chkierr(ierr)
        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(I)

        ierr = adj_create_block("LocalVelocityMassMatrix", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(I, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)
        ierr = adj_create_block("MassLocalProjection", block=P, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(P, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(P, coefficient=-1.0)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::LocalVelocity", timestep=0, iteration=0, auxiliary=.false., variable=local_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=local_u0, blocks=(/P, I/), targets=(/cartesian_u0, local_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
        call adj_chkierr(ierr)
        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(I)
        ierr = adj_destroy_block(P)
      else
        ! The equation we have to register is the derivation of u0 from geostrophic balance
        ierr = adj_create_block("Coriolis", block=L, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(L, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)
        ierr = adj_create_block("Grad", block=gC, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(gC, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(block=gC, coefficient=g)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::LocalVelocity", timestep=0, iteration=0, auxiliary=.false., variable=local_u0)
        call adj_chkierr(ierr)
        ierr = adj_create_variable("Fluid::LayerThickness", timestep=0, iteration=0, auxiliary=.false., variable=eta0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=local_u0, blocks=(/L, gC/), targets=(/local_u0, eta0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
        call adj_chkierr(ierr)
        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(L)
        ierr = adj_destroy_block(gC)

        ierr = adj_create_block("CartesianVelocityMassMatrix", block=I, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(I, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)

        ierr = adj_create_block("MassCartesianProjection", block=P, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(P, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)
        ierr = adj_block_set_coefficient(P, coefficient=-1.0)
        call adj_chkierr(ierr)

        ierr = adj_create_variable("Fluid::Velocity", timestep=0, iteration=0, auxiliary=.false., variable=cartesian_u0)
        call adj_chkierr(ierr)

        ierr = adj_create_equation(variable=cartesian_u0, blocks=(/P, I/), targets=(/local_u0, cartesian_u0/), equation=equation)
        call adj_chkierr(ierr)

        ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
        call adj_chkierr(ierr)
        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)

        ierr = adj_destroy_equation(equation)
        ierr = adj_destroy_block(P)
        ierr = adj_destroy_block(I)
      endif
#endif
    end subroutine adjoint_register_initial_u_condition

    subroutine adjoint_register_timestep(timestep, dt, states)
      integer, intent(in) :: timestep
      real, intent(in) :: dt
      type(state_type), dimension(:), intent(in) :: states
#ifdef HAVE_ADJOINT
      type(adj_block) :: Mu, minusMu, Meta, minusMeta, W, CTMC, CTML, MBCdelta, MBL, MBC, CP, CI, P
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: u, previous_cartesian_u, delta_u, eta, previous_eta, delta_eta, cartesian_u
      real :: start_time

      integer :: j, nfunctionals
      type(adj_variable), dimension(:), allocatable :: vars
      character(len=OPTION_PATH_LEN) :: buf, functional_name

      type(adj_storage_data) :: storage_u, storage_eta
      type(scalar_field), pointer :: eta_ptr
      type(vector_field), pointer :: u_ptr
      type(adj_vector) :: u_vec, eta_vec
      logical :: check_transposes

      check_transposes = have_option("/adjoint/debug/check_action_transposes")

      ! Set up adj_variables
      ierr = adj_create_variable("Fluid::LocalVelocity", timestep=timestep, iteration=0, auxiliary=.false., variable=u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::Velocity", timestep=timestep-1, iteration=0, auxiliary=.false., variable=previous_cartesian_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LocalVelocityDelta", timestep=timestep, iteration=0, auxiliary=.false., variable=delta_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::Velocity", timestep=timestep, iteration=0, auxiliary=.false., variable=cartesian_u)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LayerThickness", timestep=timestep, iteration=0, auxiliary=.false., variable=eta)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LayerThickness", timestep=timestep-1, iteration=0, auxiliary=.false., variable=previous_eta)
      call adj_chkierr(ierr)
      ierr = adj_create_variable("Fluid::LayerThicknessDelta", timestep=timestep, iteration=0, auxiliary=.false., variable=delta_eta)
      call adj_chkierr(ierr)

      ! Set up adj_blocks

      ! Blocks for delta eta equation
      ierr = adj_create_block("WaveMatrix", context=c_loc(matrices), block=W)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(W, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_create_block("DivBigMatGrad", context=c_loc(matrices), block=CTMC)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(CTMC, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=CTMC, coefficient=dt**2 * D0 * theta * g)
      call adj_chkierr(ierr)
      ierr = adj_create_block("DivMinusDivBigMatCoriolisProjection", context=c_loc(matrices), block=CTML)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(CTML, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=CTML, coefficient=-1.0)
      call adj_chkierr(ierr)

      ! Blocks for eta_n equation
      ierr = adj_create_block("LayerThicknessMassMatrix", context=c_loc(matrices), block=Meta)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(Meta, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_create_block("LayerThicknessMassMatrix", context=c_loc(matrices), block=minusMeta)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(minusMeta, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=minusMeta, coefficient=-1.0)
      call adj_chkierr(ierr)

      ! Blocks for delta u equation
      ierr = adj_create_block("MassBigMatGrad", context=c_loc(matrices), block=MBCdelta)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(MBCdelta, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=MBCdelta, coefficient=theta * dt * g)
      call adj_chkierr(ierr)
      ierr = adj_create_block("MassBigMatCoriolisProjection", context=c_loc(matrices), block=MBL)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(MBL, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=MBL, coefficient=dt)
      call adj_chkierr(ierr)
      ierr = adj_create_block("MassBigMatGrad", context=c_loc(matrices), block=MBC)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(MBC, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=MBC, coefficient=dt * g)
      call adj_chkierr(ierr)

      ! Blocks for u_n equation
      ierr = adj_create_block("LocalVelocityMassMatrix", context=c_loc(matrices), block=Mu)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(Mu, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_create_block("LocalVelocityMassMatrix", context=c_loc(matrices), block=minusMu)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(minusMu, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(block=minusMu, coefficient=-1.0)
      call adj_chkierr(ierr)
      ierr = adj_create_block("MassLocalProjection", block=P, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(P, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(P, coefficient=-1.0)
      call adj_chkierr(ierr)

      ! Blocks for embedded manifold business
      ierr = adj_create_block("MassCartesianProjection", context=c_loc(matrices), block=CP)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(CP, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)
      ierr = adj_block_set_coefficient(CP, coefficient=-1.0)
      call adj_chkierr(ierr)
      ierr = adj_create_block("CartesianVelocityMassMatrix", context=c_loc(matrices), block=CI)
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(CI, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)

      ! Ahah! Now we can register our lovely equations. 
      ierr = adj_create_equation(delta_eta, blocks=(/CTMC, CTML, W/), &
                                          & targets=(/previous_eta, previous_cartesian_u, delta_eta/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(eta, blocks=(/minusMeta, minusMeta, Meta/), &
                                    & targets=(/previous_eta, delta_eta, eta/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(delta_u, blocks=(/MBC, MBL, MBCdelta, Mu/), &
                                        & targets=(/previous_eta, previous_cartesian_u, delta_eta, delta_u/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(u, blocks=(/P, minusMu, Mu/), &
                                    & targets=(/previous_cartesian_u, delta_u, u/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(cartesian_u, blocks=(/CP, CI/), &
                                    & targets=(/u, cartesian_u/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_callback(equation, c_funloc(shallow_water_forward_source))
      call adj_chkierr(ierr)
      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)
      ierr = adj_destroy_equation(equation)
      call adj_chkierr(ierr)

      ! And now we gots to destroy some blocks
      ierr = adj_destroy_block(Mu)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(minusMu)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(Meta)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(minusMeta)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(W)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(CTMC)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(CTML)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(MBCdelta)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(MBL)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(MBC)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(CP)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(CI)
      call adj_chkierr(ierr)
      ierr = adj_destroy_block(P)
      call adj_chkierr(ierr)

      ! Set the times and functional dependencies for this timestep
      call get_option("/timestepping/current_time", start_time)
      ierr = adj_timestep_set_times(adjointer, timestep=timestep-1, start=start_time, end=start_time+dt)
      nfunctionals = option_count("/adjoint/functional")
      do j=0,nfunctionals-1
        call get_option("/adjoint/functional[" // int2str(j) // "]/functional_dependencies/algorithm", buf)
        call get_option("/adjoint/functional[" // int2str(j) // "]/name", functional_name)
        call adj_variables_from_python(adjointer, buf, start_time, start_time+dt, timestep, vars)
        ierr = adj_timestep_set_functional_dependencies(adjointer, timestep=timestep-1, functional=trim(functional_name), &
                                                      & dependencies=vars)
        call adj_chkierr(ierr)
        deallocate(vars)
        ! We also need to check if these variables will be used
        call adj_record_anything_necessary(adjointer, python_timestep=timestep, timestep_to_record=timestep, functional=trim(functional_name), states=states)
      end do

      if (have_option("/adjoint/debug/replay_forward_run")) then
        u_ptr => extract_vector_field(states(1), "Velocity")
        u_vec = field_to_adj_vector(u_ptr)
        ierr = adj_storage_memory_copy(u_vec, storage_u)
        call adj_chkierr(ierr)
        ierr = adj_storage_set_overwrite(storage_u, .true.)
        call adj_chkierr(ierr)
        ierr = adj_record_variable(adjointer, cartesian_u, storage_u)
        call adj_chkierr(ierr)
        call femtools_vec_destroy_proc(u_vec)

        eta_ptr => extract_scalar_field(states(1), "LayerThickness")
        eta_vec = field_to_adj_vector(eta_ptr)
        ierr = adj_storage_memory_copy(eta_vec, storage_eta)
        call adj_chkierr(ierr)
        ierr = adj_storage_set_overwrite(storage_eta, .true.)
        call adj_chkierr(ierr)
        ierr = adj_record_variable(adjointer, eta, storage_eta)
        call adj_chkierr(ierr)
        call femtools_vec_destroy_proc(eta_vec)
      end if

      ! And that's it!
#endif
    end subroutine adjoint_register_timestep

  end program shallow_water
