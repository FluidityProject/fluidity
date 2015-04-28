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
  program burgers_equation
    use spud
    use elements
    use fields
    use field_options
    use state_module
    use fldebug
    use vtk_interfaces
    use signal_vars
    use populate_state_module
    use write_state_module
    use timeloop_utilities
    use sparsity_patterns_meshes
    use sparse_matrices_fields
    use solvers
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use global_parameters, only: option_path_len, current_time, dt
    use memory_diagnostics
    use reserve_state_module
    use boundary_conditions_from_options
    use boundary_conditions
    use colouring
      use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
    use adjoint_functionals
    use mangle_options_tree
    use mangle_dirichlet_rows_module
    use burgers_assembly, only: assemble_advection_matrix
    use burgers_adjoint_controls
    use adjoint_controls
#ifdef HAVE_ADJOINT
    use burgers_adjoint_callbacks
    use libadjoint
    use libadjoint_data_callbacks
    use adjoint_functional_evaluation
    use adjoint_python
    use adjoint_global_variables
    use adjoint_main_loop, only: compute_adjoint
    use forward_main_loop, only: register_functional_callbacks, calculate_functional_values, compute_forward
#include "libadjoint/adj_fortran.h"
#endif


#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  implicit none
#include "petsc_legacy.h"

    type(state_type), dimension(:), pointer :: state
    type(state_type), target :: matrices ! We collect all the cached
                                         ! matrices in this state so that
                                         ! the adjoint callbacks can use
                                         ! them

    integer :: timestep
    real :: theta
    integer :: ierr
    integer :: dim

    character(len = OPTION_PATH_LEN) :: simulation_name

    integer :: stat

    logical :: stationary_problem, adjoint
    real :: change, tolerance

    integer :: dump_no=0


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


#ifdef HAVE_ADJOINT
    ierr = adj_create_adjointer(adjointer)
    ! Register the data callbacks
    call adj_register_femtools_data_callbacks(adjointer)
    ! Register the operator callbacks
    call register_burgers_operator_callbacks(adjointer)
#endif  
    call initialise_walltime

#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init()
    call read_command_line()

    call mangle_options_tree_forward
    
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
    
    timestep=0
    call populate_state(state)
    call adjoint_load_controls(timestep, dt, state)
#ifdef HAVE_ADJOINT
    call adjoint_register_initial_condition(state)
#endif

    call insert_time_in_state(state)


    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state)

    if (.not. have_option("/io/dump_period_in_timesteps/constant")) then
      call set_option('/io/dump_period_in_timesteps/constant', 1, stat=stat)
    end if

    ! No support for multiphase or multimaterial at this stage.
    if (size(state)/=1) then
       FLExit("Multiple material_phases are not supported")
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)
    call get_option("/geometry/dimension", dim)

    if (dim /= 1) then
      FLExit("Sorry! I'm only one-dimensional.")
    end if
    
    ! Set up the matrices that don't change (everything but the nonlinear term)
    call setup_matrices(state(1), matrices)

    ! Always output the initial conditions.
    call output_state(state)
    stationary_problem = have_option("/timestepping/steady_state")
    if (stationary_problem) then
      call get_option("/timestepping/steady_state/tolerance", tolerance)
    end if
    
    ! Compute the diagnostics for the initial timestep
    call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
    call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
    call adjoint_write_controls(timestep, dt, state)

    timestep_loop: do 
      timestep=timestep+1
      ewrite(1,*) "Starting timestep, current_time: ", current_time
      ! this may already have been done in populate_state, but now
      ! we evaluate at the correct "shifted" time level:
      call set_boundary_conditions_values(state, shift_time=.true.)
      
      ! evaluate prescribed fields at time = current_time+theta*dt
      call get_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/theta", theta, default=0.5)
      call set_prescribed_field_values(state, exclude_interpolated=.true., exclude_nonreprescribed=.true., time=current_time+theta*dt)
      call adjoint_load_controls(timestep, dt, state)

      call execute_timestep(state, matrices, timestep, dt, change=change)

      call set_prescribed_field_values(state, exclude_interpolated=.true., exclude_nonreprescribed=.true., time=current_time+dt)
      call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
      call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
      call adjoint_write_controls(timestep, dt, state)

      call advance_current_time(current_time, dt)
      call insert_time_in_state(state)

      if (simulation_completed(current_time, timestep)) exit timestep_loop
      if (stationary_problem) then
        if (change < tolerance) then
          exit timestep_loop
        end if
      end if
      
      if (do_write_state(current_time, timestep)) then
        call output_state(state)
      end if

#ifdef HAVE_ADJOINT
      call calculate_functional_values(timestep-1)
#endif
      call write_diagnostics(state, current_time, dt, timestep)

    end do timestep_loop

    ! One last dump
    call output_state(state)
#ifdef HAVE_ADJOINT
    call calculate_functional_values(timestep-1)
#endif
    call write_diagnostics(state, current_time, dt, timestep)

    call deallocate(state)
    call deallocate_transform_cache
    call deallocate_reserve_state    
    call close_diagnostic_files
    call uninitialise_diagnostics

    if (.not. adjoint) then
      call deallocate(matrices)
      call print_references(0)
    end if

    ! Now compute the adjoint
#ifdef HAVE_ADJOINT
    if (adjoint) then

      if (have_option("/adjoint/debug/replay_forward_run")) then
        ! Let's run the forward model through libadjoint, too, for the craic
        ewrite(1,*) "Entering forward computation through libadjoint"
        call clear_options
        call read_command_line
        call mangle_options_tree_forward
        call populate_state(state)
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
      call check_diagnostic_dependencies(state)
      call compute_matrix_transposes(matrices)

      dump_no = dump_no - 1

      call compute_adjoint(state, dump_no, burgers_adjoint_timestep_callback, c_loc(matrices))

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

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

  contains

    subroutine insert_time_in_state(state)
      type(state_type), dimension(:), intent(inout) :: state
      
      type(scalar_field) :: aux_sfield
      type(mesh_type), pointer :: x_mesh
      
      x_mesh => extract_mesh(state, "CoordinateMesh")
      call allocate(aux_sfield, x_mesh, "Time", field_type=FIELD_TYPE_CONSTANT)
      call get_option("/timestepping/current_time", current_time)
      call set(aux_sfield, current_time)
      aux_sfield%option_path = ""
      call insert(state, aux_sfield, trim(aux_sfield%name))
      call deallocate(aux_sfield)
    
    end subroutine insert_time_in_state

    subroutine execute_timestep(states, matrices, timestep, dt, change)
      implicit none
      type(state_type), dimension(:), intent(inout) :: states
      type(state_type), intent(inout) :: matrices
      integer, intent(in) :: timestep
      real, intent(in) :: dt
      real, intent(out), optional :: change

      type(scalar_field) :: iterated_velocity
      type(scalar_field) :: old_iterated_velocity, iterated_velocity_difference
      type(scalar_field) :: rhs

      type(csr_matrix) :: advection_matrix
      type(csr_matrix) :: lhs_matrix

      type(scalar_field), pointer :: u
      type(vector_field), pointer :: x
      type(csr_matrix), pointer :: mass_matrix, diffusion_matrix

      integer :: nonlinear_iterations, nit
      real :: theta
      real :: nonlin_change
      logical :: adjoint
      integer :: stat

      integer :: dummy_timestep
      type(mesh_type), pointer :: mesh_ptr

      mesh_ptr => extract_mesh(states(1), "VelocityMesh")
      dummy_timestep = timestep

      mass_matrix => extract_csr_matrix(matrices, "MassMatrix")
      diffusion_matrix => extract_csr_matrix(matrices, "DiffusionMatrix")

      u => extract_scalar_field(state, "Velocity")
      x => extract_vector_field(state, "Coordinate")

      call get_option("/timestepping/nonlinear_iterations", nonlinear_iterations, default=2)
      if (have_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/remove_advection_term")) then
        ewrite(1,*) "No advection term, setting nonlinear iterations to 1"
        call set_option("/timestepping/nonlinear_iterations", 1, stat=stat)
        nonlinear_iterations = 1
      end if

      call get_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/theta", theta, default=0.5)

      adjoint = have_option("/adjoint")

      call allocate(iterated_velocity, u%mesh, "Velocity")
      call set(iterated_velocity, u)
      call allocate(old_iterated_velocity, u%mesh, "OldIteratedVelocity")
      call allocate(iterated_velocity_difference, u%mesh, "IteratedVelocityDifference")

      call allocate(lhs_matrix, mass_matrix%sparsity, name="LeftHandSide")
      call allocate(advection_matrix, mass_matrix%sparsity, name="AdvectionMatrix")
      call allocate(rhs, u%mesh, "RightHandSide")

      nonlinear_loop: do nit=1,nonlinear_iterations
        call assemble_advection_matrix(advection_matrix, x, u, iterated_velocity)

        call assemble_left_hand_side(lhs_matrix, mass_matrix, advection_matrix, diffusion_matrix, dt, theta, u)
        call assemble_right_hand_side(rhs, mass_matrix, advection_matrix, diffusion_matrix, dt, theta, u, state(1))
        
        call mangle_dirichlet_rows(lhs_matrix, u, keep_diag=.true., rhs=rhs)
        call set_inactive_rows(lhs_matrix, u)

        call set(old_iterated_velocity, iterated_velocity)

        call petsc_solve(iterated_velocity, lhs_matrix, rhs, &
               option_path="/material_phase::Fluid/scalar_field::Velocity/")
        call compute_inactive_rows(iterated_velocity, lhs_matrix, rhs)

        call set(iterated_velocity_difference, iterated_velocity)
        call addto(iterated_velocity_difference, old_iterated_velocity, scale=-1.0)
        nonlin_change = norm2(iterated_velocity_difference, x)
        ewrite(1,*) "L2 norm of change in velocity: ", nonlin_change

        ! Tell libadjoint about the equations we are solving
#ifdef HAVE_ADJOINT
        if (adjoint) then
          call adjoint_record_velocity(timestep, nit, states, field=iterated_velocity)
        end if
#endif

        if (sig_hup .or. sig_int) then
          ewrite(1,*) "Caught signal, exiting"
          exit nonlinear_loop
        end if
      end do nonlinear_loop

      ewrite(1,*) "Out of the nonlinear loop, number of nonlinear_iterations: ", nit - 1

      if (present(change)) then
        call set(iterated_velocity_difference, u)
        call addto(iterated_velocity_difference, iterated_velocity, scale=-1.0)
        change = norm2(iterated_velocity_difference, x)
        ewrite(1,*) "L2 norm of change over timestep: ", change
      end if

      call set(u, iterated_velocity)
      call deallocate(iterated_velocity)
      call deallocate(old_iterated_velocity)
      call deallocate(iterated_velocity_difference)
      call deallocate(lhs_matrix)
      call deallocate(rhs)
      call deallocate(advection_matrix)

#ifdef HAVE_ADJOINT
      call adjoint_register_timestep(timestep, dt, nit-1, states)
#endif

    end subroutine execute_timestep

    subroutine assemble_left_hand_side(lhs_matrix, mass_matrix, advection_matrix, diffusion_matrix, dt, theta, u)
      type(csr_matrix), intent(inout) :: lhs_matrix
      type(csr_matrix), intent(in) :: mass_matrix, advection_matrix, diffusion_matrix
      type(scalar_field), intent(in) :: u
      real, intent(in) :: dt, theta

      call zero(lhs_matrix)

      if (.not. have_option(trim(u%option_path) // '/prognostic/temporal_discretisation/remove_time_term')) then
        call addto(lhs_matrix, mass_matrix, 1.0/dt)
      end if
      call addto(lhs_matrix, advection_matrix, theta)
      call addto(lhs_matrix, diffusion_matrix, theta)

    end subroutine assemble_left_hand_side

    subroutine assemble_right_hand_side(rhs, mass_matrix, advection_matrix, diffusion_matrix, dt, theta, u, state)
      type(scalar_field), intent(inout) :: rhs
      type(csr_matrix), intent(in) :: mass_matrix, advection_matrix, diffusion_matrix
      real, intent(in) :: dt, theta
      type(scalar_field), intent(in), pointer :: u
      type(state_type), intent(in) :: state

      type(scalar_field) :: tmp
      type(scalar_field), pointer :: src
      integer :: stat

      type(vector_field), pointer :: x
      type(mesh_type), pointer :: mesh
      type(csr_matrix) :: mangled_mass_matrix

      call zero(rhs)

      mesh => u%mesh
      call allocate(tmp, mesh, "TemporaryField")

      if (.not. have_option(trim(u%option_path) // '/prognostic/temporal_discretisation/remove_time_term')) then
        call mult(tmp, mass_matrix, u)
        call scale(tmp, 1.0/dt)
        call addto(rhs, tmp)
      end if

      call mult(tmp, advection_matrix, u)
      call scale(tmp, theta - 1.0)
      call addto(rhs, tmp)

      call mult(tmp, diffusion_matrix, u)
      call scale(tmp, theta - 1.0)
      call addto(rhs, tmp)

      x => extract_vector_field(state, "Coordinate")

      src => extract_scalar_field(state, "VelocitySource", stat=stat)
      if (stat == 0) then
        call allocate(mangled_mass_matrix, mass_matrix%sparsity, name="MangledMassMatrix")
        call set(mangled_mass_matrix, mass_matrix)
        call mangle_dirichlet_rows(mangled_mass_matrix, u, keep_diag=.false.)
        call mult(tmp, mangled_mass_matrix, src)
        call deallocate(mangled_mass_matrix)
        call addto(rhs, tmp)
      end if

      call deallocate(tmp)

    end subroutine assemble_right_hand_side

    subroutine advance_current_time(current_time, dt)
      implicit none
      real, intent(inout) :: current_time, dt

      ! Adaptive timestepping could go here.

      current_time=current_time + dt
      call set_option("/timestepping/current_time", current_time)

    end subroutine advance_current_time

    subroutine output_state(state, adjoint)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state
      logical, intent(in), optional :: adjoint

      call write_state(dump_no, state, adjoint)
    end subroutine output_state

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

      write (0,*) "usage: burgers_equation [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage

    subroutine setup_matrices(state, matrices)
      type(state_type), intent(inout) :: state
      type(state_type), intent(inout) :: matrices
      type(csr_matrix) :: mass_matrix, diffusion_matrix

      real :: viscosity
      integer :: ele
      type(csr_sparsity), pointer :: u_sparsity
      type(scalar_field), pointer :: u
      type(vector_field), pointer :: x


      x => extract_vector_field(state, "Coordinate")
      u => extract_scalar_field(state, "Velocity")
      u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)
      call get_option(trim(u%option_path) // "/prognostic/viscosity", viscosity)

      call allocate(mass_matrix, u_sparsity, name="MassMatrix")
      call zero(mass_matrix)
      call allocate(diffusion_matrix, u_sparsity, name="DiffusionMatrix")
      call zero(diffusion_matrix)

      do ele=1,ele_count(u)
        call setup_matrices_ele(x, u, ele, mass_matrix, diffusion_matrix)
      end do

      call scale(diffusion_matrix, viscosity)
    
      call insert(matrices, mass_matrix, "MassMatrix")
      call insert(matrices, diffusion_matrix, "DiffusionMatrix")
      call insert(matrices, u%mesh, "VelocityMesh")
      call insert(matrices, x, "Coordinate")
      call deallocate(mass_matrix)
      call deallocate(diffusion_matrix)
    end subroutine setup_matrices

    subroutine setup_matrices_ele(x, u, ele, mass_matrix, diffusion_matrix)
      type(vector_field), intent(in) :: x
      type(scalar_field), intent(in) :: u
      integer, intent(in) :: ele

      type(csr_matrix), intent(inout) :: mass_matrix, diffusion_matrix
      real, dimension(ele_ngi(u, ele)) :: detwei
      real, dimension(ele_loc(u, ele), ele_ngi(u, ele), dim) :: dshape
      real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: little_mass_matrix, little_diffusion_matrix

      call transform_to_physical(x, ele, ele_shape(u, ele), detwei=detwei, dshape=dshape)

      little_mass_matrix = shape_shape(ele_shape(u, ele), ele_shape(u, ele), detwei)
      call addto(mass_matrix, ele_nodes(u, ele), ele_nodes(u, ele), little_mass_matrix)

      little_diffusion_matrix = dshape_dot_dshape(dshape, dshape, detwei)
      call addto(diffusion_matrix, ele_nodes(u, ele), ele_nodes(u, ele), little_diffusion_matrix)

    end subroutine setup_matrices_ele

    subroutine compute_matrix_transposes(matrices)
      type(state_type), intent(inout) :: matrices
      type(csr_matrix), pointer :: mass_matrix
      type(csr_matrix) :: mass_matrix_T
      type(csr_matrix), pointer :: diffusion_matrix
      type(csr_matrix) :: diffusion_matrix_T

      mass_matrix => extract_csr_matrix(matrices, "MassMatrix")
      mass_matrix_T = transpose(mass_matrix, symmetric_sparsity=.true.)
      call insert(matrices, mass_matrix_T, "MassMatrixTranspose")
      call deallocate(mass_matrix_T)

      diffusion_matrix => extract_csr_matrix(matrices, "DiffusionMatrix")
      diffusion_matrix_T = transpose(diffusion_matrix, symmetric_sparsity=.true.)
      call insert(matrices, diffusion_matrix_T, "DiffusionMatrixTranspose")
      call deallocate(diffusion_matrix_T)
    end subroutine compute_matrix_transposes

#ifdef HAVE_ADJOINT
    subroutine adjoint_register_initial_condition(states)
      type(state_type), dimension(:), intent(in) :: states
      ! Register the initial condition.
      type(adj_block) :: I
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: u0
      real :: start_time
      real :: dt
      integer :: nfunctionals, j
      character(OPTION_PATH_LEN) :: functional_name
      type(scalar_field), pointer :: u
      logical :: check_transposes

      if (.not. adjoint) return

      check_transposes = have_option("/adjoint/debug/check_action_transposes")

      ierr = adj_create_block("VelocityIdentity", block=I, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_block_set_test_hermitian(I, check_transposes, 100, 1.0d-10)
      call adj_chkierr(ierr)

      ierr = adj_create_variable("Fluid::Velocity", timestep=0, iteration=0, auxiliary=.false., variable=u0)
      call adj_chkierr(ierr)

      ierr = adj_create_equation(variable=u0, blocks=(/I/), targets=(/u0/), equation=equation)
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
      call adj_chkierr(ierr)
      ierr = adj_equation_set_rhs_callback(equation, c_funloc(burgers_equation_forward_source))
      call adj_chkierr(ierr)

      ierr = adj_register_equation(adjointer, equation)
      call adj_chkierr(ierr)

      ierr = adj_destroy_equation(equation)
      ierr = adj_destroy_block(I)

      call get_option("/timestepping/current_time", start_time)
      call get_option("/timestepping/timestep", dt)

      ! We also may as well set the option paths for each variable now, since we know them here
      ! (in fluidity this would be done as each equation is registered)
      u => extract_scalar_field(states(1), "Velocity")
      ierr = adj_dict_set(adj_path_lookup, "Fluid::Velocity", trim(u%option_path))

      ! We may as well set the times for this timestep now
      ierr = adj_timestep_set_times(adjointer, timestep=0, start=start_time-dt, end=start_time)
      ! And we also may as well set the functional dependencies now
      nfunctionals = option_count("/adjoint/functional")
      do j=0,nfunctionals-1
        call get_option("/adjoint/functional[" // int2str(j) // "]/name", functional_name)
        call adj_record_anything_necessary(adjointer, python_timestep=1, timestep_to_record=0, functional=trim(functional_name), states=states)
      end do
      call adjoint_record_velocity(0, 1, states)
    end subroutine adjoint_register_initial_condition
#endif

#ifdef HAVE_ADJOINT
    subroutine adjoint_register_timestep(timestep, dt, niterations, states)
      integer, intent(in) :: timestep, niterations
      real, intent(in) :: dt
      type(state_type), dimension(:), intent(in) :: states
      type(adj_block) :: timestepping_block, burgers_block
      type(adj_nonlinear_block) :: burgers_advection_block, timestepping_advection_block
      integer :: ierr
      type(adj_equation) :: equation
      type(adj_variable) :: u, previous_u, iter_u
      real :: start_time

      integer :: j, nfunctionals
      type(adj_variable), dimension(:), allocatable :: vars
      character(len=OPTION_PATH_LEN) :: buf, functional_name
      integer :: iteration, niterations_prev
      real :: theta
      logical :: check_transposes, check_derivatives

      check_transposes = have_option("/adjoint/debug/check_action_transposes")
      check_derivatives = have_option("/adjoint/debug/check_action_derivative")

      if (.not. adjoint) return

      call get_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/theta", theta, default=0.5)

      do iteration=1,niterations
        ! Set up adj_variables
        ierr = adj_create_variable("Fluid::Velocity", timestep=timestep, iteration=iteration-1, auxiliary=.false., variable=u)
        call adj_chkierr(ierr)

        ! We need to find out how many nonlinear iterations we did at the previous timestep, so that we can reference it
        ierr = adj_create_variable("Fluid::Velocity", timestep=timestep-1, iteration=0, auxiliary=.false., variable=previous_u)
        call adj_chkierr(ierr)
        ierr = adj_iteration_count(adjointer, previous_u, niterations_prev)
        call adj_chkierr(ierr)
        ierr = adj_create_variable("Fluid::Velocity", timestep=timestep-1, iteration=niterations_prev-1, auxiliary=.false., variable=previous_u)
        call adj_chkierr(ierr)

        if (iteration==1) then
          ! Set up adj_blocks
          ierr = adj_create_nonlinear_block("AdvectionOperator", (/ previous_u /), context=c_loc(matrices), coefficient=theta, nblock=burgers_advection_block)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_hermitian(burgers_advection_block, check_transposes, 100, 1.0d-10)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_derivative(burgers_advection_block, check_derivatives, 5)
          call adj_chkierr(ierr)

          ierr = adj_create_nonlinear_block("AdvectionOperator", (/ previous_u /), context=c_loc(matrices), coefficient=1.0-theta, nblock=timestepping_advection_block)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_hermitian(timestepping_advection_block, check_transposes, 100, 1.0d-10)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_derivative(timestepping_advection_block, check_derivatives, 5)
          call adj_chkierr(ierr)
        else 
          ierr = adj_create_variable("Fluid::Velocity", timestep=timestep, iteration=iteration-2, auxiliary=.false., variable=iter_u) 
          call adj_chkierr(ierr)

          ierr = adj_create_nonlinear_block("AdvectionOperator", (/ previous_u, iter_u /), context=c_loc(matrices), coefficient=theta, nblock=burgers_advection_block)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_hermitian(burgers_advection_block, check_transposes, 100, 1.0d-10)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_derivative(burgers_advection_block, check_derivatives, 5)
          call adj_chkierr(ierr)

          ierr = adj_create_nonlinear_block("AdvectionOperator", (/ previous_u, iter_u /), context=c_loc(matrices), coefficient=1.0-theta, nblock=timestepping_advection_block)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_hermitian(timestepping_advection_block, check_transposes, 100, 1.0d-10)
          call adj_chkierr(ierr)
          ierr = adj_nonlinear_block_set_test_derivative(timestepping_advection_block, check_derivatives, 5)
          call adj_chkierr(ierr)
        end if

        ierr = adj_create_block("TimesteppingOperator", timestepping_advection_block, context=c_loc(matrices), block=timestepping_block)
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(timestepping_block, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)
        ierr = adj_create_block("BurgersOperator", burgers_advection_block, context=c_loc(matrices), block=burgers_block)
        call adj_chkierr(ierr)
        ierr = adj_block_set_test_hermitian(burgers_block, check_transposes, 100, 1.0d-10)
        call adj_chkierr(ierr)

        ! Now we can register our lovely equation.
        ierr = adj_create_equation(u, blocks=(/timestepping_block, burgers_block/), &
                                            & targets=(/previous_u, u/), equation=equation)
        call adj_chkierr(ierr)
        ierr = adj_equation_set_rhs_dependencies(equation, context=c_loc(matrices))
        call adj_chkierr(ierr)
        ierr = adj_equation_set_rhs_callback(equation, c_funloc(burgers_equation_forward_source))
        call adj_chkierr(ierr)

        ierr = adj_register_equation(adjointer, equation)
        call adj_chkierr(ierr)
        ierr = adj_destroy_equation(equation)
        call adj_chkierr(ierr)

        ! And now we gots to destroy some blocks
        ierr = adj_destroy_block(timestepping_block)
        call adj_chkierr(ierr)
        ierr = adj_destroy_block(burgers_block)
        call adj_chkierr(ierr)
        ierr = adj_destroy_nonlinear_block(timestepping_advection_block)
        call adj_chkierr(ierr)
        ierr = adj_destroy_nonlinear_block(burgers_advection_block)
        call adj_chkierr(ierr)
      end do

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

    end subroutine adjoint_register_timestep
#endif

#ifdef HAVE_ADJOINT
    subroutine adjoint_record_velocity(timestep, iteration, states, field)
      integer, intent(in) :: timestep, iteration
      type(state_type), dimension(:), intent(in) :: states
      type(scalar_field), intent(in), optional, target :: field
      type(adj_storage_data) :: storage
      integer :: ierr
      type(adj_vector) :: u_vec
      type(scalar_field), pointer :: u, bc_u
      type(adj_variable) :: u_var
      integer :: bc
      character(len=FIELD_NAME_LEN) :: bc_name, bc_type
      integer, dimension(:), pointer :: surface_element_list

      ierr = adj_create_variable("Fluid::Velocity", timestep=timestep, iteration=iteration-1, auxiliary=.false., variable=u_var)
      call adj_chkierr(ierr)

      if (present(field)) then
        u => field
        ! We need to record any dirichlet boundary conditions on Velocity, and put them onto the IteratedVelocity too
        bc_u => extract_scalar_field(states(1), "Velocity")
        do bc=1,get_boundary_condition_count(bc_u)
          call get_boundary_condition(bc_u, bc, type=bc_type, surface_element_list=surface_element_list, name=bc_name)
          if (trim(bc_type) == "dirichlet") then
            call add_boundary_condition_surface_elements(u, bc_name, bc_type, surface_element_list)
          end if
        end do
      else
        u => extract_scalar_field(states(1), "Velocity")
      end if
      u_vec = field_to_adj_vector(u)

      ierr = adj_storage_memory_copy(u_vec, storage)
      call adj_chkierr(ierr)
      ierr = adj_storage_set_overwrite(storage, .true.) ! In the earlier cases, they may not have the BCs
      call adj_chkierr(ierr)

      ierr = adj_record_variable(adjointer, u_var, storage);
      call femtools_vec_destroy_proc(u_vec)
      if (ierr /= ADJ_WARN_ALREADY_RECORDED) then
        call adj_chkierr(ierr)
      end if
    end subroutine adjoint_record_velocity
#endif
  end program burgers_equation
