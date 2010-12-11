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
      use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
    use adjoint_functionals
    implicit none

#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

    type(state_type), dimension(:), pointer :: state
    integer :: timestep
    integer :: ierr
    integer :: dim

    !! Matrices
    type(csr_matrix) :: mass_matrix, diffusion_matrix
    character(len = OPTION_PATH_LEN) :: simulation_name

    integer :: stat

    logical :: stationary_problem, adjoint
    real :: change, tolerance

    ! will be number of material_phases x number of timesteps
    type(state_type), dimension(:,:), pointer :: forward_state, adjoint_state

    integer :: no_timesteps

    type(scalar_field), pointer :: u, rhs
    integer :: dump_no=0


    interface
      subroutine set_global_debug_level(n)
        integer, intent(in) :: n
      end subroutine set_global_debug_level

      subroutine mpi_init(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_init

      subroutine mpi_finalize(ierr)
        integer, intent(out) :: ierr
      end subroutine mpi_finalize

      subroutine python_init
      end subroutine python_init

      subroutine petscinitialize(s, i)
        character(len=*), intent(in) :: s
        integer, intent(out) :: i
      end subroutine petscinitialize
    end interface
    
#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init()
    call read_command_line()

    call forward_options_dictionary
    call populate_state(state)
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
    call setup_matrices(state(1), mass_matrix, diffusion_matrix)
    call insert(state, mass_matrix, "MassMatrix")
    call insert(state, diffusion_matrix, "DiffusionMatrix")
    call deallocate(mass_matrix)
    call deallocate(diffusion_matrix)

    ! Allocate the space to store the forward system through time.
    adjoint = have_option("/adjoint")
    stationary_problem = have_option("/timestepping/steady_state")

    if (adjoint) then
      if (stationary_problem) then
        allocate(forward_state(1,1))
      else
        no_timesteps = get_number_of_timesteps()
        allocate(forward_state(1,1:no_timesteps+1))
      end if
    end if

    ! Always output the initial conditions.
    call output_state(state)
    if (stationary_problem) then
      call get_option("/timestepping/steady_state/tolerance", tolerance)
    end if

    if ((.not. stationary_problem) .and. adjoint) then
      ! We need a small amount of hackery. The initial Rhs should
      ! be the initial condition for u.
      u => extract_scalar_field(state(1), "Velocity")
      rhs => extract_scalar_field(state(1), "VelocityRhs", stat=stat)
      if (stat == 0) then
        call set(rhs, u)
      end if
      call insert(state(1), u, "IteratedVelocity1")
      call copy_forward_state(state, forward_state(:,1))
    end if

    timestep=0
    timestep_loop: do 
      timestep=timestep+1
      ewrite(1,*) "Starting timestep, current_time: ", current_time
      ! this may already have been done in populate_state, but now
      ! we evaluate at the correct "shifted" time level:
      call set_boundary_conditions_values(state, shift_time=.true.)
      
      ! evaluate prescribed fields at time = current_time+dt
      call set_prescribed_field_values(state, exclude_interpolated=.true., exclude_nonreprescribed=.true., time=current_time+dt)

      call execute_timestep(state(1), dt, change=change)

      call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
      call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)

      call advance_current_time(current_time, dt)
      call insert_time_in_state(state)

      if (simulation_completed(current_time, timestep)) exit timestep_loop
      if (stationary_problem) then
        if (change < tolerance) then
          exit timestep_loop
        end if
      end if

      if ((.not. stationary_problem) .and. adjoint) then
        ! Why doesn't bounds checking in the compiler do this for me?
        assert(timestep+1 <= size(forward_state, 2))
        call copy_forward_state(state, forward_state(:,timestep+1))
      end if

      if (do_write_state(current_time, timestep)) then
        call output_state(state)
      end if

      call write_diagnostics(state, current_time, dt, timestep)
      call delete_iterated_fields(state)

    end do timestep_loop

    ! One last dump
    call output_state(state)
    call write_diagnostics(state, current_time, dt, timestep)

    if (adjoint) then
      call copy_forward_state(state, forward_state(:,size(forward_state, 2)))
    end if

    call deallocate(state)
    call deallocate_transform_cache
    call deallocate_reserve_state    
    call close_diagnostic_files
    call uninitialise_diagnostics

    if (.not. adjoint) then
      call print_references(0)
    end if

    ! Now compute the adjoint
    if (adjoint) then
      ewrite(1,*) "Entering adjoint computation"

      call clear_options
      call read_command_line
      if (.not. have_option("/io/dump_period_in_timesteps/constant")) then
        call set_option('/io/dump_period_in_timesteps/constant', 1, stat=stat)
      end if
      call mangle_forward_options_paths(forward_state)
      call adjoint_options_dictionary
      call populate_state(state)

      ! Check the diagnostic field dependencies for circular dependencies
      call check_diagnostic_dependencies(state)
      call get_option('/simulation_name', simulation_name)
      call initialise_diagnostics(trim(simulation_name),state)

      allocate(adjoint_state(size(state), size(forward_state)))
      adjoint_state(:, size(forward_state)) = state
      deallocate(state)

      dump_no = dump_no - 1

      call compute_adjoint(forward_state, adjoint_state)

      call deallocate(adjoint_state)
      deallocate(adjoint_state)
      call deallocate(forward_state)
      deallocate(forward_state)

      call deallocate_transform_cache
      call deallocate_reserve_state    
      call close_diagnostic_files

      call print_references(0)
    else
      ewrite(1,*) "No adjoint specified, not entering adjoint computation"
    end if

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

    subroutine execute_timestep(state, dt, change)
      implicit none
      type(state_type), intent(inout) :: state
      real, intent(in) :: dt
      real, intent(out), optional :: change

      type(scalar_field) :: iterated_velocity
      type(scalar_field) :: old_iterated_velocity, iterated_velocity_difference, stored_iterated_velocity
      type(scalar_field) :: rhs, stored_iterated_srhs

      type(csr_matrix) :: advection_matrix
      type(csr_matrix) :: lhs_matrix

      type(scalar_field), pointer :: u, srhs
      type(vector_field), pointer :: x
      type(csr_matrix), pointer :: mass_matrix, diffusion_matrix

      integer :: nonlinear_iterations, nit
      real :: theta
      real :: nonlin_change
      logical :: adjoint, save_rhs
      integer :: stat

      mass_matrix => extract_csr_matrix(state, "MassMatrix")
      diffusion_matrix => extract_csr_matrix(state, "DiffusionMatrix")

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
      save_rhs = have_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/scalar_field::Rhs")

      call allocate(iterated_velocity, u%mesh, "IteratedVelocity")
      call set(iterated_velocity, u)
      call allocate(old_iterated_velocity, u%mesh, "OldIteratedVelocity")
      call allocate(iterated_velocity_difference, u%mesh, "IteratedVelocityDifference")

      call allocate(lhs_matrix, mass_matrix%sparsity, name="LeftHandSide")
      call allocate(advection_matrix, mass_matrix%sparsity, name="AdvectionMatrix")
      call allocate(rhs, u%mesh, "RightHandSide")
      if (save_rhs) then
       srhs => extract_scalar_field(state, "VelocityRhs")
      end if

      nonlinear_loop: do nit=1,nonlinear_iterations
        call assemble_advection_matrix(advection_matrix, x, u, iterated_velocity)

        call assemble_left_hand_side(lhs_matrix, mass_matrix, advection_matrix, diffusion_matrix, dt, theta, u)
        call assemble_right_hand_side(rhs, mass_matrix, advection_matrix, diffusion_matrix, dt, theta, u, state)
        if (save_rhs) then
          call assemble_small_right_hand_side(srhs, mass_matrix, u)
        end if
        
        call mangle_dirichlet_rows(lhs_matrix, u, keep_diag=.true., rhs=rhs)
        call apply_dirichlet_conditions(lhs_matrix, rhs, u)
        if (save_rhs) then
          call apply_dirichlet_conditions(lhs_matrix, srhs, u)  ! Put the values of the BCs on the right-hand side
        end if

        call set(old_iterated_velocity, iterated_velocity)
        call petsc_solve(iterated_velocity, lhs_matrix, rhs, &
               option_path="/material_phase::Fluid/scalar_field::Velocity/")

        call set(iterated_velocity_difference, iterated_velocity)
        call addto(iterated_velocity_difference, old_iterated_velocity, scale=-1.0)
        nonlin_change = norm2(iterated_velocity_difference, x)
        ewrite(1,*) "L2 norm of change in velocity: ", nonlin_change

        if (adjoint) then
          ! If we are adjointing, we need to store all the iterated velocities
          ! as well to propagate backwards in time
          if (nit /= nonlinear_iterations) then
            call allocate(stored_iterated_velocity, iterated_velocity%mesh, "IteratedVelocity" // int2str(nit))
            call set(stored_iterated_velocity, iterated_velocity)
            stored_iterated_velocity%option_path = u%option_path
            call insert(state, stored_iterated_velocity, trim(stored_iterated_velocity%name))
            call deallocate(stored_iterated_velocity)
          end if

          ! We might be also interested in the rhs of the forward state: For example in order to be able to compute <rhs, lambda>, which should give J(u)
          if (save_rhs) then
            call allocate(stored_iterated_srhs, srhs%mesh, "IteratedRhs" // int2str(nit))
            call set(stored_iterated_srhs, srhs)
            stored_iterated_srhs%option_path = srhs%option_path
            call insert(state, stored_iterated_srhs, trim(stored_iterated_srhs%name))
            call deallocate(stored_iterated_srhs)
          end if
        end if

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

    end subroutine execute_timestep

    subroutine assemble_advection_matrix(advection_matrix, x, u_left, u_right)
      type(csr_matrix), intent(inout) :: advection_matrix
      type(vector_field), intent(in) :: x
      type(scalar_field), intent(in), target :: u_left, u_right

      type(scalar_field) :: nu
      type(mesh_type), pointer :: mesh
      integer :: ele
      real :: itheta

      call get_option(trim(u_left%option_path) // "/prognostic/temporal_discretisation/relaxation", itheta, default=0.5)
      mesh => u_left%mesh
      call allocate(nu, mesh, "NonlinearVelocity")
      call set(nu, u_left)
      call scale(nu, (1.0 - itheta))
      call addto(nu, u_right, scale=itheta)
      nu%option_path = u_left%option_path

      call zero(advection_matrix)
      if (.not. have_option(trim(nu%option_path) // "/prognostic/remove_advection_term")) then
        do ele=1,ele_count(nu)
          call assemble_advection_matrix_ele(advection_matrix, x, nu, ele)
        end do
      end if

      call deallocate(nu)
    end subroutine assemble_advection_matrix

    subroutine assemble_advection_matrix_ele(advection_matrix, x, nu, ele)
      type(csr_matrix), intent(inout) :: advection_matrix
      type(vector_field), intent(in) :: x
      type(scalar_field), intent(in) :: nu
      integer, intent(in) :: ele

      real, dimension(ele_loc(nu, ele), ele_loc(nu, ele)) :: little_advection_matrix
      real, dimension(ele_ngi(nu, ele)) :: detwei
      real, dimension(ele_loc(nu, ele), ele_ngi(nu, ele), x%dim) :: du_t

      call transform_to_physical(x, ele, ele_shape(nu, ele), detwei=detwei, dshape=du_t)
      little_advection_matrix = shape_vector_dot_dshape(ele_shape(nu, ele), ele_val_at_quad(nu, ele), du_t, detwei)
      call addto(advection_matrix, ele_nodes(nu, ele), ele_nodes(nu, ele), little_advection_matrix)
    end subroutine assemble_advection_matrix_ele

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
      call add_neumann_boundary_conditions(rhs, x, u)

      src => extract_scalar_field(state, "VelocitySource", stat=stat)
      if (stat == 0) then
        call mult(tmp, mass_matrix, src)
        call addto(rhs, tmp)
      end if

      call deallocate(tmp)

    end subroutine assemble_right_hand_side

    
    ! The right hand side without the dependency of the old velocity field.
    subroutine assemble_small_right_hand_side(rhs, mass_matrix, u) 
      type(scalar_field), intent(inout) :: rhs
      type(csr_matrix), intent(in) :: mass_matrix
      type(scalar_field), intent(in), pointer :: u

      type(scalar_field) :: tmp
      type(scalar_field), pointer :: src
      integer :: stat

      type(vector_field), pointer :: x
      type(mesh_type), pointer :: mesh

      call zero(rhs)

      mesh => u%mesh
      call allocate(tmp, mesh, "TemporaryField")

      x => extract_vector_field(state, "Coordinate")
      call add_neumann_boundary_conditions(rhs, x, u)

      src => extract_scalar_field(state, "VelocitySource", stat=stat)
      if (stat == 0) then
        call mult(tmp, mass_matrix, src)
        call addto(rhs, tmp)
      end if

      call deallocate(tmp)

    end subroutine assemble_small_right_hand_side
    
    
    subroutine add_neumann_boundary_conditions(rhs, x, u)
      type(scalar_field), intent(inout) :: rhs
      type(scalar_field), intent(in) :: u
      type(vector_field), intent(in) :: x

      type(scalar_field) :: bc_value
      integer, dimension(surface_element_count(u)) :: bc_type_list
      integer :: sele
      real :: viscosity

      call get_option(trim(u%option_path) // "/prognostic/viscosity", viscosity)
      call get_entire_boundary_condition(u, (/"neumann"/), bc_value, bc_type_list)
      do sele=1,surface_element_count(u)
        if (bc_type_list(sele) == 0) then
          cycle
        end if

        call add_neumann_boundary_conditions_ele(rhs, x, u, sele, bc_value, viscosity)
      end do
      call deallocate(bc_value)
    end subroutine add_neumann_boundary_conditions

    subroutine add_neumann_boundary_conditions_ele(rhs, x, u, sele, bc_value, viscosity)
      type(scalar_field), intent(inout) :: rhs
      type(vector_field), intent(in) :: x
      type(scalar_field), intent(in) :: bc_value, u
      integer, intent(in) :: sele
      real, intent(in) :: viscosity

      real, dimension(face_ngi(u, sele)) :: detwei_face
      real, dimension(face_loc(u, sele), face_loc(u, sele)) :: face_mat

      call transform_facet_to_physical(x, sele, detwei_face)
      face_mat = shape_shape(face_shape(u, sele), face_shape(u, sele), detwei_face)
      call addto(rhs, face_global_nodes(u, sele), viscosity * matmul(face_mat, ele_val(bc_value, sele)))

    end subroutine add_neumann_boundary_conditions_ele

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
      integer :: increment

      if (present_and_true(adjoint)) then
        increment = -1
      else
        increment = 1
      end if

      call write_state(dump_no, state, increment)
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

    subroutine setup_matrices(state, mass_matrix, diffusion_matrix)
      type(state_type), intent(inout) :: state
      type(csr_matrix), intent(out) :: mass_matrix, diffusion_matrix

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

    subroutine compute_adjoint(forward_state, adjoint_state)
      ! material_phases x time
      type(state_type), intent(inout), dimension(:,:) :: forward_state
      type(state_type), intent(inout), dimension(:,:) :: adjoint_state


      assert(size(forward_state, 1) == 1) ! one phase
      assert(size(adjoint_state, 1) == 1) ! one phase

      if (have_option("/timestepping/steady_state")) then
        call steady_state_adjoint(forward_state, adjoint_state)
      else
        call time_dependent_adjoint(forward_state, adjoint_state)
      end if

    end subroutine compute_adjoint

    function compute_GT_ISP(sparsity, u, x, Au) result(GT)
      type(csr_sparsity), intent(in) :: sparsity
      type(scalar_field), intent(in), target :: u
      type(scalar_field), intent(in) :: Au
      type(vector_field), intent(in) :: x

      type(csr_matrix) :: GT
      type(csr_matrix) :: advection_matrix
      
      type(scalar_field) :: node_colour
      integer :: no_colours
      integer :: colour

      type(scalar_field) :: perturbed_u
      type(scalar_field) :: perturbed_Au
      real :: h

      integer :: node
      integer, dimension(:), pointer :: connected_nodes

      type(mesh_type), pointer :: mesh
      mesh => u%mesh

      call colour_graph(sparsity, mesh, node_colour, no_colours)

      call allocate(GT, sparsity, name="ForwardOperatorDerivative")
      call allocate(advection_matrix, sparsity, name="PerturbedAdvectionMatrix")
      call zero(GT)

      call allocate(perturbed_u, mesh, "PerturbedVelocity")
      perturbed_u%option_path = u%option_path
      call allocate(perturbed_Au, mesh, "PerturbedAu")

      do colour=1,no_colours
        ! Could be smarter about how to choose h.
        ! However, since A(.) is linear in u, for this problem the step size
        ! does not matter as you compute the exact derivative.
        h = 1.0

        call set(perturbed_u, u)

        ! Select the nodes of the currently-considered colour and perturb them
        do node=1,node_count(u%mesh)
          if (node_val(node_colour, node) == float(colour)) then
            call addto(perturbed_u, node, h)
          end if
        end do

        call assemble_advection_matrix(advection_matrix, x, perturbed_u, perturbed_u)
        call mult(perturbed_Au, advection_matrix, u)

        ! Subtract off Au
        call addto(perturbed_Au, Au, scale=-1.0)

        ! And divide by h
        call scale(perturbed_Au, 1.0/h)

        ! So now perturbed_Au contains the information for all of the rows of GT associated
        ! with this node. We need to split it up into the different rows.

        do node=1,node_count(u%mesh)
          if (node_val(node_colour, node) == float(colour)) then
            connected_nodes => row_m_ptr(sparsity, node)
            call set(GT, node, connected_nodes, node_val(perturbed_Au, connected_nodes))
          end if
        end do
      end do

      call deallocate(advection_matrix)
      call deallocate(perturbed_Au)
      call deallocate(perturbed_u)
      call deallocate(node_colour)

    end function compute_GT_ISP

    subroutine colour_graph(sparsity, mesh, node_colour, no_colours)
      type(csr_sparsity), intent(in) :: sparsity
      type(mesh_type), intent(inout) :: mesh

      type(scalar_field), intent(out) :: node_colour
      integer, intent(out) :: no_colours

      integer :: i
      integer, dimension(:), pointer :: row

      ! This is just to shut up the compiler about not using sparsity
      row => row_m_ptr(sparsity, 1)

      ! Here we cheat. Colouring a 1D mesh is quite easy. 
      ! We don't actually need to look at the sparsity at all.
      ! But for more general problems, you need to colour based on the sparsity
      ! pattern, so I added it to the arguments.

      
      call allocate(node_colour, mesh, "NodeColouring")

      no_colours = 3

      do i=1,node_count(mesh)
        call set(node_colour, i, float(mod(i-1, 3) + 1))
      end do
    end subroutine colour_graph

    function get_number_of_timesteps() result(no_timesteps)
      integer :: no_timesteps
      real :: start_time, dt, finish_time
      integer :: final_timestep, stat

      if (have_option("/timestepping/steady_state")) then
        no_timesteps = 1
        return
      end if

      call get_option("/timestepping/current_time", start_time)
      call get_option("/timestepping/timestep", dt)
      call get_option("/timestepping/finish_time", finish_time)
      call get_option("/timestepping/final_timestep", final_timestep, stat=stat)

      no_timesteps = ceiling((finish_time - start_time) / dt)
      if (stat == 0) then
        no_timesteps = min(no_timesteps, final_timestep)
      end if
    end function get_number_of_timesteps

    subroutine forward_options_dictionary
      ! Change the options dictionary to make it suitable
      ! for the forward model.
      integer :: i, nfields, nfields_to_delete
      character(len=OPTION_PATH_LEN) :: path, field_name
      character(len=OPTION_PATH_LEN), dimension(:), allocatable :: fields_to_delete

      nfields = option_count("/material_phase[0]/scalar_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the adjoint.
      nfields_to_delete = 0

      do i=0,nfields-1
        path = "/material_phase[0]/scalar_field["//int2str(i)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_adjoint')) then
          ! We can't actually delete it now, because if we have scalar_field[0,1,2]
          ! and delete scalar_field[1],
          ! now the numbering will be scalar_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if
      end do

      do i=1,nfields_to_delete
        path = "/material_phase[0]/scalar_field::" // trim(fields_to_delete(i))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)

    end subroutine forward_options_dictionary

    subroutine adjoint_options_dictionary
      ! Change the options dictionary to make it suitable
      ! for the adjoint model

      integer :: i, nfields, nfields_to_delete, nfields_to_move
      character(len=OPTION_PATH_LEN) :: path, field_name, simulation_name
      character(len=OPTION_PATH_LEN), dimension(:), allocatable :: fields_to_delete, fields_to_move
      integer :: stat
      real :: finish_time, current_time, dt

      nfields = option_count("/material_phase[0]/scalar_field")
      allocate(fields_to_delete(nfields))

      ! We need to delete fields that are marked as only living in the forward model.
      nfields_to_delete = 0

      do i=0,nfields-1
        path = "/material_phase[0]/scalar_field["//int2str(i)//"]"
        call get_option(trim(path)//"/name", field_name)

        if (have_option(trim(complete_field_path(trim(path))) // '/adjoint_storage/exists_in_forward')) then
          ! We can't actually delete it now, because if we have scalar_field[0,1,2]
          ! and delete scalar_field[1],
          ! now the numbering will be scalar_field[0,1]
          ! and the loop will go all wrong.
          nfields_to_delete = nfields_to_delete + 1
          fields_to_delete(nfields_to_delete) = field_name
        end if
      end do

      do i=1,nfields_to_delete
        path = "/material_phase[0]/scalar_field::" // trim(fields_to_delete(i))
        call delete_option(trim(path))
      end do

      deallocate(fields_to_delete)

      nfields = option_count("/material_phase[0]/scalar_field")
      allocate(fields_to_move(nfields))
      nfields_to_move = 0

      ! Now we need to go through and change the names of all prognostic fields
      ! to AdjointOldName.
      do i=0,nfields-1
        path = "/material_phase[0]/scalar_field["//int2str(i)//"]"
        call get_option(trim(path) // "/name", field_name)
        if (have_option(trim(path) // "/prognostic")) then
          ! We can't actually move it now, for exactly the same reason as above
          nfields_to_move = nfields_to_move + 1
          fields_to_move(nfields_to_move) = field_name
        end if
      end do

      do i=1,nfields_to_move
        path = "/material_phase[0]/scalar_field::"
        call move_option(trim(path) // trim(fields_to_move(i)), trim(path) // "Adjoint" // trim(fields_to_move(i)), stat=stat)
        assert(stat == SPUD_NO_ERROR)
        ! Delete any sources specified for the forward model -- we deal with these ourselves
        call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(i)) // "/prognostic/scalar_field::Source", stat=stat)
        ! And delete any initial conditions -- we will also specify these
        call delete_option(trim(path) // "Adjoint" // trim(fields_to_move(i)) // "/prognostic/initial_condition", stat=stat)
        call set_option(trim(path) // "Adjoint" // trim(fields_to_move(i)) // "/prognostic/initial_condition/constant", 0.0, stat=stat)
      end do

      deallocate(fields_to_move)

      ! And the simulation name
      call get_option("/simulation_name", simulation_name)
      call set_option("/simulation_name", trim(simulation_name) // "_adjoint", stat=stat)

      ! And the timestepping information
      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/finish_time", finish_time)
      call get_option("/timestepping/timestep", dt)

      call set_option("/timestepping/current_time", finish_time)
      call set_option("/timestepping/finish_time", current_time)
      call set_option("/timestepping/timestep", -dt)

    end subroutine adjoint_options_dictionary

    subroutine mangle_forward_options_paths(forward_state)
      ! We need to mangle the option_paths of the forward_state fields
      ! to change their name to AdjointWhatever
      type(state_type), dimension(:,:), intent(inout) :: forward_state

      type(scalar_field), pointer :: sfield
      type(vector_field), pointer :: vfield
      type(tensor_field), pointer :: tfield
      integer :: i, nfields, time, phase

      character(len=OPTION_PATH_LEN) :: new_path
      integer :: idx

      do phase=1,size(forward_state, 1)
        do time=1,size(forward_state, 2)
          nfields = scalar_field_count(forward_state(phase, time))
          do i=1,nfields
            sfield => extract_scalar_field(forward_state(phase, time), i)
            ! If we have a prognostic field
            if (have_option(trim(sfield%option_path) // '/prognostic')) then
              idx = index(sfield%option_path, "_field::") + len("_field::") - 1
              assert(idx /= 0)
              new_path(1:idx) = sfield%option_path(1:idx)
              new_path(idx+1:idx+8) = "Adjoint"
              new_path(idx+8:len_trim(sfield%option_path)+7) = sfield%option_path(idx+1:len_trim(sfield%option_path))
              new_path(len_trim(sfield%option_path)+8:) = " "
              sfield%option_path = new_path
            end if
          end do

          nfields = vector_field_count(forward_state(phase, time))
          do i=1,nfields
            vfield => extract_vector_field(forward_state(phase, time), i)
            ! If we have a prognostic field
            if (have_option(trim(vfield%option_path) // '/prognostic')) then
              idx = index(vfield%option_path, "_field::") + len("_field::") - 1
              assert(idx /= 0)
              new_path(1:idx) = vfield%option_path(1:idx)
              new_path(idx+1:idx+8) = "Adjoint"
              new_path(idx+8:len_trim(vfield%option_path)+7) = vfield%option_path(idx+1:len_trim(vfield%option_path))
              vfield%option_path = new_path
            end if
          end do

          nfields = tensor_field_count(forward_state(phase, time))
          do i=1,nfields
            tfield => extract_tensor_field(forward_state(phase, time), i)
            ! If we have a prognostic field
            if (have_option(trim(tfield%option_path) // '/prognostic')) then
              idx = index(tfield%option_path, "_field::") + len("_field::") - 1
              assert(idx /= 0)
              new_path(1:idx) = tfield%option_path(1:idx)
              new_path(idx+1:idx+8) = "Adjoint"
              new_path(idx+8:len_trim(tfield%option_path)+7) = tfield%option_path(idx+1:len_trim(tfield%option_path))
              tfield%option_path = new_path
            end if
          end do
        end do
      end do
    end subroutine mangle_forward_options_paths

    subroutine copy_forward_state(in_states, out_states)
      ! Copy the fields we need from in_state to out_state
      ! This routine is used to record information for the adjoint equation.
      ! Running backwards in time requires certain information about the
      ! forward state through time, and so we need to store that.

      type(state_type), dimension(:), intent(in) :: in_states
      type(state_type), dimension(size(in_states)), intent(out) :: out_states

      type(mesh_type), pointer :: mesh
      type(vector_field), pointer :: x

      integer :: i
      integer :: state, field
      type(scalar_field), pointer :: sfield
      type(scalar_field) :: cpy
      integer :: stat

      x => extract_vector_field(in_states, "Coordinate")
      if (have_option("/mesh_adaptivity")) then
        FLExit("Sorry, /adjoint and /mesh_adaptivity are currently incompatible")
      end if
      call insert(out_states, x, "Coordinate")

      do state=1,size(in_states)
        do field=1,scalar_field_count(in_states(state))
          sfield => extract_scalar_field(in_states(state), field)
          if (have_option(trim(complete_field_path(sfield%option_path, stat=stat)) // '/adjoint_storage/exists_in_both/record') .or. &
              have_option(trim(complete_field_path(sfield%option_path, stat=stat)) // '/adjoint_storage/exists_in_forward/record')) then
            call allocate(cpy, sfield%mesh, trim(sfield%name))
            call set(cpy, sfield)
            cpy%option_path = sfield%option_path
            call insert(out_states(state), cpy, trim(in_states(state)%scalar_names(field)))
            call deallocate(cpy)
          end if
        end do
      end do

      call populate_boundary_conditions(out_states)
      call set_boundary_conditions_values(out_states)

      do i=1,mesh_count(in_states(1))
        mesh => extract_mesh(in_states(1), i)
        call insert(out_states, mesh, trim(mesh%name))
      end do

      call insert(out_states, extract_csr_matrix(in_states(1), "MassMatrix"), "MassMatrix")
      call insert(out_states, extract_csr_matrix(in_states(1), "DiffusionMatrix"), "DiffusionMatrix")

      out_states%name = in_states%name
    end subroutine copy_forward_state

    subroutine mangle_dirichlet_rows(matrix, field, keep_diag, rhs)
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(in) :: field
      logical, intent(in) :: keep_diag
      type(scalar_field), intent(inout), optional :: rhs

      integer :: i, j, node
      character(len=FIELD_NAME_LEN) :: bctype
      integer, dimension(:), pointer :: node_list
      type(scalar_field), pointer :: bc_field

      real,  dimension(:), pointer :: row
      real, pointer :: diag_ptr
      real :: diag

      do i=1,get_boundary_condition_count(field)
        call get_boundary_condition(field, i, type=bctype, surface_node_list=node_list)

        if (bctype /= "dirichlet") cycle
        bc_field => extract_surface_field(field, i, "value")

        do j=1,size(node_list)
          node = node_list(j)
          diag_ptr => diag_val_ptr(matrix, node)
          diag = diag_ptr
          row  => row_val_ptr(matrix, node)
          row = 0
          if (keep_diag) then
            diag_ptr = 1.0
          end if

          if (present(rhs)) then
            call set(rhs, node, node_val(bc_field, j))
          end if
        end do
      end do
    end subroutine mangle_dirichlet_rows

    subroutine set_inactive_rows(matrix, field)
      ! Any nodes associated with Dirichlet BCs, we make them inactive
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(in) :: field

      integer :: i
      character(len=FIELD_NAME_LEN) :: bctype
      integer, dimension(:), pointer :: node_list

      do i=1,get_boundary_condition_count(field)
        call get_boundary_condition(field, i, type=bctype, surface_node_list=node_list)
        if (bctype /= "dirichlet") cycle
        call set_inactive(matrix, node_list)
      end do

    end subroutine set_inactive_rows

    subroutine compute_inactive_rows(x, A, rhs)                                                                                                                                                                              
      type(scalar_field), intent(inout):: x                                                                                                                                                                                 
      type(csr_matrix), intent(in):: A                                                                                                                                                                                      
      type(scalar_field), intent(inout):: rhs                                                                                                                                                                               

      logical, dimension(:), pointer:: inactive_mask                                                                                                                                                                        
      integer:: i                                                                                                                                                                                                           

      inactive_mask => get_inactive_mask(A)                                                                                                                                                                                 

      if (.not. associated(inactive_mask)) return                                                                                                                                                                           

      do i=1, size(A,1)                                                                                                                                                                                                     
        if (inactive_mask(i)) then                                                                                                                                                                                          

          ! we want [rhs_i - ((L+U)*x)_i]/A_ii, where L+U is off-diag part of A                                                                                                                                             
          ! by setting x_i to zero before this is the same as [rhs_i- (A*x)_i]/A_ii                                                                                                                                         
          call set(x, i, 0.0)                                                                                                                                                                                               
          call set(x, i, (node_val(rhs,i)- &                                                                                                                                                                                
          dot_product(row_val_ptr(A,i),node_val(x, row_m_ptr(A,i))) &                                                                                                                                                    
          &  )/val(A,i,i) )                                                                                                                                                                

        end if                                                                                                                                                                                                              
      end do                                                                                                                                                                                                                
    end subroutine compute_inactive_rows

    subroutine differentiate_V(x, sparsity, u_left, u_right, arg, contraction, adjoint, output)
      ! The core routine for the time-dependent adjoint equation.
      ! So, here's the arguments:
      type(vector_field), intent(in) :: x

      ! The sparsity of the advection operator
      type(csr_sparsity), intent(in) :: sparsity

      ! u_left and u_right are the arguments to the nonlinear velocity at which
      ! the advection matrix is assembled.
      ! nu = (1-itheta) * u_left + itheta * u_right
      type(scalar_field), intent(in), target :: u_left, u_right

      ! arg is an integer saying what we are differentiating with respect to.
      ! arg == -1 means differentiating with respect to the u_left argument.
      ! arg == +1 means differentiating with respect to the u_right argument.
      ! arg ==  0 means that u_left and u_right are actually the same, and we're
      !           taking the derivative with respect to that.
      integer, intent(in) :: arg

      ! So, the derivative of V is a rank-3 tensor. But we can't store them,
      ! so we immediately contract it with something else to give a rank-2
      ! tensor. What is it we're contracting it with?
      type(scalar_field), intent(in) :: contraction

      ! If you write down the form of G, we only ever multiply blocks of G
      ! by known vectors -- adjoint values previously computed. So this
      ! argument is the adjoint value to multiply it by.
      type(scalar_field), intent(in) :: adjoint

      ! The output, which is
      ! (diff(V[(1-itheta) * u_left + itheta * u_right], whatever) * contraction) * adjoint
      !                                                  ^ depends on arg
      !                                                            ^ rank-3 tensor contraction
      !                                                                           ^ plain old matrix multiplication
      type(scalar_field), intent(inout) :: output

      type(scalar_field) :: node_colour
      type(mesh_type), pointer :: mesh
      integer :: no_colours
      integer :: colour
      integer :: node
      integer, dimension(:), pointer :: connected_nodes

      type(csr_matrix) :: V ! the advection matrix, possibly perturbed, and possibly not
      type(scalar_field) :: perturbed ! some perturbed argument to assemble_advection_matrix
      type(scalar_field) :: Vc ! the unperturbed advection matrix, multiplied by contraction
      type(scalar_field) :: pVc ! the perturbed advection matrix, multiplied by contraction
      real, parameter :: h=1.0 ! the perturbation size; since our differentiation is exact regardless of h, we choose it to be 1.0

      call zero(output)
      mesh => u_left%mesh

      call allocate(V, sparsity, name="AdvectionMatrix")
      ! Firstly, assemble V unperturbed to get Vu
      call assemble_advection_matrix(V, x, u_left, u_right)

      ! We zero the rows, but not the columns.
      ! I don't know if this is the right thing to do, or 
      ! indeed what is the right thing to do.
      ! This needs a lot more thought.
      ! However, let's blunder on for now.

      call mangle_dirichlet_rows(V, u_left, keep_diag=.true.)
      call allocate(Vc, mesh, "UnperturbedVTimesContraction")
      call allocate(pVc, mesh, "PerturbedVTimesContraction")
      call mult(Vc, V, contraction)

      call colour_graph(sparsity, mesh, node_colour, no_colours)
      call allocate(perturbed, mesh, "PerturbedArgument")
      call zero(perturbed)

      do colour=1,no_colours

        if (arg == -1) then
          ! We are perturbing the left-hand argument to V
          call set(perturbed, u_left)
        else if (arg == 0) then
          ! In this case, u_left IS u_right, so it's all the same 
          call set(perturbed, u_left)
        else if (arg == 1) then
          ! Here, we're perturbing u_right
          call set(perturbed, u_right)
        else
          FLAbort("invalid value for arg")
        end if

        ! Assemble the perturbed state
        do node=1,node_count(mesh)
          if (node_val(node_colour, node) == float(colour)) then
            call addto(perturbed, node, h)
          end if
        end do

        call zero(V)
        if (arg == -1) then
          ! We are perturbing the left-hand argument to V
          call assemble_advection_matrix(V, x, perturbed, u_right)
        else if (arg == 0) then
          ! In this case, u_left IS u_right, so it's all the same 
          call assemble_advection_matrix(V, x, perturbed, perturbed)
        else if (arg == 1) then
          ! Here, we're perturbing u_right
          call assemble_advection_matrix(V, x, u_left, perturbed)
        else
          FLAbort("invalid value for arg")
        end if

        ! OK. So we've assembled the perturbed V. We multiply that by contraction:
        call mangle_dirichlet_rows(V, u_left, keep_diag=.true.)
        call mult(pVc, V, contraction)

        ! Now subtract off the unperturbed V * contraction:
        call addto(pVc, Vc, scale=-1.0)

        ! And divide by h:
        call scale(pVc, 1.0/h)

        ! So now pVc contains all the rows of GT associated with this colour.
        ! So let's extract out each row in turn, and multiply it by adjoint.
        do node=1,node_count(mesh)
          if (node_val(node_colour, node) == float(colour)) then
            connected_nodes => row_m_ptr(sparsity, node)
            ! node_val(pVc, connected_nodes) IS the row of GT associated with this node.
            call set(output, node, dot_product(node_val(pVc, connected_nodes), node_val(adjoint, connected_nodes)))
          end if
        end do
      end do

      call deallocate(V)
      call deallocate(Vc)
      call deallocate(pVc)
      call deallocate(perturbed)
      call deallocate(node_colour)
    end subroutine differentiate_V

    subroutine steady_state_adjoint(forward_state, adjoint_state)
      ! material_phases x time
      type(state_type), intent(inout), dimension(:,:) :: forward_state, adjoint_state

      type(scalar_field), pointer :: adjoint
      type(scalar_field), pointer :: u
      type(vector_field), pointer :: x
      type(scalar_field) :: Au
      type(scalar_field), pointer :: rhs
      type(state_type), dimension(size(forward_state,1)) :: derivatives

      type(csr_matrix) :: A, AT, GT, lhs, advection_matrix
      type(csr_matrix), target :: diffusion_matrix
      type(csr_matrix) :: mass_matrix
      type(csr_sparsity), pointer :: sparsity

      real :: current_time, dt

      u => extract_scalar_field(forward_state(1,1), "Velocity")
      adjoint  => extract_scalar_field(adjoint_state(1,1), "AdjointVelocity")
      x => extract_vector_field(adjoint_state(1,1), "Coordinate")

      call setup_matrices(forward_state(1,1), mass_matrix, diffusion_matrix)
      sparsity => diffusion_matrix%sparsity

      ! Assemble A again, to compute A^T
      call allocate(A, sparsity, name="ForwardOperator")
      call set(A, diffusion_matrix)

      call allocate(advection_matrix, sparsity, name="AdvectionMatrix")
      call assemble_advection_matrix(advection_matrix, x, u, u)
      call addto(A, advection_matrix)

      call allocate(Au, u%mesh, "Au")
      ! We also want (the nonlinear part of A) * u
      ! to subtract off in the formula for G^T later
      call mult(Au, advection_matrix, u)
      call deallocate(advection_matrix)

      call mangle_dirichlet_rows(A, u, keep_diag=.true.)
      AT = transpose(A, symmetric_sparsity=.true.)
      call deallocate(A)

      ! OK! Now to compute G^T.
      if (.not. have_option(trim(u%option_path) // "/prognostic/remove_advection_term")) then
        ewrite(1,*) "Equation is nonlinear, so computing derivatives of nonlinear terms"
        GT = compute_GT_ISP(sparsity, u, x, Au)
      else ! equation is linear, so G^T is zero
        ewrite(1,*) "Equation is linear, so not computing derivatives of nonlinear terms"
        call allocate(GT, sparsity, name="ForwardOperatorDerivative")
        call zero(GT)
      end if
      call deallocate(Au)

      call allocate(lhs, sparsity, name="AdjointOperator")
      call set(lhs, AT)
      call deallocate(AT)

      call addto(lhs, GT)
      call deallocate(GT)

      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/timestep", dt)
      call functional_derivative(forward_state, current_time, dt, 1, derivatives)
      rhs => extract_scalar_field(derivatives(1), "VelocityDerivative")

      call deallocate(mass_matrix)
      call deallocate(diffusion_matrix)

      call zero(adjoint)
      call set_inactive_rows(lhs, u)
      call petsc_solve(adjoint, lhs, rhs)
      call compute_inactive_rows(adjoint, lhs, rhs)
      call deallocate(lhs)
      call deallocate(derivatives)

      call calculate_diagnostic_variables(adjoint_state(:, 1), exclude_nonrecalculated = .true.)
      call calculate_diagnostic_variables_new(adjoint_state(:, 1), exclude_nonrecalculated = .true.)
      call write_diagnostics(adjoint_state(:, 1), current_time, dt, 1)

      call output_state(adjoint_state(:, 1), adjoint=.true.)
    end subroutine steady_state_adjoint

    function iteration_count(state) result(count)
      type(state_type), intent(in) :: state
      integer :: count

      integer :: i
      type(scalar_field), pointer :: iterated_velocity
      integer :: stat

      count = 0
      do
        i = count + 1
        iterated_velocity => extract_scalar_field(state, "IteratedVelocity" // int2str(i), stat=stat)
        if (stat /= 0) exit
        count = i
      end do
    end function iteration_count

    subroutine time_dependent_adjoint(forward_state, adjoint_state)
      ! material_phases x time
      type(state_type), intent(inout), dimension(:,:) :: forward_state
      type(state_type), intent(inout), dimension(:,:) :: adjoint_state
      integer :: timesteps, timestep
      real :: current_time, dt, finish_time
      real :: J
      logical :: have_functional
      integer :: iteration
      integer :: count

      timesteps = size(forward_state, 2)

      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/timestep", dt)

      timestep = timesteps
      have_functional = have_option("/adjoint/functional")

      if (have_functional) then
        J = functional(forward_state, current_time, dt, timestep)
        write(0,'(a, f0.5, a, f0.5)') "J(t=", current_time, ") == ", J
      end if

      call compute_final_adjoint_state(forward_state, adjoint_state, timestep)

      call calculate_diagnostic_variables(adjoint_state(:, timestep), exclude_nonrecalculated = .true.)
      call calculate_diagnostic_variables_new(adjoint_state(:, timestep), exclude_nonrecalculated = .true.)
      call write_diagnostics(adjoint_state(:, timestep), current_time, dt, timestep)
      call output_state(adjoint_state(:, timestep), adjoint=.true.)

      ! Rewind any nonlinear iterations, na ja?
      ! These neither cause evaluations of the functional, nor the functional derivative,
      ! nor diagnostics, nor output dumps, nor lines in the stat file, nor do they
      ! rewind time. But we still have to do them! 
      count = iteration_count(forward_state(1, timestep))
      do iteration=count-1,1,-1
        call rewind_nonlinear_iteration(forward_state, adjoint_state, timestep, iteration)
      end do

      call advance_current_time(current_time, dt)
      timestep = timestep - 1
      call setup_adjoint_state(forward_state(:, timestep), adjoint_state(:, timestep))

      do while (timestep > 1)
        if (have_functional) then
          J = functional(forward_state, current_time, dt, timestep)
          write(0,'(a, f0.5, a, f0.5)') "J(t=", current_time, ") == ", J
        end if

        FLAbort("Not implemented yet")
        ! Now compute the AdjointVelocity at the timestep itself
        !call compute_adjoint_timestep(forward_state, adjoint_state, timestep)

        call calculate_diagnostic_variables(adjoint_state(:, timestep), exclude_nonrecalculated = .true.)
        call calculate_diagnostic_variables_new(adjoint_state(:, timestep), exclude_nonrecalculated = .true.)
        call write_diagnostics(adjoint_state(:, timestep), current_time, dt, timestep)
        call output_state(adjoint_state(:, timestep), adjoint=.true.)

        ! Rewind any nonlinear iterations, na ja?
        ! These neither cause evaluations of the functional, nor the functional derivative,
        ! nor diagnostics, nor output dumps, nor lines in the stat file, nor do they
        ! rewind time. But we still have to do them! 
        count = iteration_count(forward_state(1, timestep))
        do iteration=count-1,1,-1
          call rewind_nonlinear_iteration(forward_state, adjoint_state, timestep, iteration)
        end do

        call advance_current_time(current_time, dt)
        timestep = timestep - 1
        call setup_adjoint_state(forward_state(:, timestep), adjoint_state(:, timestep))
      end do

      if (have_functional) then
        J = functional(forward_state, current_time, dt, timestep)
        write(0,'(a, f0.5, a, f0.5)') "J(t=", current_time, ") == ", J
      end if

      call compute_initial_adjoint_state(forward_state, adjoint_state)

      call calculate_diagnostic_variables(adjoint_state(:, timestep), exclude_nonrecalculated = .true.)
      call calculate_diagnostic_variables_new(adjoint_state(:, timestep), exclude_nonrecalculated = .true.)
      call write_diagnostics(adjoint_state(:, timestep), current_time, dt, timestep)
      call output_state(adjoint_state(:, timestep), adjoint=.true.)
      call get_option("/timestepping/finish_time", finish_time)
      assert(current_time == finish_time)
    end subroutine time_dependent_adjoint

    subroutine compute_final_adjoint_state(forward_state, adjoint_state, timestep)
      type(state_type), intent(inout), dimension(:,:) :: forward_state, adjoint_state
      integer, intent(in) :: timestep

      type(state_type), dimension(size(forward_state,1)) :: derivatives
      type(scalar_field), pointer :: rhs, adjoint, u, u_left, u_right
      type(csr_matrix) :: L, LT
      real :: current_time, dt, theta
      type(csr_matrix), pointer :: mass_matrix, diffusion_matrix
      type(csr_matrix) :: advection_matrix
      type(vector_field), pointer :: x
      integer :: count

      count = iteration_count(forward_state(1, timestep))
      u => extract_scalar_field(forward_state(1, timestep), "Velocity")
      adjoint => extract_scalar_field(adjoint_state(1, timestep), "AdjointVelocity")
      x => extract_vector_field(forward_state(1, timestep), "Coordinate")
      call zero(adjoint)

      call get_option(trim(u%option_path) // "/prognostic/temporal_discretisation/theta", theta, default=0.5)
      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/timestep", dt)

      mass_matrix => extract_csr_matrix(forward_state(1, timestep), "MassMatrix")
      diffusion_matrix => extract_csr_matrix(forward_state(1, timestep), "DiffusionMatrix")

      call allocate(advection_matrix, diffusion_matrix%sparsity, name="AdvectionMatrix")
      u_left => extract_scalar_field(forward_state(1, timestep-1), "Velocity")
      u_right => extract_scalar_field(forward_state(1, timestep-1), "IteratedVelocity" // int2str(count))

      call assemble_advection_matrix(advection_matrix, x, u_left, u_right)

      call allocate(L, diffusion_matrix%sparsity, name="LeftHandSide")

      call assemble_left_hand_side(L, mass_matrix, advection_matrix, diffusion_matrix, abs(dt), theta, u)
      call deallocate(advection_matrix)
      call mangle_dirichlet_rows(L, u, keep_diag=.true.)
      LT = transpose(L, symmetric_sparsity=.true.)
      call deallocate(L)

      call functional_derivative(forward_state, current_time, dt, timestep, derivatives)
      rhs => extract_scalar_field(derivatives(1), "VelocityDerivative")
      
      call set_inactive_rows(LT, u)
      call petsc_solve(adjoint, LT, rhs)
      call compute_inactive_rows(adjoint, LT, rhs)

      ! We'll also insert it as IteratedAdjointVelocityN
      call insert(adjoint_state(1, timestep), adjoint, "IteratedAdjointVelocity" // int2str(count))

      call deallocate(derivatives)
      call deallocate(LT)
    end subroutine compute_final_adjoint_state

    subroutine compute_initial_adjoint_state(forward_state, adjoint_state)
      type(state_type), intent(inout), dimension(:,:) :: forward_state, adjoint_state

      type(state_type), dimension(size(forward_state,1)) :: derivatives
      type(scalar_field), pointer :: rhs, adjoint, u, u_left, u_right
      real :: current_time, dt, theta
      type(csr_matrix), pointer :: mass_matrix, diffusion_matrix
      type(vector_field), pointer :: x
      type(scalar_field) :: output, contraction
      type(csr_sparsity), pointer :: sparsity
      type(csr_matrix) :: advection_matrix, T
      integer :: count, i

      u => extract_scalar_field(forward_state(1, 1), "Velocity")
      adjoint => extract_scalar_field(adjoint_state(1, 1), "AdjointVelocity")
      x => extract_vector_field(forward_state(1, 1), "Coordinate")
      call zero(adjoint)

      call get_option(trim(u%option_path) // "/prognostic/temporal_discretisation/theta", theta, default=0.5)
      call get_option("/timestepping/current_time", current_time)
      call get_option("/timestepping/timestep", dt)

      mass_matrix => extract_csr_matrix(forward_state(1, 1), "MassMatrix")
      diffusion_matrix => extract_csr_matrix(forward_state(1, 1), "DiffusionMatrix")
      sparsity => mass_matrix%sparsity

      ! ---------------------------------------------------------------------------------
      ! dJ/du
      ! ---------------------------------------------------------------------------------
      call functional_derivative(forward_state, current_time, dt, 1, derivatives)
      rhs => extract_scalar_field(derivatives(1), "VelocityDerivative")
      call addto(adjoint, rhs)
      call deallocate(derivatives)

      call allocate(output, u%mesh, "OutputVector")
      call allocate(contraction, u%mesh, "ContractionVector")

      ! ---------------------------------------------------------------------------------
      ! the G-blocks multiplying the adjoint RHS
      ! ---------------------------------------------------------------------------------
      ! The first G-block is a special case.
      u_left => u
      u_right => u

      call set(contraction, u)
      call scale(contraction, 1-theta)
      call addto(contraction, extract_scalar_field(forward_state(1, 2), "IteratedVelocity1"), scale=theta)

      call differentiate_V(x, sparsity, u_left, u_right, 0, contraction, &
           & extract_scalar_field(adjoint_state(1, 2), "IteratedAdjointVelocity1"), output)
      call addto(adjoint, output, scale=-1.0)

      ! Now loop over the others
      count = iteration_count(forward_state(1, 2))
      do i=2,count
        u_right => extract_scalar_field(forward_state(1, 2), "IteratedVelocity" // int2str(i-1))

        call set(contraction, u)
        call scale(contraction, 1-theta)
        call addto(contraction, extract_scalar_field(forward_state(1, 2), "IteratedVelocity" // int2str(i)), scale=theta)
        call differentiate_V(x, sparsity, u_left, u_right, -1, contraction, &
             & extract_scalar_field(adjoint_state(1, 2), "IteratedAdjointVelocity" // int2str(i)), output)
        call addto(adjoint, output, scale=-1.0)
      end do

      call deallocate(contraction)

      ! ---------------------------------------------------------------------------------
      ! the transpose of the forward timestepping operators
      ! ---------------------------------------------------------------------------------
      ! Again, the first T-block is a special case.

      ! Remember, 
      ! T = M/dt + (theta - 1) * V + (theta - 1) * D

      call allocate(advection_matrix, sparsity, name="AdvectionMatrix")
      u_left => u
      u_right => u
      call assemble_advection_matrix(advection_matrix, x, u_left, u_right)

      call allocate(T, sparsity, name="TimesteppingOperator")
      call set(T, mass_matrix)
      call scale(T, 1.0/abs(dt))

      call addto(T, diffusion_matrix, scale=theta-1.0)
      call addto(T, advection_matrix, scale=theta-1.0)
      call mangle_dirichlet_rows(T, u, keep_diag=.false.)

      call mult_T(output, T, extract_scalar_field(adjoint_state(1, 2), "IteratedAdjointVelocity1"))
      call addto(adjoint, output)

      do i=2,count
        u_right => extract_scalar_field(forward_state(1, 2), "IteratedVelocity" // int2str(i-1))
        call assemble_advection_matrix(advection_matrix, x, u_left, u_right)

        call set(T, mass_matrix)
        call scale(T, 1.0/abs(dt))
        call addto(T, diffusion_matrix, scale=theta-1.0)
        call addto(T, advection_matrix, scale=theta-1.0)
        call mangle_dirichlet_rows(T, u, keep_diag=.false.)

        call mult_T(output, T, extract_scalar_field(adjoint_state(1, 2), "IteratedAdjointVelocity" // int2str(i)))
        call addto(adjoint, output)
      end do
      
      call deallocate(advection_matrix)
      call deallocate(T)
      call deallocate(output)

      ! Note that we don't need to invert any operators, because the top-left block is
      ! just the identity matrix. That's why we've been accumulating the right-hand side
      ! into adjoint, not into an rhs variable, because the adjoint is the rhs.
      
    end subroutine compute_initial_adjoint_state

    subroutine rewind_nonlinear_iteration(forward_state, adjoint_state, timestep, iteration)
      type(state_type), intent(inout), dimension(:,:) :: forward_state, adjoint_state
      integer, intent(in) :: timestep, iteration

      type(scalar_field) :: rhs, contraction
      type(scalar_field), pointer :: adjoint, u, u_left, u_right
      real :: dt, theta
      type(csr_matrix), pointer :: mass_matrix, diffusion_matrix
      type(vector_field), pointer :: x
      type(csr_sparsity), pointer :: sparsity
      type(csr_matrix) :: advection_matrix, L, LT

      adjoint => extract_scalar_field(adjoint_state(1, timestep), "AdjointVelocity" // int2str(iteration))
      u => extract_scalar_field(forward_state(1, timestep-1), "Velocity")
      mass_matrix => extract_csr_matrix(forward_state(1, timestep), "MassMatrix")
      diffusion_matrix => extract_csr_matrix(forward_state(1, timestep), "DiffusionMatrix")
      x => extract_vector_field(adjoint_state(1, timestep), "Coordinate")
      sparsity => diffusion_matrix%sparsity

      call get_option(trim(u%option_path) // "/prognostic/temporal_discretisation/theta", theta, default=0.5)
      call get_option("/timestepping/timestep", dt)

      ! Assemble the right-hand side. This is the relevant G-block multiplied by the last value
      ! of the adjoint velocity that we've just computed.

      call allocate(rhs, adjoint%mesh, "NonlinearIterationRightHandSide")
      call zero(rhs)
      call allocate(contraction, adjoint%mesh, "ContractionVector")
      call set(contraction, u)
      call scale(contraction, (1.0 - theta))
      call addto(contraction, extract_scalar_field(forward_state(1,timestep), "IteratedVelocity" // int2str(iteration+1)), scale=theta)

      u_left => u
      u_right => extract_scalar_field(forward_state(1,timestep), "IteratedVelocity" // int2str(iteration))

      call differentiate_V(x, sparsity, u_left, u_right, 1, contraction, &
           & extract_scalar_field(adjoint_state(1, timestep), "IteratedAdjointVelocity" // int2str(iteration+1)), rhs)
      call scale(rhs, -1.0)

      call deallocate(contraction)

      ! OK. Now we have to assemble the operator. This is the relevant operator we solved
      ! for IteratedVelocity // int2str(iteration), transposed.
      call allocate(advection_matrix, sparsity, name="AdvectionMatrix")
      u_left => u
      if (iteration == 1) then
        u_right => u
      else
        u_right => extract_scalar_field(forward_state(1,timestep), "IteratedVelocity" // int2str(iteration-1))
      end if
      call assemble_advection_matrix(advection_matrix, x, u_left, u_right)

      call allocate(L, sparsity, name="LeftHandSide")
      call assemble_left_hand_side(L, mass_matrix, advection_matrix, diffusion_matrix, abs(dt), theta, u)
      call deallocate(advection_matrix)
      call mangle_dirichlet_rows(L, u, keep_diag=.true.)
      LT = transpose(L, symmetric_sparsity=.true.)
      call deallocate(L)

      call set_inactive_rows(LT, u)
      call petsc_solve(adjoint, LT, rhs)
      call compute_inactive_rows(adjoint, LT, rhs)

      call deallocate(LT)
      call deallocate(rhs)
    end subroutine rewind_nonlinear_iteration

    subroutine setup_adjoint_state(forward_state, adjoint_state)
      type(state_type), dimension(:), intent(in) :: forward_state
      type(state_type), dimension(:), intent(inout) :: adjoint_state

      type(mesh_type), pointer :: external_mesh
      type(vector_field), pointer :: x

      external_mesh => get_external_mesh(forward_state)
      x => extract_vector_field(forward_state, "Coordinate")
      call insert(adjoint_state, external_mesh, name=trim(external_mesh%name))
      call insert(adjoint_state, x, "Coordinate")
      call insert_derived_meshes(adjoint_state)
      call allocate_and_insert_fields(adjoint_state)
      call set_prescribed_field_values(adjoint_state)
    end subroutine setup_adjoint_state

    subroutine delete_iterated_fields(states)
      type(state_type), dimension(:), intent(inout) :: states
      type(scalar_field), pointer :: sfield
      integer :: field, state
      character(len=OPTION_PATH_LEN), dimension(:), allocatable :: fields_to_delete
      integer :: nfields_to_delete

      do state=1,size(states)
        allocate(fields_to_delete(scalar_field_count(states(state))))
        nfields_to_delete = 0

        do field=1,scalar_field_count(states(state))
          sfield => extract_scalar_field(states(state), field)
          if (index(trim(sfield%name), "Iterated") /= 0) then
            nfields_to_delete = nfields_to_delete + 1
            fields_to_delete(nfields_to_delete) = trim(sfield%name)
          end if
        end do

        do field=1,nfields_to_delete
          call remove_scalar_field(states(state), trim(fields_to_delete(field)))
        end do

        deallocate(fields_to_delete)
      end do
    end subroutine delete_iterated_fields
  end program burgers_equation
