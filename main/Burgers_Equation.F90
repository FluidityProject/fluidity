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
  program burgers_equation
    use spud
    use fields
    use state_module
    use FLDebug
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
    
#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call python_init()
    call read_command_line()

    call populate_state(state)
    call insert_time_in_state(state)

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state)

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

    ! Always output the initial conditions.
    call output_state(state)

    timestep=0
    timestep_loop: do 
       timestep=timestep+1
       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)
       
       ! evaluate prescribed fields at time = current_time+dt
       call set_prescribed_field_values(state, exclude_interpolated=.true., exclude_nonreprescribed=.true., time=current_time+dt)

       call execute_timestep(state(1), dt)

       call calculate_diagnostic_variables(state, exclude_nonrecalculated = .true.)
       call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)

       if (simulation_completed(current_time, timestep)) exit timestep_loop     

       call advance_current_time(current_time, dt)
       call insert_time_in_state(state)

       if (do_write_state(current_time, timestep)) then
          call output_state(state)
       end if

       call write_diagnostics(state, current_time, dt, timestep)
       
    end do timestep_loop

    ! One last dump
    call output_state(state)
    call write_diagnostics(state, current_time, dt, timestep)

    call deallocate(state)
    call deallocate_transform_cache
    call deallocate_reserve_state    
    call close_diagnostic_files

    call print_references(0)
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

    subroutine execute_timestep(state, dt)
      implicit none
      type(state_type), intent(inout) :: state
      real, intent(in) :: dt

      type(scalar_field) :: nonlinear_velocity
      type(scalar_field) :: iterated_velocity
      type(scalar_field) :: rhs

      type(csr_matrix) :: advection_matrix
      type(csr_matrix) :: lhs_matrix

      type(scalar_field), pointer :: u
      type(vector_field), pointer :: x
      type(csr_matrix), pointer :: mass_matrix, diffusion_matrix

      integer :: nonlinear_iterations, nit
      real :: itheta, theta

      mass_matrix => extract_csr_matrix(state, "MassMatrix")
      diffusion_matrix => extract_csr_matrix(state, "DiffusionMatrix")

      u => extract_scalar_field(state, "Velocity")
      x => extract_vector_field(state, "Coordinate")

      call get_option("/timestepping/nonlinear_iterations", nonlinear_iterations)
      call get_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/theta", theta)
      call get_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/temporal_discretisation/relaxation", itheta)

      call allocate(nonlinear_velocity, u%mesh, "NonlinearVelocity")
      call zero(nonlinear_velocity)
      call allocate(iterated_velocity, u%mesh, "IteratedVelocity")
      call set(iterated_velocity, u)

      call allocate(lhs_matrix, mass_matrix%sparsity, name="LeftHandSide")
      call allocate(advection_matrix, mass_matrix%sparsity, name="AdvectionMatrix")
      call allocate(rhs, u%mesh, "RightHandSide")

      do nit=1,nonlinear_iterations
        call set(nonlinear_velocity, u)
        call scale(nonlinear_velocity, itheta)
        call addto(nonlinear_velocity, iterated_velocity, scale=(1.0 - itheta))

        call assemble_advection_matrix(advection_matrix, x, nonlinear_velocity)

        call assemble_left_hand_side(lhs_matrix, mass_matrix, advection_matrix, diffusion_matrix, dt, theta)
        call assemble_right_hand_side(rhs, mass_matrix, advection_matrix, diffusion_matrix, dt, theta, u, state)

        call petsc_solve(iterated_velocity, lhs_matrix, rhs, &
               option_path="/material_phase::Fluid/scalar_field::Velocity/")
      end do

      call set(u, iterated_velocity)

      call deallocate(nonlinear_velocity)
      call deallocate(iterated_velocity)
      call deallocate(lhs_matrix)
      call deallocate(rhs)
      call deallocate(advection_matrix)

    end subroutine execute_timestep

    subroutine assemble_advection_matrix(advection_matrix, x, nu)
      type(csr_matrix), intent(inout) :: advection_matrix
      type(vector_field), intent(in) :: x
      type(scalar_field), intent(in) :: nu

      integer :: ele

      call zero(advection_matrix)

      do ele=1,ele_count(nu)
        call assemble_advection_matrix_ele(advection_matrix, x, nu, ele)
      end do
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

    subroutine assemble_left_hand_side(lhs_matrix, mass_matrix, advection_matrix, diffusion_matrix, dt, theta)
      type(csr_matrix), intent(inout) :: lhs_matrix
      type(csr_matrix), intent(in) :: mass_matrix, advection_matrix, diffusion_matrix
      real, intent(in) :: dt, theta

      call zero(lhs_matrix)

      call addto(lhs_matrix, mass_matrix, 1.0/dt)
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

      call mult(tmp, mass_matrix, u)
      call scale(tmp, 1.0/dt)
      call addto(rhs, tmp)

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

    subroutine add_neumann_boundary_conditions(rhs, x, u)
      type(scalar_field), intent(inout) :: rhs
      type(scalar_field), intent(in) :: u
      type(vector_field), intent(in) :: x

      type(scalar_field) :: bc_value
      integer, dimension(surface_element_count(u)) :: bc_type_list
      integer :: sele
      real :: viscosity

      call get_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/viscosity", viscosity)
      call get_entire_boundary_condition(u, (/"neumann"/), bc_value, bc_type_list)
      do sele=1,surface_element_count(u)
        if (bc_type_list(sele) == 0) then
          write(0,*) "Warning: surface element with no Neumann boundary condition" 
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

    subroutine output_state(state)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state

      integer, save :: dump_no=0

      call write_state(dump_no, state)
      dump_no = dump_no + 1

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

      call get_option("/material_phase::Fluid/scalar_field::Velocity/prognostic/viscosity", viscosity)

      x => extract_vector_field(state, "Coordinate")
      u => extract_scalar_field(state, "Velocity")
      u_sparsity => get_csr_sparsity_firstorder(state, u%mesh, u%mesh)

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

  end program burgers_equation
