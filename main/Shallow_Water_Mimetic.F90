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
    use global_parameters, only: option_path_len, python_func_len, current_time, dt
    use memory_diagnostics
    use reserve_state_module
    use boundary_conditions_from_options
      use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
    use iso_c_binding
    use mangle_options_tree
    use manifold_tools
    use FEFields
    use field_copies_diagnostics
    implicit none
#ifdef HAVE_PETSC
#include "finclude/petsc.h"
#endif

    ! Interface blocks for the initialisation routines we need to call
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

    type(state_type), dimension(:), pointer :: states
    character(len = OPTION_PATH_LEN) :: simulation_name

    integer :: timestep, nonlinear_iterations
    integer :: ierr
    integer, save :: dump_no=0
    integer :: stat

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

    timestep=0

    call populate_state(states)

    ! No support for multiphase or multimaterial at this stage.
    if (size(states)/=1) then
       FLExit("Multiple material_phases are not supported")
    end if

    call insert_time_in_state(states)

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(states)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),states)
    
    call setup_fields()

    call calculate_diagnostic_variables(states)
    call calculate_diagnostic_variables_new(states)
    call get_option("/timestepping/timestep", dt)
    call write_diagnostics(states, current_time, dt, timestep)

    ! Always output the initial conditions.
    call output_state(states)

    timestep_loop: do
       timestep=timestep+1
       if (simulation_completed(current_time, timestep)) exit timestep_loop
       ewrite (1,*) "SW: start of timestep ", timestep, current_time
       call execute_timestep(states(1))

       call project_local_to_cartesian(states(1))
       call calculate_diagnostic_variables(states,&
            & exclude_nonrecalculated = .true.)
       call calculate_diagnostic_variables_new(states,&
            & exclude_nonrecalculated = .true.)

       call advance_current_time(current_time, dt)
       if (simulation_completed(current_time, timestep)) exit timestep_loop

       if (do_write_state(current_time, timestep)) then
          call output_state(states)
       end if

       call write_diagnostics(states,current_time, dt, timestep)

    end do timestep_loop

    ! One last dump
    call output_state(states)
    call write_diagnostics(states,current_time, dt, timestep)

    call deallocate(states)
    call deallocate_transform_cache
    call deallocate_reserve_state
    call close_diagnostic_files
    call uninitialise_diagnostics
    ! Clean up registered diagnostics
    call destroy_registered_diagnostics 

    call print_references(0)

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

#ifdef HAVE_MPI
    call mpi_finalize(ierr)
    assert(ierr == MPI_SUCCESS)
#endif
  contains
    subroutine execute_timestep(state)
      type(state_type), intent(inout) :: state
      !
      type(vector_field), pointer :: U
      type(scalar_field), pointer :: D
      type(vector_field) :: newU
      type(scalar_field) :: newD
      integer :: nonlinear_iterations, nits

      call get_option("/timestepping/nonlinear_iterations"&
           &,nonlinear_iterations)

      !Set up iterative solutions and provide initial guess
      !From previous timestep
      U => extract_vector_field(states, "LocalVelocity", stat)
      D => extract_scalar_field(states, "LayerThickness", stat)
      call allocate(newU,U%dim,U%mesh,"NewLocalVelocity", stat)
      call allocate(newD,D%mesh,"NewLayerThickness", stat)
      call set(newD,D)
      call set(newU,U)

      !Here check for wave equation solve, if it is switched off,
      !Just advance the D and PV fields.
      do nits = 1, nonlinear_iterations
         call solve_hybridised_timestep_residual(state,newU,newD)
      end do

      call set(D,newD)
      call set(U,newU)

      call deallocate(newU)
      call deallocate(newD)
    end subroutine execute_timestep

    subroutine setup_fields()
      type(vector_field), pointer :: v_field,U,X
      type(scalar_field), pointer :: s_field,D,f_ptr
      type(mesh_type), pointer :: v_mesh
      type(vector_field) :: U_local
      character(len=PYTHON_FUNC_LEN) :: coriolis
      type(scalar_field) :: f

      X=>extract_vector_field(states, "Coordinate")
      U=>extract_vector_field(states, "Velocity")

      !SET UP LOCAL VELOCITY
      !This needs an option to switch on as we don't always want to do it.
      !   !project velocity into div-conforming space
      call allocate(U_local, mesh_dim(U), U%mesh, "LocalVelocity")
      call zero(U_local)
      call insert(states, U_local, "LocalVelocity")
      call deallocate(U_local)
      v_field => extract_vector_field(states(1), "Velocity", stat)
      call project_cartesian_to_local(states(1), v_field)
      
      !SET UP CORIOLIS FORCE
      f_ptr => extract_scalar_field(states,"Coriolis",stat=stat)
      if(stat.ne.0) then
         call allocate(f, U%mesh, "Coriolis")
         call get_option("/physical_parameters/coriolis", coriolis, stat)
         if(stat==0) then
            call set_from_python_function(f, coriolis, X, time=0.0)
         else
            call zero(f)
         end if
         call insert(states, f, "Coriolis")
         call deallocate(f)
      end if
      v_field => extract_vector_field(states(1), "LocalVelocity", stat)
      call project_to_constrained_space(states(1),v_field)

    !VARIOUS BALANCED INITIAL OPTIONS
    ! Geostrophic balanced initial condition, if required
    if(have_option("/material_phase::Fluid/vector_field::Velocity/prognostic&
         &/initial_condition::WholeMesh/balanced")) then
       call set_velocity_from_geostrophic_balance_hybridized(states(1))
    end if
    !Set velocity from commuting projection?
    if(have_option("/material_phase::Fluid/vector_field::Velocity/&
         &prognostic/initial_condition::WholeMesh/&
         &commuting_projection")) then
       call set_velocity_commuting_projection(states(1))
    end if
    if(have_option("/material_phase::Fluid/vector_field::&
         &PrescribedVelocityFromCommutingProjection")) then
       call set_velocity_commuting_projection(states(1),"PrescribedVelocityFr&
            &omCommutingProjection")
    end if

    !Set velocity from spherical components
    if(have_option("/material_phase::Fluid/vector_field::Velocity/prognost&
         &ic/initial_condition::WholeMesh/from_sphere_pullback")) then
       call set_velocity_from_sphere_pullback(states(1))
    end if
    
    if(have_option("/material_phase::Fluid/scalar_field::LayerThickness/pr&
         &ognostic/initial_condition::ProjectionFromPython")) then
       call set_layerthickness_projection(states(1))
    end if
    if(have_option("/material_phase::Fluid/scalar_field::PrescribedLayerDe&
         &pthFromProjection")) then
       call set_layerthickness_projection(states(1),&
            &"PrescribedLayerDepthFromProjection")
    end if
    end subroutine setup_fields

    subroutine advance_current_time(current_time, dt)
      implicit none
      real, intent(inout) :: current_time
      real, intent(in) :: dt

      ! Adaptive timestepping could go here.

      current_time=current_time + dt
      call set_option("/timestepping/current_time", current_time)

    end subroutine advance_current_time

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

    subroutine output_state(state)
      implicit none
      type(state_type), dimension(:), intent(inout) :: state

      ! project the local velocity to cartesian coordinates
      call project_local_to_cartesian(state(1))
      ! Now we're ready to call write_state
      call write_state(dump_no, state)
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

      write (0,*) "usage: shallow_water [-v n] <options_file>"
      write (0,*) ""
      write (0,*) "-v n sets the verbosity of debugging"
    end subroutine usage


  end program shallow_water

! TODO LIST FOR BDFM1 SWE

! Matrix-ify Helmholtz solver - CODED
! Newton iteration for linear equations - CODED
! Call to Newton iteration from main code - CODED
! Check that timestepping produces some output - DOING NOW
! Extract fluxes from DG
! Vorticity calculation
! Mass lumping for P2b
! Visualisation of P2b by mapping back to P2
! Mass mapping from P1dg to P2b
! Timestepping for PV
! Nonlinear residual calculation from PV and DG advection
! stabilisation for PV
! discontinuity capturing for PV
! stabilisation for divergence
! Improve DG timestepping
! Improve slope limiter
! Visualisation of vorticity with bubbles
