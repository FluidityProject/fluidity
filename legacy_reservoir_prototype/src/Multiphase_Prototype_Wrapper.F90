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

  subroutine multiphase_prototype_wrapper(filename, filename_len)

    use fldebug
    use state_module
    use populate_state_module
    use diagnostic_variables
    use diagnostic_fields_wrapper
    use write_state_module
    use timeloop_utilities
    use timers
    use parallel_tools
    use global_parameters, only: current_time, dt, timestep, option_path_len, &
         simulation_start_time, &
         simulation_start_cpu_time, &
         simulation_start_wall_time
    use diagnostic_fields_new, only : &
         & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
         & check_diagnostic_dependencies
    use field_priority_lists
    use spud
    use checkpoint
    use boundary_conditions_from_options
    use discrete_properties_module
    use multiphase_module
    use multimaterial_module
    use field_equations_cv, only: initialise_advection_convergence
    use memory_diagnostics
    use mp_prototype
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

    type(state_type), dimension(:), pointer :: state

    integer, intent(in) :: filename_len
    integer :: ierr, i, dump_no = 0
    integer :: ntsol, nonlinear_iterations

    character(len = filename_len), intent(in) :: filename
    character(len = option_path_len) :: simulation_name, dump_format

    real :: finish_time, nonlinear_iteration_tolerance, steady_state_tolerance
!    integer :: mpi_success

#ifdef HAVE_MPI
    call mpi_init(ierr)
    assert(ierr == MPI_SUCCESS)
#endif

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    call set_simulation_start_times()
    call initialise_walltime
    timestep = 0

#ifdef HAVE_MEMORY_STATS
    ! this is to make sure the option /io/log_output/memory_diagnostics is read
    call reset_memory_logs()
#endif

    call initialise_write_state
    ! Read state from .flml file
    call populate_state(state)

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    !--------------------------------------------------------------------------------------------------------
    ! This should be read by the Copy_Outof_Into_State subrts

    ! set the remaining timestepping options, needs to be before any diagnostics are calculated
    call get_option("/timestepping/timestep", dt)
    !  if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
    !    call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
    !    call set_option("/timestepping/timestep", dt)
    !  end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)

    call get_option('/simulation_name',simulation_name)
    call initialise_diagnostics(trim(simulation_name),state)

    ! Calculate the number of scalar fields to solve for and their correct
    ! solve order taking into account dependencies.
    call get_ntsol(ntsol)

    call initialise_field_lists_from_options(state, ntsol)


    ! For multiphase simulations, we have to call calculate_diagnostic_phase_volume_fraction *before*
    ! copy_to_stored(state,"Old") is called below. Otherwise, OldPhaseVolumeFraction (in the phase
    ! containing the diagnostic PhaseVolumeFraction) will be zero and 
    ! NonlinearPhaseVolumeFraction will be calculated incorrectly at t=0.
    if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then
       call calculate_diagnostic_phase_volume_fraction(state)
    end if

    ! set the nonlinear timestepping options, needs to be before the adapt at first timestep
    call get_option('/timestepping/nonlinear_iterations',nonlinear_iterations,&
         & default=1)
    call get_option("/timestepping/nonlinear_iterations/tolerance", &
         & nonlinear_iteration_tolerance, default=0.0)

    ! Auxilliary fields.
    call allocate_and_insert_auxilliary_fields(state)
    call copy_to_stored_values(state,"Old")
    call copy_to_stored_values(state,"Iterated")
    call relax_to_nonlinear(state)

    call enforce_discrete_properties(state)

    call run_diagnostics(state)

    !     Determine the output format.
    call get_option('/io/dump_format', dump_format)
    if(trim(dump_format) /= "vtk") then
       ewrite(-1,*) "You must specify a dump format and it must be vtk."
       FLExit("Rejig your FLML: /io/dump_format")
    end if

    ! initialise the multimaterial fields
    call initialise_diagnostic_material_properties(state)

    call calculate_diagnostic_variables(state)
    call calculate_diagnostic_variables_new(state)


    !--------------------------------------------------------------------------------------------------------


    if( &
                                ! if this is not a zero timestep simulation (otherwise, there would
                                ! be two identical dump files)
         & current_time < finish_time &
                                ! unless explicitly disabled
         & .and. .not. have_option("/io/disable_dump_at_start") &
         & ) then
       call write_state(dump_no, state)
    end if

    call initialise_convergence(filename, state)
    call initialise_steady_state(filename, state)
    call initialise_advection_convergence(state)

    if(have_option("/io/stat/output_at_start")) then
      call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
    end if

    ! ******************************
    ! *** Start of timestep loop ***
    ! ******************************
    timestep_loop: do
       timestep = timestep + 1

       ewrite(1, *) "********************"
       ewrite(1, *) "*** NEW TIMESTEP ***"
       ewrite(1, *) "********************"
       ewrite(1, *) "Current simulation time: ", current_time
       ewrite(1, *) "Timestep number: ", timestep
       ewrite(1, *) "Timestep size (dt): ", dt
       !    if(.not. allfequals(dt)) then
       !       ewrite(-1, *) "Timestep size (dt): ", dt
       !       FLAbort("The timestep is not global across all processes!")
       !    end if

       if(simulation_completed(current_time, timestep)) exit timestep_loop

       if( &
                                ! Do not dump at the start of the simulation (this is handled by write_state call earlier)
            current_time > simulation_start_time &
                                ! Do not dump at the end of the simulation (this is handled by later write_state call)
            .and. current_time < finish_time &
                                ! Test write_state conditions
            .and. do_write_state(current_time, timestep) &
            ) then
          ! Intermediate dumps
          if(do_checkpoint_simulation(dump_no)) then
             call checkpoint_simulation(state, cp_no = dump_no)
          end if
          call write_state(dump_no, state)
       end if

       ewrite(2,*)'steady_state_tolerance,nonlinear_iterations:',steady_state_tolerance,nonlinear_iterations

       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)

       call enforce_discrete_properties(state, only_prescribed=.true., &
            exclude_interpolated=.true., &
            exclude_nonreprescribed=.true.)

       ! Call the multiphase_prototype code  
       call multiphase_prototype(state, dt, current_time, finish_time, &
                                 nonlinear_iterations, nonlinear_iteration_tolerance)

       current_time = current_time + dt + 100000.

    end do timestep_loop
    ! ****************************
    ! *** END OF TIMESTEP LOOP ***
    ! ****************************

    ! Dump at end, unless explicitly disabled
    if(.not. have_option("/io/disable_dump_at_end")) then
       call write_state(dump_no, state)
    end if

    ! Deallocate state
    do i = 1, size(state)
       call deallocate(state(i))
    end do
    deallocate(state)

  contains

    subroutine set_simulation_start_times()
      !!< Set the simulation start times

      call get_option("/timestepping/current_time", simulation_start_time)

      call cpu_time(simulation_start_cpu_time)
      call allmax(simulation_start_cpu_time)

      simulation_start_wall_time = wall_time()
      call allmax(simulation_start_wall_time)

    end subroutine set_simulation_start_times

  end subroutine multiphase_prototype_wrapper
