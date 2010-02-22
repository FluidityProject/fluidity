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

module fluids_module

  use AuxilaryOptions
  use MeshDiagnostics
  use signal_vars
  use spud
  use equation_of_state
  use timers
  use adapt_state_module 
  use adapt_state_prescribed_module
  use FLDebug
  use sparse_tools
  use elements
  use fields
  use boundary_conditions_from_options
  use populate_state_module
  use reserve_state_module
  use vtk_interfaces
  use Diagnostic_variables
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
  use diagnostic_fields_wrapper
  use diagnostic_children
  use advection_diffusion_cg
  use advection_diffusion_DG
  use advection_diffusion_FV
  use field_equations_cv, only: solve_field_eqn_cv, initialise_advection_convergence, coupled_cv_field_eqn
  use volumesource
  use vertical_extrapolation_module
  use qmesh_module
  use checkpoint
  use write_state_module
  use traffic  
  use synthetic_bc
  use goals
  use adaptive_timestepping
  use conformity_measurement
  ! Use Solid-fluid coupling and ALE - Julian- 18-09-06
  use ale_module
  use adjacency_lists
  use multimaterial_module
  use parallel_tools
  use SolidConfiguration
  use MeshMovement
  use write_triangle
  use biology
  use momentum_equation
  use timeloop_utilities
  use free_surface_module
  use field_priority_lists
  use boundary_conditions
  use porous_media
  use spontaneous_potentials, only: calculate_electrical_potential
  use discrete_properties_module
  use gls
  use halos
  use memory_diagnostics
  use global_parameters, only: current_time, dt, timestep, OPTION_PATH_LEN, topology_mesh_name
  
  implicit none

  private

  public :: fluids, fluids_module_check_options

  interface
    subroutine check_options
    end subroutine check_options
  end interface

contains

  SUBROUTINE FLUIDS(filename)
    character(len = *), intent(in) :: filename
    
    INTEGER :: &
         & NTSOL,  &
         & nonlinear_iterations

    REAL :: &
         & finish_time, &
         & steady_state_tolerance
         
    real:: nonlinear_iteration_tolerance

    !     System state wrapper.
    type(state_type), dimension(:), pointer :: state => null()
    type(tensor_field) :: metric_tensor
    !     Dump index
    integer :: dump_no = 0
    !     Temporary buffer for any string options which may be required.
    character(len=OPTION_PATH_LEN) :: option_buffer

    REAL :: CHANGE,CHAOLD

    integer :: i, it, its

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !     STUFF for MEsh movement, and Solid-fluid-coupling.  ------ jem

    !     Ale mesh movement - Julian 05-02-07
    LOGICAL:: USE_ALE
    INTEGER:: fs

    !     Solid-fluid coupling - Julian 18-09-06
    !New options
    INTEGER :: ss,ph
    LOGICAL :: have_solids

    ! Pointers for scalars and velocity fields
    type(vector_field), pointer :: Velocity
    type(scalar_field), pointer :: sfield
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !     backward compatibility with new option structure - crgw 21/12/07
    logical::use_advdif=.true.  ! decide whether we enter advdif or not

    INTEGER :: adapt_count

    ! Absolute first thing: check that the options, if present, are valid.
    call check_options
    ewrite(1,*) "Options sanity check successful"

    call set_simulation_start_times()
    timestep = 0

#ifdef HAVE_MEMORY_STATS
    ! this is to make sure the option /io/log_output/memory_diagnostics is read
    call reset_memory_logs()
#endif

    call initialise_qmesh
    call initialise_write_state

    if (have_option("/geometry/disable_geometric_data_cache")) then
       ewrite(1,*) "Disabling geometric data cache"
       cache_transform_elements=.false.
    end if

    adapt_count = 0

    ! Read state from .flml file
    call populate_state(state)
    
    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)
    
    if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) then
       if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, current_time, dt)

       call adapt_state_first_timestep(state)
   
       ! Auxilliary fields.
       call allocate_and_insert_auxilliary_fields(state)
       call copy_to_stored_values(state,"Old")
       call copy_to_stored_values(state,"Iterated")
       call relax_to_nonlinear(state)
       
       call enforce_discrete_properties(state)
       if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, current_time, dt)
    else
       ! Auxilliary fields.
       call allocate_and_insert_auxilliary_fields(state)
       call copy_to_stored_values(state,"Old")
       call copy_to_stored_values(state,"Iterated")
       call relax_to_nonlinear(state)

       call enforce_discrete_properties(state)
    end if

    ! Calculate the number of scalar fields to solve for and their correct
    ! solve order taking into account dependencies.
    call get_ntsol(ntsol)

    call initialise_field_lists_from_options(state, ntsol)

    call check_old_code_path()

    ! set all timestepping options, needs to be before any diagnostics are calculated
    
    call get_option("/timestepping/timestep", dt)
    if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
       call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
       call set_option("/timestepping/timestep", dt)
    end if

    call get_option('/timestepping/nonlinear_iterations',nonlinear_iterations,&
         & default=1)
    call get_option("/timestepping/nonlinear_iterations/tolerance", &
         & nonlinear_iteration_tolerance, default=0.0)
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)

    call get_option("/timestepping/steady_state/tolerance", &
         & steady_state_tolerance, default = -666.01)

    
    call run_diagnostics(state)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     Initialise solid-fluid coupling, and ALE ----------- -Julian 17-07-2008
    have_solids=.false.
    use_ale=.false.

    !Read the amount of SolidConcentration fields
    !and save their IT numbers.
    if(have_option(trim('/imported_solids'))) then
       ss=0
       phaseloop: do ph = 1, size(state)
          write(option_buffer, '(a,i0,a)') "/material_phase[",ph-1,"]"
          if(have_option(trim(option_buffer)//"/scalar_field::SolidConcentration")) then
             ss = ph
             have_solids=.true.
          end if
       end do phaseloop
       if(ss==0) then
          FLAbort("Havent found a material phase containing solid concentration")
       end if
       EWRITE(2,*) 'solid state= ',ss
    end if

    if(have_option(trim('/mesh_adaptivity/mesh_movement/explicit_ale'))) then
       !if using ALE with two fluids, look for the prognostic Material Volume fraction field.
       !This will later be changed to be more general (i.e.: be able to do this with any field by providing
       !its name)
       use_ale=.true.
       fs=-1
       phaseloop1: do ph = 1, size(state)
          write(option_buffer, '(a,i0,a)') "/material_phase[",ph-1,"]"
          if(have_option(trim(option_buffer)//"/scalar_field::MaterialVolumeFraction/prognostic")) then
             fs = ph
          end if
       end do phaseloop1
       if (fs==-1) then
          FLAbort('No prognostic MaterialVolumeFraction was defined')
       end if
    end if
    !     end initialise solid-fluid coupling, and ALE  -Julian 17-07-2008
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    if (have_option("/mesh_adaptivity/hr_adaptivity")) then
       call allocate(metric_tensor, extract_mesh(state(1), topology_mesh_name), "ErrorMetric")
    end if

    !     Determine the output format.
    call get_option('/io/dump_format', option_buffer)
    if(trim(option_buffer) /= "vtk") then
       FLAbort("You must specify a dump format and it must be vtk")
    end if

    !        Initialisation of distance to top and bottom field
    !        Currently only needed for free surface
    if (has_scalar_field(state(1), "DistanceToTop")) then
       if (.not. have_option('/geometry/ocean_boundaries')) then
          FLAbort("There are no top and bottom boundary markers.")
       end if
       call CalculateTopBottomDistance(state(1))
    end if

    ! initialise the multimaterial fields
    call initialise_diagnostic_material_properties(state)

    call calculate_diagnostic_variables(State)
    call calculate_diagnostic_variables_new(state)
    ! This is mostly to ensure that the photosynthetic radiation
    ! has a non-zero value before the first adapt. 
    if (have_option("/ocean_biology")) then
       call calculate_biology_terms(state(1))
    end if

    ! Checkpoint at start
    if(do_checkpoint_simulation(dump_no)) call checkpoint_simulation(state, cp_no = dump_no)
    ! Dump at start  
    if( &
         ! if this is not a zero timestep simulation (otherwise, there would
         ! be two identical dump files)
         & current_time < finish_time &
         ! unless explicitly disabled
         & .and. .not. have_option("/io/disable_dump_at_start") &
         & ) then
       call write_state(dump_no, state)
    end if
    
    call initialise_diagnostics(filename, state)
    call initialise_convergence(filename, state)
    call initialise_advection_convergence(state)
    if(have_option("/io/stat/output_at_start")) call write_diagnostics(state, current_time, dt)
    
    do i=1, size(state)
       velocity => extract_vector_field(state(i), "Velocity")
       if (has_boundary_condition(velocity, "free_surface") .and. &
            & have_option('/mesh_adaptivity/mesh_movement/free_surface') .and. &
            & .not.aliased(velocity)) then
          ewrite(1,*) "Going into move_free_surface_nodes to compute surface node coordinates from initial condition"
          call move_free_surface_nodes(state(i))
       end if
    end do

    ! Initialise GLS
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/option")) then
        call gls_init(state(1))
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
       if(.not. allfequals(dt)) then
          ewrite(-1, *) "Timestep size (dt): ", dt
          FLAbort("The timestep is not global across all processes!")
       end if

       if(simulation_completed(current_time, timestep)) exit timestep_loop

       if( &
                                ! Do not dump at the start of the simulation (this is handled by write_state call earlier)
            & current_time > simulation_start_time &
                                ! Do not dump at the end of the simulation (this is handled by later write_state call)
            & .and. current_time < finish_time &
                                ! Test write_state conditions
            & .and. do_write_state(current_time, timestep) &
            & ) then

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! Regular during run state dump.
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          ! Intermediate dumps
          if(do_checkpoint_simulation(dump_no)) then
             call checkpoint_simulation(state, cp_no = dump_no)
          end if
          call write_state(dump_no, state)
       end if

       ewrite(2,*)'steady_state_tolerance,nonlinear_iterations:',steady_state_tolerance,nonlinear_iterations

       call copy_to_stored_values(state,"Old")
       if (have_option('/mesh_adaptivity/mesh_movement')) then 
          ! Coordinate isn't handled by the standard timeloop utility calls so instead
          ! we have to make sure the (Old)Coordinate field is updated with the current, most up
          ! to date Coordinates, which actually lie in IteratedCoordinate (since if the mesh
          ! is moving Coordinate is evaluated at time level n+theta).
          ! Using state(1) should be safe as they are aliased across all states.
          call set_vector_field_in_state(state(1), "OldCoordinate", "IteratedCoordinate")
          call set_vector_field_in_state(state(1), "Coordinate", "IteratedCoordinate")
          ! same for GridVelocity
          call set_vector_field_in_state(state(1), "OldGridVelocity", "IteratedGridVelocity")
          call set_vector_field_in_state(state(1), "GridVelocity", "IteratedGridVelocity")
       end if

       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)
       
       ! evaluate prescribed fields at time = current_time+dt
       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time+dt)
       call enforce_discrete_properties(state, only_prescribed=.true., &
            exclude_interpolated=.true., &
            exclude_nonreprescribed=.true.)

       ! nonlinear_iterations=maximum no of iterations within a time step
       ! NB TEMPT is T from previous iteration.
       nonlinear_iteration_loop: do  ITS=1,nonlinear_iterations

          ewrite(1,*)'###################'
          ewrite(1,*)'Start of another nonlinear iteration; ITS,nonlinear_iterations=',ITS,nonlinear_iterations
          ewrite(1,*)'###################'

          call copy_to_stored_values(state, "Iterated")
          ! relax to nonlinear has to come before copy_from_stored_values
          ! so that the up to date values don't get wiped from the field itself
          ! (this could be fixed by replacing relax_to_nonlinear on the field
          !  with a dependency on the iterated field but then copy_to_stored_values
          !  would have to come before relax_to_nonlinear)
          call relax_to_nonlinear(state)
          call copy_from_stored_values(state, "Old")

          call compute_goals(state)

          !------------------------------------------------
          ! Addition for calculating drag force ------ jem 05-06-2008
          if (have_option("/imported_solids/calculate_drag_on_surface")) then
             call drag_on_surface(state)
          end if

          !     Addition for reading solids in - jem  02-04-2008
          if(have_solids) call solid_configuration(state(ss:ss),its,nonlinear_iterations)

          !Explicit ALE ------------   jem 21/07/08         
          if (use_ale) then 
             EWRITE(0,'(A)') '----------------------------------------------'
             EWRITE(0,'(A26,E12.6)') 'Using explicit_ale. time: ',current_time       
             call explicit_ale(state,fs)
          end if
          !end explicit ale ------------  jem 21/07/08

          !-----------------------------------------------------
          ! Call to porous_media module (leading to multiphase flow in porous media)
          ! jhs - 16/01/09
          if (have_option("/porous_media")) then
             call porous_media_advection(state)
             ! compute spontaneous electrical potentials (myg - 28/10/09)
             do i=1,size(state)
                option_buffer = '/material_phase['//int2str(i-1)//']/electrical_properties/coupling_coefficients/'
                if (have_option(trim(option_buffer)//'scalar_field::Electrokinetic').or.&
                    have_option(trim(option_buffer)//'scalar_field::Thermoelectric').or.&
                    have_option(trim(option_buffer)//'scalar_field::Electrochemical')) then
                   call calculate_electrical_potential(state(i), i)
                end if
             end do
          end if
          ! End call to porous_media

          if (have_option("/ocean_biology")) then
             call calculate_biology_terms(state(1))
          end if

          field_loop: do it = 1, ntsol
             ewrite(2, "(a,i0,a,i0)") "Considering scalar field ", it, " of ", ntsol
             ewrite(1, *) "Considering scalar field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name)

             ! do we have the generic length scale vertical turbulence model?
             if( have_option("/material_phase[0]/subgridscale_parameterisations/GLS/option") ) then
                if( (trim(field_name_list(it))=="GLSTurbulentKineticEnergy")) then
                    call gls_tke(state(1))
                else if( (trim(field_name_list(it))=="GLSGenericSecondQuantity")) then
                    call gls_psi(state(1))
                endif
              end if


             call get_option(trim(field_optionpath_list(it))//&
                  '/prognostic/equation[0]/name', &
                  option_buffer, default="UnknownEquationType")
             select case(trim(option_buffer))
             case ( "AdvectionDiffusion", "ConservationOfMass", "ReducedConservationOfMass", "InternalEnergy" )
                use_advdif=.true.
             case default
                use_advdif=.false.
             end select

             IF(use_advdif)THEN

                if(starts_with(trim(field_name_list(it)), "TrafficTracer")) then
                   call traffic_tracer(trim(field_name_list(it)),state(field_state_list(it)),timestep)
                endif
                
                sfield => extract_scalar_field(state(field_state_list(it)), field_name_list(it))
                call calculate_diagnostic_children(state, field_state_list(it), sfield)


                !--------------------------------------------------
                !This addition creates a field that is a copy of
                !another to be used, i.e.: for diffusing.
                call get_copied_field(field_name_list(it), state(field_state_list(it)))
                !--------------------------------------------------

                IF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/discontinuous_galerkin")) then

                   ! Solve the DG form of the equations.
                   call solve_advection_diffusion_dg(field_name=field_name_list(it), &
                        & state=state(field_state_list(it)))

                ELSEIF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/finite_volume")) then

                   ! Solve the FV form of the equations.
                   call solve_advection_diffusion_fv(field_name=field_name_list(it), &
                        & state=state(field_state_list(it)))

                ELSEIF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/control_volumes")) then

                   ! Solve the pure control volume form of the equations
                   call solve_field_eqn_cv(field_name=trim(field_name_list(it)), &
                        state=state(field_state_list(it):field_state_list(it)), &
                        global_it=its)

                else if(have_option(trim(field_optionpath_list(it)) // &
                     & "/prognostic/spatial_discretisation/continuous_galerkin")) then

                   call solve_field_equation_cg(field_name_list(it), state(field_state_list(it)), dt)
                else

                   ewrite(2, *) "Not solving scalar field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name) //" in an advdif-like subroutine."

                end if ! End of dg/cv/cg choice.

                ! ENDOF IF((TELEDI(IT).EQ.1).AND.D3) THEN ELSE...
             ENDIF

             ewrite(1, *) "Finished field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name)
          end do field_loop

          ! Sort out the dregs of GLS after the solve on Psi (GenericSecondQuantity) has finished
          if( have_option("/material_phase[0]/subgridscale_parameterisations/GLS/option") ) then
            ! Update the diffusivity, only at the end of the loop ready for
            ! the next timestep. We do NOT want to be twiddling
            ! diffusivity/viscosity half way through a non-linear iteration
            call gls_diffusivity(state(1))
          end if


          if(option_count("/material_phase/scalar_field/prognostic/spatial_discretisation/coupled_cv")>0) then
             call coupled_cv_field_eqn(state, global_it=its)
          end if
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          !
          ! Assemble and solve N.S equations.
          !

          !-------------------------------------------------------------
          ! Call to porous_media_momentum (leading to multiphase)
          ! jhs - 16/01/09
          ! moved to here 04/02/09
          if (have_option("/porous_media")) then
             call porous_media_momentum(state)
          end if
          
          if (have_option("/traffic_model")) then
             call traffic_source(state(1),timestep)
             call traffic_density_update(state(1))
          endif
          
          ! This is where the non-legacy momentum stuff happens
          ! a loop over state (hence over phases) is incorporated into this subroutine call
          ! hence this lives outside the phase_loop
          call momentum_loop(state, at_first_timestep=((timestep==1).and.(its==1)))
                       
          if(nonlinear_iterations > 1) then
             ! Check for convergence between non linear iteration loops
             call test_and_write_convergence(state, current_time + dt, dt, its, change)
             if(its == 1) chaold = change             

             if (have_option("/timestepping/nonlinear_iterations/&
                              &tolerance")) then
               ewrite(2, *) "Nonlinear iteration change = ", change
               ewrite(2, *) "Nonlinear iteration tolerance = ", nonlinear_iteration_tolerance

               if(change < abs(nonlinear_iteration_tolerance)) then
                  ewrite(1, *) "Nonlinear iteration tolerance has been reached"
                  ewrite(1, "(a,i0,a)") "Exiting nonlinear iteration loop after ", its, " iterations"
                  exit nonlinear_iteration_loop
               endif
             end if
          end if

       end do nonlinear_iteration_loop

       if(have_option("/timestepping/nonlinear_iterations/terminate_if_not_converged")) then
          if(its >= nonlinear_iterations .and. change >= abs(nonlinear_iteration_tolerance)) then
             ewrite(0, *) "Nonlinear iteration tolerance not reached - termininating"
             exit timestep_loop
          end if
       end if

       if(have_option(trim('/mesh_adaptivity/mesh_movement/vertical_ale'))) then
          ewrite(1,*) 'Entering vertical_ale routine'
          !move the mesh and calculate the grid velocity
          call movemeshy(state(1))
       endif

       current_time=current_time+DT

       ! calculate and write diagnostics before the timestep gets changed
       call calculate_diagnostic_variables(State, exclude_nonrecalculated=.true.)
       call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)

       ! Call the modern and significantly less satanic version of study
       call write_diagnostics(state, current_time, dt)

       if(have_option("/timestepping/adaptive_timestep")) call calc_cflnumber_field_based_dt(state, dt)

       ! Update the options dictionary for the new timestep and current_time.
       call set_option("/timestepping/timestep", dt)
       call set_option("/timestepping/current_time", current_time)

       ! if strong bc or weak that overwrite then enforce the bc on the fields
       ! (should only do something for weak bcs with that options switched on)
       call set_dirichlet_consistent(state)

       if(have_option("/timestepping/steady_state")) then

          call test_steady_state(state, change)
          if(change<steady_state_tolerance) then
             ewrite(0,*)  "* Steady state has been attained, exiting the timestep loop"
             exit timestep_loop
          end if

       end if

       if(simulation_completed(current_time)) exit timestep_loop
       
       ! ******************
       ! *** Mesh adapt ***
       ! ******************
       if(have_option("/mesh_adaptivity/hr_adaptivity")) then
          if(do_adapt_mesh(current_time, timestep)) then  

             call pre_adapt_tasks()

             call qmesh(state, metric_tensor)

             if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, current_time, dt)
             call run_diagnostics(state)
             
             call adapt_state(state, metric_tensor)
             call update_state_post_adapt(state, metric_tensor, dt)
             
             if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, current_time, dt)
             call run_diagnostics(state)
          end if
       else if(have_option("/mesh_adaptivity/prescribed_adaptivity")) then
          if(do_adapt_state_prescribed(current_time)) then    
              
              call pre_adapt_tasks()
      
             if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, current_time, dt)             
             call run_diagnostics(state)
             
             call adapt_state_prescribed(state, current_time)
             call update_state_post_adapt(state, metric_tensor, dt)
             
             if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, current_time, dt)
             call run_diagnostics(state)
          end if
       end if

    end do timestep_loop
    ! ****************************
    ! *** END OF TIMESTEP LOOP ***
    ! ****************************

    ! Checkpoint at end, if enabled
    if(have_option("/io/checkpointing/checkpoint_at_end")) then
       call checkpoint_simulation(state, cp_no = dump_no)
    end if
    ! Dump at end, unless explicitly disabled
    if(.not. have_option("/io/disable_dump_at_end")) then
       call write_state(dump_no, state)
    end if

    ! cleanup GLS
    if (have_option('/material_phase[0]/subgridscale_parameterisations/GLS/')) then
        call gls_cleanup()
    end if

    ! closing .stat, .convergence and .detector files
    call close_diagnostic_files()

    ewrite(1, *) "Printing references before final deallocation"
    call print_references(1)

    ! Deallocate the metric tensor
    if(have_option("/mesh_adaptivity/hr_adaptivity")) call deallocate(metric_tensor)
    ! Deallocate state
    do i = 1, size(state)
       call deallocate(state(i))
    end do
    ! Deallocate the reserve state
    call deallocate_reserve_state()

    ! deallocate the pointer to the array of states
    deallocate(state)

    ! Delete the transform_elements cache.
    call deallocate_transform_cache

    ewrite(1, *) "Printing references after final deallocation"
    call print_references(0)

#ifdef HAVE_MEMORY_STATS
    call print_current_memory_stats(0)
#endif

  contains

    subroutine set_simulation_start_times()
      !!< Set the simulation start times

      call get_option("/timestepping/current_time", simulation_start_time)

      call cpu_time(simulation_start_cpu_time)
      call allmax(simulation_start_cpu_time)

      simulation_start_wall_time = wall_time()
      call allmax(simulation_start_wall_time)

    end subroutine set_simulation_start_times

  end subroutine fluids
 
  subroutine pre_adapt_tasks()

    ! GLS - we need to deallocate all module-level fields or the memory
    ! management system complains
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/")) then
        call gls_cleanup() ! deallocate everything
    end if

  end subroutine pre_adapt_tasks
 
  subroutine update_state_post_adapt(state, metric_tensor, dt)
    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(out) :: metric_tensor
    real, intent(inout) :: dt
    
    integer :: i
    type(vector_field), pointer :: velocity
    
    ! The adaptivity metric
    if(have_option("/mesh_adaptivity/hr_adaptivity")) then
      call allocate(metric_tensor, extract_mesh(state(1), "CoordinateMesh"), "ErrorMetric")
    end if
    
    ! Auxilliary fields.
    call allocate_and_insert_auxilliary_fields(state)
    call copy_to_stored_values(state,"Old")
    call copy_to_stored_values(state,"Iterated")
    call relax_to_nonlinear(state)

    ! Discrete properties
    call enforce_discrete_properties(state)

    ! Timestep adapt
    if(have_option("/timestepping/adaptive_timestep")) then
      call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
      call set_option("/timestepping/timestep", dt)
    end if
    
    ! Ocean boundaries
    if (has_scalar_field(state(1), "DistanceToTop")) then
      if (.not. have_option('/geometry/ocean_boundaries')) then
         FLAbort("There are no top and bottom boundary markers.")
      end if
      call CalculateTopBottomDistance(state(1))
    end if
        
    ! Diagnostic fields
    call calculate_diagnostic_variables(state)
    call calculate_diagnostic_variables_new(state)
    ! This is mostly to ensure that the photosynthetic radiation
    ! has a non-zero value before the next adapt. 
    if(have_option("/ocean_biology")) call calculate_biology_terms(state(1))
   
    ! Free surface movement
    do i=1, size(state)
      velocity => extract_vector_field(state(i), "Velocity")
      if (has_boundary_condition(velocity, "free_surface") .and. &
           & have_option('/mesh_adaptivity/mesh_movement/free_surface') .and. &
           & .not.aliased(velocity)) then
         ewrite(1,*) "Going into move_free_surface_nodes to compute surface node coordinates after adapt"
         call move_free_surface_nodes(state(i))
      end if
    end do

    ! GLS
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/")) then
        call gls_adapt_mesh(state(1))
    end if
           
  end subroutine update_state_post_adapt
  
  subroutine fluids_module_check_options
    !!< Check fluids specific options

    integer :: stat
    real :: time_limit

    ewrite(2, *) "Checking simulation completion options"

    call get_option("/timestepping/wall_time_limit", time_limit, stat)
    if(stat == SPUD_NO_ERROR) then
      if(time_limit < 0.0) then
        FLExit("Wall time limit cannot be negative")
      end if
      if(.not. wall_time_supported()) then
        FLExit("Wall time limit supplied, but wall time is not available")
      end if
    end if

    ewrite(2, *) "Finished checking simulation completion options"

  end subroutine fluids_module_check_options
    
  subroutine check_old_code_path()
  !!< checks whether any of the phases use the old code path
    character(len=OPTION_PATH_LEN):: option_path
    character(len=FIELD_NAME_LEN):: tmpstring
    logical:: use_advdif
    integer:: i, j, tmpstat
    
    ! first check for a velocity field with legacy options
    do i=1, option_count("/material_phase")
      option_path="/material_phase["//int2str(i-1)//"]/vector_field::Velocity"
      if( have_option(trim(option_path)//"/prognostic/spatial_discretisation&
                                &/legacy_continuous_galerkin") &
          .or. &
          have_option(trim(option_path)//"/prognostic/spatial_discretisation&
                                &/legacy_discretisation") &
        ) then
    
        ewrite(0,*) "ERROR: You seem to be using legacy_continuous_galerkin or "
        ewrite(0,*) "legacy_discretisation for spatial_discretisation of Velocity."
        ewrite(0,*) "This uses the old code path (navsto) which has been deleted."
        ewrite(0,*) "You should switch to continuous_galerkin."
        FLExit("The old code path is dead.")
      end if
    end do
    
    do i=1, option_count("/material_phase")
      do j=1, option_count("/material_phase["//int2str(i-1)//"]/scalar_field")
        option_path="/material_phase["//int2str(i-1)//"]/scalar_field["//int2str(j-1)//']'
        
        ! this is a copy from fluids() above:
        call get_option(trim(option_path)//&
              '/prognostic/equation[0]/name', &
              tmpstring, stat=tmpstat)
        if (tmpstat==0) then
          select case(trim(tmpstring))
          case ( "AdvectionDiffusion", "ConservationOfMass", "ReducedConservationOfMass", "InternalEnergy" )
            use_advdif=.true.
          case default
            use_advdif=.false.
          end select
        else
          use_advdif=.false.
        end if
        
        if (use_advdif .and. ( &
          have_option(trim(option_path)//&
            & "/prognostic/spatial_discretisation/legacy_continuous_galerkin").or.&
          have_option(trim(option_path)//&
            & "/prognostic/spatial_discretisation/legacy_mixed_cv_cg").or.&
          have_option(trim(option_path)//&
            & "/prognostic/spatial_discretisation/legacy_discretisation") &
          )) then
            
          ewrite(0,*) "Error: You seem to be using legacy_continuous_galerkin,"
          ewrite(0,*) "legacy_mixed_cv_cg or legacy_discretisation for the"
          ewrite(0,*) "spatial discretisation of one of your scalar fields."
          ewrite(0,*) "This uses the old code path (advdif) that has been deleted"
          ewrite(0,*) "You should switch to any of the other options."
          FLExit("The old code path is dead.")
        end if
         
      end do
    end do
    
  end subroutine check_old_code_path
  
  end module fluids_module

