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

module fluids_module

  use fldebug
  use auxilaryoptions
  use spud
  use global_parameters, only: current_time, dt, timestep, OPTION_PATH_LEN, &
       simulation_start_time, &
       simulation_start_cpu_time, &
       simulation_start_wall_time, &
       topology_mesh_name, FIELD_NAME_LEN
  use futils, only: int2str
  use reference_counting, only: print_references
  use parallel_tools
  use memory_diagnostics
  use sparse_tools
  use elements
  use adjacency_lists
  use eventcounter
  use transform_elements, only: cache_transform_elements, deallocate_transform_cache
  use meshdiagnostics
  use signal_vars
  use fields
  use state_module
  use vtk_interfaces
  use boundary_conditions
  use halos
  use equation_of_state
  use timers
  use synthetic_bc
  use k_epsilon, only: keps_advdif_diagnostics
  use tictoc
  use boundary_conditions_from_options
  use reserve_state_module
  use write_state_module
  use detector_parallel, only: sync_detector_coordinates, deallocate_detector_list_array
  use diagnostic_variables
  use populate_state_module
  use vertical_extrapolation_module
  use field_priority_lists
  use multiphase_module
  use multimaterial_module
  use free_surface_module
  use momentum_diagnostic_fields, only: calculate_densities
  use sediment_diagnostics, only: calculate_sediment_flux
  use dqmom
  use diagnostic_fields_wrapper
  use checkpoint
  use goals
  use adaptive_timestepping
  use conformity_measurement
  use timeloop_utilities
  use discrete_properties_module
  use adapt_state_module
  use adapt_state_prescribed_module
  use populate_sub_state_module
  use diagnostic_fields_new, only : &
       & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
       & check_diagnostic_dependencies
  use diagnostic_children
  use advection_diffusion_cg
  use advection_diffusion_dg
  use advection_diffusion_fv
  use field_equations_cv, only: solve_field_eqn_cv, initialise_advection_convergence, coupled_cv_field_eqn
  use qmesh_module
  use write_triangle
  use meshmovement
  use biology
  use foam_flow_module, only: calculate_potential_flow, calculate_foam_velocity
  use momentum_equation
  use gls
  use iceshelf_meltrate_surf_normal
#ifdef HAVE_HYPERLIGHT
  use hyperlight
#endif

  implicit none

  private

  public :: fluids, fluids_module_check_options

  interface
    subroutine check_options
    end subroutine check_options
  end interface

contains

  SUBROUTINE FLUIDS()
    character(len = OPTION_PATH_LEN) :: filename

    INTEGER :: &
         & NTSOL,  &
         & nonlinear_iterations,  &
         & nonlinear_iterations_adapt

    REAL :: &
         & finish_time, &
         & steady_state_tolerance

    real:: nonlinear_iteration_tolerance

    !     System state wrapper.
    type(state_type), dimension(:), pointer :: state => null()
    type(state_type), dimension(:), pointer :: sub_state => null()

    type(tensor_field) :: metric_tensor
    !     Dump index
    integer :: dump_no = 0
    !     Temporary buffer for any string options which may be required.
    character(len=OPTION_PATH_LEN) :: option_buffer
    character(len=OPTION_PATH_LEN):: option_path
    REAL :: CHANGE,CHAOLD

    integer :: i, it, its

    logical :: not_to_move_det_yet = .false.

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    ! An array of submaterials of the current phase in state(istate).
    ! Needed for k-epsilon VelocityBuoyancyDensity calculation line:~630
    ! S Parkinson 31-08-12
    type(state_type), dimension(:), pointer :: submaterials     

    ! Pointers for scalars and velocity fields
    type(scalar_field), pointer :: sfield
    type(scalar_field) :: foam_velocity_potential
    type(vector_field), pointer :: foamvel
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !     backward compatibility with new option structure - crgw 21/12/07
    logical::use_advdif=.true.  ! decide whether we enter advdif or not

    INTEGER :: adapt_count

    ! Absolute first thing: check that the options, if present, are valid.
    call check_options
    ewrite(1,*) "Options sanity check successful"

    call get_option("/simulation_name",filename)

    call set_simulation_start_times()
    call initialise_walltime
    timestep = 0

#ifdef HAVE_MEMORY_STATS
    ! this is to make sure the option /io/log_output/memory_diagnostics is read
    call reset_memory_logs()
#endif

    call initialise_qmesh
    call initialise_write_state

    ! Initialise Hyperlight
#ifdef HAVE_HYPERLIGHT
    if (have_option("ocean_biology/lagrangian_ensemble/hyperlight")) then
       if (.not.have_option("/material_phase[0]/scalar_field::Chlorophyll")) then
          FLExit("You need Chlorophyll scalar field for Hyperlight")
       else
          call hyperlight_init()       
       end if
    end if
#else
    if (have_option("ocean_biology/lagrangian_ensemble/hyperlight")) then
       ewrite(-1,*) "Hyperlight module was selected, but not compiled."
       FLExit("Please re-compile fluidity with the --enable-hyperlight option.")
    end if
#endif

    if (have_option("/geometry/disable_geometric_data_cache")) then
       ewrite(1,*) "Disabling geometric data cache"
       cache_transform_elements=.false.
    end if

    adapt_count = 0

    ! Read state from .flml file
    call populate_state(state)

    ewrite(3,*)'before have_option test'

    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)

    default_stat%zoltan_drive_call=.false.

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
    
    if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) then

       if(have_option("/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt")) then
         call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt',nonlinear_iterations_adapt)
         nonlinear_iterations = nonlinear_iterations_adapt
       end if
      
       ! set population balance initial conditions - for first adaptivity
       call dqmom_init(state)

       call adapt_state_first_timestep(state)

       ! Auxilliary fields.
       call allocate_and_insert_auxilliary_fields(state)
       call copy_to_stored_values(state,"Old")
       call copy_to_stored_values(state,"Iterated")
       call relax_to_nonlinear(state)

       call enforce_discrete_properties(state)

       ! Ensure that checkpoints do not adapt at first timestep.
       call delete_option(&
            "/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")
    else
       ! Auxilliary fields.
       call allocate_and_insert_auxilliary_fields(state)
       call copy_to_stored_values(state,"Old")
       call copy_to_stored_values(state,"Iterated")
       call relax_to_nonlinear(state)

       call enforce_discrete_properties(state)
    end if

   ! set the remaining timestepping options, needs to be before any diagnostics are calculated
    call get_option("/timestepping/timestep", dt)
    if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
       call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
       call set_option("/timestepping/timestep", dt)
    end if

    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/finish_time", finish_time)

    call get_option("/timestepping/steady_state/tolerance", &
         & steady_state_tolerance, default = -666.01)

    if(use_sub_state()) then
       call populate_sub_state(state,sub_state)
    end if

    ! set population balance initial conditions
    call dqmom_init(state)

    ! Calculate the number of scalar fields to solve for and their correct
    ! solve order taking into account dependencies.
    call get_ntsol(ntsol)

    call initialise_field_lists_from_options(state, ntsol)

    call check_old_code_path()

    !        Initialisation of distance to top and bottom field
    !        Currently only needed for free surface
    if (has_scalar_field(state(1), "DistanceToTop")) then
       if (.not. have_option('/geometry/ocean_boundaries')) then
          ewrite(-1,*) "Warning: You have a field called DistanceToTop"
          ewrite(-1,*) "but you don't have ocean_boundaries switched on."
       else
          call CalculateTopBottomDistance(state(1))
          ! Initialise the OriginalDistanceToBottom field used for wetting and drying
          if (have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")) then
             call insert_original_distance_to_bottom(state(1))
             ! Wetting and drying only works with no poisson guess ... let's check that
             call get_option("/material_phase[0]/scalar_field::Pressure/prognostic/scheme/poisson_pressure_solution", option_buffer)
             if (.not. trim(option_buffer) == "never") then 
               FLExit("Please choose 'never' under /material_phase[0]/scalar_field::Pressure/prognostic/scheme/poisson_pressure_solution when using wetting and drying")
             end if
          end if
       end if
    end if

    ! move mesh according to inital free surface:
    !    top/bottom distance needs to be up-to-date before this call, after the movement
    !    they will be updated (inside the call)
    call move_mesh_free_surface(state, initialise=.true.)

    call run_diagnostics(state)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    if (have_option("/mesh_adaptivity/hr_adaptivity")) then
       call allocate(metric_tensor, extract_mesh(state(1), topology_mesh_name), "ErrorMetric")
    end if

    !     Determine the output format.
    call get_option('/io/dump_format', option_buffer)
    if(trim(option_buffer) /= "vtk") then
       ewrite(-1,*) "You must specify a dump format and it must be vtk."
       FLExit("Rejig your FLML: /io/dump_format")
    end if

    ! Initialise multimaterial fields:
    call initialise_diagnostic_material_properties(state)

    ! Calculate diagnostic variables:
    call calculate_diagnostic_variables(state)
    call calculate_diagnostic_variables_new(state)
    
    ! This is mostly to ensure that the photosynthetic radiation
    ! has a non-zero value before the first adapt.
    if (have_option("/ocean_biology")) then
       call calculate_biology_terms(state(1))
    end if

    call initialise_diagnostics(filename, state)

    ! Initialise ice_meltrate, read constatns, allocate surface, and calculate melt rate
    if (have_option("/ocean_forcing/iceshelf_meltrate/Holland08")) then
        call melt_surf_init(state(1))
        call melt_allocate_surface(state(1))
        call melt_surf_calc(state(1))
         !BC for ice melt
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/calculate_boundaries')) then
            call melt_bc(state(1))
          endif
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

    call initialise_convergence(filename, state)
    call initialise_steady_state(filename, state)
    call initialise_advection_convergence(state)

    if(have_option("/io/stat/output_at_start")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)

    not_to_move_det_yet=.false.

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

       call tic(TICTOC_ID_TIMESTEP)

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
       if (have_option('/mesh_adaptivity/mesh_movement') .and. .not. have_option('/mesh_adaptivity/mesh_movement/free_surface')) then
          ! Coordinate isn't handled by the standard timeloop utility calls.
          ! During the nonlinear iterations of a timestep, Coordinate is
          ! evaluated at n+theta, i.e. (1-theta)*OldCoordinate+theta*IteratedCoordinate.
          ! At the end of the previous timestep however, the most up-to-date Coordinate,
          ! i.e. IteratedCoordinate, has been copied into Coordinate. This value
          ! is now used as the Coordinate at the beginning of the time step.
          ! For the free surface this is dealt with within move_mesh_free_surface() below
          call set_vector_field_in_state(state(1), "OldCoordinate", "Coordinate")
       end if
       ! if we're using an implicit (prognostic) viscous free surface then there will be a surface field stored
       ! under the boundary condition _implicit_free_surface on the FreeSurface field that we need to update
       ! - it has an old and a new timelevel and the old one needs to be set to the now old new values.
       call update_implicit_scaled_free_surface(state)

       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)
       if(use_sub_state()) call set_boundary_conditions_values(sub_state, shift_time=.true.)

       ! evaluate prescribed fields at time = current_time+dt
       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true., time=current_time+dt)

       if(use_sub_state()) call set_full_domain_prescribed_fields(state,time=current_time+dt)

       ! move the mesh according to a prescribed grid velocity
       ! NOTE: there may be a chicken and egg situation here.  This update
       ! has to come after the setting of the prescribed fields so that
       ! the grid velocity is not advected (in a lagrangian way) with the mesh
       ! after the previous timestep of mesh movement.  However, this means
       ! that other prescribed fields are actually set up according to an old
       ! Coordinate field.
       ! NOTE ALSO: this must come before the enforcement of discrete properties
       ! to ensure those properties are satisfied on the new mesh not the old one.
       call move_mesh_imposed_velocity(state)
       call move_mesh_pseudo_lagrangian(state)

       call enforce_discrete_properties(state, only_prescribed=.true., &
            exclude_interpolated=.true., &
            exclude_nonreprescribed=.true.)

#ifdef HAVE_HYPERLIGHT
       ! Calculate multispectral irradiance fields from hyperlight
       if(have_option("/ocean_biology/lagrangian_ensemble/hyperlight")) then
          call set_irradiance_from_hyperlight(state(1))
       end if
#endif

       ! nonlinear_iterations=maximum no of iterations within a time step

       nonlinear_iteration_loop: do  ITS=1,nonlinear_iterations

          ewrite(1,*)'###################'
          ewrite(1,*)'Start of another nonlinear iteration; ITS,nonlinear_iterations=',ITS,nonlinear_iterations
          ewrite(1,*)'###################'

          ! For each field, set the iterated field, if present:
          call copy_to_stored_values(state, "Iterated")
          ! For each field, set the nonlinear field, if present:
          call relax_to_nonlinear(state)
          call copy_from_stored_values(state, "Old")

          ! move the mesh according to the free surface algorithm
          ! this should not be at the end of the nonlinear iteration:
          ! if nonlinear_iterations==1:
          !    OldCoordinate is based on p^{n-1}, as it has been moved at the beginning of the previous timestep
          !    IteratedCoordinate will be moved to p^n, the current pressure achieved at the end of the previous timestep
          !    GridVelocity=(IteratedCoordinate-OldCoordinate)/dt
          !    so the coordinates and grid velocity are lagging one timestep behind the computed p, which is inevitable
          ! if nonlinear_iteration>1:
          !    In the first nonlinear iteration we use the same values of OldCoordinate (based on p^{n-1})
          !    and IteratedCoordinate (based on p at the beginning of last iteration of previous timestep, say p^n*)
          !    In subsequent iterations, we have a reasonable approximation of p^{n+1}, so we base
          !    OldCoordinate on p^n* and IteratedCoordinate on p^(n+1)*, the best approximation p^{n+1} thus
          !    far. Note that OldCoordinate should not be based on p^n (the value of p at the end of the
          !    last iteration of the previous timestep) but on p^n* (the value at the beginning of the last iteration
          !    of previous timestep).

          if (nonlinear_iterations==1) then
            call move_mesh_free_surface(state)
          else
            call move_mesh_free_surface(state, nonlinear_iteration=its)
          end if

          call compute_goals(state)

          ! Calculate source terms for population balance scalars
          call dqmom_calculate_source_terms(state, ITS)

          if (have_option("/ocean_biology")) then
             call calculate_biology_terms(state(1))
          end if

          ! Do we have the k-epsilon turbulence model?
          ! If we do then we want to calculate source terms and diffusivity for the k and epsilon 
          ! fields and also tracer field diffusivities at n + theta_nl
          do i= 1, size(state)
             if(have_option("/material_phase["//&
                  int2str(i-1)//"]/subgridscale_parameterisations/k-epsilon")) then
                if(timestep == 1 .and. its == 1 .and. have_option('/physical_parameters/gravity')) then
                   ! The very first time k-epsilon is called, VelocityBuoyancyDensity
                   ! is set to zero until calculate_densities is called in the momentum equation
                   ! solve. Calling calculate_densities here is a work-around for this problem.  
                   sfield => extract_scalar_field(state, 'VelocityBuoyancyDensity')
                   if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then 
                      call get_phase_submaterials(state, i, submaterials)
                      call calculate_densities(submaterials, buoyancy_density=sfield)
                      deallocate(submaterials)
                   else
                      call calculate_densities(state, buoyancy_density=sfield)
                   end if
                   ewrite_minmax(sfield)
                end if
                call keps_advdif_diagnostics(state(i))
             end if
          end do

          field_loop: do it = 1, ntsol
             ewrite(2, "(a,i0,a,i0)") "Considering scalar field ", it, " of ", ntsol
             ewrite(1, *) "Considering scalar field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name)

             ! do we have the generic length scale vertical turbulence model?
             if( have_option("/material_phase[0]/subgridscale_parameterisations/GLS/option")) then
                if( (trim(field_name_list(it))=="GLSTurbulentKineticEnergy")) then
                    call gls_tke(state(1))
                else if( (trim(field_name_list(it))=="GLSGenericSecondQuantity")) then
                    call gls_psi(state(1))
                end if
             end if

             ! Calculate the meltrate
             if(have_option("/ocean_forcing/iceshelf_meltrate/Holland08/") ) then
                if( (trim(field_name_list(it))=="MeltRate")) then
                   call melt_surf_calc(state(1))
                endif
             end if

             call get_option(trim(field_optionpath_list(it))//&
                  '/prognostic/equation[0]/name', &
                  option_buffer, default="UnknownEquationType")
             select case(trim(option_buffer))
             case ( "AdvectionDiffusion", "ConservationOfMass", "ReducedConservationOfMass", "InternalEnergy", "HeatTransfer", "KEpsilon" )
                use_advdif=.true.
             case default
                use_advdif=.false.
             end select

             IF(use_advdif)THEN

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
                        state=state, istate=field_state_list(it), global_it=its)

                else if(have_option(trim(field_optionpath_list(it)) // &
                     & "/prognostic/spatial_discretisation/continuous_galerkin")) then

                   call solve_field_equation_cg(field_name_list(it), state, field_state_list(it), dt)
                else

                   ewrite(2, *) "Not solving scalar field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name) //" in an advdif-like subroutine."

                end if ! End of dg/cv/cg choice.

                ! ENDOF IF((TELEDI(IT).EQ.1).AND.D3) THEN ELSE...
             ENDIF

             ewrite(1, *) "Finished field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name)
          end do field_loop

          ! Sort out the dregs of GLS after the solve on Psi (GenericSecondQuantity) has finished
          if( have_option("/material_phase[0]/subgridscale_parameterisations/GLS/option")) then
            call gls_diffusivity(state(1))
          end if
          
          !BC for ice melt
          if (have_option('/ocean_forcing/iceshelf_meltrate/Holland08/calculate_boundaries')) then
            call melt_bc(state(1))
          endif
          
          if(option_count("/material_phase/scalar_field/prognostic/spatial_discretisation/coupled_cv")>0) then
             call coupled_cv_field_eqn(state, global_it=its)
          end if
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          !
          ! Assemble and solve N.S equations.
          !

          do i=1, option_count("/material_phase")
            option_path="/material_phase["//int2str(i-1)//"]/scalar_field::FoamVelocityPotential"
            if( have_option(trim(option_path)//"/prognostic")) then
               call calculate_potential_flow(state(i), phi=foam_velocity_potential)
               call calculate_foam_velocity(state(i), foamvel=foamvel)
               ! avoid outflow bc's for velocity being zero after adapts
               call set_boundary_conditions_values(state, shift_time=.true.)
            end if
          end do 

          ! This is where the non-legacy momentum stuff happens
          ! a loop over state (hence over phases) is incorporated into this subroutine call
          ! hence this lives outside the phase_loop

          if(use_sub_state()) then
             call update_subdomain_fields(state,sub_state)
             call solve_momentum(sub_state,at_first_timestep=((timestep==1).and.(its==1)),timestep=timestep)
             call sub_state_remap_to_full_mesh(state, sub_state)
          else
             call solve_momentum(state,at_first_timestep=((timestep==1).and.(its==1)),timestep=timestep)
          end if

          ! Apply minimum weight condition on weights - population balance
          call dqmom_apply_min_weight(state)

          ! calculate abscissa in the population balance equation
          ! this must be done at the end of each non-linear iteration
          call dqmom_calculate_abscissa(state)
          do i = 1, size(state)
             call dqmom_calculate_moments(state(i))
             call dqmom_calculate_statistics(state(i))
          end do

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

       ! Calculate prognostic sediment deposit fields
       call calculate_sediment_flux(state(1))

       ! Reset the number of nonlinear iterations in case it was overwritten by nonlinear_iterations_adapt
       call get_option('/timestepping/nonlinear_iterations',nonlinear_iterations,&
         & default=1)

       if(have_option("/timestepping/nonlinear_iterations/terminate_if_not_converged")) then
          if(its >= nonlinear_iterations .and. change >= abs(nonlinear_iteration_tolerance)) then
             ewrite(0, *) "Nonlinear iteration tolerance not reached - termininating"
             exit timestep_loop
          end if
       end if

       if (have_option('/mesh_adaptivity/mesh_movement')) then
          ! During the timestep Coordinate is evaluated at n+theta, i.e.
          ! (1-theta)*OldCoordinate+theta*IteratedCoordinate. For writing
          ! the diagnostics we use the end-of-timestep n+1 coordinate however,
          ! so that we can check conservation properties.
          ! Using state(1) should be safe as they are aliased across all states.
          call set_vector_field_in_state(state(1), "Coordinate", "IteratedCoordinate")
          call IncrementEventCounter(EVENT_MESH_MOVEMENT)

          call sync_detector_coordinates(state(1))
       end if

       current_time=current_time+DT
!       ! Calculate the meltrate
!       if(have_option("/ocean_forcing/iceshelf_meltrate/Holland08/") ) then
!          call melt_surf_calc(state(1))
!       end if
       ! calculate and write diagnostics before the timestep gets changed
       call calculate_diagnostic_variables(State, exclude_nonrecalculated=.true.)
       call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)
          
       ! Call the modern and significantly less satanic version of study
       call write_diagnostics(state, current_time, dt, timestep)
       ! Work out the domain volume by integrating the water depth function over the surface if using wetting and drying
       if (have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")) then
          ewrite(1, *) "Domain volume (\int_{fs} (\eta.-b)n.n_z)): ", calculate_volume_by_surface_integral(state(1))
       end if 


       if(have_option("/timestepping/adaptive_timestep")) call calc_cflnumber_field_based_dt(state, dt)

       ! Update the options dictionary for the new timestep and current_time.
       call set_option("/timestepping/timestep", dt)
       call set_option("/timestepping/current_time", current_time)

       ! if strong bc or weak that overwrite then enforce the bc on the fields
       ! (should only do something for weak bcs with that options switched on)
       call set_dirichlet_consistent(state)

       if(have_option("/timestepping/steady_state")) then

          call test_and_write_steady_state(state, change)
          if(change<steady_state_tolerance) then
             ewrite(0,*)  "* Steady state has been attained, exiting the timestep loop"
             exit timestep_loop
          end if

       end if

       call toc(TICTOC_ID_TIMESTEP)
       call tictoc_report(2, TICTOC_ID_TIMESTEP)
       call tictoc_clear(TICTOC_ID_TIMESTEP)

       if(simulation_completed(current_time)) exit timestep_loop

       ! ******************
       ! *** Mesh adapt ***
       ! ******************

       if(have_option("/mesh_adaptivity/hr_adaptivity")) then

          if(do_adapt_mesh(current_time, timestep)) then

             call pre_adapt_tasks(sub_state)

             call qmesh(state, metric_tensor)
             if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
             call run_diagnostics(state)

             call adapt_state(state, metric_tensor)

             call update_state_post_adapt(state, metric_tensor, dt, sub_state, nonlinear_iterations, nonlinear_iterations_adapt)

             if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
             call run_diagnostics(state)
 
          end if
       else if(have_option("/mesh_adaptivity/prescribed_adaptivity")) then
          if(do_adapt_state_prescribed(current_time)) then

             call pre_adapt_tasks(sub_state)

             if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
             call run_diagnostics(state)

             call adapt_state_prescribed(state, current_time)
             call update_state_post_adapt(state, metric_tensor, dt, sub_state, nonlinear_iterations, nonlinear_iterations_adapt)

             if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, current_time, dt, timestep, not_to_move_det_yet=.true.)
             call run_diagnostics(state)

          end if

       not_to_move_det_yet=.false.

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

    ! deallocate the array of all detector lists
    call deallocate_detector_list_array()

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

    ! Deallocate sub_state:
    if(use_sub_state()) then
       do i = 1, size(sub_state)
          call deallocate(sub_state(i))
       end do
    end if

    ! deallocate the pointer to the array of states and sub-state:
    deallocate(state)
    if(use_sub_state()) deallocate(sub_state)

    ! Clean up registered diagnostics
    call destroy_registered_diagnostics 

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

  subroutine pre_adapt_tasks(sub_state)

    type(state_type), dimension(:), pointer :: sub_state

    integer :: ss

    ! GLS - we need to deallocate all module-level fields or the memory
    ! management system complains
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/")) then
        call gls_cleanup() ! deallocate everything
    end if

    ! deallocate sub-state
    if(use_sub_state()) then
      do ss = 1, size(sub_state)
        call deallocate(sub_state(ss))
      end do
      deallocate(sub_state)
    end if

  end subroutine pre_adapt_tasks

  subroutine update_state_post_adapt(state, metric_tensor, dt, sub_state, nonlinear_iterations, nonlinear_iterations_adapt)
    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(out) :: metric_tensor
    real, intent(inout) :: dt
    integer, intent(inout) :: nonlinear_iterations, nonlinear_iterations_adapt
    type(state_type), dimension(:), pointer :: sub_state

    ! Overwrite the number of nonlinear iterations if the option is switched on
    if(have_option("/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt")) then
      call get_option('/timestepping/nonlinear_iterations/nonlinear_iterations_at_adapt',nonlinear_iterations_adapt)
      nonlinear_iterations = nonlinear_iterations_adapt
    end if

    ! The adaptivity metric
    if(have_option("/mesh_adaptivity/hr_adaptivity")) then
      call allocate(metric_tensor, extract_mesh(state(1), topology_mesh_name), "ErrorMetric")
    end if

    ! Auxilliary fields.
    call allocate_and_insert_auxilliary_fields(state)
    call copy_to_stored_values(state,"Old")
    call copy_to_stored_values(state,"Iterated")
    call relax_to_nonlinear(state)

    ! Repopulate substate
    if(use_sub_state()) then
       call populate_sub_state(state,sub_state)
    end if

    ! Discrete properties
    call enforce_discrete_properties(state)

    if (have_option("/mesh_adaptivity/hr_adaptivity/adaptive_timestep_at_adapt")) then
        if (have_option("/timestepping/adaptive_timestep/minimum_timestep")) then
            call get_option("/timestepping/adaptive_timestep/minimum_timestep", dt)
            call set_option("/timestepping/timestep", dt)
        else
            ewrite(-1,*) "Warning: you have adaptive timestep adjustment after &&
                          && adapt, but have not set a minimum timestep"
        end if
    else
        ! Timestep adapt
        if(have_option("/timestepping/adaptive_timestep")) then
            call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
            call set_option("/timestepping/timestep", dt)
        end if
    end if

    ! Ocean boundaries
    if (has_scalar_field(state(1), "DistanceToTop")) then
      if (.not. have_option('/geometry/ocean_boundaries')) then
         ! Leaving this an FLAbort as there's a check up above when initilising.
         ! If we get this error after adapting, then something went very wrong!
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
