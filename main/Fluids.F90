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

  use AllSorts
  use AuxilaryOptions
  use MeshDiagnostics
  use oceansurfaceforcing, only : initialise_ocean_surface_forcing => initialise
  use signal_vars
! NOTE: VARIABLES NOT DEFINED HERE MAY BE FOUND IN Global_Parameters.F90
  use global_parameters
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
  use momentum_source
  use qmesh_module
  use checkpoint
  use write_state_module
  use traffic  
  use synthetic_bc
  use goals
  use adaptive_timestepping
  use conformity_measurement
  use EletemFeinte_module
  ! Use Solid-fluid coupling and ALE - Julian- 18-09-06
  use ale_module
  use adjacency_lists
  use solidity
  use multimaterial_module
  use solid_assembly, only: steval
  use parallel_tools
  use SolidConfiguration
  use MeshMovement
  use legacy_boundary_conditions
  use redfil_module
  use solid_update
  use write_triangle
  use biology
  use momentum_equation
  use timeloop_utilities
  use free_surface_module
  use legacy_field_lists
  use boundary_conditions
  use porous_media
  use spontaneous_potentials, only: calculate_electrical_potential
  use discrete_properties_module
  use gls
  use halos
  use memory_diagnostics
  
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

    INTEGER, PARAMETER ::MXNTSO=31
    ! MXNTSO=max no of field equations.

    ! ********************************
    ! THE VARIABLES USED IN REDFIL ...
    ! For the INTEGERS ...
    INTEGER, SAVE :: &
         & NPHASE,NTSOL,  &
         & NLOC  ,NGI   ,MLOC,&
         & ITINOI,&
         & NCOLOP,&
         & SNLOC, SNGI  ,&
         & OPTSOU,ISPHER2,&
         & MISNIT,&
         & VERSIO,&
         & NPRESS,NPROPT,&
         & RADISO,&
                                ! Phase momentum equations...
         & MULPA(MXNPHA), CGSOLQ(MXNPHA),GMRESQ(MXNPHA),&
         & MIXMAS(MXNPHA),UZAWA(MXNPHA),&
         & POISON(MXNPHA),PROJEC(MXNPHA),&
         & EQNSTA(MXNPHA),PREOPT(MXNPHA),&
                                ! Field equation...
         & DISOTT(MXNTSO),TPHASE(MXNTSO),&
         & CGSOLT(MXNTSO),&
         & GMREST(MXNTSO),&
         & IDENT(MXNTSO),&
         & TELEDI(MXNTSO),&
         & NSUBTLOC(MXNTSO),&
         & NDISOT(MXNTSO)


    ! This is for LOGICALS...
    LOGICAL, SAVE :: &
         & NAV   ,MVMESH,&
         & CMCHAN,&
         & GETTAN,ROTAT, DSPH,  BHOUT,&
         & RAD, &
         & COGRAX,COGRAY,COGRAZ,&
         & ADMESH,&
                                ! Phase momentum equations...
         & MAKSYM(MXNPHA),&
         & CONVIS(MXNPHA),&
         & CHADEN(MXNPHA),&
         & SLUMP(MXNPHA),&
         & COMPRE(MXNPHA),&
                                ! Field equation...
         & BOUSIN(MXNTSO),TLUMP(MXNTSO),&
         & SUFTEM(MXNTSO)

    ! For the REALS ...
    REAL, SAVE :: &
         & LTIME, &
         & ITHETA,ITIERR,STEDER,&
         & R0,D0,&
         & ALFST2,SPRESS,&
                                ! Phase momentum equations...
         & BETA(MXNPHA),&
         & TEMINI(MXNPHA),DENINI(MXNPHA),&
         & BSOUX(MXNPHA),BSOUY(MXNPHA),BSOUZ(MXNPHA),&
         & GAMDE2(MXNPHA),GAMDE3(MXNPHA),&
                                ! Field equation...
         & TTHETA(MXNTSO),TBETA(MXNTSO)

    ! THE FOLLOWING ARE USED IN REDSCA*********************
    INTEGER, SAVE :: NONODS,XNONOD,TOTELE,FREDOP,&
         &  NOBCU(MXNPHA),NOBCV(MXNPHA),NOBCW(MXNPHA),  &
         &  NOBCT(MXNTSO),NNODP,&
         &  NNODPP, NDPSET,&
         &  NNODRO, STOTEL

    !local memory for old style boundary conditions for NTSOL tracers
    type(real_vector), pointer, dimension(:), save :: bct1_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bct2_mem => null()
    !local memory for old style boundary conditions for velocities
    type(real_vector), pointer, dimension(:), save :: bcu1_mem => null()
    type(real_vector), pointer, dimension(:), save :: bcv1_mem => null()
    type(real_vector), pointer, dimension(:), save :: bcw1_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bcu2_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bcv2_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bcw2_mem => null()

    ! THE REMAINING VARIABLES DEFINED IN THIS SUB***********
    INTEGER, SAVE :: PROCNO


    ! **** ENDOF VARIABLES USED IN REDFIL
    ! ***********************************

    !     System state wrapper.
    type(state_type), dimension(:), pointer, save :: state => null()
    type(tensor_field) :: metric_tensor
    !     Dump index
    integer, save :: dump_no = 0
    !     Temporary buffer for any string options which may be required.
    character(len=666) :: option_buffer
    !     Status variable for option retrieval.
    integer :: option_stat

    REAL, SAVE :: CHANGE,CHAOLD


    LOGICAL, SAVE :: REMESH

    INTEGER, SAVE :: IT,ITP,ITS
    INTEGER, SAVE :: I

    ! FOR SOURCES...
    INTEGER, SAVE :: NSOUPT
    INTEGER, PARAMETER :: MXNSOU=3000
    INTEGER, SAVE :: FIESOU(MXNSOU)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !     STUFF for MEsh movement, and Solid-fluid-coupling.  ------ jem

    !     Ale mesh movement - Julian 05-02-07
    LOGICAL, save:: USE_ALE
    INTEGER, save:: fs

    !     Solid-fluid coupling - Julian 18-09-06
    !New options
    INTEGER, save :: ss,ph
    LOGICAL, SAVE :: have_solids

    ! Pointers for scalars and velocity fields
    type(vector_field) :: Velocity
    type(scalar_field), pointer :: sfield
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !     backward compatibility with new option structure - crgw 21/12/07
    logical::use_advdif=.true.  ! decide whether we enter advdif or not
    logical::use_electrical_potential_solver=.false.  !jhs 05/05/2009
    character(len=OPTION_PATH_LEN) :: tmpstring
    integer :: tmpstat

    !***********************************************************************
    !     Begin parameters for SOLIDITY            --Added by gsc 2006-03-13

    !     These options are set in "solidity_options.inp"
    !                                 (see GET_SOLIDTY_OPTIONS in Solidity.F90)
    INTEGER, save :: SOLIDS            !Solid or fluid rheology

    !     crgw - compressible schemes
    INTEGER, save :: MKCOMP
    !     mkcomp <= 0... incompressible
    !     mkcomp = 3... compressible scheme using modified gradient operator

    INTEGER, SAVE :: adapt_count

    logical, dimension(MXNPHA), save :: phase_uses_new_code_path=.false.
    ! whether any of the phases use the old code path (i.e. phase_uses_new_code_path==.false.)
    ! *or* any of the scalar fields use advdif
    logical :: uses_old_code_path=.false.

    ! Current simulation timestep
    integer, save :: timestep

    ! a counter for the total number of global non linear iterations including repeated timesteps 
    integer, save :: total_its_count

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

    DT=0.0
    REMESH=.FALSE.

    NPRESS=1
    NPROPT=1

    PROCNO = GetProcNo()

    ! Read state from .flml file
    call populate_state(state)
    
    ! Check the diagnostic field dependencies for circular dependencies
    call check_diagnostic_dependencies(state)
    
    if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) then
       if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, acctim, dt)
       call adapt_state_first_timestep(state)
       
       call allocate_and_insert_auxilliary_fields(state)
       call copy_to_stored_values(state,"Old")
       call copy_to_stored_values(state, "Iterated")
       call relax_to_nonlinear(state)
       
       if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, acctim, dt)
    else
       call allocate_and_insert_auxilliary_fields(state)
       call copy_to_stored_values(state,"Old")
       call copy_to_stored_values(state, "Iterated")
       call relax_to_nonlinear(state)
    end if

    call enforce_discrete_properties(state)
    call get_option("/timestepping/timestep", dt)
    if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
       call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
       call set_option("/timestepping/timestep", dt)
    end if

    call get_option('/timestepping/nonlinear_iterations',ITINOI,&
         & default=1)
     MVMESH=have_option("/mesh_adaptivity/mesh_movement")
     call get_option("/timestepping/current_time", ACCTIM)
     call get_option("/timestepping/finish_time", LTIME)

    ! Calculate the number of scalar fields to solve for and their correct
    ! solve order taking into account dependencies.
    call get_ntsol(ntsol)
    call get_nphase(nphase)

    call initialise_field_lists_from_options(state, ntsol)
    call initialise_state_phase_lists_from_options()

    ! We get back here after a mesh adapt with remesh = .false.
99771 CONTINUE

    call compute_uses_old_code_path(uses_old_code_path)
    if(uses_old_code_path) FLExit("The old code path is dead.")
    
    


    !     populate state or adapt_state_new_options has created a
    !     'populated' state as if read from disk
    !     so we continue as from the start:
    !     redfil() reads from state


    CALL REDFIL(&
                                ! THE VARIABLES DEFINED IN COMSCA(1ST LINES IS USED IN)
                                ! *****************************************************
         & MXNTSO,MXNPHA,&
                                ! For the INTEGERS ...
         & NPHASE,NTSOL,  &
         & NLOC,NGI,MLOC,&
         & ITINOI,&
         & NCOLOP,&
         & SNLOC, SNGI,VERSIO, &
         & NPRESS,NPROPT,&
         & RADISO,&
         & GEOBAL,OPTSOU,ISPHER2,&
         & MISNIT,&
                                ! Phase momentum equations...
         & DISOPT,MULPA, CGSOLQ,GMRESQ,&
         & MIXMAS,UZAWA,POISON,PROJEC,&
         & EQNSTA,PREOPT,NDISOP,NSUBVLOC,NSUBNVLOC,&
                                ! Field equation...
         & DISOTT,TPHASE,CGSOLT,GMREST,&
         & IDENT,&
         & TELEDI,&
         & NSUBTLOC,&
         & NDISOT,&
                                ! This is for LOGICALS...
         & D3,DCYL  ,NAV   ,MVMESH,&
         & CMCHAN,GETTAN,&
         & ROTAT,DSPH,BHOUT,RAD,COGRAX,COGRAY,COGRAZ,ADMESH,&
                                ! Phase momentum equations...
         & LUMP  ,MAKSYM,CONVIS,&
         & CHADEN, &
         & ABSLUM,SLUMP,&
         & COMPRE,&
                                ! Field equation...
         & BOUSIN,TLUMP,&
         & SUFTEM,&
                                ! For the REALS ...
         & ACCTIM,LTIME ,DT,&
         & ITHETA,ITIERR,STEDER,&
         & R0,D0,GRAVTY, &
         & ALFST2,SPRESS,&
                                ! Phase momentum equations...
         & THETA ,BETA  , &
         &DENGAM,&
         & TEMINI,DENINI,&
         & BSOUX,BSOUY,BSOUZ,&
         & GAMDE2,GAMDE3,&
                                ! Field equation...
         & TTHETA,TBETA,&
                                ! THE FOLLOWING ARE USED IN REDSCA*********************
         & NONODS,XNONOD,TOTELE,FREDOP,&
         & NOBCU,NOBCV,NOBCW,  NOBCT, NNODP,&
         & NNODPP, NDPSET,&
         & NNODRO, STOTEL, &
         & NSOUPT,MXNSOU,FIESOU,&
                                ! THE REMAINING VARIABLES DEFINED IN THIS SUB***********
         & procno, halo_tag, halo_tag_p,&
         & bcu1_mem, bcu2_mem, &
         & bcv1_mem, bcv2_mem, &
         & bcw1_mem, bcw2_mem, &
         & bct1_mem, bct2_mem, &
         &              size(state),state, uses_old_code_path)

    ! We jump to here before a mesh adapt (but after metric assembly), with
    ! remesh = .true.
9977 CONTINUE

    call run_diagnostics(state)

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !     Solidity modification -----------------------
    !        Get the last remaining old style solidity variables
    CALL GET_SOLIDITY_OPTIONS(SOLIDS,MKCOMP)

    !     -------------------------------- Added by gsc 2006-03-29

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     Initialise solid-fluid coupling, and ALE ----------- -Julian 17-07-2008
    have_solids=.false.
    use_ale=.false.

    !Read the amount of SolidConcentration fields
    !and save their IT numbers.
    if(have_option(trim('/imported_solids'))) then
       ss=0
       phaseloop: do ph = 1, size(state)
          write(tmpstring, '(a,i0,a)') "/material_phase[",ph-1,"]"
          if(have_option(trim(tmpstring)//"/scalar_field::SolidConcentration")) then
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
          write(tmpstring, '(a,i0,a)') "/material_phase[",ph-1,"]"
          if(have_option(trim(tmpstring)//"/scalar_field::MaterialVolumeFraction/prognostic")) then
             fs = ph
          end if
       end do phaseloop1
       if (fs==-1) then
          FLAbort('No prognostic MaterialVolumeFraction was defined')
       end if
    end if
    !     end initialise solid-fluid coupling, and ALE  -Julian 17-07-2008
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


    if(remesh) then
       ! ******************
       ! *** Mesh adapt ***
       ! ******************
       if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, acctim, dt)
       
       call adapt_state(state, metric_tensor)

       call allocate_and_insert_auxilliary_fields(state)
       call enforce_discrete_properties(state)
       call initialise_ocean_surface_forcing(state(1))

       ! reinitialise the auxilliary fields before the adaptive timestep check
       IF(ITINOI.GT.1) THEN
          call copy_to_stored_values(state,"Old")
          call copy_to_stored_values(state, "Iterated")
          call relax_to_nonlinear(state)
       ENDIF

       if(have_option("/timestepping/adaptive_timestep")) then
          call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
          call set_option("/timestepping/timestep", dt)
       end if
       
       if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, acctim, dt)

       remesh = .false.

       GOTO 99771
    ELSE
       if (admesh) then
          call allocate(metric_tensor, extract_mesh(state(1), "CoordinateMesh"), "ErrorMetric")
       end if

    ENDif
    ! REMESHING ************************************************
    ! **********************************************************

    call compute_phase_uses_new_code_path(state, phase_uses_new_code_path)
    if (any(.not.phase_uses_new_code_path)) FLExit("The old code path is dead")

    !     Determine the output format.
    call get_option('/io/dump_format', option_buffer, option_stat)
    if(option_stat /= 0) then
       FLAbort("You must specify a dump format and it must be vtk")
    else if(trim(option_buffer) /= "vtk") then
       FLAbort("You must specify a dump format and it must be vtk")
    end if

    !        Initialisation of distance to top and bottom field
    !        Currently only needed for free surface

    if (has_scalar_field(state(1), "DistanceToTop")) then
       if (.not. have_option('/geometry/ocean_boundaries')) then
          FLAbort("There are no top and bottom boundary markers.")
       end if
       call CalculateTopBottomDistance(state(1), cograx)
    end if

    !**********************Reduced model************************
    ! these are used only when the inital condition is inversed
    if(have_option("/model/fluids/reduced/initial")) then
       FLExit("Calling the reduced model from new options doesn't work yet")
    endif
    !******************End reduced model***********************

    ! if adaptivity hasn't just happened then initialise the multimaterial
    ! fields on the first timestep
    if(timestep == 0) then
       call initialise_diagnostic_material_properties(state)
    end if

    call calculate_diagnostic_variables(State)
    call calculate_diagnostic_variables_new(state)
    ! This is mostly to ensure that the photosynthetic radiation
    ! has a non-zero value before the first adapt. 
    if (have_option("/ocean_biology")) then
       call calculate_biology_terms(state(1))
    end if

    ! We really really only want first timestep output on the first timestep
    ! (not when leaping around with adaptivity gotos)
    if(timestep == 0) then
      ! Dump at start
      if(do_checkpoint_simulation(dump_no)) call checkpoint_simulation(state, cp_no = dump_no)
      ! Dump at start  
      if( &
             ! if this is not a zero timestep simulation (otherwise, there would
             ! be two identical dump files)
           & acctim < ltime &
             ! unless explicitly disabled
           & .and. .not. have_option("/io/disable_dump_at_start") &
           & ) then
         call write_state(dump_no, state)
      end if

      call initialise_diagnostics(filename, state)
      call initialise_convergence(filename, state)
      call initialise_advection_convergence(state)
      if(have_option("/io/stat/output_at_start")) call write_diagnostics(state, acctim, dt)
    end if

    do i=1, size(state)
       velocity = extract_vector_field(state(i), "Velocity")
       if (has_boundary_condition(velocity, "free_surface") .and. &
            & have_option('/mesh_adaptivity/mesh_movement')) then
          ewrite(1,*) "Going into move_free_surface_nodes to compute surface node coordinates from initial condition"
          call move_free_surface_nodes(state(i))
       end if
    end do

    ! ******************************
    ! *** Start of timestep loop ***
    ! ******************************
    timestep_loop: do
       timestep = timestep + 1

       ewrite(1, *) "********************"
       ewrite(1, *) "*** NEW TIMESTEP ***"
       ewrite(1, *) "********************"
       ewrite(1, *) "Current simulation time (acctim): ", acctim
       ewrite(1, *) "Timestep number: ", timestep
       ewrite(1, *) "Timestep size (dt): ", dt
       if(.not. allfequals(dt)) then
          ewrite(-1, *) "Timestep size (dt): ", dt
          FLAbort("The timestep is not global across all processes!")
       end if

       !     solidity -------------------------------------------
       ITS = 0
       !     ----------------------------------------- added by crgw 23/06/06
       total_its_count = 0

       !******************REDUCED MODEL--ffx****************
       if(have_option("/model/fluids/reduced")) then
          FLExit("Calling the reduced model from new options doesn't work yet")
       endif
       !******************END REDUCED MODEL--ffx****************

       IF(MVMESH.AND.(SOLIDS.EQ.0)) THEN
          ! Make oldcoordinate a copy of coordinate.
          call set_vector_field_in_state(state(1), "OldCoordinate", &
               "Coordinate")
       ENDIF

       if(simulation_completed(acctim, timestep)) exit timestep_loop

       if( &
                                ! Do not dump at the start of the simulation (this is handled by write_state call earlier)
            & acctim > simulation_start_time &
                                ! Do not dump at the end of the simulation (this is handled by later write_state call)
            & .and. acctim < ltime &
                                ! Test write_state conditions
            & .and. do_write_state(acctim, timestep) &
            & ) then

          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ! Regular during run state dump.
          !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          ! Intermediate dumps
          if(do_checkpoint_simulation(dump_no)) then
             call checkpoint_simulation(state, cp_no = dump_no)
          end if
          call write_state(dump_no, state)
       ENDIF

       ewrite(2,*)'STEDER,ITINOI:',STEDER,ITINOI

       IF((STEDER.GT.0).OR.(ITINOI.GT.1)) THEN
          call copy_to_stored_values(state,"Old")
       ENDIF

       ! this may already have been done in populate_state, but now
       ! we evaluate at the correct "shifted" time level:
       call set_boundary_conditions_values(state, shift_time=.true.)

       !579    CONTINUE  ! Go back to here for repeating time step.

       ! ITINOI=maximum no of iterations within a time step
       ! NB TEMPT is T from previous iteration.
       nonlinear_iteration_loop: do  ITS=1,ITINOI

          ewrite(1,*)'###################'
          ewrite(1,*)'Start of another nonlinear iteration; ITS,ITINOI=',ITS,ITINOI
          ewrite(1,*)'###################'

          total_its_count = total_its_count + 1

          IF(ITINOI.GT.1) THEN
             call copy_to_stored_values(state, "Iterated")
             call relax_to_nonlinear(state)
          ENDIF

          call compute_goals(state)

          IF(ITS.GT.1) THEN   !---CHANGED 22/03/07 CRGW
             call copy_from_stored_values(state, "Old")
          ENDIF ! if(its.gt.1)

          !------------------------------------------------
          ! Addition for calculating drag force ------ jem 05-06-2008
          if (have_option("/imported_solids/calculate_drag_on_surface")) then
             call drag_on_surface(state)
          end if

          !     Addition for reading solids in - jem  02-04-2008
          if(have_solids) call solid_configuration(state(ss:ss),its,itinoi)

          !Explicit ALE ------------   jem 21/07/08         
          if (use_ale) then 
             EWRITE(0,'(A)') '----------------------------------------------'
             EWRITE(0,'(A26,E12.6)') 'Using explicit_ale. time: ',ACCTIM       
             call explicit_ale(state,fs)
          end if
          !end explicit ale ------------  jem 21/07/08

          !-----------------------------------------------------
          ! Call to porous_media module (leading to multiphase flow in porous media)
          ! jhs - 16/01/09
          if (have_option("/porous_media")) then
             call porous_media_advection(state, nonods)
             ! compute spontaneous electrical potentials (myg - 28/10/09)
             do i=1,size(state)
                tmpstring = '/material_phase['//int2str(i-1)//']/electrical_properties/coupling_coefficients/'
                if (have_option(trim(tmpstring)//'scalar_field::Electrokinetic').or.&
                    have_option(trim(tmpstring)//'scalar_field::Thermoelectric').or.&
                    have_option(trim(tmpstring)//'scalar_field::Electrochemical')) then
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

             ITP=STATE2PHASE_INDEX(FIELD_STATE_LIST(IT))

             ! do we have the generic length scale vertical turbulence model - search in first material_phase only at the moment
             if( have_option("/material_phase[0]/subgridscale_parameterisations/GLS/option") ) then
                if( (trim(field_name_list(it))=="GLSTurbulentKineticEnergy")) then
                   !               if( (trim(field_name_list(it))=="GLSTurbulentKineticEnergy").or.(trim(field_name_list(it))=="GLSGenericSecondQuantity")) then
                   call gls_vertical_turbulence_model(state(1))            
                endif
             end if

             call get_option(trim(field_optionpath_list(it))//&
                  '/prognostic/equation[0]/name', &
                  tmpstring, stat=tmpstat)
             if (tmpstat==0) then
                select case(trim(tmpstring))
                case ( "AdvectionDiffusion", "ConservationOfMass", "ReducedConservationOfMass", "InternalEnergy" )
                   use_advdif=.true.
                case ( "ElectricalPotential" )
                   use_advdif=.false.
                   use_electrical_potential_solver=.true.
                case default
                   use_advdif=.false.
                end select
             else
                use_advdif=.false.
             end if

             IF(((ITP.GE.1).OR.(NPHASE.LE.1)).AND.(use_advdif))THEN

                ITP=MAX(1,ITP)

                if(have_option("/traffic_model/scalar_field::TrafficTracerTemplate"))then
                   call traffic_tracer(it,state(1),timestep)
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
                        & state=state(state2phase_index(field_state_list(it))))

                ELSEIF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/finite_volume")) then

                   ! Solve the FV form of the equations.
                   call solve_advection_diffusion_fv(field_name=field_name_list(it), &
                        & state=state(state2phase_index(field_state_list(it))))

                ELSEIF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/control_volumes")) then

                   ! Solve the pure control volume form of the equations
                   call solve_field_eqn_cv(field_name=trim(field_name_list(it)), &
                        state=state(field_state_list(it):field_state_list(it)), &
                        global_it=its)

                else if(have_option(trim(field_optionpath_list(it)) // &
                     & "/prognostic/spatial_discretisation/continuous_galerkin")) then

                   call solve_field_equation_cg(field_name_list(it), state(state2phase_index(field_state_list(it))), dt)
                else

                   ewrite(2, *) "Not solving scalar field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name) //" in an advdif-like subroutine."

                end if ! End of dg/cv/cg choice.

                ! ENDOF IF((TELEDI(IT).EQ.1).AND.D3) THEN ELSE...
             ENDIF

             ewrite(1, *) "Finished field " // trim(field_name_list(it)) // " in state " // trim(state(field_state_list(it))%name)
          end do field_loop

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
             call porous_media_momentum(state,nonods)
          end if
          
          if (have_option("/traffic_model")) then
             call traffic_source(state(1),timestep)
             call traffic_density_update(state(1))
          endif
          
          ! This is where the non-legacy momentum stuff happens
          ! a loop over state (hence over phases) is incorporated into this subroutine call
          ! hence this lives outside the phase_loop
          call momentum_loop(state, at_first_timestep=((timestep==1).and.(its==1)))
          
          ! Add in externally defined momentum source
          ! term. Don't even think about doing this with multiple
          ! material phases !
          call simple_source_python(state(1))             
             
          if(itinoi > 1) then
             ! Check for convergence between non linear iteration loops
             call test_and_write_convergence(state, acctim + dt, dt, its, change)
             if(its == 1) chaold = change             

             if (have_option("/timestepping/nonlinear_iterations/&
                              &tolerance")) then
               ewrite(2, *) "Nonlinear iteration change = ", change
               ewrite(2, *) "Nonlinear iteration tolerance = ", itierr

               if(change < abs(itierr)) then
                  ewrite(1, *) "Nonlinear iteration tolerance has been reached"
                  ewrite(1, "(a,i0,a)") "Exiting nonlinear iteration loop after ", its, " iterations"
                  exit nonlinear_iteration_loop
               endif
             end if
          end if

       end do nonlinear_iteration_loop

       if(have_option("/timestepping/nonlinear_iterations/terminate_if_not_converged")) then
          if(its >= itinoi .and. change >= abs(itierr)) then
             ewrite(0, *) "Nonlinear iteration tolerance not reached - termininating"
             exit timestep_loop
          end if
       end if

       if(have_option(trim('/mesh_adaptivity/mesh_movement/vertical_ale'))) then
          ewrite(1,*) 'Entering vertical_ale routine'
          !move the mesh and calculate the grid velocity
          call movemeshy(state(1))
       endif

       !     -------------------------------start of add by crgw 14/03/06
       ! solidity - move the nodes using the final grid velocities
       ! (not intermediate steps) if
       !  we are using a lagrangian mesh but itinoi>1 for material model
       !  purposes
       IF(MVMESH.AND.SOLIDS.GT.0) THEN
          call solid_coordinate_update(state(1))
       ENDIF
       !     --------------------------------------end of add by crgw 14/03/06

       ACCTIM=ACCTIM+DT

       ! calculate and write diagnostics before the timestep gets changed
       call calculate_diagnostic_variables(State, exclude_nonrecalculated=.true.)
       call calculate_diagnostic_variables_new(state, exclude_nonrecalculated = .true.)

       ! Call the modern and significantly less satanic version of study
       call write_diagnostics(state, acctim, dt)

       ! Note that other adaptive timestepping routines are deliberately disabled:
       ! - adaptive_timestepping_check_options ensures that autacc must be
       !   positive, thereby preventing non-linear iteration based adaptive
       !   timestepping
       ! - tagain is currently not available
       if(have_option("/timestepping/adaptive_timestep")) call calc_cflnumber_field_based_dt(state, dt)

       ! Update the options dictionary for the new timestep and current_time.
       call set_option("/timestepping/timestep", dt)
       call set_option("/timestepping/current_time", acctim)

       call set_prescribed_field_values(state, exclude_interpolated=.true., &
            exclude_nonreprescribed=.true.)
       call enforce_discrete_properties(state, only_prescribed=.true., &
            exclude_interpolated=.true., &
            exclude_nonreprescribed=.true.)
       ! if strong bc or weak that overwrite then enforce the bc on the fields
       ! (should only do something for weak bcs with that options switched on)
       call set_dirichlet_consistent(state)

       if(have_option("/timestepping/steady_state")) then

          call test_steady_state(state, change)
          if(change<steder) then
             ewrite(0,*)  "* Steady state has been attained, exiting the timestep loop"
             exit timestep_loop
          end if

       end if

       if(simulation_completed(acctim)) exit timestep_loop

       if(have_option("/mesh_adaptivity/hr_adaptivity")) then
          if( &
                                ! Test qmesh conditions
               & do_adapt_mesh(acctim, timestep) &
               & ) then

             call zero(metric_tensor)
             call qmesh(state, metric_tensor)

             ! Let's go adapt the mesh
             remesh = .true.
             goto 9977
          end if
       else if(have_option("/mesh_adaptivity/prescribed_adaptivity")) then
          if(do_adapt_state_prescribed(acctim)) then
             if(have_option("/io/stat/output_before_adapts")) call write_diagnostics(state, acctim, dt)
             
             call adapt_state_prescribed(state, acctim)

             call allocate_and_insert_auxilliary_fields(state)
             call enforce_discrete_properties(state)
             call initialise_ocean_surface_forcing(state(1))

             ! reinitialise the auxilliary fields before the adaptive timestep check
             if(itinoi > 1) then
                call copy_to_stored_values(state, "Old")
                call copy_to_stored_values(state, "Iterated")
                call relax_to_nonlinear(state)
             endif

             if(have_option("/timestepping/adaptive_timestep")) then
                call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
                call set_option("/timestepping/timestep", dt)
             end if

             if(have_option("/io/stat/output_after_adapts")) call write_diagnostics(state, acctim, dt)
             
             ! We now have a perfectly valid system state ready for timestepping.
             ! However we need to set up the legacy variables, so let's leap out
             ! of the timestep loop to set them up.
             remesh = .false.
             goto 99771
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
    
  subroutine compute_uses_old_code_path(uses_old_code_path)
  !!< computes whether any of the phases use the old code path (i.e. phase_uses_new_code_path==.false.)
  !!< *or* any of the scalar fields use advdif
    logical, intent(out) :: uses_old_code_path
    
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
    
        uses_old_code_path=.true.
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "** DANGER! LEGACY FUNCTIONALITY WILL BE REMOVED **"
        ewrite(0,*) "** -------------------------------------------- **"
        ewrite(0,*) "** If there is a reason why you need this       **" 
        ewrite(0,*) "** functionality, you MUST email                **"
        ewrite(0,*) "** amcg-users@imperial.ac.uk NOW!!!!!!!!!!!!!!  **"
        ewrite(0,*) "**                                              **"
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "WARNING: You seem to be using legacy_continuous_galerkin or "
        ewrite(0,*) "legacy_discretisation for spatial_discretisation of Velocity."
        ewrite(0,*) "This uses the old code path (navsto) that is going to be deleted"
        ewrite(0,*) "in the near future. You should switch to continuous_galerkin."
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        ewrite(0,*) "**************************************************" 
        return
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
            
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "** DANGER! LEGACY FUNCTIONALITY WILL BE REMOVED **"
          ewrite(0,*) "** -------------------------------------------- **"
          ewrite(0,*) "** If there is a reason why you need this       **" 
          ewrite(0,*) "** functionality, you MUST email                **"
          ewrite(0,*) "** amcg-users@imperial.ac.uk NOW!!!!!!!!!!!!!!  **"
          ewrite(0,*) "**                                              **"
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "WARNING: You seem to be using legacy_continuous_galerkin,"
          ewrite(0,*) "legacy_mixed_cv_cg or legacy_discretisation for the"
          ewrite(0,*) "spatial discretisation of one of your scalar fields."
          ewrite(0,*) "This uses the old code path (advdif) that is going to be deleted"
          ewrite(0,*) "in the near future. You should switch to any of the other options."
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          ewrite(0,*) "**************************************************" 
          
          uses_old_code_path=.true.
          return
          
        end if
         
      end do
    end do
      
    ! Hurray!!
    uses_old_code_path=.false.
    
  end subroutine compute_uses_old_code_path
  
  end module fluids_module

