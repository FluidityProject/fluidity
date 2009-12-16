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
  use AdvectionDiffusion
  use AuxilaryOptions
  use conacc_module
  use CriticalTimeStep
  use MeshDiagnostics
  use oceansurfaceforcing, only : initialise_ocean_surface_forcing => initialise
  use signal_vars
! NOTE: VARIABLES NOT DEFINED HERE MAY BE FOUND IN Global_Parameters.F90
  use global_parameters
  use spud
  use mellor_yamada
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
  use assemble_buoyancy
  use write_triangle
  use biology
  use momentum_equation
  use timeloop_utilities
  use navsto_module
  use free_surface_module
  use spaerr_module
  use legacy_field_lists
  use boundary_conditions
  use porous_media
  use spontaneous_potentials, only: calculate_electrical_potential
  use discrete_properties_module
  use flcomms_module
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

    ! KE_ident and LKE_ident are the idents for the KineticEnergy and TurbulentLengthScalexKineticEnergy fields that are used in the MY 2.5 eqn turbulence model.
    ! MY_vis_ident and MY_diff_ident are the idents for the MY 2.5 eqn turbulence model vertical viscosity and diffusivity.
    integer, parameter ::KE_ident=101,LKE_ident=102
    INTEGER, PARAMETER ::MY_vis_ident=103,MY_diff_ident=104
    INTEGER, PARAMETER ::PSIDEN=-2003
    integer, parameter:: FREES_IDENT=-29

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
         & nlevel, geostrophic_solver_option,&
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
    type(real_vector), pointer, dimension(:), save :: bct1_mem => null(), bct1w_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bct2_mem => null()
    !local memory for old style boundary conditions for velocities
    type(real_vector), pointer, dimension(:), save :: bcu1_mem => null(), bcu1w_mem => null()
    type(real_vector), pointer, dimension(:), save :: bcv1_mem => null(), bcv1w_mem => null()
    type(real_vector), pointer, dimension(:), save :: bcw1_mem => null(), bcw1w_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bcu2_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bcv2_mem => null()
    type(integer_vector), pointer, dimension(:), save :: bcw2_mem => null()

    ! THE REMAINING VARIABLES DEFINED IN THIS SUB***********
    INTEGER, SAVE :: PARA,NPROCS,PROCNO
    INTEGER, SAVE :: NSOGRASUB=0
    LOGICAL, SAVE :: GOTBOY

    INTEGER, SAVE :: GOTTRAF

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
    INTEGER, SAVE :: NCMCB

    LOGICAL, SAVE :: CMCGET,GETC12
    REAL, SAVE :: CHANGE,CHAOLD

    !     --------------------------------------start of add by crgw 31/03/06
    !     solidity - strain for plastic flow:
    REAL, ALLOCATABLE, DIMENSION(:) :: STVPXX, STVPYY, STVPZZ
    REAL, ALLOCATABLE, DIMENSION(:) :: STVPYZ, STVPXZ, STVPXY
    REAL, ALLOCATABLE, DIMENSION(:) :: YIELDF, YIELDFOLD
    !     --------------------------------------end of add by crgw 31/03/06

    ! More pointers(to do with iteration & time stepping)

    REAL, SAVE :: OLDDT
    INTEGER, SAVE :: POISO2(MXNPHA)
    INTEGER, SAVE :: IPHASE

    LOGICAL, SAVE :: REMESH

    !c Pointers used in this memory
    LOGICAL, SAVE ::  NEWMES2
    INTEGER, SAVE ::  ITFREEimp,NOBCTITFREEimp,ITS4
    ! the disott of the free surface field (ident -29) or of balanced pressure (-2003)
    ! used to binary decode into free surface options in geoeli1p()
    integer freesdisott

    INTEGER, SAVE :: IT,IT2,IP,ITP,ITS
    integer it_KE, it_LKE, it_temp, it_MY_vis, it_MY_diff
    INTEGER, SAVE :: I,II,IPN,IPF
    INTEGER, SAVE :: NOPHAS

    ! The pointers for radiation ...
    ! pointers for reals
    INTEGER, SAVE :: NRTDR

    INTEGER, SAVE :: ITKE
    LOGICAL, SAVE :: BOUINI
    LOGICAL, SAVE :: NDWISP
    ! FOR SOURCES...
    INTEGER, SAVE :: NSOUPT
    INTEGER, PARAMETER :: MXNSOU=3000
    INTEGER, SAVE :: FIESOU(MXNSOU)

    REAL, SAVE ::     YOUR_SCFACTH0 = 1.0

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     SALINITY STUFF FOR OCEAN MODELLING - cjc
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    integer, save :: it_salt
    real, save ::  ztop
    !cccccIDENTITY FOR SALTcccccccccc
    integer, parameter :: ident_sal = 42
    logical, save :: gotsal, gotvis
    !cccccccccccccccccccccccccccccccc
    logical, save :: got_top_bottom_distance
    logical, save :: got_top_bottom_markers

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     Stuff for ramped boundary conditions - cjc
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    logical, save :: copyu, copyv, copyw

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    !     STUFF for MEsh movement, and Solid-fluid-coupling.  ------ jem

    !     Ale mesh movement - Julian 05-02-07
    LOGICAL, save:: USE_ALE
    INTEGER, save:: fs

    !     Solid-fluid coupling - Julian 18-09-06
    !New options
    INTEGER, save :: NSOLIDS
    INTEGER, save ::  ITSOLID(30),ss,ph
    LOGICAL, SAVE :: have_solids

    ! Pointers for scalars and velocity fields
    type(scalar_field) :: Tracer
    type(vector_field) :: Velocity,GridVelocity
    type(scalar_field) :: VelocityX, VelocityY, VelocityZ    
    type(scalar_field) :: GridVelocityX, GridVelocityY, GridVelocityZ
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

    LOGICAL, SAVE :: OFFYIELDENV        !determines (written to by steval) whether another cycle is necesarry
    INTEGER, SAVE :: OITINOI            !remembers initial itinoi and defaults back at next timestep
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

    call initial_allocation()

    call initialise_qmesh
    call initialise_write_state

    if (have_option("/geometry/disable_geometric_data_cache")) then
       ewrite(1,*) "Disabling geometric data cache"
       cache_transform_elements=.false.
    end if

    adapt_count = 0

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     Flags to say if we need to copy u,v,w - cjc

    copyu = .true.
    copyv = .true.
    copyw = .true.
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    DT=0.0
    REMESH=.FALSE.

    NPRESS=1
    NPROPT=1

    BOUINI=.TRUE.

    ! Initialise PVM and find PROCNO, PROCNO & NPROCS ***********
    ! AND initialise some other parallel stuff.
    IF(IsParallel()) THEN
       PARA = 1
    ELSE
       PARA = 0
    END IF
    NPROCS = GetNProcs()
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
    if(have_option("/timestepping/adaptive_timestep/at_first_timestep")) then
       call get_option("/timestepping/timestep", dt)
       call calc_cflnumber_field_based_dt(state, dt, force_calculation = .true.)
       call set_option("/timestepping/timestep", dt)
    end if

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

    GOTTRAF=0
    NSOLIDS=0
    !Read the amount of SolidConcentration fields
    !and save their IT numbers.
    if(have_option(trim('/imported_solids'))) then
       ss=0
       phaseloop: do ph = 1, size(state)
          write(tmpstring, '(a,i0,a)') "/material_phase[",ph-1,"]"
          if(have_option(trim(tmpstring)//"/scalar_field::SolidConcentration")) then
             ss = ph
          end if
       end do phaseloop
       if(ss==0) then
          FLAbort("Havent found a material phase containing solid concentration")
       end if
       EWRITE(2,*) 'solid state= ',ss
       NSOLIDS=0
       DO IT=1,NTSOL
          ewrite(2,*) trim(field_name_list(it))
          if (trim(field_name_list(it))=='SolidConcentration') then
             NSOLIDS=NSOLIDS+1
             ITSOLID(NSOLIDS)=IT
          end if
       END DO
       ewrite(2,*) 'Running with imported solids',NSOLIDS
       if(NSOLIDS.GT.0) then
          have_solids=.true.
          GOTTRAF=1
       end if
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

    IF(ACCTIM.EQ.0) SCFACTH0 = SCFACTH0*YOUR_SCFACTH0 ! IF WE DON'T HAVE/NOT USING EITHER OF THESE THEN THEY WILL JUST BE ONE
    !####################################################################

    NDWISP=.FALSE.
    IF((MLOC.GE.NLOC).and.((disopt(1)<70).or.(disopt(1)>90))) NDWISP=.TRUE.

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

       GETTAN=.TRUE.
       remesh = .false.

       GOTO 99771
    ELSE
       if (admesh) then
          call allocate(metric_tensor, extract_mesh(state(1), "CoordinateMesh"), "ErrorMetric")
       end if
       got_top_bottom_markers=any(ident==frees_ident) .or.&
            &            have_option('/geometry/ocean_boundaries')

    ENDif
    ! REMESHING ************************************************
    ! **********************************************************

    ! If POISON=-1 do first time step with POISON press determin
    do  IP=1,MAX(1,NPHASE)
       POISO2(IP)=POISON(IP)
       IF(POISON(IP).EQ.-1) POISO2(IP)=1
    end do

    NOPHAS=MAX(1,NPHASE)
    call compute_phase_uses_new_code_path(state, phase_uses_new_code_path)
    if (any(.not.phase_uses_new_code_path)) FLExit("The old code path is dead")

    call nonfield_allocation()

    NEWMES2=.TRUE.

    !     Determine the output format.
    call get_option('/io/dump_format', option_buffer, option_stat)
    if(option_stat /= 0) then
       FLAbort("You must specify a dump format and it must be vtk")
    else if(trim(option_buffer) /= "vtk") then
       FLAbort("You must specify a dump format and it must be vtk")
    end if

    if (have_option('/material_phase[0]/scalar_field::FreeSurface')) then
       !solve for free surface
       geostrophic_solver_option=2
    else
       !1== use mulgridvert to solve pressure and balanced pressure
       !0== don't, which is what we want as solving is done by PETSc
       geostrophic_solver_option=0
    end if
    nlevel=0


    !        Initialisation of distance to top and bottom field
    !        Currently only needed for free surface
    got_top_bottom_distance=has_scalar_field(state(1), "DistanceToTop")

    if (got_top_bottom_distance) then
       if (.not. got_top_bottom_markers) then
          FLAbort("There are no top and bottom boundary markers.")
       end if
       call CalculateTopBottomDistance(state(1), cograx)
    end if

    NOPHAS=MAX(1,NPHASE)

    CMCGET=.TRUE.
    GETC12=.TRUE.

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

    ! itinoi loop modification - allows it to dynamically update
    OITINOI = ITINOI

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
       !     itinoi loop modification - allows it to dynamically update
       ITINOI = OITINOI
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


          ! needed for rotated boundary conditions
          !nrtdr=nonods*3
          !if (d3) nrtdr=nonods*6

          !------------------------------------------------
          ! Addition for calculating drag force ------ jem 05-06-2008
          if (have_option("/imported_solids/calculate_drag_on_surface")) then
             call drag_on_surface(state)
          end if

          !     Addition for reading solids in - jem  02-04-2008
          if(have_solids) call solid_configuration(state(ss:ss),its,itinoi)

          !Setup tests for solid forces
          !  if(have_option('/imported_solids/tests')) then
          !     call get_option('/imported_solids/tests/ntests',nstests)             
          !     DO itest=1,nstests                
          !     END DO
          !  end if

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

             ITP=TPHASE(IT)

             IF((IDENT(IT)==KE_ident).OR.(IDENT(IT)==LKE_ident)) THEN
                ! Mellor Yamada turbulence modelling

                ! Set gotsal and gotvis to false until we (maybe) find
                ! them later
                gotsal=.false.
                gotvis=.false.

                do IT2=1,NTSOL
                   ! Search for IT corresponding to prognostic
                   ! MellorYamada fields
                   IF(IDENT(IT2)==KE_ident) THEN
                      it_KE=IT2
                   else if(ident(it2)==LKE_ident) then
                      it_LKE=IT2
                   else IF(IDENT(IT2)==ident_sal) THEN
                      ! Find salt index it_salt
                      it_salt=IT2
                      GOTSAL=.TRUE.
                   else if(ident(it2)==MY_vis_ident) then
                      ! Check for presence of MellorYamada visualisation fields
                      it_MY_vis=it2
                      gotvis=.true.
                   else if(ident(it2)==MY_diff_ident) then
                      ! Check for presence of MellorYamada visualisation fields
                      it_MY_diff=it2
                      gotvis=.true.
                   ENDIF
                   ! Find temp index it_temp
                   IF(BOUSIN(IT2)) it_temp=IT2
                END DO

                ! Necessary evil
                if(.not.gotsal) it_salt=it_temp

                ! trying to fix this insanity for new options
                if (has_scalar_field(state(1), "VerticalViscosity") .or. &
                     has_scalar_field(state(1), "VerticalDiffusivity")) then
                   ewrite(0,*) "Warning: Specified diagnostic fields for Mellor Yamada,"
                   ewrite(0,*) "This doesn't currently work,fields will be ignored,"
                end if
                gotvis=.false.
                ! hopefully not overwritten with gotvis==.false.:
                it_MY_vis=1
                it_MY_diff=1

                ewrite(2,*) "index of KineticEnergy field:", it_KE
                ewrite(2,*) &
                     & "index of TurbulentLengthScalexKineticEnergy field:", &
                     & it_LKE
                ewrite(2,*) "index of Temperature field:", it_temp
                ewrite(2,*) "index of Salinity field:", it_salt
                ewrite(2,*) "index of VerticalViscosity field:", it_MY_vis
                ewrite(2,*) "index of VerticalDiffusivity field:", it_MY_diff

                II=NOBCT(it_KE)

                ewrite(1,*) 'just before calling mellor_yamada_turbulence'

                !WE NEED A TEST CASE FOR THIS BEFORE CHANGING -- cjc
                call mellor_yamada_state(state(1),&
                     DENGAM(1),DENINI(1),TEMINI(1),GRAVTY,&
                     COGRAX,BSOUX(1),BSOUY(1),BSOUZ(1),&
                     (DISOPT(1).EQ.45).OR.(DISOPT(1).EQ.46),&
                                ! Salt...
                     EQNSTA(1),GAMDE2(1),GAMDE3(1), got_top_bottom_distance, &
                                ! BCS...
                     NOBCT(it_KE),BCT1W_mem(it_KE)%ptr,BCT2_mem(it_KE)%ptr,&
                     NOBCT(it_LKE),BCT1W_mem(it_LKE)%ptr, &
                     BCT2_mem(it_LKE)%ptr,DT,&
                                ! parallel stuff...
                     nnodp,para,halo_tag)

             ENDIF

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
                        & state=state(tphase(it)))

                ELSEIF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/finite_volume")) then

                   ! Solve the FV form of the equations.
                   call solve_advection_diffusion_fv(field_name=field_name_list(it), &
                        & state=state(tphase(it)))

                ELSEIF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/control_volumes")) then

                   ! Solve the pure control volume form of the equations
                   call solve_field_eqn_cv(field_name=trim(field_name_list(it)), &
                        state=state(field_state_list(it):field_state_list(it)), &
                        global_it=its)

                else if(have_option(trim(field_optionpath_list(it)) // &
                     & "/prognostic/spatial_discretisation/continuous_galerkin")) then

                   call solve_field_equation_cg(field_name_list(it), state(tphase(it)), dt)

                ELSEIF(have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/legacy_continuous_galerkin").or.&
                     have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/legacy_mixed_cv_cg").or.&
                     have_option(trim(field_optionpath_list(it))//&
                     & "/prognostic/spatial_discretisation/legacy_discretisation")) then

                   ! Solve the equation with at least some form of legacy cg in it.

                   assert(misnit==0)

                   ewrite(1, *) "Calling advdif from fluids"
                   CALL ADVDIF(&
                        & NEWMES2,&
                        & NOBCT(IT),BCT1W_mem(IT)%ptr, &
                        & BCT2_mem(IT)%ptr,&
                        & NONODS,TOTELE,r0,d0, &
                        & NLOC,NGI,&
                        & DISOTT(IT),NDISOT(IT),TTHETA(IT),TBETA(IT),TLUMP(IT),&
                        & NONODS,NONODS,&
                        & CGSOLT(IT),GMREST(IT),&
                        & SUFTEM(IT),SNLOC,SNGI, &
                                ! sub element modelling...
                        & NSUBTLOC(IT),NSUBNVLOC(ITP),&
                        & VERSIO,ISPHER2,&
                        & metric_tensor, FREDOP,&
                        & state(field_state_list(it):field_state_list(it)), field_name_list(it))
                   ewrite(1, *) "Exited advdif"

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
          IF(NAV) THEN
             !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             !     finding out which field is salinity - cjc
             !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             it_salt=1
             GOTSAL=.FALSE.
             do IT=1,NTSOL
                IF(IDENT(IT).EQ.IDEnt_sal) THEN
                   it_salt=IT
                   GOTsal=.TRUE.
                ENDIF
             END DO

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
             
             do ip = 1, max(1,nphase)
                poiso2(ip) = max(0,poison(ip))
             end do

             ! See subroutine SOLNAV for definitions of the
             ! following from DP onwards
             GETC12=CMCHAN
             CMCGET=CMCHAN
             IF(MVMESH) THEN
                GETC12=.true.
                CMCGET=.true.
             ENDIF
             ! end of NAV if.
          ENDIF

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

       ! Calculate next timestep size
       ! two different subroutines, one for single phase, one for multi phase (although there need be no difference . . .) 
       OLDDT =DT

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

       NOPHAS=MAX(1,NPHASE)

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

    if(isparallel()) then
       ewrite(1, *) "Halos registered in halo manager at simulation end: ", halo_registered_count()
    end if

  contains

    subroutine initial_allocation()
      !!< This subroutine is called immediately after the declaration of
      !!< allocatable memory in fluids to ensure that its allocation status is defined.
      logical, save :: initialised=.false.

      if(.not.initialised) then
         allocate(STVPXX(0), STVPYY(0), STVPZZ(0), &
              STVPYZ(0), STVPXZ(0), STVPXY(0), &
              YIELDF(0), YIELDFOLD(0))
         initialised = .true.
      end if

    end subroutine initial_allocation

    subroutine bc_mem_allocation()
      !!< This subroutine allocates memory for old style boundary conditions
      !!< hopefully to be deprecated soon
      !local variables
      integer :: i

      if(associated(bct1_mem)) then
         if(ntsol>0) then
            do i = 1, ntsol
               if(associated(bct1_mem(i)%ptr)) deallocate(bct1_mem(i)%ptr)
            end do
         end if
         deallocate(bct1_mem)         
      end if

      if(associated(bct1w_mem)) then
         if(ntsol>0) then
            do i = 1, ntsol
               if(associated(bct1w_mem(i)%ptr)) deallocate(bct1w_mem(i)%ptr)
            end do
         end if
         deallocate(bct1w_mem)
      end if

      if(associated(bct2_mem)) then
         if(ntsol>0) then
            do i = 1, ntsol
               if(associated(bct2_mem(i)%ptr)) deallocate(bct2_mem(i)%ptr)
            end do
         end if
         deallocate(bct2_mem)
      end if

      allocate(bct1_mem(0))
      allocate(bct1w_mem(0))
      allocate(bct2_mem(0))

      !do i = 1, ntsol
      !   allocate( bct1_mem(i)%ptr(0) )
      !   allocate( bct1w_mem(i)%ptr(0) )
      !   allocate( bct2_mem(i)%ptr(0) )
      !end do

    end subroutine bc_mem_allocation

    subroutine bcw_mem_allocation()
      !!< This subroutine allocates memory for old style boundary conditions
      !!< in acceleration form
      !!< hopefully to be deprecated soon
      !local variables
      integer :: i

      if(associated(bct1w_mem)) deallocate(bct1w_mem)
      if(associated(bcu1w_mem)) deallocate(bcu1w_mem)
      if(associated(bcv1w_mem)) deallocate(bcv1w_mem)
      if(associated(bcw1w_mem)) deallocate(bcw1w_mem)

      allocate(bcu1w_mem(nphase))
      allocate(bcv1w_mem(nphase))
      allocate(bcw1w_mem(nphase))
      call nullify(bcu1w_mem)
      call nullify(bcv1w_mem)
      call nullify(bcw1w_mem)
      if(ntsol>0) then
         allocate(bct1w_mem(ntsol))
         call nullify(bct1w_mem)
         do i = 1, ntsol
            if(associated( bct1w_mem(i)%ptr )) deallocate(bct1w_mem(i)%ptr)
            allocate( bct1w_mem(i)%ptr(nobct(i)) )
         end do
      else
         allocate(bct1w_mem(1))
         allocate(bct1w_mem(1)%ptr(1))
      end if

      do i = 1, nphase
         if(associated( bcu1w_mem(i)%ptr )) deallocate(bcu1w_mem(i)%ptr)
         if(associated( bcv1w_mem(i)%ptr )) deallocate(bcv1w_mem(i)%ptr)
         if(associated( bcw1w_mem(i)%ptr )) deallocate(bcw1w_mem(i)%ptr)
         allocate( bcu1w_mem(i)%ptr(nobcu(i)) )
         allocate( bcv1w_mem(i)%ptr(nobcv(i)) )
         allocate( bcw1w_mem(i)%ptr(nobcw(i)) )
      end do
    end subroutine bcw_mem_allocation

    subroutine nonfield_allocation()
      !!< Allocate memory for things that aren't fields but need to exist in the fluids loop.
      !!< Allocation status must be defined so that they can be initially deallocated to prevent
      !!< memory leaks through adapts

      if(allocated(STVPXX)) deallocate(STVPXX)
      if(allocated(STVPYY)) deallocate(STVPYY)
      if(allocated(STVPZZ)) deallocate(STVPZZ)
      if(allocated(STVPYZ)) deallocate(STVPYZ)
      if(allocated(STVPXZ)) deallocate(STVPXZ)
      if(allocated(STVPXY)) deallocate(STVPXY)
      if(allocated(YIELDF)) deallocate(YIELDF)
      if(allocated(YIELDFOLD)) deallocate(YIELDFOLD)

      if(SOLIDS.GT.0) then
         allocate(STVPXX(TOTELE*NGI), STVPYY(TOTELE*NGI), STVPZZ(TOTELE*NGI), &
              STVPYZ(TOTELE*NGI), STVPXZ(TOTELE*NGI), STVPXY(TOTELE*NGI), &
              YIELDF(TOTELE*NGI), YIELDFOLD(TOTELE*NGI))
         STVPXX = 0.0
         STVPYY = 0.0
         STVPZZ = 0.0
         STVPYZ = 0.0
         STVPXZ = 0.0
         STVPXY = 0.0
         YIELDF = 0.0
         YIELDFOLD = 0.0
      else
         ALLOCATE(STVPXX(0), STVPYY(0), STVPZZ(0), &
              STVPYZ(0), STVPXZ(0), STVPXY(0), &
              YIELDF(0), YIELDFOLD(0))
      end if

    end subroutine nonfield_allocation

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

