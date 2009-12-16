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

#include "confdefs.h"
#include "fdebug.h"

module geostrophic
  use quadrature
  use fields
  use Sparse_Tools
  use allsorts
  use solvers
  use EventCounter
  use fetools, only: X_,Y_,Z_
  use FLDebug
  use state_module
  use vtk_interfaces
  use Element_Numbering, only: ilink2, silink2
  use mesh_connections
  use verticaldg_module
  use shape_transformations
  use Coordinates
  use vertical_extrapolation_module
  use Tidal_module
  use boundary_conditions
  use tr2d_module
  use global_parameters, only: new_options, OPTION_PATH_LEN
  use spud
  use position_in_matrix
  use rotated_boundary_conditions_legacy
  use diff3d_module
  use assnav_module
  use parallel_tools
  use spud
  use FLDebug
  use reduced_model
  use sml
  use free_surface_module
  use halos

  implicit none

  private

  public GEOELI1P, geostrophic_check_options

contains

  SUBROUTINE GEOELI1P(VECX,VECY,VECZ, &
       NSUBVLOC,VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
       SOUVECX,SOUVECY,SOUVECZ, &
       U,V,W, &
       INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
       SOURCX,SOURCY,SOURCZ,PSIPRE,P, &
       NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB, &
       NOBCFSimp,BCT1WFSimp,BCT2FSimp, &
       X,Y,Z,XOLD,YOLD,ZOLD,XORIG,YORIG,ZORIG, &
       UG,VG,WG,MVMESH,  &
       N2,NLX2,NLY2,NLZ2, WEIGHT2, NGI2, &
       NLOC,SNLOC, &
       TOTELE,D3,DCYL, &
       GEOBAL,&
       NDGLNO, XONDGL, &
       NONODS,XNONOD, &
       FINDRM,CENTRM,COLM,NCOLM, &
       FREDOP, &
       NLEVEL,GSOLOP,FREESDISOTT, &
       DISOPT, &
       halo_tag, NNODP, &
       scfacth0,NEWMES2, &
       GOTTOPDIS, XABSOR,YABSOR,ZABSOR, &
       GOTBOY,BSOUX,BSOUY,BSOUZ, &
       COGRAX, &
       VGRAVX,VGRAVY,VGRAVZ, &
                                ! New stuff for free surface...
       NEXTTIM,GRAVTY,DT,ACCTIM,ISPHERE,EQTD,GOT_EQTD, &
       NOFILT, &
       C1T,C2T,C3T,FINDCT,COLCT,NCT, &
                                ! the state object containing all fields of this 'phase'
       state,SOLSTAGE2,PHSOLVE, &
                                ! path to either BalancePressure or FreeSurface field options
       option_path)
    ! If SOLSTAGE2=0 then return immediatly. 
    ! If SOLSTAGE2=1 then forming bouyancy/balance contribution into VECX. 
    ! If SOLSTAGE2=2 then solve for free surface and update mesh...
    ! If SOLSTAGE2=3 then do everything
    ! If SOLSTAGE2=3 and (.NOT.NEXTTIM .OR. WETDRY=1) then set SOLSTAGE=1
    ! If PHSOLVE then solve the vertical DG eqns else use the previous 
    ! hydrostic pressure value. If PHSOLVE=.FALSE. then dont calculate C1t etc.

    ! THIS SUB SOLVES AND ELIPTIC or DG EQN FOR BALANCE.
    ! If GEOBAL=-11 or -21 then include bouyancy in forcing terms...
    ! If GOEBAL=-1? then have a linear variation in goe-hyd pressure
    ! If GEOBAL=-2? then have quadratic variation.
    ! If GEOBAL=-23 then have quadratic variation in the vertical with a DG method.
    ! If GEOBAL=-24 then have quadratic variation in the vertical with a DG method but
    ! use no b.c's when putting into the momentum eqns.
    ! If GEOBAL=-25 same as -24 but do not solve for verticaldg pressure 
    ! but use a quadratic pressure i.e. P1_(SGS or DG)-P2. 
    ! If GEOBAL=-26 same as -24 but do not solve for verticaldg pressure 
    ! & solve for p and free surface together.
    ! If GEOBAL=-27 same as -23 but solve for p and free surface together.
    ! It puts geostrophic pressure into VECX,VECY,VECZ
    ! this is converted into a source: SOUVECX,SOUVECY,SOUVECZ
    ! which is calculated from U,V,W
    ! Put the result into momentum eqns VECX=VEXC+N_i dP/dx etc
    IMPLICIT NONE
    ! Assume we are going to use a quadratic element for PG.
    INTEGER MLOC,NLEVEL,GSOLOP
    INTEGER FREESDISOTT
    INTEGER DISOPT
    ! FREEDDISOTT controls the tidal modelling >100 tidal modelling
    ! decode stuff after 100 to get the tidal constituants.
    LOGICAL GETDIAG,RESETDIA
    REAL INFINY

    PARAMETER(GETDIAG=.TRUE.,RESETDIA=.FALSE.)
    PARAMETER(INFINY=1.E+16)
    INTEGER NONODS,XNONOD
    INTEGER TOTELE,NLOC,SNLOC,NGI2
    REAL N2(NLOC,NGI2),NLX2(NLOC,NGI2),NLY2(NLOC,NGI2),NLZ2(NLOC,NGI2)
    REAL WEIGHT2(NGI2)

    REAL,target,intent(inout):: X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL XOLD(XNONOD),YOLD(XNONOD),ZOLD(XNONOD)
    REAL XORIG(XNONOD),YORIG(XNONOD),ZORIG(XNONOD)
    REAL UG(NONODS),VG(NONODS),WG(NONODS)
    LOGICAL MVMESH
    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    INTEGER NSUBVLOC
    REAL VECX_SGSADD(TOTELE*NSUBVLOC),VECY_SGSADD(TOTELE*NSUBVLOC)
    REAL VECZ_SGSADD(TOTELE*NSUBVLOC)
    REAL SOUVECX(NONODS),SOUVECY(NONODS),SOUVECZ(NONODS)
    REAL U(NONODS),V(NONODS),W(NONODS)
    INTEGER INUSUB,NSUBNVLOC
    REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
    REAL NWSUB(TOTELE*NSUBNVLOC)
  INTEGER NSOGRASUB
  REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
  REAL SOZGRASUB(TOTELE*NSOGRASUB)
    REAL SOURCX(NONODS),SOURCY(NONODS),SOURCZ(NONODS)
    REAL PSIPRE(NONODS)
    INTEGER NOBCFSimp
    REAL BCT1WFSimp(NOBCFSimp)
    INTEGER BCT2FSimp(NOBCFSimp)
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER FREDOP
    REAL P(FREDOP)
    INTEGER NCT,FINDCT(FREDOP+1),COLCT(NCT)
    INTEGER SOLSTAGE2
    LOGICAL PHSOLVE
    REAL C1T(NCT),C2T(NCT),C3T(NCT)

    LOGICAL D3,DCYL
    INTEGER GEOBAL
    INTEGER NCOLM,FINDRM(NONODS+1),CENTRM(NONODS),COLM(NCOLM)
    INTEGER halo_tag, NNODP
    real scfacth0
    LOGICAL NEWMES2
    LOGICAL GOTTOPDIS
    REAL XABSOR(NONODS),YABSOR(NONODS),ZABSOR(NONODS)
    LOGICAL GOTBOY
    REAL BSOUX,BSOUY,BSOUZ
    LOGICAL COGRAX
    REAL VGRAVX(NONODS),VGRAVY(NONODS),VGRAVZ(NONODS)
    LOGICAL NEXTTIM,GOT_EQTD
    INTEGER NOFILT,WETDRY,FULL3D
    INTEGER ELETOPBC,ELLBAL,ELLSTAGTOPBC
    REAL EQTD(NONODS)
    REAL GRAVTY,DT,ACCTIM
    INTEGER ISPHERE
    type(state_type), intent(inout):: state
    character(len=*), intent(in):: option_path

    ! Local variables...
    INTEGER NPGCOLM
    INTEGER PARA
    INTEGER I,ELE
    INTEGER NGI
    REAL L1(70),L2(70),L3(70),L4(70)

    INTEGER PGNODS, PGNODSP
    INTEGER::MeshCount=-1
    SAVE MeshCount
    INTEGER, PARAMETER::halotagT10=3
    INTEGER PROCNO, CNT, ierr
    INTEGER NOD,PPSET
   
    REAL CTHETA,ATHETA,FSTHETA
   
    LOGICAL TAKVERT,FSTAKVERT,GETFREES
    INTEGER ADJGRAV,USEPOTGRA
    
    INTEGER GOTTID,FREES3D,GEOSHORT,IMPABS,TIDALLCOM,TLOVE
    INTEGER STATID,NOFSTA
    INTEGER IDGBAL,IDGBOYONLY,ISMOOTH,IDGBCMOM
    INTEGER ELLBOYONLY,ELLCORONLY,SOLSTAGE
    INTEGER ILOC,JLOC,KLOC,IGL,JGL,PNOD
   
    INTEGER IIDUM(1),WHICH_FREE_OPT(20)
    REAL RD(1)
    LOGICAL ONLY3D,GETELEORD,SMOOTH,shared_p_fs
    INTEGER fourthfrees_diss
    ! Only the first time we enter this sub is HYDSTART=1
    INTEGER HYDSTART
    SAVE HYDSTART
    DATA HYDSTART /1/
    ! *********ALLOCATE MEMORY**************
    REAL, ALLOCATABLE, DIMENSION(:)::tempvec
    REAL, ALLOCATABLE, DIMENSION(:)::RZER
    REAL, ALLOCATABLE, DIMENSION(:)::SHORTDIFMATH
    REAL, ALLOCATABLE, DIMENSION(:)::WETDRYMAT
    REAL, ALLOCATABLE, DIMENSION(:)::SUFML
    REAL, ALLOCATABLE, DIMENSION(:)::DIFVECX
    REAL, ALLOCATABLE, DIMENSION(:)::DIFVECY
    REAL, ALLOCATABLE, DIMENSION(:)::DIFVECZ
    REAL, ALLOCATABLE, DIMENSION(:)::SHORTDIFVECV
    REAL, ALLOCATABLE, DIMENSION(:)::SHORTDP
    REAL, ALLOCATABLE, DIMENSION(:)::DIFVECV
    REAL, ALLOCATABLE, DIMENSION(:)::DIFVECN
    REAL, ALLOCATABLE, DIMENSION(:)::DIFRHS
    REAL, ALLOCATABLE, DIMENSION(:)::DIFRHSND
    REAL, ALLOCATABLE, DIMENSION(:)::EQETID
    REAL, ALLOCATABLE, DIMENSION(:)::MLND
    REAL, ALLOCATABLE, DIMENSION(:)::MLNDVERT
    INTEGER, ALLOCATABLE, DIMENSION(:)::IWORKS
    REAL, ALLOCATABLE, DIMENSION(:)::UMEAN
    REAL, ALLOCATABLE, DIMENSION(:)::VMEAN
    REAL, ALLOCATABLE, DIMENSION(:)::WMEAN
    REAL, ALLOCATABLE, DIMENSION(:)::PNORMX
    REAL, ALLOCATABLE, DIMENSION(:)::PNORMY
    REAL, ALLOCATABLE, DIMENSION(:)::PNORMZ
    REAL, ALLOCATABLE, DIMENSION(:)::DP
    REAL, ALLOCATABLE, DIMENSION(:)::NDENST

    REAL, ALLOCATABLE, DIMENSION(:)::FREESTHETA 

    REAL, ALLOCATABLE, DIMENSION(:)::ORIG_XABSOR
    REAL, ALLOCATABLE, DIMENSION(:)::ORIG_YABSOR
    REAL, ALLOCATABLE, DIMENSION(:)::ORIG_ZABSOR

    REAL, ALLOCATABLE, SAVE, DIMENSION(:)::PG
    REAL, ALLOCATABLE, SAVE, DIMENSION(:)::FREESOLD
    REAL, ALLOCATABLE, SAVE, DIMENSION(:)::PGRAVX
    REAL, ALLOCATABLE, SAVE, DIMENSION(:)::PGRAVY
    REAL, ALLOCATABLE, SAVE, DIMENSION(:)::PGRAVZ
    REAL, ALLOCATABLE, SAVE, DIMENSION(:)::POTGRA
    REAL, ALLOCATABLE, SAVE, DIMENSION(:)::DGP
    INTEGER, ALLOCATABLE, SAVE, DIMENSION(:)::ELEORD
    INTEGER, ALLOCATABLE, SAVE, DIMENSION(:)::COLELE4
    REAL, ALLOCATABLE, DIMENSION(:,:)::N
    REAL, ALLOCATABLE, DIMENSION(:,:)::NLX
    REAL, ALLOCATABLE, DIMENSION(:,:)::NLY
    REAL, ALLOCATABLE, DIMENSION(:,:)::NLZ
    REAL, ALLOCATABLE, DIMENSION(:)::WEIGHT

    Real, Allocatable :: nodfrees(:)
    Real, Allocatable :: wdgravity(:)

    !  ! Hard coded parameters!! should be turned into options
    !  real, parameter :: dd00=0.05
    real::  rnonl_coefficient

    real, allocatable :: obcflag(:)

    ! Shape functions FOR QUADRATIC TET...
    real, dimension(:), allocatable:: m, mlx, mly, mlz

    type(csr_matrix) :: pg_matrix, prolongator
    type(csr_sparsity) :: pg_sparsity
    integer, dimension(:), pointer:: PGFINDRM, PGCENTRM, PGCOLM
    real, dimension(:), allocatable:: difmath

    ! new style shape functions and meshes
    type(quadrature_type):: quad
    type(element_type):: pg_shape
    type(mesh_type):: pg_mesh
    type(scalar_field):: free_surface
    type(vector_field), pointer:: original_coordinates
    ! fields pointers pointing to fields/meshes in state:
    type(vector_field), pointer:: coordinates
    type(mesh_type), pointer:: x_mesh, velocity_mesh
    type(scalar_field), pointer:: topdis_field=>null(), botdis_field=>null()
    ! pointing to the values of these, for legacy purposes:
    real, dimension(:), pointer:: topdis, botdis
    ! pointing to free_surface scalar_field:
    real, pointer, dimension(:)::frees
    ! the quadratic pressure mesh, points at pg_mesh%ndglno
    integer, pointer, dimension(:):: pgndglno

    character(len=OPTION_PATH_LEN):: my_option_path

    logical, parameter:: DUMPMESH=.false.
 
    integer dimension
    integer :: inod, jnod

    ! the following declarations are for wetting/drying
    integer :: nloop_wetdry
    logical, save :: firstcall=.true.
    real, save :: coef_abs, ee00, dd00
    integer :: ios_wetdry

    SOLSTAGE=SOLSTAGE2
    IF(SOLSTAGE.EQ.0) RETURN
    ewrite(2, *) "SUBROUTINE GEOELI1P()"

    Allocate(wdgravity(totele))  
    wdgravity=1.

    if (gottopdis) then
       topdis_field => extract_scalar_field(state, "DistanceToTop")
       botdis_field => extract_scalar_field(state, "DistanceToBottom")
       topdis => topdis_field%val
       botdis => botdis_field%val
    else
       ! If gottopdis==.false. topdis and botdis are not used (hopefully), 
       ! but still passed as arrays to subroutines, so we need the pointers
       ! to point to a xnonod-length array. If (heaven forbid) they are used
       ! anyway we should be able to notice the difference with the following choice:
       topdis => x
       botdis => y
    end if

    x_mesh => extract_mesh(state, "CoordinateMesh")
    velocity_mesh => extract_mesh(state, "VelocityMesh")
    coordinates => extract_vector_field(state, "Coordinate")
    if (mvmesh) then
       original_coordinates => extract_vector_field(state, "OriginalCoordinate")
    end if

    allocate(obcflag(nonods))

    ! Calculate the shape functions just for this sub - it needs to 
    ! have large quadrature...

    ! this sub does not work in this place...
    !   call QuadTetShapes(n, nlx, nly, nlz, nloc, ngi)
    NGI=NGI2
    IF(NGI2.EQ.1) NGI=11

    ALLOCATE(N(NLOC,NGI))
    ALLOCATE(NLX(NLOC,NGI))
    ALLOCATE(NLY(NLOC,NGI))
    ALLOCATE(NLZ(NLOC,NGI))
    ALLOCATE(WEIGHT(NGI))
    
    shared_p_fs=(GEOBAL.EQ.-26).OR.(GEOBAL.EQ.-27)

    IF(NGI2.EQ.1) THEN
       ! Assume we are using 1 pt quad with tets...
       CALL TRIQUA(L1, L2, L3, L4, weight, .true.,NGI)
       CALL SHATRI(L1, L2, L3, L4, weight, .true., &
            nLOC,NGI, &
            n,nLX,nLY,nLZ)
    ELSE
       N=N2
       NLX=NLX2
       NLY=NLY2
       NLZ=NLZ2
       weight=weight2
    ENDIF

    PROCNO = GetProcNo()

    ewrite(2, *) 'gotboy:',gotboy

    ewrite(2, *) 'r2norm(sourcx,nonods,0):',r2norm(sourcx,nonods,0)
    ewrite(2, *) 'r2norm(sourcy,nonods,0):',r2norm(sourcy,nonods,0)
    ewrite(2, *) 'r2norm(sourcz,nonods,0):',r2norm(sourcz,nonods,0)

    FULL3D=0

    ! *********** define the options ******************************
    ! the options for the free surface:
    ! FREESDISOTT is decoded to obtain the free surface options
    ! (FREES3D,fourthfrees_diss,GEOSHORT,IMPABS,TIDALLCOM,TLOVE) (defaults are 0)
    ! FREESDISOTT is converted to a binary no to switch
    ! on and off the 1st 6 options.
    ! GOTTID switches on tidal modelling
    ! or IF((ABS(FREESDISOTT).GE.2**11).OR.(TIDALLCOM.EQ.1)) GOTTID=1
    ! Switching on individual tidal compoinent starts off at 2**11 for the
    ! 1st component + 2**12 for the 2nd etc.
    !
    ! If FREES3D=1 then form a full 3D system for the free surface.
    ! Switch on 4th order dissipation when fourthfrees_diss=1
    ! If GEOSHORT=1 then have a low order representation (linear) of the
    ! free surface - not including geostrophic balance which can still be quadratic.
    ! If IMPABS=1 then have an implicit absorption in the free surface
    ! which is good for wetting/drying
    ! TIDALLCOM=1 to switch on all tidal components.
    ! TLOVE=1 switches on a love no. of 0.3 i.e. (1.-0.3) multiplied
    ! by the tidal forcing magnitude.
    ! STATID= use a static tidal force for testing
    ! NOFSTA=Switch off free surface stabalisation term.
    ! NOFILT=Switch off pressure stabilisation term in cty eqn.
    ! WETDRY the we switch on wetting/dry (absorbtion and lumping of free s. eqns)
    ! For example suppose (option DISOTT(7) 27) then the options are:
    ! FREES3D=           1
    ! fourthfrees_diss=  0
    ! GEOSHORT=          1
    ! IMPABS=            1
    ! TIDALLCOM=         0
    ! TLOVE=             0
    ! STATID=            0
    ! NOFSTA=            0
    ! NOFILT=            0
    ! WETDRY=            0

    if (new_options) then
       my_option_path=trim(option_path)//'/prognostic/'
       ! these should all be logicals:
       frees3d=0
       if (have_option(trim(my_option_path)//'/spatial_discretisation/free_surface_3D')) frees3d=1
       fourthfrees_diss=0
       if (have_option(trim(my_option_path)//'/spatial_discretisation/fourth_order_dissipation')) fourthfrees_diss=1
       geoshort=0
       if (have_option(trim(my_option_path)//'/spatial_discretisation/low_order_free_surface')) geoshort=1
       wetdry=0
       if (have_option(trim(my_option_path)//'/spatial_discretisation/wetting_drying')) wetdry=1
       if (have_option(trim(my_option_path)//'/spatial_discretisation/default_free_surface_filter')) then
          nofsta=0
          if (wetdry==0) then
             rnonl_coefficient=0.01
          else
             rnonl_coefficient=1.0
          end if
       else if (have_option(trim(my_option_path)//'/spatial_discretisation/user_specified_free_surface_filter')) then
          nofsta=0
          call get_option(trim(my_option_path)//'/spatial_discretisation/user_specified_free_surface_filter&
               &/non_linear_filter_coefficient', rnonl_coefficient)
       else if (have_option(trim(my_option_path)//'/spatial_discretisation/switch_off_free_surface_filter')) then
          nofsta=1
       else
          nofsta=-66600
       end if
       gottid=0
       tidallcom=0
       tlove=0
       statid=0
       if (have_option(trim(my_option_path)//'/spatial_discretisation/tidal_forcing'))  then
       gottid=1
          my_option_path=trim(my_option_path)//'/spatial_discretisation/tidal_forcing'
          if (have_option(trim(my_option_path)//'/M2')) freesdisott=freesdisott+2**11
          if (have_option(trim(my_option_path)//'/S2')) freesdisott=freesdisott+2**12
          if (have_option(trim(my_option_path)//'/N2')) freesdisott=freesdisott+2**13
          if (have_option(trim(my_option_path)//'/K2')) freesdisott=freesdisott+2**14
          if (have_option(trim(my_option_path)//'/K1')) freesdisott=freesdisott+2**15 
          if (have_option(trim(my_option_path)//'/O1')) freesdisott=freesdisott+2**16
          if (have_option(trim(my_option_path)//'/P1')) freesdisott=freesdisott+2**17
          if (have_option(trim(my_option_path)//'/Q1')) freesdisott=freesdisott+2**18
          if (have_option(trim(my_option_path)//'/all_tidal_components')) tidallcom=1
          if (have_option(trim(my_option_path)//'/love_number')) tlove=1
          if (have_option(trim(my_option_path)//'/static_tidal_force')) statid=1
       end if
       if (have_option(trim(my_option_path)//'/Legacy_Free_Surface')) then
       call get_option(trim(my_option_path)//'/Legacy_Free_Surface',FREESDISOTT)
       end if

    else
       WHICH_FREE_OPT(1:12) = 0
       CALL DECODBINTID(WHICH_FREE_OPT,12,ABS(FREESDISOTT))

       ! 1st entry is empty.
       FREES3D=WHICH_FREE_OPT(1)
       fourthfrees_diss=WHICH_FREE_OPT(3)
       GEOSHORT=WHICH_FREE_OPT(4)
       IMPABS=WHICH_FREE_OPT(5)
       TIDALLCOM=WHICH_FREE_OPT(6)
       TLOVE=WHICH_FREE_OPT(7)
       STATID=WHICH_FREE_OPT(8)
       NOFSTA=WHICH_FREE_OPT(9)
       NOFILT=WHICH_FREE_OPT(10)
       WETDRY=WHICH_FREE_OPT(11)
    end if

    if(wetdry.eq.1) then
       if(firstcall) then
          open(1,file="absorption.dat",status='old',iostat=ios_wetdry)
          if(ios_wetdry.gt.0) then
             FLAbort("Error, no absorption.dat file in current directory")
          end if
          rewind(1)
          read(1,*)
          read(1,*) coef_abs, ee00,dd00
          close(1)
          firstcall=.false.
       end if
    end if

    ewrite(2,*) 'FREESDISOTT',FREESDISOTT
    ewrite(2,*) 'the options are:'
    ewrite(2,*) '2^0:FREES3D=',FREES3D
    ewrite(2,*) '2^2:fourthfrees_diss=',fourthfrees_diss
    ewrite(2,*) '2^3:GEOSHORT=',GEOSHORT
    if (.not. new_options) then
       ewrite(2,*) '2^4:IMPABS=',IMPABS
    end if
    ewrite(2,*) '2^5:TIDALLCOM=',TIDALLCOM
    ewrite(2,*) '2^6:TLOVE=',TLOVE
    ewrite(2,*) '2^7:STATID=',STATID
    ewrite(2,*) '2^8:NOFSTA=',NOFSTA
    ewrite(2,*) '2^9:NOFILT=',NOFILT
    ewrite(2,*) '2^10:WETDRY=',WETDRY

    ! GOTTID switches on tidal modelling
    GOTTID=0
    IF((ABS(FREESDISOTT).GE.2**11).OR.(TIDALLCOM.EQ.1)) GOTTID=1
    ! have a DG for vertical balanced pressure...
    IDGBAL=0
    IF((GEOBAL<=-23).AND.(GEOBAL>=-27)) IDGBAL=1
    ! Solve for free surface.
    GETFREES=(ABS(GSOLOP).GE.2)
    ! If ADJGRAV=1 then adjust gravity direction/magnetude for mesh
    ! If ADJGRAV=0 then dont adjust gravity.
    !    ADJGRAV=1
    ADJGRAV=0
    IF(.NOT.GOTBOY) ADJGRAV=0
    USEPOTGRA=1
    ! If USEPOTGRA=1 then use gravitational potential instead of gravity.
    IF(ADJGRAV.EQ.0) USEPOTGRA=0
    ! Switch on 4th order dissipation when fourthfrees_diss=1
    ! fourthfrees_diss=0
    ! If GEOSHORT=1 then have a low order representation (linear) of the
    ! free surface - not including geostrophic balance which can still be quadratic.
    ! GEOSHORT=0
    ! If IMPABS=1 then have an implicit absorbtion in the free surface
    ! which is good for wetting/drying
    !    IMPABS=0
    if (new_options) then
       my_option_path=trim(option_path)//'/prognostic/temporal_discretisation'
       call get_option(trim(my_option_path)//'/absorption_theta', atheta, default=0.0)
       ! This is the FSTHETA time stepping parameter for the free surface.
       ! (is required by schema, default is for legacy purposes):
       call get_option(trim(my_option_path)//'/theta', fstheta, default=1.0)
    else
       ATHETA=0.0
       IF(IMPABS.EQ.1) ATHETA=1.0
       ! This is the FSTHETA time stepping parameter for the free surface.
       FSTHETA=1.0
    endif
    !    FSTHETA=THETA
    ! *********** define the options ******************************

    ! *********** define GEOBAL for vertical DG and elliptic solv ***
    ! IF ELETOPBC=1 then apply b.c's to the top surface only in
    ! the elliptic solver when we have a free surface.
    ! Its dangerous to use ELETOPBC=0 with a free surface.
    ELETOPBC=1
    ELLBAL=1
    IF(IDGBAL.EQ.1) ELLBAL=0
    IDGBCMOM=1
    IDGBOYONLY=0
    ELLBOYONLY=0
    ELLCORONLY=0
    ELLSTAGTOPBC=0
    ! If ELLSTAGTOPBC=1 set the values of just the corner nodes at the surface.
    ISMOOTH=1
    ! If SMOOTH then smooth out the forces from vertical DG method.
    IF(GEOBAL.LT.-30) THEN
       WHICH_FREE_OPT(1:18) = 0
       CALL DECODBINTID(WHICH_FREE_OPT,18,ABS(GEOBAL))
       IDGBAL=WHICH_FREE_OPT(6)
       IDGBOYONLY=WHICH_FREE_OPT(7)
       ISMOOTH=WHICH_FREE_OPT(8)
       IDGBCMOM=WHICH_FREE_OPT(9)
       ELLBAL=WHICH_FREE_OPT(11)
       ELLBOYONLY=WHICH_FREE_OPT(12)
       ELLCORONLY=WHICH_FREE_OPT(13)
       ELETOPBC=WHICH_FREE_OPT(14)
       ELLSTAGTOPBC=WHICH_FREE_OPT(15)
    ENDIF

    SMOOTH=(ISMOOTH.EQ.1)
    IF(IDGBOYONLY.EQ.1) THEN
       ELLBOYONLY=0
       ELLCORONLY=1
    ENDIF
    IF((IDGBAL.EQ.1).AND.(IDGBOYONLY.EQ.1).AND.(ELLBAL.EQ.0)) THEN
       EWRITE(0,*) 'NEED TO AMEND DIFF3D TO STOP ROTATION BEING IN VECX ETC'
       FLAbort("STOP 3292")
    ENDIF
    IF((ELLBOYONLY.EQ.1).AND.(ELLBAL.EQ.1)) THEN
       EWRITE(0,*) 'NEED TO AMEND DIFF3D TO STOP ROTATION BEING IN VECX ETC'
       FLAbort("STOP 8322")
    ENDIF

    ! if GEOBAL.LT.30 use old options else use the options below.
    ! If IDGBAL=1 switch on vertical DG.
    ! If IDGBOYONLY=1 treat only bouyancy in vertical DG.
    ! If ISMOOTH=1 smooth the resulting forces from vertical DG.
    ! If IDGBCMOM=1 include element DG b.c's into VECX,VECY,VECZ.
    ! If -balnk option
    ! If ELLBAL=1 switch on elliptic eqn solver for balance.
    ! If ELLBOYONLY=1 treat bouyancy only in elliptic balance solver.
    ! If ELLCORONLY=1 treat Coriolis only in elliptic balance solver.
    ! IF ELETOPBC=1 then apply b.c's to the top surface only in
    ! the elliptic solver when we have a free surface.
    ! If ELLSTAGTOPBC=1 set the values of just the corner nodes at the surface.

    ewrite(2,*) 'the options are:'
    ewrite(2,*) '2^5:IDGBAL=',IDGBAL
    ewrite(2,*) '2^6:IDGBOYONLY=',IDGBOYONLY
    ewrite(2,*) '2^7:ISMOOTH=',ISMOOTH
    ewrite(2,*) '2^8:IDGBCMOM=',IDGBCMOM
    ewrite(2,*) '2^9:...='
    ewrite(2,*) '2^10:ELLBAL=',ELLBAL
    ewrite(2,*) '2^11:ELLBOYONLY=',ELLBOYONLY
    ewrite(2,*) '2^12:ELLCORONLY=',ELLCORONLY
    ewrite(2,*) '2^13:ELETOPBC=',ELETOPBC
    ewrite(2,*) '2^14:ELLSTAGTOPBC=',ELLSTAGTOPBC
    ! *********** define GEOBAL for vertical DG and elliptic solv ***


    ! **********initialisation of basic variables*******
    IF(GEOBAL.LE.-20) THEN
       ewrite(2,*) "quadratic pressure variation"
       MLOC=10
       IF(NLOC.EQ.8) MLOC=27
       IF(NLOC.EQ.8) THEN
          FLAbort('QUADRATIC DOES NOT EXIST YET-COMPLEX')
       ENDIF
    ELSE
       ewrite(2,*) "linear pressure variation"
       MLOC=4
       IF(NLOC.EQ.8) MLOC=8
    ENDIF

    ! Shape functions FOR QUADRATIC TET...
    allocate(m(mloc*ngi), mlx(mloc*ngi), mly(mloc*ngi), mlz(mloc*ngi))

    ASSERT(MLOC.EQ.10)

    ! The mid side nodes and corner nodes define...
    CALL GetEventCounter(EVENT_ADAPTIVITY, CNT)
    IF(MeshCount.NE.CNT) THEN
       ! after an adapt, recreate the quadratic mesh
       if (d3) then
          dimension=3
       else
          dimension=2
       end if
       quad=make_quadrature(nloc, dimension, ngi=ngi)

       ! shape and mesh of the balanced pressure
       ! (since currently only mloc==10 works, we assume quadratic):
       pg_shape=make_element_shape(nloc, dimension, 2, quad)

       ALLOCATE(PGNDGLNO(TOTELE*MLOC))
       ewrite(1,*) 'entering midsnods'
       CALL MIDSNODS(FINDRM,NONODS,PGNODS,        &
            CENTRM,COLM,NCOLM,                    &
            NDGLNO,PGNDGLNO,TOTELE,NLOC,MLOC)
       ewrite(1, *) "finished midsnods"

       pg_mesh=wrap_mesh(pgndglno, pg_shape, name="BalancePressureMesh")
       ! ugly hack, to make sure pgndglno gets deallocated 
       ! when pg_mesh is deallocated with state
       pg_mesh%wrapped=.false.
       call add_faces(pg_mesh, model=x_mesh)

       ! add this mesh to the state object
       call insert(state, pg_mesh, "BalancePressureMesh")

       if (DUMPMESH) then
          ! write the position field interpolated on the balanced pressure mesh
          call vtk_write_fields('pgmesh', cnt+1, coordinates, &
               pg_mesh, vfields=(/ coordinates /) )
       end if
       ! drop our references
       call deallocate(pg_mesh)
       call deallocate(pg_shape)
       call deallocate(quad)

#ifndef NDEBUG
       if(IsParallel()) then
          assert(halo_tag.eq.2)
          CALL flcomms_test2(halo_tag, X, 1, 1, 0, NONODS, NNODP, TOTELE, NDGLNO, ierr)
          IF(IERR.NE.0) THEN
             FLAbort("halo 2 buggered")
          END IF
       endif
#endif
       CALL flcomms_tetra4_to_tetra10(halo_tag, NNODP, NDGLNO, halotagT10, &
            PGNDGLNO, TOTELE, PGNODSP)
#ifndef NDEBUG
       CALL flcomms_test_tetra4_to_tetra10(halotagT10, NDGLNO, PGNDGLNO, &
            TOTELE, PGNODSP, PGNODS, X, ierr)
       IF(IERR.LT.0) THEN
          FLAbort("halo generated for tetra10 is foobar")
       END IF
#endif

       ! Calculate the matrix sparcity...
       IF((GEOSHORT.NE.1).OR.(ELLBAL.EQ.1)) THEN
          CALL POSINM(pg_sparsity, TOTELE,PGNODS,MLOC, &
               pgndglno, pgnods, mloc, pgndglno, &
               name="QuadraticPressureSparsity", &
               nnodp=pgnodsp, halo_tag=halotagT10)
          call insert(state, pg_sparsity, name="QuadraticPressureSparsity")
          ! deallocate our reference:
          call deallocate(pg_sparsity)
       ENDIF
       MeshCount = CNT

    END IF ! end new mesh ******************

    pg_mesh=extract_mesh(state, "BalancePressureMesh")
    pg_shape=pg_mesh%shape

    pgndglno => pg_mesh%ndglno
    pgnods=node_count(pg_mesh)

    m = reshape(pg_shape%n, (/ mloc*ngi /))
    mlx = reshape(pg_shape%dn(:,:,X_), (/ mloc*ngi /))
    mly = reshape(pg_shape%dn(:,:,Y_), (/ mloc*ngi /))
    mlz = reshape(pg_shape%dn(:,:,Z_), (/ mloc*ngi /))

    IF((GEOSHORT.NE.1).OR.(ELLBAL.EQ.1)) THEN
       pg_sparsity=extract_csr_sparsity(state, "QuadraticPressureSparsity")
       pgfindrm => pg_sparsity%findrm
       pgcolm => pg_sparsity%colm
       pgcentrm => pg_sparsity%centrm
       npgcolm=size(pgcolm)
       allocate(difmath(1:npgcolm))
       if (associated(pg_sparsity%row_halo)) then
          pgnodsp=halo_nowned_nodes(pg_sparsity%row_halo)
       else
          pgnodsp=pgnods
       end if
    ENDIF

    ! **********initialisation of basic variables*******

    ALLOCATE(tempvec(PGNODS))

    ALLOCATE(DIFVECX(PGNODS))
    ALLOCATE(DIFVECY(PGNODS))
    ALLOCATE(DIFVECZ(PGNODS))
    ALLOCATE(DIFVECV(PGNODS))
    ALLOCATE(DIFVECN(PGNODS))
    ALLOCATE(DIFRHS(PGNODS))
    ALLOCATE(DIFRHSND(NONODS))
    ALLOCATE(EQETID(PGNODS))
    ALLOCATE(MLND(NONODS))
    ALLOCATE(MLNDVERT(NONODS))
    ALLOCATE(IWORKS(NONODS))
    ALLOCATE(UMEAN(NONODS))
    ALLOCATE(VMEAN(NONODS))
    ALLOCATE(WMEAN(NONODS))
    ALLOCATE(PNORMX(PGNODS))
    ALLOCATE(PNORMY(PGNODS))
    ALLOCATE(PNORMZ(PGNODS))
    ALLOCATE(DP(PGNODS))
    ALLOCATE(NDENST(NONODS))
    ALLOCATE(FREESTHETA(PGNODS))

    IF(GOTTID.EQ.1) THEN
       ! Got equalibrium tide EQETID - only the 11th component onwards counts...
       CALL INTEQETID(EQETID,X,Y,Z,ACCTIM, &
            XNONOD,PGNODS,NLOC,TOTELE,XONDGL, &
            ISPHERE,PGNDGLNO,MLOC,INT(FREESDISOTT/2**11), &
            TIDALLCOM,TLOVE,STATID)
       IF(GOT_EQTD) CALL VISPRE(EQTD,EQETID,NONODS,PGNODS, &
            NLOC,MLOC,TOTELE, &
            NDGLNO,PGNDGLNO )
    ELSE
       EQETID(1:PGNODS) = 0.0
    ENDIF

    ! Calculate normal to domain (PNORMX,PNORMY,PNORMZ)
    ewrite(2,*) 'gotboy:',gotboy
    IF(GOTBOY) THEN
       CALL FINDDIFVECN(NONODS,XNONOD,PGNODS,TOTELE,NLOC,MLOC, &
            NDGLNO,PGNDGLNO, &
            PNORMX,PNORMY,PNORMZ, &
            COGRAX,&            
            BSOUX,BSOUY,BSOUZ, &
            VGRAVX,VGRAVY,VGRAVZ, &
            ISPHERE,X,Y,Z )
    ELSE
       PNORMX(1:PGNODS) = 0.0
       PNORMY(1:PGNODS) = 0.0
       PNORMZ(1:PGNODS) = 0.0
    ENDIF

    ! Solve for geostrophic pressure options... ************
    IF(IsParallel()) THEN
       PARA = 1
    ELSE
       PARA=0
       PGNODSP = PGNODS
    END IF

    if (.not. new_options) then ! done above for new options
       if (wetdry==1) then
          rnonl_coefficient=1.0
       else
          rnonl_coefficient=0.01
       end if
    end if

    IF(SOLSTAGE.NE.2) THEN

       IF(NEWMES2.AND.(HYDSTART.EQ.0)) THEN
          DEALLOCATE(PG)
          DEALLOCATE(FREESOLD)
          DEALLOCATE(PGRAVX)
          DEALLOCATE(PGRAVY)
          DEALLOCATE(PGRAVZ)
          DEALLOCATE(POTGRA)
          DEALLOCATE(DGP)
          DEALLOCATE(ELEORD)
          DEALLOCATE(COLELE4)
       ENDIF
       GETELEORD=.FALSE.
       IF(NEWMES2.OR.(HYDSTART.EQ.1)) THEN
          ALLOCATE(PG(PGNODS))
          ALLOCATE(FREESOLD(PGNODS))
          ALLOCATE(PGRAVX(PGNODS))
          ALLOCATE(PGRAVY(PGNODS))
          ALLOCATE(PGRAVZ(PGNODS))
          ALLOCATE(POTGRA(PGNODS*USEPOTGRA))
          ALLOCATE(DGP(MLOC*TOTELE*IDGBAL))
          ALLOCATE(ELEORD(TOTELE*IDGBAL))
          ALLOCATE(COLELE4(TOTELE*5))

          PG(1:PGNODS) = 0.0
          FREESOLD(1:PGNODS) = 0.0
          DGP(1:MLOC*TOTELE*IDGBAL) = 0.0

          if (idgbal==1) then
             geteleord=.true.
          end if

          call allocate(free_surface, pg_mesh, "QuadraticFreeSurface")
          call insert(state, free_surface, "QuadraticFreeSurface")
          call deallocate(free_surface)

          ! calculate COLELE4 points to elements surrounding an element      
          IF(SNLOC.EQ.3) CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE4,NONODS)

          ! Calculate gravity direction consitent with mesh and constant density
          IF(ADJGRAV.EQ.0) THEN
             pgravx(1:PGNODS)=-pnormx(1:PGNODS)
             pgravy(1:PGNODS)=-pnormy(1:PGNODS)
             pgravz(1:PGNODS)=-pnormz(1:PGNODS)
          ELSE
             ewrite(1,*) 'JUST GOING INTO MAPGRAV ******'
             CALL MAPGRAV(PGRAVX,PGRAVY,PGRAVZ,ADJGRAV, &
                  POTGRA,USEPOTGRA,  &
                  PNORMX,PNORMY,PNORMZ, &
                  FINDRM,COLM,NCOLM,CENTRM, &
                  PGFINDRM,PGCOLM,NPGCOLM,PGNODS,PGCENTRM,  &
                  X,Y,Z, &
                  m,mlx,mly,mlz, &
                  N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
                  TOTELE,D3,DCYL,  &
                  NDGLNO,PGNDGLNO, XONDGL, &
                  NONODS,XNONOD, &
                  halotagT10,PARA,ProcNo,PGNODSP, &
                  GETFREES,NLEVEL,ISPHERE, &
                  option_path)
             ewrite(1,*) 'OUTOF MAPGRAV'
          ENDIF

          HYDSTART=0
       ENDIF

       ! Calculate NDENST from the source SOURCX,SOURCY,SOURCZ
       ! and PGRAVX,PGRAVY,PGRAVZ (mesh adjusted gravity direction)
       CALL CALNDENST(NDENST,SOURCX,SOURCY,SOURCZ, &
            PGRAVX,PGRAVY,PGRAVZ, &
            TOTELE,NLOC,MLOC,NONODS,PGNODS, &
            NDGLNO,PGNDGLNO)
    ENDIF

    free_surface=extract_scalar_field(state, "QuadraticFreeSurface")
    frees => free_surface%val

    IF(SOLSTAGE/=2 .and. (NEWMES2.OR.(HYDSTART.EQ.1)) .and. GETFREES) THEN
       CALL GESPSI(PSIPRE,FREES,NONODS,PGNODS, &
            NLOC,MLOC,TOTELE, &
            NDGLNO,PGNDGLNO)

       FREESOLD(1:PGNODS) = FREES(1:PGNODS)
       IF(MVMESH) THEN
          ! calculate the initial X,Y,Z by perturbing X,Y,Z 
          ! by free surface height FREES.
          call CalculateNewVerticalCoordinate(state, &
               original_coordinates, free_surface, &
               cograx)

          XOLD(1:XNONOD) = X(1:XNONOD)
          YOLD(1:XNONOD) = Y(1:XNONOD)
          ZOLD(1:XNONOD) = Z(1:XNONOD)
       ENDIF
    ENDIF

    ! Next time step...
    IF(SOLSTAGE.NE.2) THEN
       IF(GETFREES) THEN
          IF(NEXTTIM) THEN
             FREESOLD(1:PGNODS) = FREES(1:PGNODS)
          ENDIF
       ENDIF
    ENDIF

    allocate(orig_XABSOR(NONODS))
    allocate(orig_yABSOR(NONODS))
    allocate(orig_zABSOR(NONODS))

    orig_XABSOR=0.0
    orig_yABSOR=0.0
    orig_zABSOR=0.0

    IF(SOLSTAGE.NE.1) THEN
       ! Solve for free surface.
       IF(GETFREES) THEN

          CTHETA=0.0

          UMEAN(1:NONODS) = U(1:NONODS)
          VMEAN(1:NONODS) = V(1:NONODS)
          WMEAN(1:NONODS) = W(1:NONODS)

          ewrite(2,*) 'r2norm(EQETID,pgnods,0):',r2norm(EQETID,pgnods,0)

             if(shared_p_fs) then
! SHORTDP is actual free surface height for shared_p_fs
               ALLOCATE(SHORTDP(NONODS))
!                SHORTDP=psipre
                SHORTDP=p/gravty
!                SHORTDP=0.0
               if (nlevel==0) then
                   call vertical_extrapolation_wrapper( &
                        shortdp, velocity_mesh, &
                        shortdp, velocity_mesh, &
                        coordinates, cograx, topdis_field)
               else 
                FLAbort("Will not work for fixed levels.")
               end if
             ! Calculate DP from SHORTDP
               DO ELE=1,TOTELE
                  DO ILOC=1,NLOC
                     IGL=NDGLNO((ELE-1)*NLOC+ILOC)
                     DO JLOC=ILOC,NLOC
                        JGL=NDGLNO((ELE-1)*NLOC+JLOC)
                        KLOC=ILINK2(ILOC,JLOC)
                        PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                        DP(PNOD)=0.5*(SHORTDP(IGL)+SHORTDP(JGL))-FREES(PNOD)
                     END DO
                  END DO
               END DO
               DEALLOCATE(SHORTDP)
             else
          ! Assemble the 3D free surface enq for free surface height...
          IF(GEOSHORT.EQ.0) THEN

             if (.not.gottopdis) then
                FLAbort("Need distance to top and bottom for longalltopsufdg, whatever that means.")
             end if

             call longalltopsufdg(difmath,difvecv, &
                  xorig,yorig,zorig, &
                  topdis,botdis,atheta,xabsor, &
                  gravty,ctheta, &
                  nloc,&
                  totele, &
                  geobal,&
                  ndglno,xondgl, &
                  nonods,xnonod, &
                  scfacth0,isphere, &
                  pgfindrm,pgcentrm,pgcolm,npgcolm,pgnods, &
                  cograx, &
                  u,v,w, &
                  INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
                  nofsta,0,dt, &
                  mloc,pgndglno,frees,freesold, &
                  fstheta, rnonl_coefficient)

             ! apply free surface b.c's ********
             CALL LONGFREESFIMPBCS(NOBCFSimp,BCT1WFSimp,BCT2FSimp,DT, &
                  NONODS, &
                  PGFINDRM,PGCENTRM,PGCOLM,NPGCOLM,PGNODS, &
                  DIFMATH,DIFVECV, &
                  FREES, FREESOLD, &
                  TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO)

             ewrite(2,*) 'r2norm(difvecv,pgnods,0):',r2norm(difvecv,pgnods,0)

             DP(1:PGNODS) = 0.0

             ! lump free surface eqn in vertical and solve...
             ! solve *quadratic* free surface equation
             ewrite(1,*) 'going into freesurface MULGRIDVERT *************'
             CALL MULGRIDVERT(DIFMATH,DIFVECV, &
                  FINDRM,COLM,CENTRM,NCOLM,NONODS, &
                  PGFINDRM,PGCOLM,PGCENTRM,NPGCOLM,PGNODS, &
                  DP, &
                  TOTELE,NLOC,MLOC, &
                  NDGLNO,PGNDGLNO, X,Y,Z,XNONOD,ISPHERE, &
                  NLEVEL, &
                  PGNODSP, &
                  halotagT10,0, &
                  .FALSE.,.FALSE.,.FALSE.,.FALSE., &
                  option_path)
             if (nlevel==0) then
                call vertical_extrapolation_wrapper( &
                     dp, pg_mesh, &
                     dp, pg_mesh, &
                     coordinates, cograx, topdis_field)
             end if

          ELSE

             ALLOCATE(SHORTDIFMATH(NCOLM))
             ALLOCATE(SHORTDIFVECV(NONODS))
             ALLOCATE(SHORTDP(NONODS))
             ALLOCATE(RZER(NONODS))
             ALLOCATE(WETDRYMAT(NCOLM*WETDRY))
             ALLOCATE(SUFML(NONODS))

             shortdifmath=0.
             shortdifvecv=0.
             shortdp=0.
             rzer=0.
             wetdrymat=0.
             sufml=0.

             wetdry_loop: do nloop_wetdry=1, 1+wetdry*5

                IF(WETDRY.EQ.1) THEN
                   allocate(nodfrees(nonods))
                   DO ELE=1,TOTELE
                      DO ILOC=1,NLOC
                         KLOC=ILINK2(ILOC,ILOC)
                         PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                         NOD=XONDGL((ELE-1)*NLOC+ILOC)
                         NODFREES(NOD)=FREES(PNOD)
                      END DO
                   END DO
                   if (.not.gottopdis) then
                      FLAbort("Need distance to top and bottom for wetting and drying.")
                   end if
                   CALL ABSWETDRY(XABSOR,YABSOR,ZABSOR,NONODS,XNONOD, &
                        NCOLM,FINDRM,COLM,TOPDIS,BOTDIS,nodfrees,shortdp,DT,&
                        coef_abs,ee00,dd00,totele,nloc,ndglno, xondgl, wdgravity)
                   ! ADD THE ORIGINAL QUADRATIC DRAG ONTO ABSORPTION
                   XABSOR=XABSOR+ORIG_XABSOR
                   YABSOR=YABSOR+ORIG_YABSOR
                   ZABSOR=ZABSOR+ORIG_ZABSOR
                   deallocate(nodfrees)
                ENDIF

                if (.not.gottopdis) then
                   FLAbort("Need distance to top and bottom for shortpgassem1free.")
                end if

                call shortpgassem1free(shortdifmath, &
                     shortdifvecv, &
                     wetdrymat,sufml, &
                     topdis,botdis,atheta,xabsor,  &
                     frees,freesold,gravty,dt,ctheta,0, &
                     xorig,yorig,zorig, umean,vmean,wmean, &
                     INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
                     nloc,mloc, &
                     totele, &
                     geobal,&
                     ndglno,pgndglno, xondgl, &
                     pgnods,nonods,xnonod,scfacth0,isphere, &
                     findrm,colm,ncolm, &
                     findrm,centrm,colm,ncolm,fredop, &
                     nofsta, &
                     cograx, &
                     fstheta,rnonl_coefficient, &
                     wetdry,colele4,wdgravity)

                ! apply free surface b.c's ********
                CALL FREESFIMPBCS(NOBCFSimp,BCT1WFSimp,BCT2FSimp,DT, &
                     NONODS,FINDRM,CENTRM,COLM,NCOLM, &
                     SHORTDIFMATH,SHORTDIFVECV, &
                     FREES,FREESOLD,PGNODS, &
                     TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO)
                ! apply free surface b.c's ********

                ewrite(2,*) 'r2norm(SHORTDIFVECV,nonods,0):',r2norm(SHORTDIFVECV,nonods,0)

                ewrite(1,*) 'reduced length of free surface           *************'
                ewrite(1,*) 'going into freesurface short MULGRIDVERT *************'

                ! solve *linear* free surface equation
                CALL MULGRIDVERT(SHORTDIFMATH,SHORTDIFVECV, &
                     FINDRM,COLM,CENTRM,NCOLM,NONODS, &
                     FINDRM,COLM,CENTRM,NCOLM,NONODS, &
                     SHORTDP, &
                     0,0,0, &
                     IIDUM,IIDUM, RD,RD,RD,0,0, &
                     NLEVEL, &
                     NNODP, &
                     halo_tag,0, &
                     .FALSE.,.FALSE.,(abs(CTHETA).gt.0.0001),.FALSE., &
                     option_path)

                if (nlevel==0) then
                   call vertical_extrapolation_wrapper( &
                        shortdp, velocity_mesh, &
                        shortdp, velocity_mesh, &
                        coordinates, cograx, topdis_field)
                end if

                ewrite(2,*) 'r2norm(Shortdp,nonods,0):',r2norm(SHORTdp,nonods,0)

                !REMOVE ewrite(2,*) 'nloop_wetdry, mindepth==',nloop_wetdry, minval(topdis+botdis+shortdp)

             end do wetdry_loop

             ewrite(1,*) 'just finished solving for free surface'

             ! Calculate DP from SHORTDP
             DO ELE=1,TOTELE
                DO ILOC=1,NLOC
                   IGL=NDGLNO((ELE-1)*NLOC+ILOC)
                   DO JLOC=ILOC,NLOC
                      JGL=NDGLNO((ELE-1)*NLOC+JLOC)
                      KLOC=ILINK2(ILOC,JLOC)
                      PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                      DP(PNOD)=0.5*(SHORTDP(IGL)+SHORTDP(JGL))
                   END DO
                END DO
             END DO
             

             ewrite(2,*) 'min&max dp====', minval(dp), maxval(dp)
             ewrite(2,*) 'min&max shortdp====', minval(shortdp), maxval(shortdp)

             DEALLOCATE(SHORTDIFMATH)
             DEALLOCATE(SHORTDIFVECV)
             DEALLOCATE(SHORTDP)
             DEALLOCATE(RZER)
! endof IF(GEOSHORT.EQ.0) THEN ELSE...
          ENDIF
          
! endof if(shared_p_fs) then else...
                endif

          DO I=1,PGNODS
             difrhs(i)=dp(i)
             FREES(I)=FREES(I)+DP(I)
             FREESTHETA(I)=FSTHETA*FREES(I)+(1.-FSTHETA)*FREESOLD(I)
          END DO

          if (mvmesh) then
             call CalculateNewVerticalCoordinate(state, &
                  original_coordinates, free_surface, cograx)
          end if

          call CalculateTopBottomDistance(state, cograx)

          ! ENDOF IF(GETFREES) THEN...
       ELSE
          if (allocated(freestheta) .and. size(freestheta) /= PGNODS) then
             deallocate(freestheta)
          end if
          if (.not. allocated(freestheta)) then
             allocate(freestheta(PGNODS))
          end if
          freestheta = 0.0
       ENDIF

       ewrite(2,*)  'wetdry-here===',wetdry
       ewrite(2,*)  'dd00===', dd00

       if(wetdry.eq.1) then

          if (.not.gottopdis) then
             FLAbort("Need distance to top and bottom for wetting and drying.")
          else if (xnonod/=nonods) then
             FLAbort("Wetting and drying doesn't work yet with periodic domains.")
          end if
          if (nlevel==0) then
             FLAbort("Wetting and drying doesn't work yet with nlevel==0")
          end if

          allocate(nodfrees(nonods))

          call alterfrees(wetdrymat,sufml,topdis,botdis,findrm, colm, &
               centrm, ncolm, xorig, yorig, zorig, nodfrees, &
               nonods, nlevel, dd00)

          do ele=1,totele
             do iloc=1,nloc
                inod=xondgl((ele-1)*nloc+iloc)
                do jloc=iloc,nloc
                   kloc=ilink2(iloc,jloc)
                   pnod=pgndglno((ele-1)*mloc+kloc)
                   jnod=xondgl((ele-1)*nloc+jloc)
                   frees(pnod)=0.5*(nodfrees(inod)+nodfrees(jnod)) 
                end do
             end do
          end do

          call CalculateNewVerticalCoordinate(state, &
               original_coordinates, free_surface, &
               cograx, minimum_depth=dd00)

          allocate(shortdp(nonods))
          shortdp=0.

          call abswetdry(xabsor,yabsor,zabsor,nonods,xnonod, &
               ncolm,findrm,colm,topdis,botdis,nodfrees,shortdp,DT,&
               coef_abs,ee00,dd00,totele,nloc,ndglno, xondgl, wdgravity)

          ! ADD THE ORIGINAL QUADRATIC DRAG ONTO ABSORPTION
          XABSOR=XABSOR+ORIG_XABSOR
          YABSOR=YABSOR+ORIG_YABSOR
          ZABSOR=ZABSOR+ORIG_ZABSOR

          deallocate(nodfrees)
          deallocate(shortdp)

       end if

       !       IF(WETDRY.EQ.1) THEN
       ! Hedong insert free surface amendment to ensure depth of 0.05m say. 
       ! WETDRYMAT is the matrix and SUFML is the surface lumped mass matrix.
       !       CALL ALTERFREES
       ! Hedong insert code that calculates the absorption. 
       !        CALL ALTERABSORB
       !       ENDIF

       ! ENDOF IF(SOLSTAGE.NE.1) THEN
    ENDIF

    ! Solve for geostrophic pressure options... ************

    IF(SOLSTAGE.NE.2) THEN
       ! Not solving for just free surface...

       IF(IDGBAL.EQ.1) THEN

          ! *******************************************************
          ! ******* DG pressure solver ****************************
          ! *******************************************************
          
          IF((TOTELE.LT.0).AND.(NLEVEL==0)) THEN
             ! Store the DG matrices for speed (only needed for unstructured meshes)...
             CALL VERTICALDGSTORE(VECX,VECY,VECZ, &
                  NSUBVLOC,VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
                  DGP,ELEORD,GETELEORD,NLEVEL, &
                  GETFREES.AND.(.NOT.SHARED_P_FS),FREESTHETA,EQETID,GRAVTY,  &
                  PNORMX,PNORMY,PNORMZ, &
                  PGRAVX,PGRAVY,PGRAVZ,NDENST, &
                  POTGRA,USEPOTGRA,  &
                  PGNODS, &
                  X,Y,Z, U,V,W, &
                  INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
                  NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,COGRAX, &
                  m,mlx,mly,mlz,&
                  N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
                  TOTELE, &
                  GEOBAL,&
                  NDGLNO,PGNDGLNO, XONDGL, &
                  NONODS,XNONOD, &
                  scfacth0,ISPHERE, &
                  FINDRM,COLM,NCOLM, &
                  IDGBCMOM,IDGBOYONLY,SMOOTH.AND.((GEOBAL.LT.-27).OR.(GEOBAL.GT.-25)), &
                  COLELE4,PHSOLVE.AND.((GEOBAL.NE.-26).AND.(GEOBAL.NE.-25)),wdgravity)
!                  IDGBCMOM,IDGBOYONLY,SMOOTH,COLELE4,PHSOLVE,wdgravity)
          ELSE
             CALL VERTICALDG(VECX,VECY,VECZ, &
                  NSUBVLOC,VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
                  DGP,ELEORD,GETELEORD,NLEVEL, &
                  GETFREES.AND.(.NOT.SHARED_P_FS),FREESTHETA,EQETID,GRAVTY,  &
                  PNORMX,PNORMY,PNORMZ, &
                  PGRAVX,PGRAVY,PGRAVZ,NDENST, &
                  POTGRA,USEPOTGRA,  &
                  PGNODS, &
                  X,Y,Z, U,V,W, &
                  INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
                  NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,COGRAX, &
                  m,mlx,mly,mlz,&
                  N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
                  TOTELE, &
                  GEOBAL,&
                  NDGLNO,PGNDGLNO, XONDGL, &
                  NONODS,XNONOD, &
                  scfacth0,ISPHERE, &
                  FINDRM,COLM,NCOLM, &
                  IDGBCMOM,IDGBOYONLY,SMOOTH.AND.((GEOBAL.LT.-27).OR.(GEOBAL.GT.-25)), &
!                  IDGBCMOM,IDGBOYONLY,.true., &
                  COLELE4,PHSOLVE.AND.((GEOBAL.NE.-26).AND.(GEOBAL.NE.-25)),wdgravity)
!                  COLELE4,.true.,wdgravity)
!                  IDGBCMOM,IDGBOYONLY,SMOOTH,COLELE4,PHSOLVE,wdgravity )
!                  IDGBCMOM,IDGBOYONLY,.false.,COLELE4,PHSOLVE,wdgravity)
!                  IDGBCMOM,IDGBOYONLY,.false.,COLELE4,.false.,wdgravity) 
          ENDIF
       ENDIF

       IF(ELLBAL.EQ.1) THEN
          ! *******************************************************
          ! ******* elliptic pressure solver **********************
          ! *******************************************************

          ! Take out the vertical presure - solve for it seperatly above...
          TAKVERT=(IDGBAL.EQ.1)
          ! Take out the vertical free surface gradient.
          FSTAKVERT=.true.

          ! Assemble the elliptic enq for geostophic pressure
          ewrite(1,*) 'just going into PGASSEM1P---'
          CALL PGASSEM1P(DIFMATH, &
               DIFVECX,DIFVECY,DIFVECZ, &
               FREESTHETA,EQETID,GETFREES,GRAVTY, &
               TAKVERT,FSTAKVERT, .FALSE.,PNORMX,PNORMY,PNORMZ,  &
               PGRAVX,PGRAVY,PGRAVZ,NDENST, &
               POTGRA,USEPOTGRA,  &
               PGFINDRM,PGCOLM,NPGCOLM,PGNODS,  &
               X,Y,Z, U,V,W, &
               m,mlx,mly,mlz, &
               N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
               TOTELE,D3,DCYL, &
               GEOBAL,&
               NDGLNO,PGNDGLNO, XONDGL, &
               NONODS,XNONOD, &
               scfacth0,ISPHERE, &
               IDGBAL,ELLBOYONLY,ELLCORONLY )
          ewrite(1,*) '********finished PGASSEM1P-'

          CALL PMINMX(DIFVECX,PGNODS,'******DIFVECX  ')
          CALL PMINMX(DIFVECY,PGNODS,'******DIFVECY  ')
          CALL PMINMX(DIFVECZ,PGNODS,'******DIFVECZ  ')

          ! Calculate DIFVECN the component in the virtical direction ****
          DO I=1,PGNODS
             DIFVECV(I)=DIFVECX(I)+DIFVECY(I)+DIFVECZ(I)
             DIFVECN(I)=DIFVECV(I)
          END DO

          ewrite(2,*) '*****************NLEVEL=',NLEVEL

          ! The pressure level to set to 0.0
          IF((PARA.EQ.0).OR.(ProcNo.EQ.1)) THEN
             PPSET=1
          ELSE
             PPSET=0
          ENDIF

          IF(GETFREES.AND.(ELETOPBC.EQ.1)) THEN
             ! Apply b.c's to the surface
             if (nlevel==0) then
                FLAbort("ELETOPBC does not yet work for non-layered meshes")
             end if
             IF(ELLSTAGTOPBC.EQ.1) THEN
                ! Just specify the corner nodes on the surface...
                CALL TOPBCDIFDIA(DIFMATH, &
                     PGFINDRM,PGCOLM,PGCENTRM,NPGCOLM,PGNODS,NONODS, &
                     TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO,DIFVECV,NLEVEL,.TRUE.)
             ELSE
                ! Apply b.c. to all nodes on the surface...
                CALL TOPBCDIFDIA(DIFMATH, &
                     PGFINDRM,PGCOLM,PGCENTRM,NPGCOLM,PGNODS,NONODS, &
                     TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO,DIFVECV,NLEVEL,.FALSE.)
             ENDIF
             PPSET=0
          ELSE
             ! Apply b.c=0 at PNOD=PPSET
             ewrite(2,*) 'set one value only'
             CALL DISSETMAVAL(DIFMATH,NPGCOLM,PGNODS,PGFINDRM,PGCOLM,PGCENTRM, &
                  PPSET,DIFVECV)
          ENDIF

          IF(NLEVEL.LT.2) THEN
             ONLY3D=.TRUE.
          ELSE
             ONLY3D=.FALSE.
          ENDIF

          pg_sparsity=extract_csr_sparsity(state, "QuadraticPressureSparsity")
          pg_matrix=wrap(pg_sparsity, val=difmath, name="QuadraticPressureMatrix")
          if (have_option(trim(option_path)// &
             "/prognostic/solver/preconditioner::mg/vertical_lumping")) then
               
             if (.not. has_csr_matrix(state, "QuadraticVerticalProlongationOperator")) then
                prolongator=vertical_prolongator_from_free_surface(state, pg_mesh)
                call insert(state, prolongator, name="QuadraticVerticalProlongationOperator")
                call deallocate(prolongator)
             end if
             
             prolongator = extract_csr_matrix(state, &
               "QuadraticVerticalProlongationOperator")
             call petsc_solve(pg, pg_matrix, difvecv, &
                option_path=option_path, &
                prolongator=prolongator)
          else
             call petsc_solve(pg, pg_matrix, difvecv, &
                option_path=option_path)
          end if
          call deallocate(pg_matrix)
          
          ! this next line is for the test case to check the geopressure norms
          ewrite(1,*) 'OUT OF MULGRIDVERT'
          ewrite(2,'(a,1pe11.4)') 'top-level-geopressure_norm: ',r2norm(pg,pgnods,0)

          DEALLOCATE(DIFVECX)
          DEALLOCATE(DIFVECY)
          DEALLOCATE(DIFVECZ)
          DEALLOCATE(DIFVECV)
          DEALLOCATE(DIFVECN)
          DEALLOCATE(DIFRHS)
          DEALLOCATE(DIFRHSND)
          DEALLOCATE(MLND)
          DEALLOCATE(MLNDVERT)
          DEALLOCATE(IWORKS)
          DEALLOCATE(UMEAN)
          DEALLOCATE(VMEAN)
          DEALLOCATE(WMEAN)
          DEALLOCATE(DP)

          deallocate(wdgravity)

          ! Diagnostic output of geostrophic pressure.
          ewrite_minmax(PG)

          ! Put the result into momentum eqns
          CALL PGMOME1P(VECX,VECY,VECZ, PG, &
               FREESTHETA,EQETID,GETFREES,GRAVTY, &
               TAKVERT,FSTAKVERT, PNORMX,PNORMY,PNORMZ,  &
               PGRAVX,PGRAVY,PGRAVZ,NDENST, &
               POTGRA,USEPOTGRA,  &
               SOURCX,SOURCY,SOURCZ,&
               X,Y,Z, U,V,W, &
               m,mlx,mly,mlz, &
               N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
               TOTELE,D3,DCYL, &
               GEOBAL,&
               NDGLNO,PGNDGLNO, XONDGL, &
               NONODS,PGNODS,XNONOD, &
               scfacth0,ISPHERE, &
               IDGBAL,ELLBOYONLY,ELLCORONLY )

          ewrite(1,*) "just out of PGMOME"

          DEALLOCATE(PNORMX)
          DEALLOCATE(PNORMY)
          DEALLOCATE(PNORMZ)

          CALL PMINMX(PG,PGNODS,'******PG  ')
          CALL PMINMX(VECX,NONODS,'******VECX  ')

          ! *********************************************************
          ! ENDOF IF(DGBAL.EQ.1) THEN ELSE...
       ENDIF

       ! its nec. to set these to zero for non-free surface flows
       SOUVECX=0.0
       SOUVECY=0.0
       SOUVECZ=0.0

       ! ENDOF IF(SOLSTAGE.NE.2) THEN...
    ENDIF

    ! copy the 1st NONODS of PG into PSIPRE
    IF(GETFREES) THEN

       if (.not.gottopdis) then
          FLAbort("Need distance to top and bottom for free surface.")
       end if

       IF(SOLSTAGE.NE.1) THEN

          IF(MVMESH) THEN
             ! Update the coordinates X,Y,Z:
             call CalculateNewVerticalCoordinate(state, &
                  original_coordinates, free_surface, &
                  cograx)
             ! Calculate the grid velocity UG,VG,WG
             call calug(x,y,z, xold,yold,zold, ug,vg,wg, &
               totele, nloc, ndglno, xondgl, dt)
             CALL PMINMX(Z,XNONOD,'******Z  ')
             CALL PMINMX(ZOLD,XNONOD,'******ZOLD  ')
             CALL PMINMX(ZORIG,XNONOD,'******ZORIG  ')
             CALL PMINMX(FREES,PGNODS,'******FREES  ')
             CALL PMINMX(UG,NONODS,'******UG  ')
             CALL PMINMX(VG,NONODS,'******VG  ')
             CALL PMINMX(WG,NONODS,'******WG  ')

             ! Recalculate TOPDIS & BOTDIS
             call CalculateTopBottomDistance(state, cograx)

             ewrite(2,*)  'max&min depth--->', maxval(topdis+botdis), minval(topdis+botdis)
             ewrite(1,*)  'just after CalculateTopBottomDistance'
          ENDIF

          ! IF(SOLSTAGE.NE.1) THEN
       ENDIF

       IF(SOLSTAGE.NE.2) THEN
          !     This subroutine caluculates C1T,C2T,C3T by NOT integrating 
          !     the pressure term by parts and integrating over the 
          !     top surface of the ocean.  

          CALL CTDGSUP(C1T,C2T,C3T, &
               XORIG,YORIG,ZORIG, &
               N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
               TOTELE,D3, &
               NDGLNO,XONDGL, &
               NONODS,XNONOD, &
               ISPHERE, &
               FINDCT,COLCT,NCT, &
               COGRAX, &
               NOBCFSimp,BCT2FSimp, &
               PGNODS,MLOC,FREES,PGNDGLNO, &
               TOPDIS,BOTDIS,COLELE4)

          ! ENDOF IF(SOLSTAGE.NE.2) THEN...
       ENDIF

       ewrite(2,*) 'SOLSTAGE,NSUBVLOC=',SOLSTAGE,NSUBVLOC
       ewrite(2,*) 'ISPHERE,DISOPT:',ISPHERE,DISOPT
       ewrite(2,*) '(((ISPHERE.NE.2).OR.(DISOPT.LE.148)).OR.(NSUBVLOC.EQ.0))', &
            (((ISPHERE.NE.2).OR.(DISOPT.LE.148)).OR.(NSUBVLOC.EQ.0))

       ! Calculate the source that corresponds to Galerkin discretised VECX etc.
       IF(SOLSTAGE.NE.2) THEN
          IF(((ISPHERE.NE.2).OR.(DISOPT.LE.148)).OR.(NSUBVLOC.NE.0)) THEN
             SOUVECX=0.0
             SOUVECY=0.0
             SOUVECZ=0.0
          ELSE
             ! If .not.Galerkin just place in SOUVECX,SOUVECY,SOUVECZ else use VECX

             ewrite(3,*) 'before VECXYZ2SOURCXYZ('
             CALL VECXYZ2SOURCXYZ(SOUVECX,SOUVECY,SOUVECZ, VECX,VECY,VECZ, &
                  X,Y,Z, NDGLNO,XONDGL, TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI, &
                  N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT,&
                  ISPHERE,FINDRM,COLM,NCOLM, &
                  halo_tag,NNODP)
          ENDIF
          ! ENDOF IF(SOLSTAGE.NE.2) THEN...
       ENDIF

       CALL VISPRE(PSIPRE,FREES,NONODS,PGNODS, &
            NLOC,MLOC,TOTELE, &
            NDGLNO,PGNDGLNO)

       ewrite(2,*) 'max,min---psipre', maxval(psipre),minval(psipre)
       ewrite(2,*) 'max,min---topdis', maxval(topdis),minval(topdis)
       ewrite(2,*) 'max,min---botdis', maxval(botdis),minval(botdis)

    ELSE
       IF(IDGBAL.EQ.1) THEN
          CALL DG2CTYNOD(PSIPRE,DGP,TOTELE,NLOC,MLOC,NONODS,NDGLNO)
       ELSE
          CALL VISPRE(PSIPRE,PG,NONODS,PGNODS, &
               NLOC,MLOC,TOTELE, &
               NDGLNO,PGNDGLNO)
       ENDIF
    ENDIF

!!!Output PG for reduced model
    if(have_option("/model/fluids/reduced")) then
       CALL IN_OUTPUT_PG(0,PG,PGNODS,TOTELE,MLOC,NGI,PGNDGLNO, &
            M,MLX,MLY,MLZ,         &
            D3,1, &
            NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM)
    endif

    IF((GEOSHORT.NE.1).OR.(ELLBAL.EQ.1)) THEN
       deallocate(difmath)
    end if

  END SUBROUTINE GEOELI1P

  SUBROUTINE DISSETMAVAL(DIFMAT,NPGCOLM,PGNODS,PGFINDRM,PGCOLM,PGCENTRM, &
       PNOD,VECRHS)
    IMPLICIT NONE
    INTEGER NPGCOLM,PGNODS,PNOD
    INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM),PGCENTRM(PGNODS)
    REAL DIFMAT(NPGCOLM),VECRHS(PGNODS)
    ! Local variables...
    INTEGER COUNT,PNOD2
    REAL SUM
    INTEGER, ALLOCATABLE, DIMENSION(:)::ISET
    REAL, ALLOCATABLE, DIMENSION(:)::DIAG

    IF(PNOD.GE.1) THEN
       ALLOCATE(ISET(PGNODS))
       ALLOCATE(DIAG(PGNODS))
       ISET(1:PGNODS) = 0
       DIAG(1:PGNODS) = 0.0
       ISET(PNOD)=1
       SUM=0.0
       DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
          SUM=SUM+ABS(DIFMAT(COUNT))
          DIFMAT(COUNT)=0.0
       END DO
       DIAG(PNOD)=SUM
       VECRHS(PNOD)=0.0

       DO PNOD2=1,PGNODS
          DO COUNT=PGFINDRM(PNOD2),PGFINDRM(PNOD2+1)-1
             IF(ISET(PGCOLM(COUNT)).EQ.1) DIFMAT(COUNT)=0.0
          END DO
       END DO
       DIFMAT(PGCENTRM(PNOD))=DIFMAT(PGCENTRM(PNOD))+DIAG(PNOD)
    ENDIF

  END SUBROUTINE DISSETMAVAL

  SUBROUTINE NDTOPOINT(IWORK, TOTELE,NLOC,NONODS,PGNODS,NCOLM, &
       NDGLNO,CENTRM,FINDRM,COLM, &
       nlevel)
    ! Calculate IWORK(PNOD)=-ve the node at the top of the domain and
    ! +ve if a node at the top of the domain.
    ! ***************************************************************
    IMPLICIT NONE
    INTEGER TOTELE,NLOC,NONODS,PGNODS,NCOLM
    INTEGER NDGLNO(TOTELE*NLOC),IWORK(PGNODS)
    INTEGER CENTRM(NONODS),FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER NLEVEL
    ! Local variables...
    INTEGER I,II,ITOPND,PNODS,ELE,ILOC,JLOC,JCOLTOP
    INTEGER INOD,JNOD,COUNT,JCOL,COUNT2
    INTEGER INODTOP,PNOD,IFOUND
    INTEGER, ALLOCATABLE, DIMENSION(:)::IWORKS
    INTEGER, ALLOCATABLE, DIMENSION(:)::DONMAT
    REAL, ALLOCATABLE, DIMENSION(:)::XP
    REAL, ALLOCATABLE, DIMENSION(:)::YP
    REAL, ALLOCATABLE, DIMENSION(:)::ZP
    ALLOCATE(XP(PGNODS))
    ALLOCATE(YP(PGNODS))
    ALLOCATE(ZP(PGNODS))
    ! Calculate IWORKS(NOD)=top surface node -similar to IWORK but for corner nodes.
    ALLOCATE(IWORKS(NONODS))
    IWORKS(1:NONODS) = 0
    ! Mark the nodes at the top of the domain
    DO I=1,NONODS
       II=MOD(I-1,NONODS/NLEVEL)+1
       ITOPND= II +(NLEVEL-1)*(NONODS/NLEVEL)
       IF(I.EQ.ITOPND) THEN
          IWORKS(I)=ITOPND
       ELSE
          IWORKS(I)=-ITOPND
       ENDIF
    END DO

    IF(NLOC.EQ.0) THEN
       CALL ICOPY(IWORK,IWORKS,NONODS)
       IF(NONODS.NE.PGNODS) THEN
          ewrite(0,*) 'THINGS HAVE GONE WRONG NONODS,PGNODS:',NONODS,PGNODS
          FLAbort("STOP 3945")
       ENDIF
    ENDIF

    IF(NLOC.NE.0) THEN
       ! Calculate DONMAT a matrix containing the PNODS
       ALLOCATE(DONMAT(NCOLM))
       DONMAT(1:NCOLM) = 0
       PNODS=0
       DO ELE=1,TOTELE
          DO ILOC=1,NLOC
             DO JLOC=1,NLOC
                INOD=NDGLNO((ELE-1)*NLOC+ILOC)
                JNOD=NDGLNO((ELE-1)*NLOC+JLOC)
                ! Find JNOD...(cONSIDER ONLY UPPER DIAGONAL PART OF DONMAT)
                DO COUNT=CENTRM(INOD),FINDRM(INOD+1)-1
                   JCOL=COLM(COUNT)
                   IF(JCOL.EQ.JNOD) THEN
                      IF(DONMAT(COUNT).EQ.0) THEN
                         PNODS=PNODS+1
                         DONMAT(COUNT)=PNODS
                      ENDIF
                   ENDIF
                ENDDO
             END DO
          END DO
       END DO

       ! Make DONMAT matrix symmetric...
       DO INOD=1,NONODS
          DO COUNT=CENTRM(INOD),FINDRM(INOD+1)-1
             JCOL=COLM(COUNT)
             ! Find the other COUNT2 for symmetry
             DO COUNT2=FINDRM(JCOL),CENTRM(JCOL)-1
                IF(COLM(COUNT2).EQ.INOD) THEN
                   IF(DONMAT(COUNT2).EQ.0) THEN
                      DONMAT(COUNT2)=DONMAT(COUNT)
                   ELSE
                      ewrite(0,*) 'DONMAT IS NOT ZERO PROBLEM'
                      FLAbort("STOP 3831")
                   ENDIF
                ENDIF
             END DO
          END DO
       END DO

       ! Calculate IWORK(PNOD) from DONMAT
       DO INOD=1,NONODS
          INODTOP=ABS(IWORKS(INOD))
          DO COUNT=CENTRM(INOD),FINDRM(INOD+1)-1
             JCOL=COLM(COUNT)
             JCOLTOP=ABS(IWORKS(JCOL))
             PNOD=DONMAT(COUNT)
             IFOUND=0
             ! On the surface calculate the node
             DO COUNT2=FINDRM(INODTOP),FINDRM(INODTOP+1)-1
                IF(COLM(COUNT2).EQ.JCOLTOP) THEN
                   IF((IWORKS(INOD).GT.0).AND.(IWORKS(JCOL).GT.0)) THEN
                      IWORK(PNOD)=DONMAT(COUNT2)
                   ELSE
                      IWORK(PNOD)=-DONMAT(COUNT2)
                   ENDIF
                   IFOUND=1
                ENDIF
             END DO
             IF(IFOUND.EQ.0) THEN
                ewrite(0,*) 'IFOUND,PNOD,INODTOP:',IFOUND,PNOD,INODTOP
                FLAbort("stop 383")
             ENDIF
          END DO
       END DO
       DEALLOCATE(DONMAT)
       ! ENDOF IF(NLOC.NE.0) THEN
    ENDIF

    DEALLOCATE(IWORKS)
    DEALLOCATE(XP)
    DEALLOCATE(YP)
    DEALLOCATE(ZP)

  END SUBROUTINE NDTOPOINT

  SUBROUTINE DIFMATLUMP(DIFMATH,DIFVECZ,IWORK,PGNODS, &
       PGFINDRM,PGCOLM,NPGCOLM, &
       SDIFMATH,SDIFVECZ,SGFINDRM,SGCOLM,SGCENTRM,NSGCOLM,SGNODS,P2SMAP)
    ! Now lump DIFMATH & DIFVECZ so that is only a horizontal system
    ! and put the result in SDIFMATH,SDIFVECV...
    ! IWORK is needed in this sub.
    INTEGER PGNODS,NPGCOLM,NSGCOLM,SGNODS
    REAL DIFMATH(NPGCOLM),DIFVECZ(PGNODS)
    REAL SDIFMATH(NSGCOLM),SDIFVECZ(SGNODS)
    INTEGER IWORK(PGNODS)
    INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    INTEGER SGFINDRM(SGNODS+1),SGCOLM(NSGCOLM),SGCENTRM(SGNODS)
    INTEGER P2SMAP(PGNODS)
    ! Local variables...
    INTEGER PNOD,PNODTOP,COUNT,PCOL,PCOLTOP,COUNT2,SNOD,SCOL
    INTEGER SNODTOP,SCOLTOP

    ewrite(1,*) 'In difmatlump'
    ! calculate SGFINDRM,SGCOLM,SGCENTRM
    SNOD=0
    COUNT2=0
    DO PNOD=1,PGNODS
       IF(IWORK(PNOD).EQ.PNOD) THEN
          SNOD=SNOD+1
          SGFINDRM(SNOD)=COUNT2+1
          DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
             PCOL=PGCOLM(COUNT)
             IF(IWORK(PCOL).EQ.PCOL) THEN
                COUNT2=COUNT2+1
                SCOL=P2SMAP(IWORK(PCOL))
                SGCOLM(COUNT2)=SCOL
             ENDIF
          END DO
          ! perform a bubble sort on the row to put in assending order...
          CALL IBUBLE(SGCOLM(SGFINDRM(SNOD)),COUNT2+1-SGFINDRM(SNOD))
       ENDIF
    END DO
    SGFINDRM(SGNODS+1)=COUNT2+1
    DO SNOD=1,SGNODS
       DO COUNT=SGFINDRM(SNOD),SGFINDRM(SNOD+1)-1
          IF(SNOD.EQ.SGCOLM(COUNT)) SGCENTRM(SNOD)=COUNT
       END DO
    END DO

    ! calculate SDIFMATH
    SDIFMATH(1:NSGCOLM) = 0.0
    DO PNOD=1,PGNODS
       PNODTOP=ABS(IWORK(PNOD))
       SNODTOP=P2SMAP(PNODTOP)
       DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
          ! Put this contribution into matrix SDIFMATH
          PCOL=PGCOLM(COUNT)
          PCOLTOP=ABS(IWORK(PCOL))
          SCOLTOP=P2SMAP(PCOLTOP)
          DO COUNT2=SGFINDRM(SNODTOP),SGFINDRM(SNODTOP+1)-1
             IF(SGCOLM(COUNT2).EQ.SCOLTOP) THEN
                SDIFMATH(COUNT2)=SDIFMATH(COUNT2)+DIFMATH(COUNT)
             ENDIF
          END DO
       END DO
    END DO

    ! Lump the vector and put unity on diagonal of matrix.
    SDIFVECZ(1:SGNODS) = 0.0
    DO PNOD=1,PGNODS
       PNODTOP=ABS(IWORK(PNOD))
       SNODTOP=P2SMAP(PNODTOP)
       SDIFVECZ(SNODTOP)=SDIFVECZ(SNODTOP)+DIFVECZ(PNOD)
    END DO


  END SUBROUTINE DIFMATLUMP

  SUBROUTINE DIFSOLVCOLNS(DIFMATH,DP,DIFVECN,IWORK,P2SMAP,PGNODS,NONODS, &
       PGFINDRM,PGCOLM,NPGCOLM,NLEVEL, &
       TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO, &
       X,Y,Z,XNONOD,ISPHERE)
    ! This sub solves for vertical colns for DP
    ! That is MATT*DP=DIFVECN.
    ! IWORK is needed in this sub.
    INTEGER PGNODS,NONODS,XNONOD,NPGCOLM,NLEVEL
    REAL DP(PGNODS),DIFVECN(PGNODS)
    REAL DIFMATH(NPGCOLM)
    INTEGER IWORK(PGNODS),P2SMAP(PGNODS)
    INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    INTEGER TOTELE,NLOC,MLOC,NDGLNO(NLOC*TOTELE),PGNDGLNO(MLOC*TOTELE)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER ISPHERE
    ! Local variables...
    REAL RINT,LX,COEF
    INTEGER PNOD,COUNT,PCOL,NOD,IW
    INTEGER INLEV
    INTEGER ELE,ILOC,JLOC,KLOC,INOD,JNOD,NPLEVEL,LOCNOD
    INTEGER CINLEV,ISOLV,IPLEVEL,JPLEVEL,LOCCOL
    INTEGER ILOOP,ISTAR,IEND,ISTEP
    LOGICAL SIMPLE
    INTEGER, ALLOCATABLE, DIMENSION(:)::CWORK
    INTEGER, ALLOCATABLE, DIMENSION(:)::HWORK
    INTEGER, ALLOCATABLE, DIMENSION(:)::LOC2GLOB
    INTEGER, ALLOCATABLE, DIMENSION(:)::GLOB2LOC
    REAL, ALLOCATABLE, DIMENSION(:,:)::MAT
    REAL, ALLOCATABLE, DIMENSION(:)::B
    REAL, ALLOCATABLE, DIMENSION(:)::XSOL
    REAL, ALLOCATABLE, DIMENSION(:)::W
    REAL, ALLOCATABLE, DIMENSION(:)::RLEV
    REAL, ALLOCATABLE, DIMENSION(:)::RADP
    REAL, ALLOCATABLE, DIMENSION(:)::RADN
    ewrite(1,*) '-IN DIFSOLVCOLNS'
    ewrite(2,*) 'totele',totele

    NPLEVEL=2*(NLEVEL-1)+1
    IF((NLOC.EQ.0).OR.(NLOC.EQ.MLOC)) NPLEVEL=NLEVEL
    CINLEV=PGNODS/NPLEVEL
    INLEV=NONODS/NLEVEL

    ALLOCATE(CWORK(PGNODS))
    ALLOCATE(HWORK(PGNODS))
    ALLOCATE(LOC2GLOB(PGNODS))
    ALLOCATE(GLOB2LOC(PGNODS))

    ALLOCATE(MAT(NPLEVEL,NPLEVEL))
    ALLOCATE(B(NPLEVEL))
    ALLOCATE(XSOL(NPLEVEL))

    DP(1:PGNODS) = 0.0

    ! CWORK is 1 at the bottom and 2*(NLEVEL-1)+1
    ewrite(2,*) 'totele,nlevel,nplevel,cinlev,inlev:',totele,nlevel,nplevel,cinlev,inlev
    IF((NLOC.EQ.0).OR.(NLOC.EQ.MLOC)) THEN
       DO NOD=1,NONODS
          CWORK(NOD)=INT((NOD-1)/INLEV)+1
          HWORK(NOD)=NOD-(CWORK(NOD)-1)*INLEV
       END DO
    ELSE
       DO ELE=1,TOTELE
          DO ILOC=1,NLOC
             INOD=NDGLNO((ELE-1)*NLOC+ILOC)
             DO JLOC=ILOC,NLOC
                JNOD=NDGLNO((ELE-1)*NLOC+JLOC)
                KLOC=ILINK2(ILOC,JLOC)
                PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                CWORK(PNOD)=INT((INOD-1)/INLEV)+1+INT((JNOD-1)/INLEV)+1 -1
             END DO
          END DO
       END DO
       DO PNOD=1,PGNODS
          HWORK(PNOD)=P2SMAP(ABS(IWORK(PNOD)))
       END DO
    ENDIF

    DO PNOD=1,PGNODS
       LOCNOD=CWORK(PNOD)  +  (HWORK(PNOD)-1)*NPLEVEL
       LOC2GLOB(LOCNOD)=PNOD
       GLOB2LOC(PNOD)=LOCNOD
    END DO

    SIMPLE=.true.
    IF(SIMPLE) THEN

       DO ILOOP=1,2
          ISTAR=1
          IEND=CINLEV
          ISTEP=1
          IF(ILOOP.EQ.2) THEN
             ISTAR=CINLEV
             IEND=1
             ISTEP=-1
          ENDIF
          ewrite(2,*) 'iloop=',iloop
          isolv_loop: DO ISOLV=ISTAR,IEND,ISTEP
             MAT = 0.0
             ! Form and solv matrix eqn...
             DO IPLEVEL=1,NPLEVEL
                LOCNOD=IPLEVEL+(ISOLV-1)*NPLEVEL
                PNOD=LOC2GLOB(LOCNOD)
                B(IPLEVEL)=DIFVECN(PNOD)
                DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
                   PCOL=PGCOLM(COUNT)
                   ! Put this contribution into matrix MAT
                   IF(HWORK(PNOD).EQ.HWORK(PCOL)) THEN
                      LOCCOL=GLOB2LOC(PCOL)
                      JPLEVEL=LOCCOL-(ISOLV-1)*NPLEVEL
                      MAT(IPLEVEL,JPLEVEL)=MAT(IPLEVEL,JPLEVEL)+DIFMATH(COUNT)
                   ENDIF
                END DO
             END DO
             ! do FBGS sweep...
             DO  IPLEVEL=1,NPLEVEL
                LOCNOD=IPLEVEL+(ISOLV-1)*NPLEVEL
                PNOD=LOC2GLOB(LOCNOD)
                DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
                   PCOL=PGCOLM(COUNT)
                   IF(HWORK(PNOD).NE.HWORK(PCOL)) then
                      B(IPLEVEL)=B(IPLEVEL)-DIFMATH(COUNT)*DP(PCOL)
                   endif
                END DO
             END DO

             CALL SMLINN(MAT,XSOL,B,NPLEVEL,NPLEVEL)
             DO IPLEVEL=1,NPLEVEL
                LOCNOD=IPLEVEL+(ISOLV-1)*NPLEVEL
                PNOD=LOC2GLOB(LOCNOD)
                DP(PNOD)=XSOL(IPLEVEL)
             END DO
          END DO isolv_loop
       END DO

    ELSE

       ! Preconditioning by assuming the ocean is 1D...
       ALLOCATE(W(NPLEVEL))
       ALLOCATE(RLEV(NPLEVEL))
       ALLOCATE(RADN(NONODS))
       ALLOCATE(RADP(PGNODS))
       ! Calculate RAD for interpolating
       DO NOD=1,NONODS
          IF(ISPHERE.EQ.2) THEN
             RADN(NOD)=SQRT(X(NOD)**2+Y(NOD)**2+Z(NOD)**2)
          ELSE
             RADN(PNOD)=Z(NOD)
          ENDIF
       END DO

       IF((NLOC.EQ.0).OR.(NLOC.EQ.MLOC)) THEN
          RADP(1:NONODS) = RADN(1:NONODS)
       ELSE
          DO ELE=1,TOTELE
             DO ILOC=1,NLOC
                INOD=NDGLNO((ELE-1)*NLOC+ILOC)
                DO JLOC=ILOC,NLOC
                   JNOD=NDGLNO((ELE-1)*NLOC+JLOC)
                   KLOC=ILINK2(ILOC,JLOC)
                   PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                   RADP(PNOD)=0.5*(RADN(INOD)+RADN(JNOD))
                END DO
             END DO
          END DO
       ENDIF

       DO ISOLV=1,CINLEV
          MAT = 0.0
          DO IPLEVEL=1,NPLEVEL
             LOCNOD=IPLEVEL+(ISOLV-1)*NPLEVEL
             PNOD=LOC2GLOB(LOCNOD)
             RLEV(IPLEVEL)=RADP(PNOD)
          END DO
          ! Form and solv matrix eqn...
          DO IPLEVEL=1,NPLEVEL
             LOCNOD=IPLEVEL+(ISOLV-1)*NPLEVEL
             PNOD=LOC2GLOB(LOCNOD)
             B(IPLEVEL)=DIFVECN(PNOD)

             DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
                PCOL=PGCOLM(COUNT)
                ! Calculate the weight W by interpolating...
                RINT=RADP(PCOL)
                ! Now interpolate...
                W(1:NPLEVEL) = 0.0
                IF(RINT.Lt.RLEV(1)) THEN
                   W(1) =0.0
                   if(hwork(pnod).eq.hwork(pcol)) w(1)=1.0
                ELSE IF(RINT.GE.RLEV(NPLEVEL)) THEN
                   W(NPLEVEL)=1.0
                ELSE
                   IF((NLOC.EQ.0).OR.(NLOC.EQ.MLOC)) THEN
                      ! Linear interpolation...
                      DO IW=1,NPLEVEL-1
                         IF((RINT.GE.RLEV(IW)).AND.(RINT.LE.RLEV(IW+1))) THEN
                            W(IW)  =(RLEV(IW+1)-RINT)/(RLEV(IW+1)-RLEV(IW))
                            W(IW+1)=1.0-W(IW)
                         ENDIF
                      END DO
                   ELSE
                      ! Quadratic interpolation...
                      DO IW=1,NPLEVEL-2,2
                         IF((RINT.GE.RLEV(IW)).AND.(RINT.LE.RLEV(IW+2))) THEN
                            LX=2.0*(RINT-RLEV(IW))/(RLEV(IW+2)-RLEV(IW)) - 1.0
                            W(IW)  =0.5*LX*(LX-1.0)
                            W(IW+1)=(1.0-LX**2)
                            W(IW+2)=0.5*LX*(LX+1.0)
                         ENDIF
                      END DO
                   ENDIF
                ENDIF

                ! Look for contribution in the same level...
                coef=1.0
                if(hwork(pnod).eq.hwork(pcol)) coef=1.0
                DO IW=1,NPLEVEL
                   MAT(IPLEVEL,IW)=MAT(IPLEVEL,IW)+coef*W(IW)*DIFMATH(COUNT)
                END DO
             END DO
          END DO
          ! Put top bc in...
          DO JPLEVEL=1,NPLEVEL
             MAT(NPLEVEL,JPLEVEL)=0.0
             MAT(JPLEVEL,NPLEVEL)=0.0
          END DO
          MAT(NPLEVEL,NPLEVEL)=1.0
          B(NPLEVEL)=0.0

          CALL SMLINN(MAT,XSOL,B,NPLEVEL,NPLEVEL)
          DO IPLEVEL=1,NPLEVEL
             LOCNOD=IPLEVEL+(ISOLV-1)*NPLEVEL
             PNOD=LOC2GLOB(LOCNOD)
             DP(PNOD)=XSOL(IPLEVEL)
          END DO

       END DO
       ! ENDOF IF(SIMPLE) THEN ELSE...
    ENDIF
    ewrite(1,*) 'just leaving DIFSOLVCOLNS'

  END SUBROUTINE DIFSOLVCOLNS

  SUBROUTINE DIFMATMULT(DIFMATH,PGV,DIFVECV, &
       PGFINDRM,PGCOLM,NPGCOLM,PGNODS)
    ! Subtract DIFVECV=-DIFMATH PGV  +  DIFVEC
    INTEGER PGNODS,NPGCOLM
    REAL DIFMATH(NPGCOLM),PGV(PGNODS),DIFVECV(PGNODS)
    INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    INTEGER PNOD,COUNT
    DO PNOD=1,PGNODS
       DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
          DIFVECV(PNOD)=DIFVECV(PNOD)-DIFMATH(COUNT)*PGV(PGCOLM(COUNT))
       END DO
    END DO

  END SUBROUTINE DIFMATMULT

  SUBROUTINE FINDDIFVECN(NONODS,XNONOD,PGNODS,TOTELE,NLOC,MLOC, &
       NDGLNO,PGNDGLNO, &
       PNORMX,PNORMY,PNORMZ, &
       COGRAX,&
       BSOUX,BSOUY,BSOUZ, &
       VGRAVX,VGRAVY,VGRAVZ, &
       ISPHERE,X,Y,Z )

    ! rotat components (DIFVECX,DIFVECY,DIFVECZ)  U'=RU **********
    !  The rotation matrix in 3-D is R=
    !   T1X   T1Y   T1Z
    !   T2X   T2Y   T2Z
    !   NORMX NORMY NORMZ
    INTEGER NONODS,XNONOD,PGNODS,TOTELE,NLOC,MLOC
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    REAL PNORMX(PGNODS),PNORMY(PGNODS),PNORMZ(PGNODS)
    logical COGRAX
    REAL BSOUX,BSOUY,BSOUZ
    REAL VGRAVX(NONODS),VGRAVY(NONODS),VGRAVZ(NONODS)
    INTEGER ISPHERE
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    ! Local variables...
    INTEGER I,ELE,ILOC,JLOC,KLOC,PNOD,INOD,JNOD
    REAL, dimension(1):: T1X,T1Y,T1Z,T2X,T2Y,T2Z
    REAL RN
    REAL, ALLOCATABLE, DIMENSION(:)::NORMX
    REAL, ALLOCATABLE, DIMENSION(:)::NORMY
    REAL, ALLOCATABLE, DIMENSION(:)::NORMZ

    IF(ISPHERE.GE.2) THEN
       DO ELE=1,TOTELE
          DO ILOC=1,NLOC
             DO JLOC=1,NLOC
                KLOC=ILINK2(ILOC,JLOC)
                INOD =NDGLNO((ELE-1)*NLOC+ILOC)
                JNOD =NDGLNO((ELE-1)*NLOC+JLOC)
                PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                PNORMX(PNOD)=0.5*(X(INOD)+X(JNOD))
                PNORMY(PNOD)=0.5*(Y(INOD)+Y(JNOD))
                PNORMZ(PNOD)=0.5*(Z(INOD)+Z(JNOD))
             END DO
          END DO
       END DO
       DO PNOD=1,PGNODS
          RN=SQRT(PNORMX(PNOD)**2+PNORMY(PNOD)**2+PNORMZ(PNOD)**2)
          PNORMX(PNOD)=PNORMX(PNOD)/RN
          PNORMY(PNOD)=PNORMY(PNOD)/RN
          PNORMZ(PNOD)=PNORMZ(PNOD)/RN
       END DO

    ELSE

       ALLOCATE(NORMX(NONODS))
       ALLOCATE(NORMY(NONODS))
       ALLOCATE(NORMZ(NONODS))

       ! Calculate normal:
       DO I=1,NONODS
          !  The rotation matrix in 3-D is R=
          !   T1X   T1Y   T1Z
          !   T2X   T2Y   T2Z
          !   NORMX NORMY NORMZ
          CALL GETROTGRAV(&
               T1X, T1Y, T1Z, &
               T2X, T2Y, T2Z, &
               NORMX(i:i), NORMY(I:i), NORMZ(I:i), &
               COGRAX,BSOUX,BSOUY,BSOUZ, &
               (/ VGRAVX(I) /), (/ VGRAVY(I) /), (/ VGRAVZ(I) /), 1)
       END DO

       ! interpolate normal to pressure nodes...
       DO ELE=1,TOTELE
          DO ILOC=1,NLOC
             DO JLOC=1,NLOC
                KLOC=ILINK2(ILOC,JLOC)
                INOD =NDGLNO((ELE-1)*NLOC+ILOC)
                JNOD =NDGLNO((ELE-1)*NLOC+JLOC)
                PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                PNORMX(PNOD)=0.5*(NORMX(INOD)+NORMX(JNOD))
                PNORMY(PNOD)=0.5*(NORMY(INOD)+NORMY(JNOD))
                PNORMZ(PNOD)=0.5*(NORMZ(INOD)+NORMZ(JNOD))
             END DO
          END DO
       END DO
       DO PNOD=1,PGNODS
          RN=SQRT(PNORMX(PNOD)**2+PNORMY(PNOD)**2+PNORMZ(PNOD)**2)
          PNORMX(PNOD)=PNORMX(PNOD)/RN
          PNORMY(PNOD)=PNORMY(PNOD)/RN
          PNORMZ(PNOD)=PNORMZ(PNOD)/RN
       END DO

       DEALLOCATE(NORMX)
       DEALLOCATE(NORMY)
       DEALLOCATE(NORMZ)
    ENDIF

  END SUBROUTINE FINDDIFVECN

  SUBROUTINE SOLVERLUM(TOTELE,NLOC,MLOC,NONODS,PGNODS,NCOLM, &
       NDGLNO,PGNDGLNO, CENTRM,FINDRM,COLM,X,Y,Z,XNONOD,ISPHERE, &
       NPGCOLM,PGFINDRM,PGCOLM, &
       NLEVEL, &
       DIFMATH,DIFVECV, &
       PGH, &
       halotagT10,CORR3D,LUMCOL,USEGMRES, &
       option_path)
    ! lump vertical eqns and solve And lump into colns and solve:
    ! If LUMCOL then lump into vertical colns as well and solve.
    INTEGER TOTELE,NLOC,MLOC,XNONOD,NONODS,PGNODS,NCOLM
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    INTEGER CENTRM(NONODS),FINDRM(NONODS+1),COLM(NCOLM)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER ISPHERE
    INTEGER NPGCOLM
    INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    INTEGER NLEVEL
    REAL DIFMATH(NPGCOLM),DIFVECV(PGNODS)
    REAL PGH(PGNODS)
    INTEGER halotagT10
    LOGICAL CORR3D,LUMCOL,USEGMRES
    ! path in options tree to field we're solving for to query solver options
    character(len=*), intent(in):: option_path
    ! Local variables...
    INTEGER I,KITS,NSGCOLM,SGNODS,PNOD,SNOD,COUNT
    
    INTEGER MINTE1,MINTE2,TLEFT1,TLEFT2,NV,NR2
    
    INTEGER,  ALLOCATABLE, DIMENSION(:)::IWORK
    INTEGER,  ALLOCATABLE, DIMENSION(:)::P2SMAP
    REAL, ALLOCATABLE, DIMENSION(:)::WORK1
    REAL, ALLOCATABLE, DIMENSION(:)::WORK2
    REAL, ALLOCATABLE, DIMENSION(:)::WORK3
    REAL, ALLOCATABLE, DIMENSION(:)::WORK4
    REAL, ALLOCATABLE, DIMENSION(:)::WORK5
    REAL, ALLOCATABLE, DIMENSION(:)::WORKP
    REAL, ALLOCATABLE, DIMENSION(:)::WORKG
    INTEGER,  ALLOCATABLE, DIMENSION(:)::SGFINDRM
    INTEGER,  ALLOCATABLE, DIMENSION(:)::SGCOLM
    INTEGER,  ALLOCATABLE, DIMENSION(:)::SGCENTRM
    REAL, ALLOCATABLE, DIMENSION(:)::SGH
    REAL, ALLOCATABLE, DIMENSION(:)::SDIFVECV
    REAL, ALLOCATABLE, DIMENSION(:)::SDIFMATH

    REAL, ALLOCATABLE, DIMENSION(:)::DIFVECN
    REAL, ALLOCATABLE, DIMENSION(:)::DP

    ALLOCATE(IWORK(PGNODS))

    ! Calculate IWORK(PNOD)=-ve the node at the top of the domain and
    ! +ve if a node at the top of the domain.
    ewrite(1,*) 'going into NDTOPOINT'
    CALL NDTOPOINT(IWORK, TOTELE,NLOC,NONODS,PGNODS,NCOLM, &
         NDGLNO,CENTRM,FINDRM,COLM, &
         NLEVEL)
    ewrite(1,*) 'calculating SGNODS,NSGCOLM'
    NSGCOLM=0
    SGNODS=0
    DO PNOD=1,PGNODS
       IF(IWORK(PNOD).EQ.PNOD) THEN
          SGNODS=SGNODS+1
          DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
             IF(IWORK(PGCOLM(COUNT)).EQ.PGCOLM(COUNT)) NSGCOLM=NSGCOLM+1
          END DO
       END IF
    END DO

    ALLOCATE(P2SMAP(PGNODS))
    ewrite(1,*) 'CALCULATING P2SMAP'

    SNOD=0
    DO PNOD=1,PGNODS
       IF(IWORK(PNOD).EQ.PNOD) THEN
          SNOD=SNOD+1
          P2SMAP(PNOD)=SNOD
       ELSE
          P2SMAP(PNOD)=0
       ENDIF
    END DO

    DO PNOD=1,PGNODS
       IF(P2SMAP(PNOD).EQ.0) P2SMAP(PNOD)=P2SMAP(ABS(IWORK(PNOD)))
    END DO

    ALLOCATE(SGFINDRM(SGNODS+1))
    ALLOCATE(SGCOLM(NSGCOLM))
    ALLOCATE(SGCENTRM(SGNODS))
    ALLOCATE(SGH(SGNODS))
    ALLOCATE(SDIFVECV(SGNODS))
    ALLOCATE(SDIFMATH(NSGCOLM))

    ! Now lump DIFMATH & DIFVECV so that is only a horizontal system...
    ewrite(1,*) 'going into DIFMATLUMP'
    CALL DIFMATLUMP(DIFMATH,DIFVECV,IWORK,PGNODS, &
         PGFINDRM,PGCOLM,NPGCOLM, &
         SDIFMATH,SDIFVECV,SGFINDRM,SGCOLM,SGCENTRM,NSGCOLM,SGNODS,P2SMAP)

    ewrite(2,*) 'pgnods,sgnods,pgnods/sgnods:',pgnods,sgnods,pgnods/sgnods

    ALLOCATE(WORK1(SGNODS))
    ALLOCATE(WORK2(SGNODS))
    ALLOCATE(WORK3(SGNODS))
    ALLOCATE(WORK4(SGNODS))
    ALLOCATE(WORK5(SGNODS))
    ALLOCATE(WORKP(4*SGNODS))
    ! INITIAL GUESS
    DO I=1,PGNODS
       SGH(P2SMAP(ABS(IWORK(I))))=PGH(I)
    END DO
    ewrite(1,*) 'solving 2D lumped system'

    IF(USEGMRES) THEN
       MINTE1=40
       MINTE2=40
       TLEFT1=1
       TLEFT2=1
       NV=(MINTE1+1)*SGNODS
       NR2=NV

       ALLOCATE(WORKG(NV))
       DEALLOCATE(WORKP)
       ALLOCATE(WORKP(NR2))
       CALL GMRES(SGH,SDIFVECV,SGNODS,SGNODS,SGNODS,.TRUE., &
            SDIFMATH,SGFINDRM,SGCOLM,NSGCOLM, NSGCOLM,&
            halotagT10,KITS, &
            option_path=option_path)
    ELSE
       CALL SOLCG(SGH,SDIFVECV,SGNODS,SGNODS,SGNODS,.TRUE.,  &
            SDIFMATH,SGFINDRM,SGCOLM,NSGCOLM, NSGCOLM,&
            halotagT10,KITS, &
            option_path=option_path)
    ENDIF


    ! Now distribute PGH in the vertical
    DO I=1,PGNODS
       PGH(I)=SGH(P2SMAP(ABS(IWORK(I))))
    END DO

    IF(LUMCOL) THEN
       IF(CORR3D) THEN
          ! *************SOLVE VERTICAL COLNS*************
          ! Now lump DIFMATH so that is only a series of vertical systems

          ALLOCATE(DIFVECN(PGNODS))
          ALLOCATE(DP(PGNODS))
          ! Subtract DIFVECN=-DIFMATH PGH  +  DIFVECV
          DO I=1,PGNODS
             DIFVECN(I)=DIFVECV(I)
          END DO
          CALL DIFMATMULT(DIFMATH,PGH,DIFVECN, &
               PGFINDRM,PGCOLM,NPGCOLM,PGNODS)

          CALL DIFSOLVCOLNS(DIFMATH,DP,DIFVECN,IWORK,P2SMAP,PGNODS,NONODS, &
               PGFINDRM,PGCOLM,NPGCOLM,NLEVEL, &
               TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO, &
               X,Y,Z,XNONOD,ISPHERE)

          DO I=1,PGNODS
             PGH(I)=DP(I)+PGH(I)
          END DO
          ! *************SOLVE VERTICAL COLNS*************
       ENDIF
    ENDIF
    ewrite(1,*) 'returning from solverlum'

  END SUBROUTINE SOLVERLUM

  SUBROUTINE TOPBCDIFDIA(DIFMATH, &
       PGFINDRM,PGCOLM,PGCENTRM,NPGCOLM,PGNODS,NONODS, &
       TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO,VECRHS,NLEVEL,CORNER)
    ! Amend matrix to apply b.c at the free surface
    ! IF CORNER then apply b.c only to corner nodes of element.
    LOGICAL CORNER
    INTEGER NPGCOLM,PGNODS,NLEVEL,nonods
    INTEGER PGCENTRM(PGNODS),PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    INTEGER TOTELE,NLOC,MLOC,NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    REAL DIFMATH(NPGCOLM),VECRHS(PGNODS)
    ! Local variables...
    INTEGER nod
    INTEGER PNOD,COL,ELE,ILOC,JLOC,KLOC,INOD,JNOD,COUNT
    REAL SUM
    REAL, ALLOCATABLE, DIMENSION(:)::WORKS
    REAL, ALLOCATABLE, DIMENSION(:)::WORKL
    REAL, ALLOCATABLE, DIMENSION(:)::DIAG

    ALLOCATE(WORKS(NONODS))
    ALLOCATE(WORKL(PGNODS))
    ALLOCATE(DIAG(PGNODS))

    WORKS(1:NONODS) = 0.0
    DO NOD=(NLEVEL-1)*NONODS/NLEVEL +1, NONODS
       WORKS(NOD)=1.0
    END DO

    WORKL(1:PGNODS) = 0.0
    DO ELE=1,TOTELE
       DO ILOC=1,NLOC
          IF(CORNER) THEN
             ! Apply b.c only to corner nodes of element..
             KLOC=ILINK2(ILOC,ILOC)
             INOD =NDGLNO((ELE-1)*NLOC+ILOC)
             PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
             WORKL(PNOD)=WORKS(INOD)
          ELSE
             ! Apply b.c to all nodes of element...
             DO JLOC=1,NLOC
                KLOC=ILINK2(ILOC,JLOC)
                INOD =NDGLNO((ELE-1)*NLOC+ILOC)
                JNOD =NDGLNO((ELE-1)*NLOC+JLOC)
                PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                WORKL(PNOD)=0.5*(WORKS(INOD)+WORKS(JNOD))
             END DO
          ENDIF
       END DO
    END DO
    COUNT=0

    ! Lump the rows:
    DIAG(1:PGNODS) = 0.0
    DO PNOD=1,PGNODS
       IF(WORKL(PNOD).GT.0.75) THEN
          SUM=0.0
          DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
             SUM=SUM+ABS(DIFMATH(COUNT))
             DIFMATH(COUNT)=0.0
          END DO
          DIAG(PNOD)=SUM
          VECRHS(PNOD)=0.0
       END IF
    END DO

    ! set the corresponding colums to zero
    DO PNOD=1,PGNODS
       DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
          COL=PGCOLM(COUNT)
          IF(WORKL(COL).GT.0.75) DIFMATH(COUNT)=0.0
       END DO
    END DO

    DO PNOD=1,PGNODS
       DIFMATH(PGCENTRM(PNOD))=DIFMATH(PGCENTRM(PNOD))+DIAG(PNOD)
    END DO

  END SUBROUTINE TOPBCDIFDIA

  SUBROUTINE LONGFREESFIMPBCS(NOBCFSimp,BCT1WFSimp,BCT2FSimp,DT, &
       NONODS, &
       PGFINDRM,PGCENTRM,PGCOLM,NPGCOLM,PGNODS, &
       LONGDIFMATH,LONGDIFVECV, &
       FREES,FREESOLD, &
       TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO)
    ! Amend matrix to apply b.c at the free surface
    ! IF CORNER then apply b.c only to corner nodes of element.
    INTEGER NOBCFSimp
    REAL BCT1WFSimp(NOBCFSimp)
    INTEGER BCT2FSimp(NOBCFSimp)
    REAL DT
    INTEGER NPGCOLM,PGNODS,NONODS
    INTEGER PGCENTRM(PGNODS),PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    INTEGER TOTELE,NLOC,MLOC
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    REAL FREES(PGNODS)
    REAL FREESOLD(PGNODS)
    REAL LONGDIFMATH(NPGCOLM),LONGDIFVECV(PGNODS)
    ! Local variables...
    INTEGER NOD, II
    INTEGER ELE,ILOC,JLOC,KLOC,INOD,JNOD,COUNT,PNOD,PCOL
    REAL RDIAG,RBCVALL

    INTEGER, ALLOCATABLE, DIMENSION(:)::WORKS
    INTEGER, ALLOCATABLE, DIMENSION(:)::WORKL
    REAL, ALLOCATABLE, DIMENSION(:)::BCVALS
    REAL, ALLOCATABLE, DIMENSION(:)::TEMPBCVALL

    ALLOCATE(WORKS(NONODS))
    ALLOCATE(WORKL(PGNODS))
    ALLOCATE(BCVALS(NONODS))
    ALLOCATE(TEMPBCVALL(PGNODS))
    ! Calculate node values of FREES i.e. NFREES

    WORKS(1:NONODS)=0
    WORKL(1:PGNODS)=0
    BCVALS(1:NONODS)=0.0
    TEMPBCVALL(1:PGNODS)=0.0
    DO II=1,NOBCFSimp
       NOD=BCT2FSimp(II)
       WORKS(NOD)=1
       BCVALS(NOD)=BCT1WFSimp(II)
    END DO

    DO ELE=1,TOTELE
       DO ILOC=1,NLOC
          DO JLOC=ILOC,NLOC
             KLOC=ILINK2(ILOC,JLOC)
             INOD=NDGLNO((ELE-1)*NLOC+ILOC)
             JNOD=NDGLNO((ELE-1)*NLOC+JLOC)
             PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
             IF(WORKS(INOD)+WORKS(JNOD).GE.2) THEN
                WORKL(PNOD)=1
                RBCVALL=0.5*(BCVALS(INOD)+BCVALS(JNOD))
                TEMPBCVALL(PNOD)=DT*RBCVALL+FREESOLD(PNOD) -FREES(PNOD)
             ENDIF
          END DO
       END DO
    END DO

    ! Reset b.c to take into account that we are solving for Delta \xi
    ! and that BCT1WFSimp contains the rate of change of free surface height.

    ! Set the corresponding colums to zero
    DO PNOD=1,PGNODS
       DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
          PCOL=PGCOLM(COUNT)
          IF(WORKL(PCOL).EQ.1) THEN
             LONGDIFVECV(PNOD)=LONGDIFVECV(PNOD)-LONGDIFMATH(COUNT)*TEMPBCVALL(PCOL)
             LONGDIFMATH(COUNT)=0.0
          ENDIF
       END DO
    END DO

    DO PNOD=1,PGNODS
       IF(WORKL(PNOD).EQ.1) THEN
          RDIAG=0.0
          DO COUNT=PGFINDRM(PNOD),PGFINDRM(PNOD+1)-1
             RDIAG=RDIAG+ABS(LONGDIFMATH(COUNT))
             LONGDIFMATH(COUNT)=0.0
          END DO
          LONGDIFMATH(PGCENTRM(PNOD))=RDIAG
          LONGDIFVECV(PNOD)=TEMPBCVALL(PNOD)*RDIAG
       ENDIF
    END DO

  END SUBROUTINE LONGFREESFIMPBCS

  SUBROUTINE FREESFIMPBCS(NOBCFSimp,BCT1WFSimp,BCT2FSimp,DT, &
       NONODS,FINDRM,CENTRM,COLM,NCOLM, &
       SHORTDIFMATH,SHORTDIFVECV, &
       FREES,FREESOLD,PGNODS, &
       TOTELE,NLOC,MLOC,NDGLNO,PGNDGLNO)
    ! Amend matrix to apply b.c at the free surface
    ! IF CORNER then apply b.c only to corner nodes of element.
    INTEGER NOBCFSimp
    REAL BCT1WFSimp(NOBCFSimp)
    INTEGER BCT2FSimp(NOBCFSimp)
    REAL DT
    INTEGER NCOLM,PGNODS,NONODS
    INTEGER CENTRM(NONODS),FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER TOTELE,NLOC,MLOC
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    REAL FREES(PGNODS)
    REAL FREESOLD(PGNODS)
    REAL SHORTDIFMATH(NCOLM),SHORTDIFVECV(NONODS)
    ! Local variables...
    INTEGER NOD, II
    INTEGER COL,ELE,ILOC,KLOC,COUNT,PNOD
    REAL RDIAG

    INTEGER, ALLOCATABLE, DIMENSION(:)::WORKL
    REAL, ALLOCATABLE, DIMENSION(:)::BCVAL
    REAL, ALLOCATABLE, DIMENSION(:)::NFREES
    REAL, ALLOCATABLE, DIMENSION(:)::NFREESOLD
    REAL, ALLOCATABLE, DIMENSION(:)::TEMPBCT1WFSimp

    ALLOCATE(WORKL(NONODS))
    ALLOCATE(BCVAL(NONODS))
    ALLOCATE(NFREES(NONODS))
    ALLOCATE(NFREESOLD(NONODS))
    ALLOCATE(TEMPBCT1WFSimp(NOBCFSimp))
    ! Calculate node values of FREES i.e. NFREES
    DO ELE=1,TOTELE
       DO ILOC=1,NLOC
          KLOC=ILINK2(ILOC,ILOC)
          NOD=NDGLNO((ELE-1)*NLOC+ILOC)
          PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
          NFREES(NOD) =FREES(PNOD)
          NFREESOLD(NOD) =FREESOLD(PNOD)
       END DO
    END DO

    ! Reset b.c to take into account that we are solving for Delta \xi
    ! and that BCT1WFSimp contains the rate of change of free surface height.
    BCVAL(1:NONODS) = 0.0
    WORKL(1:NONODS) = 0
    DO II=1,NOBCFSimp
       NOD=BCT2FSimp(II)
       IF((NOD.LT.1).OR.(NOD.GT.NONODS)) THEN
          ewrite(-1,*) 'BCT2FSimp does not contain the correct nodes'
          ewrite(-1,*) 'nonods,nod:',NONODS,NOD
          FLAbort("Bye now!")
       ENDIF
       TEMPBCT1WFSimp(II)=DT*BCT1WFSimp(II)+NFREESOLD(NOD) -NFREES(NOD)
       WORKL(NOD)=1
       BCVAL(NOD)=TEMPBCT1WFSimp(II)
    END DO

    ! Set the corresponding colums to zero
    DO NOD=1,NONODS
       DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
          COL=COLM(COUNT)
          IF(WORKL(COL).EQ.1) THEN
             SHORTDIFVECV(NOD)=SHORTDIFVECV(NOD)-SHORTDIFMATH(COUNT)*BCVAL(COL)
             SHORTDIFMATH(COUNT)=0.0
          ENDIF
       END DO
    END DO

    DO II=1,NOBCFSimp
       NOD=BCT2FSimp(II)
       RDIAG=0.0
       DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
          RDIAG=RDIAG+ABS(SHORTDIFMATH(COUNT))
          SHORTDIFMATH(COUNT)=0.0
       END DO
       SHORTDIFMATH(CENTRM(NOD))=RDIAG
       SHORTDIFVECV(NOD)=TEMPBCT1WFSimp(II)*RDIAG
    END DO

  END SUBROUTINE FREESFIMPBCS

  SUBROUTINE MULGRIDVERT(DIFMATH,DIFVECV, &
       FINDRM,COLM,CENTRM,NCOLM,NONODS, &
       PGFINDRM,PGCOLM,PGCENTRM,NPGCOLM,PGNODS, &
       PG, &
       TOTELE,NLOC,MLOC, &
       NDGLNO,PGNDGLNO, X,Y,Z,XNONOD,ISPHERE, &
       NLEVEL, &
       PGNODSP, &
       halotagT10,PPSET, &
       CORR3D, ONLY3D,USEGMRES,LUMCOL, option_path)
    ! This performs a simplified multi-grid iteration by lumping horizontal eqns.
    ! Solve DIFMATH*PG=DIFVECV
    ! If ONLY3D then solve a full 3D system only.
    ! if LUMCOL then solve a system of coln lumped eqns to help the soln as well.
    LOGICAL GETDIAG,RESETDIA
    PARAMETER(GETDIAG=.TRUE.,RESETDIA=.FALSE.)
    INTEGER NPGCOLM,PGNODS,NCOLM,NONODS
    INTEGER CENTRM(NONODS),FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER PGCENTRM(PGNODS),PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
    REAL DIFMATH(NPGCOLM),DIFVECV(PGNODS)
    REAL PG(PGNODS)
    ! Extra variables from SOLVERLUM
    INTEGER TOTELE,NLOC,MLOC,XNONOD
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER ISPHERE
    INTEGER NLEVEL
    INTEGER PGNODSP
    INTEGER halotagT10,PPSET
    LOGICAL, intent(in):: CORR3D, ONLY3D
    logical:: USEGMRES,LUMCOL
    character(len=*), intent(in):: option_path
    ! Local variables...
    INTEGER I,KITS
    LOGICAL VERLUM, checkconvergence
    logical surface_equation
    INTEGER MINTE1,MINTE2,TLEFT1,TLEFT2,NV,NR2

    REAL, ALLOCATABLE, DIMENSION(:)::DIFVECN
    REAL, ALLOCATABLE, DIMENSION(:)::DP
    REAL, ALLOCATABLE, DIMENSION(:)::WORK1
    REAL, ALLOCATABLE, DIMENSION(:)::WORK2
    REAL, ALLOCATABLE, DIMENSION(:)::WORK3
    REAL, ALLOCATABLE, DIMENSION(:)::WORK4
    REAL, ALLOCATABLE, DIMENSION(:)::WORK5
    REAL, ALLOCATABLE, DIMENSION(:)::WORKP
    REAL, ALLOCATABLE, DIMENSION(:)::WORKG
    ewrite(1,*) 'just inside MULGRIDVERT'

    ALLOCATE(DIFVECN(PGNODS))
    ALLOCATE(DP(PGNODS))
    ewrite(1,*) 'just allocated memory inside MULGRIDVERT'

    ! If we're solving a free surface equation
    ! and we can't do vertical lumping (nlevel==0)
    ! we solve with an ordinary 3D solve.
    ! Afterwards the solution has to be extrapolated from the surface
    ! to the interior of the domain.
    surface_equation= (nlevel==0) .and. .not. corr3d
    ewrite(1,*)'nlevel,corr3d,surface_equation:',nlevel,corr3d,surface_equation

    if (surface_equation) then
       ! if we're solving a free surface equation
       ! and we can't do vertical lumping
       ! just do a 3d solve, but make sure we don't have zero diagonals first
       do i=1, pgnods
          if (difmath(pgcentrm(i))==0.0) then
             difmath(pgcentrm(i))=1.0
             difvecv(i)=0.0
          end if
       end do

    end if

    ! Iterate ***************************
    ! lump vertical eqns and solve:
    VERLUM=.true.
    IF(.NOT. ONLY3D .and. .not. surface_equation) THEN
    
      IF(VERLUM) THEN
         ! Subtract DIFVECN=-DIFMATH PG  +  DIFVECV *******
         DO I=1,PGNODS
            DIFVECN(I)=DIFVECV(I)
         END DO
         ! Subtract DIFVECN=-DIFMATH PG  +  DIFVECN
         CALL DIFMATMULT(DIFMATH,PG,DIFVECN, &
              PGFINDRM,PGCOLM,NPGCOLM,PGNODS)
         IF(PPSET.GT.0) DIFVECN(PPSET)=0.0

         ! ******************SOLVING VER LUMPED then HORIZ LUMPED**
         ! Solve DIFMATH*DP=DIFVECN for DP
         CALL SOLVERLUM(TOTELE,NLOC,MLOC,NONODS,PGNODS,NCOLM, &
              NDGLNO,PGNDGLNO, CENTRM,FINDRM,COLM,X,Y,Z,XNONOD,ISPHERE, &
              NPGCOLM,PGFINDRM,PGCOLM, &
              NLEVEL, &
              DIFMATH,DIFVECN, &
              DP, &
              halotagT10,CORR3D,LUMCOL,USEGMRES, option_path)
         DO I=1,PGNODS
            PG(I)=DP(I)+PG(I)
         END DO
         ! ******************SOLVING VER LUMPED then HORIZ LUMPED**
         ! ********************************************************
         ! ENDOF IF(VERLUM) THEN
      ENDIF
      ! end of if (.not. only3d)
      ENDIF

      ewrite(2,*) 'CORR3D:',CORR3D
      IF(CORR3D .or. surface_equation) THEN
      ! Fine 3D correction not nec. for very large aspect ratio

      ! Subtract DIFVECN=-DIFMATH PG  +  DIFVECV
      DO I=1,PGNODS
         DIFVECN(I)=DIFVECV(I)
      END DO
      CALL DIFMATMULT(DIFMATH,PG,DIFVECN, &
           PGFINDRM,PGCOLM,NPGCOLM,PGNODS)
      IF(PPSET.GT.0) DIFVECN(PPSET)=0.0

      IF(ONLY3D .or. surface_equation) then
         checkconvergence=.true.
      else
         checkconvergence=.false.
      end if
      KITS=0

      ALLOCATE(WORK1(PGNODS))
      ALLOCATE(WORK2(PGNODS))
      ALLOCATE(WORK3(PGNODS))
      ALLOCATE(WORK4(PGNODS))
      ALLOCATE(WORK5(PGNODS))
      ALLOCATE(WORKP(4*PGNODS))

      IF(USEGMRES) THEN
         MINTE1=40
         MINTE2=40
         TLEFT1=1
         TLEFT2=1
         NV=(MINTE1+1)*PGNODS
         NR2=NV

         ALLOCATE(WORKG(NV))
         DEALLOCATE(WORKP)
         ALLOCATE(WORKP(NR2))
         CALL GMRES(DP, DIFVECN,PGNODS,PGNODS,PGNODSP,.TRUE.,  &
              DIFMATH,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM,&
              halotagT10,KITS, &
              option_path=option_path)
      ELSE
         CALL SOLCG(DP, DIFVECN,PGNODS,PGNODS,PGNODSP,.TRUE.,  &
              DIFMATH,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM,&
              halotagT10,KITS, &
              checkconvergence=checkconvergence, &
              option_path=option_path)
      ENDIF

      DEALLOCATE(WORK1)
      DEALLOCATE(WORK2)
      DEALLOCATE(WORK3)
      DEALLOCATE(WORK4)
      DEALLOCATE(WORK5)
      DEALLOCATE(WORKP)

      DO I=1,PGNODS
         PG(I)=DP(I)+PG(I)
      END DO
      
    END IF

    DEALLOCATE(DIFVECN)
    DEALLOCATE(DP)

  END SUBROUTINE MULGRIDVERT

  subroutine vertical_extrapolation_wrapper(&
       from_val, from_mesh, &
       to_val, to_mesh, &
       coordinates, flat_earth, topdis_field)
    ! subroutine that extrapolates a field from the free surface
    ! into the interior
    real, dimension(:), target, intent(in):: from_val
    type(mesh_type), intent(in):: from_mesh
    real, dimension(:), target, intent(out):: to_val
    type(mesh_type), intent(in):: to_mesh
    type(vector_field), target, intent(inout):: coordinates
    logical, intent(in):: flat_earth
    type(scalar_field), intent(in):: topdis_field

    type(scalar_field) from_field, to_field
    integer, dimension(:), pointer :: surface_element_list

    from_field=wrap_scalar_field(from_mesh, from_val, &
         "SomeFieldDefinedAtTheSurface")
    to_field=wrap_scalar_field(to_mesh, to_val, &
         "VerticallyExtrapolatedField")
    ! list of faces that make up the free surface
    call get_boundary_condition(topdis_field, 1, &
         surface_element_list=surface_element_list)
    call VerticalExtrapolation(from_field, to_field, coordinates, &
         flat_earth, surface_element_list)
    call deallocate(from_field)
    call deallocate(to_field)

  end subroutine vertical_extrapolation_wrapper

  SUBROUTINE CALNDENST(NDENST,SOURCX,SOURCY,SOURCZ, &
       PGRAVX,PGRAVY,PGRAVZ, &
       TOTELE,NLOC,MLOC,NONODS,PGNODS, &
       NDGLNO,PGNDGLNO)
    ! Calculate NDENST from the source SOURCX,SOURCY,SOURCZ
    ! and PNORMX,PNORMY,PNORMZ
    INTEGER TOTELE,NLOC,MLOC,NONODS,PGNODS
    INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    REAL NDENST(NONODS)
    REAL SOURCX(NONODS),SOURCY(NONODS),SOURCZ(NONODS)
    REAL PGRAVX(PGNODS),PGRAVY(PGNODS),PGRAVZ(PGNODS)
    ! Local variables...
    INTEGER ELE,ILOC,KLOC,INOD,PNOD
    REAL DIREC

    ! interpolate normal to pressure nodes...
    DO ELE=1,TOTELE
       DO ILOC=1,NLOC
          KLOC=ILINK2(ILOC,ILOC)
          INOD =NDGLNO((ELE-1)*NLOC+ILOC)
          PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
          DIREC=PGRAVX(PNOD)*SOURCX(INOD)+PGRAVY(PNOD)*SOURCY(INOD) &
               +PGRAVZ(PNOD)*SOURCZ(INOD)
          NDENST(INOD)=SIGN(1.,DIREC)*SQRT(SOURCX(INOD)**2+SOURCY(INOD)**2+SOURCZ(INOD)**2)
       END DO
    END DO

  END SUBROUTINE CALNDENST

  subroutine calug(x,y,z,xold,yold,zold,ug,vg,wg, &
    totele, nloc, ndglno,xondgl,dt)
    ! Calculate the grid velocity UG,VG,WG
    real, dimension(:), intent(in):: x, y, z
    real, dimension(:), intent(in):: xold, yold, zold
    real, dimension(:), intent(out):: ug, vg, wg
    integer, intent(in):: totele, nloc
    integer, dimension(totele*nloc), intent(in):: ndglno, xondgl
    real, intent(in):: dt

    integer nod, xnod, ele, iloc

    if (size(x)==size(ug)) then
       ! assume ug,vg,wg are on the coordinate mesh
       do nod=1, size(x)
          ug(nod)=(x(nod)-xold(nod))/dt
          vg(nod)=(y(nod)-yold(nod))/dt
          wg(nod)=(z(nod)-zold(nod))/dt
       end do
    else
       do ele=1, totele
          do iloc=1, nloc
             nod=ndglno( (ele-1)*nloc+iloc )
             xnod=xondgl( (ele-1)*nloc+iloc )
             ug(nod)=(x(xnod)-xold(xnod))/dt
             vg(nod)=(y(xnod)-yold(xnod))/dt
             wg(nod)=(z(xnod)-zold(xnod))/dt
          end do
       end do
    end if

  end subroutine calug

  SUBROUTINE INTEQETID(EQETID,X,Y,Z,ACCTI2, &
       XNONOD,PGNODS,NLOC,TOTELE,XONDGL, &
       ISPHERE,PGNDGLNO,MLOC,FREESDISOTT, &
       TIDALLCOM,TLOVE,STATID)
    ! This sub calculates the equalibrium tide EQETID by interpolating
    ! the data from a file.
    ! FREESDISOTT is decoded to obtain tidal contituants.
    ! Remember the original value passed down from GEM has been divided by 2**11.
    ! The result is converted to a binary no to switch
    ! on and off the 1st 11 tidal constituants.
    ! if STATID=1 then have a stationary tidal force.
    ! If TIDALLCOM=1 then switch on all components overriding all component options.
    IMPLICIT NONE
    INTEGER XNONOD,PGNODS,NLOC,TOTELE,MLOC
    INTEGER XONDGL(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
    INTEGER FREESDISOTT,TIDALLCOM,TLOVE
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL EQETID(PGNODS)
    REAL ACCTI2,ACCTIM
    INTEGER ISPHERE,STATID
    ! Local variables...
    INTEGER ELE,ILOC,JLOC,KLOC,PNOD,INOD,JNOD
    REAL RADI,RADJ,RADP,RAD
    ! Alocatable memory...
    REAL, ALLOCATABLE, DIMENSION(:)::XP
    REAL, ALLOCATABLE, DIMENSION(:)::YP
    REAL, ALLOCATABLE, DIMENSION(:)::ZP

    ALLOCATE(XP(PGNODS))
    ALLOCATE(YP(PGNODS))
    ALLOCATE(ZP(PGNODS))
    IF(STATID.EQ.0) THEN
       ACCTIM=ACCTI2
    ELSE
       ACCTIM=0.0
    ENDIF

    DO ELE=1,TOTELE
       DO ILOC=1,NLOC
          INOD=XONDGL((ELE-1)*NLOC+ILOC)
          DO JLOC=ILOC,NLOC
             JNOD=XONDGL((ELE-1)*NLOC+JLOC)
             KLOC=ILINK2(ILOC,JLOC)
             PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
             XP(PNOD)=0.5*(X(INOD)+X(JNOD))
             YP(PNOD)=0.5*(Y(INOD)+Y(JNOD))
             ZP(PNOD)=0.5*(Z(INOD)+Z(JNOD))
             IF(ISPHERE.GT.1) THEN
                ! Assume the mid side nodes are on the sphere
                IF(ILOC.NE.JLOC) THEN
                   RADI=SQRT(X(INOD)**2+Y(INOD)**2+Z(INOD)**2)
                   RADJ=SQRT(X(JNOD)**2+Y(JNOD)**2+Z(JNOD)**2)
                   RADP=0.5*(RADI+RADJ)
                   RAD=SQRT(XP(PNOD)**2+YP(PNOD)**2+ZP(PNOD)**2)
                   XP(PNOD)=XP(PNOD)*RADP/RAD
                   YP(PNOD)=YP(PNOD)*RADP/RAD
                   ZP(PNOD)=ZP(PNOD)*RADP/RAD
                ENDIF
             ENDIF
          END DO
       END DO
    END DO

    CALL INTSUBS(XP,YP,ZP,EQETID,PGNODS,ACCTIM,FREESDISOTT,TIDALLCOM,TLOVE)

  END SUBROUTINE INTEQETID

  SUBROUTINE INTSUBS(XP,YP,ZP,EQETID,PGNODS,ACCTIM,FREESDISOTT,TIDALLCOM,TLOVE)
    IMPLICIT NONE

    INTEGER PGNODS,FREESDISOTT,TIDALLCOM,TLOVE
    REAL ACCTIM
    REAL XP(PGNODS),YP(PGNODS),ZP(PGNODS),EQETID(PGNODS)
    ! These are the components to switch on e.g. M2 is 1st component.
    ! FREESDISOTT is decoded to obtain tidal contituants.
    ! FREESDISOTT is converted to a binary no to switch
    ! on and off the 1st 11 tidal constituants.
    ! The DISOTT value associated with the field PSIPRE now controls
    ! the tidal components. To switch on tidal modelling use
    ! DISOTT is decoded to obtain tidal contituants.
    ! (DISOTT/2**11) is converted to a binary no to switch
    ! on and off the 1st 11 tidal constituants.
    ! e.g. 1 switches on the 1st component M2.
    ! 5 is converted to
    ! binary is 101 this thus switches on the 3rd and 1st tidal
    ! components.
    ! IF TLOVE=1 then use a Love number of 0.3 has a factor of (1.-0.3)
    REAL LAT,LONG,RAD,eta,HORIZ_RESCALE,RLOVE
    INTEGER NCHI,pnod,I
    PARAMETER(NCHI=12)
    INTEGER WHICH_TIDE(11)
    REAL CHI(NCHI)
    HORIZ_RESCALE=1.0
    IF(TIDALLCOM.EQ.1) THEN
       DO I=1,11
          WHICH_TIDE(I)=1
       END DO
    ELSE
       WHICH_TIDE(1:11) = 0
       CALL DECODBINTID(WHICH_TIDE,11,FREESDISOTT)
    ENDIF

    CALL FIND_CHI(CHI,NCHI,ACCTIM,HORIZ_RESCALE)
    RLOVE=1.0
    IF(TLOVE.EQ.1) RLOVE=1.0-0.3

    DO PNOD=1,PGNODS
       CALL CONVLATLONGSPHERE(LAT,LONG,RAD,XP(PNOD),YP(PNOD),ZP(PNOD))
       ETA = equilibrium_tide(WHICH_TIDE,LAT,LONG,ACCTIM,HORIZ_RESCALE)
       EQETID(PNOD)=RLOVE*ETA
    END DO

  END SUBROUTINE INTSUBS

  SUBROUTINE CONVLATLONGSPHERE(LAT,LONG,RAD,X,Y,Z)
    IMPLICIT NONE
    REAL LAT,LONG,RAD,X,Y,Z
    RAD=SQRT(X**2+Y**2+Z**2)
    LAT=ATAN2(Z,SQRT(X**2+Y**2))
    LONG=ATAN2(Y,X)

  END SUBROUTINE CONVLATLONGSPHERE

  SUBROUTINE DECODBINTID(WHICH_TIDE,NTID,NODECOD)
    ! Calculate the 11 tidal constituants from NODECOD and
    ! put them in WHICH_TIDE.
    ! E.G. NODECOD=5 switches on 101.
    IMPLICIT NONE
    INTEGER NTID,WHICH_TIDE(NTID),NODECOD
    ! Local variables...
    INTEGER I,NODECOD2,IBIN
    NODECOD2=NODECOD
    DO I=NTID,1,-1
       IBIN=INT( NODECOD2/(2**(I-1)) )
       NODECOD2=NODECOD2 - IBIN*2**(I-1)
       WHICH_TIDE(I)=IBIN
    END DO
    ewrite(2,*) 'WHICH_TIDE:',WHICH_TIDE

  END SUBROUTINE DECODBINTID

  SUBROUTINE DG2CTYNOD(PSIPRE,DGP,TOTELE,NLOC,MLOC,NONODS,NDGLNO)
    ! CALCULATE PSIPRE from DGP for visualisation of DG pressure.
    IMPLICIT NONE
    INTEGER TOTELE,NLOC,MLOC,NONODS
    REAL PSIPRE(NONODS),DGP(MLOC*TOTELE)
    INTEGER NDGLNO(TOTELE*NLOC)
    ! Local variables
    INTEGER DGNOD,CYNOD,NOD,ELE,ILOC,KLOC
    INTEGER, ALLOCATABLE, DIMENSION(:)::ICOUNT
    ALLOCATE(ICOUNT(NONODS))

    ICOUNT(1:NONODS) = 0
    PSIPRE(1:NONODS) = 0.0

    DO ELE=1,TOTELE
       DO ILOC=1,NLOC
          KLOC=ILINK2(ILOC,ILOC)
          DGNOD=(ELE-1)*MLOC+KLOC
          CYNOD=NDGLNO((ELE-1)*NLOC+ILOC)
          ICOUNT(CYNOD)=ICOUNT(CYNOD)+1
          PSIPRE(CYNOD)=PSIPRE(CYNOD)+DGP(DGNOD)
       END DO
    END DO

    DO NOD=1,NONODS
       PSIPRE(NOD)=PSIPRE(NOD)/REAL(ICOUNT(NOD))
    END DO

  END SUBROUTINE DG2CTYNOD

  SUBROUTINE VECXYZ2SOURCXYZ(SOURCX,SOURCY,SOURCZ, VECX,VECY,VECZ, &
       X,Y,Z, NDGLNO,XONDGL, TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI, &
       N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT,&
       ISPHERE,FINDRM,COLM,NCOLM, &
       halo_tag,NNODP)
    ! This sub finds the underlying source SOURCX,SOURCY,SOURCZ
    ! from the discretised source VECX,VECY,VECZ.
    IMPLICIT NONE
    ! This sub calculates the nodes for high order tets of degree ISPHERE
    INTEGER TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,ISPHERE
    REAL SOURCX(NONODS),SOURCY(NONODS),SOURCZ(NONODS)
    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    INTEGER NDGLNO(TOTELE*NLOC),XONDGL(TOTELE*NLOC)
    ! Volume shape functions...
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
    REAL WEIGHT(NGI)
    INTEGER NCOLM
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER halo_tag,NNODP
    ! Local variables...
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI)
    REAL XD(NGI),YD(NGI),ZD(NGI)
    INTEGER ELE,KITS
    INTEGER ILOC,JLOC,GLOBI,GLOBJ,GI,POSMAT
    REAL VLM
    INTEGER IIDUM(1)
    ! *********ALOCAT MEMORY**************
    REAL, ALLOCATABLE, DIMENSION(:)::MASS

    ALLOCATE(MASS(NCOLM))
    MASS(1:NCOLM)=0.0

    DO ELE=1,TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MLOC,NGI, &
            N,NLX,NLY,NLZ, M,MLX,MLY,MLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
            NX,NY,NZ, MX,MY,MZ, &
            A11,A12,A13, A21,A22,A23, A31,A32,A33, &
            XD,YD,ZD, &
            ISPHERE)

       DO ILOC=1,NLOC
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          DO JLOC=1,NLOC
             GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
             VLM=0.0
             DO GI=1,NGI
                VLM=VLM+N(ILOC,GI)*DETWEI(GI)*N(JLOC,GI)
             END DO
             ! Find POSMAT. ***************************************
             CALL POSINMAT(POSMAT,GLOBI,GLOBJ, &
                  NONODS,FINDRM,COLM,NCOLM)
             ! Found POSMAT ***************************************
             MASS(POSMAT)=MASS(POSMAT) + VLM
          END DO
       END DO

    END DO
    ! Solve M SOURCX=VECX  & M SOURCY=VECY   & M SOURCZ=VECZ
    IIDUM(1)=0

    CALL SOLCG(SOURCX, VECX,NONODS,NONODS,NNODP,.TRUE., &
         MASS,FINDRM,COLM,NCOLM, NCOLM,&
         halo_tag, KITS)

    CALL SOLCG(SOURCY, VECY,NONODS,NONODS,NNODP,.TRUE., &
         MASS,FINDRM,COLM,NCOLM, NCOLM,&
         halo_tag, KITS)

    CALL SOLCG(SOURCZ, VECZ,NONODS,NONODS,NNODP,.TRUE., &
         MASS,FINDRM,COLM,NCOLM, NCOLM,&
         halo_tag, KITS)

  END SUBROUTINE VECXYZ2SOURCXYZ

  subroutine shortpgassem1free(shortdifmath, &
       shortdifvecx, &
       WETDRYMAT,SUFML, &
       topdis,botdis,atheta,absorb, &
       frees,freesold,gravty,dt,ctheta,velmean, &
       xorig,yorig,zorig, u,v,w,&
       INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
       nloc,mloc, &
       totele, &
       geobal,&
       ndglno,pgndglno, xondgl, &
       pgnods,nonods,xnonod,scfacth0,isphere, &
       findrm,colm,ncolm, &
       fincmc,midcmc,colcmc,ncmc,fredop, &
       nofsta, &
       cograx, &
       fstheta, rnonl_coefficient, wetdry,colele4,&
       wdgravity)
    ! put the result into momentum eqns vecx=vexc+n_i dp/dx etc
    implicit none
    logical usect,nonlin,sufint
    ! if usect then put c1t, c2t,c3t into eqn else do not.
    ! if nonlin then use nonlinear dissipation.
    parameter(usect=.true.,nonlin=.true.,sufint=.true.)
    integer nonods,pgnods,xnonod,isphere
    integer totele,nloc,mloc
    integer ncolm
    integer findrm(nonods+1),colm(ncolm)
    integer fredop,ncmc,fincmc(fredop+1),colcmc(ncmc),midcmc(fredop)
    real gravty,dt,ctheta,atheta
    integer velmean
    real topdis(xnonod),botdis(xnonod)
    real absorb(nonods)
    real shortdifmath(ncmc)
    real shortdifvecx(nonods)
   
    ! if takvert then take out the vertical component of pressure/forcing.
    real frees(pgnods),freesold(pgnods)
    real xorig(xnonod),yorig(xnonod),zorig(xnonod)
    real u(nonods),v(nonods),w(nonods)
    INTEGER INUSUB,NSUBNVLOC
    REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
    REAL NWSUB(TOTELE*NSUBNVLOC)
    integer ndglno(totele*nloc),pgndglno(totele*mloc)
    integer xondgl(totele*nloc)

    integer geobal,nofsta
    ! this is the theta time stepping parameter for the free surface.
    real fstheta
    real, intent(in):: rnonl_coefficient
    integer wetdry
    integer COLELE4(TOTELE*5)
    REAL WETDRYMAT(NCMC*WETDRY),SUFML(NONODS)


    real scfacth0
    logical cograx

    Real, intent(in) :: wdgravity(totele)

    call alltopsufdg(shortdifmath,shortdifvecx, &
         wetdrymat,sufml, &
         xorig,yorig,zorig, &
         topdis,botdis,atheta,absorb, &
         gravty,ctheta, &
         nloc, &
         totele,&
         geobal,&
         ndglno,xondgl, &
         nonods,xnonod, &
         scfacth0,isphere, &
         findrm,colm,ncolm, &
         fincmc,midcmc,colcmc,ncmc,fredop, &
         cograx, &
         u,v,w, &
         INUSUB,NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
         nofsta,nonlin,velmean,dt, &
         pgnods,mloc,pgndglno,frees,freesold, &
         fstheta,rnonl_coefficient,wetdry,colele4,wdgravity)

  end subroutine shortpgassem1free

  subroutine mapgrav(pgravx,pgravy,pgravz,adjgrav, &
       potgra,usepotgra,  &
       pnormx,pnormy,pnormz, &
       findrm,colm,ncolm,centrm, &
       pgfindrm,pgcolm,npgcolm,pgnods,pgcentrm, &
       x,y,z, &
       m,mlx,mly,mlz, &
       n,nlx,nly,nlz, weight, nloc,ngi,mloc, &
       totele,d3,dcyl,  &
       ndglno,pgndglno, xondgl, &
       nonods,xnonod, &
       halotagt10,para,procno,pgnodsp,&
       getfrees,nlevel,isphere, option_path)
    ! this sub calculates gravity (pgravx,pgravy,pgravz) taking
    ! into account the mesh from an initial guess of
    ! gravity (vgravx,vgravy,vgravz )
    ! if adjgrav=1 then adjust gravity direction/magnetude for mesh
    ! if adjgrav=0 then dont adjust gravity.
    ! if usepotgra=1 then get the gravitational potential and use that.
    implicit none
    integer adjgrav,usepotgra
    integer nonods,pgnods,xnonod
    integer totele,nloc,mloc,ngi
    integer ncolm,npgcolm
    integer centrm(nonods),findrm(nonods+1),colm(ncolm)
    integer pgfindrm(pgnods+1),pgcolm(npgcolm)
    integer pgcentrm(pgnods)
    real pgravx(pgnods),pgravy(pgnods),pgravz(pgnods)
    real pnormx(pgnods),pnormy(pgnods),pnormz(pgnods)
    real potgra(pgnods*usepotgra)
    real x(xnonod),y(xnonod),z(xnonod)
    integer ndglno(totele*nloc),pgndglno(totele*mloc)
    integer xondgl(totele*nloc)

    real m(mloc,ngi),mlx(mloc,ngi),mly(mloc,ngi),mlz(mloc,ngi)
    real n(nloc,ngi),nlx(nloc,ngi),nly(nloc,ngi),nlz(nloc,ngi)
    real weight(ngi)
    logical d3,dcyl
    integer halotagt10,para,procno,pgnodsp
    logical getfrees
    integer nlevel,isphere
    character(len=*), intent(in):: option_path

    ! local variables...
    integer iidum(1),kits,i,ppset
    logical only3d
    real, allocatable, dimension(:)::difmath
    real, allocatable, dimension(:)::vecgra1
    real, allocatable, dimension(:)::vecgra2
    real, allocatable, dimension(:)::vecgra3
    real, allocatable, dimension(:)::p
    real, allocatable, dimension(:)::work1
    real, allocatable, dimension(:)::work2
    real, allocatable, dimension(:)::work3
    real, allocatable, dimension(:)::work4
    real, allocatable, dimension(:)::work5
    real, allocatable, dimension(:)::work6

    ! initial guess for gravity direction/magnetude
    do i=1,pgnods
       pgravx(i)=-pnormx(i)
       pgravy(i)=-pnormy(i)
       pgravz(i)=-pnormz(i)
    end do
    if(adjgrav.eq.0) return

    allocate(difmath(npgcolm))
    allocate(vecgra1(pgnods))
    allocate(vecgra2(pgnods))
    allocate(vecgra3(pgnods))
    allocate(p(pgnods))
    allocate(work1(pgnods))
    allocate(work2(pgnods))
    allocate(work3(pgnods))
    allocate(work4(pgnods))
    allocate(work5(pgnods))
    allocate(work6(4*pgnods))

    iidum=0

    call mapgravdis(1,difmath,p, &
         vecgra1,vecgra2,vecgra3, &
         pgravx,pgravy,pgravz, &
         pgfindrm,pgcolm,npgcolm,pgnods, &
         x,y,z, &
         m,mlx,mly,mlz,&
         n,nlx,nly,nlz, weight, nloc,ngi,mloc, &
         totele,d3,dcyl,  &
         pgndglno, xondgl, &
         xnonod,isphere)

    ! the pressure level to set to 0.0
    if((para.eq.0).or.(procno.eq.1)) then
       ppset=1
    else
       ppset=0
    endif
    if(getfrees) then
       ! apply b.c's to the surface
       call topbcdifdia(difmath, &
            pgfindrm,pgcolm,pgcentrm,npgcolm,pgnods,nonods, &
            totele,nloc,mloc,ndglno,pgndglno,vecgra1,nlevel,.false.)
       ppset=0
    else
       ! apply b.c=0 at pnod=ppset
       call dissetmaval(difmath,npgcolm,pgnods,pgfindrm,pgcolm,pgcentrm, &
            ppset,vecgra1)
    endif

    ! set the b.c to zero at the surface in difmatv *******
    if((nlevel.lt.2).or.getfrees) then
       only3d=.true.
    else
       only3d=.false.
    endif

    p(1:pgnods) = 0.0
    ! lump free surface eqn in vertical and solve...
    call mulgridvert(difmath,vecgra1, &
         findrm,colm,centrm,ncolm,nonods, &
         pgfindrm,pgcolm,pgcentrm,npgcolm,pgnods, &
         p,&
         totele,nloc,mloc, &
         ndglno,pgndglno, x,y,z,xnonod,isphere, &
         nlevel, &
         pgnodsp, &
         halotagt10,0,  &
         .true.,only3d,.false.,.true., &
         option_path)

    ! find the pressure node ppset to set to zero...
    !       call setppset(difmath,npgcolm,pgnods,pgcentrm,infiny,ppset)
    if(usepotgra.eq.1) potgra(1:pgnods) = p(1:pgnods)


    call mapgravdis(2,difmath,p, &
         vecgra1,vecgra2,vecgra3, &
         pgravx,pgravy,pgravz, &
         pgfindrm,pgcolm,npgcolm,pgnods, &
         x,y,z, &
         m,mlx,mly,mlz, &
         n,nlx,nly,nlz, weight, nloc,ngi,mloc, &
         totele,d3,dcyl,  &
         pgndglno, xondgl, &
         xnonod,isphere)

    ! re-calculate gravity by solving a 3 mass eqns
    call solcg(pgravx, vecgra1,pgnods,pgnods,pgnodsp,.true.,  &
         difmath,pgfindrm,pgcolm,npgcolm, npgcolm,&
         halotagt10, kits)

    call solcg(pgravy, vecgra2,pgnods,pgnods,pgnodsp,.true.,  &
         difmath,pgfindrm,pgcolm,npgcolm, npgcolm,&
         halotagt10, kits)

    call solcg(pgravz, vecgra3,pgnods,pgnods,pgnodsp,.true., &
         difmath,pgfindrm,pgcolm, npgcolm, npgcolm, &
         halotagt10, kits)

  end subroutine mapgrav



  SUBROUTINE ABSWETDRY(XABSOR,YABSOR,ZABSOR,NONODS,XNONOD, &
       NCOLM,FINDRM,COLM,TOPDIS,BOTDIS,nodfrees,shortdp,DT,&
       coef_abs,ee00,dd00,totele,nloc,ndglno, xondgl, wdgravity)
    ! This sub cal;culates the absorbtion for wetting and drying...
    INTEGER,intent(in) ::  totele,NONODS,XNONOD,NCOLM
    REAL,intent(out) ::  XABSOR(NONODS),YABSOR(NONODS),ZABSOR(NONODS)
    REAL,intent(in) ::TOPDIS(XNONOD),BOTDIS(XNONOD)
    INTEGER,intent(in) ::  FINDRM(NONODS+1),COLM(NCOLM)
    Real, intent(in) :: nodfrees(nonods), shortdp(nonods)
    Integer,intent(in) :: ndglno(totele*nloc)
    integer, intent(in) :: xondgl(totele*nloc)
    Real, intent(in) :: dt,coef_abs, ee00, dd00
    REAL,intent(out) :: wdgravity(totele)

    ! Local variables...
    INTEGER NOD,COUNT,COL
    Integer :: ele, iloc, nloc, inod, ixnod
    Real :: depthtmp
    REAL RDEPTH,EXDEP,depthcol

    Real :: btlevelmax, btlevelmin
    Integer :: nbtmax, nbtmin, nbtmaxx, nbtminx

    !    real :: coef_abs=10000.
    !    real :: coef_abs=1.e5
    !    real :: ee00=16.6095    ! d=2d_0, xabsor=0.
    !    real :: ee00=9.2103    ! d=2d_0, xabsor=10.
    !    real :: ee00=6.908    ! d=2d_0, xabsor=100.

    !    real :: coef_abs=100.
    !    integer :: nalpha=2


    !   Calculate Wdgravity

    wdgravity=1.0
    Do ele=1,totele
       btlevelmax=-huge(dd00)
       btlevelmin= huge(dd00)
       
       nbtmax=ndglno((ele-1)*nloc+1)
       nbtminx=xondgl((ele-1)*nloc+1)
       nbtmin=ndglno((ele-1)*nloc+1)
       nbtmaxx=xondgl((ele-1)*nloc+1)
       
       Do iloc=1,nloc
          inod=ndglno((ele-1)*nloc+iloc)
          ixnod=xondgl((ele-1)*nloc+iloc)
          
          depthtmp=topdis(ixnod)+botdis(ixnod)+shortdp(inod)
          If(nodfrees(inod)-depthtmp.gt.btlevelmax) then
             nbtmax=inod; nbtmaxx=ixnod
          end if
          If(nodfrees(inod)-depthtmp.lt.btlevelmin) then
             nbtmin=inod; nbtminx=ixnod
          end if
       End do
       If(topdis(nbtmaxx)+botdis(nbtmaxx)+shortdp(nbtmax).le.1.1*dd00) wdgravity(ele)=0.
       If(nodfrees(nbtmin)-nodfrees(nbtmax).gt.dd00) wdgravity(ele)=1.
       !            if(wdgravity(ele).lt.0.5) then
       !              ewrite(3,*) 'ele, wdgravity, depth_btmax--->', &
       !                      & ele, wdgravity(ele),topdis(nbtmax)+botdis(nbtmax)+shortdp(nbtmax) 
       !            end if
    End do


    !         wdgravity=1.0

    !        Calculate Absorption


    if (xnonod/=nonods) then
       FLAbort("Wetting and drying does not work in periodic domains")
    end if
    
    do nod=1,nonods
       rdepth=topdis(nod)+botdis(nod)+shortdp(nod)
       depthcol=-huge(rdepth)
       DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
          COL=COLM(COUNT)
          depthcol=max(depthcol,topdis(col)+botdis(col)+shortdp(col))
          !               if(nodfrees(col).gt.nodfrees(nod)-1.0*dd00) then
          if(nodfrees(col)+shortdp(col).gt.nodfrees(nod)+shortdp(nod)) then
             rdepth=min(rdepth,topdis(col)+botdis(col)+shortdp(col))
          end if
       end do

       if(depthcol.le.dd00) rdepth=dd00

       !            if(rdepth.le.0.1*dd00) then
       !               write(0, *) 'nod, dd00, rdepth==', nod, dd00, rdepth
       !               FLAbort('zero or negative depth')
       !            end if

       EXDEP=EXPTOL(-EE00*(RDEPTH-DD00)/DD00)


       XABSOR(nod)=coef_abs*EXDEP/DT
       YABSOR(nod)=coef_abs*EXDEP/DT
       ZABSOR(nod)=coef_abs*EXDEP/DT

       !             XABSOR(nod)=alpha0/dt/rdepth
       !             YABSOR(nod)=alpha0/dt/rdepth
       !             ZABSOR(nod)=alpha0/dt/rdepth

       !             XABSOR(nod)=alpha1/(dt*(rdepth**2))
       !             YABSOR(nod)=alpha1/(dt*(rdepth**2))
       !             ZABSOR(nod)=alpha1/(dt*(rdepth**2))

       !!            XABSOR(nod)=100./dt*(dd00/rdepth)**nalpha
       !!            YABSOR(nod)=100./dt*(dd00/rdepth)**nalpha
       !!            YABSOR(nod)=100./dt*(dd00/rdepth)**nalpha

    end do

  END SUBROUTINE ABSWETDRY




  SUBROUTINE alterfrees(wetdrymat,sufml,topdis,botdis,findrm, colm, &
       centrm, ncolm, xorig, yorig, zorig, nodfrees, &
       nonods, nlevel, dd00)

    ! This sub adjust water depth in dry area

    integer, intent(in) :: ncolm, nlevel, nonods
    real, intent(inout) :: wetdrymat(ncolm), sufml(nonods)
    integer, intent(in) :: findrm(nonods+1), colm(ncolm), centrm(nonods)
    real, intent(in) :: xorig(nonods), yorig(nonods), zorig(nonods), dd00
    !! This should be xnonod:
    real, intent(inout) :: topdis(nonods), botdis(nonods)
    real, intent(inout) :: nodfrees(nonods)

    integer :: count, nod, icent, ilevel
    integer :: i, j
    integer :: nodbot, nodtop
    real :: deltdepthmax, deltdepth, originaldepth, originaltonode

    real, allocatable :: depthex(:), depthnew(:), depthold(:)

    integer :: nloop1=2000

    !! there are various places where it is assumed xnonod=nonods
    !! this should be fixed
    allocate(depthex(nonods), depthnew(nonods), depthold(nonods))

    do count=1,ncolm
       wetdrymat(count)=-wetdrymat(count)
    end do

    do nod=1,nonods
       icent=centrm(nod)
       wetdrymat(icent)=wetdrymat(icent)+sufml(nod)
    end do

    depthold=topdis+botdis
    depthnew=depthold

    !! nlevel dependency should be removed!!!
    do j=1, nloop1
       ewrite(2, *) 'depth adjustment it: j=',j
       do nod=nonods-(nonods/nlevel)+1, nonods
          depthold(nod)=depthnew(nod)
          depthex(nod)=min(depthold(nod),1.0*dd00)-dd00
       end do

       deltdepthmax=0.
       do nod=nonods-(nonods/nlevel)+1, nonods
          deltdepth=0.
          do i=findrm(nod),findrm(nod+1)-1
             deltdepth=deltdepth+wetdrymat(i)*depthex(colm(i))
             !             if(depthold(nod).lt.dd00) then
             !               ewrite(3,*) 'nod, i, colm(i), depthold(nod)==',nod, i, colm(i), depthold(nod)
             !               ewrite(3,*) 'wetdrymat(i), sufml(nod)==',wetdrymat(i), sufml(nod)
             !               ewrite(3,*) 'depthex(colm(i)), deltdepth===', depthex(colm(i)), deltdepth
             !               pause
             !             end if
          end do

          if(abs(sufml(nod)).le.1.e-8) then
             !            depthnew(nod)=depthold(nod)
             FLAbort("sufml le 1.e-8") 
          else
             depthnew(nod)=depthold(nod)-deltdepth/sufml(nod)
          endif
          if(abs(deltdepth/sufml(nod)).gt.abs(deltdepthmax)) then
             deltdepthmax=deltdepth/sufml(nod)
          end if
       end do

       ewrite(2,*) 'min depthnew max deltdepth==',&
            minval(depthnew(nonods-(nonods/nlevel)+1:nonods)),deltdepthmax
       if(minval(depthnew(nonods-(nonods/nlevel)+1:nonods)).ge.0.95*dd00) exit

    end do


    if(minval(depthnew(nonods-(nonods/nlevel)+1:nonods)).le.dd00*1.e-5) then
       !           FLAbort("wetting-drying depth diffusion failed")
       ewrite(0,*) "Caution: wetting-drying depth diffusion failed in alterfrees"
    end if

    do ilevel=1,nlevel
       do nodbot=1,nonods/nlevel
          nod=(ilevel-1)*nonods/nlevel+nodbot
          nodtop=(nlevel-1)*nonods/nlevel+nodbot
          depthnew(nod)=depthnew(nodtop)
       end do
    end do

    do ilevel=1,nlevel
       do nodbot=1,nonods/nlevel
          nod=(ilevel-1)*nonods/nlevel+nodbot
          nodtop=(nlevel-1)*nonods/nlevel+nodbot
          originaldepth=sqrt(&
               (xorig(nodtop)-xorig(nodbot))**2+&
               (yorig(nodtop)-yorig(nodbot))**2+&
               (zorig(nodtop)-zorig(nodbot))**2)
          originaltonode=sqrt(&
               (xorig(nod)-xorig(nodbot))**2+&
               (yorig(nod)-yorig(nodbot))**2+&
               (zorig(nod)-zorig(nodbot))**2)
          botdis(nod)=depthnew(nod)*originaltonode/originaldepth
          topdis(nod)=depthnew(nod)*(originaldepth-originaltonode)/originaldepth
       end do
    end do

    do ilevel=1,nlevel
       do nodbot=1,nonods/nlevel
          nod=(ilevel-1)*nonods/nlevel+nodbot
          nodtop=(nlevel-1)*nonods/nlevel+nodbot
          nodfrees(nod)=depthnew(nod)-sqrt(&
               (xorig(nodtop)-xorig(nodbot))**2+&
               (yorig(nodtop)-yorig(nodbot))**2+&
               (zorig(nodtop)-zorig(nodbot))**2)
       end do
    end do

    deallocate(depthex, depthnew, depthold)

  END SUBROUTINE alterfrees

  subroutine geostrophic_check_options

    if (have_option('/material_phase[0]/scalar_field::FreeSurface')) then

       if (option_count('/material_phase')/=1) then

          FLExit("Free surface doesn't work for multi-material/phase")

       end if

       if (have_option('/material_phase[0]/scalar_field::BalancePressure')) then

          FLExit("Can't have both a balance pressure and a free surface field.")

       end if

       if (have_option('/material_phase[0]/scalar_field::Pressure/&
            &prognostic/spatial_discretisation/continous_galerkin')) then

          if (.not. have_option('/material_phase[0]/scalar_field::Pressure/&
               &prognostic/spatial_discretisation/continous_galerkin/remove_stabilisation_term')) then

             ewrite(-1,*) 'Missing option spatial_discretisation/continous_galerkin/remove_stabilisation_term'
             ewrite(-1,*) 'under the prognostic pressure field'
             FLExit("With a free surface you have to switch off the pressure stabilisation filter.")
          end if
       end if

    end if

  end subroutine geostrophic_check_options

end module geostrophic
