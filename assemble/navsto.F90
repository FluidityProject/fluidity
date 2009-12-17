!     Copyright (C) 2006 Imperial College London and others.
!     
!     Please see the AUTHORS file in the main source directory for a full list
!     of copyright holders.
!     
!     Prof. C Pain
!     Applied Modelling and Computation Group
!     Department of Earth Science and Engineering
!     Imperial College London
!     
!     C.Pain@Imperial.ac.uk
!     
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation,
!     version 2.1 of the License.
!     
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!     Lesser General Public License for more details.
!     
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!     USA

#include "fdebug.h"

      module navsto_module

      use quadrature
      use elements
      use AdvectionDiffusion
      use AllSorts
      use flcomms_module
      use FLDebug
      use geostrophic
      use OceanSurfaceForcing
      use sparse_tools
      use fields
      use state_module
      use module_solnav
      use solmom_module
      use global_parameters, only:ISPHERE, OPTION_PATH_LEN, phase2state_index, my_new_options => new_options
      use diff3d_module
      use assemble_cmc
      use compressible_projection, only: assemble_compressible_projection
      use spud
      use drag_module
      use shape_module
      use EletemFeinte_module
      use assnav_module
      use assemble_boundary_conditions
      use futils
      use position_in_matrix
      use rotated_boundary_conditions_legacy
      use shear_stress
      use divergence_matrix_cv
      use detnlxr_module
      use ctmult_module
      use roctrt_module
      use spaerr_module
      use legacy_tensor_fields
      use legacy_boundary_conditions
      use boundary_conditions
      use sparsity_patterns_meshes
      use state_matrices_module
      use Tidal_module
      implicit none
      
      private

      public navsto, raddin,bigpha,mardmi,free_surface_diagnostics,harmonic_analysis_at_single_node

      contains

      SUBROUTINE NAVSTO&

     &     ( &
     !     --------------------------------------START OF ADD BY CRGW 31/03/06
     !          PLASTIC STRAINS:
     &     STVPXX,STVPYY,STVPZZ,&
     &     STVPYZ,STVPXZ,STVPXY,&
     !     --------------------------------------END OF ADD BY CRGW 31/03/06
     &     IGUESS, &
     &     FREESDISOTT,&
     &     NOBCFSimp,BCT1WFSimp,BCT2FSimp,&
     &     R0,D0, &
     &     MVMESH,&
     &     FREDOP,TOTELE, &
     &     NONODS,&
     &     NLOC,NGI,MLOC,DT,&
     &     DISOPT,DISOPN,THETA,BETA,LUMP,MAKSYM,&
     &     CMCGET,NDPSET,GETC12,&
     &     PREOPT,CONVIS,&
     &     SNONOD,VNONOD, &
     &     XNONOD, &
     &     DCYL,&
     &     NOBCU,BCU1,BCU2,BCU1VAL,&
     &     NOBCV,BCV1,BCV2,BCV1VAL,&
     &     NOBCW,BCW1,BCW2,BCW1VAL,&
     &     NDWISP,NCMCB,NPRESS,&
     &     NPROPT,RADISO,ALFST2,SPRESS,&
     &     ABSLUM,SLUMP,&
     &     UZAWA, MULPA, CGSOLQ,GMRESQ,&
     &     MIXMAS,POISON,PROJEC,&
     &     halo_tag,&
     &     halo_tag_p,NNODP,NNODPP,&
!     The following is for rotations only ROTAT=.TRUE.If GETTAN find tangents
     &     GETTAN, &
     &     ROTAT,NNODRO,NRTDR, &
!     New stuff for ASSNAV...
     &     DSPH, &
     &     GEOBAL,&
     &     STOTEL,SNLOC,SNGI,      &
     &     GRAVTY,&
!     The multiphase stuff...
     &     NPHASE,IPHASE,&
!     Some extra stuff MDP added for the free surface 
     &     NEWMES2, &
     &     SCFACTH0,&
     &     ITSITI,ITINOI,&
     &     GOTTRAF,&
!     sub element modelling...
     &     NSUBVLOC,NSUBNVLOC,&
!     For free surface and adjoint
     &     GOTTOPDIS, geostrophic_solver_option, nlevel, &
     &     GOTBOY,BSOUX,BSOUY,BSOUZ, &
     &     COGRAX, &
     &     VERSIO,ISPHER2,MISNIT,&
     &     state, &
!     --------------------------------------START OF ADD BY CRGW 14/03/06
!          INTEGER SOLIDITY FLAGS:     
     &     SOLIDS,&
     !     --------------------------------------END OF ADD BY CRGW 14/03/06
     !     --------------------------------------START OF ADD BY CRGW 20/03/06
     !             COMPRESSIBILITY MATRICES:
     &     MKCOMP,&
     ! for LES -Jemma
     &     metric_tensor&
     &     )

     !     --------------------------------------END OF ADD BY CRGW 20/03/06

!     branch test
!     See subroutine SOLNAV for definitions of the following from DP onwards
!     THIS SUBROUTINE ASSEBLES AND SOLVES THE NAVIER STOKES EQUATION.
!     If PLEVEL then set the pressure level. 
!     If SYM then force the discretised advec diff to be symmetric
!     PSIPRE is a pressure like variable used to deal with 
!     corriolis. 
!     
!     BCU2,BCV2,BCW2 contain the nodes at which the boundary 
!     conditions are to be applied. 
!     BCU1 caontains  the value of the specified acceleration of U 
!     that is BCU1=(Unew - Uold)/DT
!     BCU1VAL contains the actually value.
!     BIGM2 is only used in cylinderical coords
!     otherwise it can contain BIGM1. 
!     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
!     ADDSOX,ADDSOY,ADDSOZ are vec's we 
!     add to the discretised rhs of momentum equtions. 
      LOGICAL INMAT,NAV
      LOGICAL TRANSP,NTRANS
      INTEGER NBIG12,NLOC,NGI,NNODP,NNODPP
      INTEGER NRTDR,NOFILT
      INTEGER MULGRID
      REAL R0
      REAL TOLER,acctim
      PARAMETER(TRANSP=.FALSE.,NTRANS=.TRUE.)
      PARAMETER(INMAT=.TRUE.,NAV=.TRUE.,TOLER=0.0001)
!     NOFILT=1 switches off the pressure 4th order filter...
      PARAMETER(MULGRID=0)
      INTEGER DISOPT,DISOPN,GEOBAL
      INTEGER TOTELE,NONODS,SNONOD,VNONOD,XNONOD,FREDOP
      INTEGER NCMCB,NPRESS,NPROPT,RADISO
      REAL ALFST2,SPRESS
      INTEGER NOBCU,NOBCV,NOBCW
      INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
      REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
      REAL BCU1VAL(NOBCU),BCV1VAL(NOBCV),BCW1VAL(NOBCW)
!     
      INTEGER MLOC,NPHASE,NCOLM,NCT
      REAL DT,THETA,BETA
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
!     If NDWISP then we have the same pressure nodes as velocity nodes. 
      INTEGER NOBCFSimp
      REAL BCT1WFSimp(NOBCFSimp)
      INTEGER BCT2FSimp(NOBCFSimp)
      LOGICAL GOTPSI,GOT_EQTD
      LOGICAL GOT_SHEARSTRESSX,GOT_SHEARSTRESSY,GOT_SHEARSTRESSZ
      LOGICAL GOT_MEANSHEARSTRESSX,GOT_MEANSHEARSTRESSY,GOT_MEANSHEARSTRESSZ
      LOGICAL IGUESS
      INTEGER FREESDISOTT
      LOGICAL MVMESH
      real, dimension(:), allocatable ::  l_soxgra, l_soygra, l_sozgra
      INTEGER PREOPT
!     sub element modelling...
      INTEGER NSUBVLOC,NSUBNVLOC
      INTEGER NSOGRASUB

      LOGICAL GETC12
      LOGICAL LUMP,MAKSYM
      LOGICAL D3,DCYL,CONVIS 
      LOGICAL CMCGET
!     NDPSET contains the pressure nodee to be set to zero. 
!     If NDPSET=0 then no pressure node is set to zero.
      INTEGER NDPSET
!     THE FOLLOWING ARE FOR SOLNAV...
      INTEGER UZAWA,CONMAS,COUPLE,MIXMAS
      INTEGER D2HALF,MULPA,CGSOLQ,GMRESQ
      INTEGER POISON,PROJEC
      INTEGER NCMC
      LOGICAL NDWISP

!     For parallel solution...
      INTEGER  PARA,halo_tag,halo_tag_p
!     NB The normal is pointing out of the domain. 
!     The rotation matrix in 3-D is R=  
!     T1X           T1Y          T1Z
!     T2X           T2Y          T2Z
!     NORMX    NORMY    NORMZ
!     The rotation matrix in 2-D is R=  
!     T1X           T1Y   
!     NORMX    NORMY    
!     
      INTEGER NNODRO
      LOGICAL ROTAT
      LOGICAL GETTAN
!     New stuff for ASSNAV...
      LOGICAL DSPH
      INTEGER STOTEL,SNLOC,SNGI
      REAL D0 
      LOGICAL ROTMOM
      REAL     GRAVTY
!     Used in free surface routine    
      INTEGER  ITSITI,ITINOI         
!     EXTRA LOGICALS AND POINTERS
      LOGICAL  NEWMES2
      
      real :: scfacth0

!     RTDR contains RT D R   ^T where D contains DMI's and R is 
!     the rotation matrix. 
!     
!     For multiphase stuff ONLY...
!     DM1=DM1PHA(I+NONODS*(IPHASE-1)) contains b.c's of IPHASE.
!     DM2=DM2PHA(I+NONODS*(IPHASE-1)) contains b.c's of IPHASE.
!     DM3=DM3PHA(I+NONODS*(IPHASE-1)) contains b.c's of IPHASE.
!     ML=MLPHA(I+NONODS*(IPHASE-1)) contains lumped mass*density of 
!     current phase IPHASE.
!     MPDIAG contains the diagonals of the block matrix used instead 
!     of the mass matrix for the non-sym multi-phase projection method.
!     GLPHAU,GLPHAV,GLPHAW are the overall components of fluid velocity. 
!     NUPHA,NVPHA,NWPHA are the vel components of each phase to be used 
!     for advection calculations in ADVDIF or here for momentum. 
!     PHSOUX,PHSOUY,PHSOUZ are the phasic sources. 
!     NDVOLF is the volume fraction of all phases (nod-based)
!     C1TP,C2TP,C3TP contains the multi-phase r.h.d matrices acting 
!     on the volume fractions. 
      INTEGER IPHASE
      LOGICAL ABSLUM,SLUMP
      real, dimension(0):: real_length_zero_array
      real, allocatable, dimension(:)::SX,SY,SZ
      real, allocatable, dimension(:), target:: ZER
      INTEGER NOD,IDIM,NOD2,ID3,I

      !     INTEGER SOLIDITY FLAGS:
      INTEGER SOLIDS
      !     --------------------------------------END OF ADD BY CRGW 14/03/06
      !     --------------------------------------START OF ADD BY CRGW 31/03/06
      !     PLASTIC STRAINS:
      REAL STVPXX(:), STVPYY(:), STVPZZ(:)
      REAL STVPYZ(:), STVPXZ(:), STVPXY(:)
      !     --------------------------------------END OF ADD BY CRGW 31/03/06
      !     --------------------------------------START OF ADD BY CRGW 20/03/06
      !     COMPRESSIBILITY MATRICES:
      INTEGER MKCOMP
      
!**** ADJOINT MODEL***************    
      LOGICAL, intent(in):: GOTTOPDIS
      integer, intent(in):: geostrophic_solver_option, nlevel
      LOGICAL GOTBOY
      REAL BSOUX,BSOUY,BSOUZ
      LOGICAL COGRAX
      INTEGER VERSIO,ISPHER2,MISNIT
      REAL VOLELE
      INTEGER ISUBSOU,BIG_NLOC
      LOGICAL SUF_TOP_OCEAN,NEXTTIM
      INTEGER INUSUB,ISUB,NBUOY,IUNST_FREE
      INTEGER NOD1,NOD3,NOD4,ELE
      REAL, ALLOCATABLE, DIMENSION(:)::PSIPRETEM
      REAL, ALLOCATABLE, DIMENSION(:)::UOLD2
      REAL, ALLOCATABLE, DIMENSION(:)::VOLD2
      REAL, ALLOCATABLE, DIMENSION(:)::WOLD2
      REAL, ALLOCATABLE, DIMENSION(:)::VECXMOM
      REAL, ALLOCATABLE, DIMENSION(:)::VECYMOM
      REAL, ALLOCATABLE, DIMENSION(:)::VECZMOM
      REAL, ALLOCATABLE, DIMENSION(:)::XABSOR
      REAL, ALLOCATABLE, DIMENSION(:)::YABSOR
      REAL, ALLOCATABLE, DIMENSION(:)::ZABSOR
      REAL, ALLOCATABLE, DIMENSION(:)::VECX_SGSADD
      REAL, ALLOCATABLE, DIMENSION(:)::VECY_SGSADD
      REAL, ALLOCATABLE, DIMENSION(:)::VECZ_SGSADD
      REAL, ALLOCATABLE, DIMENSION(:)::UNEW
      REAL, ALLOCATABLE, DIMENSION(:)::VNEW
      REAL, ALLOCATABLE, DIMENSION(:)::WNEW
      REAL, ALLOCATABLE, DIMENSION(:)::UNEWSUB
      REAL, ALLOCATABLE, DIMENSION(:)::VNEWSUB
      REAL, ALLOCATABLE, DIMENSION(:)::WNEWSUB
      REAL, ALLOCATABLE, DIMENSION(:)::TXABSOR
      REAL, ALLOCATABLE, DIMENSION(:)::tempML
      REAL, ALLOCATABLE, DIMENSION(:)::BOY_ML
      
      real, dimension(:), pointer :: nu,nv,nw
      real, dimension(:), pointer ::  u, v, w
      real, dimension(:), pointer ::  ug, vg, wg
      real, dimension(:), pointer ::  sourcx, sourcy, sourcz
      real, dimension(:), pointer ::  soxgra, soygra, sozgra
      real, dimension(:), pointer ::  soxgrasub, soygrasub, sozgrasub
      real, dimension(:), pointer :: psipre
      real, dimension(:), pointer ::  x, y, z
      real, dimension(:), pointer ::  xold, yold, zold
      real, dimension(:), pointer ::  xorig, yorig, zorig
      real, dimension(:), pointer ::  vecx, vecy, vecz
      real, dimension(:), pointer ::  c1t, c2t, c3t
      real, dimension(:), pointer ::  c1tp, c2tp, c3tp
      real, dimension(:), pointer :: cmc, kcmc
      real, dimension(:), pointer :: ml, denpt, p, divqs
      real, dimension(:), pointer :: vgravx, vgravy, vgravz
      real, dimension(:), pointer :: eqtd
      real, dimension(:), pointer :: shearstressx, shearstressy, shearstressz
      real, dimension(:), pointer :: meanshearstressx, meanshearstressy, meanshearstressz
      real, dimension(:), pointer :: usub, vsub, wsub
      real, dimension(:), pointer :: nusub, nvsub, nwsub
      
      integer, dimension(:), pointer :: pndgln, ndglno, sondgl, vondgl, xondgl
      integer, dimension(:), pointer :: findrm, colm, centrm
      integer, dimension(:), pointer :: findct, colct
      integer, dimension(:), pointer :: fincmc, colcmc, midcmc
      
      real, dimension(:), allocatable ::  bigm1
      real, dimension(:), allocatable :: dm1, dm2, dm3
      real, dimension(:), allocatable :: dp
      real, dimension(:), allocatable :: rhs
      real, dimension(:), allocatable :: normx, normy, normz
      real, dimension(:), allocatable :: t1x, t1y, t1z
      real, dimension(:), allocatable :: t2x, t2y, t2z
      real, dimension(:), allocatable :: rtdr
     
      integer, dimension(:), allocatable :: nodrot, sndgln
      
      integer stat
      
!     Encapsulated state information containing fields. -dham
      type(state_type), dimension(:), intent(inout) :: state

!     System matrices (extracted from state)
      type(block_csr_matrix), pointer :: CT_m
      type(block_csr_matrix) :: CTP_m
      type(csr_matrix), pointer :: CMC_m, kmk_m

      ! this is now local to navsto, i.e. reallocated every call
      type(block_csr_matrix) :: big_m
      type(csr_sparsity), pointer :: u_sparsity
      integer nbigm
      
      type(scalar_field), pointer:: geoeli1p_field
      character(len=OPTION_PATH_LEN) velocity_option_path, &
     &     geoeli1p_option_path,pressure_option_path,&
     &     tmp_option_path
      integer :: istate

      type(mesh_type), pointer :: velocity_mesh, pressure_mesh, coordinate_mesh
      
      logical :: have_gravity
      real :: gravity_mag
      type(scalar_field), pointer :: buoyancy
      type(vector_field), pointer :: velocity, nonlinear_velocity, grid_velocity
      type(vector_field), pointer :: velocity_source, gravity, subgrid_source
      type(vector_field), pointer :: coordinate, old_coordinate, original_coordinate
      type(vector_field), pointer :: absorption, gravity_direction
      type(vector_field), pointer :: subgrid_velocity, nonlinear_subgrid_velocity
      type(vector_field) :: ns_rhs

      type(vector_field) :: masslump

      type(scalar_field), pointer :: pressure, density
      type(scalar_field), pointer :: equilibrium_tide
      type(scalar_field) :: compressible_projec_rhs
      type(scalar_field), pointer :: shear_stress_x, shear_stress_y, shear_stress_z
      type(scalar_field), pointer :: mean_shear_stress_x, mean_shear_stress_y, mean_shear_stress_z

      ! for LES -Jemma
      type(tensor_field) :: metric_tensor
      REAL, dimension(:), allocatable:: ML2MXX,ML2MXY,ML2MXZ
      REAL, dimension(:), allocatable:: ML2MYY,ML2MYZ,ML2MZZ
      integer :: NDIM

      ! Viscosity
      type(tensor_field) :: viscosity
      REAL, dimension(:), allocatable :: MUPTXX,MUPTXY,MUPTXZ
      REAL, dimension(:), allocatable :: MUPTYY,MUPTYZ,MUPTZZ
      
      ! shape functions, quadrature weights
      real, dimension(:), pointer :: weight, sweigh
      real, dimension(:, :), allocatable :: n, nlx, nly, nlz
      real, dimension(:, :), allocatable :: sn, snlx, snly
      real, dimension(:, :), allocatable :: m, mlx, mly, mlz

      logical :: get_grad_p_cg=.true., get_grad_p_cv = .false.

      integer:: gottraf
      type(scalar_field),pointer:: SolidPhase 
      type(vector_field),pointer:: SolidVelocity
      real,dimension(:),pointer:: utraf,vtraf,wtraf,voltraf
      real, dimension(:,:), allocatable, target:: tmp_gravity_direction

      ewrite(1, *)  "SUBROUTINE NAVSTO()"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reconstructing various variables from state and option tree
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call get_option("/timestepping/current_time", acctim)
     
      if (IsParallel()) then
        para=1
      else
        para=0
      end if
      
      call get_option("/geometry/dimension", ndim)
      select case (ndim)
      case (1)
         D3=.false.
      case (2)
        D3=.false.
      case (3)
        D3=.true.
      case default
         FLAbort("Unknown dimension")
      end select
      
      ! state with the prognostic velocity
      istate=phase2state_index(1)
      
      velocity  => extract_vector_field(state(istate), "Velocity")
! option path to velocity field for solver options
      velocity_option_path=velocity%option_path
      if (.not.have_option(trim(velocity_option_path)//"/prognostic")) then
         return
      end if
      u => velocity%val(1)%ptr
      v => velocity%val(2)%ptr
      if (D3) w => velocity%val(3)%ptr

      nonlinear_velocity  => extract_vector_field(state(istate), "NonlinearVelocity")
      nu => nonlinear_velocity%val(1)%ptr
      nv => nonlinear_velocity%val(2)%ptr
      if (D3) nw => nonlinear_velocity%val(3)%ptr
      
      allocate(zer(nonods))
      ZER = 0.0
      
      grid_velocity => extract_vector_field(state(istate), "GridVelocity", stat)
      if(stat == 0) then
        ug => grid_velocity%val(1)%ptr
        vg => grid_velocity%val(2)%ptr
        if (D3) wg => grid_velocity%val(3)%ptr
      else
        ug => ZER
        vg => ZER
        wg => ZER
      end if
      
      velocity_source => extract_vector_field(state(istate), "VelocitySource", stat)
      if(stat == 0) then
        sourcx => velocity_source%val(1)%ptr
        sourcy => velocity_source%val(2)%ptr
        if (D3) sourcz => velocity_source%val(3)%ptr
      else
        sourcx => ZER
        sourcy => ZER
        sourcz => ZER
      end if
      
      call get_option("/physical_parameters/gravity/magnitude", gravity_mag, stat=stat)
      have_gravity = stat==0
      if(have_gravity) then
         gravity => extract_vector_field(state(istate), "GravityDirection")
         buoyancy => extract_scalar_field(state(istate), "VelocityBuoyancyDensity")
         allocate(soxgra(nonods))
         allocate(soygra(nonods))
         allocate(sozgra(nonods))
         select case(gravity%field_type)
         case(FIELD_TYPE_CONSTANT)
            soxgra = gravity%val(1)%ptr(1)*gravity_mag*buoyancy%val
            soygra = gravity%val(2)%ptr(1)*gravity_mag*buoyancy%val
            if (ndim==3) sozgra = gravity%val(3)%ptr(1)*gravity_mag*buoyancy%val
         case default
            soxgra = gravity%val(1)%ptr*gravity_mag*buoyancy%val
            soygra = gravity%val(2)%ptr*gravity_mag*buoyancy%val
            if (ndim==3) sozgra = gravity%val(3)%ptr*gravity_mag*buoyancy%val
         end select
      else
         soxgra => ZER
         soygra => ZER
         sozgra => ZER
      end if
      
      nsograsub=0
      soxgrasub => ZER
      soygrasub => ZER
      sozgrasub => ZER
      
      if (have_option("/traffic_model")) then
         SolidPhase => extract_scalar_field(state(istate),"SolidPhase")
         voltraf    => SolidPhase%val
         
         SolidVelocity => extract_vector_field(state(istate), "SolidVelocity")
         utraf   => SolidVelocity%val(1)%ptr
         vtraf   => SolidVelocity%val(2)%ptr
         wtraf   => SolidVelocity%val(3)%ptr
      else
         voltraf => zer 
         utraf   => zer 
         vtraf   => zer 
         wtraf   => zer 
      endif
      
      coordinate => extract_vector_field(state(istate), "Coordinate")
      x => coordinate%val(1)%ptr
      y => coordinate%val(2)%ptr
      if (d3) then
        z => coordinate%val(3)%ptr
      else
        z => ZER
      end if
      
      old_coordinate => extract_vector_field(state(istate), "OldCoordinate", stat=stat)
      if (stat==0) then
        xold => old_coordinate%val(1)%ptr
        yold => old_coordinate%val(2)%ptr
        if (d3) then
          zold => old_coordinate%val(3)%ptr
        else
          zold => ZER
        end if
      else
        xold => x
        yold => y
        zold => z
      end if

      original_coordinate => extract_vector_field(state(istate), "OriginalCoordinate", stat=stat)
      if (stat==0) then
        xorig => original_coordinate%val(1)%ptr
        yorig => original_coordinate%val(2)%ptr
        if (d3) then
          zorig => original_coordinate%val(3)%ptr
        else
          zorig => ZER
        end if
      else
        xorig => x
        yorig => y
        zorig => z
      end if
      
      ! then allocate a rhs for the navier stokes equation
      call allocate(ns_rhs, velocity%dim, velocity%mesh, &
        name="NavierStokesRHS")
      call zero(ns_rhs)      
      vecx => ns_rhs%val(1)%ptr
      vecy => ns_rhs%val(2)%ptr
      if (d3) then
        vecz => ns_rhs%val(3)%ptr
      else
        vecz => ZER
      end if

      ! Viscosity tensor
      viscosity=extract_tensor_field(state(istate), "Viscosity", stat)
      allocate(MUPTXX(NONODS), MUPTXY(NONODS), MUPTXZ(NONODS), &
           MUPTYY(NONODS), MUPTYZ(NONODS), MUPTZZ(NONODS))
      if (stat==0) then
         call copy_tensor_to_legacy(viscosity, MUPTXX, MUPTXY, MUPTXZ, &
              MUPTYY, MUPTYZ, MUPTZZ)
      else
         MUPTXX=0.0
         MUPTXY=0.0
         MUPTXZ=0.0
         MUPTYY=0.0
         MUPTYZ=0.0
         MUPTZZ=0.0
      end if
      
      ! reconstruct meshes
      velocity_mesh => extract_mesh(state(istate), "VelocityMesh")
      pressure_mesh => extract_mesh(state(istate), "PressureMesh")
      coordinate_mesh => extract_mesh(state(istate), "CoordinateMesh")
      
      ndglno => velocity_mesh%ndglno
      pndgln => pressure_mesh%ndglno
      xondgl => coordinate_mesh%ndglno
      sondgl => ndglno
      vondgl => ndglno
      
      allocate(sndgln(1:stotel*snloc))
      call getsndgln(coordinate_mesh, sndgln)
      
      ! integration weights
      weight => velocity_mesh%shape%quadrature%weight
      sweigh => velocity_mesh%faces%shape%quadrature%weight
      
      ! reconstruct old style shape functions
      allocate(n(ele_loc(velocity_mesh, 1), ele_ngi(velocity_mesh, 1)))
      allocate(nlx(ele_loc(velocity_mesh, 1), ele_ngi(velocity_mesh, 1)))
      allocate(nly(ele_loc(velocity_mesh, 1), ele_ngi(velocity_mesh, 1)))
      allocate(nlz(ele_loc(velocity_mesh, 1), ele_ngi(velocity_mesh, 1)))
      call extract_old_element(ele_shape(velocity_mesh, 1), n, nlx, nly, nlz)
      
      allocate(sn(ele_loc(velocity_mesh, 1), face_ngi(velocity_mesh, 1)))
      allocate(snlx(ele_loc(velocity_mesh, 1), face_ngi(velocity_mesh, 1)))
      allocate(snly(ele_loc(velocity_mesh, 1), face_ngi(velocity_mesh, 1)))
      call extract_old_element(face_shape(velocity_mesh, 1), sn, snlx, snly)
      
      allocate(m(ele_loc(pressure_mesh, 1), ele_ngi(pressure_mesh, 1)))
      allocate(mlx(ele_loc(pressure_mesh, 1), ele_ngi(pressure_mesh, 1)))
      allocate(mly(ele_loc(pressure_mesh, 1), ele_ngi(pressure_mesh, 1)))
      allocate(mlz(ele_loc(pressure_mesh, 1), ele_ngi(pressure_mesh, 1)))
      call extract_old_element(ele_shape(pressure_mesh, 1), m, mlx, mly, mlz)
      
      ! the matrices
      ct_m => get_velocity_divergence_matrix(state)
      c1t => ct_m%val(1,1)%ptr
      c2t => ct_m%val(1,2)%ptr
      if (d3) then
        c3t => ct_m%val(1,3)%ptr
      else
        c3t => c1t
      end if
      findct => ct_m%sparsity%findrm
      colct => ct_m%sparsity%colm
      
      nct=size(colct)

      call allocate(ctp_m, ct_m%sparsity, blocks=ct_m%blocks, name="CompressiblePressureGradientMatrix")
      c1tp => ctp_m%val(1,1)%ptr
      c2tp => ctp_m%val(1,2)%ptr
      if (d3) then
        c3tp => ctp_m%val(1,3)%ptr
      else
        c3tp => c1tp
      end if
      
      u_sparsity => get_csr_sparsity_firstorder(state, velocity_mesh, velocity_mesh)
      findrm => u_sparsity%findrm
      colm => u_sparsity%colm
      centrm => u_sparsity%centrm
      
      ncolm=size(colm)

      nbigm=size(u_sparsity%colm)*ndim**2
      nbig12=size(u_sparsity%colm)*ndim**2

      allocate( bigm1(1:nbigm) )
      big_m=wrap(u_sparsity,blocks=(/ ndim, ndim /), val=bigm1, &
         name="CoupledMomentumMatrix")
      call zero(big_m)
      
      cmc_m => get_pressure_poisson_matrix(state)
      cmc => cmc_m%val
      fincmc => cmc_m%sparsity%findrm
      colcmc => cmc_m%sparsity%colm
      midcmc => cmc_m%sparsity%centrm
      
      ncmc=size(colcmc)

      kmk_m => get_pressure_stabilisation_matrix(state)
      kcmc => kmk_m%val
      
      call allocate(masslump, ndim, velocity_mesh, name="MassLumpedMomentumMatrix")
      ! note we're pointing to the first component only,
      ! so this doesn't work for continous_galerkin_test and "pressure_corrected" absorption
      ml => masslump%val(1)%ptr
      
      density => extract_scalar_field(state(istate), "Density")
      denpt => density%val
      pressure => extract_scalar_field(state(istate), "Pressure")
      p => pressure%val
      
      absorption => extract_vector_field(state(istate), "VelocityAbsorption", stat=stat)
      !     Set XABSOR so it is not overwriten in this sub.
      ALLOCATE(XABSOR(SNONOD))
      ALLOCATE(YABSOR(SNONOD))
      ALLOCATE(ZABSOR(SNONOD))
      if (stat==0) then
         XABSOR=absorption%val(1)%ptr
         YABSOR=absorption%val(2)%ptr
         if (D3) then
            ZABSOR=absorption%val(3)%ptr
         else
            ZABSOR=0.0
         end if
      else
         XABSOR=0.0
         YABSOR=0.0
         ZABSOR=0.0
      end if
      
      allocate(dm1(1:nonods), dm2(1:nonods), dm3(1:nonods), dp(1:fredop))
      
      ! rhs for the continuity/pressure projection equation
      allocate(rhs(1:fredop))
            
      ! then allocate a rhs for the compressible projection equation
      call allocate(compressible_projec_rhs, pressure_mesh, name="CompressibleProjectionRHS")
      divqs => compressible_projec_rhs%val
      call zero(compressible_projec_rhs)  ! in case it's not used
      
      gravity_direction => extract_vector_field(state(istate), &
         name="GravityDirection", stat=stat)
      if (stat==0) then
         if (gravity_direction%field_type==FIELD_TYPE_NORMAL) then
            vgravx => gravity_direction%val(1)%ptr
            vgravy => gravity_direction%val(2)%ptr
            if (d3) then
               vgravz => gravity_direction%val(3)%ptr
            else
               vgravz => zer
            end if
         elseif (gravity_direction%field_type==FIELD_TYPE_CONSTANT) then
            allocate( tmp_gravity_direction(1:nonods, 1:ndim) )
            do i=1, ndim
              tmp_gravity_direction(:,i)=gravity_direction%val(i)%ptr(1)
            end do
            vgravx => tmp_gravity_direction(:,1)
            vgravy => tmp_gravity_direction(:,2)
            vgravz => tmp_gravity_direction(:,3)
         else
            FLAbort("Unknown field_type for GravityDirection")
         end if
      else
         vgravx => zer
         vgravy => zer
         vgravz => zer
      endif
            
      equilibrium_tide => extract_scalar_field(state(istate), &
         name="EquilibriumTide", stat=stat)
      if (stat==0) then
        eqtd => equilibrium_tide%val
        got_eqtd=.true.
      else
        eqtd => zer
        got_eqtd=.false.
      end if
      
      shear_stress_x => extract_scalar_field(state(istate), &
         name="MaxBedShearStressX", stat=stat)
      if (stat==0) then
        shearstressx => shear_stress_x%val
        got_shearstressx=.true.
      else
        shearstressx => zer
        got_shearstressx=.false.
      end if

      shear_stress_y => extract_scalar_field(state(istate), &
         name="MaxBedShearStressy", stat=stat)
      if (stat==0) then
        shearstressy => shear_stress_y%val
        got_shearstressy=.true.
      else
        shearstressy => zer
        got_shearstressy=.false.
      end if
      
      shear_stress_z => extract_scalar_field(state(istate), &
         name="MaxBedShearStressz", stat=stat)
      if (stat==0) then
        shearstressz => shear_stress_z%val
        got_shearstressz=.true.
      else
        shearstressz => zer
        got_shearstressz=.false.
      end if

      mean_shear_stress_x => extract_scalar_field(state(istate), &
         name="MaxBedmeanshearstressX", stat=stat)
      if (stat==0) then
        meanshearstressx => mean_shear_stress_x%val
        got_meanshearstressx=.true.
      else
        meanshearstressx => zer
        got_meanshearstressx=.false.
      end if

      mean_shear_stress_y => extract_scalar_field(state(istate), &
         name="MaxBedmeanshearstressy", stat=stat)
      if (stat==0) then
        meanshearstressy => mean_shear_stress_y%val
        got_meanshearstressy=.true.
      else
        meanshearstressy => zer
        got_meanshearstressy=.false.
      end if
      
      mean_shear_stress_z => extract_scalar_field(state(istate), &
         name="MaxBedmeanshearstressz", stat=stat)
      if (stat==0) then
        meanshearstressz => mean_shear_stress_z%val
        got_meanshearstressz=.true.
      else
        meanshearstressz => zer
        got_meanshearstressz=.false.
      end if
      
      subgrid_velocity => extract_vector_field(state(istate), &
         name="VelocityInnerElement", stat=stat)
      if (stat==0) then
         usub => subgrid_velocity%val(1)%ptr
         vsub => subgrid_velocity%val(2)%ptr
         if (d3) then
            wsub => subgrid_velocity%val(3)%ptr
         else
            wsub => zer
         endif
      else
         usub => zer
         vsub => zer
         wsub => zer
      end if

      nonlinear_subgrid_velocity => extract_vector_field(state(istate), &
         name="NonlinearVelocityInnerElement", stat=stat)
      if (stat==0) then
         nusub => nonlinear_subgrid_velocity%val(1)%ptr
         nvsub => nonlinear_subgrid_velocity%val(2)%ptr
         if (d3) then
            nwsub => nonlinear_subgrid_velocity%val(3)%ptr
         else
            nwsub => zer
         endif
      else
         nusub => zer
         nvsub => zer
         nwsub => zer
      end if

      allocate(nodrot(0))
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ewrite_minmax(denpt)
#if defined(HAVE_NETCDF)&&defined(HAVE_LIBUDUNITS)&&defined(DOUBLEP)
      if(have_option("/environmental_data/ERA40")) then
         call fluxes_settimeseconds(acctim)
      end if
#endif
      

      
      ALLOCATE(VECXMOM(NONODS)) ; VECXMOM = 0.0
      ALLOCATE(VECYMOM(NONODS)) ; VECYMOM = 0.0
      ALLOCATE(VECZMOM(NONODS)) ; VECZMOM = 0.0
      
      ALLOCATE(UOLD2(NONODS)) ; UOLD2 = 0.0
      ALLOCATE(VOLD2(NONODS)) ; VOLD2 = 0.0
      ALLOCATE(WOLD2(NONODS)) ; WOLD2 = 0.0
      
!     work out which field is solved in geoeli1p and its option_path
      if (geobal<=-10) then
         if (has_scalar_field(state(istate), "BalancePressure")) then
            geoeli1p_field => extract_scalar_field(state(istate), &
     &           "BalancePressure")
            geoeli1p_option_path=geoeli1p_field%option_path
         else if (has_scalar_field(state(istate), "FreeSurface")) then
            geoeli1p_field => extract_scalar_field(state(istate),&
     &           "FreeSurface")
            geoeli1p_option_path=geoeli1p_field%option_path
         else if (my_new_options) then
!     with new options the field must be declared otherwise
!     we don't have options for them
            ewrite(0,*) "geobal=", geobal
            ewrite(0,*) "gotpsi=", gotpsi
            ewrite(0,*) "With new options, for geoeli1p: "
            FLAbort("You need either a BalancePressure or FreeSurface field")
         end if
         psipre => geoeli1p_field%val
         gotpsi=.true.
      else
         psipre => u
         gotpsi=.false.
      end if
        
      NOFILT = 0
      
      NBUOY=0
      pressure_option_path='/material_phase['//int2str(istate-1)//&
     &     ']/scalar_field::Pressure'
      tmp_option_path="/prognostic/spatial_discretisation"// &
     &     "/continuous_galerkin/remove_stabilisation_term"
      if(have_option(trim(pressure_option_path)//trim(tmp_option_path))) then
         NOFILT=1
      end if      
      ewrite(2,*) 'pressure_option_path NOFILT = ',NOFILT
      
      ISPHERE=ISPHER2
      IF((.NOT.COGRAX).AND.(GEOBAL.LE.-10)) THEN
!     ISPHERE controls the order of supermarametric mapping 
!     for the spherical Earth
         ISPHERE=2
         ISPHER2=2
      ENDIF
      
      IF(DISOPT.EQ.999) RETURN
            
      ! allocate a local copy of the sources
      allocate(sx(nonods), sy(nonods), sz(nonods))
      ! add in the prescribed sources
      SX = SOURCX
      SY = SOURCY
      SZ = SOURCZ

      ! if the hydrostatic solver and free surface aren't being used then
      ! add the buoyancy terms into the local source vector sx, sy, sz
      ! this will then contain both buoyancy and prescribed sources
      if (my_new_options.and.(GEOBAL>-10)) then
         SX = SX + SOXGRA
         SY = SY + SOYGRA
         SZ = SZ + SOZGRA
      end if

      ! create a local copy of the buoyancy sources for the hydrostatic pressure
      ! and free surface solvers to use
      allocate(l_soxgra(nonods), l_soygra(nonods), l_sozgra(nonods))
      l_soxgra = soxgra
      l_soygra = soygra
      l_sozgra = sozgra

      
!**********************************************
!     SET UP MEM FOR DISCRETISED SGS SOURCE
      ALLOCATE(VECX_SGSADD(TOTELE*NSUBVLOC))
      ALLOCATE(VECY_SGSADD(TOTELE*NSUBVLOC))
      ALLOCATE(VECZ_SGSADD(TOTELE*NSUBVLOC))
      VECX_SGSADD=0.0
      VECY_SGSADD=0.0
      VECZ_SGSADD=0.0

      IF(GEOBAL.LT.0) THEN
!     Solve an elliptic eqn for geostrophic balance and 
!     put the resulting pressure in VECX,VECY,VECZ.

         IF(GEOBAL.LE.-10) THEN
!     If GEOBAL=-11 or -21 then include bouyancy in forcing terms...
!     If GOEBAL=-1? then have a linear variation in goe-hyd pressure
!     If GEOBAL=-2? then have quadratic variation.
!     If GEOBAL.LT.-30 then use a full blown option list
!     solve for quadratic pressure and put into momentum equations...
            ALLOCATE(PSIPRETEM(NONODS))
            IF(GOTPSI) THEN
               PSIPRETEM(1:NONODS)=PSIPRE(1:NONODS)
            ELSE
               PSIPRETEM(1:NONODS)=0.0
            ENDIF
            VECX=0.0
            VECY=0.0
            VECZ=0.0
            
            CALL GEOELI1P(VECX,VECY,VECZ,&
     &           NSUBNVLOC,VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
     &           SX,SY,SZ, &
     &           NU,NV,NW, &
     &           1,NSUBVLOC,NUSUB,NVSUB,NWSUB,&
     &           l_SOXGRA,l_SOYGRA,l_SOZGRA,PSIPRETEM,P,&
     &           NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
     &           NOBCFSimp,BCT1WFSimp,BCT2FSimp,&
     &           X,Y,Z, XOLD,YOLD,ZOLD,XORIG,YORIG,ZORIG, &
     &           UG,VG,WG,MVMESH, &
     &           N,NLX,NLY,NLZ, WEIGHT, NGI,&
     &           NLOC,SNLOC,&
     &           TOTELE,D3,DCYL, &
     &           GEOBAL,&
     &           NDGLNO,XONDGL,&
     &           NONODS,XNONOD,&
     &           FINDRM,CENTRM,COLM,NCOLM,&
     &           FREDOP,&
     &           NLEVEL,geostrophic_solver_option,FREESDISOTT,&
     &           DISOPT,&
     &           halo_tag_p, NNODP,&
     &           scfacth0,NEWMES2,&
     &           GOTTOPDIS, XABSOR,YABSOR,ZABSOR,&
     &           GOTBOY,BSOUX,BSOUY,BSOUZ, &
     &           COGRAX, &
     &           VGRAVX,VGRAVY,VGRAVZ,&
     &           (.NOT.IGUESS),GRAVTY,DT,ACCTIM,ISPHERE,EQTD,GOT_EQTD,&
     &           NOFILT,&
     &           C1T,C2T,C3T,FINDCT,COLCT,NCT, &
!     state contains all fields of this phase...
     &           state(istate),3,((.NOT.IGUESS).OR.(MISNIT.EQ.0)),&
     &           geoeli1p_option_path)
            
            ! only use NBUOY=1 if we go into DIFF3DSIMP            
            IF((ISPHERE.EQ.2).OR.((DISOPT.ge.144).AND.(DISOPT.LE.154))) NBUOY=1
! switch off implicit buoyancy on the sphere until we sort out bug...
            IF(ISPHERE.EQ.2) NBUOY=0
            
            if((geobal<0).and.(geostrophic_solver_option>=2)) then
               VECXMOM=VECX
               VECYMOM=VECY
               VECZMOM=VECZ
            ENDIF
            
!     We have just calculated C1T     
            IF(geostrophic_solver_option>=2) GETC12=.FALSE.
            
            SX(1:NONODS)=SX(1:NONODS)+SOURCX(1:NONODS)
            SY(1:NONODS)=SY(1:NONODS)+SOURCY(1:NONODS)
            SZ(1:NONODS)=SZ(1:NONODS)+SOURCZ(1:NONODS)
                        
            UOLD2(1:NONODS)=U(1:NONODS)
            VOLD2(1:NONODS)=V(1:NONODS)
            WOLD2(1:NONODS)=W(1:NONODS)

            IF(GOTPSI) THEN
               PSIPRE(1:NONODS)=PSIPRETEM(1:NONODS)
            ENDIF
            DEALLOCATE(PSIPRETEM)
            
            if ((itsiti==itinoi).and.(my_new_options)) then
               call free_surface_diagnostics(state(istate))
            endif
            
            
!     
         ELSE IF(GEOBAL.LE.-4) THEN

            FLAbort("With new options -10<GEOBAL<=-4 won't work. -Stephan")
            
         ENDIF
      ENDIF

      call quadratic_drag(state(istate), xabsor, yabsor, zabsor)

      IF(NSUBVLOC.NE.0) THEN
         ISUBSOU=0
         BIG_NLOC=3*NLOC
         SUF_TOP_OCEAN=(GEOBAL.LE.-10)
         INUSUB=1
         ISUB  =1
         IF(NSUBVLOC.EQ.0) STOP 6754
         
         ALLOCATE(UNEW(NONODS))
         ALLOCATE(VNEW(NONODS))
         ALLOCATE(WNEW(NONODS))
         ALLOCATE(UNEWSUB(TOTELE*NSUBVLOC))
         ALLOCATE(VNEWSUB(TOTELE*NSUBVLOC))
         ALLOCATE(WNEWSUB(TOTELE*NSUBVLOC))
         UNEW=U+(NU-U)/MAX(THETA,1.E-20)
         VNEW=V+(NV-V)/MAX(THETA,1.E-20)
         WNEW=W+(NW-W)/MAX(THETA,1.E-20)
         UNEWSUB=USUB+(NUSUB-USUB)/MAX(THETA,1.E-20)
         VNEWSUB=VSUB+(NVSUB-VSUB)/MAX(THETA,1.E-20)
         WNEWSUB=WSUB+(NWSUB-WSUB)/MAX(THETA,1.E-20)
         
         ALLOCATE(TXABSOR(NONODS))

         if (has_boundary_condition(velocity,"wind_forcing")) then
            call wind_forcing(state(istate), ns_rhs)
         end if
         
         if(have_option("/environmental_data/ERA40")) then
!     calculate ML
            ALLOCATE(tempML(NONODS))
            tempML=0.0
            VECX=0.0
            VECY=0.0
            VECZ=0.0
            CALL surface_stress(U,V,W,NU,NV,NW,X,Y,Z,&
     &           snloc, sngi, sn, snlx, snly, sweigh,&
     &           theta,dt,centrm,tempML,0,&
     &           BIGM1, VECX, VECY, VECZ)
!     Calculate ML...
            ML=0.0
      do ELE=1,TOTELE
               NOD1=NDGLNO((ELE-1)*NLOC+1)
               NOD2=NDGLNO((ELE-1)*NLOC+2)
               NOD3=NDGLNO((ELE-1)*NLOC+3)
               NOD4=NDGLNO((ELE-1)*NLOC+4)
               VOLELE=ABS(TetVolume( &
     &              x(NOD1), y(NOD1), z(NOD1),         &
     &              x(NOD2), y(NOD2), z(NOD2),         &
     &              x(NOD3), y(NOD3), z(NOD3),         &
     &              x(NOD4), y(NOD4), z(NOD4)))
               ML(NOD1)=ML(NOD1)+0.25*VOLELE
               ML(NOD2)=ML(NOD2)+0.25*VOLELE
               ML(NOD3)=ML(NOD3)+0.25*VOLELE
               ML(NOD4)=ML(NOD4)+0.25*VOLELE
            END DO
            
            TXABSOR=XABSOR+tempML/(ML*DT)
            SX=SX+VECX/ML
            SY=SY+VECY/ML
            SZ=SZ+VECZ/ML
         ELSE
            TXABSOR=XABSOR
         end if
         
         NEXTTIM=(.NOT.IGUESS)
         IUNST_FREE=0
         IF((GEOBAL.EQ.-26).OR.(GEOBAL.EQ.-27)) IUNST_FREE=1
         
         ! Chris: I've hard-coded ISUBSOU to 0 below as the interface
         ! is broken, we should not *ever* pass in arrays 
         ! that are not conforming in shape and size
         CALL SGSPROJ(U,V,W, UNEW,VNEW,WNEW,&
     &        NU,NV,NW,UG,VG,WG,&
     &        SX,SY,SZ,X,Y,Z,D0,&
     &        real_length_zero_array,real_length_zero_array, real_length_zero_array, 0,&
     &        VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA,&
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
     &        BIGM1, ML,&
     &        FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &        FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &        M,MLX,MLY,MLZ,&
     &        N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC,BIG_NLOC,&
     &        DISOPT,DISOPN,DT,THETA,BETA,LUMP,PREOPT,&
     &        ABSLUM,SLUMP,&
     &        NDGLNO,PNDGLN,&
     &        SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &        DENPT,&
     &        MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &        TXABSOR,VERSIO,ISPHERE,&
     &        NNODP,&
     &        GEOBAL,SCFACTH0,&
!     for the SGS...
     &        NLOC,&
     &        NUSUB,NVSUB,NWSUB,INUSUB,&
     &        USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &        NBIG12,COGRAX,NOBCFSimp,BCT2FSimp,P, &
!     &     NBIG12,COGRAX,P, 
     &        SUF_TOP_OCEAN,&
     &        NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,BCU1val,BCV1val,BCW1val,&
!     GMRES solver...
     &        PARA,halo_tag,halo_tag_p,&
!     form DGCMC
     &        NCMC,FINCMC,COLCMC,MIDCMC,&
     &        NDPSET,D3,NEWMES2, velocity_option_path, &
     &        NEXTTIM,THETA,GRAVTY,IUNST_FREE) 
     

            IF (GOT_SHEARSTRESSX .or. GOT_MEANSHEARSTRESSX) THEN
               ewrite(1,*) 'Calling calculate_shear_stress'
               CALL calculate_shear_stress(ACCTIM,U,V,W,SCFACTH0,&
     &              VECX,VECY,VECZ, NONODS,STOTEL,SNLOC,SNGI,&
     &              D3,DCYL,DSPH, &
     &              SNDGLN,  &
     &              SN,SNLX,SNLY,SWEIGH, &
     &              X,Y,Z,R0,D0, &
     &              BIGM1,NBIGM,THETA,&
     &              FINDRM,COLM,NCOLM,CENTRM,ITSITI,ITINOI,&
!     make implicit in pressure as well
     &              DT,ML,SHEARSTRESSX,GOT_SHEARSTRESSX,SHEARSTRESSY,&
     &              GOT_SHEARSTRESSY,SHEARSTRESSZ,GOT_SHEARSTRESSZ,&
     &              MEANSHEARSTRESSX,GOT_MEANSHEARSTRESSX,MEANSHEARSTRESSY,&
     &              GOT_MEANSHEARSTRESSY,MEANSHEARSTRESSZ,GOT_MEANSHEARSTRESSZ)
            ENDIF

         ! Added to solve memory leak issues with the Inner Element method
            call deallocate(ns_rhs)
            call deallocate(masslump)
            call deallocate(compressible_projec_rhs)
            call deallocate(big_m)

         RETURN
      ENDIF
      
      ALLOCATE(BOY_ML(NBUOY*9*NONODS))

      ! IF LES THEN PUT INVERSE OF METRIC INTO MUPTXX,MUPTXY ETC
      !
      ! If DISOPT=42, constant length scale
      ! If DISOPT=43, isotropic length scale
      ! If DISOPT=44, inverted metric used
      ! If DISOPT=45 or 46, modified inverted metric using rotation tensors
  
      allocate(ml2mxx(snonod), ml2mxy(snonod), ml2mxz(snonod))
      allocate(ml2myy(snonod), ml2myz(snonod), ml2mzz(snonod))
      ml2mxx=0.0
      ml2mxy=0.0
      ml2mxz=0.0
      ml2myy=0.0
      ml2myz=0.0
      ml2mzz=0.0

      if(disopt>=42 .and. disopt<=50) then

        call eletens(x, y, z, &
          ml2mxx, ml2mxy, ml2mxz, ml2myy, ml2myz, ml2mzz,&
          fredop, snonod, xnonod, totele, nloc, disopt, ndglno)

        IF(DISOPT.EQ.42) THEN
            CALL LESIMEC(SNONOD,NDIM,&
              & ML2MXX,ML2MXY,ML2MXZ,&
              & ML2MYY,ML2MYZ,ML2MZZ)
        ENDIF

        IF(DISOPT.EQ.43) THEN
            CALL LESIME0(SNONOD,NDIM, metric_tensor%val,&
              & ML2MXX,ML2MXY,ML2MXZ,&
              & ML2MYY,ML2MYZ,ML2MZZ)
        ENDIF

        IF(DISOPT.EQ.44) THEN
            CALL LESIME(SNONOD,NDIM, metric_tensor%val,&
              &  ML2MXX,ML2MXY,ML2MXZ,&
              &  ML2MYY,ML2MYZ,ML2MZZ)
        ENDIF

        IF((DISOPT.EQ.45).OR.(DISOPT.EQ.46)) THEN
            ewrite(3,*) 'MUPTXX,MUPTXY,MUPTXZ:',&
                & MUPTXX,MUPTXY,MUPTXZ
            ewrite(3,*) 'MUPTYY,MUPTYZ,MUPTZZ:',&
                & MUPTYY,MUPTYZ,MUPTZZ
            ewrite(3,*) 'JUST BEFORE LESIME2'

            CALL LESIME2(SNONOD,NDIM, metric_tensor%val,&
              & ML2MXX,ML2MXY,ML2MXZ,&
              & ML2MYY,ML2MYZ,ML2MZZ)

        ENDIF
      end if
  
        ewrite(1, *) 'about to enter assnav NPHASE',NPHASE
  
        CALL ASSNAV&
      &       (acctim,U,V,W,&
      !     --------------------------------------START OF ADD BY CRGW 31/03/06
      !          PLASTIC STRAINS:
      &      STVPXX,STVPYY,STVPZZ,&
      &      STVPYZ,STVPXZ,STVPXY,&
      !     --------------------------------------END OF ADD BY CRGW 31/03/06
      &      NU,NV,NW,UG,VG,WG,&
      &      SX,SY,SZ,X,Y,Z,R0,D0,&
      &      NBUOY,BOY_ML,SOXGRA,SOYGRA,SOZGRA,&
      &      NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
      &      VECX,VECY,VECZ, &
      &      C1T,C2T, BIGM1, NBIG12, &
      &      NCMC,&
      &      SOURCZ,ML,&
      &      FINDCT,COLCT,NCT,FREDOP,TOTELE,&
      &      FINDRM,COLM,NCOLM,NONODS,CENTRM, &
      &      M,&
      &      N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, &
      &      DISOPT,DISOPN,DT,THETA,BETA,&
      &      GEOBAL,&
      &      LUMP,&
      &      MAKSYM,&
      &      GETC12,&
      &      CONVIS, &
      &      ABSLUM,SLUMP,&
      &      NDGLNO,PNDGLN, &
      &      SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
      &      DENPT,&
      &      MUPTXX,MUPTXY,MUPTXZ, MUPTYY,MUPTYZ,MUPTZZ,&
      &      ML2MXX,ML2MXY,ML2MXZ, ML2MYY,ML2MYZ,ML2MZZ, &
      &      XABSOR,YABSOR,ZABSOR, ZER,&
      &      D3,DCYL,DSPH,&
      &      STOTEL,SNLOC,SNGI,SNDGLN,  &
      &      SN,SNLX,SNLY,SWEIGH, &
      &      SCFACTH0,&
      &      GOTTRAF, &
      &      VERSIO,ISPHERE,&
      &      nnodp,para,halo_tag,rotat,nodrot,nnodro,&
      &      state, istate, &
      &      big_m, CT_m, CTP_m, CMC_m,&
      !     --------------------------------------START OF ADD BY CRGW 14/03/06
      !           INTEGER SOLIDS FLAG:
      &      SOLIDS, MKCOMP, &
      &      pressure_option_path &
      !     --------------------------------------END OF ADD BY CRGW 14/03/06
      &      )

      IF(MVMESH.AND.SOLIDS.EQ.0) THEN

!     Amend the r.h.s vector for a matrix eqn with a moving mesh. 
         CALL AMENDMVMESH(3,VECX,VECY,VECZ,&
     &        U,V,W,DT,&
     &        X,Y,Z, XOLD,YOLD,ZOLD, &
     &        N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
     &        TOTELE,D3,DCYL,  &
     &        XONDGL,NDGLNO, NONODS,XNONOD)
      ENDIF
      
!     FOR TRAFFIC *******************************
      IF(GOTTRAF.EQ.1) THEN

         do NOD=1,NONODS
            ML(NOD)=ML(NOD)/MAX(1.E-10,1.-VOLTRAF(NOD))
         END DO
         ewrite(3,*) 'IN NAVSTO.f:After Line 548:AFTER TRAFFIC EQUATION'
         
         CALL PMINMX(VOLTRAF,NONODS,'******VOLTRAF  ')
         CALL PMINMX(UTRAF,NONODS,'******UTRAF  ')
         CALL PMINMX(VTRAF,NONODS,'******VTRAF  ')
         CALL PMINMX(WTRAF,NONODS,'******WRAF  ')
      ENDIF
      !     FOR TRAFFIC *******************************
      !     
      !     
      !     call to CVCTS (now ASSEMBLE_CTS_CV, which used to be via CVCTS2)
      !     has moved to SUBROUTINE CMC_WRAPPER IN ASSEMBLE_CMC.F90 (CRGW, 03/09/06)
      !     
      ewrite(3,*) "DISOPN,GEOBAL = ",DISOPN,GEOBAL
      !     Add in the 1st order spatial derivative -transport. ******
      IF(DISOPN.NE.0) THEN
         !     Use and extended halo in the high resolution method...
         IF(D3.AND.(NLOC.EQ.4)) THEN
            IF(PARA.EQ.1) THEN
               CALL HALGET(U,NONODS,NONODS,NNODP,halo_tag_p)
               CALL HALGET(V,NONODS,NONODS,NNODP,halo_tag_p)
               IF(D3) CALL HALGET(W,NONODS,NONODS,NNODP,halo_tag_p)
               CALL HALGET(NU,NONODS,NONODS,NNODP,halo_tag_p)
               CALL HALGET(NV,NONODS,NONODS,NNODP,halo_tag_p)
               IF(D3) CALL HALGET(NW,NONODS,NONODS,NNODP,halo_tag_p)
               CALL HALGET(UG,NONODS,NONODS,NNODP,halo_tag_p)
               CALL HALGET(VG,NONODS,NONODS,NNODP,halo_tag_p)
               IF(D3) CALL HALGET(WG,NONODS,NONODS,NNODP,halo_tag_p)
            ENDIF
         ENDIF
         
         CALL ASSEMBLE_MOMENTUM_ADVECTION_CV(&
     &        NONODS,VNONOD,XNONOD,TOTELE,NLOC,NGI,&
     &        SNLOC, STOTEL,&
     &        NDGLNO,XONDGL,VONDGL,&
     &        SNDGLN,&
     &        FINDRM,COLM,NCOLM,CENTRM,&
     &        NBIGM,BIGM1, VECX,VECY,VECZ,&
     &        X,Y,Z,&
     &        N,NLX,NLY,NLZ, WEIGHT,&
     &        NU,NV,NW, UG,VG,WG, &
     &        U,V,W, U,V,W,&
     &        DENPT,DENPT,&
     &        DISOPN,DT,THETA,BETA, D3,DCYL,&
     &        PARA,halo_tag,NNODP,&
     &        NOBCU,BCU1,BCU2,&
     &        NOBCV,BCV1,BCV2,&
     &        NOBCW,BCW1,BCW2, velocity_option_path)
      ENDIF
!     Add in the 1st order spatial derivative -transport. ******
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     hydrostatic pressure
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      if (has_boundary_condition(velocity,"wind_forcing")) then
         call wind_forcing(state(istate), ns_rhs)
      end if
      
      if(have_option("/environmental_data/ERA40")) then
         CALL surface_stress(U,V,W,NU,NV,NW,X,Y,Z,&
     &        snloc, sngi, sn, snlx, snly, sweigh,&
     &        theta,dt,centrm,ML,ncolm,&
     &        BIGM1, VECX, VECY, VECZ)
      end if
      
      
      IDIM=2
      IF(D3) IDIM=3

      ! allocation of arrays for rotated bcs (if not used nnodro is 0)
      ! ( rtdr is a rotation matrix used in solnav() )
      if(allocated(nodrot)) deallocate(nodrot)
      allocate( normx(1:nnodro), normy(1:nnodro), normz(1:nnodro), &
         t1x(1:nnodro), t1y(1:nnodro), t1z(1:nnodro), &
         t2x(1:nnodro), t2y(1:nnodro), t2z(1:nnodro), &
         nodrot(1:nnodro), rtdr(1:nrtdr))
         
      IF(ROTAT) THEN
         call getnodrot(velocity, nodrot, gettan, &
              normx, normy, normz, &
              t1x, t1y, t1z, &
              t2x, t2y, t2z)
         IF(GETTAN) THEN
            ewrite(3,*)'about to go into getrot NNODRO',NNODRO
            CALL GETNOR( NORMX,NORMY,NORMZ, &
     &           NODROT,NNODRO,D3, &
     &           NONODS,FREDOP,&
     &           C1T,C2T,C3T,FINDCT,COLCT,NCT,&
     &           PARA,halo_tag,&
     &           NNODP)

            IF((GEOBAL.LE.-4)) THEN
!     Rotate for ocean modelling. 
               
               CALL ROTOCE( GETTAN, NODROT,&
     &              NOBCU,BCU2,NOBCV,BCV2,NOBCW,BCW2,&
     &              NORMX,NORMY,NORMZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &              NNODRO,D3)
            ELSE
!     Rotate for other things.
               CALL GETROT( GETTAN, &
     &              NORMX,NORMY,NORMZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &              NNODRO,D3)

            ENDIF

            ! store in new structures
            call copy_back_tangents(velocity, &
              normx, normy, normz, &
              t1x, t1y, t1z, &
              t2x, t2y, t2z)
            ! don't recalculate next time
            GETTAN=.FALSE.
         ENDIF
         ewrite(3,*) 'GETC12,nloc:',GETC12,nloc

!     ROTAT VECX,VECY,VECZ before we apply the b.c's
!     Rotat  vec=R vec - to unrotated system. 
         CALL ROTUVW(VECX,VECY,VECZ,NONODS,&
     &        NORMX,NORMY,NORMZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z, &
     &        NODROT,NNODRO,D3,NTRANS)
      ENDIF
      
      ewrite(1,*)'going into boucon ROTAT,GETTAN',ROTAT,GETTAN
      CALL BOUCON(IDIM, NONODS, &
     &     DM1,NOBCU,BCU1,BCU2,VECX,&
     &     DM2,NOBCV,BCV1,BCV2,VECY,&
     &     DM3,NOBCW,BCW1,BCW2,VECZ    )
!     
      ROTMOM=ROTAT
      IF(ROTAT) THEN
         IF(NPHASE.LE.1) THEN
            ewrite(1,*)'going into RBIGRT'
            CALL RBIGRT(BIGM1,CENTRM,FINDRM,COLM,&
     &           NCOLM,NBIGM,NONODS,&
     &           NORMX,NORMY,NORMZ, &
     &           T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &           NODROT,NNODRO,D3)
         ENDIF
         ROTMOM=.FALSE.
!     If ROTMOM=.FALSE. then momentum equations have already been rotated. 
      ENDIF

      IF((.NOT.ROTAT).OR.(.NOT.ROTMOM)) THEN
         IF(NPHASE.LE.1) THEN
            ewrite(1,*)'going into bouVIS'
            CALL BOUVIS(IDIM, NONODS, DM1,DM2,DM3,&
     &           BIGM1,NBIGM,NCOLM,CENTRM)
         ENDIF
      ENDIF
!     Get the pressure matrix.  
!     NB If we are starting to use reduced quadrature on the mass matrix 
!     then we probably need to have CMCGET=.TRUE. once, to obain 
!     another CMC matrix.
!     
!     
      IF(ROTAT.AND.GETC12) THEN
!     rotat the CiT matrices before calling GETCMC. 
         ewrite(3,*)'about to go into ROCTRT'
         CALL ROCTRT(C1T,C2T,C3T,FINDCT,COLCT,&
     &        NCT,FREDOP,NONODS,&
     &        NORMX,NORMY,NORMZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &        NODROT,NNODRO,D3)
      ENDIF
      
      IF(NPHASE.GT.1) CMCGET=.FALSE.
!     
!     
!     calls to GETCMC amd CMCFILT have moved to ASSEMBLE_CMC.F90 (CRGW, 03/09/06)
!     
      ewrite(2,*) 'pressure_option_path NOFILT = ',NOFILT
      
      CALL CMC_WRAPPER(&
     &     NONODS, XNONOD, TOTELE&
     &     , GETC12, NCT, FINDCT, COLCT&
     &     , CMCGET, FREDOP, FINCMC, COLCMC, MIDCMC&
     &     , NDGLNO, PNDGLN, XONDGL, SONDGL&
     &     , NCMC, C1T, C2T, C3T, C1TP, C2TP, C3TP&
     &     , PREOPT, X, Y, Z&
     &     , D3, MKCOMP&
     &     , NBUOY,BOY_ML,NDWISP&
     &     , NLOC, NGI, N, NLX, NLY, NLZ&
     &     , MLOC, M, MLX, MLY, MLZ, WEIGHT&
     &     , DM1, DM2, DM3&
     &     , NDPSET, ML&
     &     , CMC, KCMC&
     &     , ROTAT, NNODRO, NODROT&
     &     , NORMX, NORMY, NORMZ, T1X, T1Y, T1Z, T2X, T2Y, T2Z&
     &     , ISPHERE&
!     consistency additions:
     &     , NOFILT&
     &     , state)
!     --------------------------------------start of add by crgw 20/03/06

      if(have_option(trim(pressure_option_path)//"/prognostic/scheme/use_compressible_projection_method")) THEN

        call assemble_compressible_projection(state, cmc_m, compressible_projec_rhs, dt, cmcget)

      end if
      
      !     solve the Navier Stokes equations
!++++++++++++++++++++++++++++++++++++++++++++++
!     
!     Convert variables for Navier Stokes solution. 
      ID3=0
      IF(D3) ID3=1
      CONMAS=0
      IF(.NOT.LUMP) CONMAS=1
      IF(MIXMAS.EQ.1) CONMAS=0
      D2HALF=0
      COUPLE=1
      ewrite(3,*)'INSIDE NAVSTO:'
      ewrite(3,*)'poison,projec,uzawa,mulpa,cgsolq,gmresq',&
     &     poison,projec,uzawa,mulpa,cgsolq,gmresq
      ewrite(3,*) 'iphase:',iphase
!     
      ewrite(1,*) 'entering solnav NCMC,NCMCB,NPRESS,NPROPT',&
     &     NCMC,NCMCB,NPRESS,NPROPT
         
      IF((GEOBAL.LT.0).AND.(geostrophic_solver_option>=2)) THEN
      !     Work out rhs of momentum not including geobal        
        VECXMOM=VECX-VECXMOM
        VECYMOM=VECY-VECYMOM
        VECZMOM=VECZ-VECZMOM
      ENDIF

      !     NOW SOLVE FOR A FREE SURFACE...
      if((geobal<0).and.(geostrophic_solver_option>=2)) then
           UOLD2(1:NONODS)=U(1:NONODS)
           VOLD2(1:NONODS)=V(1:NONODS)
           WOLD2(1:NONODS)=W(1:NONODS)
      ENDIF
        
      ewrite(1,*) 'entering solnav'
      CALL SOLNAV(U,V,W,NU,NV,NW,IGUESS,P,DP,DT,&
      &              NONODS,FREDOP,&
      &              C1T,C2T,C3T,FINDCT,COLCT,NCT, &
      &              ML,BIGM1,BIGM1,BIGM1,&
      &              FINDRM,COLM,CENTRM,NCOLM,NBIGM, &
      &              NDWISP,KCMC,&
      &              CMC,FINCMC,COLCMC,NCMC,&
      &              NCMCB,NPRESS,NPROPT,RADISO,ALFST2,SPRESS,&
      &              DM1,DM2,DM3,VECX,VECY,VECZ,&
      &              RHS,&
      &              UZAWA,CONMAS,COUPLE,D2HALF,ID3, &
      &              MULPA,CGSOLQ,GMRESQ,&
      &              POISON,PROJEC, &
      &              NDPSET,&
      &              PARA,halo_tag,&
      &              halo_tag_p,NNODP,NNODPP,&
      !     The following is for rotations only ROTAT=.TRUE.
      &              ROTAT,NNODRO,NRTDR,   NODROT,&
      &              NORMX,NORMY,NORMZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
      &              RTDR,ROTMOM,&
      !     The multiphase stuff...
      &              NPHASE,IPHASE,&
      &              C1TP,C2TP,C3TP, &
      &              NBUOY,BOY_ML, &
      !     For traffic 
      &              VOLTRAF,UTRAF,VTRAF,WTRAF,GOTTRAF,&
      !     --------------------------------------START OF ADD BY CRGW 14/03/06
      !     SINGLE PHASE SOLID COMPRESSIBILITY:
      &              MKCOMP, divqs, &
      &              state)
           !     ---------------------------END OF ADD BY CRGW 14/03/06
      ewrite(1,*)'just outside solnav'
      
      ! Copy viscosity back in case it's changed.
      if (has_tensor_field(state(1),"Viscosity")) then
         call copy_tensor_to_legacy(viscosity, MUPTXX, MUPTXY, MUPTXZ, &
              MUPTYY, MUPTYZ, MUPTZZ)
      end if

      call deallocate(ns_rhs)
      call deallocate(masslump)
      call deallocate(compressible_projec_rhs)
      call deallocate(big_m)
      call deallocate(ctp_m)
      
      if(have_gravity) then
        deallocate(soxgra)
        deallocate(soygra)
        deallocate(sozgra)
      end if
            
    END SUBROUTINE NAVSTO
      
      
    SUBROUTINE BIGPHA(PHAMAT,BIGTEM,NPHAMA,NCOLM,&
     &     MPDIAG,NMPDIA, MLPHA, CENTRM,&
     &     NPHASE,NONODS,D3, &
     &     ADDSOX,ADDSOY,ADDSOZ, VECX,VECY,VECZ)
!     This sub compiles the large matrix equation that results 
!     from 2 smaller matrices AND THE 4 matrices in ADDMAT.
!     WE MUST TAKE CARE AS THE MATRICES BIGM1,BIGM2,ADDMAT 
!     ARE ENTIRELY CONTAINED IN THE STORAGE SPACE ALOCATED TO PHAMAT. 
!     NB NADDMA=NCOLM*NPHASE*NPHASE. 
!     NB NMPDIA=NONODS*NPHASE*NPHASE 
!     Add MLPHA to MPDIAG(NB. MLPHA is already in matrix PHAMAT)
!     BIGTEM is a working array.
      use FLDebug
      IMPLICIT NONE
      INTEGER NPHAMA,NCOLM,NPHASE,NONODS
      REAL PHAMAT(NPHAMA),BIGTEM(4*NCOLM)
      REAL MLPHA(NPHASE*NONODS)
      INTEGER NMPDIA
      REAL MPDIAG(NMPDIA)
      INTEGER CENTRM(NONODS)
      REAL ADDSOX(NONODS*NPHASE),ADDSOY(NONODS*NPHASE)
      REAL ADDSOZ(NONODS*NPHASE)
      REAL VECX(NONODS*NPHASE),VECY(NONODS*NPHASE)
      REAL VECZ(NONODS*NPHASE)
      LOGICAL D3
!     Local variables...
      INTEGER NDIM,IMAT,I,NORD,INOD
      INTEGER IM,JM,ICENT,NOD,IP,ID

      ewrite(1,*) 'In bigpha'
!     Form the BIG matrix from all the other matrices -stored 
!     temporarily in PHAMAT.(ADDMAT,BIGM1,BIGM2 are in PHAMAT). 
      CALL MAKBIG(PHAMAT,BIGTEM,NPHAMA,NCOLM,D3)
      
      ewrite(1,*) 'After makbig'

!     Copy the diagonals in MPDIAG to PHAMAT...
      NDIM=2
      IF(D3) NDIM=3
      NORD=NPHASE*NDIM
      do IM=1,NORD
         do JM=1,NORD
            do NOD=1,NONODS
               ICENT=CENTRM(NOD)
               IMAT=ICENT + ((IM-1)*NORD+JM-1)*NCOLM
               INOD=NOD   + ((IM-1)*NORD+JM-1)*NONODS
               PHAMAT(IMAT)=PHAMAT(IMAT)+MPDIAG(INOD)
            END DO
         END DO
      END DO

!     Add ADDSOX,ADDSOY into rhs of momentum eqn's.
      do I=1,NONODS*NPHASE
         VECX(I)=VECX(I)+ADDSOX(I)
         VECY(I)=VECY(I)+ADDSOY(I)
         IF(D3) VECZ(I)=VECZ(I)+ADDSOZ(I)
      END DO
  
!     Add MLPHA to MPDIAG(NB. MLPHA is already in matrix PHAMAT)
      do IP=1,NPHASE
         do ID=1,NDIM
            IM=(IP-1)*NDIM+ID
            do NOD=1,NONODS
               INOD=NOD   + ((IM-1)*NORD+IM-1)*NONODS
               MPDIAG(INOD)=MPDIAG(INOD)+MLPHA(NOD+(IP-1)*NONODS)
            END DO
         END DO
      END DO

      ewrite(1,*)'Leaving bigpha'

    END SUBROUTINE BIGPHA

    SUBROUTINE MAKBIG(PHAMAT,BIGTEM,NPHAMA,NCOLM,D3)
      use FLDebug
      IMPLICIT NONE
      INTEGER NPHAMA,NCOLM,I
      REAL PHAMAT(NPHAMA),BIGTEM(4*NCOLM)
      LOGICAL D3
!     Form the BIG matrix from all the other matrices - stored.  
!     temporarily in PHAMAT.
      ewrite(1,*) 'In makbig, D3',D3
      ewrite(2,*) 'NPHAMA',NPHAMA
!     
      IF(.NOT.D3) THEN
!     Before manipulating matrix store contribution for gas.
         BIGTEM(1:1 + 2*NCOLM - 1) = PHAMAT(4*NCOLM+1:4*NCOLM+1 + 2*NCOLM - 1)
         BIGTEM(2*NCOLM+1:2*NCOLM+1 + 2*NCOLM - 1) = PHAMAT(10*NCOLM+1:10*NCOLM+1 + 2*NCOLM - 1)

!     PUT THE 2 BIG MATRICES INTO PHAMAT ****************
!     PUT TOP LEFT MATIRX TOGETHER. 
         PHAMAT(4*NCOLM+1:4*NCOLM+1 + 2*NCOLM - 1) = PHAMAT(2*NCOLM+1:2*NCOLM+1 + 2*NCOLM - 1)
         PHAMAT(2*NCOLM+1:2*NCOLM+1 + 2*NCOLM - 1) = 0.0
!     PUT BOTTOM RIGHT MATRIX TOGETHER. 
         PHAMAT(10*NCOLM+1:10*NCOLM+1 + 2*NCOLM - 1) = PHAMAT(12*NCOLM+1:12*NCOLM+1 + 2*NCOLM - 1)
         PHAMAT(12*NCOLM+1:12*NCOLM+1 + 2*NCOLM - 1) = 0.0
!     ****************************************************

!     Place ADDMAT matrices over the block diagonals of PHAMAT. *****
!     TOP LEFT BLOCKS...
         CALL RADDIN(PHAMAT(        1),PHAMAT(6*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(5*NCOLM+1),PHAMAT(6*NCOLM+1),NCOLM)
         PHAMAT(6*NCOLM+1:6*NCOLM+1 + NCOLM - 1) = 0.0
!     BOTTOM RHIGHT BLOCKS...
         CALL RADDIN(PHAMAT(10*NCOLM+1),PHAMAT(9*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(15*NCOLM+1),PHAMAT(9*NCOLM+1),NCOLM)
         PHAMAT(9 *NCOLM+1:9 *NCOLM+1 + NCOLM - 1) = 0.0
!     ****************************************************************

!     NOW PUT OFF DIAGONAL BLOCKS IN***
         PHAMAT(2 *NCOLM+1:2 *NCOLM+1 + NCOLM - 1) = PHAMAT(7*NCOLM+1:7*NCOLM+1 + NCOLM - 1)
         PHAMAT(13*NCOLM+1:13*NCOLM+1 + NCOLM - 1) = PHAMAT(8*NCOLM+1:8*NCOLM+1 + NCOLM - 1)
!     *********************************

!     Add contibution to PHAMAT stored in BIGTEM
         CALL RADDIN(PHAMAT(10*NCOLM+1),BIGTEM(1),        2*NCOLM)
         CALL RADDIN(PHAMAT(14*NCOLM+1),BIGTEM(2*NCOLM+1),2*NCOLM)
      ENDIF

      IF(D3) THEN
!     Copy memory used in MATPHA into bottom of PHAMAT
         CALL RADDIN(PHAMAT(27*NCOLM+1),PHAMAT(15*NCOLM+1),9*NCOLM)
!     Construct the bottom right matrix.
         ewrite(3,*) 'here 1'
         do I=1,3*NCOLM
            PHAMAT(21*NCOLM+I) = PHAMAT(27*NCOLM+I)
         ENDDO
         do I=1,3*NCOLM
            PHAMAT(27*NCOLM+I) = PHAMAT(30*NCOLM+I)
         ENDDO

!     Clear bottom left matrix.
         ewrite(3,*) 'here 3'
         PHAMAT(18*NCOLM+1:18*NCOLM+1 + 3*NCOLM - 1) = 0.0
         PHAMAT(24*NCOLM+1:24*NCOLM+1 + 3*NCOLM - 1) = 0.0
         PHAMAT(30*NCOLM+1:30*NCOLM+1 + 3*NCOLM - 1) = 0.0
!     Finish off bottom of matrix by adding ADDMAT...
         CALL RADDIN(PHAMAT(18*NCOLM+1),PHAMAT(11*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(25*NCOLM+1),PHAMAT(11*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(32*NCOLM+1),PHAMAT(11*NCOLM+1),NCOLM)
         PHAMAT(11*NCOLM+1:11*NCOLM+1 + NCOLM - 1) = 0.0
         CALL RADDIN(PHAMAT(21*NCOLM+1),PHAMAT(12*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(28*NCOLM+1),PHAMAT(12*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(35*NCOLM+1),PHAMAT(12*NCOLM+1),NCOLM)
         PHAMAT(12*NCOLM+1:12*NCOLM+1 + NCOLM - 1) = 0.0
!     
!     Construct top left matrix...
         do I=1,3*NCOLM
            PHAMAT(12*NCOLM+I) = PHAMAT(6*NCOLM+I)
         ENDDO

         do I=1,3*NCOLM
            PHAMAT(6*NCOLM+I) = PHAMAT(3*NCOLM+I)
         ENDDO

!     Add ADDMAT into top left matrix
         CALL RADDIN(PHAMAT(0*NCOLM+1), PHAMAT(9*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(7*NCOLM+1), PHAMAT(9*NCOLM+1),NCOLM)
         CALL RADDIN(PHAMAT(14*NCOLM+1),PHAMAT(9*NCOLM+1),NCOLM)
!     Top right matrix construct...
         PHAMAT(3*NCOLM+1:3*NCOLM+1 + 3*NCOLM - 1) = 0.0
         PHAMAT(9*NCOLM+1:9*NCOLM+1 + 1*NCOLM - 1) = 0.0
         PHAMAT(11*NCOLM+1:11*NCOLM+1 + 1*NCOLM - 1) = 0.0
         PHAMAT(15*NCOLM+1:15*NCOLM+1 + 3*NCOLM - 1) = 0.0
         PHAMAT(3*NCOLM+1:3*NCOLM+1 + NCOLM - 1) = PHAMAT(10*NCOLM+1:10*NCOLM+1 + NCOLM - 1)
         PHAMAT(17*NCOLM+1:17*NCOLM+1 + NCOLM - 1) = PHAMAT(10*NCOLM+1:10*NCOLM+1 + NCOLM - 1)
      ENDIF
      ewrite(1,*)'Leaving makbig'

    END SUBROUTINE MAKBIG

    SUBROUTINE MARDMI(BIGM1,NBIGM, &
     &     CENTRM,FINDRM,COLM,&
     &     NCOLM,NONODS,NPHASE,NDIM,&
     &     DM1PHA,DM2PHA,DM3PHA,&
     &     NORMX,NORMY,NORMZ, &
     &     T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &     NODROT,NNODRO,ROTAT,D3, &
     &     MAWORK)
!     Rotate the matrix BIGM1 then add in the b.c's (DM1PHA,DM2PHA,DM3PHA)...
      IMPLICIT NONE
      INTEGER NBIGM,NCOLM,NPHASE,NONODS,NDIM
      INTEGER FINDRM(NONODS+1),COLM(NCOLM)
      INTEGER CENTRM(NONODS)
      REAL BIGM1(NBIGM)
      REAL DM1PHA(NPHASE*NONODS),DM2PHA(NPHASE*NONODS)
      REAL DM3PHA(NPHASE*NONODS)
      LOGICAL ROTAT,D3
      REAL MAWORK(NCOLM*NDIM*NDIM)
      INTEGER NNODRO
      REAL NORMX(NNODRO),NORMY(NNODRO),NORMZ(NNODRO)
      REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
      REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
      INTEGER NODROT(NNODRO)
!     Local variables...
      INTEGER IP,JP,II

      do IP=1,NPHASE
         do JP=1,NPHASE

         IF(ROTAT.OR.(IP.EQ.JP)) THEN
!     Take MAWORK out of BIGM1...
            do II=1,NDIM
               MAWORK(1+(II-1)*NDIM*NCOLM:1+(II-1)*NDIM*NCOLM +&
     &                 NCOLM*NDIM - 1) = BIGM1(1+ (IP-1)*NCOLM*NDIM*NDIM*NPHASE&
     &                 +(II-1)*NCOLM*NDIM*NPHASE+ (JP-1)*NCOLM*NDIM :1 + &
     &                 (IP-1)*NCOLM*NDIM*NDIM*NPHASE +(II-1)*NCOLM*NDIM*NPHASE + &
     &                 (JP-1)*NCOLM*NDIM  + NCOLM*NDIM - 1)
            END DO
!     Rotate each quadrant of the matrix BIGM1...
            IF(ROTAT) THEN
               CALL RBIGRT(MAWORK,CENTRM,FINDRM,COLM,&
     &                 NCOLM,NCOLM*NDIM*NDIM,NONODS,&
     &                 NORMX,NORMY,NORMZ, &
     &                 T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &                 NODROT,NNODRO,D3)
            ENDIF

            IF(IP.EQ.JP) THEN
!     Add in the b.c's
               CALL BOUVIS(NDIM, NONODS, &
     &               DM1PHA(1+(IP-1)*NONODS),DM2PHA(1+(IP-1)*NONODS),&
     &               DM3PHA(1+(IP-1)*NONODS),&
     &               MAWORK,NCOLM*NDIM*NDIM,NCOLM,CENTRM)
            ENDIF

!     Put MAWORK back into BIGM1...
            do II=1,NDIM
               BIGM1(1+ (IP-1)*NCOLM*NDIM*NDIM*NPHASE &
     &               +(II-1)*NCOLM*NDIM*NPHASE+ (JP-1)*NCOLM*NDIM :1 &
     &               + (IP-1)*NCOLM*NDIM*NDIM*NPHASE + &
     &               (II-1)*NCOLM*NDIM*NPHASE+ (JP-1)*NCOLM*NDIM + &
     &               NCOLM*NDIM - 1) = &
     &               MAWORK(1+(II-1)*NDIM*NCOLM:1+(II-1)*NDIM*NCOLM&
     &               + NCOLM*NDIM - 1)
            END DO

!     ENDOF IF(ROTAT.OR.(IP.EQ.JP)) THEN...
          ENDIF
       END DO
      END DO

    END SUBROUTINE MARDMI
        
    subroutine free_surface_diagnostics(state)
!
    type(state_type), intent(in) :: state
    type(scalar_field), pointer :: free_surface_field
    type(scalar_field), pointer :: max_over_time_free_surface_field
    type(scalar_field), pointer :: min_over_time_free_surface_field
    type(scalar_field), pointer :: range_free_surface_field
    type(scalar_field), pointer :: harmonic_current    
    character(len=OPTION_PATH_LEN) :: option_path  
    integer :: stat_free, stat_range, stat_max, stat_min, stat, stride, M, N, ii
    real :: spin_up_time, current_time
    integer, save :: harmonic_counter=0, total_counter=0
    real, save :: harmonic_times(50) !
    real :: sigma(11)
    character(len=3), dimension(11) :: constituents=(/"M2 ","S2 ","N2 ","K2 ","K1 ","O1 ","P1 ","Q1 ","Mf ","Mm ","SSa"/)            
    logical :: which_constituents(11)
!    
    ewrite(3,*),'in free_surface_diagnostics'
!
    free_surface_field => extract_scalar_field(state,"FreeSurface",stat_free)
    max_over_time_free_surface_field => extract_scalar_field(state, "MaxFreeSurface",stat_max)       
    min_over_time_free_surface_field => extract_scalar_field(state, "MinFreeSurface",stat_min) 
    range_free_surface_field => extract_scalar_field(state, "TidalRange",stat_range)
!
    if(stat_free /= 0) then
      ewrite(1,*)'I do not have a FreeSurface field so can not calculate diagnostics on it'
      return
    end if
    if(stat_max /= 0) then
      ewrite(1,*)'You have not defined a MaxFreeSurface field so jumping out of free_surface_diagnostics'
      return
    end if
    if(stat_min /= 0) then
      ewrite(1,*)'You have not defined a MinFreeSurface field so jumping out of free_surface_diagnostics'
      return
    end if
    if(stat_range /= 0) then
      ewrite(1,*)'You have not defined a TidalRange field so jumping out of free_surface_diagnostics'
      return
    end if
!           
    option_path = range_free_surface_field%option_path
    call get_option(trim(option_path)//'/diagnostic/spin_up_time', spin_up_time)
    ewrite(3,*)'spin_up_time:',spin_up_time
    call get_option("/timestepping/current_time", current_time)
    if(current_time<spin_up_time) then
       max_over_time_free_surface_field%val = 0.0
       min_over_time_free_surface_field%val = 0.0
       range_free_surface_field%val = 0.0
    else         
       max_over_time_free_surface_field%val = max( max_over_time_free_surface_field%val, free_surface_field%val )   
       min_over_time_free_surface_field%val = min( min_over_time_free_surface_field%val, free_surface_field%val )      
       range_free_surface_field%val = max_over_time_free_surface_field%val - min_over_time_free_surface_field%val 
    endif
!    

! Harmonic analysis stuff
    sigma(1:11)=0  ! eventually we should set the frequencies via diamond, assume for now we have less than 11 consituents.
    which_constituents(1:11) = .false.
    M = 0   ! how many consituents are we going to analyse?
    N = 50
    do ii = 1,11
      if (    has_scalar_field(state,'HarmonicAmplitude'//constituents(ii) ) &
         .or. has_scalar_field(state,'HarmonicPhase'//constituents(ii)) ) then
        M = M+1 
        sigma(M) = get_tidal_frequency(constituents(ii))/(2.*3.141592654)
        which_constituents(ii) = .true.
      endif      
    end do 

    ewrite(3,*)'Found this many constituents to analyse:', M, which_constituents, sigma
        
    if ( M .gt. 0 ) then
       total_counter = total_counter + 1  ! don't do this every time step so keep a counter
       stride = 5 ! add data value every stride timesteps ! this number will be an option from diamond, hard code to a reasonable value for now.
       if( mod(total_counter,stride)==0 ) then
          harmonic_counter = harmonic_counter+1

          harmonic_current => extract_scalar_field(state,'harmonic'//int2str(mod(harmonic_counter-1,N)+1),stat)
          call set(harmonic_current,free_surface_field)
          harmonic_times(mod(harmonic_counter-1,N)+1) = current_time

          if( harmonic_counter .ge. N ) then ! wait until we have filled up the time_series
             call harmonic_analysis(state, harmonic_times, harmonic_counter, sigma, M, which_constituents, node_count(free_surface_field))
          end if

       end if     
    end if
    

    end subroutine free_surface_diagnostics

!!!!!

    subroutine harmonic_analysis(state, harmonic_times, harmonic_counter, sigma, M, which_constituents, nonods)
    type(state_type), intent(in) :: state
    real, intent(in) :: harmonic_times(50), sigma(11)
    integer, intent(in) :: harmonic_counter, M, nonods
    logical, intent(in) :: which_constituents(11)
    type(scalar_field), pointer :: HarmonicAmplitude, HarmonicPhase, harmonic_current   
    integer :: stat, N, i, j, k, node, ii, MM
    real :: phase    
    real, dimension(:,:), allocatable :: harmonic_A  ! for solving Ax=b system
    real, dimension(:), allocatable :: harmonic_x, harmonic_b, harmonic_time_series_vals_at_node, harmonic_times_reordered
    character(len=3), dimension(11) :: constituents=(/"M2 ","S2 ","N2 ","K2 ","K1 ","O1 ","P1 ","Q1 ","Mf ","Mm ","SSa"/)            


! length of time series    
    N = 50     
    
    allocate(harmonic_A(2*M+1,2*M+1))
    allocate(harmonic_x(2*M+1))
    allocate(harmonic_b(2*M+1))
    allocate(harmonic_time_series_vals_at_node(N))
    allocate(harmonic_times_reordered(N))

!Extract the times (do a reordering)
      do i = 1,N
        harmonic_times_reordered(i) = harmonic_times(mod( (mod(harmonic_counter-1,N)+1)+i-1, N)+1)
      end do        

! Need to analyse for evary node - this can be optimised to only do at the surface      
    do node = 1,nonods


!Extract the free surface elevations at the current node
      do i = 1,N
        harmonic_current => extract_scalar_field(state,'harmonic'//int2str( mod( (mod(harmonic_counter-1,N)+1)+i-1, N)+1 ),stat)
        harmonic_time_series_vals_at_node(i) = node_val(harmonic_current,node) 
      end do    

! Form and invert the least squares system
      call harmonic_analysis_at_single_node(N,harmonic_times_reordered,harmonic_time_series_vals_at_node,M,sigma,&
                                                harmonic_A,harmonic_x,harmonic_b)

! stick the amplitude and phase into something that will be output
       MM = 0
       do ii = 1,11
         if ( which_constituents(ii) ) then
           MM = MM+1    
           if ( has_scalar_field(state,'HarmonicAmplitude'//constituents(ii)) ) then
             HarmonicAmplitude => extract_scalar_field(state,'HarmonicAmplitude'//constituents(ii),stat)
             if (stat == 0) call set( HarmonicAmplitude, node, sqrt( harmonic_x(MM+1)**2 + harmonic_x(MM+1+M)**2 ) )        
           end if
           if ( has_scalar_field(state,'HarmonicPhase'//constituents(ii)) ) then
             HarmonicPhase => extract_scalar_field(state,'HarmonicPhase'//constituents(ii),stat)
             phase = atan2(harmonic_x(MM+1+M),harmonic_x(MM+1))
             !*180.0/pi
             !if (phase < 0.0) phase = phase + 360.0
             if (stat == 0) call set( HarmonicPhase,     node, phase )       
           end if                               
         end if      
       end do 

    end do ! node loop
    
    deallocate(harmonic_A)
    deallocate(harmonic_x)
    deallocate(harmonic_b)    
    deallocate(harmonic_time_series_vals_at_node)
    deallocate(harmonic_times_reordered)
        
    end subroutine harmonic_analysis

!!!!!!!!
    subroutine harmonic_analysis_at_single_node(N,harmonic_times_reordered,harmonic_time_series_vals_at_node,M,sigma,&
                                                harmonic_A,harmonic_x,harmonic_b)
    real, intent(in) :: harmonic_times_reordered(:),harmonic_time_series_vals_at_node(:),sigma(:)
    real, intent(inout) :: harmonic_A(:,:),harmonic_x(:),harmonic_b(:)
    integer, intent(in) :: M,N
    real :: C_k, S_k, CC_jk, SS_jk, SC_jk, CS_kj, phase
    real :: pi = 3.141592654
    integer :: i, j, k, stat
    

! For the least squares system    
       harmonic_A(1,1) = N

! Need the C_k and S_k for first row/column
       j = 1
       do k = 1,M    
          C_k = 0.0
          S_k = 0.0
          do i = 1,N
             C_k = C_k + cos(2.*pi*sigma(k)*harmonic_times_reordered(i))
             S_k = S_k + sin(2.*pi*sigma(k)*harmonic_times_reordered(i))
          end do
          harmonic_A(1,k+1)   = C_k
          harmonic_A(1,k+1+M) = S_k
          harmonic_A(k+1,1)   = C_k
          harmonic_A(k+1+M,1) = S_k
       end do
! rest of the rows and columns of the matrix
       do j = 1,M
         do k = 1,M
          CC_jk = 0.0
          SS_jk = 0.0
          SC_jk = 0.0
          do i = 1,N
             CC_jk = CC_jk + cos(2.*pi*sigma(k)*harmonic_times_reordered(i))*cos(2.*pi*sigma(j)*harmonic_times_reordered(i))
             SS_jk = SS_jk + sin(2.*pi*sigma(k)*harmonic_times_reordered(i))*sin(2.*pi*sigma(j)*harmonic_times_reordered(i))
             SC_jk = SC_jk + cos(2.*pi*sigma(k)*harmonic_times_reordered(i))*sin(2.*pi*sigma(j)*harmonic_times_reordered(i))
          end do 
          CS_kj = SC_jk         
          harmonic_A(j+1,k+1)     = CC_jk ! top left quadrant
          harmonic_A(j+1+M,k+1+M) = SS_jk ! bottom right quadrant
          harmonic_A(j+1+M,k+1)   = SC_jk ! bottom left quadrant 
          harmonic_A(k+1,j+1+M)   = SC_jk ! top right quadrant  (swap order we fill up matrix as CS_kj.eq.SC_jk, CS_kj.ne.SC_kj)               
         end do
       end do
! now the rhs vector
       harmonic_b(1:2*M+1) = 0.0

       do i = 1,N
          harmonic_b(1) = harmonic_b(1) + harmonic_time_series_vals_at_node(i) 
       end do
   
       do j = 1,M
          do i = 1,N
             harmonic_b(j+1)   = harmonic_b(j+1)   + harmonic_time_series_vals_at_node(i)*cos(2.*pi*sigma(j)*harmonic_times_reordered(i))
             harmonic_b(j+1+M) = harmonic_b(j+1+M) + harmonic_time_series_vals_at_node(i)*sin(2.*pi*sigma(j)*harmonic_times_reordered(i))
          end do    
       end do
! solve the system
       call solve(harmonic_A, harmonic_b, stat) ! Solve Ax=b, note that b will be overwritten
       harmonic_x = harmonic_b


    end subroutine harmonic_analysis_at_single_node



  end module navsto_module
      
      
