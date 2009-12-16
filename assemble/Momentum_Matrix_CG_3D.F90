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

  module diff3d_module
    
  use sml
  use AllSorts
  use Coordinates
  use Element_Numbering
  use FETools
  use FLDebug
  use Global_Parameters, only: OPTION_PATH_LEN
  use hart3d_allsorts
  use legacy_cv_numbering
  use legacy_cv_shape_functions
  use limit_gradient
  use mesh_connections
  use position_in_matrix
  use shape_transformations
  use sparse_tools
  use Solvers
  use tr2d_module
  use rotated_boundary_conditions_legacy
  use spud
  use flcomms_module
  use detnlxr_module
  use spaerr_module
  use coriolis_module
 
  implicit none

  private

  public :: ano_ocean_scal, get_quad_loc, half_lesvis_hart,&
       &  high_visc, lesderivs, lesela, lesvis, lesvis2, oneeletens,&
       &  sgsproj, sizegieletens, solv_field_sgsdg, &
       &  diff3d, diff3dsimp

contains

  SUBROUTINE DIFF3D &
  &         (acctim,U,V,W,NU,NV,NW,UG,VG,WG,&
  &          SOURCX,SOURCY,SOURCZ,X,Y,Z, D0,&
  &          VECX,VECY,VECZ, &
  &          BIGM,NBIGM, &
  &          ML,&
  &          TOTELE, &
  &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
  &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
  &          DISOPT,DISOPN,DT,THETA1,BETA,&
  &          GEOBAL,&
  &          LUMP,MAKSYM,CONVIS, &
  &          ABSLUM,SLUMP,&
  &          NDGLNO, &
  &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
  &          DENPT,&
  &          MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
  &          ML2MXX,ML2MXY,ML2MXZ,ML2MYY,ML2MYZ,ML2MZZ, &
  &          XABSOR,YABSOR,ZABSOR, &
  &          CVMVEC,&
  &          LES4TH,&
  &          LESNOD,&
  &          SCFACTH0,&
!     stuff we need for LANS
  &          nnodp,para,halo_tag,&
!     traffic stuff
  &          GOTTRAF,VERSIO,ISPHERE)
! This subroutine forms the discretised momentum equations. 
! The discretisation scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
! All methods use a balancing diffusion term and THETA time stepping. 
! The following are absolute values of DISOPT...
! DISOPT=1 - balancing diffusion based on (x,y) space.
! DISOPT=2 - Laxwendrof balancing diffusion.
! DISOPT=3 - (x,y,t) -balancing diffusion. 
! DISOPT=4 - No balancing diffusion.
! DISOPT=5 - nonlinear streamline and cross stream diffusion.
! DISOPT=6 - nonlinear upwind in steapest direction.
! DISOPT=7 - nonlinear streamline+ cross stream diffusion(but restricted)

! DISOPT=42- LES option using constant length scale.
! DISOPT=43- LES option using isotropic length scale.
! DISOPT=44- LES option which uses no balancing diffusion.
! DISOPT=45- LES option which uses no balancing diffusion.
! DISOPT=46- same as 45 but with 4th order dissipation.
! DISOPT=47 -LES but in tensor form like hart3d.
! DISOPT=48 -LES 4th order version of 47.

! DISOPT=125 - NO balancing diffusion(DISOPT=4)and take out non-linear terms.

! if DISOPN.ne.0 then set advection to zero and treat advection 
! using the high res method.
! If LES4th=1 then wehn DISOPT=45 USE A 4TH order Smagorinsky. 
!    with this option it is best to use THETA=0.5, ITHETA=0.5 and ITINOI>=2
! SORN contains the div q field NOT N_i * div q. (IS NOT USED NOW)
! GEOBAL=option to use a Geostrophic balance solver instead of 
! putting corriolis in here.
! IF SYM then make BIGM symmetrix .
! IF LUMP then lump the mass matrix else dont. 
! BETA controls the conservative discretisation 
! BETA=0. is the usual advection form, BETA=1. is divergence form. 
! BETA=0.5 is an average of the divergence and advection form.(recommended).
! BETA\in[-2,-1] then use the integrated by parts form of the advection term.
! this must be used with PREOPT=1 because we do not want pressure  
! BETA=-1 is non-conservative form
! BETA=-2 is conservative form (recommended)
! BETA=-1.5 is a compromise.
! ABSORBX =absorbtion of energy of x-momentum equation. 
! CVMVEC=virtual mass vector (is 0 in non 2-phase flow). 
! MUPTXX,MUPTYY =components of diffusivities. 
! ML2MXX,ML2MXY etc contain the length scales for LES. 
! NU,NV,NW are for the non-linear terms ordinarily NU=U,NV=V,NW=W.
! UG,VG,WG are the grid velocities. 
! R0=(NOT NEEDED IN CARTEASEAN COORDS) 
! the minimum distance from ocean bed to centre of earth.
! D0=H/r0 where H is the height of ocean above R0 (H not used here)
! For dimensional run D0=0.
! NB the centrifugal force is \row*(distance from line (x,y)=(0,0)
! NDGLNO=element pter for unknowns. 
! SONDGL=element pter for materials and sources. 
! VONDGL=element pter for advection NU,UG.
! XONDGL=element pter for coordinates(x,y,z)
!
! If CONVIS then assume the conservation of viscosity due to div q term, 
! if div q is not zero then this must be set to .TRUE. 
! NU,NV,NW are for the non-linear terms ordinarily NU=U,NV=V,NW=W.
! UG,VG,WG are the grid velocities. 
! NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
! THE BLOCKS OF THE MATRIX BIGM ARE ARRANGED AS FOLLOWS
!   WHEN BLOCK SYMMETRIC ...
!   1 * *
!   2 3 *
!   4 5 6 
!   WHEN **NOT** BLOCK SYMMETRIC ...
!   1 2 3
!   4 5 6
!   7 8 9                     *- MATRIX NOT STORED.
    REAL RWIND
    PARAMETER(RWIND=0.5)
    INTEGER SNONOD,VNONOD,XNONOD
    INTEGER DISOPT,DISOPN,GEOBAL,DISCRA
    INTEGER NGI,NLOC,NBIGM,NCOLM
    REAL DENPT(SNONOD)
    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
    REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
    REAL ML2MXX(SNONOD),ML2MXY(SNONOD),ML2MXZ(SNONOD)
    REAL ML2MYY(SNONOD),ML2MYZ(SNONOD),ML2MZZ(SNONOD)
    REAL XABSOR(SNONOD),YABSOR(SNONOD),ZABSOR(SNONOD)
    REAL CVMVEC(SNONOD)
    REAL DT,THETA,BETA,acctim
    INTEGER TOTELE,NONODS
! If THETA=0.5 then Crank-Nickolson is used in Navier Stokes equ's.
! If THETA=1.0 then bakward Euler is used in NS equns.
! Similarly for the temperature equation. 
    REAL, intent(out) :: U(NONODS),V(NONODS),W(NONODS)
    REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
    REAL UD(NGI),VD(NGI),WD(NGI)
    REAL UDx(NGI),VDx(NGI),WDx(NGI)
    REAL UDy(NGI),VDy(NGI),WDy(NGI)
    REAL UDz(NGI),VDz(NGI),WDz(NGI)
    REAL HOVERQ
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL BECTWE(NGI)
    REAL DDETWE(NGI)
    REAL SCFACTH0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     stuff we need for LANS -- cjc
    integer, intent(in) :: nnodp,para,halo_tag
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    INTEGER, intent(in):: ISPHERE
      
! HX,HY,HZ-characteristic length scales in x,y,z directions.
    REAL HXGI,HYGI,HZGI
    REAL GAMMA(NGI)
    REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
    REAL MYZDTW(NGI),MZZDTW(NGI)
    REAL LXXDTW(NGI),LXYDTW(NGI),LXZDTW(NGI),LYYDTW(NGI)
    REAL LYZDTW(NGI),LZZDTW(NGI)
    REAL XABDGI(NGI),YABDGI(NGI),ZABDGI(NGI)
    REAL MUGIXX,MUGIXY,MUGIXZ,MUGIYY,MUGIYZ,MUGIZZ
    REAL L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ
    REAL DENGI(NGI)
    REAL E(6,6,NGI),OP(3,3),OP4TH(3,3),EE(6,6),EEN(6,6)

    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL BIGM(NBIGM)
    REAL ML(NONODS)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
    REAL MC1,HL1,HL2,HL3,VLND,VLM,VLMD,VLMOME
    REAL VLKXX,VLKXY,VLKXZ
    REAL VLKYX,VLKYY,VLKYZ
    REAL VLKZX,VLKZY,VLKZZ 
    REAL VXVXX,VXVXY,VXVXZ,VYVYX,VYVYY,VYVYZ
    REAL VZVZX,VZVZY,VZVZZ
    REAL VABSX,VABSY,VABSZ
    REAL VLKTEN
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER ELE,ILOC,JLOC,L
    INTEGER GI,COUNT,I,GLOBI,GLOBJ,GLOBIS,GLOBJS
    INTEGER LOWER,UPPER
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER CENTRM(NONODS)

    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
! For 4th order LES (LES4TH=1)
    INTEGER LES4TH
! if LESNOD=1 then have a node-wise material 
! for LES - nec. for proper 4th order LES.
    INTEGER LESNOD
    REAL, dimension(:,:,:), allocatable :: EV
    REAL UUD(NGI),UVD(NGI),UWD(NGI)
    REAL VUD(NGI),VVD(NGI),VWD(NGI)
    REAL WUD(NGI),WVD(NGI),WWD(NGI)
    REAL XD(NGI),YD(NGI),ZD(NGI)
    REAL GAMMAU(NGI),GAMMAV(NGI),GAMMAW(NGI)

    REAL TWOTHI
    REAL ABLUMP
    LOGICAL ABSLUM,SLUMP
    LOGICAL LUMP,MAKSYM,CONVIS,UPWIND,XTPETV
    LOGICAL BLKSYM
! Ocean stuff
    REAL D0
    LOGICAL NODIM
    REAL INROSB,ROSBY,INROSS,INDD0,DD0
    REAL COLUMP,RLUMP,RSLUMP,RSYM,RGAMMA,R,R1,R2,RNN
    REAL maxr,maxvln,maxhl,UPNAX,UPNAY,UPNAZ,UPNAA
    INTEGER::IGLS,IGLV=-1,IGLX,ICENT,INUM
    INTEGER IBL11,IBL12,IBL13
    INTEGER IBL21,IBL22,IBL23
    INTEGER IBL31,IBL32,IBL33
    LOGICAL NONLIN
    REAL R1U,R2U, R1V,R2V, R1W,R2W 
    REAL UU
    REAL HH2,ROTRHS
    REAL LOVERQ,BETA1,BETA2,ADV1,ADV2
    LOGICAL LES,LESTWO,TENSOR

!     local stuff for LANS-alpha model

    real, allocatable, dimension(:) :: smooth_u,smooth_v,smooth_w
    real, allocatable, dimension(:) :: smooth_nu,smooth_nv,smooth_nw
    real::LANS_nonlin_XX=0, LANS_nonlin_XY=0, LANS_nonlin_XZ=0
    real::LANS_nonlin_YX=0, LANS_nonlin_YY=0, LANS_nonlin_YZ=0
    real::LANS_nonlin_ZX=0, LANS_nonlin_ZY=0, LANS_nonlin_ZZ=0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! END: some stuff for LANS-alpha model - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    INTEGER GOTTRAF,II,JJ,KLOC,KNOD
    INTEGER VERSIO
    REAL ABTHETA
!     ADJOINT
    REAL THETA1
    LOGICAL TENELEWIS,TEN4TH
    REAL TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI) 
! Local coords...
    INTEGER NCLOC
    REAL, ALLOCATABLE, DIMENSION(:,:)::NC
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCLX
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCLY
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCLZ
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCX
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCY
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDZ
    logical, parameter :: using_lans = .false.
    integer, parameter :: nonlinearity_option = 0
    real:: f0

    THETA=THETA1

    ewrite(1,*) 'in subroutine diff3d' 
    ewrite(2,*) 'cjc maxval(sourcx)', maxval(abs(sourcx))

    allocate(EV(6,6,LESNOD*NONODS))
       
! Hvae an element-wise tensor calculated in each element
    TENELEWIS=.false.
    NCLOC=4
    if (NLOC == 8) NCLOC = 8

    ALLOCATE(NC(NCLOC,NGI))
    ALLOCATE(NCLX(NCLOC,NGI))
    ALLOCATE(NCLY(NCLOC,NGI))
    ALLOCATE(NCLZ(NCLOC,NGI))
    ALLOCATE(NCX(NCLOC,NGI))
    ALLOCATE(NCY(NCLOC,NGI))
    ALLOCATE(NCZ(NCLOC,NGI))
    IF(ISPHERE==0) THEN
       NC = N
       NCLX = NLX
       NCLY = NLY
       NCLZ = NLZ
    else
!        don't know what this is supposed to do 
!        (isphere==1 the same as ==0? No! it isn't!)
       FLAbort("Unexpected value for ISPHERE in DIFF3D")
    ENDIF

    EWRITE_MINMAX(ML2MXX)
    EWRITE_MINMAX(ML2MXY)
    EWRITE_MINMAX(ML2MXZ)
    EWRITE_MINMAX(ML2MYY)
    EWRITE_MINMAX(ML2MYZ)
    EWRITE_MINMAX(ML2MZZ)
      
    ewrite_minmax(denpt)
    ewrite(2,*)  "R2NORM(denpt, nonods, 1) = ", R2NORM(denpt, nonods, 1)
    ewrite(2,*)  "R2NORM(nu, nonods, 1) = ", R2NORM(nu, nonods, 1)
    ewrite(2,*)  "R2NORM(nv, nonods, 1) = ", R2NORM(nv, nonods, 1)
    ewrite(2,*)  "R2NORM(nw, nonods, 1) = ", R2NORM(nw, nonods, 1)
    ewrite(2,*) 'les4th,lesnod:',les4th,lesnod
       
    TEN4TH=(DISOPT.EQ.48)
       
    IF((LES4TH.EQ.1).OR.(TEN4TH)) THEN
! calculate the derivatives of NU,NV,NW wrt X,Y,Z
       ALLOCATE(VEDUDX(NONODS))
       ALLOCATE(VEDUDY(NONODS))
       ALLOCATE(VEDUDZ(NONODS))
       ALLOCATE(VEDVDX(NONODS))
       ALLOCATE(VEDVDY(NONODS))
       ALLOCATE(VEDVDZ(NONODS))
       ALLOCATE(VEDWDX(NONODS))
       ALLOCATE(VEDWDY(NONODS))
       ALLOCATE(VEDWDZ(NONODS))
       ewrite(1,*) 'just before LESDERIVS'
       CALL LESDERIVS(&
     &          VEDUDX,VEDUDY,VEDUDZ,&
     &          VEDVDX,VEDVDY,VEDVDZ,&
     &          VEDWDX,VEDWDY,VEDWDZ,&
     &          ML,&
     &          NU,NV,NW,&
     &          X,Y,Z, XONDGL, TOTELE,NONODS,NLOC,NGI, &
     &          N,NLX,NLY,NLZ, WEIGHT, DETWEI,DDETWE, .TRUE.,.FALSE.,&
     &          NX,NY,NZ, &
     &          NDGLNO, nnodp,para,halo_tag,VECX,VECY,&
     &          LESNOD,EV,&
     &          ML2MXX,ML2MXY,ML2MXZ,ML2MYY,ML2MYZ,ML2MZZ, &
     &          UG,VG,WG, &
     &          VONDGL)
       ewrite(1,*) 'just after LESDERIVS'
       CALL PMINMX(VEDUDX,NONODS,'******VEDUDX  ')
       CALL PMINMX(VEDUDY,NONODS,'******VEDUDY  ')
       CALL PMINMX(VEDUDZ,NONODS,'******VEDUDZ  ')
         
       CALL PMINMX(VEDVDX,NONODS,'******VEDVDX  ')
       CALL PMINMX(VEDVDY,NONODS,'******VEDVDY  ')
       CALL PMINMX(VEDVDZ,NONODS,'******VEDVDZ  ')
          
       CALL PMINMX(VEDWDX,NONODS,'******VEDWDX  ')
       CALL PMINMX(VEDWDY,NONODS,'******VEDWDY  ')
       CALL PMINMX(VEDWDZ,NONODS,'******VEDWDZ  ')
    ENDIF
    ewrite(2,*) 'using_lans:',using_lans
       
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     smoothing operator to nonlinear velocities for LANS - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
    if(using_lans.and.nonlinearity_option<3) then
       allocate( smooth_u(nonods), smooth_v(nonods), smooth_w(nonods) )
       allocate( smooth_nu(nonods), smooth_nv(nonods), smooth_nw(nonods) )
         
       smooth_nu = NU(1:nonods)
       smooth_nv = Nv(1:nonods)
       smooth_nw = Nw(1:nonods)
       smooth_u = U(1:nonods)
       smooth_v = v(1:nonods)
       smooth_w = w(1:nonods)
         
       ewrite(2,*) size(X), XNONOD
    end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     END: smoothing operator to nonlinear velocities for LANS - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! This subroutine forms a contabution to the Right Hand Side
! of Poissons pressure equation, as well as  F1 & F2.
! THIS SUBROUTINE FORMS THE MATRICES BIGM, AND C1T,C2T.

! THE STRUCTURE OF THE MATRIX ...
    IF(NBIGM.EQ.9*NCOLM) THEN
       BLKSYM=.FALSE.
    ELSE
       BLKSYM=.TRUE.
    ENDIF
    IF(BLKSYM) THEN
       IBL11=0
       IBL12=1*NCOLM
       IBL13=3*NCOLM
       IBL21=1*NCOLM
       IBL22=2*NCOLM
       IBL23=4*NCOLM
       IBL31=3*NCOLM
       IBL32=4*NCOLM
       IBL33=5*NCOLM
    ELSE
       IBL11=0
       IBL12=1*NCOLM
       IBL13=2*NCOLM
       IBL21=3*NCOLM
       IBL22=4*NCOLM
       IBL23=5*NCOLM
       IBL31=6*NCOLM
       IBL32=7*NCOLM
       IBL33=8*NCOLM
    ENDIF

    ML=0.0

    IF(GEOBAL.EQ.0) THEN
! Dont set to zero if Galerkin and visited GEOELI1P or other 
! geostrophic pressure solvers
       VECX=0.0
       VECY=0.0
       VECZ=0.0
    ENDIF
    BIGM(1:NBIGM)=0.0
      
    COLUMP=0.0
    RLUMP =0.
    ABLUMP=0.
    RSLUMP=0.
    RSYM  =1.
    TWOTHI=0.
    IF(LUMP)   RLUMP=1.
    IF(ABSLUM) ABLUMP=1.
    IF(SLUMP)  RSLUMP=1.
    IF(MAKSYM) RSYM=0.
    IF(CONVIS) TWOTHI=2./3.
! FORM OD THE ADVECTION TERM...BETA1,BETA2,ADV1,ADV2
    IF(BETA.GT.-0.5) THEN
! Standard form...
       BETA1=BETA
       BETA2=0.0
       ADV1=1.0
       ADV2=0.0
    ELSE
! have integrated advection terms by parts...
       BETA1=0.
       BETA2=ABS(BETA+1.0)
       ADV1=0.0
       ADV2=1.0
    ENDIF

    ABTHETA=THETA
    IF(GOTTRAF.EQ.1) ABTHETA=1.0

! Ocean scaling.
    NODIM=.TRUE.
    IF(D0.EQ.0) NODIM=.FALSE.
    ! workout omega:
    if (have_option("/physical_parameters/coriolis/sine_of_latitude")) then
      call get_option("/physical_parameters/coriolis/sine_of_latitude/omega", f0)
    else if (have_option("/physical_parameters/coriolis/on_sphere"))&
         & then
      call get_option("/physical_parameters/coriolis/on_sphere/&
            &omega", f0)
    else
      f0=0.0
    end if
    f0=f0/2.0 ! f=2omega
    INROSB=f0
    ROSBY =1.
    IF(f0.NE.0) ROSBY =1./INROSB
    INROSS=1.
    INDD0=1.
    DD0=1.
    IF(NODIM) THEN
       INROSS=INROSB
       INDD0=1./D0
       DD0=D0
    ENDIF
    IF(INROSS.EQ.0) INROSS=1.
    ROTRHS=1.0
    IF(GEOBAL.LE.-3) THEN
       ROTRHS=0.
       COLUMP=1.
    ENDIF

! START OF OPTIONS ******************************
    UPWIND=.FALSE.
    NONLIN=.FALSE.
! OPWIND=optimal upwinding method. 
    XTPETV=.FALSE.
    RGAMMA=0.0
! DISOPT=discretisation option 
    DISCRA=ABS(DISOPT)
    ewrite(2,*) 'DISOPT, DISCRA = ', DISOPT, DISCRA
    LES=(DISOPT.GE.42).AND.(DISOPT.LE.44)
    LESTWO=(DISOPT.GE.45).AND.(DISOPT.LE.48)
    TENSOR=LESTWO

    UUD(1:NGI) = 0.0
    UVD(1:NGI) = 0.0
    UWD(1:NGI) = 0.0
    VUD(1:NGI) = 0.0
    VVD(1:NGI) = 0.0
    VWD(1:NGI) = 0.0
    WUD(1:NGI) = 0.0
    WVD(1:NGI) = 0.0
    WWD(1:NGI) = 0.0
    XD(1:NGI) = 0.0
    YD(1:NGI) = 0.0
    ZD(1:NGI) = 0.0
    GAMMAU(1:NGI) = 0.0
    GAMMAV(1:NGI) = 0.0
    GAMMAW(1:NGI) = 0.0

    IF(DISCRA.EQ.1) THEN
! balancing diffusion based on (x,y) space.
       UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.2) THEN
! Laxwendrof balancing diffusion.
       UPWIND=.FALSE.
       RGAMMA=0.5*DT
    ENDIF
    IF(DISCRA.EQ.3) THEN
! (x,y,t) -balancing diffusion.
       UPWIND=.TRUE.
       XTPETV=.TRUE.
    ENDIF
    IF(DISCRA.EQ.4) THEN
! No balancing diffusion.
    ENDIF
    IF((DISCRA.GE.5).AND.(DISCRA.LE.7)) THEN
! nonlinear streamline and cross stream diffusion
       UPWIND=.TRUE.
       NONLIN=.TRUE.
    ENDIF

    do GI=1,NGI
       GAMMA(GI)=RGAMMA
    end do

    maxr=-1.e+20
    maxvln=-1.e+20
    maxhl=-1.e+20

! END OF OPTIONS *******************************

    do  ELE=1,TOTELE
! Get determinant and derivatives...
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
              & N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
              & NX,NY,NZ, NCX,NCY,NCZ,&
              & A11,A12,A13, A21,A22,A23, A31,A32,A33,&
              & XD,YD,ZD,&
              & ISPHERE) 
       IF(TENELEWIS) THEN
          CALL ONEELETENS(ELE,X,Y,Z,DISOPT,&
                 & TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ, &
                 & TOTELE,DETWEI,NGI,NLOC, &
                 & XONDGL,XNONOD)
       ENDIF

       do  GI=1,NGI
          UD(GI)=0.
          VD(GI)=0.
          WD(GI)=0.
            
          udx(gi) = 0.
          udy(gi) = 0.
          udz(gi) = 0.
          vdx(gi) = 0.
          vdy(gi) = 0.
          vdz(gi) = 0.
          wdx(gi) = 0.
          wdy(gi) = 0.
          wdz(gi) = 0.

          DENGI(GI)=0.
            
          MUGIXX=0.
          MUGIXY=0.
          MUGIXZ=0.
          MUGIYY=0.
          MUGIYZ=0.
          MUGIZZ=0.
            
          L2GIXX=0.
          L2GIXY=0.
          L2GIXZ=0.
          L2GIYY=0.
          L2GIYZ=0.
          L2GIZZ=0.

          XABDGI(GI)=0.
          YABDGI(GI)=0.
          ZABDGI(GI)=0.

          do  L=1,NLOC
             IGLS=SONDGL((ELE-1)*NLOC+L)
             IGLV=VONDGL((ELE-1)*NLOC+L)
             IGLX=XONDGL((ELE-1)*NLOC+L)
! NB R0 does not appear here although the z-coord might be Z+R0. 
! Introduce grid vels in non-linear terms.
             IF(DISOPN.EQ.0) THEN
                IF(DISOPT.NE.125) THEN
                   if(using_lans.and.nonlinearity_option.eq.2) then
                      UD(GI)=UD(GI) + N(L,GI)*(smooth_NU(IGLV)-UG(IGLV))
                      VD(GI)=VD(GI) + N(L,GI)*(smooth_NV(IGLV)-VG(IGLV))
                      WD(GI)=WD(GI) + N(L,GI)*(smooth_NW(IGLV)-WG(IGLV))
                   else
                      UD(GI)=UD(GI) + N(L,GI)*(NU(IGLV)-UG(IGLV))
                      VD(GI)=VD(GI) + N(L,GI)*(NV(IGLV)-VG(IGLV))
                      WD(GI)=WD(GI) + N(L,GI)*(NW(IGLV)-WG(IGLV))
                   end if
                ENDIF
             ENDIF
             DENGI(GI)=DENGI(GI)+ N(L,GI)*DENPT(IGLS)

             MUGIXX=MUGIXX+ N(L,GI)*MUPTXX(IGLS)
             MUGIXY=MUGIXY+ N(L,GI)*MUPTXY(IGLS)
             MUGIXZ=MUGIXZ+ N(L,GI)*MUPTXZ(IGLS)
             MUGIYY=MUGIYY+ N(L,GI)*MUPTYY(IGLS)
             MUGIYZ=MUGIYZ+ N(L,GI)*MUPTYZ(IGLS)
             MUGIZZ=MUGIZZ+ N(L,GI)*MUPTZZ(IGLS)

! ***********************************
! For LES...
             L2GIXX=L2GIXX+ N(L,GI)*ML2MXX(IGLS)
             L2GIXY=L2GIXY+ N(L,GI)*ML2MXY(IGLS)
             L2GIXZ=L2GIXZ+ N(L,GI)*ML2MXZ(IGLS)
             L2GIYY=L2GIYY+ N(L,GI)*ML2MYY(IGLS)
             L2GIYZ=L2GIYZ+ N(L,GI)*ML2MYZ(IGLS)
             L2GIZZ=L2GIZZ+ N(L,GI)*ML2MZZ(IGLS)

             XABDGI(GI)=XABDGI(GI) + N(L,GI)*XABSOR(IGLS)
             YABDGI(GI)=YABDGI(GI) + N(L,GI)*YABSOR(IGLS)
             ZABDGI(GI)=ZABDGI(GI) + N(L,GI)*ZABSOR(IGLS)

          end do

          IF((DISOPT.EQ.47).OR.(DISOPT.EQ.48)) THEN
         ! this has a QUADRATURE POINT varying element size and shape
             CALL SIZEGIELETENS(NX,NY,NZ,NLOC,NGI,GI,&
                  & L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ )
          ENDIF

          DDETWE(GI)=DENGI(GI)*DETWEI(GI)
          MXXDTW(GI)=MUGIXX*DETWEI(GI)
          MXYDTW(GI)=MUGIXY*DETWEI(GI)
          MXZDTW(GI)=MUGIXZ*DETWEI(GI)
          MYYDTW(GI)=MUGIYY*DETWEI(GI)
          MYZDTW(GI)=MUGIYZ*DETWEI(GI)
          MZZDTW(GI)=MUGIZZ*DETWEI(GI)
! for LES...
          LXXDTW(GI)=L2GIXX*DETWEI(GI)
          LXYDTW(GI)=L2GIXY*DETWEI(GI)
          LXZDTW(GI)=L2GIXZ*DETWEI(GI)
          LYYDTW(GI)=L2GIYY*DETWEI(GI)
          LYZDTW(GI)=L2GIYZ*DETWEI(GI)
          LZZDTW(GI)=L2GIZZ*DETWEI(GI)
          R=0.

          do  L=1,NLOC
             IGLV=VONDGL((ELE-1)*NLOC+L)
             if(using_lans.and.(nonlinearity_option<3)) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     use smooth velocities - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                R=R+NX(L,GI)*smooth_NU(IGLV)+NY(L,GI)&
                 &      *smooth_NV(IGLV)+ NZ(L,GI)*smooth_NW(IGLV)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     END: use smooth velocities - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             else
!     standard Navier-stokes
                R=R+NX(L,GI)*NU(IGLV)+NY(L,GI)*NV(IGLV)&
                     & + NZ(L,GI)*NW(IGLV)
             endif
          end do

          IF(LES.OR.(DISOPT.EQ.47).OR.(DISOPT.EQ.48)) THEN
! Amend viscosity (LXXDTW) for LES
! LXXDTW = contains the inverse of the metric 
! * DETWEI( determinant * quadrature weight )
             CALL LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             LXXDTW,LXYDTW,LXZDTW, &
     &             LYYDTW,LYZDTW,LZZDTW) 
! Add in the normal viscosity...
             MXXDTW(GI) = MXXDTW(GI)+LXXDTW(GI)
             MXYDTW(GI) = MXYDTW(GI)+LXYDTW(GI)
             MXZDTW(GI) = MXZDTW(GI)+LXZDTW(GI)
             MYYDTW(GI) = MYYDTW(GI)+LYYDTW(GI)
             MYZDTW(GI) = MYZDTW(GI)+LYZDTW(GI)
             MZZDTW(GI) = MZZDTW(GI)+LZZDTW(GI)
          ENDIF

          IF(LESTWO.AND.(.NOT.((DISOPT.EQ.47).OR.(DISOPT.EQ.48)))) THEN
             IF(LESNOD.NE.1) THEN
! Recalculate tensor for oceans...
                IF(TENELEWIS) THEN
                   L2GIXX=TENSXX
                   L2GIXY=TENSXY
                   L2GIXZ=TENSXZ
                   L2GIYY=TENSYY
                   L2GIYZ=TENSYZ
                   L2GIZZ=TENSZZ
                ENDIF

! Calculate E here
! E is the matrix of length-scales appropriate for each element

                CALL LESELA(NONODS,&
     &             L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ, &
     &             ELE,GI,NGI,NLOC,E)
         
! Mulitply E by 4*C_S^2*Sbars to give matrix of viscosities
                CALL LESVIS2(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             E)
             ENDIF
          ENDIF
         
          IF(DISOPN.EQ.0) THEN
             BECTWE(GI)=R*DETWEI(GI)
          ELSE
             BECTWE(GI)=0.0
          ENDIF
! Work out GAMMA at each Gauss pt(OR centre of element)
          IF(UPWIND) THEN
             HXGI=ABS(A11(GI)*UD(GI)+A12(GI)*VD(GI)+A13(GI)*WD(GI))
             HYGI=ABS(A21(GI)*UD(GI)+A22(GI)*VD(GI)+A23(GI)*WD(GI))
             HZGI=ABS(A31(GI)*UD(GI)+A32(GI)*VD(GI)+A33(GI)*WD(GI))

             IF(NONLIN) THEN

            ! (HH2 streamline direction)
                HH2=2./MAX(HXGI,HYGI,HZGI,1.E-7)
            ! limit for U...
                CALL LIMGRA(.TRUE.,NGI,NLOC,NONODS,TOTELE, &
     &               ELE,GI,&
     &               NDGLNO, UD,VD,WD, U, HH2, &
     &               N,NX,NY,NZ, (DISCRA.EQ.7),&
! 2-d Inverse of Jacobian...
     &               A11(GI),  A11(GI),A11(GI),&
     &               A11(GI),A11(GI),&
! 3-d inverse of Jacobian...
     &               A11(GI),A12(GI),A13(GI),&
     &               A21(GI),A22(GI),A23(GI),&
     &               A31(GI),A32(GI),A33(GI),&
! THE LIMITING VARIABLEs
     &               UUD(GI),UVD(GI),UWD(GI),LOVERQ )
                GAMMAU(GI)=RWIND*LOVERQ

! limit for V...
                CALL LIMGRA(.TRUE.,NGI,NLOC,NONODS,TOTELE, &
     &               ELE,GI,&
     &               NDGLNO, UD,VD,WD, V, HH2, &
     &               N,NX,NY,NZ, (DISCRA.EQ.7),&
! 2-d Inverse of Jacobian...
     &               A11(GI),  A11(GI),A11(GI),&
     &               A11(GI),A11(GI),&
! 3-d inverse of Jacobian...
     &               A11(GI),A12(GI),A13(GI),&
     &               A21(GI),A22(GI),A23(GI),&
     &               A31(GI),A32(GI),A33(GI),&
! THE LIMITING VARIABLEs
     &               VUD(GI),VVD(GI),VWD(GI),LOVERQ )
                GAMMAV(GI)=RWIND*LOVERQ

! limit for V...
                CALL LIMGRA(.TRUE.,NGI,NLOC,NONODS,TOTELE, &
     &               ELE,GI,&
     &               NDGLNO, UD,VD,WD, W, HH2, &
     &               N,NX,NY,NZ, (DISCRA.EQ.7),&
! 2-d Inverse of Jacobian...
     &               A11(GI),  A11(GI),A11(GI),&
     &               A11(GI),A11(GI),&
! 3-d inverse of Jacobian...
     &               A11(GI),A12(GI),A13(GI),&
     &               A21(GI),A22(GI),A23(GI),&
     &               A31(GI),A32(GI),A33(GI),&
! THE LIMITING VARIABLEs
     &               WUD(GI),WVD(GI),WWD(GI),LOVERQ )
                GAMMAW(GI)=RWIND*LOVERQ
               
             ENDIF

! HX,HY are the characteristic length scales in x,y directions. 
! If  XTPETV  then (x,y,t) Petrov Galerkin Method  else (x,y).
             IF(XTPETV) THEN
                HOVERQ=MIN(2./MAX(HXGI,HYGI,HZGI,1.E-7), DT)
             ELSE
                HOVERQ=2./MAX(HXGI,HYGI,HZGI,1.E-7)
             ENDIF
             GAMMA(GI)=RWIND*HOVERQ
          ENDIF
       end do

       do ILOC=1,NLOC
          GLOBI =NDGLNO((ELE-1)*NLOC+ILOC)
          GLOBIS=SONDGL((ELE-1)*NLOC+ILOC)
          do  JLOC=1,NLOC
             GLOBJ =NDGLNO((ELE-1)*NLOC+JLOC)
             GLOBJS=SONDGL((ELE-1)*NLOC+JLOC)

             HL1 =0.
             HL2 =0.
             HL3 =0.
             VLND=0.
             VLM =0.
             VLMD=0.
             VLMOME=0.
            
             VLKXX=0.
             VLKXY=0.
             VLKXZ=0.
             VLKYX=0.
             VLKYY=0.     
             VLKYZ=0. 
             VLKZX=0.
             VLKZY=0.
             VLKZZ=0.
            
             VXVXX=0.
             VXVXY=0.
             VXVXZ=0.
             VYVYX=0.
             VYVYY=0.
             VYVYZ=0.
             VZVZX=0.
             VZVZY=0.
             VZVZZ=0.
            
             VLKTEN=0.
            
             VABSX=0.
             VABSY=0.
             VABSZ=0.

             UPNAA=0.
             UPNAX=0.
             UPNAY=0.
             UPNAZ=0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     extra non-linear term for LANS-alpha - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

             if(using_LANS) then
                LANS_nonlin_XX=0.
                LANS_nonlin_XY=0.
                LANS_nonlin_XZ=0.
                LANS_nonlin_YX=0.
                LANS_nonlin_YY=0.
                LANS_nonlin_YZ=0.
                LANS_nonlin_ZX=0.
                LANS_nonlin_ZY=0.
                LANS_nonlin_ZZ=0.
             end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! END: extra non-linear term for LANS-alpha - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

             OP = 0.0
             OP4TH = 0.0
             do  GI=1,NGI
                R1= UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI) &
     &           + WD(GI)*NZ(ILOC,GI)
                R2=(UD(GI)*NX(JLOC,GI)+VD(GI)*NY(JLOC,GI)&
     &            +WD(GI)*NZ(JLOC,GI))*DETWEI(GI) 

                R1U= UUD(GI)*NX(ILOC,GI)+UVD(GI)*NY(ILOC,GI) &
     &           + UWD(GI)*NZ(ILOC,GI)   
                R2U=(UUD(GI)*NX(JLOC,GI)+UVD(GI)*NY(JLOC,GI)&
     &            +UWD(GI)*NZ(JLOC,GI))*DETWEI(GI) 

                R1V= VUD(GI)*NX(ILOC,GI)+VVD(GI)*NY(ILOC,GI) &
     &           + VWD(GI)*NZ(ILOC,GI)   
                R2V=(VUD(GI)*NX(JLOC,GI)+VVD(GI)*NY(JLOC,GI)&
     &            +VWD(GI)*NZ(JLOC,GI))*DETWEI(GI) 

                R1W= WUD(GI)*NX(ILOC,GI)+WVD(GI)*NY(ILOC,GI) &
     &           + WWD(GI)*NZ(ILOC,GI)   
                R2W=(WUD(GI)*NX(JLOC,GI)+WVD(GI)*NY(JLOC,GI)&
     &            +WWD(GI)*NZ(JLOC,GI))*DETWEI(GI) 
                VLND=VLND+ADV1*DENGI(GI)*N(ILOC,GI)*(R2&
     &               + BETA1*BECTWE(GI)*N(JLOC,GI)  )
                VLND=VLND+ADV2*DENGI(GI)*(-R1*N(JLOC,GI)*DETWEI(GI)&
     &               - N(ILOC,GI)*(1.-BETA2)*BECTWE(GI)*N(JLOC,GI)  )

                RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
!  N.B.   BECTWE(GI) = CTY(GI)*DETWEI(GI).
                UU=GAMMA(GI) *DENGI(GI)*R1 *R2

                HL1=HL1 + UU + GAMMAU(GI)*DENGI(GI)*R1U*R2U
                HL2=HL2 + UU + GAMMAV(GI)*DENGI(GI)*R1V*R2V
                HL3=HL3 + UU + GAMMAW(GI)*DENGI(GI)*R1W*R2W

                VABSX=VABSX + XABDGI(GI)*RNN
                VABSY=VABSY + YABDGI(GI)*RNN
                VABSZ=VABSZ + ZABDGI(GI)*RNN

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     calculate extra nonlinear term for LANS-alpha - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               
                if(using_lans) then
                   if(nonlinearity_option.eq.2) then
                      assert(iglv/=-1) ! Assert not uninitialised
                      udx(gi) = udx(gi) + NX(iloc,gi)*smooth_Nu(iglv)
                      udy(gi) = udy(gi) + Ny(iloc,gi)*smooth_NU(iglv)
                      udz(gi) = udz(gi) + Nz(iloc,gi)*smooth_NU(iglv)
                      vdx(gi) = vdx(gi) + NX(iloc,gi)*smooth_Nv(iglv)
                      vdy(gi) = vdy(gi) + Ny(iloc,gi)*smooth_Nv(iglv)
                      vdz(gi) = vdz(gi) + Nz(iloc,gi)*smooth_Nv(iglv)
                      wdx(gi) = wdx(gi) + NX(iloc,gi)*smooth_Nw(iglv)
                      wdy(gi) = wdy(gi) + Ny(iloc,gi)*smooth_Nw(iglv)
                      wdz(gi) = wdz(gi) + Nz(iloc,gi)*smooth_Nw(iglv)
                  
                      LANS_nonlin_XX = LANS_nonlin_XX + RNN*UDx(gi)
                      LANS_nonlin_XY = LANS_nonlin_XY + RNN*vDx(gi)
                      LANS_nonlin_XZ = LANS_nonlin_XZ + RNN*wDx(gi)
                      LANS_nonlin_YX = LANS_nonlin_YX + RNN*uDy(gi)
                      LANS_nonlin_YY = LANS_nonlin_YY + RNN*vDy(gi)
                      LANS_nonlin_YZ = LANS_nonlin_YZ + RNN*wDy(gi)
                      LANS_nonlin_ZX = LANS_nonlin_ZX + RNN*uDz(gi)
                      LANS_nonlin_ZY = LANS_nonlin_ZY + RNN*vDz(gi)
                      LANS_nonlin_ZZ = LANS_nonlin_ZZ + RNN*wDZ(gi)
                   end if
                endif
               
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! END: calculate extra nolinear term for LANS-alpha - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

                VLM =VLM +RNN
                VLMD=VLMD+RNN*DENGI(GI)
                VLMOME=VLMOME+RNN*DENGI(GI)&
     &               *FUNOME(XD(GI),YD(GI),ZD(GI))

                VLKXX=VLKXX+NX(ILOC,GI)*NX(JLOC,GI)*MXXDTW(GI)
                VLKXY=VLKXY+NX(ILOC,GI)*NY(JLOC,GI)*MXYDTW(GI)
                VLKXZ=VLKXZ+NX(ILOC,GI)*NZ(JLOC,GI)*MXZDTW(GI)
               
                VLKYX=VLKYX+NY(ILOC,GI)*NX(JLOC,GI)*MXYDTW(GI)
                VLKYY=VLKYY+NY(ILOC,GI)*NY(JLOC,GI)*MYYDTW(GI)
                VLKYZ=VLKYZ+NY(ILOC,GI)*NZ(JLOC,GI)*MYZDTW(GI)
               
                VLKZX=VLKZX+NZ(ILOC,GI)*NX(JLOC,GI)*MXZDTW(GI)
                VLKZY=VLKZY+NZ(ILOC,GI)*NY(JLOC,GI)*MYZDTW(GI)
                VLKZZ=VLKZZ+NZ(ILOC,GI)*NZ(JLOC,GI)*MZZDTW(GI)
! cty terms...
                VXVXY=VXVXY+NX(ILOC,GI)*NY(JLOC,GI)*MXXDTW(GI)
                VXVXZ=VXVXZ+NX(ILOC,GI)*NZ(JLOC,GI)*MXXDTW(GI)

                VYVYX=VYVYX+NY(ILOC,GI)*NX(JLOC,GI)*MYYDTW(GI)
                VYVYZ=VYVYZ+NY(ILOC,GI)*NZ(JLOC,GI)*MYYDTW(GI)

                VZVZX=VZVZX+NZ(ILOC,GI)*NX(JLOC,GI)*MZZDTW(GI)
                VZVZY=VZVZY+NZ(ILOC,GI)*NY(JLOC,GI)*MZZDTW(GI)

                IF(TENSOR) THEN
                   VLKTEN=VLKTEN+ NX(ILOC,GI)*(NX(JLOC,GI)*MXXDTW(GI)&
     &                     +NY(JLOC,GI)*MXYDTW(GI)+NZ(JLOC,GI)*MXZDTW(GI))&
     &                     +NY(ILOC,GI)*(NX(JLOC,GI)*MXYDTW(GI)&
     &                     +NY(JLOC,GI)*MYYDTW(GI)+NZ(JLOC,GI)*MYZDTW(GI))&
     &                     +NZ(ILOC,GI)*(NX(JLOC,GI)*MXZDTW(GI)&
     &                     +NY(JLOC,GI)*MYZDTW(GI)+NZ(JLOC,GI)*MZZDTW(GI))
                   IF(VERSIO.LE.20) VLKTEN=0.0
                   IF(TEN4TH) THEN
! subtract out extended stencil 2nd order operator to get a 4th order operator...
                      VECX(GLOBI)=VECX(GLOBI) + NX(ILOC,GI)*N(JLOC,GI)*(&
     &  MXXDTW(GI)*VEDUDX(GLOBJ) +MXYDTW(GI)*VEDUDY(GLOBJ) +MXZDTW(GI)*VEDUDZ(GLOBJ) )  &
     &                     + NY(ILOC,GI)*N(JLOC,GI)*(&
     &                  MXYDTW(GI)*VEDUDX(GLOBJ) +MYYDTW(GI)*VEDUDY(GLOBJ) +MYZDTW(GI)*VEDUDZ(GLOBJ) ) &
     &                     + NZ(ILOC,GI)*N(JLOC,GI)*( &
     &                  MXZDTW(GI)*VEDUDX(GLOBJ) +MYZDTW(GI)*VEDUDY(GLOBJ) +MZZDTW(GI)*VEDUDZ(GLOBJ) ) 
     
                      VECY(GLOBI)=VECY(GLOBI) + NX(ILOC,GI)*N(JLOC,GI)*( &
     &                  MXXDTW(GI)*VEDVDX(GLOBJ) +MXYDTW(GI)*VEDVDY(GLOBJ) +MXZDTW(GI)*VEDVDZ(GLOBJ) ) &
     &                     + NY(ILOC,GI)*N(JLOC,GI)*( &
     &                  MXYDTW(GI)*VEDVDX(GLOBJ) +MYYDTW(GI)*VEDVDY(GLOBJ) +MYZDTW(GI)*VEDVDZ(GLOBJ) ) &
     &                     + NZ(ILOC,GI)*N(JLOC,GI)*( &
     &                  MXZDTW(GI)*VEDVDX(GLOBJ) +MYZDTW(GI)*VEDVDY(GLOBJ) +MZZDTW(GI)*VEDVDZ(GLOBJ) ) 
     
                      VECZ(GLOBI)=VECZ(GLOBI) + NX(ILOC,GI)*N(JLOC,GI)*( &
     &                  MXXDTW(GI)*VEDWDX(GLOBJ) +MXYDTW(GI)*VEDWDY(GLOBJ) +MXZDTW(GI)*VEDWDZ(GLOBJ) ) &
     &                     + NY(ILOC,GI)*N(JLOC,GI)*( &
     &                  MXYDTW(GI)*VEDWDX(GLOBJ) +MYYDTW(GI)*VEDWDY(GLOBJ) +MYZDTW(GI)*VEDWDZ(GLOBJ) ) &
     &                     + NZ(ILOC,GI)*N(JLOC,GI)*( &
     &                  MXZDTW(GI)*VEDWDX(GLOBJ) +MYZDTW(GI)*VEDWDY(GLOBJ) +MZZDTW(GI)*VEDWDZ(GLOBJ) )
                   ENDIF
                ENDIF
               
                IF(LESTWO.AND.(.NOT.((DISOPT.EQ.47).OR.(DISOPT.EQ.48)))) THEN
                   do II=1,6
                      do JJ=1,6
                         IF(LESNOD.EQ.1) THEN
! For proper 4th order dissipation.
                            EE(II,JJ)=0.0
                            do KLOC=1,NLOC
                               KNOD=NDGLNO((ELE-1)*NLOC+KLOC)
                               EE(II,JJ)=EE(II,JJ) + N(KLOC,GI)*EV(II,JJ,KNOD)*DENPT(KNOD)
                            END DO
                            EEN(II,JJ)=EV(II,JJ,GLOBJ)
                         ELSE
                            EE(II,JJ)=E(II,JJ,GI)*DENGI(GI)
                            EEN(II,JJ)=E(II,JJ,GI)*DENGI(GI)
                         ENDIF
                      END DO
                   END DO

! First operator (on u for x-momentum equation)
                   OP(1,1) = OP(1,1) +&
     &           DETWEI(GI)*( &
     &           EE(1,1)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(4,4)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(5,5)*NZ(ILOC,GI)*NZ(JLOC,GI) + &
     &          (EE(1,4)+EE(4,1))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(1,5)+EE(5,1))*NX(ILOC,GI)*NZ(JLOC,GI) + &
     &          (EE(4,5)+EE(5,4))*NY(ILOC,GI)*NZ(JLOC,GI))

! Second operator (on v for x-momentum equation)
                   OP(1,2) = OP(1,2) +&
     &           DETWEI(GI)*(&
     &           EE(1,4)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(4,2)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(5,6)*NZ(ILOC,GI)*NZ(JLOC,GI) + &
     &          (EE(1,2)+EE(4,4))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(1,6)+EE(5,4))*NX(ILOC,GI)*NZ(JLOC,GI) + &
     &          (EE(4,6)+EE(5,2))*NY(ILOC,GI)*NZ(JLOC,GI))

! Third operator (on w for x-momentum equation)
                   OP(1,3) = OP(1,3) +&
     &           DETWEI(GI)*(&
     &           EE(1,5)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(4,6)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(5,3)*NZ(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(1,6)+EE(4,5))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(1,3)+EE(5,5))*NX(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(4,3)+EE(5,6))*NY(ILOC,GI)*NZ(JLOC,GI))

! Fourth operator (on u for y-momentum equation)
                   OP(2,1) = OP(2,1) +&
     &           DETWEI(GI)*(&
     &           EE(4,1)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(2,4)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(6,5)*NZ(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(2,1)+EE(4,4))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(4,5)+EE(6,1))*NX(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(2,5)+EE(6,4))*NY(ILOC,GI)*NZ(JLOC,GI))

! Fifth operator (on v for y-momentum equation)
                   OP(2,2) = OP(2,2) +&
     &           DETWEI(GI)*(&
     &           EE(4,4)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(2,2)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(6,6)*NZ(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(2,4)+EE(4,2))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(4,6)+EE(6,4))*NX(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(2,6)+EE(6,2))*NY(ILOC,GI)*NZ(JLOC,GI))

! Sixth operator (on w for y-momentum equation)
                   OP(2,3) = OP(2,3) +&
     &           DETWEI(GI)*(&
     &           EE(4,5)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(2,6)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(6,3)*NZ(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(2,5)+EE(4,6))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(4,3)+EE(6,5))*NX(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(2,3)+EE(6,6))*NY(ILOC,GI)*NZ(JLOC,GI))

! Seventh operator (on u for z-momentum equation)
                   OP(3,1) = OP(3,1) +&
     &           DETWEI(GI)*(&
     &           EE(5,1)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(6,4)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(3,5)*NZ(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(5,4)+EE(6,1))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(3,1)+EE(5,5))*NX(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(3,4)+EE(6,5))*NY(ILOC,GI)*NZ(JLOC,GI))

! Eighth operator (on v for z-momentum equation)
                   OP(3,2) = OP(3,2) +&
     &           DETWEI(GI)*(&
     &           EE(5,4)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(6,2)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(3,6)*NZ(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(5,2)+EE(6,4))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(3,4)+EE(5,6))*NX(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(3,2)+EE(6,6))*NY(ILOC,GI)*NZ(JLOC,GI))

! Ninth operator (on w for z-momentum equation)
                   OP(3,3) = OP(3,3) +&
     &           DETWEI(GI)*(&
     &           EE(5,5)*NX(ILOC,GI)*NX(JLOC,GI) +&
     &           EE(6,6)*NY(ILOC,GI)*NY(JLOC,GI) +&
     &           EE(3,3)*NZ(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(5,6)+EE(6,5))*NX(ILOC,GI)*NY(JLOC,GI) +&
     &          (EE(3,5)+EE(5,3))*NX(ILOC,GI)*NZ(JLOC,GI) +&
     &          (EE(3,6)+EE(6,3))*NY(ILOC,GI)*NZ(JLOC,GI))
     
! For the enlarged stencil viscouse operator used to find 4th order operator...
                   IF(LES4TH.EQ.1) THEN
! First operator (on u for x-momentum equation)
                      OP4TH(1,1) = OP4TH(1,1) +&
     &           DETWEI(GI)*( &
     &           EEN(1,1)*NX(ILOC,GI)*N(JLOC,GI)*VEDUDX(GLOBJ) +&
     &           EEN(4,4)*NY(ILOC,GI)*N(JLOC,GI)*VEDUDY(GLOBJ) +&
     &           EEN(5,5)*NZ(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ) + &
     &          (EEN(1,4)+EEN(4,1))*NX(ILOC,GI)*N(JLOC,GI)*VEDUDY(GLOBJ) +&
     &          (EEN(1,5)+EEN(5,1))*NX(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ) + &
     &          (EEN(4,5)+EEN(5,4))*NY(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ))

! Second operator (on v for x-momentum equation)
                      OP4TH(1,2) = OP4TH(1,2) +&
     &           DETWEI(GI)*(&
     &           EEN(1,4)*NX(ILOC,GI)*N(JLOC,GI)*VEDVDX(GLOBJ) +&
     &           EEN(4,2)*NY(ILOC,GI)*N(JLOC,GI)*VEDVDY(GLOBJ) +&
     &           EEN(5,6)*NZ(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ) + &
     &          (EEN(1,2)+EEN(4,4))*NX(ILOC,GI)*N(JLOC,GI)*VEDVDY(GLOBJ) +&
     &          (EEN(1,6)+EEN(5,4))*NX(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ) + &
     &          (EEN(4,6)+EEN(5,2))*NY(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ))

! Third operator (on w for x-momentum equation)
                      OP4TH(1,3) = OP4TH(1,3) +&
     &           DETWEI(GI)*(&
     &           EEN(1,5)*NX(ILOC,GI)*N(JLOC,GI)*VEDWDX(GLOBJ) +&
     &           EEN(4,6)*NY(ILOC,GI)*N(JLOC,GI)*VEDWDY(GLOBJ) +&
     &           EEN(5,3)*NZ(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ) +&
     &          (EEN(1,6)+EEN(4,5))*NX(ILOC,GI)*N(JLOC,GI)*VEDWDY(GLOBJ) +&
     &          (EEN(1,3)+EEN(5,5))*NX(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ) +&
     &          (EEN(4,3)+EEN(5,6))*NY(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ))

! Fourth operator (on u for y-momentum equation)
                      OP4TH(2,1) = OP4TH(2,1) +&
     &           DETWEI(GI)*(&
     &           EEN(4,1)*NX(ILOC,GI)*N(JLOC,GI)*VEDUDX(GLOBJ) +&
     &           EEN(2,4)*NY(ILOC,GI)*N(JLOC,GI)*VEDUDY(GLOBJ) +&
     &           EEN(6,5)*NZ(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ) +&
     &          (EEN(2,1)+EEN(4,4))*NX(ILOC,GI)*N(JLOC,GI)*VEDUDY(GLOBJ) +&
     &          (EEN(4,5)+EEN(6,1))*NX(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ) +&
     &          (EEN(2,5)+EEN(6,4))*NY(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ))

! Fifth operator (on v for y-momentum equation)
                      OP4TH(2,2) = OP4TH(2,2) +&
     &           DETWEI(GI)*(&
     &           EEN(4,4)*NX(ILOC,GI)*N(JLOC,GI)*VEDVDX(GLOBJ) +&
     &           EEN(2,2)*NY(ILOC,GI)*N(JLOC,GI)*VEDVDY(GLOBJ) +&
     &           EEN(6,6)*NZ(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ) +&
     &          (EEN(2,4)+EEN(4,2))*NX(ILOC,GI)*N(JLOC,GI)*VEDVDY(GLOBJ) +&
     &          (EEN(4,6)+EEN(6,4))*NX(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ) +&
     &          (EEN(2,6)+EEN(6,2))*NY(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ))

! Sixth operator (on w for y-momentum equation)
                      OP4TH(2,3) = OP4TH(2,3) +&
     &           DETWEI(GI)*(&
     &           EEN(4,5)*NX(ILOC,GI)*N(JLOC,GI)*VEDWDX(GLOBJ) +&
     &           EEN(2,6)*NY(ILOC,GI)*N(JLOC,GI)*VEDWDY(GLOBJ) +&
     &           EEN(6,3)*NZ(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ) +&
     &          (EEN(2,5)+EEN(4,6))*NX(ILOC,GI)*N(JLOC,GI)*VEDWDY(GLOBJ) +&
     &          (EEN(4,3)+EEN(6,5))*NX(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ) +&
     &          (EEN(2,3)+EEN(6,6))*NY(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ))

! Seventh operator (on u for z-momentum equation)
                      OP4TH(3,1) = OP4TH(3,1) +&
     &           DETWEI(GI)*(&
     &           EEN(5,1)*NX(ILOC,GI)*N(JLOC,GI)*VEDUDX(GLOBJ) +&
     &           EEN(6,4)*NY(ILOC,GI)*N(JLOC,GI)*VEDUDY(GLOBJ) +&
     &           EEN(3,5)*NZ(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ) +&
     &          (EEN(5,4)+EEN(6,1))*NX(ILOC,GI)*N(JLOC,GI)*VEDUDY(GLOBJ) +&
     &          (EEN(3,1)+EEN(5,5))*NX(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ) +&
     &          (EEN(3,4)+EEN(6,5))*NY(ILOC,GI)*N(JLOC,GI)*VEDUDZ(GLOBJ))

! Eighth operator (on v for z-momentum equation)
                      OP4TH(3,2) = OP4TH(3,2) +&
     &           DETWEI(GI)*(&
     &           EEN(5,4)*NX(ILOC,GI)*N(JLOC,GI)*VEDVDX(GLOBJ) +&
     &           EEN(6,2)*NY(ILOC,GI)*N(JLOC,GI)*VEDVDY(GLOBJ) +&
     &           EEN(3,6)*NZ(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ) +&
     &          (EEN(5,2)+EEN(6,4))*NX(ILOC,GI)*N(JLOC,GI)*VEDVDY(GLOBJ) +&
     &          (EEN(3,4)+EEN(5,6))*NX(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ) +&
     &          (EEN(3,2)+EEN(6,6))*NY(ILOC,GI)*N(JLOC,GI)*VEDVDZ(GLOBJ))

! Ninth operator (on w for z-momentum equation)
                      OP4TH(3,3) = OP4TH(3,3) +&
     &           DETWEI(GI)*(&
     &           EEN(5,5)*NX(ILOC,GI)*N(JLOC,GI)*VEDWDX(GLOBJ) +&
     &           EEN(6,6)*NY(ILOC,GI)*N(JLOC,GI)*VEDWDY(GLOBJ) +&
     &           EEN(3,3)*NZ(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ) +&
     &          (EEN(5,6)+EEN(6,5))*NX(ILOC,GI)*N(JLOC,GI)*VEDWDY(GLOBJ) +&
     &          (EEN(3,5)+EEN(5,3))*NX(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ) +&
     &          (EEN(3,6)+EEN(6,3))*NY(ILOC,GI)*N(JLOC,GI)*VEDWDZ(GLOBJ))
! endof IF(LES4TH.EQ.1) THEN...
                   ENDIF
                ENDIF
             end do

! The 3 terms below are not calculated because we know them...
             VXVXX=VLKXX
             VYVYY=VLKYY
             VZVZZ=VLKZZ

! Make the lumped absorption consistent with a lumped source. 
             VABSX=VABSX*(1.-ABLUMP) &
     &           + VLM*XABSOR(GLOBIS)*ABLUMP
             VABSY=VABSY*(1.-ABLUMP) &
     &           + VLM*YABSOR(GLOBIS)*ABLUMP
             VABSZ=VABSZ*(1.-ABLUMP) &
     &           + VLM*ZABSOR(GLOBIS)*ABLUMP

! lump the mass matrix.(LUMP=1. if we are to LUMP the mass matrix)
             ICENT=CENTRM(GLOBI)
             BIGM(ICENT+IBL11)=BIGM(ICENT+IBL11)+VLMD*RLUMP&
     &                        + CVMVEC(GLOBIS)*VLM &
     &                        +DT*ABTHETA*VABSX*ABLUMP
             BIGM(ICENT+IBL22)=BIGM(ICENT+IBL22)+VLMD*RLUMP&
     &                        + CVMVEC(GLOBIS)*VLM &
     &                        +DT*ABTHETA*VABSY*ABLUMP
             BIGM(ICENT+IBL33)=BIGM(ICENT+IBL33)+VLMD*RLUMP&
     &                        + CVMVEC(GLOBIS)*VLM &
     &                        +DT*ABTHETA*VABSZ*ABLUMP
! Presumably we lump the absorption to put in to pressure.
             ML(GLOBI)=ML(GLOBI) + VLMD + CVMVEC(GLOBIS)*VLM &
     &                + DT*ABTHETA*VABSX*ABLUMP

             R=VLND
             IF(LESTWO) THEN
                VECX(GLOBI)=VECX(GLOBI) -(R+HL1+VABSX*(1.-ABLUMP)+OP(1,1))&
     &            *U(GLOBJ) -OP(1,2)*V(GLOBJ) -OP(1,3)*W(GLOBJ) +2.*VLMOME&
     &            *V(GLOBJ)*(1.-COLUMP)*ROTRHS +2. *VLMOME*V(GLOBI)*COLUMP&
     &            *ROTRHS +vlm*sourcx(globjS)  *(1. -RSLUMP)
                VECY(GLOBI)=VECY(GLOBI)&
     &            -(R+HL2+VABSY*(1.-ABLUMP)+OP(2,2))*V(GLOBJ)&
     &            -OP(2,1)*U(GLOBJ)&
     &            -OP(2,3)*W(GLOBJ)&
     &            -2.*VLMOME*U(GLOBJ)*(1.-COLUMP)*ROTRHS&
     &            -2.*VLMOME*U(GLOBI)*COLUMP*ROTRHS&
     &            +vlm*sourcy(globjS)  *(1.-RSLUMP)
! We are using the transpose.
                VECZ(GLOBI)=VECZ(GLOBI)&
     &            -(R+HL3+VABSZ*(1.-ABLUMP)+OP(3,3))*W(GLOBJ)&
     &            -OP(3,1)*U(GLOBJ)&
     &            -OP(3,2)*V(GLOBJ)&
     &            +vlm*sourcz(globjS)  *(1.-RSLUMP)
                IF(LES4TH.EQ.1) THEN
! nb change the sign so we subtract out the enlarged stencil...
                   VECX(GLOBI)=VECX(GLOBI) +OP4TH(1,1) +OP4TH(1,2) +OP4TH(1,3)
                   VECY(GLOBI)=VECY(GLOBI) +OP4TH(2,1) +OP4TH(2,2) +OP4TH(2,3)
                   VECZ(GLOBI)=VECZ(GLOBI) +OP4TH(3,1) +OP4TH(3,2) +OP4TH(3,3)
                ENDIF
! put in the viscous terms...
                IF(.NOT.TENSOR) THEN
                   VECX(GLOBI)=VECX(GLOBI) &
     &               -(2.*VLKXX+VLKYY+VLKZZ)*U(GLOBJ)&
     &               -VLKYX*V(GLOBJ)&
     &               -VLKZX*W(GLOBJ)&
     &               +TWOTHI*(VXVXX*U(GLOBJ)+VXVXY*V(GLOBJ)+VXVXZ*W(GLOBJ))

                   VECY(GLOBI)=VECY(GLOBI)&
     &               -(VLKXX+VLKZZ+2.*VLKYY)*V(GLOBJ)&
     &               -VLKXY*U(GLOBJ)&
     &               -VLKZY*W(GLOBJ)&
     &               +TWOTHI*(VYVYX*U(GLOBJ)+VYVYY*V(GLOBJ)+VYVYZ*W(GLOBJ))

                   VECZ(GLOBI)=VECZ(GLOBI)&
     &               -(VLKXX+VLKYY+2.*VLKZZ)*W(GLOBJ)&
     &               -VLKXZ*U(GLOBJ)&
     &               -VLKYZ*V(GLOBJ)&
     &               +TWOTHI*(VZVZX*U(GLOBJ)+VZVZY*V(GLOBJ)+VZVZZ*W(GLOBJ))
                ELSE
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXICY OF STRESS FORM...
                   VECX(GLOBI)=VECX(GLOBI) - VLKTEN*U(GLOBJ)
                   VECY(GLOBI)=VECY(GLOBI) - VLKTEN*V(GLOBJ)
                   VECZ(GLOBI)=VECZ(GLOBI) - VLKTEN*W(GLOBJ)
                ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     formulate source terms in LANS-alpha form - cjc 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

             ELSEif(using_lans) then
                if(nonlinearity_option<3) then
                   VECX(GLOBI)=VECX(GLOBI) -(R+HL1+VABSX*(1.-ABLUMP)+2.*VLKXX&
     &                 +VLKYY+VLKZZ)*U(GLOBJ) -VLKYX*V(GLOBJ) -VLKZX*W(GLOBJ)&
     &                 +TWOTHI*(VXVXX*U(GLOBJ)+VXVXY*V(GLOBJ)+VXVXZ*W(GLOBJ))+&
     &                 lans_nonlin_xx*U(globj)+ lans_nonlin_yx*V(globj) +&
     &                 lans_nonlin_zx*W(globj)&
     &                 +2.*VLMOME*smooth_nv(globj) &
     &                 *(1.-COLUMP)*ROTRHS&
     &                 +2.*VLMOME*smooth_nv(globi)&
     &                 *COLUMP*ROTRHS&
     &                 +vlm *sourcx(globjS)  *(1. -RSLUMP)
                  
                   VECY(GLOBI)=VECY(GLOBI) -(R+HL2+VABSY*(1.-ABLUMP)+VLKXX+VLKZZ&
     &                 +2.*VLKYY)*V(GLOBJ) -VLKXY*U(GLOBJ) -VLKZY*W(GLOBJ)&
     &                 +TWOTHI*(VYVYX*U(GLOBJ)+VYVYY*V(GLOBJ)+VYVYZ*W(GLOBJ))&
     &                 -2.*VLMOME*smooth_nu(globj) &
     &                 *(1.-COLUMP)*ROTRHS&
     &                 -2.*VLMOME*smooth_nu(globi)&
     &                 *COLUMP*ROTRHS&
     &                 +lans_nonlin_xy*U(globj)+ lans_nonlin_yy*V(globj) +&
     &                 lans_nonlin_zy*W(globj) +vlm *sourcy(globjS)  *(1.-RSLUMP)
                   VECZ(GLOBI)=VECZ(GLOBI) -(R+HL3+VABSZ*(1.-ABLUMP)+VLKXX+VLKYY&
     &                 +2.*VLKZZ)*W(GLOBJ) -VLKXZ*U(GLOBJ) -VLKYZ*V(GLOBJ)&
     &                 +TWOTHI*(VZVZX*U(GLOBJ)+VZVZY*V(GLOBJ)+VZVZZ*W(GLOBJ))&
     &                 +lans_nonlin_xz*U(globj)+ lans_nonlin_yz*V(globj) +&
     &                 lans_nonlin_zz*W(globj) +vlm *sourcz(globjS)  *(1.&
     &                 -RSLUMP)
                else
                   VECX(GLOBI)=VECX(GLOBI) -(R+HL1+VABSX*(1.-ABLUMP)+2.&
     &                 *VLKXX+VLKYY+VLKZZ)*U(GLOBJ) -VLKYX*V(GLOBJ)&
     &                 -VLKZX*W(GLOBJ)+TWOTHI*(VXVXX*U(GLOBJ)+VXVXY&
     &                 *V(GLOBJ)+VXVXZ*W(GLOBJ))+vlm *sourcx(globjS)&
     &                 *(1. -RSLUMP)
                  
                   VECY(GLOBI)=VECY(GLOBI) -(R+HL2+VABSY*(1.-ABLUMP)&
     &                 +VLKXX+VLKZZ+2.*VLKYY)*V(GLOBJ) -VLKXY*U(GLOBJ)&
     &                 -VLKZY*W(GLOBJ)+TWOTHI*(VYVYX*U(GLOBJ)+VYVYY&
     &                 *V(GLOBJ)+VYVYZ*W(GLOBJ))+vlm *sourcy(globjS)&
     &                 *(1.-RSLUMP)
                  
                   VECZ(GLOBI)=VECZ(GLOBI) -(R+HL3+VABSZ*(1.-ABLUMP)&
     &                 +VLKXX+VLKYY+2.*VLKZZ)*W(GLOBJ) -VLKXZ*U(GLOBJ)&
     &                 -VLKYZ*V(GLOBJ)+TWOTHI*(VZVZX*U(GLOBJ)+VZVZY&
     &                 *V(GLOBJ)+VZVZZ*W(GLOBJ))+vlm *sourcz(globjS)&
     &                 *(1.-RSLUMP)
                end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! END: formulate source terms in LANS-alpha form - cjc 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

             else     
                VECX(GLOBI)=VECX(GLOBI)&
     &            -(R+HL1+VABSX*(1.-ABLUMP)+2.*VLKXX+VLKYY+VLKZZ)*U(GLOBJ)&
     &            -VLKYX*V(GLOBJ)&
     &            -VLKZX*W(GLOBJ)&
     &            +TWOTHI*(VXVXX*U(GLOBJ)+VXVXY*V(GLOBJ)+VXVXZ*W(GLOBJ))&
     &            +2.*VLMOME*V(GLOBJ)*(1.-COLUMP)*ROTRHS&
     &            +2.*VLMOME*V(GLOBI)*COLUMP*ROTRHS&
     &            +vlm*sourcx(globjS)  *(1.-RSLUMP)

                VECY(GLOBI)=VECY(GLOBI)&
     &            -(R+HL2+VABSY*(1.-ABLUMP)+VLKXX+VLKZZ+2.*VLKYY)*V(GLOBJ)&
     &            -VLKXY*U(GLOBJ)&
     &            -VLKZY*W(GLOBJ)&
     &            +TWOTHI*(VYVYX*U(GLOBJ)+VYVYY*V(GLOBJ)+VYVYZ*W(GLOBJ))&
     &            -2.*VLMOME*U(GLOBJ)*(1.-COLUMP)*ROTRHS&
     &            -2.*VLMOME*U(GLOBI)*COLUMP*ROTRHS&
     &            +vlm*sourcy(globjS)  *(1.-RSLUMP)
! We are using the transpose. 

                VECZ(GLOBI)=VECZ(GLOBI)&
     &            -(R+HL3+VABSZ*(1.-ABLUMP)+VLKXX+VLKYY+2.*VLKZZ)*W(GLOBJ)&
     &            -VLKXZ*U(GLOBJ)&
     &            -VLKYZ*V(GLOBJ)&
     &            +TWOTHI*(VZVZX*U(GLOBJ)+VZVZY*V(GLOBJ)+VZVZZ*W(GLOBJ))&
     &            +vlm*sourcz(globjS)  *(1.-RSLUMP)
             ENDIF

! Lumped source & absorption.
             VECX(GLOBI)=VECX(GLOBI)-VABSX*ABLUMP*U(GLOBI)&
     &                            +VLM*sourcx(GLOBIS) *RSLUMP
             VECY(GLOBI)=VECY(GLOBI)-VABSY*ABLUMP*V(GLOBI)&
     &                            +VLM*sourcy(GLOBIS) *RSLUMP
             VECZ(GLOBI)=VECZ(GLOBI)-VABSZ*ABLUMP*W(GLOBI)&
     &                            +VLM*sourcz(GLOBIS) *RSLUMP

! we will now place contributions into the matrices BIGM and BIGT
! these are symm matrices so we only store the upper diagonal
! parts of these matrices.
             MC1=(VLND*RSYM)*THETA*DT + VLMD*(1.-RLUMP)

             LOWER=FINDRM(GLOBI) 
             UPPER=FINDRM(GLOBI+1)-1
7000         CONTINUE

             INUM=LOWER+(UPPER-LOWER+1)/2
             IF(GLOBJ.GE.COLM(INUM) )  THEN
                LOWER=INUM
             ELSE
                UPPER=INUM
             ENDIF
             IF(UPPER-LOWER.LE.1) THEN
                IF(GLOBJ.EQ.COLM(LOWER)) THEN
                   COUNT=LOWER
                ELSE
                   COUNT=UPPER
                ENDIF
                GOTO 9000
             ENDIF
             GOTO 7000
9000         CONTINUE

! Before embarking on this bit, call subroutine SFS2 to calculate 
! differential operators.

             IF(LESTWO) THEN
                BIGM(COUNT+IBL11)=BIGM(COUNT+IBL11)&
     &            +MC1 +DT*THETA*(HL1&
     &            +OP(1,1))  + DT*ABTHETA*VABSX*(1.-ABLUMP)
                BIGM(COUNT+IBL22)=BIGM(COUNT+IBL22)&
     &            +MC1 +DT*THETA*(HL2&
     &            +OP(2,2))  +DT*ABTHETA*VABSY*(1.-ABLUMP)
                BIGM(COUNT+IBL33)=BIGM(COUNT+IBL33)&
     &            +MC1 +DT*THETA*(HL3&
     &            +OP(3,3))   +DT*ABTHETA*VABSZ*(1.-ABLUMP)

                BIGM(COUNT+IBL21)=BIGM(COUNT+IBL21)&
     &            +DT*THETA*(OP(2,1)  + 2.*VLMOME*(1.-COLUMP))
                BIGM(ICENT+IBL21)=BIGM(ICENT+IBL21)&
     &            +DT*THETA*2.*VLMOME*COLUMP

                BIGM(COUNT+IBL31)=BIGM(COUNT+IBL31)&
     &            +DT*THETA*OP(3,1)
                BIGM(COUNT+IBL32)=BIGM(COUNT+IBL32)&
     &            +DT*THETA*OP(3,2) 

                BIGM(COUNT+IBL12)=BIGM(COUNT+IBL12)&
     &            +DT*THETA*(OP(1,2)  - 2.*VLMOME*(1.-COLUMP))
                BIGM(ICENT+IBL12)=BIGM(ICENT+IBL12)&
     &            -DT*THETA*2.*VLMOME*COLUMP
                BIGM(COUNT+IBL13)=BIGM(COUNT+IBL13)&
     &            +DT*THETA*OP(1,3)

                BIGM(COUNT+IBL23)=BIGM(COUNT+IBL23)&
     &            +DT*THETA*OP(2,3)

! Put in the viscous terms...
                IF(.NOT.TENSOR) THEN
                   BIGM(COUNT+IBL11)=BIGM(COUNT+IBL11)  +DT*THETA*(&
     &               +2.*VLKXX+VLKYY+VLKZZ - TWOTHI&
     &               *VXVXX)  
                   BIGM(COUNT+IBL22)=BIGM(COUNT+IBL22)  +DT*THETA*(&
     &               +VLKXX+2.*VLKYY+VLKZZ - TWOTHI&
     &               *VYVYY)
                   BIGM(COUNT+IBL33)=BIGM(COUNT+IBL33)&
     &               +DT*THETA*(&
     &               +VLKXX+VLKYY+2.*VLKZZ - TWOTHI*VZVZZ) 

                   BIGM(COUNT+IBL21)=BIGM(COUNT+IBL21) +DT*THETA*(VLKXY -&
     &               TWOTHI*VYVYX )

                   BIGM(COUNT+IBL31)=BIGM(COUNT+IBL31)&
     &               +DT*THETA*(VLKXZ - TWOTHI*VZVZX)
                   BIGM(COUNT+IBL32)=BIGM(COUNT+IBL32)&
     &               +DT*THETA*(VLKYZ - TWOTHI*VZVZY)

                   BIGM(COUNT+IBL12)=BIGM(COUNT+IBL12) +DT*THETA*(VLKYX -&
     &               TWOTHI*VXVXY )
                   BIGM(COUNT+IBL13)=BIGM(COUNT+IBL13)&
     &               +DT*THETA*(VLKZX - TWOTHI*VXVXZ)

                   BIGM(COUNT+IBL23)=BIGM(COUNT+IBL23)&
     &               +DT*THETA*(VLKZY - TWOTHI*VYVYZ) 
                ELSE
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXICY OF STRESS FORM.
                   BIGM(COUNT+IBL11)=BIGM(COUNT+IBL11)  +DT*THETA*VLKTEN
                   BIGM(COUNT+IBL22)=BIGM(COUNT+IBL22)  +DT*THETA*VLKTEN
                   BIGM(COUNT+IBL33)=BIGM(COUNT+IBL33)  +DT*THETA*VLKTEN
                ENDIF
               
             ELSEif(using_lans) then
           
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     set up the matrix for LANS-alpha model - cjc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                if(blksym) then
                   ewrite(-1,*)  'Cannot use BLKSYM=true'
                   ewrite(-1,*)  'with LANS-alpha'
                   ewrite(-1,*)  'as the matrix is not symmetric,'
                   ewrite(-1,*)  'so stopping. - cjc'
                   FLAbort('cannot use BLKSYM=true with LANS-alpha')
                endif

                BIGM(COUNT+IBL11)=BIGM(COUNT+IBL11) +MC1 +DT*THETA*(HL1&
     &            +2.*VLKXX+VLKYY+VLKZZ - TWOTHI&
     &            *VXVXX +lans_nonlin_xx)  +DT*ABTHETA*VABSX*(1.-ABLUMP)

                BIGM(COUNT+IBL22)=BIGM(COUNT+IBL22) +MC1 +DT*THETA*(HL2&
     &            +VLKXX+2.*VLKYY+VLKZZ - TWOTHI&
     &            *VYVYY +lans_nonlin_yy) + DT*ABTHETA*VABSY*(1.-ABLUMP)

                BIGM(COUNT+IBL33)=BIGM(COUNT+IBL33) +MC1 +DT*THETA*(HL3&
     &            +VLKXX+VLKYY+2.*VLKZZ - TWOTHI&
     &            *VZVZZ +lans_nonlin_zz)  + DT*ABTHETA*VABSZ*(1.-ABLUMP) 

                BIGM(COUNT+IBL21)=BIGM(COUNT+IBL21) +DT*THETA*(VLKXY -&
     &            TWOTHI*VYVYX + &
     &            lans_nonlin_xy)

                BIGM(COUNT+IBL31)=BIGM(COUNT+IBL31) +DT*THETA*(VLKXZ -&
     &            TWOTHI*VZVZX + lans_nonlin_xz)
              
                BIGM(COUNT+IBL32)=BIGM(COUNT+IBL32) +DT*THETA*(VLKYZ -&
     &            TWOTHI*VZVZY+ lans_nonlin_yz) 

                BIGM(COUNT+IBL12)=BIGM(COUNT+IBL12) +DT*THETA*(VLKYX -&
     &            TWOTHI*VXVXY  +&
     &            lans_nonlin_yx)

                BIGM(COUNT+IBL13)=BIGM(COUNT+IBL13) +DT*THETA*(VLKZX -&
     &            TWOTHI*VXVXZ + lans_nonlin_zx)

                BIGM(COUNT+IBL23)=BIGM(COUNT+IBL23) +DT*THETA*(VLKZY -&
     &            TWOTHI*VYVYZ + lans_nonlin_zy) 
             else
                BIGM(COUNT+IBL11)=BIGM(COUNT+IBL11) +MC1 +DT*THETA*(HL1&
     &            +2.*VLKXX+VLKYY+VLKZZ - TWOTHI&
     &            *VXVXX)  + DT*ABTHETA*VABSX*(1.-ABLUMP)
                BIGM(COUNT+IBL22)=BIGM(COUNT+IBL22) +MC1 +DT*THETA*(HL2&
     &            +VLKXX+2.*VLKYY+VLKZZ - TWOTHI&
     &            *VYVYY)  +DT*ABTHETA*VABSY*(1.-ABLUMP)
                BIGM(COUNT+IBL33)=BIGM(COUNT+IBL33)&
     &            +MC1 +DT*THETA*(HL3&
     &            +VLKXX+VLKYY+2.*VLKZZ - TWOTHI*VZVZZ) +DT*ABTHETA*VABSZ*(1.-ABLUMP)

                BIGM(COUNT+IBL21)=BIGM(COUNT+IBL21) +DT*THETA*(VLKXY -&
     &            TWOTHI*VYVYX  + 2.*VLMOME*(1.-COLUMP))
                BIGM(ICENT+IBL21)=BIGM(ICENT+IBL21)&
     &            +DT*THETA*2.*VLMOME*COLUMP

                BIGM(COUNT+IBL31)=BIGM(COUNT+IBL31)&
     &            +DT*THETA*(VLKXZ - TWOTHI*VZVZX)
                BIGM(COUNT+IBL32)=BIGM(COUNT+IBL32)&
     &            +DT*THETA*(VLKYZ - TWOTHI*VZVZY)

                IF(.NOT.BLKSYM) THEN
                   BIGM(COUNT+IBL12)=BIGM(COUNT+IBL12) +DT*THETA*(VLKYX -&
     &               TWOTHI*VXVXY  - 2.*VLMOME*(1.-COLUMP))
                   BIGM(ICENT+IBL12)=BIGM(ICENT+IBL12)&
     &               -DT*THETA*2.*VLMOME*COLUMP
                   BIGM(COUNT+IBL13)=BIGM(COUNT+IBL13)&
     &               +DT*THETA*(VLKZX - TWOTHI*VXVXZ)

                   BIGM(COUNT+IBL23)=BIGM(COUNT+IBL23)&
     &               +DT*THETA*(VLKZY - TWOTHI*VYVYZ)
                ENDIF
             ENDIF

          end do
       end do

    end do

    ewrite(2,*) 'NONODS,XNONOD:',NONODS,XNONOD
    ewrite(2,*) 'abslum:',abslum
    CALL PMINMX(ML,NONODS,'******ML  ')
    CALL PMINMX(XABSOR,NONODS,'******XABSOR  ')
    CALL PMINMX(YABSOR,NONODS,'******YABSOR  ')
    CALL PMINMX(ZABSOR,NONODS,'******ZABSOR  ')
    CALL PMINMX(DENPT,NONODS,'******DENPT  ')
    CALL PMINMX(sourcx,NONODS,'******sourcx  ')
    CALL PMINMX(sourcy,NONODS,'******sourcy  ')
    CALL PMINMX(sourcz,NONODS,'******sourcz  ')
    CALL PMINMX(MUPTXX,nonods,'******MUPTXX  ')
    CALL PMINMX(MUPTXY,nonods,'******MUPTXY  ')
    CALL PMINMX(MUPTXZ,nonods,'******MUPTXZ  ')
    CALL PMINMX(MUPTYY,nonods,'******MUPTYY  ')
    CALL PMINMX(MUPTYZ,nonods,'******MUPTYZ  ')
    CALL PMINMX(MUPTZZ,nonods,'******MUPTZZ  ')
    ewrite(2,*)  "+R2NORM(vecx, nonods, 0) = ", R2NORM(vecx, nonods, 0)
    ewrite(2,*)  "+R2NORM(vecy, nonods, 0) = ", R2NORM(vecy, nonods, 0)
    ewrite(2,*)  "+R2NORM(vecz, nonods, 0) = ", R2NORM(vecz, nonods, 0)
    ewrite_minmax(vecx)
    ewrite_minmax(vecy)
    ewrite_minmax(vecz)
    
    ewrite(2,*) 'BETA,D0,f0,twothi:',BETA,D0,f0,twothi
     
    do I=1,-NONODS
       IF(ABS(x(I))+abs(Y(I)).lt.0.1) then
          ewrite(2,*) 'i,z(i),SOURCZ(i),denpt(i):',&
     &                      i,z(i),SOURCZ(i),denpt(i)
       endif
    ENDDO
    ewrite(2,*) 'horizontal:'
    do I=1,-NONODS
       IF(ABS(Z(I)-86.58).lt.0.5) then
          ewrite(2,*) 'x(i),y(i),SOURCx(i),SOURCy(i),SOURCZ(i),denpt(i):',&
     &                   x(i),y(i),SOURCx(i),SOURCy(i),SOURCZ(i),denpt(i)
       endif
    ENDDO
    RETURN

    deallocate(NC)
    deallocate(NCLX)
    deallocate(NCLY)
    deallocate(NCLZ)
    deallocate(NCX)
    deallocate(NCY)
    deallocate(NCZ)

    deallocate(VEDUDX)
    deallocate(VEDUDY)
    deallocate(VEDUDZ)
    deallocate(VEDVDX)
    deallocate(VEDVDY)
    deallocate(VEDVDZ)
    deallocate(VEDWDX)
    deallocate(VEDWDY)
    deallocate(VEDWDZ)

  END SUBROUTINE DIFF3D

  SUBROUTINE DIFF3DSIMP&
     &         (acctim,U,V,W,NU,NV,NW,UG,VG,WG,&
     &          SOURCX,SOURCY,SOURCZ,X,Y,Z, D0,&
     &          NBUOY,BOY_ML,SOXGRA,SOYGRA,SOZGRA,&
     &          VECX,VECY,VECZ, &
     &          BIGM,NBIGM, &
     &          ML,&
     &          TOTELE, &
     &          FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &          N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI, &
     &          DISOPT2,DISOPN,DT,THETA1,BETA,&
     &          GEOBAL,&
     &          LUMP,CONVIS, &
     &          ABSLUM,&
     &          NDGLNO, &
     &          SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &          DENPT,&
     &          MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &          ML2MXX,ML2MXY,ML2MXZ,ML2MYY,ML2MYZ,ML2MZZ, &
     &          XABSOR,YABSOR,ZABSOR, &
     &          CVMVEC,&
     &          LES4TH,&
     &          LESNOD,&
     &          SCFACTH0,&
     &          nnodp, para, halo_tag,&
     &      ISPHERE)
! This subroutine forms the discretised momentum equations. 
! the discretisation is based on Galerkin Least squares of the mom-eqns only. 
! The discretisation scheme is controlled through DISOPT(DISCRETISATION OPTION) 
! All methods use a balancing diffusion term and THETA time stepping. 
! The following are absolute values of DISOPT...
! Galerkin with no vertical derivatives: 
! DISOPT=144 Galerkin LES but no derivatives in the vertical.
! DISOPT=145 Galerkin LES but no derivatives in the vertical and 4th order LES
! Galerkin: 
! DISOPT=146- no LES
! DISOPT=147 -LES but in tensor form like hart3d.
! DISOPT=148 -LES 4th order version of 147.
! Least Squares without time term in the Least Squares weighting: 
! DISOPT=149- no LES
! DISOPT=150 -LES but in tensor form like hart3d.
! DISOPT=151 -LES 4th order version of 147.
! Least Squares weighting: 
! DISOPT=152- no LES
! DISOPT=153 -LES but in tensor form like hart3d.
! DISOPT=154 -LES 4th order version of 147.
!
!
! if DISOPN.ne.0 then set advection to zero and treat advection 
! using the high res method.
! If LES4th=1 then wehn DISOPT=45 USE A 4TH order Smagorinsky. 
!    with this option it is best to use THETA=0.5, ITHETA=0.5 and ITINOI>=2
! SORN contains the div q field NOT N_i * div q. (IS NOT USED NOW)
! GEOBAL=option to use a Geostophic balance solver instead of 
! putting corriolis in here.
! IF SYM then make BIGM symmetrix .
! IF LUMP then lump the mass matrix else dont. 
! BETA controls the conservative discretisation 
! BETA=0. is the usual advection form, BETA=1. is divergence form. 
! BETA=0.5 is an average of the divergence and advection form.(recommended).
! BETA\in[-2,-1] then use the integrated by parts form of the advection term.
! this must be used with PREOPT=1 because we do not want pressure  
! BETA=-1 is non-conservative form
! BETA=-2 is conservative form (recommended)
! BETA=-1.5 is a compromise.
! ABSORBX =absorbtion of energy of x-commentum equation. 
! CVMVEC=virtual mass vector (is 0 in non 2-phase flow). 
! MUPTXX,MUPTYY =components of diffusivities. 
! ML2MXX,ML2MXY etc contain the length scales for LES. 
! NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
! UG,VG,WG are the grid velocities. 
! R0=(NOT NEEDED IN CARTEASEAN COORDS) 
! the minimum distance from ocean bed to centre of earth.
! D0=H/r0 where H is the height of ocean above R0 (H not used here)
! For dimensional run D0=0.
! NB the centrifugal force is \row*(distance from line (x,y)=(0,0)
! NDGLNO=element pter for unknowns. 
! SONDGL=element pter for materials and sources. 
! VONDGL=element pter for advection NU,UG.
! XONDGL=element pter for coordinates(x,y,z)
!
! If CONVIS then assume the conservation of viscosity due to div q term, 
! if div q is not zero then this must be set to .TRUE. 
! NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
! UG,VG,WG are the grid velocities. 
! NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
! THE BLOCKS OF THE MATRIX BIGM ARE ARRANGED AS FOLLOWS
!   WHEN BLOCK SYMMETRIC ...
!   1 * *
!   2 3 *
!   4 5 6 
!   WHEN **NOT** BLOCK SYMMETRIC ...
!   1 2 3
!   4 5 6
!   7 8 9                     *- MATRIX NOT STORED.
    REAL RWIND
       LOGICAL VERTICAL_BOY
    PARAMETER(RWIND=0.5,VERTICAL_BOY=.true.)
! if VERTICAL_BOY only add in a vertical bouyancy term.. which acts like absorption. 
    INTEGER SNONOD,VNONOD,XNONOD
    INTEGER DISOPT2,DISOPT,DISOPN,GEOBAL
    INTEGER NGI,NLOC,NBIGM,NCOLM
    REAL DENPT(SNONOD)
    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
    REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
    REAL ML2MXX(SNONOD),ML2MXY(SNONOD),ML2MXZ(SNONOD)
    REAL ML2MYY(SNONOD),ML2MYZ(SNONOD),ML2MZZ(SNONOD)
    REAL XABSOR(SNONOD),YABSOR(SNONOD),ZABSOR(SNONOD)
    REAL CVMVEC(SNONOD)
    REAL DT,THETA,BETA,acctim
    INTEGER TOTELE,NONODS
! If THETA=0.5 then Crank-Nickolson is used in Navier Stokes equ's.
! If THETA=1.0 then bakward Euler is used in NS equns.
! Similarly for the temperature equation. 
    REAL U(NONODS),V(NONODS),W(NONODS)
    REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)

    REAL UD(NGI),VD(NGI),WD(NGI)
    REAL UDDX(NGI),VDDY(NGI),WDDZ(NGI)

    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL DDETWE(NGI)
    REAL SCFACTH0
    INTEGER ISPHERE
! HX,HY,HZ-characteristic length scales in x,y,z directions.
    REAL HXGI,HYGI,HZGI
    REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
    REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
    REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
    REAL MYZDTW(NGI),MZZDTW(NGI)
    REAL MUGIXX,MUGIXY,MUGIXZ,MUGIYY,MUGIYZ,MUGIZZ
    REAL L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ
    REAL DENGI(NGI)

    REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
    REAL BIGM(NBIGM)
    REAL ML(NONODS)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
       INTEGER NBUOY
       REAL BOY_ML(NBUOY*9*NONODS)
! IF NBUOY=1 treat buoyancy implcicitly else =0
! soxgra,soygra,sozgra is the direction of gravity force when geobal.le.-10.or.(HYDROS.NE.0)
       REAL SOXGRA(SNONOD),SOYGRA(SNONOD),SOZGRA(SNONOD)
       
    REAL VLND,VLM
    REAL HINX,HINY,HINZ
    REAL VLKTEN
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER ELE,ILOC,JLOC,L
    INTEGER GI,COUNT,I,GLOBI,GLOBJ,GLOBIS,GLOBJS
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER CENTRM(NONODS)

    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
! For 4th order LES (LES4TH=1)
    INTEGER LES4TH
! if LESNOD=1 then have a node-wise material 
! for LES - nec. for proper 4th order LES.
    INTEGER LESNOD
    REAL, dimension(:,:,:), allocatable :: EV

    REAL XD(NGI),YD(NGI),ZD(NGI)

    REAL ABLUMP
    LOGICAL ABSLUM
    LOGICAL LUMP,MAKSYM,CONVIS
    LOGICAL BLKSYM
! Ocean stuff
    REAL D0
    REAL COLUMP,R2,RNN
    INTEGER IGLV,ICENT,INUM
    INTEGER IBL11,IBL12,IBL13
    INTEGER IBL21,IBL22,IBL23
    INTEGER IBL31,IBL32,IBL33
    REAL LOVERQ,BETA1,BETA2,ADV1,ADV2
        
    REAL::OMEGAGI,TEN,ADVILOC=0.0,ADVJLOC,RDIFF
    LOGICAL BALDIF,BALGLS
    REAL RLHNLNNC11,RLHNLNNC12
    REAL RLHNLNNC21,RLHNLNNC22
    REAL RLHNLNNC33
    REAL LHNLNNC11,LHNLNNC12,LHNLNNC13
    REAL LHNLNNC21,LHNLNNC22,LHNLNNC23
    REAL LHNLNNC31,LHNLNNC32,LHNLNNC33
    REAL LHNLN11,LHNLN12,LHNLN13
    REAL LHNLN21,LHNLN22,LHNLN23
    REAL LHNLN31,LHNLN32,LHNLN33
         
    REAL RLHNN11,RLHNN12
    REAL RLHNN21,RLHNN22
    REAL RLHNN33
    REAL LHNN11,LHNN12,LHNN13
    REAL LHNN21,LHNN22,LHNN23
    REAL LHNN31,LHNN32,LHNN33
         
    REAL NLN11,NLN12,NLN13
    REAL NLN21,NLN22,NLN23
    REAL NLN31,NLN32,NLN33
    REAL NLNNC11,NLNNC12,NLNNC13
    REAL NLNNC21,NLNNC22,NLNNC23
    REAL NLNNC31,NLNNC32,NLNNC33
    REAL TAKMASX,TAKMASY,TAKMASZ
    REAL VLMXAB,VLMYAB,VLMZAB
    real::BETBAL=-1.0
    REAL XNDMUILOC,YNDMUILOC,ZNDMUILOC
    REAL XNDMUJLOC,YNDMUJLOC,ZNDMUJLOC
    REAL MUPTXXGI,MUPTXYGI,MUPTXZGI,MUPTYYGI,MUPTYZGI,MUPTZZGI
    REAL VLN,VLMOME
    REAL LHNN11C2,LHNN22C2,VLMC2
    REAL RLES,RN
         REAL VLMABXX,VLMABXY,VLMABXZ
         REAL VLMABYX,VLMABYY,VLMABYZ
         REAL VLMABZX,VLMABZY,VLMABZZ
    integer counter

    integer, intent(in) :: nnodp,para,halo_tag

    INTEGER II,JJ,KLOC,KNOD
    REAL ABTHETA,DTRABS
!     ADJOINT
    REAL THETA1
    REAL TENXAB,TENYAB,TENZAB
    REAL XABTEN,YABTEN,ZABTEN
    REAL XABGITOT,YABGITOT,ZABGITOT
    REAL L1(NGI),L2(NGI),L3(NGI),L4(NGI),PGWEIG(NGI) 
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI) 
    REAL XABGI(NGI),YABGI(NGI),ZABGI(NGI)
    REAL XABSDX(NGI),XABSDY(NGI),XABSDZ(NGI)
    REAL YABSDX(NGI),YABSDY(NGI),YABSDZ(NGI)
    REAL ZABSDX(NGI),ZABSDY(NGI),ZABSDZ(NGI)
    REAL GAMB(NGI),GAMBAL(NGI)
    REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
    REAL NLOCIDENT(NLOC,NLOC)
           REAL BOY_ABSXX(NGI),BOY_ABSXY(NGI),BOY_ABSXZ(NGI)
           REAL BOY_ABSYX(NGI),BOY_ABSYY(NGI),BOY_ABSYZ(NGI)
           REAL BOY_ABSZX(NGI),BOY_ABSZY(NGI),BOY_ABSZZ(NGI)
           REAL RSOXGRAGI,RSOYGRAGI,RSOZGRAGI
           REAL GRAV_DEN_NOD_NX,GRAV_DEN_NOD_NY,GRAV_DEN_NOD_NZ
           REAL GRAV_DEN_NOD,RMID,GRA_SIG_CHAN
           REAL OVER_RELAX_BOY
! Local coords...
    INTEGER NCLOC
    REAL, ALLOCATABLE, DIMENSION(:,:)::NC
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCLX
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCLY
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCLZ
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCX
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCY
    REAL, ALLOCATABLE, DIMENSION(:,:)::NCZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDUDZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDVDZ
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDX
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDY
    REAL, ALLOCATABLE, DIMENSION(:)::VEDWDZ
    logical, parameter :: using_lans = .false.
         
         BOY_ML=0.0
    
    ewrite(1,*) 'in subroutine diff3dsimp' 

    allocate(EV(6,6,LESNOD*NONODS))

    NLOCIDENT(1:NLOC,1:NLOC)=0.0
    do I=1,NLOC
       NLOCIDENT(I,I)=1.0
    END DO
       
    THETA=THETA1
    DISOPT=DISOPT2
    IF(DISOPT2.EQ.47) DISOPT=147

!    IF(DISOPT2.EQ.48) DISOPT=147
!    IF(DISOPT2.EQ.148) DISOPT=147
    if(disopt2==48 .or. disopt2==148) then
       ewrite(-1,*) "You have DISOPT=48 or DISOPT=148."
       ewrite(-1,*) "This option has not been coded yet."
       ewrite(-1,*) "If this option used to work for you, that is "
       ewrite(-1,*) "because your DISOPT used to be reset to 147 "
       ewrite(-1,*) "(see the beginning of diff3dsimp)."
       ewrite(-1,*) "I have replaced this hardcoded resetting of "
       ewrite(-1,*) "your options with a warning because I thought "
       ewrite(-1,*) "you ought to know about it."
       ewrite(-1,*) "Best wishes, Jemma"
       ewrite(-1,*) "P.S. To continue as before, just set DISOPT to 147."
       FLAbort('Option not coded yet.')
    end if

    if((DISOPT.GE.149).AND. ( &
          have_option("/physical_parameters/coriolis/sine_of_latitude") .or. &
          have_option("/physical_parameters/coriolis/on_sphere"))) then
       FLAbort('sub does not work for non-uniform omega yet - subtraction')
    endif
    if(using_lans) then
       FLAbort('LANS not in diff3dsimp yet, contact Colin to request')
    end if

    NCLOC=4
    IF(ISPHERE.EQ.2) NCLOC=10
    ALLOCATE(NC(NCLOC,NGI))
    ALLOCATE(NCLX(NCLOC,NGI))
    ALLOCATE(NCLY(NCLOC,NGI))
    ALLOCATE(NCLZ(NCLOC,NGI))
    ALLOCATE(NCX(NCLOC,NGI))
    ALLOCATE(NCY(NCLOC,NGI))
    ALLOCATE(NCZ(NCLOC,NGI))
    IF(ISPHERE.NE.0) THEN
       CALL TRIQUA(L1, L2, L3, L4, PGWEIG, .true.,NGI)
! Work out the shape functions and there derivatives...
       CALL SHATRI(L1, L2, L3, L4, PGWEIG, .true., &
     &            NCLOC,NGI,                          &
     &            NC,NCLX,NCLY,NCLZ) 
    ELSE
       NC = N
       NCLX = NLX
       NCLY = NLY
       NCLZ = NLZ
    ENDIF
       
!    IF((DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154)) THEN
       IF((DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154).or.(DISOPT.EQ.145)) THEN
! calculate the derivatives of NU,NV,NW wrt X,Y,Z
       ewrite(-1,*) 'just before LESDERIVS -option not ready yet'
       stop 4933
       ALLOCATE(VEDUDX(NONODS))
       ALLOCATE(VEDUDY(NONODS))
       ALLOCATE(VEDUDZ(NONODS))
       ALLOCATE(VEDVDX(NONODS))
       ALLOCATE(VEDVDY(NONODS))
       ALLOCATE(VEDVDZ(NONODS))
       ALLOCATE(VEDWDX(NONODS))
       ALLOCATE(VEDWDY(NONODS))
       ALLOCATE(VEDWDZ(NONODS))
       CALL LESDERIVS(&
     &          VEDUDX,VEDUDY,VEDUDZ,&
     &          VEDVDX,VEDVDY,VEDVDZ,&
     &          VEDWDX,VEDWDY,VEDWDZ,&
     &          ML,&
     &          NU,NV,NW,&
     &          X,Y,Z, XONDGL, TOTELE,NONODS,NLOC,NGI, &
     &          N,NLX,NLY,NLZ, WEIGHT, DETWEI,DDETWE, .TRUE.,.FALSE.,&
     &          NX,NY,NZ, &
     &          NDGLNO, nnodp,para,halo_tag,VECX,VECY,&
     &          LESNOD,EV,&
     &          ML2MXX,ML2MXY,ML2MXZ,ML2MYY,ML2MYZ,ML2MZZ, &
     &          UG,VG,WG, &
     &          VONDGL)
       ewrite(1,*) 'just after LESDERIVS'
    ENDIF

    CALL PMINMX(ML2MXX,NONODS,'******ML2MXX  ')
    CALL PMINMX(ML2MXY,NONODS,'******ML2MXY  ')
    CALL PMINMX(ML2MXZ,NONODS,'******ML2MXZ  ')
    CALL PMINMX(ML2MYY,NONODS,'******ML2MYY  ')
    CALL PMINMX(ML2MYZ,NONODS,'******ML2MYZ  ')
    CALL PMINMX(ML2MZZ,NONODS,'******ML2MZZ  ')
          
    IBL11=0
    IBL12=1*NCOLM
    IBL13=2*NCOLM
    IBL21=3*NCOLM
    IBL22=4*NCOLM
    IBL23=5*NCOLM
    IBL31=6*NCOLM
    IBL32=7*NCOLM
    IBL33=8*NCOLM

    ML=0.0
!    IF(.NOT.((GEOBAL.LE.-10).AND.(DISOPT.GE.146)&
!     &                             .AND.(DISOPT.LE.148))) THEN
           IF(.NOT.((GEOBAL.LE.-10).AND.(DISOPT.GE.144) &
      &                            .AND.(DISOPT.LE.148))) THEN
! Dont set to zero if Galerkin and visited GEOELI1P
       VECX=0.0
       VECY=0.0
       VECZ=0.0
    ENDIF
    BIGM(1:NBIGM)=0.0

    do  ELE=1,TOTELE! Was loop 340
! Get determinant and derivatives...
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
     &               N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
     &               NX,NY,NZ, NCX,NCY,NCZ,&
     &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
     &               XD,YD,ZD,&
     &               ISPHERE) 
     
      do  GI=1,NGI! Was loop 331
      
             IF(ISPHERE.LE.1) THEN
               DIRX(GI)=0.0
               DIRY(GI)=0.0
               DIRZ(GI)=1.0
             ELSE
               DIRX(GI)=XD(GI)
               DIRY(GI)=YD(GI)
               DIRZ(GI)=ZD(GI)
               RN=SQRT(DIRX(GI)**2+DIRY(GI)**2+DIRZ(GI)**2)
               DIRX(GI)=DIRX(GI)/RN
               DIRY(GI)=DIRY(GI)/RN
               DIRZ(GI)=DIRZ(GI)/RN
             ENDIF

         UD(GI)=0.    
         VD(GI)=0.
         WD(GI)=0.
         UDDX(GI)=0.    
         VDDY(GI)=0.
         WDDZ(GI)=0.
             
         XABGI(GI)=0.0
         YABGI(GI)=0.0
         ZABGI(GI)=0.0
             
         XABSDX(GI)=0.0
         XABSDY(GI)=0.0
         XABSDZ(GI)=0.0
             
         YABSDX(GI)=0.0
         YABSDY(GI)=0.0
         YABSDZ(GI)=0.0
             
         ZABSDX(GI)=0.0
         ZABSDY(GI)=0.0
         ZABSDZ(GI)=0.0
             
         MUPTXXGI=0.0
         MUPTXYGI=0.0
         MUPTXZGI=0.0
         MUPTYYGI=0.0
         MUPTYZGI=0.0
         MUPTZZGI=0.0
             
               RSOXGRAGI=0.0
               RSOYGRAGI=0.0
               RSOZGRAGI=0.0
               
               GRAV_DEN_NOD_NX=0.0
               GRAV_DEN_NOD_NY=0.0
               GRAV_DEN_NOD_NZ=0.0

         do  L=1,NLOC! Was loop 79

            IGLV=VONDGL((ELE-1)*NLOC+L)

            UD(GI)=UD(GI) + N(L,GI)*(NU(IGLV)-UG(IGLV))
            VD(GI)=VD(GI) + N(L,GI)*(NV(IGLV)-VG(IGLV))
            WD(GI)=WD(GI) + N(L,GI)*(NW(IGLV)-WG(IGLV))
        
            UDDX(GI)=UDDX(GI) + NX(L,GI)*(NU(IGLV)-UG(IGLV))
            VDDY(GI)=VDDY(GI) + NY(L,GI)*(NV(IGLV)-VG(IGLV))
            WDDZ(GI)=WDDZ(GI) + NZ(L,GI)*(NW(IGLV)-WG(IGLV))
              
            IF(.NOT.ABSLUM) THEN
               XABGI(GI)=XABGI(GI)+N(L,GI)*XABSOR(IGLV)
               YABGI(GI)=YABGI(GI)+N(L,GI)*YABSOR(IGLV)
               ZABGI(GI)=ZABGI(GI)+N(L,GI)*ZABSOR(IGLV)
               
               XABSDX(GI)=XABSDX(GI)+NX(L,GI)*XABSOR(IGLV)
               XABSDY(GI)=XABSDY(GI)+NY(L,GI)*XABSOR(IGLV)
               XABSDZ(GI)=XABSDZ(GI)+NZ(L,GI)*XABSOR(IGLV)
              
               YABSDX(GI)=YABSDX(GI)+NX(L,GI)*YABSOR(IGLV)
               YABSDY(GI)=YABSDY(GI)+NY(L,GI)*YABSOR(IGLV)
               YABSDZ(GI)=YABSDZ(GI)+NZ(L,GI)*YABSOR(IGLV)
              
               ZABSDX(GI)=ZABSDX(GI)+NX(L,GI)*ZABSOR(IGLV)
               ZABSDY(GI)=ZABSDY(GI)+NY(L,GI)*ZABSOR(IGLV)
               ZABSDZ(GI)=ZABSDZ(GI)+NZ(L,GI)*ZABSOR(IGLV)
            ENDIF
              
            MUPTXXGI=MUPTXXGI+N(L,GI)*MUPTXX(IGLV)
            MUPTXYGI=MUPTXYGI+N(L,GI)*MUPTXY(IGLV)
            MUPTXZGI=MUPTXZGI+N(L,GI)*MUPTXZ(IGLV)
            MUPTYYGI=MUPTYYGI+N(L,GI)*MUPTYY(IGLV)
            MUPTYZGI=MUPTYZGI+N(L,GI)*MUPTYZ(IGLV)
            MUPTZZGI=MUPTZZGI+N(L,GI)*MUPTZZ(IGLV)
              IF(NBUOY.NE.0) THEN
               RSOXGRAGI=RSOXGRAGI+N(L,GI)*SOXGRA(IGLV)
               RSOYGRAGI=RSOYGRAGI+N(L,GI)*SOYGRA(IGLV)
               RSOZGRAGI=RSOZGRAGI+N(L,GI)*SOZGRA(IGLV)
               
               GRAV_DEN_NOD=SQRT(SOXGRA(IGLV)**2+SOYGRA(IGLV)**2+SOZGRA(IGLV)**2)
               
               GRAV_DEN_NOD_NX=GRAV_DEN_NOD_NX + NX(L,GI)*GRAV_DEN_NOD
               GRAV_DEN_NOD_NY=GRAV_DEN_NOD_NY + NY(L,GI)*GRAV_DEN_NOD
               GRAV_DEN_NOD_NZ=GRAV_DEN_NOD_NZ + NZ(L,GI)*GRAV_DEN_NOD
              ENDIF
         end do ! Was loop 79

              IF(NBUOY.NE.0) THEN
! Calculate the absorbtion BOY_ABSXX etc associated with vertical buoyancy...               
               GRA_SIG_CHAN= &
     &                    SIGN(1.0,DIRX(GI)*RSOXGRAGI+DIRY(GI)*RSOYGRAGI &
     &                     +DIRZ(GI)*RSOZGRAGI)   
               
               OVER_RELAX_BOY=1.0
               GRAV_DEN_NOD_NX=GRA_SIG_CHAN*GRAV_DEN_NOD_NX
               GRAV_DEN_NOD_NY=GRA_SIG_CHAN*GRAV_DEN_NOD_NY
               GRAV_DEN_NOD_NZ=GRA_SIG_CHAN*GRAV_DEN_NOD_NZ 
               RMID=DT* MAX(0.0, DIRX(GI)*GRAV_DEN_NOD_NX+DIRY(GI)*GRAV_DEN_NOD_NY &
     &                     +DIRZ(GI)*GRAV_DEN_NOD_NZ )*OVER_RELAX_BOY
!               RMID=DT* abs((DIRX(GI)*GRAV_DEN_NOD_NX+DIRY(GI)*GRAV_DEN_NOD_NY &
!     &                     +DIRZ(GI)*GRAV_DEN_NOD_NZ )*OVER_RELAX_BOY)
              ELSE
               RMID=0.0
              ENDIF
           IF(NBUOY.NE.0) THEN
             IF(VERTICAL_BOY) THEN
! only add in a vertical bouyancy term.. which acts like absorption. 
               BOY_ABSXX(GI)=RMID*DIRX(GI)*DIRX(GI)
               BOY_ABSXY(GI)=RMID*DIRX(GI)*DIRY(GI)
               BOY_ABSXZ(GI)=RMID*DIRX(GI)*DIRZ(GI)
               
               BOY_ABSYX(GI)=RMID*DIRY(GI)*DIRX(GI)
               BOY_ABSYY(GI)=RMID*DIRY(GI)*DIRY(GI)
               BOY_ABSYZ(GI)=RMID*DIRY(GI)*DIRZ(GI)
               
               BOY_ABSZX(GI)=RMID*DIRZ(GI)*DIRX(GI)
               BOY_ABSZY(GI)=RMID*DIRZ(GI)*DIRY(GI)
               BOY_ABSZZ(GI)=RMID*DIRZ(GI)*DIRZ(GI)
             ELSE
               RMID=DT
               BOY_ABSXX(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NX
               BOY_ABSXY(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NY
               BOY_ABSXZ(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NZ
               
               BOY_ABSYX(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NX
               BOY_ABSYY(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NY
               BOY_ABSYZ(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NZ
               
               BOY_ABSZX(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NX
               BOY_ABSZY(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NY
               BOY_ABSZZ(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NZ
             ENDIF
           ELSE
               BOY_ABSXX(GI)=0.0
               BOY_ABSXY(GI)=0.0
               BOY_ABSXZ(GI)=0.0
               
               BOY_ABSYX(GI)=0.0
               BOY_ABSYY(GI)=0.0
               BOY_ABSYZ(GI)=0.0
               
               BOY_ABSZX(GI)=0.0
               BOY_ABSZY(GI)=0.0
               BOY_ABSZZ(GI)=0.0
           ENDIF
!
! calculate tensor for oceans:
          
!         IF((DISOPT.EQ.147).or.(DISOPT.EQ.150).or.(DISOPT.EQ.153)&
!     & .or.(DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154)) THEN
     
        IF((DISOPT.EQ.147).or.(DISOPT.EQ.150).or.(DISOPT.EQ.153) &
     & .or.(DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154) &
     & .or.(DISOPT.EQ.144).or.(DISOPT.EQ.145)) THEN

! based on anisotropic length scales...
            CALL ANO_OCEAN_SCAL(&
     &    HIMXXDTW(GI),HIMXYDTW(GI),HIMXZDTW(GI),&
     &    HIMYYDTW(GI),HIMYZDTW(GI),HIMZZDTW(GI),&
     &    MXXDTW(GI),MXYDTW(GI),MXZDTW(GI),&
     &    MYYDTW(GI),MYZDTW(GI),MZZDTW(GI),&
     &    XD(GI),YD(GI),ZD(GI),&
     &    ELE,NLOC,TOTELE,XONDGL,XNONOD,X,Y,Z,ISPHERE) 
     
            CALL HALF_LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             HIMXXDTW,HIMXYDTW,HIMXZDTW,&
     &             HIMYYDTW,HIMYZDTW,HIMZZDTW,&
     &             MXXDTW,MXYDTW,MXZDTW,&
     &             MYYDTW,MYZDTW,MZZDTW, &
     &             (DISOPT.EQ.144).OR.(DISOPT.EQ.145),DIRX,DIRY,DIRZ)

            HIMXXDTW(GI)=HIMXXDTW(GI)*DETWEI(GI)
            HIMXYDTW(GI)=HIMXYDTW(GI)*DETWEI(GI)
            HIMXZDTW(GI)=HIMXZDTW(GI)*DETWEI(GI)

            HIMYYDTW(GI)=HIMYYDTW(GI)*DETWEI(GI)
            HIMYZDTW(GI)=HIMYZDTW(GI)*DETWEI(GI)

            HIMZZDTW(GI)=HIMZZDTW(GI)*DETWEI(GI)
        
            MUGIXX=MXXDTW(GI)
            MUGIXY=MXYDTW(GI)
            MUGIXZ=MXZDTW(GI)
            MUGIYY=MYYDTW(GI)
            MUGIYZ=MYZDTW(GI)
            MUGIZZ=MZZDTW(GI)
         ELSE

            MUGIXX=0.0
            MUGIXY=0.0
            MUGIXZ=0.0
            MUGIYY=0.0
            MUGIYZ=0.0
            MUGIZZ=0.0
         ENDIF

         RLES=0.1
         MXXDTW(GI)=(RLES*MUGIXX+MUPTXXGI)*DETWEI(GI)
         MXYDTW(GI)=(RLES*MUGIXY+MUPTXYGI)*DETWEI(GI)
         MXZDTW(GI)=(RLES*MUGIXZ+MUPTXZGI)*DETWEI(GI)

         MYYDTW(GI)=(RLES*MUGIYY+MUPTYYGI)*DETWEI(GI)
         MYZDTW(GI)=(RLES*MUGIYZ+MUPTYZGI)*DETWEI(GI)

         MZZDTW(GI)=(RLES*MUGIZZ+MUPTZZGI)*DETWEI(GI)
! for LES...
      end do ! Was loop 331

      do ILOC=1,NLOC! Was loop 350
         GLOBI =NDGLNO((ELE-1)*NLOC+ILOC)
         GLOBIS=SONDGL((ELE-1)*NLOC+ILOC)

         do JLOC=1,NLOC! Was loop 360
            GLOBJ =NDGLNO((ELE-1)*NLOC+JLOC)
            GLOBJS=SONDGL((ELE-1)*NLOC+JLOC)

            VLM=0. 
             
            HINX=0.0
            HINY=0.0
            HINZ=0.0
        
            VLKTEN=0.0
        
            LHNLNNC11=0.0
            LHNLNNC12=0.0
            LHNLNNC13=0.0
             
            LHNLNNC21=0.0
            LHNLNNC22=0.0
            LHNLNNC23=0.0
             
            LHNLNNC31=0.0
            LHNLNNC32=0.0
            LHNLNNC33=0.0
             
            LHNLN11=0.0
            LHNLN12=0.0
            LHNLN13=0.0
             
            LHNLN21=0.0
            LHNLN22=0.0
            LHNLN23=0.0
             
            LHNLN31=0.0
            LHNLN32=0.0
            LHNLN33=0.0
             
            LHNN11=0.0
            LHNN12=0.0
            LHNN13=0.0
             
            LHNN21=0.0
            LHNN22=0.0
            LHNN23=0.0
             
            LHNN31=0.0
            LHNN32=0.0
            LHNN33=0.0
             
            NLN11=0.0
            NLN12=0.0
            NLN13=0.0
            NLN21=0.0
            NLN22=0.0
            NLN23=0.0
            NLN31=0.0
            NLN32=0.0
            NLN33=0.0
              
            TAKMASX=0.0
            TAKMASY=0.0
            TAKMASZ=0.0
              
            VLND  =0.0
            VLMOME=0.0
             
            VLMXAB=0.0
            VLMYAB=0.0
            VLMZAB=0.0
             
            LHNN11C2=0.0
            LHNN22C2=0.0
             
            VLMC2=0.0
             
             VLMABXX=0.0
             VLMABXY=0.0
             VLMABXZ=0.0
             
             VLMABYX=0.0
             VLMABYY=0.0
             VLMABYZ=0.0
             
             VLMABZX=0.0
             VLMABZY=0.0
             VLMABZZ=0.0
        
            do GI=1,NGI! Was loop 458
           
               OMEGAGI=FUNOME(XD(GI),YD(GI),ZD(GI))
           
               BALDIF=.TRUE.
               BALGLS=.TRUE.
             
               IF((DISOPT.GE.146).AND.(DISOPT.LE.148)) THEN
! Galerkin....
                  BETBAL=0.0
                  GAMB(GI)=0.0
                  BALDIF=.FALSE.
                  BALGLS=.FALSE.
               ELSE IF((DISOPT.GE.149).AND.(DISOPT.LE.151)) THEN
! Galerkin Least Squares BETBAL=0.0....
                  HXGI=ABS(A11(GI)*UD(GI)+A12(GI)*VD(GI)+A13(GI)*WD(GI))
                  HYGI=ABS(A21(GI)*UD(GI)+A22(GI)*VD(GI)+A23(GI)*WD(GI))
                  HZGI=ABS(A31(GI)*UD(GI)+A32(GI)*VD(GI)+A33(GI)*WD(GI))
                  BETBAL=0.0

                  GAMB(GI)=MIN(2./MAX(HXGI,HYGI,HZGI,1.E-7), DT, 1./max(2.*OMEGAGI,1.e-9))
               ELSE IF((DISOPT.GE.152).AND.(DISOPT.LE.154)) THEN
! Galerkin Least Squares 
                  HXGI=ABS(A11(GI)*UD(GI)+A12(GI)*VD(GI)+A13(GI)*WD(GI))
                  HYGI=ABS(A21(GI)*UD(GI)+A22(GI)*VD(GI)+A23(GI)*WD(GI))
                  HZGI=ABS(A31(GI)*UD(GI)+A32(GI)*VD(GI)+A33(GI)*WD(GI))
                  BETBAL=1.0

                  GAMB(GI)=MIN(2./MAX(HXGI,HYGI,HZGI,1.E-7), DT, 1./max(2.*OMEGAGI,1.e-9))
               ENDIF
             
               GAMBAL(GI)=GAMB(GI)  / (1.0+GAMB(GI)*BETBAL/DT)
             
               XABGITOT=XABGI(GI) 
               YABGITOT=YABGI(GI) 
               ZABGITOT=ZABGI(GI) 
          
               VLN=+N(ILOC,GI)*(UD(GI)*NX(JLOC,GI)+VD(GI)*NY(JLOC,GI)&
     &                       +WD(GI)*NZ(JLOC,GI))*DETWEI(GI) 

               VLND=VLND+VLN

               RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
          
               VLM=VLM+RNN 
               
             VLMABXX=VLMABXX+RNN * BOY_ABSXX(GI)
             VLMABXY=VLMABXY+RNN * BOY_ABSXY(GI)
             VLMABXZ=VLMABXZ+RNN * BOY_ABSXZ(GI)
             
             VLMABYX=VLMABYX+RNN * BOY_ABSYX(GI)
             VLMABYY=VLMABYY+RNN * BOY_ABSYY(GI)
             VLMABYZ=VLMABYZ+RNN * BOY_ABSYZ(GI)
             
             VLMABZX=VLMABZX+RNN * BOY_ABSZX(GI)
             VLMABZY=VLMABZY+RNN * BOY_ABSZY(GI)
             VLMABZZ=VLMABZZ+RNN * BOY_ABSZZ(GI)
     
               XNDMUJLOC=(NX(JLOC,GI)*MXXDTW(GI)&
     &                 +NY(JLOC,GI)*MXYDTW(GI)+NZ(JLOC,GI)*MXZDTW(GI))
               YNDMUJLOC=(NX(JLOC,GI)*MXYDTW(GI)&
     &                 +NY(JLOC,GI)*MYYDTW(GI)+NZ(JLOC,GI)*MYZDTW(GI))
               ZNDMUJLOC=(NX(JLOC,GI)*MXZDTW(GI)&
     &                 +NY(JLOC,GI)*MYZDTW(GI)+NZ(JLOC,GI)*MZZDTW(GI))
     
               XNDMUILOC=(NX(ILOC,GI)*MXXDTW(GI)&
     &                 +NY(ILOC,GI)*MXYDTW(GI)+NZ(ILOC,GI)*MXZDTW(GI))
               YNDMUILOC=(NX(JLOC,GI)*MXYDTW(GI)&
     &                 +NY(ILOC,GI)*MYYDTW(GI)+NZ(ILOC,GI)*MYZDTW(GI))
               ZNDMUILOC=(NX(ILOC,GI)*MXZDTW(GI)&
     &                 +NY(ILOC,GI)*MYZDTW(GI)+NZ(ILOC,GI)*MZZDTW(GI))
     
               TEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC+NZ(ILOC,GI)*ZNDMUJLOC

               TENXAB=  +N(ILOC,GI)*(XABSDX(GI)*XNDMUJLOC+XABSDY(GI)*YNDMUJLOC&
     &                              +XABSDZ(GI)*ZNDMUJLOC)  +XABGITOT*TEN
               TENYAB=  +N(ILOC,GI)*(YABSDX(GI)*XNDMUJLOC+YABSDY(GI)*YNDMUJLOC&
     &                              +YABSDZ(GI)*ZNDMUJLOC)  +YABGITOT*TEN
               TENZAB=  +N(ILOC,GI)*(ZABSDX(GI)*XNDMUJLOC+ZABSDY(GI)*YNDMUJLOC&
     &                              +ZABSDZ(GI)*ZNDMUJLOC)  +ZABGITOT*TEN
     
               XABTEN=  +(XABSDX(GI)*XNDMUILOC+XABSDY(GI)*YNDMUILOC&
     &                   +XABSDZ(GI)*ZNDMUILOC)*N(JLOC,GI)  +XABGITOT*TEN
               YABTEN=  +(YABSDX(GI)*XNDMUILOC+YABSDY(GI)*YNDMUILOC&
     &                   +YABSDZ(GI)*ZNDMUILOC)*N(JLOC,GI)  +YABGITOT*TEN
               ZABTEN=  +(ZABSDX(GI)*XNDMUILOC+ZABSDY(GI)*YNDMUILOC&
     &                   +ZABSDZ(GI)*ZNDMUILOC)*N(JLOC,GI)  +ZABGITOT*TEN
          
        
               VLKTEN=VLKTEN+TEN 
           
          IF((DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154).or.(DISOPT.EQ.145)) THEN
!               IF((DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154)) THEN
! Subtract out extended stencil 2nd order operator to get a 4th order operator...
                  VECX(GLOBI)=VECX(GLOBI) + NX(ILOC,GI)*N(JLOC,GI)*( &
     &     MXXDTW(GI)*VEDUDX(GLOBJ) +MXYDTW(GI)*VEDUDY(GLOBJ) +MXZDTW(GI)*VEDUDZ(GLOBJ) )  &
     &                            + NY(ILOC,GI)*N(JLOC,GI)*( &
     &     MXYDTW(GI)*VEDUDX(GLOBJ) +MYYDTW(GI)*VEDUDY(GLOBJ) +MYZDTW(GI)*VEDUDZ(GLOBJ) ) &
     &                            + NZ(ILOC,GI)*N(JLOC,GI)*( &
     &     MXZDTW(GI)*VEDUDX(GLOBJ) +MYZDTW(GI)*VEDUDY(GLOBJ) +MZZDTW(GI)*VEDUDZ(GLOBJ) ) 
     
                  VECY(GLOBI)=VECY(GLOBI) + NX(ILOC,GI)*N(JLOC,GI)*( &
     &     MXXDTW(GI)*VEDVDX(GLOBJ) +MXYDTW(GI)*VEDVDY(GLOBJ) +MXZDTW(GI)*VEDVDZ(GLOBJ) ) &
     &                            + NY(ILOC,GI)*N(JLOC,GI)*( &
     &     MXYDTW(GI)*VEDVDX(GLOBJ) +MYYDTW(GI)*VEDVDY(GLOBJ) +MYZDTW(GI)*VEDVDZ(GLOBJ) ) &
     &                            + NZ(ILOC,GI)*N(JLOC,GI)*( &
     &     MXZDTW(GI)*VEDVDX(GLOBJ) +MYZDTW(GI)*VEDVDY(GLOBJ) +MZZDTW(GI)*VEDVDZ(GLOBJ) ) 
     
                  VECZ(GLOBI)=VECZ(GLOBI) + NX(ILOC,GI)*N(JLOC,GI)*( &
     &     MXXDTW(GI)*VEDWDX(GLOBJ) +MXYDTW(GI)*VEDWDY(GLOBJ) +MXZDTW(GI)*VEDWDZ(GLOBJ) ) &
     &                            + NY(ILOC,GI)*N(JLOC,GI)*( &
     &     MXYDTW(GI)*VEDWDX(GLOBJ) +MYYDTW(GI)*VEDWDY(GLOBJ) +MYZDTW(GI)*VEDWDZ(GLOBJ) ) &
     &                            + NZ(ILOC,GI)*N(JLOC,GI)*( &
     &     MXZDTW(GI)*VEDWDX(GLOBJ) +MYZDTW(GI)*VEDWDY(GLOBJ) +MZZDTW(GI)*VEDWDZ(GLOBJ) ) 
               ENDIF
         
               IF(BALDIF) THEN
                  ADVILOC=UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI)+WD(GI)*NZ(ILOC,GI)
                  ADVJLOC=UD(GI)*NX(JLOC,GI)+VD(GI)*NY(JLOC,GI)+WD(GI)*NZ(JLOC,GI)
                  RDIFF=ADVILOC*ADVJLOC*DETWEI(GI) 
        
                  RLHNLNNC11=GAMBAL(GI)*RDIFF  &
     &          +GAMBAL(GI)*XABGITOT*( XABGITOT*RNN + N(ILOC,GI)*ADVJLOC*DETWEI(GI) &
     &                                          + ADVILOC*N(JLOC,GI)*DETWEI(GI)  )&
     &                      +GAMBAL(GI)*(XABTEN +TENXAB)        
                  RLHNLNNC12=GAMBAL(GI)*( -2.*OMEGAGI*N(ILOC,GI)*ADVJLOC*DETWEI(GI) )&
     &                 +GAMBAL(GI)*( -2.*OMEGAGI*YABGITOT*RNN)&
     &                 +GAMBAL(GI)*( -2.*OMEGAGI*TEN )
                  RLHNLNNC21=GAMBAL(GI)*( +2.*OMEGAGI*N(ILOC,GI)*ADVJLOC*DETWEI(GI) )&
     &                 +GAMBAL(GI)*( +2.*OMEGAGI*XABGITOT*RNN)&
     &                 +GAMBAL(GI)*( +2.*OMEGAGI*TEN )
                  RLHNLNNC22=GAMBAL(GI)*RDIFF  &
     &          +GAMBAL(GI)*YABGITOT*( YABGITOT*RNN + N(ILOC,GI)*ADVJLOC*DETWEI(GI) &
     &                                          + ADVILOC*N(JLOC,GI)*DETWEI(GI)  )&
     &                      +GAMBAL(GI)*(YABTEN +TENYAB)        
                  RLHNLNNC33=GAMBAL(GI)*RDIFF &
     &          +GAMBAL(GI)*ZABGITOT*( ZABGITOT*RNN + N(ILOC,GI)*ADVJLOC*DETWEI(GI) &
     &                                          + ADVILOC*N(JLOC,GI)*DETWEI(GI)  )&
     &                      +GAMBAL(GI)*(ZABTEN +TENZAB)  
        
                  LHNLNNC11=LHNLNNC11+RLHNLNNC11   
                  LHNLNNC12=LHNLNNC12+RLHNLNNC12 
                  LHNLNNC21=LHNLNNC21+RLHNLNNC21 
                  LHNLNNC22=LHNLNNC22+RLHNLNNC22    
                  LHNLNNC33=LHNLNNC33+RLHNLNNC33
        
                  LHNLN11=LHNLN11+RLHNLNNC11&
     &         - GAMBAL(GI)*(4.*OMEGAGI*OMEGAGI*RNN )
                  LHNLN12=LHNLN12+RLHNLNNC12&
     &              +GAMBAL(GI)*(-2.*OMEGAGI*ADVILOC*N(JLOC,GI)*DETWEI(GI) )&
     &              +GAMBAL(GI)*(-2.*OMEGAGI*XABGITOT*RNN)&
     &              +GAMBAL(GI)*(-2.*OMEGAGI*TEN )
                  LHNLN21=LHNLN21+RLHNLNNC21&
     &              +GAMBAL(GI)*(+2.*OMEGAGI*ADVILOC*N(JLOC,GI)*DETWEI(GI) )&
     &              +GAMBAL(GI)*(+2.*OMEGAGI*YABGITOT*RNN)&
     &              +GAMBAL(GI)*(+2.*OMEGAGI*TEN )
                  LHNLN22=LHNLN22+RLHNLNNC22&
     &          -GAMBAL(GI)*(4.*OMEGAGI*OMEGAGI*RNN )
                  LHNLN33=LHNLN33+RLHNLNNC33 
               ENDIF
             
               IF(BALGLS) THEN
                  ! If this assert fails when ADVILOC would be used
                  ! uninitialised.
                  assert(BALDIF)
                  RLHNN11=GAMBAL(GI)*(ADVILOC*N(JLOC,GI)*DETWEI(GI) +TEN +XABGITOT*RNN)
                  RLHNN12=-GAMBAL(GI)*2.*OMEGAGI*RNN
                  RLHNN21=+GAMBAL(GI)*2.*OMEGAGI*RNN
                  RLHNN22=GAMBAL(GI)*(ADVILOC*N(JLOC,GI)*DETWEI(GI) +TEN +YABGITOT*RNN)
                  RLHNN33=GAMBAL(GI)*(ADVILOC*N(JLOC,GI)*DETWEI(GI) +TEN +ZABGITOT*RNN)
     
                  LHNN11=LHNN11+RLHNN11
                  LHNN12=LHNN12+RLHNN12
                  LHNN21=LHNN21+RLHNN21
                  LHNN22=LHNN22+RLHNN22
                  LHNN33=LHNN33+RLHNN33
             
                  LHNN11C2=LHNN11C2+RLHNN11*2.*OMEGAGI
                  LHNN22C2=LHNN22C2+RLHNN22*2.*OMEGAGI
               ENDIF
         
               NLN11=NLN11 + (VLN+TEN +XABGI(GI)*RNN) 
               NLN12=NLN12 - (2.*OMEGAGI*RNN)         
               NLN21=NLN21 + (2.*OMEGAGI*RNN)        
               NLN22=NLN22 + (VLN+TEN +YABGI(GI)*RNN) 
               NLN33=NLN33 + (VLN+TEN +ZABGI(GI)*RNN) 
             
               VLMOME=VLMOME+OMEGAGI*RNN
             
               VLMC2=VLMC2+RNN*2.*OMEGAGI
             
               VLMXAB=VLMXAB+XABGI(GI)*RNN
               VLMYAB=VLMYAB+YABGI(GI)*RNN
               VLMZAB=VLMZAB+ZABGI(GI)*RNN

            end do ! Was loop 458
             
            NLNNC11=NLN11
            NLNNC12=0.0
            NLNNC13=0.0
            NLNNC21=0.0
            NLNNC22=NLN22
            NLNNC23=0.0
            NLNNC31=0.0
            NLNNC32=0.0
            NLNNC33=NLN33

! Find count. ***************************************
            CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
! lump the mass matrix.
            ICENT=CENTRM(GLOBI)
           
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...
            VECX(GLOBI)=VECX(GLOBI) &
     &        -(NLNNC11+THETA*LHNLNNC11)*U(GLOBJ)&
     &        -(NLNNC12+THETA*LHNLNNC12)*V(GLOBJ)&
     &        -(NLNNC13+THETA*LHNLNNC13)*W(GLOBJ)&
     & +(VLM+THETA*LHNN11)*SOURCX(GLOBJ)      +THETA*LHNN12*SOURCY(GLOBJ)      +THETA*LHNN13*SOURCZ(GLOBJ)&
! Buoyancy stabilisation term...
     & +VLMABXX*(NU(GLOBI)-U(GLOBI))+VLMABXY*(NV(GLOBI)-V(GLOBI))+VLMABXZ*(NW(GLOBI)-W(GLOBI)) &
! added term to get good time accuracy of Coriolis...
     & +(VLMC2+THETA*LHNN11C2)*(V(GLOBJ)-NV(GLOBJ))    
            VECY(GLOBI)=VECY(GLOBI) &
     &        -(NLNNC21+THETA*LHNLNNC21)*U(GLOBJ)&
     &        -(NLNNC22+THETA*LHNLNNC22)*V(GLOBJ)&
     &        -(NLNNC23+THETA*LHNLNNC23)*W(GLOBJ)&
     &       +THETA*LHNN21*SOURCX(GLOBJ)+(VLM+THETA*LHNN22)*SOURCY(GLOBJ)      +THETA*LHNN23*SOURCZ(GLOBJ)&
! Buoyancy stabilisation term...
     & +VLMABYX*(NU(GLOBI)-U(GLOBI))+VLMABYY*(NV(GLOBI)-V(GLOBI))+VLMABYZ*(NW(GLOBI)-W(GLOBI)) &
! added term to get good time accuracy of Coriolis...
     &                                  -(VLMC2+THETA*LHNN22C2)*(U(GLOBJ)-NU(GLOBJ)) 
            VECZ(GLOBI)=VECZ(GLOBI) &
     &        -(NLNNC31+THETA*LHNLNNC31)*U(GLOBJ)&
     &        -(NLNNC32+THETA*LHNLNNC32)*V(GLOBJ)&
     &        -(NLNNC33+THETA*LHNLNNC33)*W(GLOBJ)&
     &   +THETA*LHNN31*SOURCX(GLOBJ)          +THETA*LHNN32*SOURCY(GLOBJ)+(VLM+THETA*LHNN33)*SOURCZ(GLOBJ)&
! Buoyancy stabilisation term...
     & +VLMABZX*(NU(GLOBI)-U(GLOBI))+VLMABZY*(NV(GLOBI)-V(GLOBI))+VLMABZZ*(NW(GLOBI)-W(GLOBI))
          
          IF(NBUOY.NE.0) THEN
            boy_ml(globi)         =boy_ml(globi)         +VLMABXX*THETA*DT
            boy_ml(globi+NONODS)  =boy_ml(globi+NONODS)  +VLMABXY*THETA*DT
            boy_ml(globi+2*NONODS)=boy_ml(globi+2*NONODS)+VLMABXZ*THETA*DT
          
            boy_ml(globi+3*NONODS)=boy_ml(globi+3*NONODS)+VLMABYX*THETA*DT
            boy_ml(globi+4*NONODS)=boy_ml(globi+4*NONODS)+VLMABYY*THETA*DT
            boy_ml(globi+5*NONODS)=boy_ml(globi+5*NONODS)+VLMABYZ*THETA*DT
          
            boy_ml(globi+6*NONODS)=boy_ml(globi+6*NONODS)+VLMABZX*THETA*DT
            boy_ml(globi+7*NONODS)=boy_ml(globi+7*NONODS)+VLMABZY*THETA*DT
            boy_ml(globi+8*NONODS)=boy_ml(globi+8*NONODS)+VLMABZZ*THETA*DT
          ENDIF
          
! lump the mass matrix.(LUMP=1. if we are to LUMP the mass matrix)
            BIGM(ICENT+IBL11)=BIGM(ICENT+IBL11)  + VLM 
            BIGM(ICENT+IBL22)=BIGM(ICENT+IBL22)  + VLM 
            BIGM(ICENT+IBL33)=BIGM(ICENT+IBL33)  + VLM 
! Buoyancy stabilisation term...
         BIGM(ICENT+IBL11)=BIGM(ICENT+IBL11)  + DT*THETA*VLMABXX
         BIGM(ICENT+IBL12)=BIGM(ICENT+IBL12)  + DT*THETA*VLMABXY
         BIGM(ICENT+IBL13)=BIGM(ICENT+IBL13)  + DT*THETA*VLMABXZ
         
         BIGM(ICENT+IBL21)=BIGM(ICENT+IBL21)  + DT*THETA*VLMABYX
         BIGM(ICENT+IBL22)=BIGM(ICENT+IBL22)  + DT*THETA*VLMABYY
         BIGM(ICENT+IBL23)=BIGM(ICENT+IBL23)  + DT*THETA*VLMABYZ
         
         BIGM(ICENT+IBL31)=BIGM(ICENT+IBL31)  + DT*THETA*VLMABZX
         BIGM(ICENT+IBL32)=BIGM(ICENT+IBL32)  + DT*THETA*VLMABZY
         BIGM(ICENT+IBL33)=BIGM(ICENT+IBL33)  + DT*THETA*VLMABZZ
         
! Balancing tensor form...
            BIGM(COUNT+IBL11)=BIGM(COUNT+IBL11)  +THETA*(DT*NLN11+LHNN11+DT*THETA*LHNLN11)
            BIGM(COUNT+IBL12)=BIGM(COUNT+IBL12)  +THETA*(DT*NLN12+LHNN12+DT*THETA*LHNLN12)
            BIGM(COUNT+IBL13)=BIGM(COUNT+IBL13)  +THETA*(DT*NLN13+LHNN13+DT*THETA*LHNLN13)
              
            BIGM(COUNT+IBL21)=BIGM(COUNT+IBL21)  +THETA*(DT*NLN21+LHNN21+DT*THETA*LHNLN21)
            BIGM(COUNT+IBL22)=BIGM(COUNT+IBL22)  +THETA*(DT*NLN22+LHNN22+DT*THETA*LHNLN22)
            BIGM(COUNT+IBL23)=BIGM(COUNT+IBL23)  +THETA*(DT*NLN23+LHNN23+DT*THETA*LHNLN23)
              
            BIGM(COUNT+IBL31)=BIGM(COUNT+IBL31)  +THETA*(DT*NLN31+LHNN31+DT*THETA*LHNLN31)
            BIGM(COUNT+IBL32)=BIGM(COUNT+IBL32)  +THETA*(DT*NLN32+LHNN32+DT*THETA*LHNLN32)
            BIGM(COUNT+IBL33)=BIGM(COUNT+IBL33)  +THETA*(DT*NLN33+LHNN33+DT*THETA*LHNLN33)
         
! Lump the absorption to put in to pressure. 
            ML(GLOBI)=ML(GLOBI) + VLM 
! put diagonal of viscosity in here...
            ML(GLOBI)=ML(GLOBI) + VLKTEN*NLOCIDENT(ILOC,JLOC)
             
            IF(ABSLUM) THEN 
               ML(GLOBI)=ML(GLOBI)+DT*VLM*(XABSOR(GLOBI)+YABSOR(GLOBI)+ZABSOR(GLOBI))/3.
               BIGM(ICENT+IBL11)=BIGM(ICENT+IBL11)  +DT*VLM*XABSOR(GLOBI)
               BIGM(ICENT+IBL22)=BIGM(ICENT+IBL22)  +DT*VLM*YABSOR(GLOBI)
               BIGM(ICENT+IBL33)=BIGM(ICENT+IBL33)  +DT*VLM*ZABSOR(GLOBI)
               VECX(GLOBI)=VECX(GLOBI) - VLM*XABSOR(GLOBI)*U(GLOBI)
               VECY(GLOBI)=VECY(GLOBI) - VLM*YABSOR(GLOBI)*V(GLOBI)
               VECZ(GLOBI)=VECZ(GLOBI) - VLM*ZABSOR(GLOBI)*W(GLOBI)
            ENDIF
             
         end do ! Was loop 360

      end do ! Was loop 350
   end do ! Was loop 340

        IF(NBUOY.NE.0) THEN
          DO I=1,NONODS
            BOY_ML(I)         =BOY_ML(I)         +ml(I)
            BOY_ML(I+4*NONODS)=BOY_ML(I+4*NONODS)+ml(I)
            BOY_ML(I+8*NONODS)=BOY_ML(I+8*NONODS)+ml(I)
          END DO
        ENDIF

   ewrite(2,*) 'NONODS,XNONOD:',NONODS,XNONOD
   ewrite(2,*) 'abslum:',abslum
   CALL PMINMX(ML,NONODS,'******ML  ')
   CALL PMINMX(XABSOR,NONODS,'******XABSOR  ')
   CALL PMINMX(DENPT,NONODS,'******DENPT  ')
   CALL PMINMX(sourcx,NONODS,'******sourcx  ')
   CALL PMINMX(sourcy,NONODS,'******sourcy  ')
   CALL PMINMX(sourcz,NONODS,'******sourcz  ')
   CALL PMINMX(MUPTXX,nonods,'******MUPTXX  ')
   CALL PMINMX(MUPTXY,nonods,'******MUPTXY  ')
   CALL PMINMX(MUPTXZ,nonods,'******MUPTXZ  ')
   CALL PMINMX(MUPTYY,nonods,'******MUPTYY  ')
   CALL PMINMX(MUPTYZ,nonods,'******MUPTYZ  ')
   CALL PMINMX(MUPTZZ,nonods,'******MUPTZZ  ')

 END SUBROUTINE DIFF3DSIMP
 
 SUBROUTINE LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             MXXDTW,MXYDTW,MXZDTW, &
     &             MYYDTW,MYZDTW,MZZDTW) 
   INTEGER NONODS,TOTELE,NLOC,NGI, ELE,GI
   INTEGER VONDGL(TOTELE*NLOC)
   REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
   REAL NU(NONODS),NV(NONODS),NW(NONODS)
   REAL UG(NONODS),VG(NONODS),WG(NONODS)
   REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
   REAL MYZDTW(NGI),MZZDTW(NGI)
! Local variables...
   REAL SXXGI,SXYGI,SXZGI,SYYGI,SYZGI,SZZGI,SYXGI,SZXGI,SZYGI
   REAL VIS,CS2,FOURCS
   REAL VELUD,VELVD,VELWD
   INTEGER L,IGLV
! WORK OUT STRAIN RATES
   !ewrite(1,*) 'in subroutine lesvis' 
   SXXGI=0.
   SXYGI=0.
   SXZGI=0.
   SYYGI=0.
   SYZGI=0.
   SZZGI=0.
   do L=1,NLOC
      IGLV=VONDGL((ELE-1)*NLOC+L)
      VELUD=NU(IGLV)-UG(IGLV)
      VELVD=NV(IGLV)-VG(IGLV)
      VELWD=NW(IGLV)-WG(IGLV)
      SXXGI = SXXGI + NX(L,GI)*VELUD
      SXYGI = SXYGI + 0.5*(NY(L,GI)*VELUD+NX(L,GI)*VELVD)
      SXZGI = SXZGI + 0.5*(NZ(L,GI)*VELUD+NX(L,GI)*VELWD)
      SYYGI = SYYGI + NY(L,GI)*VELVD
      SYZGI = SYZGI + 0.5*(NZ(L,GI)*VELVD+NY(L,GI)*VELWD)
      SZZGI = SZZGI + NZ(L,GI)*VELWD
   end do
   SYXGI=SXYGI
   SZXGI=SXZGI
   SZYGI=SYZGI

   VIS=SQRT(2.* (SXXGI*SXXGI + SXYGI*SXYGI + SXZGI*SXZGI&
     &                +  SYXGI*SYXGI + SYYGI*SYYGI + SYZGI*SYZGI&
     &                +  SZXGI*SZXGI + SZYGI*SZYGI + SZZGI*SZZGI))

! THEN FIND TURBULENT 'VISCOSITIES'
   CS2=0.1**2

   FOURCS=4.*CS2

   MXXDTW(GI) = FOURCS*VIS*MXXDTW(GI)
   MXYDTW(GI) = FOURCS*VIS*MXYDTW(GI)
   MXZDTW(GI) = FOURCS*VIS*MXZDTW(GI)
   MYYDTW(GI) = FOURCS*VIS*MYYDTW(GI)
   MYZDTW(GI) = FOURCS*VIS*MYZDTW(GI)
   MZZDTW(GI) = FOURCS*VIS*MZZDTW(GI)

 END SUBROUTINE LESVIS




 SUBROUTINE DG_LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             DGNU,DGNV,DGNW, UG,VG,WG, &
     &             MXXDTW,MXYDTW,MXZDTW, &
     &             MYYDTW,MYZDTW,MZZDTW, &
     &             JUSTH,DIRX,DIRY,DIRZ) 
   INTEGER NONODS,TOTELE,NLOC,NGI, ELE,GI
   INTEGER VONDGL(TOTELE*NLOC)
   REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
   REAL DGNU(TOTELE*NLOC),DGNV(TOTELE*NLOC),DGNW(TOTELE*NLOC)
   REAL UG(NONODS),VG(NONODS),WG(NONODS)
   REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
   REAL MYZDTW(NGI),MZZDTW(NGI)
         LOGICAL JUSTH
         REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
! If JUSTH then do not have vertical derivatives in calculating LES viscocity. 
! Local variables...
   REAL SXXGI,SXYGI,SXZGI,SYYGI,SYZGI,SZZGI,SYXGI,SZXGI,SZYGI
   REAL RNX,RNY,RNZ
   REAL VIS,CS2,FOURCS
   REAL VELUD,VELVD,VELWD
   INTEGER L,IGLV
! WORK OUT STRAIN RATES
   ewrite(1,*) 'in subroutine dg_lesvis' 
           SXXGI=0.
           SXYGI=0.
           SXZGI=0.
           SYYGI=0.
           SYZGI=0.
           SZZGI=0.
           DO 374 L=1,NLOC
           IGLV=VONDGL((ELE-1)*NLOC+L)
           IF(JUSTH) THEN
             RNX=NX(L,GI)*(1.-ABS(DIRX(GI)))
             RNY=NY(L,GI)*(1.-ABS(DIRY(GI)))
             RNZ=NZ(L,GI)*(1.-ABS(DIRZ(GI)))
           ELSE
             RNX=NX(L,GI)
             RNY=NY(L,GI)
             RNZ=NZ(L,GI)
           ENDIF
!             VELUD=DGNU((ELE-1)*NLOC+L)-UG(IGLV)
!             VELVD=DGNV((ELE-1)*NLOC+L)-VG(IGLV)
!             VELWD=DGNW((ELE-1)*NLOC+L)-WG(IGLV)
             VELUD=DGNU((ELE-1)*NLOC+L)
             VELVD=DGNV((ELE-1)*NLOC+L)
             VELWD=DGNW((ELE-1)*NLOC+L)
             SXXGI = SXXGI + RNX*VELUD
             SXYGI = SXYGI + 0.5*(RNY*VELUD+RNX*VELVD)
             SXZGI = SXZGI + 0.5*(RNZ*VELUD+RNX*VELWD)
             SYYGI = SYYGI + RNY*VELVD
             SYZGI = SYZGI + 0.5*(RNZ*VELVD+RNY*VELWD)
             SZZGI = SZZGI + RNZ*VELWD
374        CONTINUE
           SYXGI=SXYGI
           SZXGI=SXZGI
           SZYGI=SYZGI
!
           VIS=SQRT(2.* (SXXGI*SXXGI + SXYGI*SXYGI + SXZGI*SXZGI &
     &                +  SYXGI*SYXGI + SYYGI*SYYGI + SYZGI*SYZGI &
     &                +  SZXGI*SZXGI + SZYGI*SZYGI + SZZGI*SZZGI))
!
! THEN FIND TURBULENT 'VISCOSITIES'
           CS2=0.1**2
!           CS2=0.2**2
!           CS2=0.4**2
!           CS2=0.8**2
           FOURCS=10.*4.*CS2
!
         IF(JUSTH) THEN
           MXXDTW(GI) = FOURCS*VIS*MXXDTW(GI)*(1.-ABS(DIRX(GI)*DIRX(GI)))
           MXYDTW(GI) = FOURCS*VIS*MXYDTW(GI)*(1.-ABS(DIRX(GI)*DIRY(GI)))
           MXZDTW(GI) = FOURCS*VIS*MXZDTW(GI)*(1.-ABS(DIRX(GI)*DIRZ(GI)))
           MYYDTW(GI) = FOURCS*VIS*MYYDTW(GI)*(1.-ABS(DIRY(GI)*DIRY(GI)))
           MYZDTW(GI) = FOURCS*VIS*MYZDTW(GI)*(1.-ABS(DIRY(GI)*DIRZ(GI)))
           MZZDTW(GI) = FOURCS*VIS*MZZDTW(GI)*(1.-ABS(DIRZ(GI)*DIRZ(GI)))
         ELSE
           MXXDTW(GI) = FOURCS*VIS*MXXDTW(GI)
           MXYDTW(GI) = FOURCS*VIS*MXYDTW(GI)
           MXZDTW(GI) = FOURCS*VIS*MXZDTW(GI)
           MYYDTW(GI) = FOURCS*VIS*MYYDTW(GI)
           MYZDTW(GI) = FOURCS*VIS*MYZDTW(GI)
           MZZDTW(GI) = FOURCS*VIS*MZZDTW(GI)
         ENDIF
         
 END SUBROUTINE DG_LESVIS

 SUBROUTINE DG_DG_LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             N,NX,NY,NZ, &
     &             DGNU,DGNV,DGNW,&
     &             MXXDTW,MXYDTW,MXZDTW, &
     &             MYYDTW,MYZDTW,MZZDTW,&
     &             LOCDGNUDX,LOCDGNUDY,LOCDGNUDZ,&
     &             LOCDGNVDX,LOCDGNVDY,LOCDGNVDZ,&
     &             LOCDGNWDX,LOCDGNWDY,LOCDGNWDZ,&
     &             JUSTH,DIRX,DIRY,DIRZ)
   INTEGER NONODS,TOTELE,NLOC,NGI, ELE,GI
   INTEGER VONDGL(TOTELE*NLOC)
   REAL N(NLOC,NGI),NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
   REAL DGNU(TOTELE*NLOC),DGNV(TOTELE*NLOC),DGNW(TOTELE*NLOC)
   REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
   REAL MYZDTW(NGI),MZZDTW(NGI)
   REAL LOCDGNUDX(NLOC),LOCDGNUDY(NLOC),LOCDGNUDZ(NLOC)
   REAL LOCDGNVDX(NLOC),LOCDGNVDY(NLOC),LOCDGNVDZ(NLOC)
   REAL LOCDGNWDX(NLOC),LOCDGNWDY(NLOC),LOCDGNWDZ(NLOC)
   LOGICAL JUSTH
   REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
! Local variables...
   REAL SXXGI,SXYGI,SXZGI,SYYGI,SYZGI,SZZGI,SYXGI,SZXGI,SZYGI
   REAL RLOCX,RLOCY,RLOCZ
   REAL VIS,CS2,FOURCS
   INTEGER L,IGLV
! WORK OUT STRAIN RATES
!   ewrite(1,*) 'in subroutine dg_dg_lesvis' 
   SXXGI=0.
   SXYGI=0.
   SXZGI=0.
   SYYGI=0.
   SYZGI=0.
   SZZGI=0.
   do L=1,NLOC
!      SXXGI = SXXGI + N(L,GI)*LOCDGNUDX(L)
!      SXYGI = SXYGI + 0.5*(N(L,GI)*LOCDGNUDY(L)+N(L,GI)*LOCDGNVDX(L))
!      SXZGI = SXZGI + 0.5*(N(L,GI)*LOCDGNUDZ(L)+N(L,GI)*LOCDGNWDX(L))
!      SYYGI = SYYGI + N(L,GI)*LOCDGNVDY(L)
!      SYZGI = SYZGI + 0.5*(N(L,GI)*LOCDGNVDZ(L)+N(L,GI)*LOCDGNWDY(L))
!      SZZGI = SZZGI + N(L,GI)*LOCDGNWDZ(L)
           IF(JUSTH) THEN
             RLOCX=1.-ABS(DIRX(GI))
             RLOCY=1.-ABS(DIRY(GI))
             RLOCZ=1.-ABS(DIRZ(GI))
           ELSE
             RLOCX=1.0
             RLOCY=1.0
             RLOCZ=1.0
           ENDIF

             SXXGI = SXXGI + N(L,GI)*LOCDGNUDX(L)*RLOCX
             SXYGI = SXYGI + 0.5*(N(L,GI)*LOCDGNUDY(L)*RLOCY+N(L,GI)*LOCDGNVDX(L)*RLOCX)
             SXZGI = SXZGI + 0.5*(N(L,GI)*LOCDGNUDZ(L)*RLOCZ+N(L,GI)*LOCDGNWDX(L)*RLOCX)
             SYYGI = SYYGI + N(L,GI)*LOCDGNVDY(L)*RLOCY
             SYZGI = SYZGI + 0.5*(N(L,GI)*LOCDGNVDZ(L)*RLOCZ+N(L,GI)*LOCDGNWDY(L)*RLOCY)
             SZZGI = SZZGI + N(L,GI)*LOCDGNWDZ(L)*RLOCZ
   end do
   SYXGI=SXYGI
   SZXGI=SXZGI
   SZYGI=SYZGI

   VIS=SQRT(2.* (SXXGI*SXXGI + SXYGI*SXYGI + SXZGI*SXZGI&
     &                +  SYXGI*SYXGI + SYYGI*SYYGI + SYZGI*SYZGI&
     &                +  SZXGI*SZXGI + SZYGI*SZYGI + SZZGI*SZZGI))

! THEN FIND TURBULENT 'VISCOSITIES'
   CS2=0.1**2
   FOURCS=10.*4.*CS2

         IF(JUSTH) THEN
           MXXDTW(GI) = FOURCS*VIS*MXXDTW(GI)*(1.-ABS(DIRX(GI)*DIRX(GI)))
           MXYDTW(GI) = FOURCS*VIS*MXYDTW(GI)*(1.-ABS(DIRX(GI)*DIRY(GI)))
           MXZDTW(GI) = FOURCS*VIS*MXZDTW(GI)*(1.-ABS(DIRX(GI)*DIRZ(GI)))
           MYYDTW(GI) = FOURCS*VIS*MYYDTW(GI)*(1.-ABS(DIRY(GI)*DIRY(GI)))
           MYZDTW(GI) = FOURCS*VIS*MYZDTW(GI)*(1.-ABS(DIRY(GI)*DIRZ(GI)))
           MZZDTW(GI) = FOURCS*VIS*MZZDTW(GI)*(1.-ABS(DIRZ(GI)*DIRZ(GI)))
         ELSE
           MXXDTW(GI) = FOURCS*VIS*MXXDTW(GI)
           MXYDTW(GI) = FOURCS*VIS*MXYDTW(GI)
           MXZDTW(GI) = FOURCS*VIS*MXZDTW(GI)
           MYYDTW(GI) = FOURCS*VIS*MYYDTW(GI)
           MYZDTW(GI) = FOURCS*VIS*MYZDTW(GI)
           MZZDTW(GI) = FOURCS*VIS*MZZDTW(GI)
         ENDIF 

 END SUBROUTINE DG_DG_LESVIS
 
 
 

 SUBROUTINE HALF_LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             HIMXXDTW,HIMXYDTW,HIMXZDTW, &
     &             HIMYYDTW,HIMYZDTW,HIMZZDTW,&
     &             MXXDTW,MXYDTW,MXZDTW, &
     &             MYYDTW,MYZDTW,MZZDTW, &
     &             JUSTH,DIRX,DIRY,DIRZ)
   INTEGER NONODS,TOTELE,NLOC,NGI, ELE,GI
   INTEGER VONDGL(TOTELE*NLOC)
   REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
   REAL NU(NONODS),NV(NONODS),NW(NONODS)
   REAL UG(NONODS),VG(NONODS),WG(NONODS)
   REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
   REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
   REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
   REAL MYZDTW(NGI),MZZDTW(NGI)
   LOGICAL JUSTH
   REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
! If JUSTH then do not include vertical derivatives in calculating viscocity.
! Local variables...
   REAL SXXGI,SXYGI,SXZGI,SYYGI,SYZGI,SZZGI,SYXGI,SZXGI,SZYGI
   REAL VIS,CS2,FOURCS
   REAL VELUD,VELVD,VELWD,PROD,SQRTPROD
   REAL RNX,RNY,RNZ
   INTEGER L,IGLV
! WORK OUT STRAIN RATES
   !ewrite(1,*) 'in subroutine half_lesvis' 
   SXXGI=0.
   SXYGI=0.
   SXZGI=0.
   SYYGI=0.
   SYZGI=0.
   SZZGI=0.
   do  L=1,NLOC
           IGLV=VONDGL((ELE-1)*NLOC+L)
           IF(JUSTH) THEN
             RNX=NX(L,GI)*(1.-ABS(DIRX(GI)))
             RNY=NY(L,GI)*(1.-ABS(DIRY(GI)))
             RNZ=NZ(L,GI)*(1.-ABS(DIRZ(GI)))
           ELSE
             RNX=NX(L,GI)
             RNY=NY(L,GI)
             RNZ=NZ(L,GI)
           ENDIF
             VELUD=NU(IGLV)-UG(IGLV)
             VELVD=NV(IGLV)-VG(IGLV)
             VELWD=NW(IGLV)-WG(IGLV)
             SXXGI = SXXGI + RNX*VELUD
             SXYGI = SXYGI + 0.5*(RNY*VELUD+RNX*VELVD)
             SXZGI = SXZGI + 0.5*(RNZ*VELUD+RNX*VELWD)
             SYYGI = SYYGI + RNY*VELVD
             SYZGI = SYZGI + 0.5*(RNZ*VELVD+RNY*VELWD)
             SZZGI = SZZGI + RNZ*VELWD
   end do
   SYXGI=SXYGI
   SZXGI=SXZGI
   SZYGI=SYZGI

   VIS=SQRT(2.* (SXXGI*SXXGI + SXYGI*SXYGI + SXZGI*SXZGI&
     &                +  SYXGI*SYXGI + SYYGI*SYYGI + SYZGI*SYZGI&
     &                +  SZXGI*SZXGI + SZYGI*SZYGI + SZZGI*SZZGI))

! THEN FIND TURBULENT 'VISCOSITIES'
   CS2=0.1**2

   FOURCS=10.*4.*CS2 
           
   PROD=FOURCS*VIS
   SQRTPROD=SQRT(PROD)

   IF(JUSTH) THEN
     HIMXXDTW(GI) = SQRTPROD*HIMXXDTW(GI)*(1.-ABS(DIRX(GI)*DIRX(GI)))
     HIMXYDTW(GI) = SQRTPROD*HIMXYDTW(GI)*(1.-ABS(DIRX(GI)*DIRY(GI)))
     HIMXZDTW(GI) = SQRTPROD*HIMXZDTW(GI)*(1.-ABS(DIRX(GI)*DIRZ(GI)))
     HIMYYDTW(GI) = SQRTPROD*HIMYYDTW(GI)*(1.-ABS(DIRY(GI)*DIRY(GI)))
     HIMYZDTW(GI) = SQRTPROD*HIMYZDTW(GI)*(1.-ABS(DIRY(GI)*DIRZ(GI)))
     HIMZZDTW(GI) = SQRTPROD*HIMZZDTW(GI)*(1.-ABS(DIRZ(GI)*DIRZ(GI)))

     MXXDTW(GI) = PROD*MXXDTW(GI)*(1.-ABS(DIRX(GI)*DIRX(GI)))
     MXYDTW(GI) = PROD*MXYDTW(GI)*(1.-ABS(DIRX(GI)*DIRY(GI)))
     MXZDTW(GI) = PROD*MXZDTW(GI)*(1.-ABS(DIRX(GI)*DIRZ(GI)))
     MYYDTW(GI) = PROD*MYYDTW(GI)*(1.-ABS(DIRY(GI)*DIRY(GI)))
     MYZDTW(GI) = PROD*MYZDTW(GI)*(1.-ABS(DIRY(GI)*DIRZ(GI)))
     MZZDTW(GI) = PROD*MZZDTW(GI)*(1.-ABS(DIRZ(GI)*DIRZ(GI)))
   ELSE
     HIMXXDTW(GI) = SQRTPROD*HIMXXDTW(GI)
     HIMXYDTW(GI) = SQRTPROD*HIMXYDTW(GI)
     HIMXZDTW(GI) = SQRTPROD*HIMXZDTW(GI)
     HIMYYDTW(GI) = SQRTPROD*HIMYYDTW(GI)
     HIMYZDTW(GI) = SQRTPROD*HIMYZDTW(GI)
     HIMZZDTW(GI) = SQRTPROD*HIMZZDTW(GI)

     MXXDTW(GI) = PROD*MXXDTW(GI)
     MXYDTW(GI) = PROD*MXYDTW(GI)
     MXZDTW(GI) = PROD*MXZDTW(GI)
     MYYDTW(GI) = PROD*MYYDTW(GI)
     MYZDTW(GI) = PROD*MYZDTW(GI)
     MZZDTW(GI) = PROD*MZZDTW(GI)
   ENDIF

 END SUBROUTINE HALF_LESVIS
 
 
 

 SUBROUTINE HALF_LESVIS_HART(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             HIMXXDTW,HIMXYDTW,HIMXZDTW, &
     &             HIMYYDTW,HIMYZDTW,HIMZZDTW,&
     &             MXXDTW,MXYDTW,MXZDTW, &
     &             MYYDTW,MYZDTW,MZZDTW) 
   INTEGER NONODS,TOTELE,NLOC,NGI, ELE,GI
   INTEGER VONDGL(TOTELE*NLOC)
   REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
   REAL NU(NONODS),NV(NONODS),NW(NONODS)
   REAL UG(NONODS),VG(NONODS),WG(NONODS)
   REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
   REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
   REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
   REAL MYZDTW(NGI),MZZDTW(NGI)
! Local variables...
   REAL SXXGI,SXYGI,SXZGI,SYYGI,SYZGI,SZZGI,SYXGI,SZXGI,SZYGI
   REAL VIS,CS2,FOURCS
   REAL VELUD,VELVD,VELWD,PROD,SQRTPROD
   INTEGER L,IGLV
! WORK OUT STRAIN RATES
   ewrite(1,*) 'in subroutine half_lesvis_hart' 
   SXXGI=0.
   SXYGI=0.
   SXZGI=0.
   SYYGI=0.
   SYZGI=0.
   SZZGI=0.
   do L=1,NLOC
      IGLV=VONDGL((ELE-1)*NLOC+L)
      VELUD=NU(IGLV)-UG(IGLV)
      VELVD=NV(IGLV)-VG(IGLV)
      VELWD=NW(IGLV)-WG(IGLV)
      SXXGI = SXXGI + NX(L,GI)*VELUD
      SXYGI = SXYGI + 0.5*(NY(L,GI)*VELUD+NX(L,GI)*VELVD)
      SXZGI = SXZGI + 0.5*(NZ(L,GI)*VELUD+NX(L,GI)*VELWD)
      SYYGI = SYYGI + NY(L,GI)*VELVD
      SYZGI = SYZGI + 0.5*(NZ(L,GI)*VELVD+NY(L,GI)*VELWD)
      SZZGI = SZZGI + NZ(L,GI)*VELWD
   end do
   SYXGI=SXYGI
   SZXGI=SXZGI
   SZYGI=SYZGI

   VIS=SQRT(2.* (SXXGI*SXXGI + SXYGI*SXYGI + SXZGI*SXZGI&
     &                +  SYXGI*SYXGI + SYYGI*SYYGI + SYZGI*SYZGI&
     &                +  SZXGI*SZXGI + SZYGI*SZYGI + SZZGI*SZZGI))

! THEN FIND TURBULENT 'VISCOSITIES'
   CS2=0.1**2
   FOURCS=1.*4.*CS2 

   PROD=FOURCS*VIS
   SQRTPROD=SQRT(PROD)

   HIMXXDTW(GI) = SQRTPROD*HIMXXDTW(GI)
   HIMXYDTW(GI) = SQRTPROD*HIMXYDTW(GI)
   HIMXZDTW(GI) = SQRTPROD*HIMXZDTW(GI)
   HIMYYDTW(GI) = SQRTPROD*HIMYYDTW(GI)
   HIMYZDTW(GI) = SQRTPROD*HIMYZDTW(GI)
   HIMZZDTW(GI) = SQRTPROD*HIMZZDTW(GI)

   MXXDTW(GI) = PROD*MXXDTW(GI)
   MXYDTW(GI) = PROD*MXYDTW(GI)
   MXZDTW(GI) = PROD*MXZDTW(GI)
   MYYDTW(GI) = PROD*MYYDTW(GI)
   MYZDTW(GI) = PROD*MYZDTW(GI)
   MZZDTW(GI) = PROD*MZZDTW(GI)
 END SUBROUTINE HALF_LESVIS_HART

 SUBROUTINE LESVIS2(NONODS,TOTELE,NLOC,VONDGL,NGI,ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             E) 
! This subroutine multiplies the length-scale matrix E by 4C_s^2 * |S|
! to give the viscosity matrix, necessary for SFS2, the 
! anisotropic, inhomogeneous SFS model which solves the filtered NS
! equations in stress form.

   INTEGER NONODS,TOTELE,NLOC,NGI,ELE,GI
   INTEGER VONDGL(TOTELE*NLOC)
   REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
   REAL NU(NONODS),NV(NONODS),NW(NONODS)
   REAL UG(NONODS),VG(NONODS),WG(NONODS)
   REAL E(6,6,NGI)
! Local variables...
   REAL SXXGI,SXYGI,SXZGI,SYYGI,SYZGI,SZZGI,SYXGI,SZXGI,SZYGI
   REAL VIS,CS2,FOURCS
   REAL VELUD,VELVD,VELWD
   INTEGER L,IGLV,IDIM,JDIM
! WORK OUT STRAIN RATES
   ewrite(1,*) 'in subroutine lesvis2' 
   SXXGI=0.
   SXYGI=0.
   SXZGI=0.
   SYYGI=0.
   SYZGI=0.
   SZZGI=0.
   do L=1,NLOC
      IGLV=VONDGL((ELE-1)*NLOC+L)
      VELUD=NU(IGLV)-UG(IGLV)
      VELVD=NV(IGLV)-VG(IGLV)
      VELWD=NW(IGLV)-WG(IGLV)
      SXXGI = SXXGI + NX(L,GI)*VELUD
      SXYGI = SXYGI + 0.5*(NY(L,GI)*VELUD+NX(L,GI)*VELVD)
      SXZGI = SXZGI + 0.5*(NZ(L,GI)*VELUD+NX(L,GI)*VELWD)
      SYYGI = SYYGI + NY(L,GI)*VELVD
      SYZGI = SYZGI + 0.5*(NZ(L,GI)*VELVD+NY(L,GI)*VELWD)
      SZZGI = SZZGI + NZ(L,GI)*VELWD
   end do
   SYXGI=SXYGI
   SZXGI=SXZGI
   SZYGI=SYZGI

   VIS=SQRT(2.* (SXXGI*SXXGI + SXYGI*SXYGI + SXZGI*SXZGI&
     &                +  SYXGI*SYXGI + SYYGI*SYYGI + SYZGI*SYZGI&
     &                +  SZXGI*SZXGI + SZYGI*SZYGI + SZZGI*SZZGI))

! THEN FIND TURBULENT 'VISCOSITIES'
   CS2=0.1**2
   FOURCS=4.*CS2

! Put a bit in here which multiplies E by FOURCS*VIS           
   do IDIM=1,6
      do JDIM=1,6
         E(IDIM,JDIM,GI) = FOURCS*VIS*E(IDIM,JDIM,GI)
      END DO
   END DO
 END SUBROUTINE LESVIS2

 SUBROUTINE LESELA(NONODS,&
     &             MUGIXX,MUGIXY,MUGIXZ,MUGIYY,MUGIYZ,MUGIZZ, &
     &             ELE,GI,NGI,NLOC,E)

! This subroutine calculates the length-scale (material 
! properties) matrix in local co-ordinates, then rotates 
! it, in order to solve the Navier-Stokes equations in 
! stress form. It returns E - the six-by-six matrix of
! material properties which is multiplied by the stress 
! vector in DIFF3D.
! MUGIXX etc contain the matrix at the Gauss points...
! Te is the rotation matrix.

   INTEGER NONODS,NDIM,NLOC
   REAL TOLER
   PARAMETER(NDIM=3,TOLER=1.E-14)
   INTEGER GI,NGI,ELE
   REAL MUGIXX,MUGIXY,MUGIXZ,MUGIYY,MUGIYZ,MUGIZZ
   REAL E(6,6,NGI)
! Local variables...
   INTEGER NOD,IDIM,JDIM,NDIM2,L,IGLV
   INTEGER I,J,M,NN,IDSIP
   REAL V(3,3),D(3),A(3,3)
   REAL Te(6,6),EDASH(6,6),AA(6,6)
   REAL FERRGI(3,3)
   REAL RSUM

   ewrite(1,*) 'in subroutine lesela' 

   NDIM2=NDIM**2

! Calculate the eigen values and vectors at the Gauss pt GI..
   FERRGI(1,1)=MUGIXX
   FERRGI(1,2)=MUGIXY
   FERRGI(1,3)=MUGIXZ

   FERRGI(2,1)=MUGIXY
   FERRGI(2,2)=MUGIYY
   FERRGI(2,3)=MUGIYZ

   FERRGI(3,1)=MUGIXZ
   FERRGI(3,2)=MUGIYZ
   FERRGI(3,3)=MUGIZZ

! This sub performs Jacobi rotations of a symmetric matrix in order to 
! find the eigen-vectors V and the eigen values D so 

   CALL JACDIA(FERRGI,V,D,NDIM,A)

   do IDIM=1,NDIM
      D(IDIM)=1./SQRT(MAX(TOLER,D(IDIM)))
   END DO
! Then find rotation matirx Te

! First do top-left quarter 
   do IDIM=1,NDIM
      do JDIM=1,NDIM
         Te(IDIM,JDIM)=V(IDIM,JDIM)*V(IDIM,JDIM)
      END DO
   END DO

! Then top-right quarter
   do IDIM=1,NDIM
      do JDIM=1,2
         M =JDIM+1
         NN=JDIM+3
         Te(IDIM,NN)=V(IDIM,JDIM)*V(IDIM,M)
      END DO
      Te(IDIM,6)=V(IDIM,3)*V(IDIM,1)
   END DO

! Then bottom-left quarter
   do JDIM=1,NDIM
      do IDIM=1,2                
         M =IDIM+1
         NN=IDIM+3
         Te(NN,JDIM)=2*V(IDIM,JDIM)*V(M,JDIM)
      END DO
      Te(6,JDIM)=2*V(3,JDIM)*V(1,JDIM)
   END DO

! Then bottom-right quarter
   Te(4,4)=V(1,1)*V(2,2)+V(1,2)*V(2,1)
   Te(4,5)=V(1,2)*V(2,3)+V(1,3)*V(2,2)
   Te(4,6)=V(1,3)*V(2,1)+V(1,1)*V(2,3)
   Te(5,4)=V(2,1)*V(3,2)+V(2,2)*V(3,1)
   Te(5,5)=V(2,2)*V(3,3)+V(2,3)*V(3,2)
   Te(5,6)=V(2,3)*V(3,1)+V(2,1)*V(3,3)
   Te(6,4)=V(3,1)*V(1,2)+V(3,2)*V(1,1)
   Te(6,5)=V(3,2)*V(1,3)+V(3,3)*V(1,2)
   Te(6,6)=V(3,3)*V(1,1)+V(3,1)*V(1,3)

! Next find length-scale matrix E' or EDASH
! set matrix to zero initially...
   EDASH = 0.0

! First do top-left quarter
   do IDIM=1,NDIM
      do JDIM=1,NDIM
         IF(IDIM.EQ.JDIM) THEN
            EDASH(IDIM,JDIM)=(4/3)*D(IDIM)*D(JDIM)
         ELSE
            EDASH(IDIM,JDIM)=-(2/3)*D(IDIM)*D(JDIM)
         END IF
      END DO
   END DO

! Top-right and bottom-left quarters (all zeros)

! Then do bottom-right quarter
   EDASH(4,4)=D(1)*D(2)
   EDASH(5,5)=D(1)*D(3)
   EDASH(6,6)=D(2)*D(3)
! Finally, multiply up to get [E] = [Te]T [E'] [Te]

! First do [AA] = [E'] [Te]
   do IDIM=1,6
      do JDIM=1,6
         RSUM=0.
         do M=1,6
            RSUM=RSUM+EDASH(IDIM,M)*Te(M,JDIM)
         END DO
         AA(IDIM,JDIM)=RSUM
      END DO
   END DO

! Then do [E] = [Te]T [AA]
   do IDIM=1,6
      do JDIM=1,6
         RSUM=0.
         do M=1,6
            RSUM=RSUM+Te(M,IDIM)*AA(M,JDIM)
         END DO
         E(IDIM,JDIM,GI)=RSUM
      END DO
   END DO
 END SUBROUTINE LESELA

 SUBROUTINE LESDERIVS(&
     &          VEDUDX,VEDUDY,VEDUDZ,&
     &          VEDVDX,VEDVDY,VEDVDZ,&
     &          VEDWDX,VEDWDY,VEDWDZ,&
     &          ML,&
     &          NU,NV,NW,&
     &          X,Y,Z, XONDGL, TOTELE,NONODS,NLOC,NGI, &
     &          N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA, D3,DCYL,&
     &          NX,NY,NZ,&
     &          NDGLNO, nnodp,para,halo_tag,WORK,WORK2,&
     &          LESNOD,EV,&
     &          MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ, &
     &          UG,VG,WG, &
     &          VONDGL)
! calculate the derivatives of NU,NV,NW wrt X,Y,Z
! the derivatives are the 9 variables listed below...
   INTEGER NONODS
   REAL VEDUDX(NONODS),VEDUDY(NONODS),VEDUDZ(NONODS)
   REAL VEDVDX(NONODS),VEDVDY(NONODS),VEDVDZ(NONODS)
   REAL VEDWDX(NONODS),VEDWDY(NONODS),VEDWDZ(NONODS)
   REAL ML(NONODS)
   REAL NU(NONODS),NV(NONODS),NW(NONODS)
   REAL X(NONODS),Y(NONODS),Z(NONODS)
   INTEGER TOTELE,NLOC,NGI
   LOGICAL D3,DCYL
   INTEGER XONDGL(TOTELE*NLOC),NDGLNO(TOTELE*NLOC)
   REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
   REAL WEIGHT(NGI)
   REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
   REAL DETWEI(NGI),RA(NGI)
   INTEGER nnodp,para,halo_tag
   REAL WORK(NONODS),WORK2(NONODS)
   INTEGER LESNOD
   REAL EV(6,6,LESNOD*NONODS)
   REAL MUPTXX(NONODS),MUPTXY(NONODS),MUPTXZ(NONODS)
   REAL MUPTYY(NONODS),MUPTYZ(NONODS),MUPTZZ(NONODS)
   REAL UG(NONODS),VG(NONODS),WG(NONODS)
   INTEGER VONDGL(TOTELE*NLOC)
   ! Local variables...
   REAL VOLUME
   INTEGER IGL,JGL,ILOC,JLOC,GI,ELE,IDIM,JDIM
   REAL EET(6,6,1)
   REAL SXXGI,SXYGI,SXZGI,SYYGI,SYZGI,SZZGI,SYXGI,SZXGI,SZYGI
   REAL VIS,CS2,FOURCS
   ewrite(1,*) 'in subroutine lesderivs' 
     
   ML(1:NONODS) = 0.0
   VEDUDX(1:NONODS) = 0.0
   VEDUDY(1:NONODS) = 0.0
   VEDUDZ(1:NONODS) = 0.0
   VEDVDX(1:NONODS) = 0.0
   VEDVDY(1:NONODS) = 0.0
   VEDVDZ(1:NONODS) = 0.0
   VEDWDX(1:NONODS) = 0.0
   VEDWDY(1:NONODS) = 0.0
   VEDWDZ(1:NONODS) = 0.0
       
   do  ELE=1,TOTELE
! Calculate DETWEI,RA,NX,NY,NZ for element ELE
      CALL DETNLXR(ELE, X,Y,Z, XONDGL, TOTELE,NONODS,NLOC,NGI, &
     &                 N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, D3,DCYL,&
     &                 NX,NY,NZ)

      do ILOC=1,NLOC
         IGL=NDGLNO((ELE-1)*NLOC+ILOC)
         do JLOC=1,NLOC
            JGL=NDGLNO((ELE-1)*NLOC+JLOC)
            do GI=1,NGI
               ML(IGL)=ML(IGL)+DETWEI(GI)*N(ILOC,GI)*N(JLOC,GI)
               VEDUDX(IGL)=VEDUDX(IGL)+DETWEI(GI)*N(ILOC,GI)*NX(JLOC,GI)*NU(JGL)
               VEDUDY(IGL)=VEDUDY(IGL)+DETWEI(GI)*N(ILOC,GI)*NY(JLOC,GI)*NU(JGL)
               VEDUDZ(IGL)=VEDUDZ(IGL)+DETWEI(GI)*N(ILOC,GI)*NZ(JLOC,GI)*NU(JGL)
                
               VEDVDX(IGL)=VEDVDX(IGL)+DETWEI(GI)*N(ILOC,GI)*NX(JLOC,GI)*NV(JGL)
               VEDVDY(IGL)=VEDVDY(IGL)+DETWEI(GI)*N(ILOC,GI)*NY(JLOC,GI)*NV(JGL)
               VEDVDZ(IGL)=VEDVDZ(IGL)+DETWEI(GI)*N(ILOC,GI)*NZ(JLOC,GI)*NV(JGL)
                
               VEDWDX(IGL)=VEDWDX(IGL)+DETWEI(GI)*N(ILOC,GI)*NX(JLOC,GI)*NW(JGL)
               VEDWDY(IGL)=VEDWDY(IGL)+DETWEI(GI)*N(ILOC,GI)*NY(JLOC,GI)*NW(JGL)
               VEDWDZ(IGL)=VEDWDZ(IGL)+DETWEI(GI)*N(ILOC,GI)*NZ(JLOC,GI)*NW(JGL)
                
            END DO
         END DO
      END DO
   end do

   do IGL=1,NONODS
      VEDUDX(IGL)=VEDUDX(IGL)/ML(IGL)
      VEDUDY(IGL)=VEDUDY(IGL)/ML(IGL)
      VEDUDZ(IGL)=VEDUDZ(IGL)/ML(IGL)
                
      VEDVDX(IGL)=VEDVDX(IGL)/ML(IGL)
      VEDVDY(IGL)=VEDVDY(IGL)/ML(IGL)
      VEDVDZ(IGL)=VEDVDZ(IGL)/ML(IGL)
                
      VEDWDX(IGL)=VEDWDX(IGL)/ML(IGL)
      VEDWDY(IGL)=VEDWDY(IGL)/ML(IGL)
      VEDWDZ(IGL)=VEDWDZ(IGL)/ML(IGL)
   END DO

! Collect halo's for the derivatives...VEDUDX etc.
   IF(PARA.EQ.1) THEN
      CALL HALGET(VEDUDX,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(VEDUDY,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(VEDUDZ,NONODS,NONODS,NNODP,halo_tag)
         
      CALL HALGET(VEDVDX,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(VEDVDY,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(VEDVDZ,NONODS,NONODS,NNODP,halo_tag)
         
      CALL HALGET(VEDWDX,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(VEDWDY,NONODS,NONODS,NNODP,halo_tag)
      CALL HALGET(VEDWDZ,NONODS,NONODS,NNODP,halo_tag)
   ENDIF
       
   IF(LESNOD.EQ.1) THEN
! For proper 4th order discretisation...
! Calculate E here
! E is the matrix of length-scales appropriate for each element

      do IGL=1,NONODS
                 
         CALL LESELA(NONODS,&
     & MUPTXX(IGL),MUPTXY(IGL),MUPTXZ(IGL),MUPTYY(IGL),MUPTYZ(IGL),MUPTZZ(IGL), &
     &             1,1,1,NLOC,EET)
!
! Mulitply E by 4*C_S^2*Sbars to give matrix of viscosities
         SXXGI = VEDUDX(IGL)
         SXYGI = 0.5*(VEDUDY(IGL)+VEDVDX(IGL)) 
         SXZGI = 0.5*(VEDUDZ(IGL)+VEDWDX(IGL)) 
         SYYGI = VEDVDY(IGL)  
         SYZGI = 0.5*(VEDVDZ(IGL)+VEDWDY(IGL)) 
         SZZGI = VEDWDZ(IGL) 
         SYXGI=SXYGI
         SZXGI=SXZGI
         SZYGI=SYZGI

         VIS=SQRT(2.* (SXXGI*SXXGI + SXYGI*SXYGI + SXZGI*SXZGI&
     &                +  SYXGI*SYXGI + SYYGI*SYYGI + SYZGI*SYZGI&
     &                +  SZXGI*SZXGI + SZYGI*SZYGI + SZZGI*SZZGI))
! THEN FIND TURBULENT 'VISCOSITIES'
         CS2=0.1**2
         FOURCS=4.*CS2

! Put a bit in here which multiplies E by FOURCS*VIS           
         do IDIM=1,6
            do JDIM=1,6
               EV(IDIM,JDIM,IGL)= FOURCS*VIS*EET(IDIM,JDIM,1)
            END DO
         END DO
           
      END DO
   ENDIF
         
 END SUBROUTINE LESDERIVS

 SUBROUTINE ANO_OCEAN_SCAL(&
     &  HIMXXDTW,HIMXYDTW,HIMXZDTW,&
     &  HIMYYDTW,HIMYZDTW,HIMZZDTW,&
     &  MXXDTW,MXYDTW,MXZDTW,&
     &  MYYDTW,MYZDTW,MZZDTW,&
     &  XD,YD,ZD,&
     &  ELE,NLOC,TOTELE,XONDGL,XNONOD,X,Y,Z,ISPHERE) 
! fIND HALF LENGTH SCALES HIMXXDTW and full 
! length scales MXXDTW based on anisotropic length scales...
   LOGICAL JUST_H_DIF

   PARAMETER(JUST_H_DIF=.false.)
! IF(JUST_H_DIF) Just horizontal diffusion (based on max dist between nodes)...         
   REAL HIMXXDTW,HIMXYDTW,HIMXZDTW,&
     &  HIMYYDTW,HIMYZDTW,HIMZZDTW,&
     &  MXXDTW,MXYDTW,MXZDTW,&
     &  MYYDTW,MYZDTW,MZZDTW,&
     &  XD,YD,ZD
   INTEGER ELE,NLOC,TOTELE,XONDGL(NLOC*TOTELE)
   INTEGER XNONOD
   REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
   INTEGER ISPHERE
! IF ISPHERE.GE.1 then on the sphere else z- is vertical.
! Local variables...
   REAL HIMYXDTW,HIMZXDTW,HIMZYDTW
   REAL MYXDTW,MZXDTW,MZYDTW
   REAL R1,R2,R3,DXMIN2,DXMAX2,DXMIN,DXMAX,DIST2
   INTEGER ILOC,JLOC,INOD,JNOD
   REAL RN,NORMX,NORMY,NORMZ
   REAL RVERTXX,RVERTXY,RVERTXZ
   REAL RVERTYX,RVERTYY,RVERTYZ
   REAL RVERTZX,RVERTZY,RVERTZZ
   REAL HORXX,HORXY,HORXZ 
   REAL HORYX,HORYY,HORYZ 
   REAL HORZX,HORZY,HORZZ 
       
   DXMIN2=+1.E+30
   DXMAX2=-1.E+30
   do ILOC=1,NLOC
      INOD=XONDGL((ELE-1)*NLOC+ILOC)
      do JLOC=ILOC+1,NLOC,1
         JNOD=XONDGL((ELE-1)*NLOC+JLOC)
         DIST2=(X(INOD)-X(JNOD))**2+(Y(INOD)-Y(JNOD))**2&
     &            +(Z(INOD)-Z(JNOD))**2
         DXMIN2=MIN(DXMIN2,DIST2)
         DXMAX2=MAX(DXMAX2,DIST2)
      END DO
   END DO
   DXMIN=SQRT(DXMIN2)
   DXMAX=SQRT(DXMAX2)
             
   IF(ISPHERE.GE.1) THEN
      RN=SQRT(XD**2+YD**2+ZD**2)
      NORMX=XD/RN
      NORMY=YD/RN
      NORMZ=ZD/RN
   ELSE
      NORMX=0.0
      NORMY=0.0
      NORMZ=1.0
   ENDIF
         
   RVERTXX=NORMX*NORMX
   RVERTXY=NORMX*NORMY
   RVERTXZ=NORMX*NORMZ
         
   RVERTYX=NORMY*NORMX
   RVERTYY=NORMY*NORMY
   RVERTYZ=NORMY*NORMZ
         
   RVERTZX=NORMZ*NORMX
   RVERTZY=NORMZ*NORMY
   RVERTZZ=NORMZ*NORMZ
         
   HORXX=-RVERTXX  + 1.0
   HORXY=-RVERTXY  
   HORXZ=-RVERTXZ  
         
   HORYX=-RVERTYX  
   HORYY=-RVERTYY  + 1.0
   HORYZ=-RVERTYZ 
         
   HORZX=-RVERTZX  
   HORZY=-RVERTZY  
   HORZZ=-RVERTZZ  + 1.0
         
   IF(JUST_H_DIF) THEN
! Just horizontal diffusion (based on max dist between nodes)...         
      HIMXXDTW=DXMAX*HORXX 
      HIMXYDTW=DXMAX*HORXY 
      HIMXZDTW=DXMAX*HORXZ 
         
      HIMYXDTW=DXMAX*HORYX 
      HIMYYDTW=DXMAX*HORYY 
      HIMYZDTW=DXMAX*HORYZ 
         
      HIMZXDTW=DXMAX*HORZX 
      HIMZYDTW=DXMAX*HORZY 
      HIMZZDTW=DXMAX*HORZZ 
   ELSE
! Both horizontal (max. dist between nodes) and vertical (min dist) diffusion.
      HIMXXDTW=DXMAX*HORXX + DXMIN*RVERTXX
      HIMXYDTW=DXMAX*HORXY + DXMIN*RVERTXY
      HIMXZDTW=DXMAX*HORXZ + DXMIN*RVERTXZ
         
      HIMYXDTW=DXMAX*HORYX + DXMIN*RVERTYX
      HIMYYDTW=DXMAX*HORYY + DXMIN*RVERTYY
      HIMYZDTW=DXMAX*HORYZ + DXMIN*RVERTYZ
         
      HIMZXDTW=DXMAX*HORZX + DXMIN*RVERTZX
      HIMZYDTW=DXMAX*HORZY + DXMIN*RVERTZY
      HIMZZDTW=DXMAX*HORZZ + DXMIN*RVERTZZ
   ENDIF
         
! MXXDTW contains L^2 and HIMXXDTW contains L.
   CALL SMLMA3( &
     &            MXXDTW,    MXYDTW,    MXZDTW,&
     &            R1,        MYYDTW,    MYZDTW,&
     &            R2,        R3,        MZZDTW,&
! =
     &            HIMXXDTW,    HIMXYDTW,    HIMXZDTW,&
     &            HIMXYDTW,    HIMYYDTW,    HIMYZDTW,&
     &            HIMXZDTW,    HIMYZDTW,    HIMZZDTW,&
! *          
     &            HIMXXDTW,    HIMXYDTW,    HIMXZDTW,&
     &            HIMXYDTW,    HIMYYDTW,    HIMYZDTW,&
     &            HIMXZDTW,    HIMYZDTW,    HIMZZDTW)

 END SUBROUTINE ANO_OCEAN_SCAL

 SUBROUTINE SIZEGIELETENS(NX,NY,NZ,NLOC,NGI,GI,&
     &     TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ )
! This sub caclulates the tensor for use in LES at quadrature point GI.
! The tensor contains L^2. 
   INTEGER NLOC,NGI
   INTEGER GI

   REAL TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ

   REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
      
!     Local variables...
   REAL RN
   REAL UDL,VDL,WDL

   INTEGER ILOC

   TENSXX=0.0
   TENSXY=0.0
   TENSXZ=0.0
   TENSYY=0.0
   TENSYZ=0.0
   TENSZZ=0.0

!     nb we want L^2 - at the moment we have L^2 on the diagonal.
   do  ILOC=1,NLOC

      RN=NX(ILOC,GI)**2+NY(ILOC,GI)**2+NZ(ILOC,GI)**2
      UDL=1.*NX(ILOC,GI)/RN
      VDL=1.*NY(ILOC,GI)/RN 
      WDL=1.*NZ(ILOC,GI)/RN
                  
      TENSXX=TENSXX + UDL*UDL
      TENSXY=TENSXY + UDL*VDL
      TENSXZ=TENSXZ + UDL*WDL
      TENSYY=TENSYY + VDL*VDL
      TENSYZ=TENSYZ + VDL*WDL
      TENSZZ=TENSZZ + WDL*WDL

   end do
 END SUBROUTINE SIZEGIELETENS
 
 SUBROUTINE ONEELETENS(ELE,X,Y,Z,DISOPT,&
     &     TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ, &
     &     TOTELE,DETWEI,NGI,NLOC, &
     &     XONDGL,XNONOD)
!     This sub calculates the NODE-WISE TENSOR TENS?? THAT 
!     REPRESENTS THE SIZE AND SHAPE OF THE SURROUNDING ELEMENTS.
!     DISOPT=LES option.
   INTEGER XNONOD
   INTEGER NLOC,NGI
   INTEGER TOTELE,DISOPT

   REAL TENSXX,TENSXY,TENSXZ
   REAL TENSYY,TENSYZ,TENSZZ
   REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
   INTEGER XONDGL(TOTELE*NLOC)

   REAL DETWEI(NGI)
!     HX,HY-characteristic length scales in x,y directions.
!     Local variables...
   REAL RN
   REAL XC,YC,ZC
   INTEGER NDIM
   PARAMETER(NDIM=3)
   REAL AA(NDIM,NDIM),V(NDIM,NDIM),D(NDIM),A(NDIM,NDIM)

   REAL UDL(NLOC*NLOC,NGI),VDL(NLOC*NLOC,NGI)
   REAL WDL(NLOC*NLOC,NGI)
   REAL GAMMA(NLOC*NLOC,NGI)

   INTEGER ELE,ILOC,L,L1,L2,IGLX1,IGLX2,IGLX,ID,NID
   INTEGER GI

   REAL HOVERQ
   REAL RWIND
   REAL RFACT,RT1,RT2,RT3,D1,D2,D3,VOLUME

   RWIND =1./REAL(6)
   NID=NLOC*NLOC

   TENSXX=0.0
   TENSXY=0.0
   TENSXZ=0.0
   TENSYY=0.0
   TENSYZ=0.0
   TENSZZ=0.0

!     This subroutine forms a contabution to the Right Hand Side
!     of Poissons pressure equation, as well as  F1 & F2.

   VOLUME=0.0
   do GI=1,NGI
      VOLUME=VOLUME+DETWEI(GI)

!     C The first is the old filter term, the second the new one MDP getting
!     c different results and stabiltiy for tidal applications ????
      RWIND =1./REAL(6)
      NID=NLOC*NLOC
!     **********calculate normalised velocitys across element...  
      ID=0
      do L1=1,NLOC
         IGLX1=XONDGL((ELE-1)*NLOC+L1) 
         do L2=1,NLOC
            IGLX2=XONDGL((ELE-1)*NLOC+L2) 
            ID=ID+1
            if(l1.eq.l2) then
               UDL(ID,GI)=0.0
               VDL(ID,GI)=0.0
               WDL(ID,GI)=0.0
               GAMMA(ID,GI)=0.0
            else
               UDL(ID,GI)=X(IGLX1)-X(IGLX2)
               VDL(ID,GI)=Y(IGLX1)-Y(IGLX2)
               WDL(ID,GI)=Z(IGLX1)-Z(IGLX2)

!     Normalise 
               RN=SQRT(UDL(ID,GI)**2+VDL(ID,GI)**2+WDL(ID,GI)**2)
               UDL(ID,GI)=UDL(ID,GI)/RN
               VDL(ID,GI)=VDL(ID,GI)/RN
               WDL(ID,GI)=WDL(ID,GI)/RN
!     HX,HY are the characteristic length scales in x,y directions. 
               HOVERQ=RN
               GAMMA(ID,GI)=RWIND*HOVERQ 
            endif
         END DO
      END DO
!     **********calculate normalised velocitys across element... 

   end do ! Was loop 331

   do  GI=1,NGI
      do  ID=1,NID

         RFACT=DETWEI(GI)*GAMMA(ID,GI)/VOLUME
                  
         TENSXX=TENSXX + RFACT*UDL(ID,GI)*UDL(ID,GI)
         TENSXY=TENSXY + RFACT*UDL(ID,GI)*VDL(ID,GI)
         TENSXZ=TENSXZ + RFACT*UDL(ID,GI)*WDL(ID,GI)
         TENSYY=TENSYY + RFACT*VDL(ID,GI)*VDL(ID,GI)
         TENSYZ=TENSYZ + RFACT*VDL(ID,GI)*WDL(ID,GI)
         TENSZZ=TENSZZ + RFACT*WDL(ID,GI)*WDL(ID,GI)
                  
!     USE THE COMPONENT OF DIFLIN THE X,Y & Z-DIRECTIONS 
!     RESPECTIVELY FOR C1T,C2T,C3T.
      end do
   end do
 
!     nb we want 1/L^2 - at the moment we have L on the diagonal.
!     Make sure the eigen-values are positive...
   AA(1,1)=TENSXX
   AA(1,2)=TENSXY
   AA(1,3)=TENSXZ
         
   AA(2,1)=TENSXY
   AA(2,2)=TENSYY
   AA(2,3)=TENSYZ
         
   AA(3,1)=TENSXZ
   AA(3,2)=TENSYZ
   AA(3,3)=TENSZZ
   CALL JACDIA(AA,V,D,NDIM,A)

   IF((DISOPT.EQ.45).OR.(DISOPT.EQ.46)) THEN
!     SET to metric which has 1/h^2 in it...
      D1=1./MAX(1.E-16,D(1)**2)
      D2=1./MAX(1.E-16,D(2)**2)
      D3=1./MAX(1.E-16,D(3)**2)
   ELSE IF((DISOPT.EQ.43).OR.(DISOPT.EQ.44)) THEN
!     set to inverse of metric which is a multiple of the tensor
      D1=MAX(1.E-16,D(1)**2)
      D2=MAX(1.E-16,D(2)**2)
      D3=MAX(1.E-16,D(3)**2)
   ELSE
!            ERROR("NOT A VALID OPTION FOR LES ASSEMBLED EQNS")
      STOP 9331
   ENDIF
   CALL SMLMA3( &
     &        TENSXX,       TENSXY,       TENSXZ,&
     &        RT1,          TENSYY,       TENSYZ,&
     &        RT2,          RT3,          TENSZZ, &
!     =
     &        V(1,1),   V(2,1),   V(3,1),  &
     &        V(1,2),   V(2,2),   V(3,2),  &
     &        V(1,3),   V(2,3),   V(3,3),             &
!     *          
     &        D1*V(1,1),   D1*V(1,2),   D1*V(1,3),  &
     &        D2*V(2,1),   D2*V(2,2),   D2*V(2,3),  &
     &        D3*V(3,1),   D3*V(3,2),   D3*V(3,3) ) 
 
 END SUBROUTINE ONEELETENS

 LOGICAL FUNCTION SAMLEVEL(NOD,COL, NONODS,NLEVEL)
   INTEGER NOD,COL, NONODS,NLEVEL
   integer ndslev
! See if nodes NOD and COL are on the same level of the structured mesh
   NDSLEV=NONODS/NLEVEL
   IF(INT((NOD-1)/NDSLEV).EQ.INT((COL-1)/NDSLEV)) THEN
      SAMLEVEL=.TRUE.
   ELSE
      SAMLEVEL=.FALSE.
   ENDIF
 END FUNCTION SAMLEVEL

 SUBROUTINE HIGH_VISC(VECX,VECY,VECZ, NU,NV,NW,&
     &     HIMATX,HIMATY,HIMATZ,HIMASS, &
     &     NONODS,FINDRM,COLM,CENTRM,NCOLM, &
     &     PARA,halo_tag,NNODP)
!     This sub forms VECX=VECX+KMAT*NU, VECY=VECY+KMAT*NV, VECZ=VECZ+KMAT*NW
!     in which KMAT is the high order diffusion matrix.
!     The discretised source is VECX,VECY,VECZ. 
!     This sub calculates the nodes for high order tets of degree ISPHERE 
   INTEGER NONODS
   REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
   REAL NU(NONODS),NV(NONODS),NW(NONODS)
   INTEGER NCOLM
   INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
   REAL HIMATX(NCOLM),HIMATY(NCOLM),HIMATZ(NCOLM),HIMASS(NCOLM)
   INTEGER PARA,halo_tag,NNODP
!     Local variables...
   INTEGER IMATST
   REAL TRELAX
   LOGICAL MOMSYM
   REAL TEROR1
   INTEGER TNOIT1,KITS
   INTEGER NOD,COUNT,COL
   INTEGER IIDUM(1)
   REAL RD(1)
!     *********ALOCAT MEMORY**************
   REAL, ALLOCATABLE, DIMENSION(:)::VECXU
   REAL, ALLOCATABLE, DIMENSION(:)::VECYU
   REAL, ALLOCATABLE, DIMENSION(:)::VECZU
      
   REAL, ALLOCATABLE, DIMENSION(:)::VECXV
   REAL, ALLOCATABLE, DIMENSION(:)::VECYV
   REAL, ALLOCATABLE, DIMENSION(:)::VECZV
      
   REAL, ALLOCATABLE, DIMENSION(:)::VECXW
   REAL, ALLOCATABLE, DIMENSION(:)::VECYW
   REAL, ALLOCATABLE, DIMENSION(:)::VECZW
      
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECXU
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECYU
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECZU
      
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECXV
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECYV
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECZV
      
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECXW
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECYW
   REAL, ALLOCATABLE, DIMENSION(:)::MIVECZW
      
      
   ALLOCATE(VECXU(NONODS))
   ALLOCATE(VECYU(NONODS))
   ALLOCATE(VECZU(NONODS))
      
   ALLOCATE(VECXV(NONODS))
   ALLOCATE(VECYV(NONODS))
   ALLOCATE(VECZV(NONODS))
      
   ALLOCATE(VECXW(NONODS))
   ALLOCATE(VECYW(NONODS))
   ALLOCATE(VECZW(NONODS))
      
   ALLOCATE(MIVECXU(NONODS))
   ALLOCATE(MIVECYU(NONODS))
   ALLOCATE(MIVECZU(NONODS))
      
   ALLOCATE(MIVECXV(NONODS))
   ALLOCATE(MIVECYV(NONODS))
   ALLOCATE(MIVECZV(NONODS))
      
   ALLOCATE(MIVECXW(NONODS))
   ALLOCATE(MIVECYW(NONODS))
   ALLOCATE(MIVECZW(NONODS))
      
      
   VECXU(1:NONODS)=0.0
   VECYU(1:NONODS)=0.0
   VECZU(1:NONODS)=0.0
      
   VECXV(1:NONODS)=0.0
   VECYV(1:NONODS)=0.0
   VECZV(1:NONODS)=0.0
      
   VECXW(1:NONODS)=0.0
   VECYW(1:NONODS)=0.0
   VECZW(1:NONODS)=0.0
   do NOD=1,NONODS
      do COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
         COL=COLM(COUNT)
         VECXU(NOD)=VECXU(NOD)+HIMATX(COUNT)*NU(COL)
         VECYU(NOD)=VECYU(NOD)+HIMATY(COUNT)*NU(COL)
         VECZU(NOD)=VECZU(NOD)+HIMATZ(COUNT)*NU(COL)
            
         VECXV(NOD)=VECXV(NOD)+HIMATX(COUNT)*NV(COL)
         VECYV(NOD)=VECYV(NOD)+HIMATY(COUNT)*NV(COL)
         VECZV(NOD)=VECZV(NOD)+HIMATZ(COUNT)*NV(COL)
            
         VECXW(NOD)=VECXW(NOD)+HIMATX(COUNT)*NW(COL)
         VECYW(NOD)=VECYW(NOD)+HIMATY(COUNT)*NW(COL)
         VECZW(NOD)=VECZW(NOD)+HIMATZ(COUNT)*NW(COL)
      END DO
   END DO
      
!     Solve M MIVECXU=VECXU  & M MIVECYU=VECYU  & M MIVECZU=VECZU 
!     Solve M MIVECXV=VECXV  & M MIVECYV=VECYV  & M MIVECZV=VECZV
!     Solve M MIVECXW=VECXW  & M MIVECYW=VECYW  & M MIVECZW=VECZW  
   IMATST=0
   IIDUM(1)=0

!     Solve M MIVECXU=VECXU  & M MIVECYU=VECYU  & M MIVECZU=VECZU 
   CALL SOLCG(MIVECXU, VECXU,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
   CALL SOLCG(MIVECYU, VECYU,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
   CALL SOLCG(MIVECZU, VECZU,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
      
!     Solve M MIVECXV=VECXV  & M MIVECYV=VECYV  & M MIVECZV=VECZV
   CALL SOLCG(MIVECXV, VECXV,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
   CALL SOLCG(MIVECYV, VECYV,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
   CALL SOLCG(MIVECZV, VECZV,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM, &
     &     halo_tag, KITS)
      
!     Solve M MIVECXW=VECXW  & M MIVECYW=VECYW  & M MIVECZW=VECZW 
   CALL SOLCG(MIVECXW, VECXW,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
   CALL SOLCG(MIVECYW, VECYW,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
   CALL SOLCG(MIVECZW, VECZW,NONODS,NONODS,NNODP,.TRUE., &
     &     HIMASS,FINDRM,COLM,NCOLM, NCOLM,&
     &     halo_tag, KITS)
      
!     Now the transpose part...    
   do NOD=1,NONODS
      do COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
         COL=COLM(COUNT)
         VECX(COL)=VECX(COL)-( HIMATX(COUNT)*MIVECXU(NOD) &
     &           +HIMATY(COUNT)*MIVECYU(NOD) &
     &           +HIMATZ(COUNT)*MIVECZU(NOD) )
         VECY(COL)=VECY(COL)-( HIMATX(COUNT)*MIVECXV(NOD) &
     &           +HIMATY(COUNT)*MIVECYV(NOD) &
     &           +HIMATZ(COUNT)*MIVECZV(NOD) )
         VECZ(COL)=VECZ(COL)-( HIMATX(COUNT)*MIVECXW(NOD) &
     &           +HIMATY(COUNT)*MIVECYW(NOD) &
     &           +HIMATZ(COUNT)*MIVECZW(NOD) )
      END DO
   END DO
      
 END SUBROUTINE HIGH_VISC

 SUBROUTINE SGSPROJ(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA, &
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB, &
     &     BIGM, ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC,BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT3,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     NNODP,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp,P, &
     &     SUF_TOP_OCEAN,&
     &     NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,BCU1val,BCV1val,BCW1val,&
     &     PARA,halo_tag,halo_tag_p,&
! form DGCMC
     &     NCMC,FINCMC,COLCMC,MIDCMC,&
     &     NDPSET,D3,NEWMES2, velocity_option_path, &
     &     NEXTTIM,FTHETA,GRAVTY,IUNST_FREE)
! if PREOPT3.GE.1000 then have a quadratic pressure else linear. 
! if PREOPT3.GE.10000 then use fully blown DG.  
! This code solves and SGS momentum and cty eqns...
   INTEGER MATSTR,IMP_P_ROT,IMP_P_ROT_UDG
   PARAMETER(MATSTR=0)
   INTEGER SNONOD,VNONOD,XNONOD,NCT,NCOLM,NLOC,NGI
   INTEGER MLOC,DISOPT2,DISOPT,DISOPN,DISCRA
   REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
   REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
   REAL DENPT(SNONOD)
   REAL ABSORB(SNONOD)
   REAL  DT,THETA,BETA
   LOGICAL NEXTTIM
   REAL FTHETA,GRAVTY
   INTEGER IUNST_FREE
! If IUNST_FREE=1 then solve for free surface height as part of pressure. 
   INTEGER TOTELE,NONODS,FREDOP,NBIGM,VERSIO,ISPHERE
   INTEGER BIG_NLOC
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
   REAL U(NONODS),V(NONODS),W(NONODS)
   REAL UNEW(NONODS),VNEW(NONODS),WNEW(NONODS)
   INTEGER ISUB,INUSUB,NSUBVLOC
   REAL NUSUB(TOTELE*NSUBVLOC),NVSUB(TOTELE*NSUBVLOC)
   REAL NWSUB(TOTELE*NSUBVLOC)
   REAL USUB(TOTELE*NSUBVLOC),VSUB(TOTELE*NSUBVLOC)
   REAL WSUB(TOTELE*NSUBVLOC)
   REAL UNEWSUB(TOTELE*NSUBVLOC),VNEWSUB(TOTELE*NSUBVLOC)
   REAL WNEWSUB(TOTELE*NSUBVLOC)
   REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
   REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
! This is what we add into the SGS model - a discretised source.      
   REAL VECX_SGSADD(TOTELE*NSUBVLOC),VECY_SGSADD(TOTELE*NSUBVLOC)
   REAL VECZ_SGSADD(TOTELE*NSUBVLOC)
       INTEGER NBUOY,NSOGRASUB
! soxgra,soygra,sozgra is the direction of gravity force when geobal.le.-10.or.(HYDROS.NE.0)
       REAL SOXGRA(SNONOD),SOYGRA(SNONOD),SOZGRA(SNONOD)
       REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
       REAL SOZGRASUB(TOTELE*NSOGRASUB)
   REAL BIGM(NBIGM)
   LOGICAL COGRAX,SUF_TOP_OCEAN
   INTEGER NOBCFSimp
   INTEGER BCT2FSimp(NOBCFSimp)
   INTEGER NOBCU,NOBCV,NOBCW
   INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
   REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
   REAL BCU1val(NOBCU),BCV1val(NOBCV),BCW1val(NOBCW)
   REAL ML(NONODS)
   REAL X(XNONOD),Y(XNONOD),Z(XNONOD),D0
   REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
   INTEGER ISUBSOU
   REAL SUBSOUX(ISUBSOU*TOTELE*NLOC),SUBSOUY(ISUBSOU*TOTELE*NLOC)
   REAL SUBSOUZ(ISUBSOU*TOTELE*NLOC)
   INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
   INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
   INTEGER XONDGL(TOTELE*NLOC)
      
   INTEGER FINDRM(NONODS+1),COLM(NCOLM)
   INTEGER CENTRM(NONODS)
   INTEGER FINDCT(FREDOP+1),COLCT(NCT)
   REAL P(FREDOP)

   REAL M(MLOC,NGI)
   REAL MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
   REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI)
   REAL, target:: NLY(NLOC,NGI),NLZ(NLOC,NGI)
   REAL WEIGHT(NGI)
   LOGICAL LUMP,UPWIND,OPWIND,XTPETV
   INTEGER PREOPT3
   LOGICAL ABSLUM,SLUMP
   INTEGER NNODP
   INTEGER GEOBAL
   REAL SCFACTH0
! GMRES solver...
   INTEGER PARA,halo_tag
   REAL MOMER1,MOMER2
   INTEGER MOMNO1,MOMNO2
! Pressure matrix soln:
   INTEGER PREPRE
   INTEGER PMISOU,halo_tag_p
   REAL PREERR
   INTEGER PRENOI,PINNOI,NNODPP
! form DGCMC
   INTEGER NCMC,FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)
   INTEGER NDPSET
   LOGICAL D3,NEWMES2
   character(len=*), intent(in):: velocity_option_path
! Local vartiables...
   INTEGER IDUM(1),KITS,NR2,ILENG,PREOPT
   INTEGER ELE,ILOC,NODDG,NOD,I,count
   integer preopt2
   REAL RDUM(1),RMIN
   LOGICAL ROTATDG

! Only the first time we enter this sub is SGSPROJ_START=1
   INTEGER SGSPROJ_START
   SAVE SGSPROJ_START
   DATA SGSPROJ_START /1/
   INTEGER, PARAMETER::halotagT10=3
       
   real, dimension(:), allocatable:: m2, m2lx, m2ly, m2lz
       
   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:)::PGNDGLNO
   type(csr_matrix), save:: pgM2
       
   REAL, ALLOCATABLE, SAVE, DIMENSION(:)::PG,PGOLD,POLD
   INTEGER PGNODS,PGNODSP
   SAVE PGNODS
   SAVE PGNODSP
       
   integer, dimension(:), pointer:: PGFINDRM, PGCENTRM, PGCOLM
   real, dimension(:), pointer:: difmath
   integer npgcolm
       
   INTEGER MLOC2
   ewrite(1,*) 'in subroutine sgsproj' 
       
   PREOPT2=PREOPT3
       
   IF(MOD(PREOPT3,10000).GE.1000) THEN
! Have a quadratic free surface...
      
      mloc2=10
      allocate(m2(mloc2*ngi), m2lx(mloc2*ngi), m2ly(mloc2*ngi), m2lz(mloc2*ngi))

! shape and mesh of the balanced pressure
! (since currently only mloc==10 works, we assume quadratic):
      call QuadTetShapes(M2, M2lx, M2ly, M2lz, Mloc2, ngi)
        
      IF(NEWMES2.AND.(SGSPROJ_START.EQ.0)) THEN
         DEALLOCATE(PGNDGLNO)
         DEALLOCATE(PG)
         DEALLOCATE(PGOLD)
      ENDIF
        
      IF(NEWMES2.OR.(SGSPROJ_START.EQ.1)) THEN
         ALLOCATE(PGNDGLNO(TOTELE*MLOC2))
! calculate PGNDGLNO...
         CALL MIDSNODS_DIFF(FINDRM,NONODS,PGNODS,        &
     &               CENTRM,COLM,NCOLM,                    &
     &               NDGLNO,PGNDGLNO,TOTELE,NLOC,MLOC2)
         PGNODSP=PGNODS
             
         CALL POSINM(pgM2%sparsity, TOTELE,PGNODS,MLOC2, &
     &                pgndglno, pgnods, mloc2, pgndglno, &
     &                name='SGSPressureSparsity')
         ALLOCATE(PG(PGNODS))
         ALLOCATE(PGOLD(PGNODS*IUNST_FREE))
! Interpolate P to get PG...
         CALL GESPSI_DIFF(P,PG,NONODS,PGNODS, &
     &                NLOC,MLOC2,TOTELE, &
     &                NDGLNO,PGNDGLNO)
         IF(IUNST_FREE.EQ.1) PGOLD=PG
         SGSPROJ_START=0
      ENDIF

! Prepare for next time step...      
      IF((IUNST_FREE.EQ.1).AND.NEXTTIM) PGOLD=PG

      pgfindrm => pgM2%sparsity%findrm
      pgcolm => pgM2%sparsity%colm
      pgcentrm => pgM2%sparsity%centrm
      npgcolm=size(pgcolm)
      difmath => pgM2%val
      
      CALL SGSPROJ2(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA,&
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
     &     BIGM, ML,&
     &     pgfindrm,pgcolm,npgcolm,PGNODS,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M2,M2LX,M2LY,M2LZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC2,BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT2,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PGNDGLNO,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     NNODP,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp,PG,PGOLD, &
     &     SUF_TOP_OCEAN,&
     &     NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,BCU1val,BCV1val,BCW1val,&
! GMRES solver...
     &     PARA,halo_tag,&
! Pressure matrix soln:
     &     halotagt10,&
     &     PGNODSP,&
! form DGCMC
     &     npgcolm,pgfindrm,pgcolm,pgcentrm,&
     &     NDPSET,D3,NEWMES2, velocity_option_path, &
     &     FTHETA,GRAVTY,IUNST_FREE)
     
! convert PG to P for plotting 
! and also for working out free surface height when IUNST_FREE=1... 
      CALL VISPRE_DIFF(P,PG,NONODS,PGNODS, &
     &           NLOC,MLOC2,TOTELE, &
     &           NDGLNO,PGNDGLNO)   
     
   ELSE
! Have a linear pressure...
      IF(NEWMES2.AND.(SGSPROJ_START.EQ.0)) THEN
         DEALLOCATE(POLD)
      ENDIF
        
      IF(NEWMES2.OR.(SGSPROJ_START.EQ.1)) THEN
         ALLOCATE(POLD(FREDOP*IUNST_FREE))
         IF(IUNST_FREE.EQ.1) POLD=P
         SGSPROJ_START=0
      ENDIF
! Prepare for next time step...      
      IF((IUNST_FREE.EQ.1).AND.NEXTTIM) POLD=P
      
      CALL SGSPROJ2(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA,&
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
     &     BIGM, ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC,BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT2,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     NNODP,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp,P,POLD, &
     &     SUF_TOP_OCEAN,&
     &     NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,BCU1val,BCV1val,BCW1val,&
! GMRES solver...
     &     PARA,halo_tag,&
! Pressure matrix soln:
     &     halo_tag_p, &
     &     NNODPP,&
! form DGCMC
     &     NCMC,FINCMC,COLCMC,MIDCMC,&
     &     NDPSET,D3,NEWMES2, velocity_option_path, &
     &     FTHETA,GRAVTY,IUNST_FREE)
   ENDIF
      
 END SUBROUTINE SGSPROJ

 SUBROUTINE SGSPROJ2(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA,&
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
     &     BIGM, ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC,BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT2,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     NNODP,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp,P,POLD, &
     &     SUF_TOP_OCEAN,&
     &     NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,BCU1val,BCV1val,BCW1val,&
! GMRES solver...
     &     PARA,halo_tag,&
! Pressure matrix soln:
     &     halo_tag_p,&
     &     nnodpp,&
! form DGCMC
     &     NCMC,FINCMC,COLCMC,MIDCMC,&
     &     NDPSET,D3,NEWMES2, velocity_option_path, &
     &     FTHETA,GRAVTY,IUNST_FREE)
! This code solves and SGS momentum and cty eqns...
   INTEGER MATSTR,IMP_P_ROT,IMP_P_ROT_UDG
   PARAMETER(MATSTR=0)
   INTEGER SNONOD,VNONOD,XNONOD,NCT,NCOLM,NLOC,NGI
   INTEGER MLOC,DISOPT2,DISOPT,DISOPN,DISCRA
   REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
   REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
   REAL DENPT(SNONOD)
   REAL ABSORB(SNONOD)
   REAL  DT,THETA,BETA
   REAL FTHETA,GRAVTY
   INTEGER IUNST_FREE
   INTEGER TOTELE,NONODS,FREDOP,NBIGM,VERSIO,ISPHERE
   INTEGER BIG_NLOC
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
   REAL U(NONODS),V(NONODS),W(NONODS)
   REAL UNEW(NONODS),VNEW(NONODS),WNEW(NONODS)
   INTEGER ISUB,INUSUB,NSUBVLOC
   REAL NUSUB(TOTELE*NSUBVLOC),NVSUB(TOTELE*NSUBVLOC)
   REAL NWSUB(TOTELE*NSUBVLOC)
   REAL USUB(TOTELE*NSUBVLOC),VSUB(TOTELE*NSUBVLOC)
   REAL WSUB(TOTELE*NSUBVLOC)
   REAL UNEWSUB(TOTELE*NSUBVLOC),VNEWSUB(TOTELE*NSUBVLOC)
   REAL WNEWSUB(TOTELE*NSUBVLOC)
   REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
   REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
! This is what we add into the SGS model - a discretised source.      
   REAL VECX_SGSADD(TOTELE*NSUBVLOC),VECY_SGSADD(TOTELE*NSUBVLOC)
   REAL VECZ_SGSADD(TOTELE*NSUBVLOC)
       INTEGER NBUOY,NSOGRASUB
! soxgra,soygra,sozgra is the direction of gravity force when geobal.le.-10.or.(HYDROS.NE.0)
       REAL SOXGRA(SNONOD),SOYGRA(SNONOD),SOZGRA(SNONOD)
       REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
       REAL SOZGRASUB(TOTELE*NSOGRASUB)
   REAL BIGM(NBIGM)
   LOGICAL COGRAX,SUF_TOP_OCEAN
   INTEGER NOBCFSimp
   INTEGER BCT2FSimp(NOBCFSimp)
   INTEGER NOBCU,NOBCV,NOBCW
   INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
   REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
   REAL BCU1val(NOBCU),BCV1val(NOBCV),BCW1val(NOBCW)
   REAL ML(NONODS)
   REAL X(XNONOD),Y(XNONOD),Z(XNONOD),D0
   REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
   INTEGER ISUBSOU
   REAL SUBSOUX(ISUBSOU*TOTELE*NLOC),SUBSOUY(ISUBSOU*TOTELE*NLOC)
   REAL SUBSOUZ(ISUBSOU*TOTELE*NLOC)
   INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
   INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
   INTEGER XONDGL(TOTELE*NLOC)
      
   INTEGER FINDRM(NONODS+1),COLM(NCOLM)
   INTEGER CENTRM(NONODS)
   INTEGER FINDCT(FREDOP+1),COLCT(NCT)
   REAL P(FREDOP),POLD(FREDOP*IUNST_FREE)

   REAL M(MLOC,NGI)
   REAL MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
   REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI)
   REAL, target:: NLY(NLOC,NGI),NLZ(NLOC,NGI)
   REAL WEIGHT(NGI)
   LOGICAL LUMP,UPWIND,OPWIND,XTPETV
   INTEGER PREOPT2
   LOGICAL ABSLUM,SLUMP
   INTEGER NNODP
   INTEGER GEOBAL
   REAL SCFACTH0
! GMRES solver...
   INTEGER PARA,halo_tag
! Pressure matrix soln:
   INTEGER halo_tag_p
   INTEGER NNODPP
! form DGCMC
   INTEGER NCMC,FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)
   INTEGER NDPSET
   LOGICAL D3,NEWMES2
   character(len=*), intent(in):: velocity_option_path
! Local vartiables...
   INTEGER IDUM(1),KITS,NR2,ILENG,PREOPT
   INTEGER ELE,ILOC,NODDG,NOD,I,j,count
   REAL RDUM(1),RMIN
   LOGICAL ROTATDG
      
   REAL, ALLOCATABLE, DIMENSION(:)::VECX
   REAL, ALLOCATABLE, DIMENSION(:)::VECY
   REAL, ALLOCATABLE, DIMENSION(:)::VECZ
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::CMAT_STORE
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::DMAT_STORE
   REAL, ALLOCATABLE, DIMENSION(:,:)::DINVSOUSUBU_STORE
   REAL, ALLOCATABLE, DIMENSION(:,:)::DINVSOUSUBV_STORE
   REAL, ALLOCATABLE, DIMENSION(:,:)::DINVSOUSUBW_STORE
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT11
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT12
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT13
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT21
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT22
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT23
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT31
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT32
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT33
       
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA11
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA12
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA13
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA21
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA22
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA23
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA31
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA32
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA33
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::MASSDGELE
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::C1TDGELE
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::C2TDGELE
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::C3TDGELE
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::G1TDGELE
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::G2TDGELE
   REAL, ALLOCATABLE, DIMENSION(:,:,:)::G3TDGELE
   REAL, ALLOCATABLE, DIMENSION(:)::UDG
   REAL, ALLOCATABLE, DIMENSION(:)::VDG
   REAL, ALLOCATABLE, DIMENSION(:)::WDG
   REAL, ALLOCATABLE, DIMENSION(:)::UDGOLD
   REAL, ALLOCATABLE, DIMENSION(:)::VDGOLD
   REAL, ALLOCATABLE, DIMENSION(:)::WDGOLD
   REAL, ALLOCATABLE, DIMENSION(:)::DGCMC
   REAL, ALLOCATABLE, DIMENSION(:)::DP
   REAL, ALLOCATABLE, DIMENSION(:)::RHSVEC
   REAL, ALLOCATABLE, DIMENSION(:)::UVW
   REAL, ALLOCATABLE, DIMENSION(:)::VECXYZ
   REAL, ALLOCATABLE, DIMENSION(:)::DGU_THETA,DGV_THETA,DGW_THETA
   LOGICAL, ALLOCATABLE, DIMENSION(:)::ROT_THIS_DGNOD
   REAL, ALLOCATABLE, DIMENSION(:)::DGNX,DGNY,DGNZ, &
     &                                  DGT1X,DGT1Y,DGT1Z, &
     &                                  DGT2X,DGT2Y,DGT2Z, &
     &                                  DGDM1,DGDM2,DGDM3
   REAL, ALLOCATABLE, DIMENSION(:)::PRES_RHS
   character(len=OPTION_PATH_LEN) :: pressure_option_path
       
      ! warning needs to be fixed for multi-material/phase:
   ewrite(1,*) 'in subroutine sgsproj2' 
   pressure_option_path='/material_phase[0]/scalar_field::Pressure'
       
   IF(PREOPT2.GE.10000) THEN
! real DG...
      PREOPT=PREOPT2
      CALL DGPROJ(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, &
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA,&
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
     &     BIGM, ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC,BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT2,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     NNODP,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp,P,POLD, &
     &     SUF_TOP_OCEAN,&
     &     NOBCU,NOBCV,NOBCW,BCU1val,BCV1val,BCW1val,BCU2,BCV2,BCW2,&
! GMRES solver...
     &     PARA,halo_tag,&
! Pressure matrix soln:
     &     halo_tag_p,&
     &     NNODPP,&
! form DGCMC
     &     NCMC,FINCMC,COLCMC,MIDCMC,&
     &     NDPSET,D3,NEWMES2,velocity_option_path, &
     &     FTHETA,GRAVTY,IUNST_FREE)
      RETURN
   ENDIF
       
   ILENG=3*NONODS
   ALLOCATE(VECX(NONODS))
   ALLOCATE(VECY(NONODS))
   ALLOCATE(VECZ(NONODS))
   ALLOCATE(ROT_THIS_DGNOD(TOTELE*NLOC))
   ALLOCATE(DGNX(TOTELE*NLOC))
   ALLOCATE(DGNY(TOTELE*NLOC))
   ALLOCATE(DGNZ(TOTELE*NLOC))
   ALLOCATE(DGT1X(TOTELE*NLOC))
   ALLOCATE(DGT1Y(TOTELE*NLOC))
   ALLOCATE(DGT1Z(TOTELE*NLOC))
   ALLOCATE(DGT2X(TOTELE*NLOC))
   ALLOCATE(DGT2Y(TOTELE*NLOC))
   ALLOCATE(DGT2Z(TOTELE*NLOC))
   ALLOCATE(DGDM1(TOTELE*NLOC))
   ALLOCATE(DGDM2(TOTELE*NLOC))
   ALLOCATE(DGDM3(TOTELE*NLOC))
   
   ROTATDG=.FALSE.
! IF PREOPT2.GE.10 then put Coriolis into pressure calc. 
! IF PREOPT2.GE.100 then also put Coriolis in distribution of velocity 
! between SGS and global solns.  
! IF PREOPT=2 then treat every boundary like a free surface - need 
! to do this for open boundaries. 
   IMP_P_ROT=0
   IMP_P_ROT_UDG=0
   IF(MOD(PREOPT2,100).GE.10)   IMP_P_ROT=1
   IF(MOD(PREOPT2,1000).GE.100) IMP_P_ROT_UDG=1
   PREOPT=MOD(MOD(PREOPT2,10),100)
     
   ALLOCATE(CMAT_STORE(TOTELE,NLOC,NLOC))
   ALLOCATE(DMAT_STORE(TOTELE,NLOC,NLOC))
   ALLOCATE(DINVSOUSUBU_STORE(TOTELE,NLOC))
   ALLOCATE(DINVSOUSUBV_STORE(TOTELE,NLOC))
   ALLOCATE(DINVSOUSUBW_STORE(TOTELE,NLOC))
       
       ALLOCATE(COMAT11(TOTELE,NLOC))
       ALLOCATE(COMAT12(TOTELE,NLOC))
       ALLOCATE(COMAT13(TOTELE,NLOC))
       ALLOCATE(COMAT21(TOTELE,NLOC))
       ALLOCATE(COMAT22(TOTELE,NLOC))
       ALLOCATE(COMAT23(TOTELE,NLOC))
       ALLOCATE(COMAT31(TOTELE,NLOC))
       ALLOCATE(COMAT32(TOTELE,NLOC))
       ALLOCATE(COMAT33(TOTELE,NLOC))
       
       ALLOCATE(BOYMA11(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA12(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA13(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA21(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA22(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA23(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA31(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA32(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA33(NBUOY*TOTELE,NLOC))
       
   ALLOCATE(MASSDGELE(TOTELE,NLOC,NLOC))
   ALLOCATE(C1TDGELE(TOTELE,MLOC,NLOC))
   ALLOCATE(C2TDGELE(TOTELE,MLOC,NLOC))
   ALLOCATE(C3TDGELE(TOTELE,MLOC,NLOC))
   ALLOCATE(G1TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC))
   ALLOCATE(G2TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC))
   ALLOCATE(G3TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC))
   ALLOCATE(PRES_RHS(FREDOP))
       
   ALLOCATE(DGCMC(NCMC*IUNST_FREE))
   
! Get the DG soln UDGOLD from SGS soln...
   ALLOCATE(UDGOLD(TOTELE*NLOC))
   ALLOCATE(VDGOLD(TOTELE*NLOC))
   ALLOCATE(WDGOLD(TOTELE*NLOC))
   CALL GETDG_FROM_SUB(UDGOLD,VDGOLD,WDGOLD,U,V,W,USUB,VSUB,WSUB,&
     &                     NONODS,NLOC,TOTELE,NDGLNO)

! form momentum eqns and put pressure into rhs,
! also put b.c's into matrix...
   CALL DIFF3DSGS(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD,&
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA,&
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB, &
     &     COMAT11,COMAT12,COMAT13,COMAT21,&
     &     COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &     BOYMA11,BOYMA12,BOYMA13,BOYMA21,&
     &     BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33, &
     &     VECX,VECY,VECZ, &
     &     BIGM, ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT,&
     &     ABSLUM,SLUMP, &
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     nnodp,para,halo_tag,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBVLOC,NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp, &
     &     CMAT_STORE,DMAT_STORE,&
     &     DINVSOUSUBU_STORE,DINVSOUSUBV_STORE,&
     &     DINVSOUSUBW_STORE,&
     &     MASSDGELE,C1TDGELE,C2TDGELE,C3TDGELE,P,&
     &     FTHETA,GRAVTY,IUNST_FREE,G1TDGELE,G2TDGELE,G3TDGELE,POLD,&
     &     DGCMC,NCMC,FINCMC,COLCMC,MIDCMC, &
     &     SUF_TOP_OCEAN,&
     &     ROTATDG,ROT_THIS_DGNOD,&
     &     DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &     DGDM1,DGDM2,DGDM3,&
     &     PRES_RHS, velocity_option_path )

! Apply standard b.c.'s
   CALL APPLY_STAN_BC(BIGM,CENTRM,NONODS,NCOLM,NBIGM,&
     &                     U,V,W,DT,&
     &                     VECX,VECY,VECZ,&
     &                     NOBCU,NOBCV,NOBCW,&
     &                     BCU1,BCV1,BCW1,&
     &                     BCU2,BCV2,BCW2, &
     &                     TOTELE,NLOC,NDGLNO,&
     &                     DGDM1,DGDM2,DGDM3)
     
   CALL PMINMX(U,       NONODS,'******U  ')
   CALL PMINMX(V,       NONODS,'******V  ')
   CALL PMINMX(W,       NONODS,'******W  ')
   CALL PMINMX(Usub,       totele*nloc,'******Usub  ')
   CALL PMINMX(Vsub,       totele*nloc,'******Vsub  ')
   CALL PMINMX(Wsub,       totele*nloc,'******Wsub  ')
      
! Solve momentum eqns... 
   ILENG=3*NONODS
   NR2=3*ILENG
   ALLOCATE(UVW(ILENG))
   UVW=0.0
   ALLOCATE(VECXYZ(ILENG))
   UVW(1:NONODS)=UNEW(1:NONODS)
   UVW(NONODS+1:2*NONODS)=VNEW(1:NONODS)
   UVW(2*NONODS+1:3*NONODS)=WNEW(1:NONODS)
   VECXYZ(1:NONODS)=VECX(1:NONODS)
   VECXYZ(NONODS+1:2*NONODS)=VECY(1:NONODS)
   VECXYZ(2*NONODS+1:3*NONODS)=VECZ(1:NONODS)
   CALL GMRES(UVW,VECXYZ,NONODS,ILENG,NNODP,.FALSE.,&
     &       BIGM,FINDRM,COLM,NCOLM,NBIGM,&
     & halo_tag,KITS,&
! Where to find new solver options     
     &       option_path=velocity_option_path)

! UVW contains the acceleration.
   U(1:NONODS)=UVW(1:NONODS)           
   V(1:NONODS)=UVW(NONODS+1:2*NONODS)
   W(1:NONODS)=UVW(2*NONODS+1:3*NONODS)
   DEALLOCATE(UVW)
   DEALLOCATE(VECXYZ)
   CALL PMINMX(U,       NONODS,'******U  ')
   CALL PMINMX(V,       NONODS,'******V  ')
   CALL PMINMX(W,       NONODS,'******W  ')
       
   CALL SOLV_FOR_SUB_UVW(USUB,VSUB,WSUB, U,V,W,&
     &               NLOC,BIG_NLOC,TOTELE,NONODS,NDGLNO,&
     &               CMAT_STORE,DMAT_STORE,&
     &               DINVSOUSUBU_STORE,DINVSOUSUBV_STORE,&
     &               DINVSOUSUBW_STORE,&
     &               COMAT11,COMAT12,COMAT13,COMAT21,&
     &               COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)

   CALL PMINMX(Usub,       totele*nloc,'******Usub  ')
   CALL PMINMX(Vsub,       totele*nloc,'******Vsub  ')
   CALL PMINMX(Wsub,       totele*nloc,'******Wsub  ')

! *********free up memory**********
   DEALLOCATE(CMAT_STORE)
   DEALLOCATE(DMAT_STORE)
   DEALLOCATE(DINVSOUSUBU_STORE)
   DEALLOCATE(DINVSOUSUBV_STORE)
   DEALLOCATE(DINVSOUSUBW_STORE)
! *********free up memory**********
     
! Rotate the DGCT matrices...
   CALL ROTDGCT(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,&
     &         NONODS,FREDOP,&
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3) 
     
      IF(IUNST_FREE.EQ.1) THEN
! Rotate the DGgT matrices FOR FREE SURFACE...
       CALL ROTDGCT(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,&
     &         G1TDGELE,G2TDGELE,G3TDGELE,&
     &         NONODS,FREDOP,&
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3) 
      ENDIF   


! Get the DG soln UDG from SGS soln...
   ALLOCATE(UDG(TOTELE*NLOC))
   ALLOCATE(VDG(TOTELE*NLOC))
   ALLOCATE(WDG(TOTELE*NLOC))
   CALL GETDG_FROM_SUB(UDG,VDG,WDG,U,V,W,USUB,VSUB,WSUB,&
     &                     NONODS,NLOC,TOTELE,NDGLNO)

! form r.h.s of pressure eqn...
       ALLOCATE(RHSVEC(FREDOP))
       IF(IUNST_FREE.EQ.0) THEN
         CALL DGCTMULT(RHSVEC,.TRUE.,PNDGLN,MLOC,NLOC,TOTELE,FREDOP,&
     &                 DT,C1TDGELE,C2TDGELE,C3TDGELE,&
     &                 UDG,VDG,WDG)
       ELSE
         ewrite(1,*) 'before rhsvec'
         ALLOCATE(DGU_THETA(TOTELE*NLOC))
         ALLOCATE(DGV_THETA(TOTELE*NLOC))
         ALLOCATE(DGW_THETA(TOTELE*NLOC))
         DGU_THETA=UDG+((1.-FTHETA)/FTHETA)*UDGOLD
         DGV_THETA=VDG+((1.-FTHETA)/FTHETA)*VDGOLD
         DGW_THETA=WDG+((1.-FTHETA)/FTHETA)*WDGOLD 
         CALL DGCTMULT(RHSVEC,.TRUE.,PNDGLN,MLOC,NLOC,TOTELE,FREDOP,&
     &                 DT,C1TDGELE,C2TDGELE,C3TDGELE,&
     &                 DGU_THETA,DGV_THETA,DGW_THETA)
         DEALLOCATE(DGU_THETA)
         DEALLOCATE(DGV_THETA)
         DEALLOCATE(DGW_THETA)
! RHSVEC=RHSVEC-(M_S/(GRAVTY*DT*DT))*(P-POLD)
         DO I=1,FREDOP
           DO COUNT=FINCMC(I),FINCMC(I+1)-1
             J=COLCMC(COUNT)
             RHSVEC(I)=RHSVEC(I)-DGCMC(COUNT)*(P(J)-POLD(J))
           END DO
         END DO
       ENDIF
       RHSVEC=RHSVEC-PRES_RHS/DT
   
   
! form DGCMC...
      IF(IUNST_FREE.EQ.0) THEN
        DEALLOCATE(DGCMC)
        ALLOCATE(DGCMC(NCMC))
        DGCMC=0.0
      ENDIF

      IF((NBUOY.NE.0).AND.(IMP_P_ROT.EQ.0)) THEN
! Treating just buoyancy implcitly in pressure...
       CALL GETDGCMC(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,DGCMC,.FALSE.,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,FREDOP,NCMC,&
     &         FINCMC,COLCMC,MIDCMC,&
     &         DT,NBUOY,&
     &         BOYMA11,BOYMA12,BOYMA13,BOYMA21,&
     &         BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33, & 
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3 )
      ELSE
! Treating Coriolis and buoyancy implcitly in pressure (no symm matrix)...
       CALL GETDGCMC(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,DGCMC,.FALSE.,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,FREDOP,NCMC,&
     &         FINCMC,COLCMC,MIDCMC,&
     &         DT,IMP_P_ROT,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3 )
      ENDIF


! Solve pressure eqn...
! region has become isolated. 
   ALLOCATE(DP(FREDOP))

   RMIN=1.0E+20

   do count=1,ncmc
      if((colcmc(count).lt.1).or.(colcmc(count).gt.fredop)) then
         ewrite(-1,*)'count,colcmc(count):',count,colcmc(count)
         FLAbort('Oh dear')
      endif
   end do

   do I=1,FREDOP
      RMIN=MIN(RMIN,ABS(DGCMC(MIDCMC(I))))
   END DO

   ewrite(2,*)'RMIN=',RMIN
   DP=0.0

   IF(IMP_P_ROT.EQ.0) THEN
      CALL SOLCG(DP,RHSVEC,FREDOP,FREDOP,NNODPP,.FALSE., &
     &            DGCMC,FINCMC,COLCMC,NCMC,NCMC,&
     &            halo_tag_p,KITS, &
     &            option_path=pressure_option_path)
   ELSE
      CALL GMRES(DP,RHSVEC,FREDOP,FREDOP,NNODPP,.FALSE.,&
     &       DGCMC,FINCMC,COLCMC,NCMC,NCMC,&
     &       halo_tag_p,KITS, &
     &       option_path=pressure_option_path)
   ENDIF

   P=P+DP

! Solve DG mass matrix eqn for velocity update M_DG dU_DG=dt*C*dP
! then choose how to distribute dU_DG between global and SGS...


      IF((NBUOY.NE.0).AND.(IMP_P_ROT.EQ.0)) THEN
       CALL GET_SGS_PROJ_VELS(U,V,W,USUB,VSUB,WSUB, DP,&
     &         BIG_NLOC,NLOC,TOTELE,FREDOP,PNDGLN,MLOC,&
     &         DT,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,&
     &         NBUOY,&
     &         IMP_P_ROT_UDG,NBIGM,BIGM,&
     &         BOYMA11,BOYMA12,BOYMA13,BOYMA21,&
     &         BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33, &
! Th following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,&
! Extra variables used in GET_USUB_FROM_UDG...
     &         DGDM1,DGDM2,DGDM3,&
     &         NDGLNO,&
     &         FINDRM,COLM,CENTRM,NCOLM,&
     &         NNODP,PARA,halo_tag )
      ELSE
       CALL GET_SGS_PROJ_VELS(U,V,W,USUB,VSUB,WSUB, DP,&
     &         BIG_NLOC,NLOC,TOTELE,FREDOP,PNDGLN,MLOC,&
     &         DT,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,&
     &         IMP_P_ROT,&
     &         IMP_P_ROT_UDG,NBIGM,BIGM,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
! Th following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,&
! Extra variables used in GET_USUB_FROM_UDG...
     &         DGDM1,DGDM2,DGDM3,&
     &         NDGLNO,&
     &         FINDRM,COLM,CENTRM,NCOLM,&
     &         NNODP,PARA,halo_tag )
      ENDIF

 END SUBROUTINE SGSPROJ2

 SUBROUTINE ROTDGCT(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,&
     &         NONODS,FREDOP,&
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3)
   INTEGER BIG_NLOC,NLOC,TOTELE,NDGLNO(TOTELE*NLOC),NDPSET
   INTEGER NONODS,FREDOP
   INTEGER MLOC,PNDGLN(TOTELE*MLOC)
   LOGICAL D3
   REAL C1TDGELE(TOTELE,MLOC,NLOC),C2TDGELE(TOTELE,MLOC,NLOC)
   REAL C3TDGELE(TOTELE,MLOC,NLOC)
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY     
! 
   REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
   REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
   REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
   LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
   REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Local variables...
   INTEGER ELE,ILOC,JLOC,ILOC3,JLOC3,NODDG
   REAL C1T2(MLOC,NLOC),C1T2TRAN(NLOC,MLOC)
   REAL C2T2(MLOC,NLOC),C2T2TRAN(NLOC,MLOC)
   REAL C3T2(MLOC,NLOC),C3T2TRAN(NLOC,MLOC)
   REAL CT(MLOC,BIG_NLOC),CTRT(MLOC,BIG_NLOC)
   REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
   LOGICAL ROT_THIS_ELE
       
   IF(.NOT.ROTATDG) RETURN

   do ELE=1,TOTELE
       
      ROT_THIS_ELE=.FALSE.
      do ILOC=1,NLOC
         NODDG=(ELE-1)*NLOC+ILOC
         IF(ROT_THIS_DGNOD(NODDG)) ROT_THIS_ELE=.TRUE.
      END DO
      IF(.NOT.ROT_THIS_ELE) THEN
         C1T2(1:MLOC,1:NLOC)=C1TDGELE(ELE,1:MLOC,1:NLOC)
         C2T2(1:MLOC,1:NLOC)=C2TDGELE(ELE,1:MLOC,1:NLOC)
         C3T2(1:MLOC,1:NLOC)=C3TDGELE(ELE,1:MLOC,1:NLOC)

! Now for the rotation approach...
         RLOC=0.0
         do ILOC=1,NLOC
            NODDG=(ELE-1)*NLOC+ILOC
            IF(ROT_THIS_DGNOD(NODDG)) THEN
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
               RLOC(ILOC,ILOC)       =DGT1X(NODDG)
               RLOC(ILOC,ILOC+NLOC)  =DGT1Y(NODDG)
               RLOC(ILOC,ILOC+2*NLOC)=DGT1Z(NODDG)
               RLOC(ILOC+NLOC,ILOC)       =DGT2X(NODDG)
               RLOC(ILOC+NLOC,ILOC+NLOC)  =DGT2Y(NODDG)
               RLOC(ILOC+NLOC,ILOC+2*NLOC)=DGT2Z(NODDG)
               RLOC(ILOC+2*NLOC,ILOC)       =DGNX(NODDG)
               RLOC(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(NODDG)
               RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(NODDG)
            ELSE
               RLOC(ILOC,ILOC)              =1.0
               RLOC(ILOC+NLOC,ILOC+NLOC)    =1.0
               RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=1.0
            ENDIF
         END DO
         do ILOC3=1,BIG_NLOC
            do JLOC3=1,BIG_NLOC
               RTLOC(ILOC3,JLOC3)=RLOC(JLOC3,ILOC3)
            END DO
         END DO
! Calculate C^T R^T...
         CT(1:MLOC,1:NLOC)         =C1T2(1:MLOC,1:NLOC)
         CT(1:MLOC,1+NLOC:2*NLOC)  =C2T2(1:MLOC,1:NLOC)
         CT(1:MLOC,1+2*NLOC:3*NLOC)=C3T2(1:MLOC,1:NLOC)
         CALL ABMATRIXMUL(CTRT,   CT,MLOC,BIG_NLOC, &
     &                              RTLOC,BIG_NLOC,BIG_NLOC)
! Re-define C1TDGELE,C2TDGELE,C3TDGELE
         C1TDGELE(ELE,1:MLOC,1:NLOC)=CTRT(1:MLOC,1:NLOC)
         C2TDGELE(ELE,1:MLOC,1:NLOC)=CTRT(1:MLOC,NLOC+1:2*NLOC)
         C2TDGELE(ELE,1:MLOC,1:NLOC)=CTRT(1:MLOC,2*NLOC+1:3*NLOC)
! ENDOF IF(.NOT.ROT_THIS_ELE) THEN...
      ENDIF

   END DO
       
 END SUBROUTINE ROTDGCT

 SUBROUTINE GETDGCMC(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,DGCMC,ZERO_DGCMC, &
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,FREDOP,NCMC,&
     &         FINCMC,COLCMC,MIDCMC,&
     &         DT,IMP_P_ROT,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33,  &
! Th following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3 )
! form DGCMC
! ALFST2= max volume fraction of solids. 
! SPRESS=2nd pressure scaling parameter.
   INTEGER BIG_NLOC,NLOC,TOTELE
   INTEGER NDGLNO(TOTELE*NLOC)
   INTEGER NONODS,FREDOP
   INTEGER NCMC
   INTEGER MLOC,PNDGLN(TOTELE*MLOC)
   INTEGER FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)
! If ROTAT then it is assumed that C1T etc already conbtain the 
! rotation matrices via C1T=R^T*C1T
! This sub obtains the pressure matrices CMC & CMCP(which is non-sym). 
   REAL C1TDGELE(TOTELE,MLOC,NLOC),C2TDGELE(TOTELE,MLOC,NLOC)
   REAL C3TDGELE(TOTELE,MLOC,NLOC)
   REAL MASSDGELE(TOTELE,NLOC,NLOC)
   REAL DGCMC(NCMC)
   LOGICAL ZERO_DGCMC
   REAL DT
   INTEGER IMP_P_ROT
         REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC),&
     &        COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC),&
     &        COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
! NDVOLF=volume frac of all phases. 
! Put volume fraction of this phase into RMEM(1) & make a no not 
! too close to 0.0 if not multi-phase then set volume frac=1.
! NB NMPDIA=NPHASE*NPHASE*NONODS*NDIM*NDIM.
! NMXBD=maximum value of NDIM*NPHASE
! IF GETINV then get the inverse of MPDIAG with b.c's in. 
   REAL INFINY
   PARAMETER(INFINY=1.e+20)

   LOGICAL D3
   INTEGER NDPSET
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY     

   REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
   REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
   REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
   LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
   REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! A COLOURIN TECHNIQUE WORK BE MORE EFFICIENT 
! HERE(SEE BOTTON OF SUB). 
! The following forms the matrix C1T RML C1 + C2T RML C2.
! Often there is only one diagonal so that ML=ML1,ML2,ML3. 
! If D3 then suppose the flow is 3-D else suppose it is 2-D. 
! May have to alter this to make it more efficient. 
! If a multi-phase solution then accumulate the 
! pressure matrix
! Local mem...
   REAL MASS2(NLOC,NLOC),MASS2INV(NLOC,NLOC)
   REAL C1T2(MLOC,NLOC),C2T2(MLOC,NLOC),C3T2(MLOC,NLOC)
   REAL C1T2TRAN(NLOC,MLOC),C2T2TRAN(NLOC,MLOC)
   REAL C3T2TRAN(NLOC,MLOC)
   REAL CMCLOC(MLOC,MLOC)
   REAL CMCLOC1(MLOC,MLOC),CMCLOC2(MLOC,MLOC),CMCLOC3(MLOC,MLOC)
   REAL MINVC1(NLOC,MLOC),MINVC2(NLOC,MLOC),MINVC3(NLOC,MLOC)
   REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
   REAL MASSBLOC(BIG_NLOC,BIG_NLOC)
   REAL RTMR(BIG_NLOC,BIG_NLOC),MRLOC(BIG_NLOC,BIG_NLOC)
   REAL RTMINVR(BIG_NLOC,BIG_NLOC)
   REAL CT(MLOC,BIG_NLOC),CTTRAN(BIG_NLOC,MLOC)
   REAL RTMINVRC(BIG_NLOC,MLOC)
   INTEGER ELE,ILOC,JLOC,ILOC3,JLOC3,GLOBI,GLOBJ,NODDG
   INTEGER COUNT
   LOGICAL ROT_THIS_ELE
   real rmin,rmax,rmin2,rmax2,rmin3,rmax3

   ewrite(1,*) 'inside GETDGCMC'

   IF(ZERO_DGCMC) DGCMC(1:NCMC)=0.0

   do ELE=1,TOTELE
       
      MASS2(1:NLOC,1:NLOC)=MASSDGELE(ELE,1:NLOC,1:NLOC)
      C1T2(1:MLOC,1:NLOC)=C1TDGELE(ELE,1:MLOC,1:NLOC)
      C2T2(1:MLOC,1:NLOC)=C2TDGELE(ELE,1:MLOC,1:NLOC)
      C3T2(1:MLOC,1:NLOC)=C3TDGELE(ELE,1:MLOC,1:NLOC)
      do ILOC=1,MLOC
         do JLOC=1,NLOC
            C1T2TRAN(JLOC,ILOC)=C1T2(ILOC,JLOC)
            C2T2TRAN(JLOC,ILOC)=C2T2(ILOC,JLOC)
            C3T2TRAN(JLOC,ILOC)=C3T2(ILOC,JLOC)
         END DO
      END DO
      ROT_THIS_ELE=.FALSE.
      do ILOC=1,NLOC
         NODDG=(ELE-1)*NLOC+ILOC
         IF(ROT_THIS_DGNOD(NODDG)) ROT_THIS_ELE=.TRUE.
      END DO
      IF((.NOT.ROT_THIS_ELE).AND.(IMP_P_ROT.EQ.0)) THEN

! calculate DMATINV...
         do ILOC=1,NLOC
            NODDG=(ELE-1)*NLOC+ILOC
            MASS2(ILOC,ILOC)=MASS2(ILOC,ILOC)+DGDM1(NODDG)
         END DO
           CALL MATDMATINV(MASS2,MASS2INV,NLOC)
           ! calculate M^{-1} C
           CALL ABMATRIXMUL(MINVC1,MASS2INV,NLOC,NLOC,&
     &                             C1T2TRAN,NLOC,MLOC)
! calculate C^T M^{-1} C
           CALL ABMATRIXMUL(CMCLOC1,C1T2,MLOC,NLOC,&
     &                              MINVC1,NLOC,MLOC)
     
           do ILOC=1,NLOC
              NODDG=(ELE-1)*NLOC+ILOC
              MASS2(ILOC,ILOC)=MASS2(ILOC,ILOC)+DGDM2(NODDG)
           END DO
           CALL MATDMATINV(MASS2,MASS2INV,NLOC)
! calculate M^{-1} C
           CALL ABMATRIXMUL(MINVC2,MASS2INV,NLOC,NLOC,&
     &                             C2T2TRAN,NLOC,MLOC)
! calculate C^T M^{-1} C
           CALL ABMATRIXMUL(CMCLOC2,C2T2,MLOC,NLOC,&
     &                              MINVC2,NLOC,MLOC)

           do ILOC=1,NLOC
              NODDG=(ELE-1)*NLOC+ILOC
              MASS2(ILOC,ILOC)=MASS2(ILOC,ILOC)+DGDM3(NODDG)
           END DO
           CALL MATDMATINV(MASS2,MASS2INV,NLOC)
! calculate M^{-1} C
           CALL ABMATRIXMUL(MINVC3,MASS2INV,NLOC,NLOC,&
     &                             C3T2TRAN,NLOC,MLOC)
! calculate C^T M^{-1} C
           CALL ABMATRIXMUL(CMCLOC3,C3T2,MLOC,NLOC,&
     &                              MINVC3,NLOC,MLOC)
           CMCLOC=CMCLOC1+CMCLOC2+CMCLOC3
        ELSE
! now for the rotation approach...
           RLOC=0.0
           do ILOC=1,NLOC
              NODDG=(ELE-1)*NLOC+ILOC
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
              IF(ROT_THIS_DGNOD(NODDG)) THEN
                 RLOC(ILOC,ILOC)              =DGT1X(NODDG)
                 RLOC(ILOC,ILOC+NLOC)         =DGT1Y(NODDG)
                 RLOC(ILOC,ILOC+2*NLOC)       =DGT1Z(NODDG)
                 RLOC(ILOC+NLOC,ILOC)         =DGT2X(NODDG)
                 RLOC(ILOC+NLOC,ILOC+NLOC)    =DGT2Y(NODDG)
                 RLOC(ILOC+NLOC,ILOC+2*NLOC)  =DGT2Z(NODDG)
                 RLOC(ILOC+2*NLOC,ILOC)       =DGNX(NODDG)
                 RLOC(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(NODDG)
                 RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(NODDG)
              ELSE
                 RLOC(ILOC,ILOC)              =1.0
                 RLOC(ILOC+NLOC,ILOC+NLOC)    =1.0
                 RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=1.0
              ENDIF
           END DO
           do ILOC3=1,BIG_NLOC
              do JLOC3=1,BIG_NLOC
                 RTLOC(ILOC3,JLOC3)=RLOC(JLOC3,ILOC3)
              END DO
           END DO
           MASSBLOC=0.0
           MASSBLOC(1:NLOC,1:NLOC)                  =MASS2(1:NLOC,1:NLOC)
           MASSBLOC(NLOC+1:2*NLOC,NLOC+1:2*NLOC)    =MASS2(1:NLOC,1:NLOC)
           MASSBLOC(2*NLOC+1:3*NLOC,2*NLOC+1:3*NLOC)=MASS2(1:NLOC,1:NLOC)
           IF(IMP_P_ROT.EQ.1) THEN
! Add in Corriolis...
             DO ILOC=1,NLOC
               MASSBLOC(ILOC,       ILOC)        =MASSBLOC(ILOC,         ILOC)     + DT*COMAT11(ELE,ILOC)
               MASSBLOC(ILOC+NLOC,  ILOC+NLOC)   =MASSBLOC(ILOC+NLOC,ILOC+NLOC)    + DT*COMAT22(ELE,ILOC)
               MASSBLOC(ILOC+2*NLOC,ILOC+2*NLOC) =MASSBLOC(ILOC+2*NLOC,ILOC+2*NLOC)+ DT*COMAT33(ELE,ILOC)
               MASSBLOC(ILOC,       ILOC+NLOC)   =DT*COMAT12(ELE,ILOC)
               MASSBLOC(ILOC,       ILOC+2*NLOC) =DT*COMAT13(ELE,ILOC)
               MASSBLOC(ILOC+NLOC,  ILOC)        =DT*COMAT21(ELE,ILOC)
               MASSBLOC(ILOC+NLOC,  ILOC+2*NLOC) =DT*COMAT23(ELE,ILOC)
               MASSBLOC(ILOC+2*NLOC,ILOC)      =DT*COMAT31(ELE,ILOC)
               MASSBLOC(ILOC+2*NLOC,ILOC+NLOC) =DT*COMAT32(ELE,ILOC)
             END DO

           ENDIF
! calculate M^{-1} R
           CALL ABMATRIXMULTI(MRLOC,MASSBLOC,RLOC,BIG_NLOC)
! calculate R^T M^{-1} R
           CALL ABMATRIXMULTI(RTMR,RTLOC,MRLOC,BIG_NLOC)
! RTMINVR...
           do ILOC=1,NLOC
              NODDG=(ELE-1)*NLOC+ILOC
              RTMR(ILOC,ILOC)=RTMR(ILOC,ILOC)+DGDM1(NODDG)
              RTMR(ILOC+NLOC,ILOC+NLOC)=RTMR(ILOC+NLOC,ILOC+NLOC)+DGDM2(NODDG)
              RTMR(ILOC+2*NLOC,ILOC+2*NLOC)=RTMR(ILOC+2*NLOC,ILOC+2*NLOC)+DGDM3(NODDG)
           END DO
! calculate RTMINVR...
           CALL MATDMATINV(RTMR,RTMINVR,BIG_NLOC)
! calculate C^T RTMINVR C
           CT(1:MLOC,1:NLOC)         =C1T2(1:MLOC,1:NLOC)
           CT(1:MLOC,1+NLOC:2*NLOC)  =C2T2(1:MLOC,1:NLOC)
           CT(1:MLOC,1+2*NLOC:3*NLOC)=C3T2(1:MLOC,1:NLOC)
           
           CTTRAN(1:NLOC,1:MLOC)         =C1T2TRAN(1:NLOC,1:MLOC)
           CTTRAN(1+NLOC:2*NLOC,1:MLOC)  =C2T2TRAN(1:NLOC,1:MLOC)
           CTTRAN(1+2*NLOC:3*NLOC,1:MLOC)=C3T2TRAN(1:NLOC,1:MLOC)
           
           CALL ABMATRIXMUL(RTMINVRC, RTMINVR,BIG_NLOC,BIG_NLOC, &
     &                                CTTRAN,BIG_NLOC,MLOC)
           CALL ABMATRIXMUL(CMCLOC,   CT,MLOC,BIG_NLOC,&
     &                                RTMINVRC,BIG_NLOC,MLOC)
        ENDIF

! Place into matrix - new CMC has same sparcity as BIGM.         
        do ILOC=1,MLOC
           GLOBI=PNDGLN((ELE-1)*MLOC+ILOC)
           do JLOC=1,MLOC
              GLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
              CALL POSINMAT(COUNT,GLOBI,GLOBJ,FREDOP,FINCMC,COLCMC,NCMC)

              DGCMC(COUNT)=DGCMC(COUNT)+CMCLOC(ILOC,JLOC)

           END DO
        END DO

     END DO

     IF((NDPSET.GE.1).AND.(NDPSET.LE.FREDOP)) THEN 
! SET PRSEEURE NODE NDPSET TO ZERO.
        DGCMC(MIDCMC(NDPSET)) =INFINY
     ENDIF

     ewrite(1,*)'exiting GETDGCMC'

   END SUBROUTINE GETDGCMC

   SUBROUTINE GETDG_FROM_SUB(UDG,VDG,WDG,U,V,W,USUB,VSUB,WSUB,&
     &                           NONODS,NLOC,TOTELE,NDGLNO)
! Calculate DG vels UDG,VDG,WDG from SGS soln
     INTEGER NONODS,NLOC,TOTELE,NDGLNO(NLOC*TOTELE)
     REAL UDG(TOTELE*NLOC),VDG(TOTELE*NLOC),WDG(TOTELE*NLOC)
     REAL U(NONODS),V(NONODS),W(NONODS)
     REAL USUB(TOTELE*NLOC),VSUB(TOTELE*NLOC),WSUB(TOTELE*NLOC)
! Local variables...
     INTEGER ELE,ILOC,NODI,DGNODI
     do ELE=1,TOTELE
        do ILOC=1,NLOC
           NODI=NDGLNO((ELE-1)*NLOC+ILOC)
           DGNODI=(ELE-1)*NLOC+ILOC
           UDG(DGNODI)=U(NODI)+USUB(DGNODI)
           VDG(DGNODI)=V(NODI)+VSUB(DGNODI)
           WDG(DGNODI)=W(NODI)+WSUB(DGNODI)
        END DO
     END DO

   END SUBROUTINE GETDG_FROM_SUB

   SUBROUTINE DGCTMULT(RHSVEC,ZERO_RHS,PNDGLN,MLOC,NLOC,TOTELE,FREDOP,&
     &             DT,C1TDGELE,C2TDGELE,C3TDGELE,&
     &             UDG,VDG,WDG)
! form r.h.s of pressure eqn...
!     RHS=RHS-(C1TU+C2TV+C3TW)/DT
     INTEGER MLOC,NLOC,TOTELE,FREDOP,PNDGLN(TOTELE*MLOC)
     REAL RHSVEC(FREDOP)
     LOGICAL ZERO_RHS
     REAL DT
     REAL C1TDGELE(TOTELE,MLOC,NLOC),C2TDGELE(TOTELE,MLOC,NLOC)
     REAL C3TDGELE(TOTELE,MLOC,NLOC)
     REAL UDG(TOTELE*NLOC),VDG(TOTELE*NLOC),WDG(TOTELE*NLOC)
! Local variables...
     INTEGER ELE,ILOC,JLOC,PNODI
     IF(ZERO_RHS) RHSVEC(1:FREDOP)=0.0
     do ELE=1,TOTELE
        do ILOC=1,MLOC
           PNODI=PNDGLN((ELE-1)*MLOC+ILOC)
           do JLOC=1,NLOC
              RHSVEC(PNODI)=RHSVEC(PNODI) &
     &      -    (C1TDGELE(ELE,ILOC,JLOC)*UDG((ELE-1)*NLOC+JLOC)&
     &          + C2TDGELE(ELE,ILOC,JLOC)*VDG((ELE-1)*NLOC+JLOC)&
     &          + C3TDGELE(ELE,ILOC,JLOC)*WDG((ELE-1)*NLOC+JLOC) )/DT
           END DO
        END DO
     END DO

   END SUBROUTINE DGCTMULT

   SUBROUTINE GET_SGS_PROJ_VELS(U,V,W,USUB,VSUB,WSUB, DP,&
     &         BIG_NLOC,NLOC,TOTELE,FREDOP,PNDGLN,MLOC,&
     &         DT,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,&
     &         IMP_P_ROT,&
     &         IMP_P_ROT_UDG,NBIGM,BIGM,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,&
! Extra variables used in GET_USUB_FROM_UDG...
     &         DGDM1,DGDM2,DGDM3,&
     &         NDGLNO,&
     &         FINDRM,COLM,CENTRM,NCOLM,&
     &         NNODP,PARA,halo_tag )
! Update the velocity U,V,W with the corrected pressure...
! Solve DG mass matrix eqn for velocity update M_DG dU_DG=dt*C*dP
! then choose how to distribute dU_DG between global and SGS...
     INTEGER BIG_NLOC,NLOC,TOTELE,MLOC
     INTEGER PNDGLN(TOTELE*MLOC)
     INTEGER NONODS,FREDOP
     REAL U(NONODS),V(NONODS),W(NONODS)
     REAL USUB(TOTELE*NLOC),VSUB(TOTELE*NLOC),WSUB(TOTELE*NLOC)
     REAL DP(FREDOP)
! If ROTAT then it is assumed that C1T etc already conbtain the 
! rotation matrices via C1T=R^T*C1T
! This sub obtains the pressure matrices CMC & CMCP(which is non-sym). 
     REAL C1TDGELE(TOTELE,MLOC,NLOC),C2TDGELE(TOTELE,MLOC,NLOC)
     REAL C3TDGELE(TOTELE,MLOC,NLOC)
     REAL MASSDGELE(TOTELE,NLOC,NLOC)
     INTEGER IMP_P_ROT
     INTEGER IMP_P_ROT_UDG,NBIGM
     REAL BIGM(IMP_P_ROT_UDG*NBIGM)
         REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC),&
     &        COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC),&
     &        COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
! NDVOLF=volume frac of all phases. 
! Put volume fraction of this phase into RMEM(1) & make a no not 
! too close to 0.0 if not multi-phase then set volume frac=1.
! NB NMPDIA=NPHASE*NPHASE*NONODS*NDIM*NDIM.
! NMXBD=maximum value of NDIM*NPHASE
! IF GETINV then get the inverse of MPDIAG with b.c's in. 

     REAL DT
     LOGICAL D3
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY     

     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     LOGICAL ROT_THIS_DGNOD(TOTELE*NLOC)
     LOGICAL ROTATDG
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
     INTEGER NDGLNO(TOTELE*NLOC)
     INTEGER NOBCU,NOBCV,NOBCW
     INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
     REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
! Extra variables used in GET_USUB_FROM_UDG...
     INTEGER NCOLM
     INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
     INTEGER NNODP,PARA,MISOUT,halo_tag
! Local variables...
     REAL, ALLOCATABLE, DIMENSION(:)::COR_DGU
     REAL, ALLOCATABLE, DIMENSION(:)::COR_DGV
     REAL, ALLOCATABLE, DIMENSION(:)::COR_DGW
       
     ALLOCATE(COR_DGU(TOTELE*NLOC))
     ALLOCATE(COR_DGV(TOTELE*NLOC))
     ALLOCATE(COR_DGW(TOTELE*NLOC))
       
! Does U+USUB satisfy cty. 
     ewrite(1,*) 'in subroutine get_sgs_proj_vels' 

     CALL CHECK_SUBDIVQ(U,V,W,USUB,VSUB,WSUB,&
     &         MLOC,NLOC,TOTELE,NDGLNO,PNDGLN,NONODS,FREDOP,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE)

! STARTSolve M_DG dDG =dt*C*dP ************************************
     CALL GET_DG_VEL_FROM_DP(COR_DGU,COR_DGV,COR_DGW, DP,&
     &         BIG_NLOC,NLOC,TOTELE,FREDOP,PNDGLN,MLOC,&
     &         DT,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,&
     &         IMP_P_ROT,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
! Th following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3 )
! END  Solve M_DG dUDG =dt*C*dP ************************************

! Choose how to distribute dU_DG between global and SGS ********
     CALL GET_USUB_FROM_UDG(U,V,W,USUB,VSUB,WSUB, &
     &         COR_DGU,COR_DGV,COR_DGW,&
     &         NLOC,TOTELE,NDGLNO,&
     &         NONODS,FINDRM,COLM,CENTRM,NCOLM,&
     &         DT,D3,&
     &         MASSDGELE,&
     &         IMP_P_ROT_UDG,NBIGM,BIGM,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &         NNODP,PARA,halo_tag,&
     &         NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2)
     
! Does U+USUB satisfy cty. 

     CALL CHECK_SUBDIVQ(U,V,W,USUB,VSUB,WSUB,&
     &         MLOC,NLOC,TOTELE,NDGLNO,PNDGLN,NONODS,FREDOP,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE)

   END SUBROUTINE GET_SGS_PROJ_VELS

   SUBROUTINE CHECK_SUBDIVQ(U,V,W,USUB,VSUB,WSUB,&
     &         MLOC,NLOC,TOTELE,NDGLNO,PNDGLN,NONODS,FREDOP,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE) 
! Does U+USUB satisfy cty - PRINT VALUE
     INTEGER BIG_NLOC,NLOC,TOTELE,MLOC
     INTEGER PNDGLN(TOTELE*MLOC)
     INTEGER NONODS,FREDOP
     REAL U(NONODS),V(NONODS),W(NONODS)
     REAL USUB(TOTELE*NLOC),VSUB(TOTELE*NLOC),WSUB(TOTELE*NLOC)
! If ROTAT then it is assumed that C1T etc already conbtain the 
! rotation matrices via C1T=R^T*C1T
! This sub obtains the pressure matrices CMC & CMCP(which is non-sym). 
     REAL C1TDGELE(TOTELE,MLOC,NLOC),C2TDGELE(TOTELE,MLOC,NLOC)
     REAL C3TDGELE(TOTELE,MLOC,NLOC)
     REAL MASSDGELE(TOTELE,NLOC,NLOC)
     INTEGER NDGLNO(TOTELE*NLOC)
! Does U+USUB satisfy cty. 
     INTEGER ELE,ILOC,GLOBI,PGLOBI,JLOC,DGNODI,DGNODJ
     REAL, ALLOCATABLE, DIMENSION(:)::UDG
     REAL, ALLOCATABLE, DIMENSION(:)::VDG
     REAL, ALLOCATABLE, DIMENSION(:)::WDG
     REAL, ALLOCATABLE, DIMENSION(:)::DISDIVQ
     REAL, ALLOCATABLE, DIMENSION(:)::PML
     ewrite(1,*) 'in subroutine check_subdivq' 
     ALLOCATE(UDG(TOTELE*NLOC))
     ALLOCATE(VDG(TOTELE*NLOC))
     ALLOCATE(WDG(TOTELE*NLOC))
     do ELE=1,TOTELE
        do ILOC=1,NLOC
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           DGNODI=(ELE-1)*NLOC+ILOC
           UDG(DGNODI)=USUB(DGNODI)+U(GLOBI)
           VDG(DGNODI)=VSUB(DGNODI)+V(GLOBI)
           WDG(DGNODI)=WSUB(DGNODI)+W(GLOBI)
        END DO
     END DO
     ALLOCATE(DISDIVQ(FREDOP))
     ALLOCATE(PML(FREDOP))
     DISDIVQ=0.0
     PML=0.0
     do ELE=1,TOTELE
        do ILOC=1,MLOC
           PGLOBI=PNDGLN((ELE-1)*MLOC+ILOC)
           do JLOC=1,NLOC
              DGNODJ=(ELE-1)*NLOC+JLOC
              DISDIVQ(PGLOBI)=DISDIVQ(PGLOBI)&
     &          +C1TDGELE(ELE,ILOC,JLOC)*UDG(DGNODJ)&
     &          +C2TDGELE(ELE,ILOC,JLOC)*VDG(DGNODJ)&
     &          +C3TDGELE(ELE,ILOC,JLOC)*WDG(DGNODJ)
              IF(FREDOP.EQ.NONODS) PML(PGLOBI)=PML(PGLOBI)+MASSDGELE(ELE,ILOC,JLOC)

           END DO
        END DO
     END DO
     IF(FREDOP.EQ.NONODS) THEN
        do PGLOBI=1,FREDOP
           DISDIVQ(PGLOBI)=DISDIVQ(PGLOBI)/PML(PGLOBI)
        END DO
     ELSE
        ewrite(1,*) 'DISCRETISED DIVQ:'
     ENDIF
     CALL PMINMX(DISDIVQ,FREDOP,'DISDIVQ*******  ')

   END SUBROUTINE CHECK_SUBDIVQ

   SUBROUTINE GET_DG_VEL_FROM_DP(COR_DGU,COR_DGV,COR_DGW, DP,&
     &         BIG_NLOC,NLOC,TOTELE,FREDOP,PNDGLN,MLOC,&
     &         DT,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE,&
     &         NONODS,&
     &         IMP_P_ROT,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
! Th following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3 )
! STARTSolve M_DG dDG =dt*C*dP ************************************
     INTEGER BIG_NLOC,NLOC,TOTELE,MLOC
     INTEGER PNDGLN(TOTELE*MLOC)
     INTEGER FREDOP
     INTEGER NONODS
     INTEGER IMP_P_ROT
         REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC),&
     &        COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC),&
     &        COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
     REAL COR_DGU(TOTELE*NLOC),COR_DGV(TOTELE*NLOC),COR_DGW(TOTELE*NLOC)
     REAL DP(FREDOP)
! If ROTAT then it is assumed that C1T etc already conbtain the 
! rotation matrices via C1T=R^T*C1T
! This sub obtains the pressure matrices CMC & CMCP(which is non-sym). 
     REAL C1TDGELE(TOTELE,MLOC,NLOC),C2TDGELE(TOTELE,MLOC,NLOC)
     REAL C3TDGELE(TOTELE,MLOC,NLOC)
     REAL MASSDGELE(TOTELE,NLOC,NLOC)
! NDVOLF=volume frac of all phases. 
! Put volume fraction of this phase into RMEM(1) & make a no not 
! too close to 0.0 if not multi-phase then set volume frac=1.
! NB NMPDIA=NPHASE*NPHASE*NONODS*NDIM*NDIM.
! NMXBD=maximum value of NDIM*NPHASE
! IF GETINV then get the inverse of MPDIAG with b.c's in. 

     REAL DT
     LOGICAL D3
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY     

     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     LOGICAL ROT_THIS_DGNOD(TOTELE*NLOC)
     LOGICAL ROTATDG
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Local variables...
     INTEGER ELE,ILOC,JLOC,PNODJ,INODDG,JNODDG,NODDG,ILOC3,JLOC3
     REAL MASS2(NLOC,NLOC),MASS2INV(NLOC,NLOC)
     REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
     REAL MASSBLOC(BIG_NLOC,BIG_NLOC)
     REAL MASSBLOCR(BIG_NLOC,BIG_NLOC),RTMR(BIG_NLOC,BIG_NLOC)
     REAL RTMINVR(BIG_NLOC,BIG_NLOC)
     REAL MASS2X(NLOC,NLOC),MASS2INVX(NLOC,NLOC)
     REAL MASS2Y(NLOC,NLOC),MASS2INVY(NLOC,NLOC)
     REAL MASS2Z(NLOC,NLOC),MASS2INVZ(NLOC,NLOC)
     LOGICAL ROT_THIS_ELE
     REAL, ALLOCATABLE, DIMENSION(:)::RHSMASU
     REAL, ALLOCATABLE, DIMENSION(:)::RHSMASV
     REAL, ALLOCATABLE, DIMENSION(:)::RHSMASW
     REAL, ALLOCATABLE, DIMENSION(:)::COR_UDG
     REAL, ALLOCATABLE, DIMENSION(:)::COR_VDG
     REAL, ALLOCATABLE, DIMENSION(:)::COR_WDG
     
     ewrite(1,*) 'in subroutine get_dg_vel_from_dp' 
     ALLOCATE(RHSMASU(TOTELE*NLOC))
     ALLOCATE(RHSMASV(TOTELE*NLOC))
     ALLOCATE(RHSMASW(TOTELE*NLOC))
       
     ALLOCATE(COR_UDG(TOTELE*NLOC))
     ALLOCATE(COR_VDG(TOTELE*NLOC))
     ALLOCATE(COR_WDG(TOTELE*NLOC))
     RHSMASU(1:TOTELE*NLOC)=0.0
     RHSMASV(1:TOTELE*NLOC)=0.0
     RHSMASW(1:TOTELE*NLOC)=0.0
     do ELE=1,TOTELE
        do JLOC=1,MLOC
           PNODJ=PNDGLN((ELE-1)*MLOC+JLOC)
           do ILOC=1,NLOC
              INODDG=(ELE-1)*NLOC+ILOC
              RHSMASU(INODDG)=RHSMASU(INODDG) &
     &        + DT*C1TDGELE(ELE,JLOC,ILOC)*DP(PNODJ)
              RHSMASV(INODDG)=RHSMASV(INODDG) &
     &        + DT*C2TDGELE(ELE,JLOC,ILOC)*DP(PNODJ)
              RHSMASW(INODDG)=RHSMASW(INODDG) &
     &        + DT*C3TDGELE(ELE,JLOC,ILOC)*DP(PNODJ)
           END DO
        END DO
     END DO
       
! solve M_DG COR_UDG=RHSMASU
     do ELE=1,TOTELE
        ROT_THIS_ELE=.FALSE.
        do ILOC=1,NLOC
           INODDG=(ELE-1)*NLOC+ILOC
           IF(ROT_THIS_DGNOD(INODDG)) ROT_THIS_ELE=.TRUE.
        END DO
        IF(ROT_THIS_ELE.OR.(IMP_P_ROT.EQ.1)) THEN
! now for the rotation approach...
           RLOC=0.0
           do ILOC=1,NLOC
              NODDG=(ELE-1)*NLOC+ILOC
              IF(ROT_THIS_DGNOD(NODDG)) THEN
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
                 RLOC(ILOC,ILOC)       =DGT1X(NODDG)
                 RLOC(ILOC,ILOC+NLOC)  =DGT1Y(NODDG)
                 RLOC(ILOC,ILOC+2*NLOC)=DGT1Z(NODDG)
                 RLOC(ILOC+NLOC,ILOC)       =DGT2X(NODDG)
                 RLOC(ILOC+NLOC,ILOC+NLOC)  =DGT2Y(NODDG)
                 RLOC(ILOC+NLOC,ILOC+2*NLOC)=DGT2Z(NODDG)
                 RLOC(ILOC+2*NLOC,ILOC)       =DGNX(NODDG)
                 RLOC(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(NODDG)
                 RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(NODDG)
              ELSE
                 RLOC(ILOC,ILOC)              =1.0
                 RLOC(ILOC+NLOC,ILOC+NLOC)    =1.0
                 RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=1.0
              ENDIF
           END DO
           do ILOC3=1,BIG_NLOC
              do JLOC3=1,BIG_NLOC
                 RTLOC(ILOC3,JLOC3)=RLOC(JLOC3,ILOC3)
              END DO
           END DO
! calculate DMATINV...
           MASS2(1:NLOC,1:NLOC)=MASSDGELE(ELE,1:NLOC,1:NLOC)
           MASSBLOC=0.0
           MASSBLOC(1:NLOC,1:NLOC)                  =MASS2(1:NLOC,1:NLOC)
           MASSBLOC(NLOC+1:2*NLOC,NLOC+1:2*NLOC)    =MASS2(1:NLOC,1:NLOC)
           MASSBLOC(2*NLOC+1:3*NLOC,2*NLOC+1:3*NLOC)=MASS2(1:NLOC,1:NLOC)
           IF(IMP_P_ROT.EQ.1) THEN
! Add in Corriolis...
             DO ILOC=1,NLOC
               MASSBLOC(ILOC,       ILOC)      =MASSBLOC(ILOC,         ILOC)          + DT*COMAT11(ELE,ILOC)
               MASSBLOC(ILOC+NLOC,  ILOC+NLOC) =MASSBLOC(ILOC+NLOC,ILOC+NLOC)         + DT*COMAT22(ELE,ILOC)
               MASSBLOC(ILOC+2*NLOC,  ILOC+2*NLOC) =MASSBLOC(ILOC+2*NLOC,ILOC+2*NLOC) + DT*COMAT33(ELE,ILOC)
               MASSBLOC(ILOC,       ILOC+NLOC)   =DT*COMAT12(ELE,ILOC)
               MASSBLOC(ILOC,       ILOC+2*NLOC) =DT*COMAT13(ELE,ILOC)
               MASSBLOC(ILOC+NLOC,  ILOC)        =DT*COMAT21(ELE,ILOC)
               MASSBLOC(ILOC+NLOC,  ILOC+2*NLOC) =DT*COMAT23(ELE,ILOC)
               MASSBLOC(ILOC+2*NLOC,ILOC)      =DT*COMAT31(ELE,ILOC)
               MASSBLOC(ILOC+2*NLOC,ILOC+NLOC) =DT*COMAT32(ELE,ILOC)
             END DO
           ENDIF
! calculate M^{-1} R
           CALL ABMATRIXMULTI(MASSBLOCR,MASSBLOC,RLOC,BIG_NLOC)
! calculate R^T M^{-1} R
           CALL ABMATRIXMULTI(RTMR,RTLOC,MASSBLOCR,BIG_NLOC)
! RTMINVR...
           do ILOC=1,NLOC
              NODDG=(ELE-1)*NLOC+ILOC
              RTMR(ILOC,ILOC)=RTMR(ILOC,ILOC)+DGDM1(NODDG)
              RTMR(ILOC+NLOC,ILOC+NLOC)=RTMR(ILOC+NLOC,ILOC+NLOC)+DGDM2(NODDG)
              RTMR(ILOC+2*NLOC,ILOC+2*NLOC)=RTMR(ILOC+2*NLOC,ILOC+2*NLOC)+DGDM3(NODDG)
           END DO
! calculate DMATINV...
           CALL MATDMATINV(RTMR,RTMINVR,BIG_NLOC)
! Calculate COR_DGU,...
           do ILOC=1,NLOC
              INODDG=(ELE-1)*NLOC+ILOC
              COR_DGU(INODDG)=0.0
              COR_DGV(INODDG)=0.0
              COR_DGW(INODDG)=0.0
              do JLOC=1,NLOC
                 JNODDG=(ELE-1)*NLOC+JLOC
                 COR_DGU(INODDG)=COR_DGU(INODDG)&
     &                          +RTMINVR(ILOC,JLOC)*RHSMASU(JNODDG)&
     &                          +RTMINVR(ILOC,JLOC+NLOC)*RHSMASV(JNODDG)&
     &                          +RTMINVR(ILOC,JLOC+2*NLOC)*RHSMASW(JNODDG)
                 COR_DGV(INODDG)=COR_DGV(INODDG)&
     &                          +RTMINVR(ILOC+NLOC,JLOC)*RHSMASU(JNODDG)&
     &                          +RTMINVR(ILOC+NLOC,JLOC+NLOC)*RHSMASV(JNODDG)&
     &                          +RTMINVR(ILOC+NLOC,JLOC+2*NLOC)*RHSMASW(JNODDG)
                 COR_DGW(INODDG)=COR_DGW(INODDG)&
     &                          +RTMINVR(ILOC+2*NLOC,JLOC)*RHSMASU(JNODDG)&
     &                          +RTMINVR(ILOC+2*NLOC,JLOC+NLOC)*RHSMASV(JNODDG)&
     &                          +RTMINVR(ILOC+2*NLOC,JLOC+2*NLOC)*RHSMASW(JNODDG)
              END DO
           END DO
        ELSE
! No rotations needed...
           MASS2X(1:NLOC,1:NLOC)=MASSDGELE(ELE,1:NLOC,1:NLOC)
           MASS2Y=MASS2X
           MASS2Z=MASS2X
           do ILOC=1,NLOC
              MASS2X(ILOC,ILOC)=MASS2X(ILOC,ILOC)+DGDM1((ELE-1)*NLOC+ILOC)
              MASS2Y(ILOC,ILOC)=MASS2Y(ILOC,ILOC)+DGDM2((ELE-1)*NLOC+ILOC)
              MASS2Z(ILOC,ILOC)=MASS2Z(ILOC,ILOC)+DGDM3((ELE-1)*NLOC+ILOC)
           END DO
           CALL MATDMATINV(MASS2X,MASS2INVX,NLOC)
           CALL MATDMATINV(MASS2Y,MASS2INVY,NLOC)
           CALL MATDMATINV(MASS2Z,MASS2INVZ,NLOC)
           do ILOC=1,NLOC
              INODDG=(ELE-1)*NLOC+ILOC
              COR_DGU(INODDG)=0.0
              COR_DGV(INODDG)=0.0
              COR_DGW(INODDG)=0.0
              do JLOC=1,NLOC
                 JNODDG=(ELE-1)*NLOC+JLOC
                 COR_DGU(INODDG)=COR_DGU(INODDG)&
     &                          +MASS2INVX(ILOC,JLOC)*RHSMASU(JNODDG)
                 COR_DGV(INODDG)=COR_DGV(INODDG)&
     &                          +MASS2INVY(ILOC,JLOC)*RHSMASV(JNODDG)
                 COR_DGW(INODDG)=COR_DGW(INODDG)&
     &                          +MASS2INVZ(ILOC,JLOC)*RHSMASW(JNODDG)
              END DO
           END DO
        ENDIF
     END DO

! END  Solve M_DG dUDG =dt*C*dP ************************************
   END SUBROUTINE GET_DG_VEL_FROM_DP

   SUBROUTINE GET_USUB_FROM_UDG(U,V,W,USUB,VSUB,WSUB, &
     &         COR_DGU,COR_DGV,COR_DGW,&
     &         NLOC,TOTELE,NDGLNO,&
     &         NONODS,FINDRM,COLM,CENTRM,NCOLM,&
     &         DT,D3,&
     &         MASSDGELE,&
     &         IMP_P_ROT_UDG,NBIGM,BIGM,&
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &         NNODP,PARA,halo_tag,&
     &         NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2)
! choose how to distribute dU_DG between global and SGS ********
     INTEGER MATSTR
     REAL INFINY
     PARAMETER(MATSTR=0,INFINY=1.0E+20)
     INTEGER NLOC,TOTELE
     INTEGER NDGLNO(TOTELE*NLOC)
     INTEGER NONODS,NCOLM
     INTEGER FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
     REAL U(NONODS),V(NONODS),W(NONODS)
     REAL USUB(TOTELE*NLOC),VSUB(TOTELE*NLOC),WSUB(TOTELE*NLOC)
     REAL COR_DGU(TOTELE*NLOC),COR_DGV(TOTELE*NLOC),COR_DGW(TOTELE*NLOC)
     REAL MASSDGELE(TOTELE,NLOC,NLOC)
     INTEGER IMP_P_ROT_UDG,NBIGM
     REAL BIGM(IMP_P_ROT_UDG*NBIGM)
         REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC),&
     &        COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC),&
     &        COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
     REAL DT
     LOGICAL D3
     INTEGER NNODP,PARA,halo_tag
     INTEGER NOBCU,NOBCV,NOBCW
     INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
     REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
! Local variables...
     INTEGER IDUM(1),COUNT,GLOBI,GLOBJ,ILOC,JLOC,ELE,DGNODI,DGNODJ
     INTEGER II,NOD,ILENG,NR2
     INTEGER KITS,IZER
     REAL RDUM(1)
     REAL VLM
     REAL, ALLOCATABLE, DIMENSION(:)::COR_U
     REAL, ALLOCATABLE, DIMENSION(:)::COR_V
     REAL, ALLOCATABLE, DIMENSION(:)::COR_W
     REAL, ALLOCATABLE, DIMENSION(:)::COR_USUB
     REAL, ALLOCATABLE, DIMENSION(:)::COR_VSUB
     REAL, ALLOCATABLE, DIMENSION(:)::COR_WSUB
     REAL, ALLOCATABLE, DIMENSION(:)::MASMAT
     REAL, ALLOCATABLE, DIMENSION(:)::VECSECU
     REAL, ALLOCATABLE, DIMENSION(:)::VECSECV
     REAL, ALLOCATABLE, DIMENSION(:)::VECSECW
     REAL, ALLOCATABLE, DIMENSION(:)::MAS_DIA
     REAL, ALLOCATABLE, DIMENSION(:)::UVW
     REAL, ALLOCATABLE, DIMENSION(:)::VECXYZ
       
     character(len=OPTION_PATH_LEN) :: velocity_option_path
     ewrite(1,*) 'in subroutine get_usub_from_udg' 
       
       ! won't work with multi-phase/material
     velocity_option_path='/material_phase[0]/vector_field::Velocity'
       
     ALLOCATE(COR_U(NONODS))
     ALLOCATE(COR_V(NONODS))
     ALLOCATE(COR_W(NONODS))
     COR_U=0.0
     COR_V=0.0
     COR_W=0.0

     ALLOCATE(COR_USUB(TOTELE*NLOC))
     ALLOCATE(COR_VSUB(TOTELE*NLOC))
     ALLOCATE(COR_WSUB(TOTELE*NLOC))

     ALLOCATE(VECSECU(NONODS))
     VECSECU=0.0
     ALLOCATE(VECSECV(NONODS))
     VECSECV=0.0
     ALLOCATE(VECSECW(NONODS))
     VECSECW=0.0
       
     IZER=0
       
     IF(IMP_P_ROT_UDG.EQ.1) THEN
          
        BIGM=0.0
        do ELE=1,TOTELE
           do ILOC=1,NLOC
              GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
              DGNODI=(ELE-1)*NLOC+ILOC
              do JLOC=1,NLOC
                 GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
                 DGNODJ=(ELE-1)*NLOC+JLOC
                 VLM=MASSDGELE(ELE,ILOC,JLOC)
             
                 VECSECU(GLOBI)=VECSECU(GLOBI)+VLM*COR_DGU(DGNODJ)
                 VECSECV(GLOBI)=VECSECV(GLOBI)+VLM*COR_DGV(DGNODJ)
                 VECSECW(GLOBI)=VECSECW(GLOBI)+VLM*COR_DGW(DGNODJ)
             
                 CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
                 BIGM(COUNT)=BIGM(COUNT)                 + VLM
                 BIGM(COUNT+4*NCOLM)=BIGM(COUNT+4*NCOLM) + VLM
                 BIGM(COUNT+8*NCOLM)=BIGM(COUNT+8*NCOLM) + VLM
              END DO
! Add Coriolis in... 
           VECSECU(GLOBI)=VECSECU(GLOBI) &
     &       +DT*COMAT11(ELE,ILOC)*COR_DGU(DGNODI)+DT*COMAT12(ELE,ILOC)*COR_DGV(DGNODI)&
     &       +DT*COMAT13(ELE,ILOC)*COR_DGW(DGNODI)
           VECSECV(GLOBI)=VECSECV(GLOBI) &
     &       +DT*COMAT21(ELE,ILOC)*COR_DGU(DGNODI)+DT*COMAT22(ELE,ILOC)*COR_DGV(DGNODI)&
     &       +DT*COMAT23(ELE,ILOC)*COR_DGW(DGNODI)
           VECSECW(GLOBI)=VECSECW(GLOBI) &
     &       +DT*COMAT31(ELE,ILOC)*COR_DGU(DGNODI)+DT*COMAT32(ELE,ILOC)*COR_DGV(DGNODI)&
     &       +DT*COMAT33(ELE,ILOC)*COR_DGW(DGNODI)
           COUNT=CENTRM(GLOBI)
           BIGM(COUNT)        =BIGM(COUNT)         + DT*COMAT11(ELE,ILOC)
           BIGM(COUNT+NCOLM)  =BIGM(COUNT+NCOLM)   + DT*COMAT12(ELE,ILOC)
           BIGM(COUNT+2*NCOLM)=BIGM(COUNT+2*NCOLM) + DT*COMAT13(ELE,ILOC)
           
           BIGM(COUNT+3*NCOLM)=BIGM(COUNT+3*NCOLM) + DT*COMAT21(ELE,ILOC)
           BIGM(COUNT+4*NCOLM)=BIGM(COUNT+4*NCOLM) + DT*COMAT22(ELE,ILOC)
           BIGM(COUNT+5*NCOLM)=BIGM(COUNT+5*NCOLM) + DT*COMAT23(ELE,ILOC)
           
           BIGM(COUNT+6*NCOLM)=BIGM(COUNT+6*NCOLM) + DT*COMAT31(ELE,ILOC)
           BIGM(COUNT+7*NCOLM)=BIGM(COUNT+7*NCOLM) + DT*COMAT32(ELE,ILOC)
           BIGM(COUNT+8*NCOLM)=BIGM(COUNT+8*NCOLM) + DT*COMAT33(ELE,ILOC)
           END DO
        END DO
! put b.c's in:
        do II=1,NOBCU
           NOD=BCU2(II)
           BIGM(CENTRM(NOD))        =INFINY
        END DO
        do II=1,NOBCV
           NOD=BCV2(II)
           BIGM(CENTRM(NOD)+4*NCOLM)=INFINY
        END DO
        do II=1,NOBCW
           NOD=BCW2(II)
           BIGM(CENTRM(NOD)+8*NCOLM)=INFINY
        END DO

! Solve amended mass matrix with Corriolis in... 
        ILENG=3*NONODS
        ALLOCATE(UVW(ILENG))
        UVW=0.0
        ALLOCATE(VECXYZ(ILENG))
        VECXYZ(1:NONODS)=VECSECU(1:NONODS)
        VECXYZ(NONODS+1:2*NONODS)=VECSECV(1:NONODS)
        VECXYZ(2*NONODS+1:3*NONODS)=VECSECW(1:NONODS)
        CALL GMRES(UVW,VECXYZ,NONODS,ILENG,NNODP,.TRUE.,&
     &       BIGM,FINDRM,COLM,NCOLM,NBIGM,&
     & halo_tag,KITS)
! UVW contains the acceleration.
        COR_U(1:NONODS)=UVW(1:NONODS)           
        COR_V(1:NONODS)=UVW(NONODS+1:2*NONODS)
        COR_W(1:NONODS)=UVW(2*NONODS+1:3*NONODS)
        DEALLOCATE(UVW)
        DEALLOCATE(VECXYZ)

! ENDOF (FOR ELSE) IF(IMP_P_ROT_UDG.EQ.1) THEN ...
     ELSE
       
        ALLOCATE(MASMAT(NCOLM))
        MASMAT=0.0
       
        do ELE=1,TOTELE
           do ILOC=1,NLOC
              GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
              do JLOC=1,NLOC
                 GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
                 DGNODJ=(ELE-1)*NLOC+JLOC
                 VLM=MASSDGELE(ELE,ILOC,JLOC)
                 VECSECU(GLOBI)=VECSECU(GLOBI)+VLM*COR_DGU(DGNODJ)
                 VECSECV(GLOBI)=VECSECV(GLOBI)+VLM*COR_DGV(DGNODJ)
                 VECSECW(GLOBI)=VECSECW(GLOBI)+VLM*COR_DGW(DGNODJ)
                 CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
                 MASMAT(COUNT)=MASMAT(COUNT)+VLM
              END DO
           END DO
        END DO
              
        ALLOCATE(MAS_DIA(NONODS))
       
        do NOD=1,NONODS
           MAS_DIA(NOD)=MASMAT(CENTRM(NOD))
        END DO
        do II=1,NOBCU
           NOD=BCU2(II)
           MASMAT(CENTRM(NOD))=INFINY
        END DO

        CALL SOLCG(COR_U,VECSECU,NONODS,NONODS,NNODP,.TRUE.,&
     &                 MASMAT,FINDRM,COLM,NCOLM, NCOLM,&
     &                 halo_tag, KITS, &
     &                 option_path=velocity_option_path)
     
        do NOD=1,NONODS
           MASMAT(CENTRM(NOD))=MAS_DIA(NOD)
        END DO
        do II=1,NOBCV
           NOD=BCV2(II)
           MASMAT(CENTRM(NOD))=INFINY
        END DO

        CALL SOLCG(COR_V,VECSECV,NONODS,NONODS,NNODP,.TRUE.,&
     &                 MASMAT,FINDRM,COLM,NCOLM, NCOLM,&
     &                 halo_tag, KITS, &
     &                 option_path=velocity_option_path)
     
        do NOD=1,NONODS
           MASMAT(CENTRM(NOD))=MAS_DIA(NOD)
        END DO
        do II=1,NOBCW
           NOD=BCW2(II)
           MASMAT(CENTRM(NOD))=INFINY
        END DO

        CALL SOLCG(COR_W,VECSECW,NONODS,NONODS,NNODP,.TRUE.,&
     &                 MASMAT,FINDRM,COLM,NCOLM, NCOLM,&
     &                 halo_tag, KITS, &
     &                 option_path=velocity_option_path)
     
! ENDOF IF(IMP_P_ROT_UDG.EQ.1) THEN ELSE...
     ENDIF
              
     do ELE=1,TOTELE
        do ILOC=1,NLOC
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           DGNODI=(ELE-1)*NLOC+ILOC
           COR_USUB(DGNODI)=COR_DGU(DGNODI)-COR_U(GLOBI)
           COR_VSUB(DGNODI)=COR_DGV(DGNODI)-COR_V(GLOBI)
           COR_WSUB(DGNODI)=COR_DGW(DGNODI)-COR_W(GLOBI)
        END DO
     END DO
! Update velocities...  
     U=U+COR_U 
     V=V+COR_V 
     W=W+COR_W 
     USUB=USUB+COR_USUB 
     VSUB=VSUB+COR_VSUB 
     WSUB=WSUB+COR_WSUB 
     CALL PMINMX(COR_U,NONODS,'COR_U*******  ')
     CALL PMINMX(COR_V,NONODS,'COR_V*******  ')
     CALL PMINMX(COR_W,NONODS,'COR_W*******  ')
         
     CALL PMINMX(COR_DGU,TOTELE*NLOC,'COR_DGU*******  ')
     CALL PMINMX(COR_DGV,TOTELE*NLOC,'COR_DGV*******  ')
     CALL PMINMX(COR_DGW,TOTELE*NLOC,'COR_DGW*******  ')
         
     CALL PMINMX(COR_USUB,TOTELE*NLOC,'COR_USUB*******  ')
     CALL PMINMX(COR_VSUB,TOTELE*NLOC,'COR_VSUB*******  ')
     CALL PMINMX(COR_WSUB,TOTELE*NLOC,'COR_WSUB*******  ')

   END SUBROUTINE GET_USUB_FROM_UDG




   SUBROUTINE DIFF3DSGS(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD,&
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA,&
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB,&
     &     COMAT11,COMAT12,COMAT13,COMAT21,&
     &     COMAT22,COMAT23,COMAT31,COMAT32,COMAT33,  &
     &     BOYMA11,BOYMA12,BOYMA13,BOYMA21,&
     &     BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33, &
     &     VECX,VECY,VECZ, &
     &     BIGM, ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     nnodp,para,halo_tag,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBTLOC,NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp, &
     &     CMAT_STORE,DMAT_STORE,&
     &     DINVSOUSUBU_STORE,DINVSOUSUBV_STORE,&
     &     DINVSOUSUBW_STORE,&
     &     MASSDGELE,C1TDGELE,C2TDGELE,C3TDGELE,P,&
     &     FTHETA,GRAVTY,IUNST_FREE,G1TDGELE,G2TDGELE,G3TDGELE,POLD,&
     &     DGCMC,NCMC,FINCMC,COLCMC,MIDCMC, &
     &     SUF_TOP_OCEAN,&
     &     ROTATDG,ROT_THIS_DGNOD,&
     &     DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &     DGDM1,DGDM2,DGDM3,&
     &     PRES_RHS, velocity_option_path)
      
!     This subroutine discretises a field equation using Galerkin Least squares. 
!     The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
!     
! All methods use a balancing diffusion term and THETA time stepping. 
! The following are absolute values of DISOPT...
! Galerkin (146 is defaulted to if non of the below are used): 
! 154 same as 147 and no added viscocity.
! 153 same as 152 and no added viscocity.
! 152 same as 147 but reduced magnitude by a factor of 8.
! 151 is LES applied to the global system. 
! 150 is LES applied to the local Residual system (involving pertation soln) only RECOMMENDED
! 149 is LES applied to the local system only RECOMMENDED
! 148 is LES applied to the global system only 
! 147 is LES applied to every thing. 
! 145 PG applied to just the global part and using global velocity.
! 144 PG applied to just the global part and using the full velocity.
      
!     if DISOPN.ne.0 then set advection to zero and treat advection 
!     using the high res method.
!     IF LUMP then lump the mass matrix else dont. 
!     BETA controls the conservative discretisation 
!     BETA=0. is the usual advection form, BETA=1. is divergence form. 
!     BETA=0.5 is an average of the divergence and advection form.(recommended).
!     ABSORB =absorbtion of energy etc. 
!     MUPTXX,MUPTYY =components of diffusivities. 
!     ML2MXX,ML2MXY etc contain the length scales for LES. 
!     NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
!     UG,VG,WG are the grid velocities. 
!     NDGLNO=element pter for unknowns. 
!     SONDGL=element pter for materials and sources. 
!     VONDGL=element pter for advection NU,UG.
!     XONDGL=element pter for coordinates(x,y,z)
!     If INMAT then return matricie(s). 
!     R0=the minimum distance from ocean bed to centre of earth.
!     D0=H/r0 where H is the height of ocean above R0 (H not used here)
!     For dimensional run D0=0.(R0 is not used used here)
!     NB This sub is only used for field equations in ocean modelling. 
!     ------------------------------------------------------------------------
!     If we are solving for another variable like temperature 
!     or chemical species then NU,NV,NW will be the velocities 
!     and U=TEMPERATURE OR CHEMICAL SPECIES. 
!     ------------------------------------------------------------------------
!     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
     REAL RWIND,INFINY
     PARAMETER(RWIND=0.5,INFINY=1.0E+20)
      LOGICAL VERTICAL_BOY 
      PARAMETER(VERTICAL_BOY=.true.)
! if VERTICAL_BOY only add in a vertical bouyancy term.. which acts like absorption. 
     INTEGER SNONOD,VNONOD,XNONOD,NCT,NCOLM,NLOC,BIG_NLOC,NGI
     INTEGER MLOC,DISOPT2,DISOPT,DISOPN,DISCRA
     REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
     REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
     REAL DENPT(SNONOD)
     REAL ABSORB(SNONOD)
     REAL  DT,THETA,BETA
     INTEGER TOTELE,NONODS,FREDOP
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
     REAL U(NONODS),V(NONODS),W(NONODS)
     REAL UNEW(NONODS),VNEW(NONODS),WNEW(NONODS)
     INTEGER ISUB,INUSUB,NSUBTLOC,NSUBVLOC
     REAL NUSUB(TOTELE*NSUBVLOC),NVSUB(TOTELE*NSUBVLOC)
     REAL NWSUB(TOTELE*NSUBVLOC)
     REAL USUB(TOTELE*NSUBTLOC),VSUB(TOTELE*NSUBTLOC),WSUB(TOTELE*NSUBTLOC)
     REAL UNEWSUB(TOTELE*NSUBTLOC),VNEWSUB(TOTELE*NSUBTLOC),WNEWSUB(TOTELE*NSUBTLOC)
     REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
     REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
      
! This is what we add into the SGS model - a discretised source.      
     REAL VECX_SGSADD(TOTELE*NSUBTLOC),VECY_SGSADD(TOTELE*NSUBTLOC)
     REAL VECZ_SGSADD(TOTELE*NSUBTLOC)
       INTEGER NBUOY,NSOGRASUB
! soxgra,soygra,sozgra is the direction of gravity force when geobal.le.-10.or.(HYDROS.NE.0)
       REAL SOXGRA(SNONOD),SOYGRA(SNONOD),SOZGRA(SNONOD)
       REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
       REAL SOZGRASUB(TOTELE*NSOGRASUB)
       REAL COMAT11(TOTELE,NSUBTLOC),COMAT12(TOTELE,NSUBTLOC),COMAT13(TOTELE,NSUBTLOC)
       REAL COMAT21(TOTELE,NSUBTLOC),COMAT22(TOTELE,NSUBTLOC),COMAT23(TOTELE,NSUBTLOC)
       REAL COMAT31(TOTELE,NSUBTLOC),COMAT32(TOTELE,NSUBTLOC),COMAT33(TOTELE,NSUBTLOC)
       REAL BOYMA11(NBUOY*TOTELE,NLOC),BOYMA12(NBUOY*TOTELE,NLOC),BOYMA13(NBUOY*TOTELE,NLOC)
       REAL BOYMA21(NBUOY*TOTELE,NLOC),BOYMA22(NBUOY*TOTELE,NLOC),BOYMA23(NBUOY*TOTELE,NLOC)
       REAL BOYMA31(NBUOY*TOTELE,NLOC),BOYMA32(NBUOY*TOTELE,NLOC),BOYMA33(NBUOY*TOTELE,NLOC)
     REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
     INTEGER NBIGM
     REAL BIGM(NBIGM)
     LOGICAL COGRAX
     INTEGER NOBCFSimp
     INTEGER BCT2FSimp(NOBCFSimp)
     REAL ML(NONODS)
     REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
     REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
     INTEGER ISUBSOU
     REAL SUBSOUX(ISUBSOU*TOTELE*NLOC),SUBSOUY(ISUBSOU*TOTELE*NLOC)
     REAL SUBSOUZ(ISUBSOU*TOTELE*NLOC)
     INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
     INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
     INTEGER XONDGL(TOTELE*NLOC)
      
     INTEGER FINDRM(NONODS+1),COLM(NCOLM)
     INTEGER CENTRM(NONODS)
     INTEGER FINDCT(FREDOP+1),COLCT(NCT)
     REAL CMAT_STORE(TOTELE,NLOC,NLOC),DMAT_STORE(TOTELE,NLOC,NLOC)
     REAL DINVSOUSUBU_STORE(TOTELE,NLOC),DINVSOUSUBV_STORE(TOTELE,NLOC)
     REAL DINVSOUSUBW_STORE(TOTELE,NLOC)
     REAL MASSDGELE(TOTELE,NLOC,NLOC)
     REAL C1TDGELE(TOTELE,MLOC,NLOC)
     REAL C2TDGELE(TOTELE,MLOC,NLOC),C3TDGELE(TOTELE,MLOC,NLOC)
     REAL P(FREDOP)
      INTEGER IUNST_FREE
      REAL FTHETA,GRAVTY
      REAL G1TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC)
      REAL G2TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC),G3TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC)
      REAL POLD(FREDOP*IUNST_FREE)
      INTEGER NCMC
      REAL DGCMC(NCMC*IUNST_FREE)
      INTEGER FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)
     LOGICAL SUF_TOP_OCEAN
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY     

     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
     REAL PRES_RHS(FREDOP)
     character(len=*), intent(in):: velocity_option_path

     REAL M(MLOC,NGI)
     REAL MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
     REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
     REAL WEIGHT(NGI)
     LOGICAL LUMP,UPWIND,OPWIND,XTPETV
     INTEGER PREOPT
! IF PREOPT=2 then treat every boundary like a free surface - need 
! to do this for open boundaries. 
     LOGICAL ABSLUM,SLUMP
     LOGICAL INMAT
     INTEGER VERSIO
     INTEGER nnodp,para,halo_tag
     INTEGER GEOBAL
     REAL SCFACTH0
     INTEGER SNSUBTLOC,ISUBNOD
      
     REAL MASS,MATRIX,SOUMAT,VLMGI,VLKGI,VABS
     REAL::VLM=0.0
     REAL VF1M,VF2M,VF2,VF3M
     REAL THETA1,THETA2,THETA3,AHAT1,ALPHA,ALPHA2,ALPHA3,THTWDT,INDT,TWOTHI
     REAL PE,HOVERQ,PWEIGI
     REAL GALERK,LEASQR
     REAL UD(NGI),VD(NGI),WD(NGI)
     REAL TUD(NGI),TVD(NGI),TWD(NGI)
     REAL UDDX(NGI),VDDY(NGI),WDDZ(NGI)
     REAL UDDX2(NGI),VDDY2(NGI),WDDZ2(NGI)
     REAL UGDX(NGI),VGDY(NGI),WGDZ(NGI)
     REAL MUGIXX,MUGIXY,MUGIXZ
     REAL MUGIYY,MUGIYZ,MUGIZZ
     REAL L2GIXX,L2GIXY,L2GIXZ
     REAL L2GIYY,L2GIYZ,L2GIZZ
     REAL DENGI(NGI),DENGI2(NGI),L1GI(NGI)
     REAL ABDGI(NGI)
     REAL AGI,BGI,CGI,DGI
     REAL EGI,FGI,GGI,HGI,KGI
     REAL DETJ,DETWEI(NGI)
     REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
     REAL GIVOLX,GIVOLY,GIVOLZ
!     HX,HY-characteristic length scales in x,y directions.
     REAL HXGI,HYGI,HZGI
     REAL AHAT(NGI),GAMMA(NGI),GAMMAG(NGI),MASCON(NGI)
     REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
     REAL BECTWE(NGI)

     INTEGER ELE,ILOC,JLOC,L
     INTEGER GI,COUNT,I,GLOBI,GLOBJ,GLOBIS,GLOBJS
     INTEGER POS,LOWER,UPPER,PGLOBJ

     REAL D0,DD0,INDD0
     LOGICAL NODIM
     LOGICAL GALWIN
     LOGICAL CTWIND,LES
     REAL ABLUMP,RSLUMP
     REAL R,RUD,RVD,RWD,R1,R2,RR,RMXD,RMXMXX,RMXMYY,RMXMZZ,RMXABS
     REAL::RNN=0.0
     REAL RMXSX,RMXSY,RMXSZ,RMID,RMIMXX,RMIMYY,RMIMZZ,RMIABS
     REAL RMISX,RMISY,RMISZ,RLUMP,RSYM,RCONSI,RGAMMA
     INTEGER IGLS,IGLV,IGLX,INUM
     INTEGER ISPHERE,ICENT
     REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
     REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
     REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
     REAL MYZDTW(NGI),MZZDTW(NGI)
     REAL AMXXDTW(NGI),AMXYDTW(NGI),AMXZDTW(NGI),AMYYDTW(NGI)
     REAL AMYZDTW(NGI),AMZZDTW(NGI)
     REAL MUGIXX2,MUGIXY2,MUGIXZ2
     REAL MUGIYY2,MUGIYZ2,MUGIZZ2
     REAL MUPTXXGI,MUPTXYGI,MUPTXZGI,MUPTYYGI,MUPTYZGI,MUPTZZGI
     REAL XD(NGI),YD(NGI),ZD(NGI)
     REAL VLND,VLN2,VLND2,VLN3,VLND3
     REAL TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ
     REAL L2GIXX2,L2GIXY2,L2GIXZ2,L2GIYY2,L2GIYZ2,L2GIZZ2
!  for the extra terms in the adjoint model, local variables....
     REAL VLND_UX_ADJ,VLND_UY_ADJ,VLND_UZ_ADJ,    &
     &       VLND_VX_ADJ,VLND_VY_ADJ,VLND_VZ_ADJ, &
     &       VLND_WX_ADJ,VLND_WY_ADJ,VLND_WZ_ADJ  
     REAL HINX,HINY,HINZ,VLKTEN,NLN,VLMAB
     REAL XNDMUJLOC,YNDMUJLOC,ZNDMUJLOC,XNDMUILOC,YNDMUILOC,ZNDMUILOC
     REAL TEN,TENAB,ABTEN
     REAL VLN,ADVILOC,ADVJLOC,RDIFF
     REAL UATTHETA,VATTHETA,WATTHETA
     LOGICAL BALDIF,BALGLS
     REAL L1(NGI),L2(NGI),L3(NGI),L4(NGI),PGWEIG(NGI) 
     REAL A11(NGI),A12(NGI),A13(NGI)
     REAL A21(NGI),A22(NGI),A23(NGI)
     REAL A31(NGI),A32(NGI),A33(NGI) 
     REAL GAMB(NGI),GAMBAL(NGI)
     REAL ABGI(NGI),ABSDX(NGI),ABSDY(NGI),ABSDZ(NGI),DDETWE(NGI)
     REAL AMAT(NLOC,NLOC),BMAT(NLOC,NSUBTLOC)
     REAL CMAT(NSUBTLOC,NLOC),DMAT(NSUBTLOC,NSUBTLOC)
     REAL STDMAT(NSUBTLOC,NSUBTLOC)
     REAL OMAT(NLOC,NLOC),LOMAT(NLOC,NLOC),COMAT(NLOC,NLOC)
     REAL DMATINV(NSUBTLOC,NSUBTLOC),DMATINVC(NSUBTLOC,NLOC)
     REAL BDINVC(NLOC,NLOC),BDINVMAT(NLOC,NSUBTLOC)
     REAL MATN1XN2(NLOC,NSUBTLOC),MATN1YN2(NLOC,NSUBTLOC)
     REAL MATN1ZN2(NLOC,NSUBTLOC)
     REAL MATN1XN1(NLOC,NLOC),MATN1YN1(NLOC,NLOC)
     REAL MATN1ZN1(NLOC,NLOC)
     REAL MATN2XN2(NSUBTLOC,NSUBTLOC),MATN2YN2(NSUBTLOC,NSUBTLOC)
     REAL MATN2ZN2(NSUBTLOC,NSUBTLOC)
     REAL MATN2N2X(NSUBTLOC,NSUBTLOC),MATN2N2Y(NSUBTLOC,NSUBTLOC)
     REAL MATN2N2Z(NSUBTLOC,NSUBTLOC)
     REAL MASSM1M1(NLOC,NLOC),MASSM1M2(NLOC,NSUBTLOC)
     REAL MASSM2M1(NSUBTLOC,NLOC),MASSM2M2(NSUBTLOC,NSUBTLOC)
     REAL MASSINVM2M2(NSUBTLOC,NSUBTLOC)
     REAL MASSINVM1M1(NLOC,NLOC)
     REAL KMAT22X(NSUBTLOC,NSUBTLOC),KMAT22Y(NSUBTLOC,NSUBTLOC)
     REAL KMAT22Z(NSUBTLOC,NSUBTLOC)
     REAL DIFMAT(NSUBTLOC,NSUBTLOC),DIFMATB(NLOC,NSUBTLOC)
     REAL DIFMATX(NSUBTLOC,NSUBTLOC)
     REAL DIFMATY(NSUBTLOC,NSUBTLOC),DIFMATZ(NSUBTLOC,NSUBTLOC)
     REAL DIFMATX2(NSUBTLOC,NSUBTLOC)
     REAL DIFMATY2(NSUBTLOC,NSUBTLOC),DIFMATZ2(NSUBTLOC,NSUBTLOC)
     REAL DIFMATX2B(NLOC,NSUBTLOC)
     REAL DIFMATY2B(NLOC,NSUBTLOC),DIFMATZ2B(NLOC,NSUBTLOC)
     REAL MASSLUMM1M1(NLOC,NLOC),MASSLUMM2M2(NSUBTLOC,NSUBTLOC)
     REAL SOUSUBU(NSUBTLOC),SOUSUBV(NSUBTLOC),SOUSUBW(NSUBTLOC)
     REAL DERSUBNUX(NSUBVLOC),DERSUBNVY(NSUBVLOC),DERSUBNWZ(NSUBVLOC)
     REAL SOUSUBNU(NSUBVLOC),SOUSUBNV(NSUBVLOC),SOUSUBNW(NSUBVLOC)
     REAL DMATINVCSTORE(BIG_NLOC,BIG_NLOC)
     REAL DINVSOUSUBUSTORE(BIG_NLOC),DINVSOUSUBVSTORE(BIG_NLOC)
     REAL DINVSOUSUBWSTORE(BIG_NLOC)
     REAL BIG_AMAT(BIG_NLOC,BIG_NLOC),BIG_BMAT(BIG_NLOC,BIG_NLOC)
     REAL BIG_CMAT(BIG_NLOC,BIG_NLOC),BIG_DMAT(BIG_NLOC,BIG_NLOC)
     REAL BIG_DMATINVC(BIG_NLOC,BIG_NLOC),BIG_DMATINV(BIG_NLOC,BIG_NLOC)
     REAL BIG_BDINVC(BIG_NLOC,BIG_NLOC),BIG_BDINVMAT(BIG_NLOC,BIG_NLOC)
     REAL C1TDGELELOC(MLOC,NLOC)
     REAL C2TDGELELOC(MLOC,NLOC),C3TDGELELOC(MLOC,NLOC)
     REAL G1TDGELELOC(MLOC,NLOC)
     REAL G2TDGELELOC(MLOC,NLOC),G3TDGELELOC(MLOC,NLOC)
     REAL TOTSOUX(NLOC),TOTSOUY(NLOC),TOTSOUZ(NLOC)
     INTEGER LESNOD,JNODSUB
     REAL RINMAT,RHSMAT,SRHSMAT
     REAL VLNXN,VLNYN,VLNZN
     REAL VLNNX,VLNNY,VLNNZ
     REAL VLMXN,VLMYN,VLMZN
     REAL C1TCONU,C2TCONV,C3TCONW
! Local variables...
     INTEGER SNLOC,SNGI,SNCLOC
     PARAMETER(SNLOC=3,SNCLOC=6,SNGI=7)
     REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
     REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
     REAL SNX(SNLOC,SNGI),SNY(SNLOC,SNGI),SNZ(SNLOC,SNGI)
     REAL SWEIGH(SNGI)
     REAL SNC(SNCLOC,SNGI),SNCLX(SNCLOC,SNGI),SNCLY(SNCLOC,SNGI)
     REAL SNCX(SNCLOC,SNGI),SNCY(SNCLOC,SNGI),SNCZ(SNCLOC,SNGI)
     REAL XSL(SNCLOC),YSL(SNCLOC),ZSL(SNCLOC)
     REAL SUD(SNGI),SVD(SNGI),SWD(SNGI)
     REAL SUDGLOB(SNGI),SVDGLOB(SNGI),SWDGLOB(SNGI)
     REAL NDOTQ(SNGI),NDOTQGLOB(SNGI)
     REAL SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI)
     REAL SDETWE(SNGI)
     REAL SDIRX(SNGI),SDIRY(SNGI),SDIRZ(SNGI)
     REAL SGS_UPWIND_FAC(NGI)
     INTEGER IFACE,NFACE
     PARAMETER(NFACE=4)
     INTEGER LOCLIST(NFACE,3),SILOC2ILOC(SNLOC)
     INTEGER SILOC,SJLOC,SKLOC,IGL,COL,INOD,JNOD,SGI,KLOC
     INTEGER BIG_ILOC,BIG_JLOC,ELE2,IJDISP,ID,JD,IXNOD,NODDG
     INTEGER IDGNOD,NOD,II,ICOUNT,PGLOBI
     REAL RADI,RADJ,RADP,RAD,RNU,RNV,RNW,SAREA,VOL
     REAL NORMX,NORMY,NORMZ,VLNDOTQ,allVLNDOTQ,STABINEL,UPWIN,TDIVU
     REAL VLOMEGA,OMEGAGI,STABOME
     REAL A11MAT,A12MAT,A21MAT,A22MAT,O11MAT,O12MAT,O22MAT,O21MAT
     REAL RO11MAT,RO12MAT,RO22MAT,RO21MAT
     REAL RNUSUBGI,RNVSUBGI,RNWSUBGI
     REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
     REAL NSUB(NSUBTLOC,NGI)
     REAL NSUBLX(NSUBTLOC,NGI),NSUBLY(NSUBTLOC,NGI),NSUBLZ(NSUBTLOC,NGI)
     REAL NSUBX(NSUBTLOC,NGI),NSUBY(NSUBTLOC,NGI),NSUBZ(NSUBTLOC,NGI)
     REAL AVE_SNORMXN,AVE_SNORMYN,AVE_SNORMZN,ONE_OR_ZERO,ROCEAN
     REAL LOCPNORMX,LOCPNORMY,LOCPNORMZ
     REAL XC,YC,ZC,RN,sum,stVLKTEN,stTEN
     REAL ALPHA11,ALPHA12,ALPHA21,ALPHA22,SUMP,SUMU
     REAL RLES,RHAVEVIS
     LOGICAL DMATLUM,SUF_TOP,STREAM_UPWIND,BNDSUF
           REAL BOY_ABSXX(NGI),BOY_ABSXY(NGI),BOY_ABSXZ(NGI)
           REAL BOY_ABSYX(NGI),BOY_ABSYY(NGI),BOY_ABSYZ(NGI)
           REAL BOY_ABSZX(NGI),BOY_ABSZY(NGI),BOY_ABSZZ(NGI)
           REAL RSOXGRAGI,RSOYGRAGI,RSOZGRAGI
           REAL GRAV_DEN_NOD_NX,GRAV_DEN_NOD_NY,GRAV_DEN_NOD_NZ
           REAL GRAV_DEN_NOD,GRA_SIG_CHAN
           REAL OVER_RELAX_BOY
           REAL RPRES
! Local coords...
     INTEGER NCLOC,SMLOC
     REAL, ALLOCATABLE, DIMENSION(:,:)::NC
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCLX
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCLY
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCLZ
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCX
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCY
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCZ
     REAL, ALLOCATABLE, DIMENSION(:)::VEDUDX
     REAL, ALLOCATABLE, DIMENSION(:)::VEDUDY
     REAL, ALLOCATABLE, DIMENSION(:)::VEDUDZ
     REAL, ALLOCATABLE, DIMENSION(:)::VEDVDX
     REAL, ALLOCATABLE, DIMENSION(:)::VEDVDY
     REAL, ALLOCATABLE, DIMENSION(:)::VEDVDZ
     REAL, ALLOCATABLE, DIMENSION(:)::VEDWDX
     REAL, ALLOCATABLE, DIMENSION(:)::VEDWDY
     REAL, ALLOCATABLE, DIMENSION(:)::VEDWDZ

     REAL, ALLOCATABLE, DIMENSION(:)::HIMATX
     REAL, ALLOCATABLE, DIMENSION(:)::HIMATY
     REAL, ALLOCATABLE, DIMENSION(:)::HIMATZ
     REAL, ALLOCATABLE, DIMENSION(:)::HIMASS
 
     REAL, ALLOCATABLE, DIMENSION(:)::VECU
     REAL, ALLOCATABLE, DIMENSION(:)::VECV
     REAL, ALLOCATABLE, DIMENSION(:)::VECW

     REAL, ALLOCATABLE, DIMENSION(:)::VECYT

     REAL, ALLOCATABLE, DIMENSION(:,:,:)::EV
     INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
     INTEGER, ALLOCATABLE, DIMENSION(:)::COLELE4
         
     REAL, ALLOCATABLE, DIMENSION(:)::MASSBIGM
     REAL, ALLOCATABLE, DIMENSION(:)::MASSVEC
     
     REAL, ALLOCATABLE, DIMENSION(:,:)::SNSUB
     INTEGER, ALLOCATABLE, DIMENSION(:)::MLOCSILOC2ILOC 
     INTEGER, ALLOCATABLE, DIMENSION(:)::P_MLOCSILOC2ILOC
     REAL, ALLOCATABLE, DIMENSION(:,:)::SM
     REAL, ALLOCATABLE, DIMENSION(:,:)::ctemp
     INTEGER, ALLOCATABLE, DIMENSION(:)::BNDCON
! BNDCON marks the boundary condition nodes...
     ALLOCATE(BNDCON(NONODS))

     ewrite(1, *) "in subroutine diff3dsgs"
     ewrite(2, *) "isphere=",isphere

     BNDCON(1:NONODS)=0
     do II=1,NOBCFSimp
        GLOBI=BCT2FSimp(II)
        BNDCON(GLOBI)=1
     END DO
 
     DISOPT=DISOPT2
     IF(DISOPT2.EQ.47) DISOPT=147
     IF(DISOPT2.EQ.48) DISOPT=148
      
     NCLOC=4
     IF(ISPHERE.EQ.2) NCLOC=10

     ALLOCATE(NC(NCLOC,NGI))
     ALLOCATE(NCLX(NCLOC,NGI))
     ALLOCATE(NCLY(NCLOC,NGI))
     ALLOCATE(NCLZ(NCLOC,NGI))
     ALLOCATE(NCX(NCLOC,NGI))
     ALLOCATE(NCY(NCLOC,NGI))
     ALLOCATE(NCZ(NCLOC,NGI))
       
     IF(ISPHERE.NE.0) THEN
        ASSERT(ncloc==10)
          
        call QuadTetShapes(nc, nclx, ncly, nclz, ncloc, ngi)
         
     ELSE
        NC = N
        NCLX = NLX
        NCLY = NLY
        NCLZ = NLZ
     ENDIF

      COMAT11=0.0
      COMAT12=0.0
      COMAT13=0.0
      COMAT21=0.0
      COMAT22=0.0
      COMAT23=0.0
      COMAT31=0.0
      COMAT32=0.0
      COMAT33=0.0
      
      BOYMA11=0.0
      BOYMA12=0.0
      BOYMA13=0.0
      BOYMA21=0.0
      BOYMA22=0.0
      BOYMA23=0.0
      BOYMA31=0.0
      BOYMA32=0.0
      BOYMA33=0.0
     do I=1,NONODS
        ML(I)=0.
        VECX(I)=0.
        VECY(I)=0.
        VECZ(I)=0.
     end do
     BIGM(1:NBIGM) = 0.0
       
     PRES_RHS(1:FREDOP)=0.0

     RLUMP=1.
     IF(LUMP) RLUMP=1.
     RCONSI=1.-RLUMP
     ABLUMP=1.
     RSLUMP=1.
     IF(ABSLUM) ABLUMP=1.
     IF(SLUMP)  RSLUMP=1.
! ALPHA contains the inner element stabilisation that ensures 
! a non-singular system...

     ALPHA=0.01 


     ewrite(3,*) "InnerElement: Filter section"
     if (have_option(trim(velocity_option_path) // "/prognostic/&
          &spatial_discretisation/inner_element")) then
       ewrite(3,*) "InnerElement: Filter section - detected model"
       if (have_option(trim(velocity_option_path) // "/prognostic/&
            &spatial_discretisation/inner_element/&
            &use_filter")) then
          ewrite(3,*) "InnerElement: SGS filtering turned on"
          call get_option(trim(velocity_option_path) // "/prognostic/&
            &spatial_discretisation/inner_element/&
            &use_filter/strength", alpha, default=0.0) 
       else
          ewrite(3,*) "InnerElement: SGS filtering turned off"
          alpha=0.0
       end if
     end if
     
    ewrite(1,*) "InnerElement: Filter section, alpha[",alpha,']'


! alpha=1.0 is good or DMATLUM=.TRUE. this way it works without 
! transport
     DMATLUM=.false.

! alpha2 control works for void but others do not.
     ALPHA2=0.0

     ALPHA3=0.0
     ROCEAN=0.0
     IF(SUF_TOP_OCEAN) ROCEAN=1.0

! This is for streamline upwinding only on the global part
! and not a PG method.       

     STREAM_UPWIND=(DISOPT.EQ.144).OR.(DISOPT.EQ.145)
     ALPHA11=1.0
     ALPHA12=1.0
     ALPHA21=1.0
     ALPHA22=1.0
! DISOPT=151 is LES applied to the global system. 
! DISOPT=150 is LES applied to the local Residual system (involving pertation soln) only RECOMMENDED
! DISOPT=149 is LES applied to the local system only RECOMMENDED
! 150 and 149 produce the same results.
     IF(DISOPT.EQ.149) THEN 
        ALPHA11=0.0
        ALPHA12=0.0
     ENDIF
     IF(DISOPT.EQ.150) THEN
        ALPHA11=0.0
        ALPHA21=0.0
     ENDIF
     IF(DISOPT.EQ.151) THEN
        ALPHA21=0.0
        ALPHA22=0.0
     ENDIF

     RLES=1.0
     IF((DISOPT.EQ.152).OR.(DISOPT.EQ.153)) RLES=0.125
     RHAVEVIS=1.0
     IF((DISOPT.EQ.153).OR.(DISOPT.EQ.154)) RHAVEVIS=0.0

! Calculate surface shape functions SN, SM etc...
     CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
     CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
     &             SNLOC,SNGI, &
     &             SN,SNLX,SNLY,SNLX)

     CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
     &             SNCLOC,SNGI, &
     &             SNC,SNCLX,SNCLY,SNCLX)

     IF(NSUBTLOC.EQ.10) THEN
        call QuadTetShapes(NSUB, NSUBLX, NSUBLY, NSUBLZ, NSUBTLOC, NGI)
        SNSUBTLOC=SNCLOC
        ALLOCATE(SNSUB(SNSUBTLOC,SNGI))
        SNSUB=SNC
     ELSE
        NSUB = N
        NSUBLX = NLX
        NSUBLY = NLY
        NSUBLZ = NLZ
        SNSUBTLOC=SNLOC
        ALLOCATE(SNSUB(SNSUBTLOC,SNGI))
        SNSUB=SN
     ENDIF
     ALLOCATE(MLOCSILOC2ILOC(SNSUBTLOC))
       
     IF(MLOC.EQ.10) THEN
! Quadratic...
        SMLOC=6
        ALLOCATE(SM(SMLOC,SNGI))
        SM=SNC
     ELSE
! linear...
        SMLOC=3
        ALLOCATE(SM(SMLOC,SNGI))
        SM=SN
     ENDIF
     ALLOCATE(P_MLOCSILOC2ILOC(SMLOC))

     ALLOCATE(FINELE(TOTELE+1))
     ALLOCATE(COLELE4(TOTELE*5))
     ewrite(1,*) 'going into GETFINELE4'
     CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE4,NONODS)
      
     ROT_THIS_DGNOD=.FALSE.
     DGDM1=0.0
     DGDM2=0.0
     DGDM3=0.0
        
     DGT1X=1.0
     DGT1Y=0.0
     DGT1Z=0.0
     DGT2X=0.0
     DGT2Y=1.0
     DGT2Z=0.0 
     DGNX =0.0
     DGNY =0.0
     DGNZ =1.0

! The local cords are in anti-clockwise order
     IFACE=1
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=3
     IFACE=2
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=4
     IFACE=3
     LOCLIST(IFACE,1)=3
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=4
     IFACE=4
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=3
     LOCLIST(IFACE,3)=4
     
      IF(IUNST_FREE.EQ.1) THEN
        DGCMC=0.0
      ENDIF
      
     do  ELE=1,TOTELE! Was loop 340

! Get determinant and derivatives...
        CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
     &               N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
     &               NX,NY,NZ, NCX,NCY,NCZ,&
     &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
     &               XD,YD,ZD,&
     &               ISPHERE) 
     
        do ILOC=1,NSUBTLOC
           do GI=1,NGI
              NSUBX(ILOC,GI)= A11(GI)*NSUBLX(ILOC,GI)+A12(GI)*NSUBLY(ILOC,GI)+A13(GI)*NSUBLZ(ILOC,GI)
              NSUBY(ILOC,GI)= A21(GI)*NSUBLX(ILOC,GI)+A22(GI)*NSUBLY(ILOC,GI)+A23(GI)*NSUBLZ(ILOC,GI)
              NSUBZ(ILOC,GI)= A31(GI)*NSUBLX(ILOC,GI)+A32(GI)*NSUBLY(ILOC,GI)+A33(GI)*NSUBLZ(ILOC,GI)
           END DO
        END DO
        do ILOC=1,MLOC
           do GI=1,NGI
              MX(ILOC,GI)= A11(GI)*MLX(ILOC,GI)+A12(GI)*MLY(ILOC,GI)+A13(GI)*MLZ(ILOC,GI)
              MY(ILOC,GI)= A21(GI)*MLX(ILOC,GI)+A22(GI)*MLY(ILOC,GI)+A23(GI)*MLZ(ILOC,GI)
              MZ(ILOC,GI)= A31(GI)*MLX(ILOC,GI)+A32(GI)*MLY(ILOC,GI)+A33(GI)*MLZ(ILOC,GI)
           END DO
        END DO
! ***PUT ALL SOURCES IN TOTSOU*****
        do ILOC=1,NLOC
           IGL=SONDGL((ELE-1)*NLOC+ILOC)
           IF(ISUBSOU.EQ.0) THEN
              TOTSOUX(ILOC)=SOURCX(IGL)
              TOTSOUY(ILOC)=SOURCY(IGL)
              TOTSOUZ(ILOC)=SOURCZ(IGL)
           ELSE
              NODDG=(ELE-1)*NLOC+ILOC
              TOTSOUX(ILOC)=SOURCX(IGL)+SUBSOUX(NODDG)
              TOTSOUY(ILOC)=SOURCY(IGL)+SUBSOUY(NODDG)
              TOTSOUZ(ILOC)=SOURCZ(IGL)+SUBSOUZ(NODDG)
           ENDIF
        END DO

! ***START CALCULATING DERIVATIVES OF SUBU AT NODES**********
        IF(NSUBVLOC.NE.0) THEN
           do ILOC=1,NSUBVLOC
              do JLOC=1,NSUBVLOC
                 VLM=0.0
                 VLNXN=0.0
                 VLNYN=0.0
                 VLNZN=0.0
                 do GI=1,NGI
                    RNN=NSUB(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)
                    VLM=VLM+RNN 
                    VLNXN=VLNXN+NX(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)   
                    VLNYN=VLNYN+NY(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) 
                    VLNZN=VLNZN+NZ(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) 
                 END DO

                 MATN1XN1(ILOC,JLOC)=VLNXN
                 MATN1YN1(ILOC,JLOC)=VLNYN
                 MATN1ZN1(ILOC,JLOC)=VLNZN
                 MASSM1M1(ILOC,JLOC)=VLM
              END DO
           END DO
           CALL MATDMATINV(MASSM1M1,MASSINVM1M1,NSUBVLOC)
           do ILOC=1,NSUBVLOC
              SOUSUBNU(ILOC)=0.0
              SOUSUBNV(ILOC)=0.0
              SOUSUBNW(ILOC)=0.0
              do JLOC=1,NSUBVLOC
                 JNODSUB=(ELE-1)*NSUBVLOC+JLOC
                 SOUSUBNU(ILOC)=SOUSUBNU(ILOC)-MATN1XN1(ILOC,JLOC)*NUSUB(JNODSUB)
                 SOUSUBNV(ILOC)=SOUSUBNV(ILOC)-MATN1YN1(ILOC,JLOC)*NVSUB(JNODSUB)
                 SOUSUBNW(ILOC)=SOUSUBNW(ILOC)-MATN1ZN1(ILOC,JLOC)*NWSUB(JNODSUB)
              END DO
           END DO
           do ILOC=1,NSUBVLOC
              DERSUBNUX(ILOC)=0.0
              DERSUBNVY(ILOC)=0.0
              DERSUBNWZ(ILOC)=0.0
              do JLOC=1,NSUBVLOC
                 DERSUBNUX(ILOC)=DERSUBNUX(ILOC)+MASSINVM1M1(ILOC,JLOC)*SOUSUBNU(JLOC)
                 DERSUBNVY(ILOC)=DERSUBNVY(ILOC)+MASSINVM1M1(ILOC,JLOC)*SOUSUBNV(JLOC)
                 DERSUBNWZ(ILOC)=DERSUBNWZ(ILOC)+MASSINVM1M1(ILOC,JLOC)*SOUSUBNW(JLOC)
              END DO
           END DO
        ENDIF
! *****END CALCULATING DERIVATIVES OF SUBU AT NODES**********
         
        do  GI=1,NGI! Was loop 331

           UD(GI)=0.    
           VD(GI)=0.
           WD(GI)=0.
           TUD(GI)=0.    
           TVD(GI)=0.
           TWD(GI)=0.
           UDDX(GI)=0.    
           VDDY(GI)=0.
           WDDZ(GI)=0.
           UDDX2(GI)=0.    
           VDDY2(GI)=0.
           WDDZ2(GI)=0.
             
           UGDX(GI)=0.0
           VGDY(GI)=0.0
           WGDZ(GI)=0.0
           
           ABGI(GI)=0.0
             
           ABSDX(GI)=0.0
           ABSDY(GI)=0.0
           ABSDZ(GI)=0.0
             
           MUPTXXGI=0.0
           MUPTXYGI=0.0
           MUPTXZGI=0.0
           MUPTYYGI=0.0
           MUPTYZGI=0.0
           MUPTZZGI=0.0
             
           do ILOC=1,NSUBVLOC
              UD(GI)=UD(GI)+NSUB(ILOC,GI)*NUSUB((ELE-1)*NSUBVLOC+ILOC)
              VD(GI)=VD(GI)+NSUB(ILOC,GI)*NVSUB((ELE-1)*NSUBVLOC+ILOC)
              WD(GI)=WD(GI)+NSUB(ILOC,GI)*NWSUB((ELE-1)*NSUBVLOC+ILOC)
               
              UDDX(GI)=UDDX(GI) + NSUB(ILOC,GI)*DERSUBNUX(ILOC)
              VDDY(GI)=VDDY(GI) + NSUB(ILOC,GI)*DERSUBNVY(ILOC)
              WDDZ(GI)=WDDZ(GI) + NSUB(ILOC,GI)*DERSUBNWZ(ILOC)
                
              UDDX2(GI)=UDDX2(GI) + NSUB(ILOC,GI)*DERSUBNUX(ILOC)
              VDDY2(GI)=VDDY2(GI) + NSUB(ILOC,GI)*DERSUBNVY(ILOC)
              WDDZ2(GI)=WDDZ2(GI) + NSUB(ILOC,GI)*DERSUBNWZ(ILOC)
           END DO
           IF(COGRAX) THEN
              DIRX(GI)=0.0
              DIRY(GI)=0.0
              DIRZ(GI)=1.0
           ELSE
              DIRX(GI)=XD(GI)
              DIRY(GI)=YD(GI)
              DIRZ(GI)=ZD(GI)
              RN=SQRT(DIRX(GI)**2+DIRY(GI)**2+DIRZ(GI)**2)
              DIRX(GI)=DIRX(GI)/RN
              DIRY(GI)=DIRY(GI)/RN
              DIRZ(GI)=DIRZ(GI)/RN
           ENDIF

               RSOXGRAGI=0.0
               RSOYGRAGI=0.0
               RSOZGRAGI=0.0
               
               GRAV_DEN_NOD_NX=0.0
               GRAV_DEN_NOD_NY=0.0
               GRAV_DEN_NOD_NZ=0.0
               
           do  L=1,NLOC
              IGLV=VONDGL((ELE-1)*NLOC+L)
        
              IF(DISOPN.EQ.0) THEN
                 UD(GI)=UD(GI) + N(L,GI)*(NU(IGLV)-UG(IGLV))
                 VD(GI)=VD(GI) + N(L,GI)*(NV(IGLV)-VG(IGLV))
                 WD(GI)=WD(GI) + N(L,GI)*(NW(IGLV)-WG(IGLV))
        
                 UDDX(GI)=UDDX(GI) + NX(L,GI)*(NU(IGLV)-UG(IGLV))
                 VDDY(GI)=VDDY(GI) + NY(L,GI)*(NV(IGLV)-VG(IGLV))
                 WDDZ(GI)=WDDZ(GI) + NZ(L,GI)*(NW(IGLV)-WG(IGLV))
        
                 UDDX2(GI)=UDDX2(GI) + NX(L,GI)*NU(IGLV)
                 VDDY2(GI)=VDDY2(GI) + NY(L,GI)*NV(IGLV)
                 WDDZ2(GI)=WDDZ2(GI) + NZ(L,GI)*NW(IGLV)
        
                 UGDX(GI)=UGDX(GI) + NX(L,GI)*UG(IGLV)
                 VGDY(GI)=VGDY(GI) + NY(L,GI)*VG(IGLV)
                 WGDZ(GI)=WGDZ(GI) + NZ(L,GI)*WG(IGLV)
              ENDIF
           
              ABGI(GI)=ABGI(GI)+N(L,GI)*ABSORB(IGLV)
      
              ABSDX(GI)=ABSDX(GI)+NX(L,GI)*ABSORB(IGLV)
              ABSDY(GI)=ABSDY(GI)+NY(L,GI)*ABSORB(IGLV)
              ABSDZ(GI)=ABSDZ(GI)+NZ(L,GI)*ABSORB(IGLV)
      
              MUPTXXGI=MUPTXXGI+N(L,GI)*MUPTXX(IGLV)
              MUPTXYGI=MUPTXYGI+N(L,GI)*MUPTXY(IGLV)
              MUPTXZGI=MUPTXZGI+N(L,GI)*MUPTXZ(IGLV)
              MUPTYYGI=MUPTYYGI+N(L,GI)*MUPTYY(IGLV)
              MUPTYZGI=MUPTYZGI+N(L,GI)*MUPTYZ(IGLV)
              MUPTZZGI=MUPTZZGI+N(L,GI)*MUPTZZ(IGLV)
              IF(NBUOY.NE.0) THEN
               RSOXGRAGI=RSOXGRAGI+N(L,GI)*SOXGRA(IGLV)
               RSOYGRAGI=RSOYGRAGI+N(L,GI)*SOYGRA(IGLV)
               RSOZGRAGI=RSOZGRAGI+N(L,GI)*SOZGRA(IGLV)
               
               GRAV_DEN_NOD=SQRT(SOXGRA(IGLV)**2+SOYGRA(IGLV)**2+SOZGRA(IGLV)**2)
               
               GRAV_DEN_NOD_NX=GRAV_DEN_NOD_NX + NX(L,GI)*GRAV_DEN_NOD
               GRAV_DEN_NOD_NY=GRAV_DEN_NOD_NY + NY(L,GI)*GRAV_DEN_NOD
               GRAV_DEN_NOD_NZ=GRAV_DEN_NOD_NZ + NZ(L,GI)*GRAV_DEN_NOD
              ENDIF
           end do

              IF(NBUOY.NE.0) THEN
! Calculate the absorbtion BOY_ABSXX etc associated with vertical buoyancy...               
               GRA_SIG_CHAN= &
     &                    SIGN(1.0,DIRX(GI)*RSOXGRAGI+DIRY(GI)*RSOYGRAGI &
     &                     +DIRZ(GI)*RSOZGRAGI)   
               
               OVER_RELAX_BOY=1.0
               GRAV_DEN_NOD_NX=GRA_SIG_CHAN*GRAV_DEN_NOD_NX
               GRAV_DEN_NOD_NY=GRA_SIG_CHAN*GRAV_DEN_NOD_NY
               GRAV_DEN_NOD_NZ=GRA_SIG_CHAN*GRAV_DEN_NOD_NZ 
               RMID=DT* MAX(0.0, DIRX(GI)*GRAV_DEN_NOD_NX+DIRY(GI)*GRAV_DEN_NOD_NY &
     &                     +DIRZ(GI)*GRAV_DEN_NOD_NZ )*OVER_RELAX_BOY
              ELSE
               RMID=0.0
              ENDIF
           IF(NBUOY.NE.0) THEN
             IF(VERTICAL_BOY) THEN
! only add in a vertical bouyancy term.. which acts like absorption. 
               BOY_ABSXX(GI)=RMID*DIRX(GI)*DIRX(GI)
               BOY_ABSXY(GI)=RMID*DIRX(GI)*DIRY(GI)
               BOY_ABSXZ(GI)=RMID*DIRX(GI)*DIRZ(GI)
               
               BOY_ABSYX(GI)=RMID*DIRY(GI)*DIRX(GI)
               BOY_ABSYY(GI)=RMID*DIRY(GI)*DIRY(GI)
               BOY_ABSYZ(GI)=RMID*DIRY(GI)*DIRZ(GI)
               
               BOY_ABSZX(GI)=RMID*DIRZ(GI)*DIRX(GI)
               BOY_ABSZY(GI)=RMID*DIRZ(GI)*DIRY(GI)
               BOY_ABSZZ(GI)=RMID*DIRZ(GI)*DIRZ(GI)
             ELSE
               RMID=DT
               BOY_ABSXX(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NX
               BOY_ABSXY(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NY
               BOY_ABSXZ(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NZ
               
               BOY_ABSYX(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NX
               BOY_ABSYY(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NY
               BOY_ABSYZ(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NZ
               
               BOY_ABSZX(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NX
               BOY_ABSZY(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NY
               BOY_ABSZZ(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NZ
             ENDIF
           ELSE
               BOY_ABSXX(GI)=0.0
               BOY_ABSXY(GI)=0.0
               BOY_ABSXZ(GI)=0.0
               
               BOY_ABSYX(GI)=0.0
               BOY_ABSYY(GI)=0.0
               BOY_ABSYZ(GI)=0.0
               
               BOY_ABSZX(GI)=0.0
               BOY_ABSZY(GI)=0.0
               BOY_ABSZZ(GI)=0.0
           ENDIF
           IF(STREAM_UPWIND) THEN
              do L=1,NLOC
                 IGLV=VONDGL((ELE-1)*NLOC+L)
                 IDGNOD=(ELE-1)*NLOC+L
                 IF(DISOPT.EQ.145) THEN
                    RN=sqrt(NU(IGLV)**2+NV(IGLV)**2+NW(IGLV)**2)
                 ELSE
! DISOPT=144
                    RN=sqrt((NU(IGLV)+NUSUB(IDGNOD))**2&
                         &                 +(NV(IGLV)+NVSUB(IDGNOD))**2&
                         &                 +(NW(IGLV)+NWSUB(IDGNOD))**2)
                 ENDIF
                 TUD(GI)=TUD(GI)+NX(L,GI)*RN
                 TVD(GI)=TVD(GI)+NY(L,GI)*RN
                 TWD(GI)=TWD(GI)+NZ(L,GI)*RN
              END DO
              RN=max(SQRT(TUD(GI)**2+TVD(GI)**2+TWD(GI)**2),1.e-14)
              TUD(GI)=TUD(GI)/RN
              TVD(GI)=TVD(GI)/RN
              TWD(GI)=TWD(GI)/RN
              HXGI=ABS(A11(GI)*TUD(GI)+A12(GI)*TVD(GI)+A13(GI)*TWD(GI))
              HYGI=ABS(A21(GI)*TUD(GI)+A22(GI)*TVD(GI)+A23(GI)*TWD(GI))
              HZGI=ABS(A31(GI)*TUD(GI)+A32(GI)*TVD(GI)+A33(GI)*TWD(GI))

              HOVERQ=0.5/MAX(HXGI,HYGI,HZGI,1.0E-9)

              SGS_UPWIND_FAC(GI)=0.5*HOVERQ*sqrt(UD(GI)**2+VD(GI)**2+WD(GI)**2)
           ELSE
              SGS_UPWIND_FAC(GI)=0.0
           ENDIF
! calculate tensor for oceans:
! this has a QUADRATURE POINT varying element size and shape
           CALL SIZEGIELETENS(NX,NY,NZ,NLOC,NGI,GI,&
     &     L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ )

           IF((DISOPT.GE.147).AND.(DISOPT.LE.154)) THEN

              MXXDTW(GI)=L2GIXX
              MXYDTW(GI)=L2GIXY
              MXZDTW(GI)=L2GIXZ
              MYYDTW(GI)=L2GIYY
              MYZDTW(GI)=L2GIYZ
              MZZDTW(GI)=L2GIZZ
              CALL LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             NX,NY,NZ, &
     &             NU,NV,NW, UG,VG,WG, &
     &             MXXDTW,MXYDTW,MXZDTW,&
     &             MYYDTW,MYZDTW,MZZDTW)

              MUGIXX=MXXDTW(GI)
              MUGIXY=MXYDTW(GI)
              MUGIXZ=MXZDTW(GI)
              MUGIYY=MYYDTW(GI)
              MUGIYZ=MYZDTW(GI)
              MUGIZZ=MZZDTW(GI)
           ELSE

              MUGIXX=0.0
              MUGIXY=0.0
              MUGIXZ=0.0
              MUGIYY=0.0
              MUGIYZ=0.0
              MUGIZZ=0.0
           ENDIF
       
! 0.25 is used to scale down LES dissipation from ocean dissipation
! factor of 4 is used to reduce the length scale. 

           AMXXDTW(GI)=(RLES*MUGIXX+RHAVEVIS*MUPTXXGI)*DETWEI(GI)
           AMXYDTW(GI)=(RLES*MUGIXY+RHAVEVIS*MUPTXYGI)*DETWEI(GI)
           AMXZDTW(GI)=(RLES*MUGIXZ+RHAVEVIS*MUPTXZGI)*DETWEI(GI)

           AMYYDTW(GI)=(RLES*MUGIYY+RHAVEVIS*MUPTYYGI)*DETWEI(GI)
           AMYZDTW(GI)=(RLES*MUGIYZ+RHAVEVIS*MUPTYZGI)*DETWEI(GI)

           AMZZDTW(GI)=(RLES*MUGIZZ+RHAVEVIS*MUPTZZGI)*DETWEI(GI) 

           IF(DISOPT.EQ.148) THEN
! Put LES only in the matrix AMAT.
              MXXDTW(GI)=(MUPTXXGI)*DETWEI(GI)
              MXYDTW(GI)=(MUPTXYGI)*DETWEI(GI)
              MXZDTW(GI)=(MUPTXZGI)*DETWEI(GI)

              MYYDTW(GI)=(MUPTYYGI)*DETWEI(GI)
              MYZDTW(GI)=(MUPTYZGI)*DETWEI(GI)

              MZZDTW(GI)=(MUPTZZGI)*DETWEI(GI)
           ELSE
              MXXDTW(GI)=(RLES*MUGIXX+RHAVEVIS*MUPTXXGI)*DETWEI(GI)
              MXYDTW(GI)=(RLES*MUGIXY+RHAVEVIS*MUPTXYGI)*DETWEI(GI)
              MXZDTW(GI)=(RLES*MUGIXZ+RHAVEVIS*MUPTXZGI)*DETWEI(GI)

              MYYDTW(GI)=(RLES*MUGIYY+RHAVEVIS*MUPTYYGI)*DETWEI(GI)
              MYZDTW(GI)=(RLES*MUGIYZ+RHAVEVIS*MUPTYZGI)*DETWEI(GI)

              MZZDTW(GI)=(RLES*MUGIZZ+RHAVEVIS*MUPTZZGI)*DETWEI(GI)
           ENDIF
! for LES...
        end do ! Was loop 331
   
        MASSLUMM1M1=0.0
        MASSLUMM2M2=0.0
        COMAT=0.0
        OMAT=0.0
        LOMAT=0.0

! AMAT**********************
        do  ILOC=1,NLOC! Was loop 3501
           do  JLOC=1,NLOC! Was loop 3601

              VLM=0. 
        
              VLKTEN=0.0
      
              VLND  =0.0
              VLND2 =0.0
    
              VLMAB=0.0
              VLOMEGA=0.0
              stVLKTEN=0.0

              do  GI=1,NGI! Was loop 4581
                 OMEGAGI=2.*FUNOME(XD(GI),YD(GI),ZD(GI))

                 VLN =-(NX(ILOC,GI)*UD(GI) &
     &                 +NY(ILOC,GI)*VD(GI) &
     &                 +NZ(ILOC,GI)*WD(GI) &
     &           -N(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*N(JLOC,GI)*DETWEI(GI) 
     
                 VLN2 =N(ILOC,GI)*(NX(JLOC,GI)*UD(GI) &
     &                 +NY(JLOC,GI)*VD(GI) &
     &                 +NZ(JLOC,GI)*WD(GI) )*DETWEI(GI)
             
                 VLND=VLND+VLN
                 VLND2=VLND2+VLN2
             
                 RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
! Coriolis terms lumped...

                 COMAT11(ELE,ILOC)=COMAT11(ELE,ILOC)+RNN*(0.0 - DIRX(GI)*DIRY(GI))*OMEGAGI
                 COMAT12(ELE,ILOC)=COMAT12(ELE,ILOC)+RNN*(-1.0+ DIRX(GI)*DIRX(GI))*OMEGAGI
                 COMAT21(ELE,ILOC)=COMAT21(ELE,ILOC)+RNN*(1.0 - DIRY(GI)*DIRY(GI))*OMEGAGI
                 COMAT22(ELE,ILOC)=COMAT22(ELE,ILOC)+RNN*(0.0 + DIRY(GI)*DIRX(GI))*OMEGAGI
                 COMAT31(ELE,ILOC)=COMAT31(ELE,ILOC)+RNN*(0.0 - DIRZ(GI)*DIRY(GI))*OMEGAGI
                 COMAT32(ELE,ILOC)=COMAT32(ELE,ILOC)+RNN*(0.0 + DIRZ(GI)*DIRX(GI))*OMEGAGI
          
           IF(NBUOY.NE.0) THEN
                 COMAT11(ELE,ILOC)=COMAT11(ELE,ILOC)+RNN*BOY_ABSXX(GI)
                 COMAT12(ELE,ILOC)=COMAT12(ELE,ILOC)+RNN*BOY_ABSXY(GI)
                 COMAT13(ELE,ILOC)=COMAT13(ELE,ILOC)+RNN*BOY_ABSXZ(GI)
                 
                 COMAT21(ELE,ILOC)=COMAT21(ELE,ILOC)+RNN*BOY_ABSYX(GI)
                 COMAT22(ELE,ILOC)=COMAT22(ELE,ILOC)+RNN*BOY_ABSYY(GI)
                 COMAT23(ELE,ILOC)=COMAT23(ELE,ILOC)+RNN*BOY_ABSYZ(GI)
                 
                 COMAT31(ELE,ILOC)=COMAT31(ELE,ILOC)+RNN*BOY_ABSZX(GI)
                 COMAT32(ELE,ILOC)=COMAT32(ELE,ILOC)+RNN*BOY_ABSZY(GI)
                 COMAT33(ELE,ILOC)=COMAT33(ELE,ILOC)+RNN*BOY_ABSZZ(GI)
                 
                 BOYMA11(ELE,ILOC)=BOYMA11(ELE,ILOC)+RNN*BOY_ABSXX(GI)
                 BOYMA12(ELE,ILOC)=BOYMA12(ELE,ILOC)+RNN*BOY_ABSXY(GI)
                 BOYMA13(ELE,ILOC)=BOYMA13(ELE,ILOC)+RNN*BOY_ABSXZ(GI)
                 
                 BOYMA21(ELE,ILOC)=BOYMA21(ELE,ILOC)+RNN*BOY_ABSYX(GI)
                 BOYMA22(ELE,ILOC)=BOYMA22(ELE,ILOC)+RNN*BOY_ABSYY(GI)
                 BOYMA23(ELE,ILOC)=BOYMA23(ELE,ILOC)+RNN*BOY_ABSYZ(GI)
                 
                 BOYMA31(ELE,ILOC)=BOYMA31(ELE,ILOC)+RNN*BOY_ABSZX(GI)
                 BOYMA32(ELE,ILOC)=BOYMA32(ELE,ILOC)+RNN*BOY_ABSZY(GI)
                 BOYMA33(ELE,ILOC)=BOYMA33(ELE,ILOC)+RNN*BOY_ABSZZ(GI)
           ENDIF
           
                 VLM=VLM+RNN 
                 VLOMEGA=VLOMEGA+RNN*OMEGAGI
     
                 XNDMUJLOC=(NX(JLOC,GI)*AMXXDTW(GI)&
     &                 +NY(JLOC,GI)*AMXYDTW(GI)+NZ(JLOC,GI)*AMXZDTW(GI))
                 YNDMUJLOC=(NX(JLOC,GI)*AMXYDTW(GI)&
     &                 +NY(JLOC,GI)*AMYYDTW(GI)+NZ(JLOC,GI)*AMYZDTW(GI))
                 ZNDMUJLOC=(NX(JLOC,GI)*AMXZDTW(GI)&
     &                 +NY(JLOC,GI)*AMYZDTW(GI)+NZ(JLOC,GI)*AMZZDTW(GI))
     
                 XNDMUILOC=(NX(ILOC,GI)*AMXXDTW(GI)&
     &                 +NY(ILOC,GI)*AMXYDTW(GI)+NZ(ILOC,GI)*AMXZDTW(GI))
                 YNDMUILOC=(NX(JLOC,GI)*AMXYDTW(GI)&
     &                 +NY(ILOC,GI)*AMYYDTW(GI)+NZ(ILOC,GI)*AMYZDTW(GI))
                 ZNDMUILOC=(NX(ILOC,GI)*AMXZDTW(GI)&
     &                 +NY(ILOC,GI)*AMYZDTW(GI)+NZ(ILOC,GI)*AMZZDTW(GI))
     
                 TEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC+NZ(ILOC,GI)*ZNDMUJLOC
                 VLKTEN=VLKTEN+TEN
             
                 VLMAB=VLMAB+ABGI(GI)*RNN 

                 if(STREAM_UPWIND) then
                    XNDMUJLOC=NX(JLOC,GI)*TUD(GI)*TUD(GI)&
                         &              +NY(JLOC,GI)*TUD(GI)*TVD(GI)&
                         &              +NZ(JLOC,GI)*TUD(GI)*TWD(GI)
                    YNDMUJLOC=NX(JLOC,GI)*TVD(GI)*TUD(GI)&
                         &              +NY(JLOC,GI)*TVD(GI)*TVD(GI)&
                         &              +NZ(JLOC,GI)*TVD(GI)*TWD(GI)
                    ZNDMUJLOC=NX(JLOC,GI)*TWD(GI)*TUD(GI)&
                         &              +NY(JLOC,GI)*TWD(GI)*TWD(GI)&
                         &              +NZ(JLOC,GI)*TWD(GI)*TWD(GI)
     
                    stTEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC+NZ(ILOC,GI)*ZNDMUJLOC
                    stVLKTEN=stVLKTEN+stTEN*DETWEI(GI)*SGS_UPWIND_FAC(GI)
                 endif

              end do ! Was loop 4581
         
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
! conservative form BETA=1...
! non-conservative form BETA=0...
              AMAT(ILOC,JLOC)=BETA*VLND      +(1.-BETA)*VLND2 +VLMAB +VLKTEN*ALPHA11 + stVLKTEN

              COMAT(ILOC,JLOC)=0.0
              OMAT(ILOC,JLOC)=0.0
              lOMAT(ILOC,JLOC)=0.0
             
              MASSM1M1(ILOC,JLOC)=VLM
              MASSLUMM1M1(ILOC,ILOC)=MASSLUMM1M1(ILOC,ILOC)+VLM
         
           end do ! Was loop 3601
        end do ! Was loop 3501
 
! BMAT**********************
        do  ILOC=1,NLOC! Was loop 3502
           do  JLOC=1,NSUBTLOC! Was loop 3602

              VLM=0. 
      
              VLND  =0.0
              VLND3 =0.0
    
              VLMAB=0.0
              VLNXN=0.0
              VLNYN=0.0
              VLNZN=0.0

              do  GI=1,NGI! Was loop 4582
     
                 VLN =-(NX(ILOC,GI)*UD(GI) &
     &                 +NY(ILOC,GI)*VD(GI) &
     &                 +NZ(ILOC,GI)*WD(GI) &
     &           -N(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI) 
     
                 VLN3 =-(NX(ILOC,GI)*UD(GI) &
     &                  +NY(ILOC,GI)*VD(GI) &
     &                  +NZ(ILOC,GI)*WD(GI) &
     &           +N(ILOC,GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI)
             
                 VLND=VLND+VLN
                 VLND3=VLND3+VLN3
             
                 RNN=N(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)
          
                 VLM=VLM+RNN 
             
                 VLMAB=VLMAB+ABGI(GI)*RNN   
                 VLNXN=VLNXN+NX(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)   
                 VLNYN=VLNYN+NY(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI) 
                 VLNZN=VLNZN+NZ(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI) 

              end do ! Was loop 4582
         
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
! conservative form BETA=1...
! non-conservative form BETA=0... 
              BMAT(ILOC,JLOC)=BETA*VLND          +(1.-BETA)*VLND3 +VLMAB
             
              MATN1XN2(ILOC,JLOC)=VLNXN
              MATN1YN2(ILOC,JLOC)=VLNYN
              MATN1ZN2(ILOC,JLOC)=VLNZN
             
              MASSM1M2(ILOC,JLOC)=VLM
             
           end do ! Was loop 3602
        end do ! Was loop 3502
 
! CMAT**********************
        do  ILOC=1,NSUBTLOC! Was loop 3503
           do  JLOC=1,NLOC! Was loop 3603

              VLM=0. 
               
              VLND2 =0.0
    
              VLMAB=0.0
               
              TDIVU=0.0

              do  GI=1,NGI
     
                 VLN2 =NSUB(ILOC,GI)*(NX(JLOC,GI)*UD(GI) &
     &                 +NY(JLOC,GI)*VD(GI) &
     &                 +NZ(JLOC,GI)*WD(GI) )*DETWEI(GI)
             
                 VLND2=VLND2+VLN2
                 TDIVU=TDIVU+NSUB(ILOC,GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI))*N(JLOC,GI)*DETWEI(GI) 
             
                 RNN=NSUB(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
          
                 VLM=VLM+RNN 
             
                 VLMAB=VLMAB+ABGI(GI)*RNN   

              end do
         
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
! conservative form BETA=1...
! non-conservative form BETA=0...
              CMAT(ILOC,JLOC)=BETA*(VLND2 +TDIVU)+(1.-BETA)*VLND2 +VLMAB
             
              MASSM2M1(ILOC,JLOC)=VLM
             
           end do ! Was loop 3603
        end do ! Was loop 3503
 
! DMAT**********************
        do  ILOC=1,NSUBTLOC! Was loop 3504
           do  JLOC=1,NSUBTLOC! Was loop 3604

              VLM=0. 
      
              VLND  =0.0
              VLND3 =0.0
    
              VLMAB=0.0
              VLNXN=0.0
              VLNYN=0.0
              VLNZN=0.0
             
              VLNNX=0.0  
              VLNNY=0.0
              VLNNZ=0.0
              STABINEL=0.0
              STABOME=0.0

              do  GI=1,NGI! Was loop 4584
                 OMEGAGI=2.*FUNOME(XD(GI),YD(GI),ZD(GI))
     
                 VLN =-(NSUBX(ILOC,GI)*UD(GI) &
     &                 +NSUBY(ILOC,GI)*VD(GI) &
     &                 +NSUBZ(ILOC,GI)*WD(GI) &
     &           -NSUB(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI) 
             
                 VLN3 =-(NSUBX(ILOC,GI)*UD(GI) &
     &                  +NSUBY(ILOC,GI)*VD(GI) &
     &                  +NSUBZ(ILOC,GI)*WD(GI) &
     &           +NSUB(ILOC,GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI)
             
                 VLND=VLND+VLN
                 VLND3=VLND3+VLN3
             
                 RNN=NSUB(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)
          
                 VLM=VLM+RNN 
             
                 VLMAB=VLMAB+ABGI(GI)*RNN   
                 VLNXN=VLNXN+NSUBX(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)   
                 VLNYN=VLNYN+NSUBY(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI) 
                 VLNZN=VLNZN+NSUBZ(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)  
             
                 VLNNX=VLNNX+NSUB(ILOC,GI)*NSUBX(JLOC,GI)*DETWEI(GI)   
                 VLNNY=VLNNY+NSUB(ILOC,GI)*NSUBY(JLOC,GI)*DETWEI(GI) 
                 VLNNZ=VLNNZ+NSUB(ILOC,GI)*NSUBZ(JLOC,GI)*DETWEI(GI)  
             
                 STABINEL=STABINEL+(NSUBX(ILOC,GI)*(NSUBX(JLOC,GI)*L2GIXX&
     &                +NSUBY(JLOC,GI)*L2GIXY+NSUBZ(JLOC,GI)*L2GIXZ)&
     &                            +NSUBY(ILOC,GI)*(NSUBX(JLOC,GI)*L2GIXY&
     &                +NSUBY(JLOC,GI)*L2GIYY+NSUBZ(JLOC,GI)*L2GIYZ)&
     &                            +NSUBZ(ILOC,GI)*(NSUBX(JLOC,GI)*L2GIXZ&
     &                +NSUBY(JLOC,GI)*L2GIYZ+NSUBZ(JLOC,GI)*L2GIZZ)&
     &                              )*DETWEI(GI)

              end do ! Was loop 4584
         
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
! conservative form BETA=1...
! non-conservative form BETA=0...
              DMAT(ILOC,JLOC)=BETA*VLND          +(1.-BETA)*VLND3 +VLMAB 

              STDMAT(ILOC,JLOC)= (ALPHA/DT)*STABINEL + STABOME*100.0 +0.*vlm/dt
             
              MATN2XN2(ILOC,JLOC)=VLNXN
              MATN2YN2(ILOC,JLOC)=VLNYN
              MATN2ZN2(ILOC,JLOC)=VLNZN
             
              MATN2N2X(ILOC,JLOC)=VLNNX
              MATN2N2Y(ILOC,JLOC)=VLNNY
              MATN2N2Z(ILOC,JLOC)=VLNNZ
             
              MASSM2M2(ILOC,JLOC)=VLM
              MASSLUMM2M2(ILOC,ILOC)=MASSLUMM2M2(ILOC,ILOC)+VLM
              MASSDGELE(ELE,ILOC,JLOC)=VLM
             
           end do ! Was loop 3604
        end do ! Was loop 3504
 
! Put surface integrals into DMAT *********************
! C1TDGELELOC**********************
        do  ILOC=1,MLOC
           do  JLOC=1,NLOC

              VLMXN=0.0
              VLMYN=0.0
              VLMZN=0.0
              do GI=1,NGI
                 VLMXN=VLMXN+MX(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)   
                 VLMYN=VLMYN+MY(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) 
                 VLMZN=VLMZN+MZ(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) 
              end do
               IF(IUNST_FREE.EQ.1) THEN
                 C1TDGELELOC(ILOC,JLOC)=-FTHETA*VLMXN
                 C2TDGELELOC(ILOC,JLOC)=-FTHETA*VLMYN
                 C3TDGELELOC(ILOC,JLOC)=-FTHETA*VLMZN
                 G1TDGELELOC(ILOC,JLOC)=0.0
                 G2TDGELELOC(ILOC,JLOC)=0.0
                 G3TDGELELOC(ILOC,JLOC)=0.0
               ELSE
                 C1TDGELELOC(ILOC,JLOC)=-VLMXN
                 C2TDGELELOC(ILOC,JLOC)=-VLMYN
                 C3TDGELELOC(ILOC,JLOC)=-VLMZN
               ENDIF
             
           end do
        end do
 
! calculate SNX,SNY,SNZ***
! Surface element of domain...
        do  IFACE=1,NFACE! Was loop 3344
         
           do SILOC=1,SNLOC
              SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
           END DO
           IF(MLOC.EQ.10) THEN
! Calculate P_MLOCSILOC2ILOC
              do SILOC=1,SNLOC
                 ILOC=SILOC2ILOC(SILOC)
                 do SJLOC=1,SNLOC
                    JLOC=SILOC2ILOC(SJLOC)
                    SKLOC=SILINK2(SILOC,SJLOC)
                    KLOC=ILINK2(ILOC,JLOC)
                    P_MLOCSILOC2ILOC(SKLOC)=KLOC
                 END DO
              END DO
           ELSE
              P_MLOCSILOC2ILOC=SILOC2ILOC
           ENDIF
         
           do SILOC=1,SNLOC
              ILOC=SILOC2ILOC(SILOC)
              INOD=XONDGL((ELE-1)*NLOC+ILOC)
              do SJLOC=SILOC,SNLOC,1
                 JLOC=SILOC2ILOC(SJLOC)
                 JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                 SKLOC=SILINK2(SILOC,SJLOC)
                 XSL(SKLOC)=0.5*(X(INOD)+X(JNOD))
                 YSL(SKLOC)=0.5*(Y(INOD)+Y(JNOD))
                 ZSL(SKLOC)=0.5*(Z(INOD)+Z(JNOD))
! Assume the mide side nodes are on the sphere
                 IF(ISPHERE.GE.2) THEN
                    IF(ILOC.NE.JLOC) THEN
                       RADI=SQRT(X(INOD)**2+Y(INOD)**2+Z(INOD)**2)
                       RADJ=SQRT(X(JNOD)**2+Y(JNOD)**2+Z(JNOD)**2)
                       RADP=0.5*(RADI+RADJ)
                       RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                       XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                       YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                       ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
                    ENDIF
                 ENDIF
              END DO
           END DO
       
! Form approximate surface normal (NORMX,NORMY,NORMZ)
           CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
     &                     X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)
! Recalculate the normal...
           CALL DGSDETNXLOC2(SNCLOC,SNGI, &
     &        XSL,YSL,ZSL, &
     &        SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
     &        SNORMXN,SNORMYN,SNORMZN, &
     &        NORMX,NORMY,NORMZ)

           SUD=0.0
           SVD=0.0
           SWD=0.0
           SUDGLOB=0.0
           SVDGLOB=0.0
           SWDGLOB=0.0
           do SILOC=1,SNLOC
              ILOC=SILOC2ILOC(SILOC)
              IGL =NDGLNO((ELE-1)*NLOC+ILOC)
              ISUBNOD=(ELE-1)*NSUBVLOC+ILOC

              do SGI=1,SNGI
                 IF(NSUBVLOC.EQ.0) THEN
                    SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*(NU(IGL)-UG(IGL))
                    SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*(NV(IGL)-VG(IGL))
                    SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*(NW(IGL)-WG(IGL))
                 ELSE
                    SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*(NU(IGL)+NUSUB(ISUBNOD)-UG(IGL))
                    SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*(NV(IGL)+NVSUB(ISUBNOD)-VG(IGL))
                    SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*(NW(IGL)+NWSUB(ISUBNOD)-WG(IGL))
                 ENDIF
                 SUDGLOB(SGI)=SUDGLOB(SGI) + SN(SILOC,SGI)*(NU(IGL)-UG(IGL))
                 SVDGLOB(SGI)=SVDGLOB(SGI) + SN(SILOC,SGI)*(NV(IGL)-VG(IGL))
                 SWDGLOB(SGI)=SWDGLOB(SGI) + SN(SILOC,SGI)*(NW(IGL)-WG(IGL))
              END DO
           END DO

           do SGI=1,SNGI
              NDOTQ(SGI)=SUD(SGI)*SNORMXN(SGI)+SVD(SGI)*SNORMYN(SGI) &
     &                 +SWD(SGI)*SNORMZN(SGI)
              NDOTQGLOB(SGI)=SUDGLOB(SGI)*SNORMXN(SGI)+SVDGLOB(SGI)*SNORMYN(SGI) &
                   &                     +SWDGLOB(SGI)*SNORMZN(SGI)
           END DO

! ************************
! Perform surface integration...
! NB colele4 is ordered in terms of faces.
           ELE2=COLELE4((ELE-1)*5+IFACE)
           ONE_OR_ZERO=0.0
           BNDSUF=.false.
           SUF_TOP=.FALSE.
           IF(ELE2.EQ.0) THEN

! BNDSUF=.true. if we are on a boundary condition surface i.e. inlet-outlet
! for the free surface.
              ICOUNT=0
              do SILOC=1,SNLOC
                 ILOC=SILOC2ILOC(SILOC)
                 GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                 ICOUNT=ICOUNT+BNDCON(GLOBI)
              ENDDO
              BNDSUF=(ICOUNT.EQ.SNLOC)
              IF(SUF_TOP_OCEAN.OR.(PREOPT.EQ.2).or.BNDSUF) THEN
! Only put surface term in if at top of ocean...
                 IF(COGRAX) THEN
                    LOCPNORMX=0.0
                    LOCPNORMY=0.0
                    LOCPNORMZ=1.0
                 ELSE
                    XC=0.0
                    YC=0.0
                    ZC=0.0
                    do ILOC=1,NLOC
                       IXNOD=XONDGL((ELE-1)*NLOC+ILOC)
                       XC=XC+X(IXNOD)/NLOC
                       YC=YC+Y(IXNOD)/NLOC
                       ZC=ZC+Z(IXNOD)/NLOC
                    END DO
                    RN=SQRT(XC**2+YC**2+ZC**2)
                    LOCPNORMX=XC/RN
                    LOCPNORMY=YC/RN
                    LOCPNORMZ=ZC/RN
                 ENDIF
                 AVE_SNORMXN=0.0
                 AVE_SNORMYN=0.0
                 AVE_SNORMZN=0.0
                 do SGI=1,SNGI
                    AVE_SNORMXN=AVE_SNORMXN+SNORMXN(SGI)/SNGI
                    AVE_SNORMYN=AVE_SNORMYN+SNORMYN(SGI)/SNGI
                    AVE_SNORMZN=AVE_SNORMZN+SNORMZN(SGI)/SNGI
                 END DO
                 IF(AVE_SNORMXN*LOCPNORMX&
                      &       + AVE_SNORMYN*LOCPNORMY&
                      &       + AVE_SNORMZN*LOCPNORMZ.GT.0.8) THEN
                    ONE_OR_ZERO=1.0
                    SUF_TOP=.TRUE.
                 ENDIF
                 IF((PREOPT.EQ.2).or.BNDSUF) ONE_OR_ZERO=1.0
                 IF(ONE_OR_ZERO.NE.0.0) THEN

                    do SILOC=1,SMLOC
                       ILOC=P_MLOCSILOC2ILOC(SILOC)

                       do SJLOC=1,SNSUBTLOC
                          JLOC=SILOC2ILOC(SJLOC)
! Have a surface integral on outlet boundary... 

                          C1TCONU=0.0
                          C2TCONV=0.0
                          C3TCONW=0.0
                          do SGI=1,SNGI
                             RNN=SDETWE(SGI)*SM(SILOC,SGI)*SNSUB(SJLOC,SGI)
                             C1TCONU=C1TCONU+SNORMXN(SGI)*RNN
                             C2TCONV=C2TCONV+SNORMYN(SGI)*RNN
                             C3TCONW=C3TCONW+SNORMZN(SGI)*RNN
                          END DO
                           IF(IUNST_FREE.EQ.1) THEN
                             G1TDGELELOC(ILOC,JLOC)=G1TDGELELOC(ILOC,JLOC)+FTHETA*ONE_OR_ZERO*C1TCONU
                             G2TDGELELOC(ILOC,JLOC)=G2TDGELELOC(ILOC,JLOC)+FTHETA*ONE_OR_ZERO*C2TCONV
                             G3TDGELELOC(ILOC,JLOC)=G3TDGELELOC(ILOC,JLOC)+FTHETA*ONE_OR_ZERO*C3TCONW
                           ELSE
                             C1TDGELELOC(ILOC,JLOC)=C1TDGELELOC(ILOC,JLOC)+ONE_OR_ZERO*C1TCONU
                             C2TDGELELOC(ILOC,JLOC)=C2TDGELELOC(ILOC,JLOC)+ONE_OR_ZERO*C2TCONV
                             C3TDGELELOC(ILOC,JLOC)=C3TDGELELOC(ILOC,JLOC)+ONE_OR_ZERO*C3TCONW
                           ENDIF
                       END DO
                    END DO
                    
                   if(IUNST_FREE.eq.1) then
! Put in the surface mass matrix into DGCMC which is a compressive term.
                     IF(COGRAX) THEN
                       SDIRX=0.0
                       SDIRY=0.0
                       SDIRZ=1.0
                     ELSE
                       SDIRX=0.0
                       SDIRY=0.0
                       SDIRZ=0.0
                       DO SGI=1,SNGI
                         DO SILOC=1,SNCLOC
                           SDIRX(SGI)=SDIRX(SGI)+SNC(SILOC,SGI)*XSL(SILOC)
                           SDIRY(SGI)=SDIRY(SGI)+SNC(SILOC,SGI)*YSL(SILOC)
                           SDIRZ(SGI)=SDIRZ(SGI)+SNC(SILOC,SGI)*ZSL(SILOC)
                         END DO
                         RNN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                         SDIRX(SGI)=SDIRX(SGI)/RNN
                         SDIRY(SGI)=SDIRY(SGI)/RNN
                         SDIRZ(SGI)=SDIRZ(SGI)/RNN
                       END DO
                     ENDIF
                     do SILOC=1,SMLOC
                        ILOC=P_MLOCSILOC2ILOC(SILOC)
                        PGLOBI=PNDGLN((ELE-1)*MLOC+ILOC)
                        do SJLOC=1,SMLOC
                           JLOC=P_MLOCSILOC2ILOC(SJLOC)
                           PGLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
! Have a surface integral on outlet boundary... 
                           RNN=0.0
                           do SGI=1,SNGI
                              RNN=RNN+SDETWE(SGI)*SM(SILOC,SGI)*SM(SJLOC,SGI) &
           *(SNORMXN(SGI)*SDIRX(SGI) +SNORMYN(SGI)*SDIRY(SGI)  +SNORMZN(SGI)*SDIRZ(SGI)   )
                           END DO
                           CALL POSINMAT(POS,PGLOBI,PGLOBJ,FREDOP,FINCMC,COLCMC,NCMC)
                           DGCMC(POS)=DGCMC(POS)+RNN/(DT*DT*GRAVTY)
                        END DO
                     END DO
                   endif

                 ENDIF
              ENDIF
           ENDIF

           suf_top=suf_top.or.(PREOPT.EQ.2).or.BNDSUF   
           do  SILOC=1,SNSUBTLOC
              ILOC=SILOC2ILOC(SILOC)
              do  SJLOC=1,SNSUBTLOC
                 JLOC=SILOC2ILOC(SJLOC)
! Have a surface integral on outlet boundary...                
                 VLNDOTQ=0.0               
                 allVLNDOTQ=0.0
                 do  SGI=1,SNGI
                    RNN=SDETWE(SGI)*SNSUB(SILOC,SGI)*SNSUB(SJLOC,SGI)

                    VLNDOTQ=VLNDOTQ+MAX(NDOTQ(SGI),0.0)*RNN
                    allVLNDOTQ=allVLNDOTQ+NDOTQ(SGI)*RNN
                 END DO
                 IF(SUF_TOP) THEN
! avoid applying any sort of b.c at the free surface
                    AMAT(ILOC,JLOC)=AMAT(ILOC,JLOC)+allVLNDOTQ*BETA
                    BMAT(ILOC,JLOC)=BMAT(ILOC,JLOC)+allVLNDOTQ
                    DMAT(ILOC,JLOC)=DMAT(ILOC,JLOC)+allVLNDOTQ
                 ELSE IF(ELE2.EQ.0) THEN
! Bottom bathymetry (ZERO BC)...

                 ELSE
                    DMAT(ILOC,JLOC)=DMAT(ILOC,JLOC)+VLNDOTQ
                 ENDIF
              END DO
           END DO

! Rotate if on sides or bottom of ocean only...
           IF(SUF_TOP_OCEAN) THEN

              IF((ELE2.EQ.0).AND.(.NOT.SUF_TOP).AND.(.NOT.BNDSUF)) THEN
! Place into rotations...
                 CALL PLAC_IN_SUF_ROT_DG(ELE,NLOC,TOTELE,SNSUBTLOC,SNGI,&
     &               SDETWE,SNSUB, SNORMXN,SNORMYN,SNORMZN, SILOC2ILOC,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)
              ENDIF
           ENDIF

        end do ! Was loop 3344

! END OF Put surface integrals into DMAT *********************
          
! Calculate the diffusion terms*****
! calculate MASS2INV...
        CALL MATDMATINV(MASSM2M2,MASSINVM2M2,NSUBTLOC)
        CALL ABMATRIXMULTI(KMAT22X, MASSINVM2M2, MATN2XN2,NSUBTLOC)
        KMAT22X=-KMAT22X
        CALL ABMATRIXMULTI(KMAT22Y, MASSINVM2M2, MATN2YN2,NSUBTLOC)
        KMAT22Y=-KMAT22Y
        CALL ABMATRIXMULTI(KMAT22Z, MASSINVM2M2, MATN2ZN2,NSUBTLOC)
        KMAT22Z=-KMAT22Z
         
        MUGIXX=0.0
        MUGIXY=0.0
        MUGIXZ=0.0
        MUGIYY=0.0
        MUGIYZ=0.0
        MUGIZZ=0.0
        VOL=0.0
        do GI=1,NGI
           MUGIXX=MUGIXX+MXXDTW(GI)
           MUGIXY=MUGIXY+MXYDTW(GI)
           MUGIXZ=MUGIXZ+MXZDTW(GI)
           MUGIYY=MUGIYY+MYYDTW(GI)
           MUGIYZ=MUGIYZ+MYZDTW(GI)
           MUGIZZ=MUGIZZ+MZZDTW(GI)
           VOL=VOL+DETWEI(GI)
        END DO
        MUGIXX=MUGIXX/VOL
        MUGIXY=MUGIXY/VOL
        MUGIXZ=MUGIXZ/VOL
        MUGIYY=MUGIYY/VOL
        MUGIYZ=MUGIYZ/VOL
        MUGIZZ=MUGIZZ/VOL
         
        DIFMATX=MUGIXX*KMAT22X+MUGIXY*KMAT22Y+MUGIXZ*KMAT22Z
        DIFMATY=MUGIXY*KMAT22X+MUGIYY*KMAT22Y+MUGIYZ*KMAT22Z
        DIFMATZ=MUGIXZ*KMAT22X+MUGIYZ*KMAT22Y+MUGIZZ*KMAT22Z
         
        CALL ABMATRIXMULTI(DIFMATX2, MATN2N2X, DIFMATX,NSUBTLOC)
        CALL ABMATRIXMULTI(DIFMATY2, MATN2N2Y, DIFMATY,NSUBTLOC)
        CALL ABMATRIXMULTI(DIFMATZ2, MATN2N2Z, DIFMATZ,NSUBTLOC)
        DIFMAT=DIFMATX2+DIFMATY2+DIFMATZ2
         
        DIFMATX=MUGIXX*KMAT22X+MUGIXY*KMAT22Y+MUGIXZ*KMAT22Z
        DIFMATY=MUGIXY*KMAT22X+MUGIYY*KMAT22Y+MUGIYZ*KMAT22Z
        DIFMATZ=MUGIXZ*KMAT22X+MUGIYZ*KMAT22Y+MUGIZZ*KMAT22Z
        
        CALL ABMATRIXMUL(DIFMATX2B, MATN1XN2,NLOC,NSUBTLOC, &
     &                               DIFMATX,NSUBTLOC,NSUBTLOC)
        CALL ABMATRIXMUL(DIFMATY2B, MATN1YN2,NLOC,NSUBTLOC, &
     &                               DIFMATY,NSUBTLOC,NSUBTLOC)
        CALL ABMATRIXMUL(DIFMATZ2B, MATN1ZN2,NLOC,NSUBTLOC, &
     &                               DIFMATZ,NSUBTLOC,NSUBTLOC)
        DIFMATB=DIFMATX2B+DIFMATY2B+DIFMATZ2B
         
        BMAT=BMAT+DIFMATB*ALPHA12
        DMAT=DMAT-DIFMAT*ALPHA22
! C matrix...         
        do ILOC=1,NSUBTLOC
           do JLOC=1,NLOC
              CMAT(ILOC,JLOC)=CMAT(ILOC,JLOC)+DIFMATB(JLOC,ILOC)*ALPHA21
           END DO
        END DO
! Calculate the diffusion terms***** 

! Calculate the discretised sources...
        SOUSUBU=0.0
        SOUSUBV=0.0
        SOUSUBW=0.0
! ******add in the pre-discretised source*****
        do ILOC=1,NLOC   
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           IDGNOD=(ELE-1)*NLOC+ILOC
           VECX(GLOBI)=VECX(GLOBI) + VECX_SGSADD(IDGNOD)
           VECY(GLOBI)=VECY(GLOBI) + VECY_SGSADD(IDGNOD)
           VECZ(GLOBI)=VECZ(GLOBI) + VECZ_SGSADD(IDGNOD)
           SOUSUBU(ILOC)=SOUSUBU(ILOC)+VECX_SGSADD(IDGNOD)
           SOUSUBV(ILOC)=SOUSUBV(ILOC)+VECY_SGSADD(IDGNOD)
           SOUSUBW(ILOC)=SOUSUBW(ILOC)+VECZ_SGSADD(IDGNOD)
        END DO

! AMAT:
        do ILOC=1,NLOC   
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           do JLOC=1,NLOC
              GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
              A11MAT=MASSM1M1(ILOC,JLOC)/DT -(1.-THETA)*AMAT(ILOC,JLOC)
              O11MAT=(1.)*OMAT(ILOC,JLOC)*0.0
              RO11MAT=ROCEAN*COMAT(ILOC,JLOC)
              VECX(GLOBI)=VECX(GLOBI) &
     &                  +A11MAT*U(GLOBJ) +MASSM1M1(ILOC,JLOC)*TOTSOUX(JLOC)&
     &                  +(O11MAT*V(GLOBJ)-RO11MAT*NV(GLOBJ))
              VECY(GLOBI)=VECY(GLOBI)  &
     &                  +A11MAT*V(GLOBJ) +MASSM1M1(ILOC,JLOC)*TOTSOUY(JLOC)&
     &                  -(O11MAT*U(GLOBJ)-RO11MAT*NU(GLOBJ))
              VECZ(GLOBI)=VECZ(GLOBI)   &
     &                  +A11MAT*W(GLOBJ) +MASSM1M1(ILOC,JLOC)*TOTSOUZ(JLOC)
! Lump the absorption to put in to pressure. 
              ML(GLOBI)=ML(GLOBI) + MASSM1M1(ILOC,JLOC)
           END DO

           VECX(GLOBI)=VECX(GLOBI) +COMAT11(ELE,ILOC)*UNEW(GLOBI)+COMAT12(ELE,ILOC)*VNEW(GLOBI) &
     &                             +COMAT13(ELE,ILOC)*WNEW(GLOBI)
           VECY(GLOBI)=VECY(GLOBI) +COMAT21(ELE,ILOC)*UNEW(GLOBI)+COMAT22(ELE,ILOC)*VNEW(GLOBI) &
     &                             +COMAT23(ELE,ILOC)*WNEW(GLOBI)
           VECZ(GLOBI)=VECZ(GLOBI) +COMAT31(ELE,ILOC)*UNEW(GLOBI)+COMAT32(ELE,ILOC)*VNEW(GLOBI) &
     &                             +COMAT33(ELE,ILOC)*WNEW(GLOBI)
! BMAT:
           do JLOC=1,NSUBTLOC
              A12MAT=MASSM1M2(ILOC,JLOC)/DT -(1.-THETA)*BMAT(ILOC,JLOC)
              O12MAT=(1.)*OMAT(ILOC,JLOC)*0.0
              RO12MAT=ROCEAN*COMAT(ILOC,JLOC)
              VECX(GLOBI)=VECX(GLOBI) &
     &                  +A12MAT*USUB((ELE-1)*NSUBTLOC+JLOC)&
     &                  +(O12MAT*VSUB((ELE-1)*NSUBTLOC+JLOC)-RO12MAT*NVSUB((ELE-1)*NSUBTLOC+JLOC))
              VECY(GLOBI)=VECY(GLOBI)  &
     &                  +A12MAT*VSUB((ELE-1)*NSUBTLOC+JLOC)&
     &                  -(O12MAT*USUB((ELE-1)*NSUBTLOC+JLOC)-RO12MAT*NUSUB((ELE-1)*NSUBTLOC+JLOC))
              VECZ(GLOBI)=VECZ(GLOBI)   &
     &                  +A12MAT*WSUB((ELE-1)*NSUBTLOC+JLOC)
           END DO
           VECX(GLOBI)=VECX(GLOBI) &
     & +COMAT11(ELE,ILOC)*UNEWSUB((ELE-1)*NSUBTLOC+ILOC)+COMAT12(ELE,ILOC)*VNEWSUB((ELE-1)*NSUBTLOC+ILOC) &
     & +COMAT13(ELE,ILOC)*WNEWSUB((ELE-1)*NSUBTLOC+ILOC)

           VECY(GLOBI)=VECY(GLOBI)  &
     & +COMAT21(ELE,ILOC)*UNEWSUB((ELE-1)*NSUBTLOC+ILOC)+COMAT22(ELE,ILOC)*VNEWSUB((ELE-1)*NSUBTLOC+ILOC) &
     & +COMAT23(ELE,ILOC)*WNEWSUB((ELE-1)*NSUBTLOC+ILOC)

           VECZ(GLOBI)=VECZ(GLOBI)   &
     & +COMAT31(ELE,ILOC)*UNEWSUB((ELE-1)*NSUBTLOC+ILOC)+COMAT32(ELE,ILOC)*VNEWSUB((ELE-1)*NSUBTLOC+ILOC) &
     & +COMAT33(ELE,ILOC)*WNEWSUB((ELE-1)*NSUBTLOC+ILOC)

        END DO
         
! CMAT:
        do ILOC=1,NSUBTLOC   
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           do JLOC=1,NLOC
              GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
              A21MAT=MASSM2M1(ILOC,JLOC)/DT -(1.-THETA)*CMAT(ILOC,JLOC)
              O21MAT=(1)*OMAT(ILOC,JLOC)*0.0
              RO21MAT=ROCEAN*COMAT(ILOC,JLOC)

! Put the SGS source in: 
              SOUSUBU(ILOC)=SOUSUBU(ILOC)&
     &                 +A21MAT*U(GLOBJ) +MASSM2M1(ILOC,JLOC)*TOTSOUX(JLOC)&
     &                 +(O21MAT*V(GLOBJ)-RO21MAT*NV(GLOBJ))
              SOUSUBV(ILOC)=SOUSUBV(ILOC)&
     &                 -(O21MAT*U(GLOBJ)-RO21MAT*NU(GLOBJ))&
     &                 +A21MAT*V(GLOBJ) +MASSM2M1(ILOC,JLOC)*TOTSOUY(JLOC)
              SOUSUBW(ILOC)=SOUSUBW(ILOC)&
     &                 +A21MAT*W(GLOBJ) +MASSM2M1(ILOC,JLOC)*TOTSOUZ(JLOC)
           END DO
           SOUSUBU(ILOC)=SOUSUBU(ILOC)+COMAT11(ELE,ILOC)*UNEW(GLOBI)+COMAT12(ELE,ILOC)*VNEW(GLOBI) &
     &                                +COMAT13(ELE,ILOC)*WNEW(GLOBI)
           SOUSUBV(ILOC)=SOUSUBV(ILOC)+COMAT21(ELE,ILOC)*UNEW(GLOBI)+COMAT22(ELE,ILOC)*VNEW(GLOBI) &
     &                                +COMAT23(ELE,ILOC)*WNEW(GLOBI)
           SOUSUBW(ILOC)=SOUSUBW(ILOC)+COMAT31(ELE,ILOC)*UNEW(GLOBI)+COMAT32(ELE,ILOC)*VNEW(GLOBI) &
     &                                +COMAT33(ELE,ILOC)*WNEW(GLOBI)
         
! DMAT:
           do JLOC=1,NSUBTLOC
              IF(DMATLUM) THEN
                 A22MAT=MASSlumM2M2(ILOC,JLOC)/DT -(1.-THETA)*DMAT(ILOC,JLOC)
              ELSE
                 A22MAT=MASSM2M2(ILOC,JLOC)/DT -(1.-THETA)*DMAT(ILOC,JLOC)
              ENDIF
              O22MAT=(1.)*LOMAT(ILOC,JLOC)*0.0
              RO22MAT=ROCEAN*COMAT(ILOC,JLOC)

! Put the SGS source in: 
              SOUSUBU(ILOC)=SOUSUBU(ILOC)&
     &                 +A22MAT*USUB((ELE-1)*NSUBTLOC+JLOC)&
     &                 +(O22MAT*VSUB((ELE-1)*NSUBTLOC+JLOC)-RO22MAT*NVSUB((ELE-1)*NSUBTLOC+JLOC))
              SOUSUBV(ILOC)=SOUSUBV(ILOC)&
     &                 -(O22MAT*USUB((ELE-1)*NSUBTLOC+JLOC)-RO22MAT*NUSUB((ELE-1)*NSUBTLOC+JLOC))&
     &                 +A22MAT*VSUB((ELE-1)*NSUBTLOC+JLOC)
              SOUSUBW(ILOC)=SOUSUBW(ILOC)&
     &                 +A22MAT*WSUB((ELE-1)*NSUBTLOC+JLOC)
           END DO
           SOUSUBU(ILOC)=SOUSUBU(ILOC)&
     & +COMAT11(ELE,ILOC)*UNEWSUB((ELE-1)*NSUBTLOC+ILOC)+COMAT12(ELE,ILOC)*VNEWSUB((ELE-1)*NSUBTLOC+ILOC) &
     & +COMAT13(ELE,ILOC)*WNEWSUB((ELE-1)*NSUBTLOC+ILOC)

           SOUSUBV(ILOC)=SOUSUBV(ILOC)&
     & +COMAT21(ELE,ILOC)*UNEWSUB((ELE-1)*NSUBTLOC+ILOC)+COMAT22(ELE,ILOC)*VNEWSUB((ELE-1)*NSUBTLOC+ILOC) &
     & +COMAT23(ELE,ILOC)*WNEWSUB((ELE-1)*NSUBTLOC+ILOC)

           SOUSUBW(ILOC)=SOUSUBW(ILOC)&
     & +COMAT31(ELE,ILOC)*UNEWSUB((ELE-1)*NSUBTLOC+ILOC)+COMAT32(ELE,ILOC)*VNEWSUB((ELE-1)*NSUBTLOC+ILOC) &
     & +COMAT33(ELE,ILOC)*WNEWSUB((ELE-1)*NSUBTLOC+ILOC)

        END DO
         
! Put pressure terms in...         
        do ILOC=1,NLOC  
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC) 
           do JLOC=1,MLOC
              PGLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
                 IF(IUNST_FREE.EQ.1) THEN
                   RPRES=P(PGLOBJ)+((1.-FTHETA)/FTHETA)*POLD(PGLOBJ)
                   VECX(GLOBI)=VECX(GLOBI) +C1TDGELELOC(JLOC,ILOC)*RPRES
                   VECY(GLOBI)=VECY(GLOBI) +C2TDGELELOC(JLOC,ILOC)*RPRES
                   VECZ(GLOBI)=VECZ(GLOBI) +C3TDGELELOC(JLOC,ILOC)*RPRES
                   SOUSUBU(ILOC)=SOUSUBU(ILOC) +C1TDGELELOC(JLOC,ILOC)*RPRES
                   SOUSUBV(ILOC)=SOUSUBV(ILOC) +C2TDGELELOC(JLOC,ILOC)*RPRES
                   SOUSUBW(ILOC)=SOUSUBW(ILOC) +C3TDGELELOC(JLOC,ILOC)*RPRES
                 ELSE
                   VECX(GLOBI)=VECX(GLOBI) +C1TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                   VECY(GLOBI)=VECY(GLOBI) +C2TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                   VECZ(GLOBI)=VECZ(GLOBI) +C3TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                   SOUSUBU(ILOC)=SOUSUBU(ILOC) +C1TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                   SOUSUBV(ILOC)=SOUSUBV(ILOC) +C2TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                   SOUSUBW(ILOC)=SOUSUBW(ILOC) +C3TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                 ENDIF
                 
           END DO
        END DO

        omat=0.0
        comat=0.0
        lomat=0.0

! Amend AMAT,BMAT,CMAT,DMAT to take into account time dep.
        AMAT=(1./DT)*MASSM1M1+THETA*AMAT
        BMAT=(1./DT)*MASSM1M2+THETA*BMAT
        CMAT=(1./DT)*MASSM2M1+THETA*CMAT
        IF(DMATLUM) THEN
           DMAT=(1./DT)*MASSLUMM2M2+THETA*DMAT + STDMAT
        ELSE
           DMAT=(1./DT)*MASSM2M2+THETA*DMAT    + STDMAT
        ENDIF
        
        BIG_AMAT=0.0
        BIG_BMAT=0.0
        BIG_CMAT=0.0
        BIG_DMAT=0.0
        
        BIG_AMAT(1:NLOC,1:NLOC)                  =AMAT(1:NLOC,1:NLOC)
        BIG_AMAT(1+NLOC:2*NLOC,1+NLOC:2*NLOC)    =AMAT(1:NLOC,1:NLOC)
        BIG_AMAT(1+2*NLOC:3*NLOC,1+2*NLOC:3*NLOC)=AMAT(1:NLOC,1:NLOC)
        
        BIG_BMAT(1:NLOC,1:NLOC)                  =BMAT(1:NLOC,1:NLOC)
        BIG_BMAT(1+NLOC:2*NLOC,1+NLOC:2*NLOC)    =BMAT(1:NLOC,1:NLOC)
        BIG_BMAT(1+2*NLOC:3*NLOC,1+2*NLOC:3*NLOC)=BMAT(1:NLOC,1:NLOC)
        
        BIG_CMAT(1:NLOC,1:NLOC)                  =CMAT(1:NLOC,1:NLOC)
        BIG_CMAT(1+NLOC:2*NLOC,1+NLOC:2*NLOC)    =CMAT(1:NLOC,1:NLOC)
        BIG_CMAT(1+2*NLOC:3*NLOC,1+2*NLOC:3*NLOC)=CMAT(1:NLOC,1:NLOC)
        
        BIG_DMAT(1:NLOC,1:NLOC)                  =DMAT(1:NLOC,1:NLOC)
        BIG_DMAT(1+NLOC:2*NLOC,1+NLOC:2*NLOC)    =DMAT(1:NLOC,1:NLOC)
        BIG_DMAT(1+2*NLOC:3*NLOC,1+2*NLOC:3*NLOC)=DMAT(1:NLOC,1:NLOC)
        
! This is implicit Coriolis lumped...
        do ILOC=1,NLOC
          BIG_AMAT(ILOC,       ILOC)      =BIG_AMAT(ILOC,         ILOC)  + COMAT11(ELE,ILOC)
          BIG_AMAT(ILOC+NLOC,  ILOC+NLOC) =BIG_AMAT(ILOC+NLOC,ILOC+NLOC) + COMAT22(ELE,ILOC)
          BIG_AMAT(ILOC+2*NLOC,ILOC+2*NLOC) =BIG_AMAT(ILOC+2*NLOC,ILOC+2*NLOC) + COMAT33(ELE,ILOC)
          BIG_AMAT(ILOC,       ILOC+NLOC) =COMAT12(ELE,ILOC)
          BIG_AMAT(ILOC,       ILOC+2*NLOC) =COMAT13(ELE,ILOC)
          BIG_AMAT(ILOC+NLOC,  ILOC)      =COMAT21(ELE,ILOC)
          BIG_AMAT(ILOC+NLOC,  ILOC+2*NLOC) =COMAT23(ELE,ILOC)
          BIG_AMAT(ILOC+2*NLOC,ILOC)      =COMAT31(ELE,ILOC)
          BIG_AMAT(ILOC+2*NLOC,ILOC+NLOC) =COMAT32(ELE,ILOC)
        
          BIG_BMAT(ILOC,       ILOC)      =BIG_BMAT(ILOC,         ILOC)  + COMAT11(ELE,ILOC)
          BIG_BMAT(ILOC+NLOC,  ILOC+NLOC) =BIG_BMAT(ILOC+NLOC,ILOC+NLOC) + COMAT22(ELE,ILOC)
          BIG_BMAT(ILOC+2*NLOC,ILOC+2*NLOC) =BIG_BMAT(ILOC+2*NLOC,ILOC+2*NLOC) + COMAT33(ELE,ILOC)
          BIG_BMAT(ILOC,       ILOC+NLOC) =COMAT12(ELE,ILOC)
          BIG_BMAT(ILOC,       ILOC+2*NLOC) =COMAT13(ELE,ILOC)
          BIG_BMAT(ILOC+NLOC,  ILOC)      =COMAT21(ELE,ILOC)
          BIG_BMAT(ILOC+NLOC,  ILOC+2*NLOC) =COMAT23(ELE,ILOC)
          BIG_BMAT(ILOC+2*NLOC,ILOC)      =COMAT31(ELE,ILOC)
          BIG_BMAT(ILOC+2*NLOC,ILOC+NLOC) =COMAT32(ELE,ILOC)
        
          BIG_CMAT(ILOC,       ILOC)      =BIG_CMAT(ILOC,         ILOC)  + COMAT11(ELE,ILOC)
          BIG_CMAT(ILOC+NLOC,  ILOC+NLOC) =BIG_CMAT(ILOC+NLOC,ILOC+NLOC) + COMAT22(ELE,ILOC)
          BIG_CMAT(ILOC+2*NLOC,ILOC+2*NLOC) =BIG_CMAT(ILOC+2*NLOC,ILOC+2*NLOC) + COMAT33(ELE,ILOC)
          BIG_CMAT(ILOC,       ILOC+NLOC) =COMAT12(ELE,ILOC)
          BIG_CMAT(ILOC,       ILOC+2*NLOC) =COMAT13(ELE,ILOC)
          BIG_CMAT(ILOC+NLOC,  ILOC)      =COMAT21(ELE,ILOC)
          BIG_CMAT(ILOC+NLOC,  ILOC+2*NLOC) =COMAT23(ELE,ILOC)
          BIG_CMAT(ILOC+2*NLOC,ILOC)      =COMAT31(ELE,ILOC)
          BIG_CMAT(ILOC+2*NLOC,ILOC+NLOC) =COMAT32(ELE,ILOC)
        
          BIG_DMAT(ILOC,       ILOC)      =BIG_DMAT(ILOC,         ILOC)  + COMAT11(ELE,ILOC)
          BIG_DMAT(ILOC+NLOC,  ILOC+NLOC) =BIG_DMAT(ILOC+NLOC,ILOC+NLOC) + COMAT22(ELE,ILOC)
          BIG_DMAT(ILOC+2*NLOC,ILOC+2*NLOC) =BIG_DMAT(ILOC+2*NLOC,ILOC+2*NLOC) + COMAT33(ELE,ILOC)
          BIG_DMAT(ILOC,       ILOC+NLOC) =COMAT12(ELE,ILOC)
          BIG_DMAT(ILOC,       ILOC+2*NLOC) =COMAT13(ELE,ILOC)
          BIG_DMAT(ILOC+NLOC,  ILOC)      =COMAT21(ELE,ILOC)
          BIG_DMAT(ILOC+NLOC,  ILOC+2*NLOC) =COMAT23(ELE,ILOC)
          BIG_DMAT(ILOC+2*NLOC,ILOC)      =COMAT31(ELE,ILOC)
          BIG_DMAT(ILOC+2*NLOC,ILOC+NLOC) =COMAT32(ELE,ILOC)
        END DO
        
        CMAT_STORE(ELE,1:NLOC,1:NLOC)=CMAT(1:NLOC,1:NLOC)
        DMAT_STORE(ELE,1:NLOC,1:NLOC)=DMAT(1:NLOC,1:NLOC)

! Now make bigger matrix for coupled vels

! Apply b.c's to submatrices BIG_CMAT,BIG_DMAT:
        CALL APPLY_BCS_ROT_DG(ELE,NLOC,BIG_NLOC,TOTELE,&
     &               BIG_CMAT,BIG_DMAT,&
     &               SOUSUBU,SOUSUBV,SOUSUBW,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)

! calculate DMATINV...
        CALL MATDMATINV(BIG_DMAT,BIG_DMATINV,BIG_NLOC)
! calculate D^{-1} C
        CALL ABMATRIXMUL(BIG_DMATINVC,BIG_DMATINV,BIG_NLOC,BIG_NLOC,&
     &                                 BIG_CMAT,   BIG_NLOC,BIG_NLOC)
! BDINVC=B D^{-1} C
        CALL ABMATRIXMUL(BIG_BDINVC,BIG_BMAT,    BIG_NLOC,BIG_NLOC,&
     &                               BIG_DMATINVC,BIG_NLOC,BIG_NLOC)
! BDINVMAT=B D^{-1} 
        CALL ABMATRIXMUL(BIG_BDINVMAT,BIG_BMAT,   BIG_NLOC,BIG_NLOC,&
     &                                 BIG_DMATINV,BIG_NLOC,BIG_NLOC)

! store critical matrices to calculate new SGS soln UNEWSUB...
        do ILOC=1,NSUBTLOC
           DINVSOUSUBU_STORE(ELE,ILOC)=0.0
           DINVSOUSUBV_STORE(ELE,ILOC)=0.0
           DINVSOUSUBW_STORE(ELE,ILOC)=0.0
           do JLOC=1,NSUBTLOC
              DINVSOUSUBU_STORE(ELE,ILOC)=DINVSOUSUBU_STORE(ELE,ILOC)&
     &             +BIG_DMATINV(ILOC,JLOC)              *SOUSUBU(JLOC)&
     &             +BIG_DMATINV(ILOC,JLOC+NLOC)         *SOUSUBV(JLOC)&
     &             +BIG_DMATINV(ILOC,JLOC+2*NLOC)       *SOUSUBW(JLOC)
              DINVSOUSUBV_STORE(ELE,ILOC)=DINVSOUSUBV_STORE(ELE,ILOC)&
     &             +BIG_DMATINV(ILOC+NLOC,JLOC)         *SOUSUBU(JLOC)&
     &             +BIG_DMATINV(ILOC+NLOC,JLOC+NLOC)    *SOUSUBV(JLOC)&
     &             +BIG_DMATINV(ILOC+NLOC,JLOC+2*NLOC)  *SOUSUBW(JLOC)
              DINVSOUSUBW_STORE(ELE,ILOC)=DINVSOUSUBW_STORE(ELE,ILOC)&
     &             +BIG_DMATINV(ILOC+2*NLOC,JLOC)       *SOUSUBU(JLOC)&
     &             +BIG_DMATINV(ILOC+2*NLOC,JLOC+NLOC)  *SOUSUBV(JLOC)&
     &             +BIG_DMATINV(ILOC+2*NLOC,JLOC+2*NLOC)*SOUSUBW(JLOC)
           END DO
        END DO

! amend VECX to take into account -B D^{-1}*SOUSUBU:
        do ILOC=1,NLOC  
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC) 
           do JLOC=1,NSUBTLOC
              VECX(GLOBI)=VECX(GLOBI)&
     &             -BIG_BDINVMAT(ILOC,JLOC)              *SOUSUBU(JLOC)&
     &             -BIG_BDINVMAT(ILOC,JLOC+NLOC)         *SOUSUBV(JLOC)&
     &             -BIG_BDINVMAT(ILOC,JLOC+2*NLOC)       *SOUSUBW(JLOC)
              VECY(GLOBI)=VECY(GLOBI)&
     &             -BIG_BDINVMAT(ILOC+NLOC,JLOC)         *SOUSUBU(JLOC)&
     &             -BIG_BDINVMAT(ILOC+NLOC,JLOC+NLOC)    *SOUSUBV(JLOC)&
     &             -BIG_BDINVMAT(ILOC+NLOC,JLOC+2*NLOC)  *SOUSUBW(JLOC)
              VECZ(GLOBI)=VECZ(GLOBI)&
     &             -BIG_BDINVMAT(ILOC+2*NLOC,JLOC)       *SOUSUBU(JLOC)&
     &             -BIG_BDINVMAT(ILOC+2*NLOC,JLOC+NLOC)  *SOUSUBV(JLOC)&
     &             -BIG_BDINVMAT(ILOC+2*NLOC,JLOC+2*NLOC)*SOUSUBW(JLOC)

           END DO
        END DO
         
        do ILOC=1,NLOC
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           do JLOC=1,NLOC
              GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
! Put into matrix...
! Find count. ***************************************
              CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
              do ID=1,3
                 BIG_ILOC=(ID-1)*NLOC+ILOC
                 do JD=1,3
                    BIG_JLOC=(JD-1)*NLOC+JLOC
                    IJDISP=(ID-1)*3*NCOLM + (JD-1)*NCOLM
                    BIGM(COUNT+IJDISP)=BIGM(COUNT+IJDISP)&
     &           +BIG_AMAT(BIG_ILOC,BIG_JLOC)-BIG_BDINVC(BIG_ILOC,BIG_JLOC)

                 END DO
              END DO
           END DO
        END DO
         
! Store pressure matrix...     
      IF(IUNST_FREE.EQ.1) THEN   
        G1TDGELE(ELE,1:MLOC,1:NLOC)=G1TDGELELOC(1:MLOC,1:NLOC)
        G2TDGELE(ELE,1:MLOC,1:NLOC)=G2TDGELELOC(1:MLOC,1:NLOC)
        G3TDGELE(ELE,1:MLOC,1:NLOC)=G3TDGELELOC(1:MLOC,1:NLOC)
      ENDIF       
        C1TDGELE(ELE,1:MLOC,1:NLOC)=C1TDGELELOC(1:MLOC,1:NLOC)
        C2TDGELE(ELE,1:MLOC,1:NLOC)=C2TDGELELOC(1:MLOC,1:NLOC)
        C3TDGELE(ELE,1:MLOC,1:NLOC)=C3TDGELELOC(1:MLOC,1:NLOC)
          
     end do ! Was loop 340

     CALL PMINMX(usub,totele*nSUBTloc,'usub*******  ')
     CALL PMINMX(MUPTXX,NONODS,'MUPTXX*******  ')
     CALL PMINMX(MUPTXY,NONODS,'MUPTXY*******  ')
     CALL PMINMX(MUPTXZ,NONODS,'MUPTXZ*******  ')
     CALL PMINMX(MUPTYY,NONODS,'MUPTYY*******  ') 
     CALL PMINMX(MUPTYZ,NONODS,'MUPTYZ*******  ')
     CALL PMINMX(MUPTZZ,NONODS,'MUPTZZ*******  ')

   END SUBROUTINE DIFF3DSGS

   SUBROUTINE PLAC_IN_SUF_ROT_DG(ELE,NLOC,TOTELE,SNSUBTLOC,SNGI,&
     &               SDETWE,SNSUB, SNORMXN,SNORMYN,SNORMZN, SILOC2ILOC,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)
     REAL INFINY
     PARAMETER(INFINY=1.0E+20)
! Calculate rotations matrices and b.c's.
     INTEGER ELE,NLOC,TOTELE,SNSUBTLOC,SNGI
     REAL SDETWE(SNGI),SNSUB(SNSUBTLOC,SNGI)
     REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
     INTEGER SILOC2ILOC(SNSUBTLOC)
! place into rotations matrices:
! DGT1X,DGT1Y,DGT1Z, 
! DGT2X,DGT2Y,DGT2Z, 
! DGNX,DGNY,DGNZ,
! and b.c's               DGDM1,DGDM2,DGDM3
     LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
    
! Local variables...
     INTEGER SILOC,SJLOC,ILOC,SGI,DGNODI
     REAL SUMDMI,C1TCONU,C2TCONV,C3TCONW,RNN,RN
     REAL SLIPDIRX,SLIPDIRY,SLIPDIRZ
     REAL SUFLOCVECX(SNSUBTLOC),SUFLOCVECY(SNSUBTLOC)
     REAL SUFLOCVECZ(SNSUBTLOC)
     REAL SUFMASS(SNSUBTLOC,SNSUBTLOC)
     REAL SUFMASSINV(SNSUBTLOC,SNSUBTLOC)
     REAL SNODNORMX(SNSUBTLOC),SNODNORMY(SNSUBTLOC)
     REAL SNODNORMZ(SNSUBTLOC)

     IF(ROTATDG) THEN
! Calculate the normals to the surface at the nodes...
        do SILOC=1,SNSUBTLOC
           ILOC=SILOC2ILOC(SILOC)
! Have a surface integral on bottom boundary... 
           C1TCONU=0.0
           C2TCONV=0.0
           C3TCONW=0.0
           do SGI=1,SNGI
              RNN=SDETWE(SGI)*SNSUB(SILOC,SGI)
              C1TCONU=C1TCONU+SNORMXN(SGI)*RNN
              C2TCONV=C2TCONV+SNORMYN(SGI)*RNN
              C3TCONW=C3TCONW+SNORMZN(SGI)*RNN
           END DO
           SUFLOCVECX(SILOC)=C1TCONU
           SUFLOCVECY(SILOC)=C2TCONV
           SUFLOCVECZ(SILOC)=C3TCONW
        END DO
        SUFMASS=0.0
        do SILOC=1,SNSUBTLOC
           do SJLOC=1,SNSUBTLOC
              do SGI=1,SNGI
                 RNN=SDETWE(SGI)*SNSUB(SILOC,SGI)*SNSUB(SJLOC,SGI)
                 SUFMASS(SILOC,SJLOC)=SUFMASS(SILOC,SJLOC)+RNN
              END DO
           END DO
        END DO
        CALL MATDMATINV(SUFMASS,SUFMASSINV,SNSUBTLOC)
        do SILOC=1,SNSUBTLOC
           SNODNORMX(SILOC)=0.0
           SNODNORMY(SILOC)=0.0
           SNODNORMZ(SILOC)=0.0
           do SJLOC=1,SNSUBTLOC
              SNODNORMX(SILOC)=SNODNORMX(SILOC)&
     &                           +SUFMASSINV(SILOC,SJLOC)*SUFLOCVECX(SJLOC)
              SNODNORMY(SILOC)=SNODNORMY(SILOC)&
     &                           +SUFMASSINV(SILOC,SJLOC)*SUFLOCVECY(SJLOC)
              SNODNORMZ(SILOC)=SNODNORMZ(SILOC)&
     &                           +SUFMASSINV(SILOC,SJLOC)*SUFLOCVECZ(SJLOC)
           END DO
           RN=SQRT(SNODNORMX(SILOC)**2+SNODNORMY(SILOC)**2&
     &                +SNODNORMZ(SILOC)**2)
           SNODNORMX(SILOC)=SNODNORMX(SILOC)/RN
           SNODNORMY(SILOC)=SNODNORMY(SILOC)/RN
           SNODNORMZ(SILOC)=SNODNORMZ(SILOC)/RN
        END DO
! place into rotations...
        do SILOC=1,SNSUBTLOC
           ILOC=SILOC2ILOC(SILOC)
           DGNODI=(ELE-1)*NLOC+ILOC
           ROT_THIS_DGNOD(DGNODI)=.TRUE.
           SUMDMI=DGDM1(DGNODI)+DGDM2(DGNODI)+DGDM3(DGNODI)
           IF(SUMDMI.EQ.0.0) THEN
! 1st time we rotate this node...
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NORMX NORMY NORMZ
! Calculate T1X,T1Y,T1Z, T2X,T2Y,T2Z from NORMX,NORMY,NORMZ
              DGNX(DGNODI)=SNODNORMX(SILOC)
              DGNY(DGNODI)=SNODNORMY(SILOC)
              DGNZ(DGNODI)=SNODNORMZ(SILOC)
              CALL GETROTSIM(.TRUE., &
     &                    DGT1X(DGNODI),DGT1Y(DGNODI),DGT1Z(DGNODI), &
     &                    DGT2X(DGNODI),DGT2Y(DGNODI),DGT2Z(DGNODI), &
     &                    DGNX(DGNODI),DGNY(DGNODI),DGNZ(DGNODI))
              DGDM3(DGNODI)=INFINY
           ELSE IF(SUMDMI.LT.1.5*INFINY) THEN
! 2nd time we rotate this node...
! Form cross product of new and previous noramals. 
              CALL XPROD(SLIPDIRX,SLIPDIRY,SLIPDIRZ, &
     &                    SNODNORMX(SILOC),SNODNORMY(SILOC),SNODNORMZ(SILOC),&
     &                    DGNX(DGNODI),DGNY(DGNODI),DGNZ(DGNODI))

              DGNX(DGNODI)=SLIPDIRX
              DGNY(DGNODI)=SLIPDIRY
              DGNZ(DGNODI)=SLIPDIRZ
              CALL GETROTSIM(.TRUE., &
     &                    DGT1X(DGNODI),DGT1Y(DGNODI),DGT1Z(DGNODI), &
     &                    DGT2X(DGNODI),DGT2Y(DGNODI),DGT2Z(DGNODI), &
     &                    DGNX(DGNODI),DGNY(DGNODI),DGNZ(DGNODI))
              DGDM1(DGNODI)=INFINY
              DGDM2(DGNODI)=INFINY
              DGDM3(DGNODI)=0.0

           ELSE
! 3rd time we rotate this node (apply b.c. to all directions)...
              DGT1X(DGNODI)=1.0
              DGT1Y(DGNODI)=0.0
              DGT1Z(DGNODI)=0.0
              DGT2X(DGNODI)=0.0
              DGT2Y(DGNODI)=1.0
              DGT2Z(DGNODI)=0.0
              DGNX(DGNODI)=0.0
              DGNY(DGNODI)=0.0
              DGNZ(DGNODI)=1.0
              DGDM1(DGNODI)=INFINY
              DGDM2(DGNODI)=INFINY
              DGDM3(DGNODI)=INFINY
           ENDIF

        END DO
! ENDOF IF(ROTATDG) THEN...
     ENDIF
             
   END SUBROUTINE PLAC_IN_SUF_ROT_DG

   SUBROUTINE APPLY_BCS_ROT_DG(ELE,NLOC,BIG_NLOC,TOTELE,&
     &               BIG_CMAT,BIG_DMAT,&
     &               SOUSUBU,SOUSUBV,SOUSUBW,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)
! apply b.c's to submatrices BIG_BMAT,BIG_CMAT,BIG_DMAT:
     INTEGER ELE,NLOC,TOTELE
     INTEGER BIG_NLOC
     REAL BIG_CMAT(BIG_NLOC,BIG_NLOC)
     REAL BIG_DMAT(BIG_NLOC,BIG_NLOC)
     REAL SOUSUBU(NLOC),SOUSUBV(NLOC),SOUSUBW(NLOC)
! place into rotations matrices:
! DGT1X,DGT1Y,DGT1Z, 
! DGT2X,DGT2Y,DGT2Z, 
! DGNX,DGNY,DGNZ,
! and b.c's               DGDM1,DGDM2,DGDM3
     LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Local variables...
     INTEGER SILOC,ILOC,SGI,DGNODI,ILOC3,JLOC3,NODDG
     REAL RDIAG
     REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
     REAL RSOUSUBU(NLOC),RSOUSUBV(NLOC),RSOUSUBW(NLOC)
     REAL CMATRT(BIG_NLOC,BIG_NLOC),RCMATRT(BIG_NLOC,BIG_NLOC)
     REAL DMATRT(BIG_NLOC,BIG_NLOC),RDMATRT(BIG_NLOC,BIG_NLOC)
     LOGICAL ROT
         
     IF(.NOT.ROTATDG) RETURN
     ROT=.FALSE.
     do ILOC=1,NLOC
        IF(ROT_THIS_DGNOD((ELE-1)*NLOC+ILOC)) ROT=.TRUE.
     END DO
     IF(ROT) THEN
        RLOC=0.0
        do ILOC=1,NLOC
           NODDG=(ELE-1)*NLOC+ILOC
           IF(ROT_THIS_DGNOD(NODDG)) THEN
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
              RLOC(ILOC,ILOC)       =DGT1X(NODDG)
              RLOC(ILOC,ILOC+NLOC)  =DGT1Y(NODDG)
              RLOC(ILOC,ILOC+2*NLOC)=DGT1Z(NODDG)
              RLOC(ILOC+NLOC,ILOC)       =DGT2X(NODDG)
              RLOC(ILOC+NLOC,ILOC+NLOC)  =DGT2Y(NODDG)
              RLOC(ILOC+NLOC,ILOC+2*NLOC)=DGT2Z(NODDG)
              RLOC(ILOC+2*NLOC,ILOC)       =DGNX(NODDG)
              RLOC(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(NODDG)
              RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(NODDG)

           ELSE
              RLOC(ILOC,ILOC)              =1.0
              RLOC(ILOC+NLOC,ILOC+NLOC)    =1.0
              RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=1.0
           ENDIF
        END DO
        do ILOC3=1,BIG_NLOC
           do JLOC3=1,BIG_NLOC
              RTLOC(ILOC3,JLOC3)=RLOC(JLOC3,ILOC3)
           END DO
        END DO
! Calculate rotated b.c's...
        CALL ABMATRIXMUL(CMATRT, BIG_CMAT,BIG_NLOC,BIG_NLOC, &
     &                              RTLOC,   BIG_NLOC,BIG_NLOC)
        CALL ABMATRIXMUL(RCMATRT,RLOC,    BIG_NLOC,BIG_NLOC, &
     &                              CMATRT,  BIG_NLOC,BIG_NLOC)
        CALL ABMATRIXMUL(DMATRT, BIG_DMAT,BIG_NLOC,BIG_NLOC, &
     &                              RTLOC,   BIG_NLOC,BIG_NLOC)
        CALL ABMATRIXMUL(RDMATRT,RLOC,    BIG_NLOC,BIG_NLOC, &
     &                              DMATRT,  BIG_NLOC,BIG_NLOC)
! Define b.c's...
        do ILOC=1,NLOC
           NODDG=(ELE-1)*NLOC+ILOC
           IF(ROT_THIS_DGNOD(NODDG)) THEN
! Rotate rhs vecs and set to zero (for zero bc)...
              RSOUSUBU(ILOC)=RLOC(ILOC,ILOC)*SOUSUBU(ILOC)&
     &                       +RLOC(ILOC,ILOC+NLOC)*SOUSUBV(ILOC)&
     &                       +RLOC(ILOC,ILOC+2*NLOC)*SOUSUBW(ILOC)   
              RSOUSUBV(ILOC)=RLOC(ILOC+NLOC,ILOC)*SOUSUBU(ILOC)&
     &                       +RLOC(ILOC+NLOC,ILOC+NLOC)*SOUSUBV(ILOC)&
     &                       +RLOC(ILOC+NLOC,ILOC+2*NLOC)*SOUSUBW(ILOC)
              RSOUSUBW(ILOC)=RLOC(ILOC+2*NLOC,ILOC)*SOUSUBU(ILOC)&
     &                       +RLOC(ILOC+2*NLOC,ILOC+NLOC)*SOUSUBV(ILOC)&
     &                       +RLOC(ILOC+2*NLOC,ILOC+2*NLOC)*SOUSUBW(ILOC)
              IF(DGDM1(NODDG).NE.0.0) THEN
                 RDIAG=ABS(BIG_DMAT(ILOC,ILOC))
                 RSOUSUBU(ILOC)=0.0
                 RDMATRT(ILOC,1:BIG_NLOC)=0.0
                 RDMATRT(ILOC,ILOC)=RDIAG
                 RCMATRT(ILOC,1:BIG_NLOC)=0.0
                 RCMATRT(ILOC,ILOC)=RDIAG
              ENDIF
              IF(DGDM2(NODDG).NE.0.0) THEN
                 RDIAG=ABS(BIG_DMAT(ILOC+NLOC,ILOC+NLOC))
                 RSOUSUBV(ILOC)=0.0
                 RDMATRT(ILOC+NLOC,1:BIG_NLOC)=0.0
                 RDMATRT(ILOC+NLOC,ILOC+NLOC)=RDIAG
                 RCMATRT(ILOC+NLOC,1:BIG_NLOC)=0.0
                 RCMATRT(ILOC+NLOC,ILOC+NLOC)=RDIAG
              ENDIF
              IF(DGDM3(NODDG).NE.0.0) THEN
                 RDIAG=ABS(BIG_DMAT(ILOC+2*NLOC,ILOC+2*NLOC))
                 RSOUSUBW(ILOC)=0.0
                 RDMATRT(ILOC+2*NLOC,1:BIG_NLOC)=0.0
                 RDMATRT(ILOC+2*NLOC,ILOC+2*NLOC)=RDIAG
                 RCMATRT(ILOC+2*NLOC,1:BIG_NLOC)=0.0
                 RCMATRT(ILOC+2*NLOC,ILOC+2*NLOC)=RDIAG
              ENDIF
! Now rotate source back...
              SOUSUBU(ILOC)=RTLOC(ILOC,ILOC)*RSOUSUBU(ILOC)&
     &                      +RTLOC(ILOC,ILOC+NLOC)*RSOUSUBV(ILOC)&
     &                      +RTLOC(ILOC,ILOC+2*NLOC)*RSOUSUBW(ILOC)   
              SOUSUBV(ILOC)=RTLOC(ILOC+NLOC,ILOC)*RSOUSUBU(ILOC)&
     &                      +RTLOC(ILOC+NLOC,ILOC+NLOC)*RSOUSUBV(ILOC)&
     &                      +RTLOC(ILOC+NLOC,ILOC+2*NLOC)*RSOUSUBW(ILOC)
              SOUSUBW(ILOC)=RTLOC(ILOC+2*NLOC,ILOC)*RSOUSUBU(ILOC)&
     &                      +RTLOC(ILOC+2*NLOC,ILOC+NLOC)*RSOUSUBV(ILOC)&
     &                      +RTLOC(ILOC+2*NLOC,ILOC+2*NLOC)*RSOUSUBW(ILOC)
           ENDIF
        END DO
           
! Now rotate RCMATRT & RDMATRT back 
        CALL ABMATRIXMUL(CMATRT, RTLOC,     BIG_NLOC,BIG_NLOC, &
     &                              RCMATRT,   BIG_NLOC,BIG_NLOC)
        CALL ABMATRIXMUL(BIG_CMAT,CMATRT,   BIG_NLOC,BIG_NLOC, &
     &                               RLOC,     BIG_NLOC,BIG_NLOC)
        CALL ABMATRIXMUL(DMATRT, RTLOC,     BIG_NLOC,BIG_NLOC, &
     &                              RDMATRT,   BIG_NLOC,BIG_NLOC)
        CALL ABMATRIXMUL(BIG_DMAT,DMATRT,   BIG_NLOC,BIG_NLOC, &
     &                               RLOC,     BIG_NLOC,BIG_NLOC)
     ENDIF

   END SUBROUTINE APPLY_BCS_ROT_DG

   SUBROUTINE APPLY_STAN_BC(BIGM,CENTRM,NONODS,NCOLM,NBIGM,&
     &                     U,V,W,DT,&
     &                     VECX,VECY,VECZ,&
     &                     NOBCU,NOBCV,NOBCW,&
     &                     BCU1,BCV1,BCW1,&
     &                     BCU2,BCV2,BCW2,&
     &                     TOTELE,NLOC,NDGLNO,&
     &                     DGDM1,DGDM2,DGDM3)
! Apply standard b.c.'s
     REAL INFINY
     PARAMETER(INFINY=1.0E+20)
     INTEGER NONODS,NCOLM,NBIGM,CENTRM(NONODS)
     REAL BIGM(NBIGM)
     REAL U(NONODS),V(NONODS),W(NONODS)
     REAL DT
     REAL VECX(NONODS),VECY(NONODS),VECZ(NONODS)
     INTEGER NOBCU,NOBCV,NOBCW
     INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
     REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
     INTEGER TOTELE,NLOC,NDGLNO(TOTELE*NLOC)
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Local variables...
     INTEGER II,NOD,ELE,ILOC,NODDG,GLOBI
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKU
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKV
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKW
     ALLOCATE(MARKU(NONODS))
     ALLOCATE(MARKV(NONODS))
     ALLOCATE(MARKW(NONODS))
     MARKU=.FALSE.

     do II=1,NOBCU
        NOD=BCU2(II)
        VECX(NOD)=VECX(NOD)+INFINY*(U(NOD)+DT*BCU1(II))
        BIGM(CENTRM(NOD))=BIGM(CENTRM(NOD))+INFINY
        MARKU(NOD)=.TRUE.
     END DO
        
     MARKV=.FALSE.
     do II=1,NOBCV
        NOD=BCV2(II)
        VECY(NOD)=VECY(NOD)+INFINY*(V(NOD)+DT*BCV1(II))
        BIGM(CENTRM(NOD)+4*NCOLM)=BIGM(CENTRM(NOD)+4*NCOLM)+INFINY
        MARKV(NOD)=.TRUE.
     END DO
        
     MARKW=.FALSE.
     do II=1,NOBCW
        NOD=BCW2(II)
        VECZ(NOD)=VECZ(NOD)+INFINY*(W(NOD)+DT*BCW1(II))
        BIGM(CENTRM(NOD)+8*NCOLM)=BIGM(CENTRM(NOD)+8*NCOLM)+INFINY
        MARKW(NOD)=.TRUE.
     END DO
        
   END SUBROUTINE APPLY_STAN_BC

   SUBROUTINE SOLV_FOR_SUB_UVW(USUB,VSUB,WSUB, U,V,W,&
     &               NLOC,BIG_NLOC,TOTELE,NONODS,NDGLNO,&
     &               CMAT_STORE,DMAT_STORE,&
     &               DINVSOUSUBU_STORE,DINVSOUSUBV_STORE,&
     &               DINVSOUSUBW_STORE,&
     &               COMAT11,COMAT12,COMAT13,COMAT21, &
     &               COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)
! Solve for the sub-grid-scale vels from stored matrices...
     INTEGER NLOC,TOTELE,NONODS,NDGLNO(TOTELE*NLOC)
     INTEGER BIG_NLOC
     REAL USUB(TOTELE*NLOC),VSUB(TOTELE*NLOC),WSUB(TOTELE*NLOC)
     REAL U(NONODS),V(NONODS),W(NONODS)
     REAL CMAT_STORE(TOTELE,NLOC,NLOC),DMAT_STORE(TOTELE,NLOC,NLOC)
     REAL DINVSOUSUBU_STORE(TOTELE,NLOC),DINVSOUSUBV_STORE(TOTELE,NLOC)
     REAL DINVSOUSUBW_STORE(TOTELE,NLOC)
        REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC)
        REAL COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC)
        REAL COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
! place into rotations matrices:
! DGT1X,DGT1Y,DGT1Z, 
! DGT2X,DGT2Y,DGT2Z, 
! DGNX,DGNY,DGNZ,
! and b.c's               DGDM1,DGDM2,DGDM3
     LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Local variables...
     INTEGER ILOC,JLOC,NODDG,ELE,GLOBJ
     REAL BIG_CMAT(BIG_NLOC,BIG_NLOC)
     REAL BIG_DMAT(BIG_NLOC,BIG_NLOC)
     REAL SOUSUBU(NLOC),SOUSUBV(NLOC),SOUSUBW(NLOC)
     REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
     REAL CMATRT(BIG_NLOC,BIG_NLOC),RCMATRT(BIG_NLOC,BIG_NLOC)
     REAL DMATRT(BIG_NLOC,BIG_NLOC),RDMATRT(BIG_NLOC,BIG_NLOC)
     REAL BIG_DMATINV(BIG_NLOC,BIG_NLOC),BIG_DMATINVC(BIG_NLOC,BIG_NLOC)
     REAL LOMAT(NLOC,NLOC)
         
     SOUSUBU=0.0
     SOUSUBV=0.0
     SOUSUBW=0.0
     do ELE=1,TOTELE
        LOMAT=0.0
        
        BIG_CMAT=0.0
        BIG_CMAT(1:NLOC,         1:NLOC)         =CMAT_STORE(ELE,1:NLOC,1:NLOC)
        BIG_CMAT(1+NLOC:2*NLOC,  1+NLOC:2*NLOC)  =CMAT_STORE(ELE,1:NLOC,1:NLOC)
        BIG_CMAT(1+2*NLOC:3*NLOC,1+2*NLOC:3*NLOC)=CMAT_STORE(ELE,1:NLOC,1:NLOC)
        BIG_DMAT=0.0
        BIG_DMAT(1:NLOC,         1:NLOC)         =DMAT_STORE(ELE,1:NLOC,1:NLOC)
        BIG_DMAT(1+NLOC:2*NLOC,  1+NLOC:2*NLOC)  =DMAT_STORE(ELE,1:NLOC,1:NLOC)
        BIG_DMAT(1+2*NLOC:3*NLOC,1+2*NLOC:3*NLOC)=DMAT_STORE(ELE,1:NLOC,1:NLOC)

! This is implicit Coriolis lumped...
        do ILOC=1,NLOC
          BIG_CMAT(ILOC,       ILOC)          =BIG_CMAT(ILOC,         ILOC)      + COMAT11(ELE,ILOC)
          BIG_CMAT(ILOC+NLOC,  ILOC+NLOC)     =BIG_CMAT(ILOC+NLOC,ILOC+NLOC)     + COMAT22(ELE,ILOC)
          BIG_CMAT(ILOC+2*NLOC,  ILOC+2*NLOC) =BIG_CMAT(ILOC+2*NLOC,ILOC+2*NLOC) + COMAT33(ELE,ILOC)
          BIG_CMAT(ILOC,       ILOC+NLOC)   =COMAT12(ELE,ILOC)
          BIG_CMAT(ILOC,       ILOC+2*NLOC) =COMAT13(ELE,ILOC)
          BIG_CMAT(ILOC+NLOC,  ILOC)        =COMAT21(ELE,ILOC)
          BIG_CMAT(ILOC+NLOC,  ILOC+2*NLOC) =COMAT23(ELE,ILOC)
          BIG_CMAT(ILOC+2*NLOC,ILOC)        =COMAT31(ELE,ILOC)
          BIG_CMAT(ILOC+2*NLOC,ILOC+NLOC)   =COMAT32(ELE,ILOC)
        
          BIG_DMAT(ILOC,       ILOC)          =BIG_DMAT(ILOC,         ILOC)      + COMAT11(ELE,ILOC)
          BIG_DMAT(ILOC+NLOC,  ILOC+NLOC)     =BIG_DMAT(ILOC+NLOC,ILOC+NLOC)     + COMAT22(ELE,ILOC)
          BIG_DMAT(ILOC+2*NLOC,  ILOC+2*NLOC) =BIG_DMAT(ILOC+2*NLOC,ILOC+2*NLOC) + COMAT33(ELE,ILOC)
          BIG_DMAT(ILOC,       ILOC+NLOC)   =COMAT12(ELE,ILOC)
          BIG_DMAT(ILOC,       ILOC+2*NLOC) =COMAT13(ELE,ILOC)
          BIG_DMAT(ILOC+NLOC,  ILOC)        =COMAT21(ELE,ILOC)
          BIG_DMAT(ILOC+NLOC,  ILOC+2*NLOC) =COMAT23(ELE,ILOC)
          BIG_DMAT(ILOC+2*NLOC,ILOC)        =COMAT31(ELE,ILOC)
          BIG_DMAT(ILOC+2*NLOC,ILOC+NLOC)   =COMAT32(ELE,ILOC)
        END DO

! Apply b.c's to submatrices BIG_CMAT,BIG_DMAT:
        CALL APPLY_BCS_ROT_DG(ELE,NLOC,BIG_NLOC,TOTELE,&
     &               BIG_CMAT,BIG_DMAT,&
     &               SOUSUBU,SOUSUBV,SOUSUBW,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)

! calculate DMATINV...
        CALL MATDMATINV(BIG_DMAT,BIG_DMATINV,BIG_NLOC)
! calculate D^{-1} C
        CALL ABMATRIXMUL(BIG_DMATINVC,BIG_DMATINV,BIG_NLOC,BIG_NLOC,&
     &                                   BIG_CMAT,   BIG_NLOC,BIG_NLOC)
! solve for USUB...
        do ILOC=1,NLOC
           USUB((ELE-1)*NLOC+ILOC)=DINVSOUSUBU_STORE(ELE,ILOC)
           VSUB((ELE-1)*NLOC+ILOC)=DINVSOUSUBV_STORE(ELE,ILOC)
           WSUB((ELE-1)*NLOC+ILOC)=DINVSOUSUBW_STORE(ELE,ILOC)
        END DO
           
        do JLOC=1,NLOC
           GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
           do ILOC=1,NLOC
              USUB((ELE-1)*NLOC+ILOC)=USUB((ELE-1)*NLOC+ILOC)&
     &                -BIG_DMATINVC(ILOC,JLOC)       *U(GLOBJ)&
     &                -BIG_DMATINVC(ILOC,JLOC+NLOC)  *V(GLOBJ)&
     &                -BIG_DMATINVC(ILOC,JLOC+2*NLOC)*W(GLOBJ)
              VSUB((ELE-1)*NLOC+ILOC)=VSUB((ELE-1)*NLOC+ILOC)&
     &                -BIG_DMATINVC(ILOC+NLOC,JLOC)       *U(GLOBJ)&
     &                -BIG_DMATINVC(ILOC+NLOC,JLOC+NLOC)  *V(GLOBJ)&
     &                -BIG_DMATINVC(ILOC+NLOC,JLOC+2*NLOC)*W(GLOBJ)
              WSUB((ELE-1)*NLOC+ILOC)=WSUB((ELE-1)*NLOC+ILOC)&
     &                -BIG_DMATINVC(ILOC+2*NLOC,JLOC)       *U(GLOBJ)&
     &                -BIG_DMATINVC(ILOC+2*NLOC,JLOC+NLOC)  *V(GLOBJ)&
     &                -BIG_DMATINVC(ILOC+2*NLOC,JLOC+2*NLOC)*W(GLOBJ)
! DMATINVCSTORE is not a square matrix if NSUBTLOC.ne.NLOC...
           END DO
        END DO
     END DO

   END SUBROUTINE SOLV_FOR_SUB_UVW

   SUBROUTINE MIDSNODS_DIFF(FINDRM,NONODS,FREDOP,&
     &                CENTRM,COLM,NCOLM,&
     &                NDGLNO,PNDGLN,TOTELE,NLOC,MLOC) 
     INTEGER NONODS,NCOLM,NLOC,MLOC,TOTELE,FREDOP
     INTEGER FINDRM(NONODS+1),CENTRM(NONODS)
     INTEGER COLM(NCOLM)
     INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
! FUNCTION
     INTEGER COUNT,I,ELE,PNODS,ILOC,JLOC,JCOL,INOD,JNOD,KLOC
! WORK SPACE...
     INTEGER, ALLOCATABLE, DIMENSION(:)::DONMAT
     ALLOCATE(DONMAT(NCOLM))

     do I=1,NCOLM
        DONMAT(I)=0
     END DO
     do ELE=1,TOTELE
        do KLOC=1,MLOC
           PNDGLN((ELE-1)*MLOC+KLOC)=-1
        END DO
     END DO

     PNODS=0
     do ELE=1,TOTELE
        do ILOC=1,NLOC
           do JLOC=1,NLOC
              INOD=NDGLNO((ELE-1)*NLOC+ILOC)
              JNOD=NDGLNO((ELE-1)*NLOC+JLOC)
              KLOC=ILINK2(ILOC,JLOC)
! Find JNOD...(cONSIDER ONLY UPPER DIAGONAL PART OF DONMAT)
              do COUNT=CENTRM(INOD),FINDRM(INOD+1)-1
                 JCOL=COLM(COUNT)
                 IF(JCOL.EQ.JNOD) THEN 
                    IF(DONMAT(COUNT).EQ.0) THEN
                       PNODS=PNODS+1
                       PNDGLN((ELE-1)*MLOC+KLOC)=PNODS
                       DONMAT(COUNT)=PNODS
                    ELSE
                       PNDGLN((ELE-1)*MLOC+KLOC)=DONMAT(COUNT)
                    ENDIF
                 ENDIF
              ENDDO

           END DO
        END DO
     END DO

     IF(MLOC.EQ.11) THEN
! Add a bubble function...
        do ELE=1,TOTELE
           PNDGLN((ELE-1)*MLOC+MLOC)=PNODS + ELE
        END DO
        PNODS=PNODS+TOTELE
     ENDIF

     FREDOP=PNODS

!       Sanity checks
     do ELE=1,TOTELE
        do KLOC=1,MLOC
           IF((PNDGLN((ELE-1)*MLOC+KLOC).LE.0).OR.&
     &          (PNDGLN((ELE-1)*MLOC+KLOC).GT.FREDOP)) THEN
              ewrite(2,*) 'ELE,KLOC,PNDGLN((ELE-1)*MLOC+KLOC):',&
     &                   ELE,KLOC,PNDGLN((ELE-1)*MLOC+KLOC)
              do ILOC=1,4
                 INOD=NDGLNO((ELE-1)*NLOC+ILOC)
                 ewrite(2,*) 'ILOC,INOD:',ILOC,INOD
              END DO
                FLAbort('oh dear')
             ENDIF
          END DO
       END DO

     END SUBROUTINE MIDSNODS_DIFF

     SUBROUTINE VISPRE_DIFF(PSIPRE,PG,NONODS,PGNODS,&
     &               NLOC,MLOC,TOTELE,&
     &               NDGLNO,PGNDGLNO )
! This sub maps the quadratic pressure to linear pressure...
       INTEGER NONODS,PGNODS,NLOC,MLOC,TOTELE
       REAL PSIPRE(NONODS),PG(PGNODS)
       INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
! Local variables...
       INTEGER PNOD,NOD,ILOC,KLOC,ELE
! Function..
       IF((MLOC.LE.5).OR.(MLOC.EQ.8)) THEN
          PSIPRE(1:NONODS) = PG(1:NONODS)
       ELSE
          do ELE=1,TOTELE
             do ILOC=1,NLOC
                KLOC=ILINK2(ILOC,ILOC)
                NOD =NDGLNO((ELE-1)*NLOC+ILOC)
                PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                PSIPRE(NOD)=PG(PNOD)
             END DO
          END DO
       ENDIF

     END SUBROUTINE VISPRE_DIFF

     SUBROUTINE GESPSI_DIFF(PSIPRE,PG,NONODS,PGNODS,&
     &               NLOC,MLOC,TOTELE,&
     &               NDGLNO,PGNDGLNO )
! THIS SUB MAPS LINEAR PRESSURE TO QUADRATIC FOR INITIAL GUESS...
       INTEGER NONODS,PGNODS,NLOC,MLOC,TOTELE
       REAL PSIPRE(NONODS),PG(PGNODS)
       INTEGER NDGLNO(TOTELE*NLOC),PGNDGLNO(TOTELE*MLOC)
! Local variables...
       INTEGER PNOD,INOD,JNOD,ILOC,JLOC,KLOC,ELE
! Function..
       ewrite(2,*) 'NONODS,PGNODS,NLOC,MLOC,TOTELE:',&
     &            NONODS,PGNODS,NLOC,MLOC,TOTELE
       IF((MLOC.LE.5).OR.(MLOC.EQ.8)) THEN
          PG(1:NONODS) = PSIPRE(1:NONODS)
          IF(MLOC.EQ.5) PG(NONODS+1:NONODS+1 + TOTELE - 1) = 0.0
       ELSE
          do ELE=1,TOTELE
             do ILOC=1,NLOC
                do JLOC=1,NLOC
                   KLOC=ILINK2(ILOC,JLOC)
                   INOD =NDGLNO((ELE-1)*NLOC+ILOC)
                   JNOD =NDGLNO((ELE-1)*NLOC+JLOC)
                   PNOD=PGNDGLNO((ELE-1)*MLOC+KLOC)
                   PG(PNOD)=0.5*(PSIPRE(INOD)+PSIPRE(JNOD))
                END DO
             END DO
          END DO
          IF(MLOC.EQ.11) PG(PGNODS-TOTELE+1:PGNODS-TOTELE+1 + TOTELE - 1) = 0.0
       ENDIF

     END SUBROUTINE GESPSI_DIFF

     SUBROUTINE DGPROJ(U,V,W,UNEW,VNEW,WNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD, & 
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA, &
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB, &
     &     BIGM, ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC,BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT2,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     NNODP,&
     &     GEOBAL,SCFACTH0,&
! for the SGS...
     &     NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     USUB,VSUB,WSUB,UNEWSUB,VNEWSUB,WNEWSUB,ISUB,&
     &     NBIGM,COGRAX,NOBCFSimp,BCT2FSimp,P,POLD, &
     &     SUF_TOP_OCEAN,&
     &     NOBCU,NOBCV,NOBCW,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2,&
     &     PARA,halo_tag,&     
! Pressure matrix soln:
     &     halo_tag_p,&
     &     NNODPP, &
! form DGCMC
     &     NCMC,FINCMC,COLCMC,MIDCMC,&
     &     NDPSET,D3,NEWMES2, velocity_option_path, &
     &     FTHETA,GRAVTY,IUNST_FREE)
! This code solves and SGS momentum and cty eqns...
       INTEGER MATSTR,IMP_P_ROT,IMP_P_ROT_UDG
       LOGICAL SMALL_DGMAT_MEM,WEAK_BC
       PARAMETER(MATSTR=0,SMALL_DGMAT_MEM=.TRUE.,WEAK_BC=.TRUE.)
       INTEGER SNONOD,VNONOD,XNONOD,NCT,NCOLM,NLOC,NGI
       INTEGER MLOC,DISOPT2,DISOPT,DISOPN,DISCRA
       REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
       REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
       REAL DENPT(SNONOD)
       REAL ABSORB(SNONOD)
       REAL  DT,THETA,BETA
       REAL FTHETA,GRAVTY
       INTEGER IUNST_FREE
! If IUNST_FREE=1 (else=0) then calulcate free surface 
! height as part of pressure. 
! FTHETA is the free surface theta when IUNST_FREE=1 is used
       INTEGER TOTELE,NONODS,FREDOP,NBIGM,VERSIO,ISPHERE
       INTEGER BIG_NLOC
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
       REAL U(NONODS),V(NONODS),W(NONODS)
       REAL UNEW(NONODS),VNEW(NONODS),WNEW(NONODS)
       INTEGER ISUB,INUSUB,NSUBVLOC
       REAL NUSUB(TOTELE*NSUBVLOC),NVSUB(TOTELE*NSUBVLOC)
       REAL NWSUB(TOTELE*NSUBVLOC)
       REAL USUB(TOTELE*NSUBVLOC),VSUB(TOTELE*NSUBVLOC)
       REAL WSUB(TOTELE*NSUBVLOC)
       REAL UNEWSUB(TOTELE*NSUBVLOC),VNEWSUB(TOTELE*NSUBVLOC)
       REAL WNEWSUB(TOTELE*NSUBVLOC)
       REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
       REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
! This is what we add into the SGS model - a discretised source.      
       REAL VECX_SGSADD(TOTELE*NSUBVLOC),VECY_SGSADD(TOTELE*NSUBVLOC)
       REAL VECZ_SGSADD(TOTELE*NSUBVLOC)
       INTEGER NBUOY,NSOGRASUB
! soxgra,soygra,sozgra is the direction of gravity force when geobal.le.-10.or.(HYDROS.NE.0)
       REAL SOXGRA(SNONOD),SOYGRA(SNONOD),SOZGRA(SNONOD)
       REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
       REAL SOZGRASUB(TOTELE*NSOGRASUB)
       REAL BIGM(NBIGM)
       LOGICAL COGRAX,SUF_TOP_OCEAN
       INTEGER NOBCFSimp
       INTEGER BCT2FSimp(NOBCFSimp)
       INTEGER NOBCU,NOBCV,NOBCW
       INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
       REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
       REAL ML(NONODS)
       REAL X(XNONOD),Y(XNONOD),Z(XNONOD),D0
       REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
       INTEGER ISUBSOU
       REAL SUBSOUX(ISUBSOU*TOTELE*NLOC),SUBSOUY(ISUBSOU*TOTELE*NLOC)
       REAL SUBSOUZ(ISUBSOU*TOTELE*NLOC)
       INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
       INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
       INTEGER XONDGL(TOTELE*NLOC)
      
       INTEGER FINDRM(NONODS+1),COLM(NCOLM)
       INTEGER CENTRM(NONODS)
       INTEGER FINDCT(FREDOP+1),COLCT(NCT)
       REAL P(FREDOP),POLD(FREDOP*IUNST_FREE)

       REAL M(MLOC,NGI)
       REAL MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
       REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI)
       REAL, target:: NLY(NLOC,NGI),NLZ(NLOC,NGI)
       REAL WEIGHT(NGI)
       LOGICAL LUMP,UPWIND,OPWIND,XTPETV
       INTEGER PREOPT2
       LOGICAL ABSLUM,SLUMP
       INTEGER NNODP
       INTEGER GEOBAL
       REAL SCFACTH0
       INTEGER PARA,halo_tag
! Pressure matrix soln:
       INTEGER halo_tag_p
       INTEGER NNODPP
! form DGCMC
       INTEGER NCMC,FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)
       INTEGER NDPSET
       LOGICAL D3,NEWMES2
       character(len=*), intent(in):: velocity_option_path
! Local vartiables...
       INTEGER IDUM(1),KITS,NR2,ILENG,PREOPT
       INTEGER ELE,ILOC,NODDG,NOD,I,j,count,ILE,GLOBI,IDG,MATSTRVEL
       INTEGER NCOLMELE,NBIGMELE
       REAL RDUM(1),RMIN
       LOGICAL ROTATDG,NO_IN_FLO
! Function to decode PREPRE for free surface
       INTEGER DECOPREPRE
      
       REAL, ALLOCATABLE, DIMENSION(:)::VECX
       REAL, ALLOCATABLE, DIMENSION(:)::VECY
       REAL, ALLOCATABLE, DIMENSION(:)::VECZ
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT11
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT12
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT13
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT21
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT22
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT23
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT31
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT32
       REAL, ALLOCATABLE, DIMENSION(:,:)::COMAT33
       
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA11
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA12
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA13
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA21
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA22
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA23
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA31
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA32
       REAL, ALLOCATABLE, DIMENSION(:,:)::BOYMA33
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::MASSDGELE
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::C1TDGELE
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::C2TDGELE
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::C3TDGELE
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::G1TDGELE
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::G2TDGELE
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::G3TDGELE
       REAL, ALLOCATABLE, DIMENSION(:)::UDG
       REAL, ALLOCATABLE, DIMENSION(:)::VDG
       REAL, ALLOCATABLE, DIMENSION(:)::WDG
       REAL, ALLOCATABLE, DIMENSION(:)::DGCMC
       REAL, ALLOCATABLE, DIMENSION(:)::DP
       REAL, ALLOCATABLE, DIMENSION(:)::RHSVEC
       REAL, ALLOCATABLE, DIMENSION(:)::DGUVW
       REAL, ALLOCATABLE, DIMENSION(:)::VECXYZ
       LOGICAL, ALLOCATABLE, DIMENSION(:)::ROT_THIS_DGNOD
       REAL, ALLOCATABLE, DIMENSION(:)::DGNX,DGNY,DGNZ, &
     &                                  DGT1X,DGT1Y,DGT1Z, &
     &                                  DGT2X,DGT2Y,DGT2Z, &
     &                                  DGDM1,DGDM2,DGDM3
       REAL, ALLOCATABLE, DIMENSION(:)::PRES_RHS
       REAL, ALLOCATABLE, DIMENSION(:)::COR_DGU
       REAL, ALLOCATABLE, DIMENSION(:)::COR_DGV
       REAL, ALLOCATABLE, DIMENSION(:)::COR_DGW
       INTEGER, ALLOCATABLE, DIMENSION(:)::FINDRMELE
       INTEGER, ALLOCATABLE, DIMENSION(:)::COLMELE
       INTEGER, ALLOCATABLE, DIMENSION(:)::CENTRMELE
       REAL, ALLOCATABLE, DIMENSION(:)::DGU,DGV,DGW,DGUNEW,DGVNEW,DGWNEW,DGNU,DGNV,DGNW

       REAL, ALLOCATABLE, DIMENSION(:)::UG2,VG2,WG2,DGNU2,DGNV2,DGNW2
       REAL, ALLOCATABLE, DIMENSION(:)::BIGMELE,MLELE,DGZER,DGONE
       REAL, ALLOCATABLE, DIMENSION(:)::DGU_THETA,DGV_THETA,DGW_THETA
       REAL, ALLOCATABLE, DIMENSION(:)::DGUOLD,DGVOLD,DGWOLD
       ewrite(1,*) 'in subroutine dgproj' 
       
       ILENG=3*TOTELE*NLOC
       ALLOCATE(VECX(TOTELE*NLOC))
       ALLOCATE(VECY(TOTELE*NLOC))
       ALLOCATE(VECZ(TOTELE*NLOC))
       ALLOCATE(ROT_THIS_DGNOD(TOTELE*NLOC))
       ALLOCATE(DGNX(TOTELE*NLOC))
       ALLOCATE(DGNY(TOTELE*NLOC))
       ALLOCATE(DGNZ(TOTELE*NLOC))
       ALLOCATE(DGT1X(TOTELE*NLOC))
       ALLOCATE(DGT1Y(TOTELE*NLOC))
       ALLOCATE(DGT1Z(TOTELE*NLOC))
       ALLOCATE(DGT2X(TOTELE*NLOC))
       ALLOCATE(DGT2Y(TOTELE*NLOC))
       ALLOCATE(DGT2Z(TOTELE*NLOC))
       ALLOCATE(DGDM1(TOTELE*NLOC))
       ALLOCATE(DGDM2(TOTELE*NLOC))
       ALLOCATE(DGDM3(TOTELE*NLOC))

       ROTATDG=.FALSE.
! IF PREOPT2.GE.10 then put Coriolis into pressure calc. 
! IF PREOPT2.GE.100 then also put Coriolis in distribution of velocity 
! between SGS and global solns.  
! IF PREOPT=2 then treat every boundary like a free surface - need 
! to do this for open boundaries. 
       IMP_P_ROT=0
       IMP_P_ROT_UDG=0
       IF(MOD(PREOPT2,100).GE.10)   IMP_P_ROT=1
       PREOPT=MOD(MOD(PREOPT2,10),100)

       ALLOCATE(COMAT11(TOTELE,NLOC))
       ALLOCATE(COMAT12(TOTELE,NLOC))
       ALLOCATE(COMAT13(TOTELE,NLOC))
       ALLOCATE(COMAT21(TOTELE,NLOC))
       ALLOCATE(COMAT22(TOTELE,NLOC))
       ALLOCATE(COMAT23(TOTELE,NLOC))
       ALLOCATE(COMAT31(TOTELE,NLOC))
       ALLOCATE(COMAT32(TOTELE,NLOC))
       ALLOCATE(COMAT33(TOTELE,NLOC))
       
       ALLOCATE(BOYMA11(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA12(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA13(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA21(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA22(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA23(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA31(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA32(NBUOY*TOTELE,NLOC))
       ALLOCATE(BOYMA33(NBUOY*TOTELE,NLOC))
       
       ALLOCATE(MASSDGELE(TOTELE,NLOC,NLOC))
       ALLOCATE(C1TDGELE(TOTELE,MLOC,NLOC))
       ALLOCATE(C2TDGELE(TOTELE,MLOC,NLOC))
       ALLOCATE(C3TDGELE(TOTELE,MLOC,NLOC))
       ALLOCATE(G1TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC))
       ALLOCATE(G2TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC))
       ALLOCATE(G3TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC))
       ALLOCATE(PRES_RHS(FREDOP))
! Form matrix structure...
       ALLOCATE(FINDRMELE(TOTELE+1))
       ALLOCATE(COLMELE(TOTELE*5))
       ALLOCATE(CENTRMELE(TOTELE))
       CALL SIMP_DG_MAT_STRU(FINDRMELE,COLMELE,NCOLMELE,TOTELE,NLOC,&
     &                       NDGLNO,CENTRMELE,NONODS)
! Map variables to DG variles...
       ALLOCATE(DGU(TOTELE*NLOC))
       ALLOCATE(DGV(TOTELE*NLOC))
       ALLOCATE(DGW(TOTELE*NLOC))
       ALLOCATE(DGUNEW(TOTELE*NLOC))
       ALLOCATE(DGVNEW(TOTELE*NLOC))
       ALLOCATE(DGWNEW(TOTELE*NLOC))
       ALLOCATE(DGNU(TOTELE*NLOC))
       ALLOCATE(DGNV(TOTELE*NLOC))
       ALLOCATE(DGNW(TOTELE*NLOC))
       ALLOCATE(MLELE(TOTELE*NLOC))
       do ELE=1,TOTELE
          do ILOC=1,NLOC
             GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
             IDG=(ELE-1)*NLOC+ILOC
             DGU(IDG)=U(GLOBI)+USUB(IDG)
             DGV(IDG)=V(GLOBI)+VSUB(IDG)
             DGW(IDG)=W(GLOBI)+WSUB(IDG)
             DGUNEW(IDG)=UNEW(GLOBI)+UNEWSUB(IDG)
             DGVNEW(IDG)=VNEW(GLOBI)+VNEWSUB(IDG)
             DGWNEW(IDG)=WNEW(GLOBI)+WNEWSUB(IDG)
             DGNU(IDG)=NU(GLOBI)+NUSUB(IDG)
             DGNV(IDG)=NV(GLOBI)+NVSUB(IDG)
             DGNW(IDG)=NW(GLOBI)+NWSUB(IDG)
          END DO
       END DO
       ALLOCATE(DGUOLD(TOTELE*NLOC))
       ALLOCATE(DGVOLD(TOTELE*NLOC))
       ALLOCATE(DGWOLD(TOTELE*NLOC))
       
       DGUOLD=DGU
       DGVOLD=DGV
       DGWOLD=DGW
       ALLOCATE(DGNU2(TOTELE*NLOC))
       ALLOCATE(DGNV2(TOTELE*NLOC))
       ALLOCATE(DGNW2(TOTELE*NLOC))
       ALLOCATE(UG2(NONODS))
       ALLOCATE(VG2(NONODS))
       ALLOCATE(WG2(NONODS))
       NO_IN_FLO=.true.
       IF(NO_IN_FLO) THEN
          UG2=0.0
          VG2=0.0
          WG2=0.0
          do ELE=1,TOTELE
             do ILOC=1,NLOC
                GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                IDG=(ELE-1)*NLOC+ILOC
                DGNU2(IDG)=DGNU(IDG)-UG(GLOBI)
                DGNV2(IDG)=DGNV(IDG)-VG(GLOBI)
                DGNW2(IDG)=DGNW(IDG)-WG(GLOBI)
             END DO
          END DO
          CALL MAKE_NO_INFLOW(TOTELE,NLOC,XONDGL,3,XNONOD,&
     &                       DGNU2,DGNV2,DGNW2,X,Y,Z)
       ELSE
          UG2=UG
          VG2=VG
          WG2=WG
          DGNU2=DGNU
          DGNV2=DGNV
          DGNW2=DGNW
       ENDIF
       IF(SMALL_DGMAT_MEM) THEN
          NBIGMELE=NLOC*NLOC * NCOLMELE + 21*TOTELE*NLOC
       ELSE
          NBIGMELE=NLOC*3*NLOC*3 * NCOLMELE
       ENDIF
       ALLOCATE(BIGMELE(NBIGMELE))
       
       ALLOCATE(DGCMC(NCMC*IUNST_FREE))

! form momentum eqns and put pressure into rhs,
! also put b.c's into matrix...
       CALL DIFF3DSGS_DG(DGU,DGV,DGW,DGUNEW,DGVNEW,DGWNEW,&
     &     DGNU2,DGNV2,DGNW2,UG2,VG2,WG2,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD,&
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA, &
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB, &
     &     COMAT11,COMAT12,COMAT13,COMAT21,COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &     BOYMA11,BOYMA12,BOYMA13,BOYMA21,BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33, &
     &     VECX,VECY,VECZ, &
     &     BIGMELE,SMALL_DGMAT_MEM, MLELE,ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRMELE,COLMELE,NCOLMELE,NONODS,CENTRMELE, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT,&
     &     ABSLUM,SLUMP, &
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     nnodp,para,halo_tag,&
     &     GEOBAL,SCFACTH0,&
! for the DG...
     &     NBIGMELE,COGRAX,NOBCFSimp,BCT2FSimp, &
     &     MASSDGELE,C1TDGELE,C2TDGELE,C3TDGELE,P,&
     &     FTHETA,GRAVTY,IUNST_FREE,G1TDGELE,G2TDGELE,G3TDGELE,POLD,&
     &     DGCMC,NCMC,FINCMC,COLCMC,MIDCMC, &
     &     SUF_TOP_OCEAN,&
     &     ROTATDG,ROT_THIS_DGNOD,&
     &     DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &     DGDM1,DGDM2,DGDM3,&
     &     PRES_RHS,WEAK_BC,&
     &     NOBCU,NOBCV,NOBCW,&
     &     BCU1,BCV1,BCW1,&
     &     BCU2,BCV2,BCW2, NDPSET, velocity_option_path)

! Apply DG to R VEC and form this R VEC

       IF(.NOT.WEAK_BC) THEN
          ewrite(1,*) 'going into DGDM123_DG_ROT_BC'
          CALL DGDM123_DG_ROT_BC(&
     &                     NLOC,BIG_NLOC,TOTELE,NDGLNO,NONODS,&
     &                     VECX,VECY,VECZ,&
     &                     DGU,DGV,DGW,DT,&
     &                     NOBCU,NOBCV,NOBCW,&
     &                     BCU1,BCV1,BCW1,&
     &                     BCU2,BCV2,BCW2, &
     &                     ROTATDG,ROT_THIS_DGNOD,&
     &                     DGT1X,DGT1Y,DGT1Z, &
     &                     DGT2X,DGT2Y,DGT2Z, &
     &                     DGNX,DGNY,DGNZ,&
     &                     DGDM1,DGDM2,DGDM3)
          ewrite(1,*) 'finished DGDM123_DG_ROT_BC SMALL_DGMAT_MEM=',SMALL_DGMAT_MEM
       ENDIF

       IF(SMALL_DGMAT_MEM) THEN
! Apply key variables to the end of BIGMELE which are retrieved in 
! the matrix solver...
          CALL PACK_BIGMELE_VECS(TOTELE,NLOC,NBIGMELE,BIGMELE,&
     &                     COMAT11,COMAT12,COMAT13,COMAT21, &
     &                     COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &                     DGT1X,DGT1Y,DGT1Z, &
     &                     DGT2X,DGT2Y,DGT2Z, &
     &                     DGNX,DGNY,DGNZ,&
     &                     DGDM1,DGDM2,DGDM3)
       ELSE
! Apply DG b.c.'s also rotate R A R^T to matrix 
          CALL APPLY_DG_ROT_BC_MAT(BIGMELE,CENTRMELE,FINDRMELE,&
     &                     COLMELE,NCOLMELE,NBIGMELE,&
     &                     NLOC,BIG_NLOC,TOTELE,&
     &                     ROTATDG,ROT_THIS_DGNOD,&
     &                     DGT1X,DGT1Y,DGT1Z, &
     &                     DGT2X,DGT2Y,DGT2Z, &
     &                     DGNX,DGNY,DGNZ,&
     &                     DGDM1,DGDM2,DGDM3)
       ENDIF
     
       CALL PMINMX(DGU,       TOTELE*NLOC,'******DGU  ')
       CALL PMINMX(DGV,       TOTELE*NLOC,'******DGV  ')
       CALL PMINMX(DGW,       TOTELE*NLOC,'******DGW  ')
! Rotate U=R U...
       CALL UVW_ROT_DG(NLOC,BIG_NLOC,TOTELE,&
     &               .FALSE.,DGU,DGV,DGW,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ)
! Solve momentum eqns... 
       ILENG=3*TOTELE*NLOC
       ILE=TOTELE
       NR2=3*ILENG
       ALLOCATE(DGUVW(ILENG))
       DGUVW=0.0
       ALLOCATE(VECXYZ(ILENG))
       do ELE=1,TOTELE
          do ILOC=1,NLOC
             DGUVW((ELE-1)*3*NLOC+ILOC)       =DGUNEW((ELE-1)*NLOC+ILOC)
             DGUVW((ELE-1)*3*NLOC+ILOC+NLOC)  =DGVNEW((ELE-1)*NLOC+ILOC)
             DGUVW((ELE-1)*3*NLOC+ILOC+2*NLOC)=DGWNEW((ELE-1)*NLOC+ILOC)
             VECXYZ((ELE-1)*3*NLOC+ILOC)       =VECX((ELE-1)*NLOC+ILOC)
             VECXYZ((ELE-1)*3*NLOC+ILOC+NLOC)  =VECY((ELE-1)*NLOC+ILOC)
             VECXYZ((ELE-1)*3*NLOC+ILOC+2*NLOC)=VECZ((ELE-1)*NLOC+ILOC)
          END DO
       END DO

       CALL PMINMX(DGU,TOTELE*NLOC,'******DGU  ')
       CALL PMINMX(DGV,TOTELE*NLOC,'******DGV  ')
       CALL PMINMX(DGW,TOTELE*NLOC,'******DGW  ')

       MATSTRVEL=3
       IF(SMALL_DGMAT_MEM) MATSTRVEL=6
       CALL GMRES(DGUVW,VECXYZ,ILE,ILENG,ILE,.FALSE.,&
     &       BIGMELE,FINDRMELE,COLMELE,NCOLMELE,NBIGMELE,&
     & halo_tag,KITS)

       do ELE=1,TOTELE
          do ILOC=1,NLOC
             DGU((ELE-1)*NLOC+ILOC)=DGUVW((ELE-1)*3*NLOC+ILOC)           
             DGV((ELE-1)*NLOC+ILOC)=DGUVW((ELE-1)*3*NLOC+ILOC+NLOC)
             DGW((ELE-1)*NLOC+ILOC)=DGUVW((ELE-1)*3*NLOC+ILOC+2*NLOC)
          END DO
       END DO

       CALL PMINMX(DGU,TOTELE*NLOC,'******DGU  ')
       CALL PMINMX(DGV,TOTELE*NLOC,'******DGV  ')
       CALL PMINMX(DGW,TOTELE*NLOC,'******DGW  ')
       DEALLOCATE(BIGMELE)
       DEALLOCATE(DGUVW)
       DEALLOCATE(VECXYZ)

! Rotate the DGCT matrices...
       CALL ROTDGCT(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,&
     &         C1TDGELE,C2TDGELE,C3TDGELE,&
     &         NONODS,FREDOP,&
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3) 
     
      IF(IUNST_FREE.EQ.1) THEN
! Rotate the DGgT matrices FOR FREE SURFACE...
       CALL ROTDGCT(BIG_NLOC,NLOC,TOTELE,NDGLNO,&
     &         MLOC,PNDGLN,&
     &         NDPSET,D3,&
     &         G1TDGELE,G2TDGELE,G3TDGELE,&
     &         NONODS,FREDOP,&
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD,&
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &         DGDM1,DGDM2,DGDM3) 
      ENDIF
! form r.h.s of pressure eqn...
       ALLOCATE(RHSVEC(FREDOP))
       IF(IUNST_FREE.EQ.0) THEN
         CALL DGCTMULT(RHSVEC,.TRUE.,PNDGLN,MLOC,NLOC,TOTELE,FREDOP,&
     &                 DT,C1TDGELE,C2TDGELE,C3TDGELE,&
     &                 DGU,DGV,DGW)
       ELSE
         ewrite(1,*)'before rhsvec'
         ALLOCATE(DGU_THETA(TOTELE*NLOC))
         ALLOCATE(DGV_THETA(TOTELE*NLOC))
         ALLOCATE(DGW_THETA(TOTELE*NLOC))
         DGU_THETA=DGU+((1.-FTHETA)/FTHETA)*DGUOLD
         DGV_THETA=DGV+((1.-FTHETA)/FTHETA)*DGVOLD
         DGW_THETA=DGW+((1.-FTHETA)/FTHETA)*DGWOLD 
         CALL DGCTMULT(RHSVEC,.TRUE.,PNDGLN,MLOC,NLOC,TOTELE,FREDOP,&
     &                 DT,C1TDGELE,C2TDGELE,C3TDGELE,&
     &                 DGU_THETA,DGV_THETA,DGW_THETA)
         DEALLOCATE(DGU_THETA)
         DEALLOCATE(DGV_THETA)
         DEALLOCATE(DGW_THETA)
! RHSVEC=RHSVEC-(M_S/(GRAVTY*DT*DT))*(P-POLD)
         DO I=1,FREDOP
           DO COUNT=FINCMC(I),FINCMC(I+1)-1
             J=COLCMC(COUNT)
             RHSVEC(I)=RHSVEC(I)-DGCMC(COUNT)*(P(J)-POLD(J))
           END DO
         END DO
       ENDIF
!      IF(IUNST_FREE.EQ.1) THEN
       RHSVEC=RHSVEC-PRES_RHS/DT
         ewrite(1,*) 'after rhsvec'
!
! form DGCMC...
      IF(IUNST_FREE.EQ.0) THEN
        DEALLOCATE(DGCMC)
        ALLOCATE(DGCMC(NCMC))
        DGCMC=0.0
      ENDIF
       ewrite(1,*) 'before getdgcmc fredop,nonods,mloc:',&
     &                          fredop,nonods,mloc
      IF((NBUOY.NE.0).AND.(IMP_P_ROT.EQ.0)) THEN
       CALL GETDGCMC(BIG_NLOC,NLOC,TOTELE,NDGLNO, &
     &         MLOC,PNDGLN, &
     &         NDPSET,D3,DGCMC,.FALSE., &
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE, &
     &         NONODS,FREDOP,NCMC, &
     &         FINCMC,COLCMC,MIDCMC, &
     &         DT,NBUOY, &
     &         BOYMA11,BOYMA12,BOYMA13,BOYMA21, &
     &         BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33,  & 
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD, &
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z, &
     &         DGDM1,DGDM2,DGDM3 )
      ELSE
       CALL GETDGCMC(BIG_NLOC,NLOC,TOTELE,NDGLNO, &
     &         MLOC,PNDGLN, &
     &         NDPSET,D3,DGCMC,.FALSE., &
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE, &
     &         NONODS,FREDOP,NCMC, &
     &         FINCMC,COLCMC,MIDCMC, &
     &         DT,IMP_P_ROT, &
     &         COMAT11,COMAT12,COMAT13,COMAT21,&
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
! The following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD, &
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z, &
     &         DGDM1,DGDM2,DGDM3 )
      ENDIF

! Rotate back U=R^T U...
       CALL UVW_ROT_DG(NLOC,BIG_NLOC,TOTELE,&
     &               .true.,DGU,DGV,DGW,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ)

! Solve pressure eqn...
! region has become isolated. 
       ALLOCATE(DP(FREDOP))
       ewrite(2,*) 'R2NORM(RHSVEC,FREDOP,0):',R2NORM(RHSVEC,FREDOP,0)
       ewrite(2,*) 'R2NORM(DGCMC,NCMC,0):',R2NORM(DGCMC,NCMC,0)
       RMIN=1.0E+20
       ewrite(2,*) 'ncmc,fredop,nnodpp,nonods=',ncmc,fredop,nnodpp,nonods
       do count=1,ncmc
          if((colcmc(count).lt.1).or.(colcmc(count).gt.fredop)) then
             ewrite(-1,*) 'count,colcmc(count):',count,colcmc(count)
             FLAbort('Oops')
          endif
       end do
       do I=1,FREDOP
          RMIN=MIN(RMIN,ABS(DGCMC(MIDCMC(I))))
       END DO

       DP=0.0
       IF(IMP_P_ROT.EQ.0) THEN
          CALL SOLCG(DP,RHSVEC,FREDOP,FREDOP,NNODPP,.FALSE., &
     &            DGCMC,FINCMC,COLCMC,NCMC,NCMC, &
     & halo_tag_p, KITS)
       ELSE
          CALL GMRES(DP,RHSVEC,FREDOP,FREDOP,NNODPP,.FALSE.,&
     &       DGCMC,FINCMC,COLCMC,NCMC,NCMC,&
     &       halo_tag_p,KITS)
       ENDIF
       P=P+DP

! Solve DG mass matrix eqn for velocity update M_DG dU_DG=dt*C*dP
! then choose how to distribute dU_DG between global and SGS...             

! STARTSolve M_DG dDG =dt*C*dP ************************************
       ALLOCATE(COR_DGU(TOTELE*NLOC))
       ALLOCATE(COR_DGV(TOTELE*NLOC))
       ALLOCATE(COR_DGW(TOTELE*NLOC))
      IF((NBUOY.NE.0).AND.(IMP_P_ROT.EQ.0)) THEN
       CALL GET_DG_VEL_FROM_DP(COR_DGU,COR_DGV,COR_DGW, DP, &
     &         BIG_NLOC,NLOC,TOTELE,FREDOP,PNDGLN,MLOC, &
     &         DT,D3, &
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE, &
     &         NONODS, &
     &         NBUOY, &
     &         BOYMA11,BOYMA12,BOYMA13,BOYMA21, &
     &         BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33, &
! Th following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD, &
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z, &
     &         DGDM1,DGDM2,DGDM3 )
      ELSE
       CALL GET_DG_VEL_FROM_DP(COR_DGU,COR_DGV,COR_DGW, DP, &
     &         BIG_NLOC,NLOC,TOTELE,FREDOP,PNDGLN,MLOC, &
     &         DT,D3, &
     &         C1TDGELE,C2TDGELE,C3TDGELE,MASSDGELE, &
     &         NONODS, &
     &         IMP_P_ROT, &
     &         COMAT11,COMAT12,COMAT13,COMAT21, &
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
! Th following is for rotations only ROTAT=.TRUE.
     &         ROTATDG,ROT_THIS_DGNOD, &
     &         DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z, &
     &         DGDM1,DGDM2,DGDM3 )
      ENDIF
       DGU=DGU+COR_DGU
       DGV=DGV+COR_DGV
       DGW=DGW+COR_DGW
! END  Solve M_DG dUDG =dt*C*dP ************************************

! Choose how to distribute dU_DG between global and SGS ********
       U=0.0
       V=0.0
       W=0.0
       USUB=0.0
       VSUB=0.0
       WSUB=0.0
       CALL GET_USUB_FROM_UDG(U,V,W,USUB,VSUB,WSUB, &
     &         DGU,DGV,DGW,&
     &         NLOC,TOTELE,NDGLNO,&
     &         NONODS,FINDRM,COLM,CENTRM,NCOLM,&
     &         DT,D3,&
     &         MASSDGELE,&
     &         0,0,BIGM,&
     &         COMAT11,COMAT12,COMAT13,COMAT21, &
     &         COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &         NNODP,PARA,halo_tag,&
     &         0,0,0,BCU1,BCV1,BCW1,BCU2,BCV2,BCW2)

     END SUBROUTINE DGPROJ

     SUBROUTINE SIMP_DG_MAT_STRU(FINDRMELE,COLMELE,NCOLMELE,TOTELE,NLOC,&
          NDGLNO,CENTRMELE,NONODS)
! Calculate the structure of the element connectivity list...
       INTEGER SNLOC
       PARAMETER(SNLOC=3)
       INTEGER TOTELE,NLOC,NDGLNO(TOTELE*NLOC),NONODS
       INTEGER FINDRMELE(TOTELE+1)
       INTEGER COLMELE(TOTELE*5),CENTRMELE(TOTELE),NCOLMELE
! Local variables...
       INTEGER IELE,JELE,KELE,II,JJ,COUNT,ELE,IL
! Calculate COLMELE points to elements surrounding an element   
       CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLMELE,NONODS)
      
! Place in assending order by bubble sort...
       do IELE=1,TOTELE
          CALL IBUBLE(COLMELE((IELE-1)*5+1),5)
       END DO

! Ignore blanks 0's in COLMELE
       NCOLMELE=0
       do IELE=1,TOTELE
          FINDRMELE(IELE)=NCOLMELE+1
          do II=1,5
             JELE=COLMELE((IELE-1)*5+II)
             IF(JELE.NE.0) THEN
                NCOLMELE=NCOLMELE+1
                COLMELE(NCOLMELE)=JELE
                IF(IELE.EQ.JELE) CENTRMELE(IELE)=NCOLMELE
             ENDIF
          END DO
       END DO
       FINDRMELE(TOTELE+1)=NCOLMELE+1

     END SUBROUTINE SIMP_DG_MAT_STRU
     
     
     

     SUBROUTINE DIFF3DSGS_DG(DGU,DGV,DGW,DGUNEW,DGVNEW,DGWNEW,&
     &     DGNU,DGNV,DGNW,UG,VG,WG,&
     &     SOURCX,SOURCY,SOURCZ,X,Y,Z,D0,&
     &     SUBSOUX,SUBSOUY,SUBSOUZ,ISUBSOU,&
     &     VECX_SGSADD,VECY_SGSADD,VECZ_SGSADD,&
     &     NBUOY,SOXGRA,SOYGRA,SOZGRA, &
     &     NSOGRASUB,SOXGRASUB,SOYGRASUB,SOZGRASUB, &
     &     COMAT11,COMAT12,COMAT13,COMAT21,COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &     BOYMA11,BOYMA12,BOYMA13,BOYMA21,BOYMA22,BOYMA23,BOYMA31,BOYMA32,BOYMA33, &
     &     VECX,VECY,VECZ, &
     &     BIGMELE,SMALL_DGMAT_MEM, MLELE,ML,&
     &     FINDCT,COLCT,NCT,FREDOP,TOTELE, &
     &     FINDRMELE,COLMELE,NCOLMELE,NONODS,CENTRMELE, &
     &     M,MLX,MLY,MLZ,&
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC, BIG_NLOC,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,PNDGLN,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     nnodp,para,halo_tag,&
     &     GEOBAL,SCFACTH0,&
! for the DG...
     &     NBIGMELE,COGRAX,NOBCFSimp,BCT2FSimp, &
     &     MASSDGELE,C1TDGELE,C2TDGELE,C3TDGELE,P,&
     &     FTHETA,GRAVTY,IUNST_FREE,G1TDGELE,G2TDGELE,G3TDGELE,POLD,&
     &     DGCMC,NCMC,FINCMC,COLCMC,MIDCMC, &
     &     SUF_TOP_OCEAN,&
     &     ROTATDG,ROT_THIS_DGNOD,&
     &     DGNX,DGNY,DGNZ, DGT1X,DGT1Y,DGT1Z, DGT2X,DGT2Y,DGT2Z,&
     &     DGDM1,DGDM2,DGDM3,&
     &     PRES_RHS,WEAK_BC,&
     &     NOBCU,NOBCV,NOBCW,&
     &     BCU1,BCV1,BCW1,&
     &     BCU2,BCV2,BCW2, NDPSET, velocity_option_path)
      
!     This subroutine discretises a field equation using Galerkin Least squares. 
!     The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 

! All methods use a balancing diffusion term and THETA time stepping. 
! The following are absolute values of DISOPT...
! Galerkin (146 is defaulted to if non of the below are used): 
! 156 Galerkin LES but no derivatives in the vertical.
! 155 Galerkin LES but no derivatives in the vertical -reduced magnitude by a factor of 10.
! 154 same as 147 and no added viscocity.
! 153 same as 152 and no added viscocity.
! 152 same as 147 but reduced magnitude by a factor of 10.
! 151 is LES applied to the global system. 
! 150 is LES applied to the local Residual system (involving pertation soln) only RECOMMENDED
! 149 is LES applied to the local system only RECOMMENDED
! 148 is LES applied to the global system only 
! 147 is LES applied to every thing. 
! 145 PG applied to just the global part and using global velocity.
! 144 PG applied to just the global part and using the full velocity.
      
!     if DISOPN.ne.0 then set advection to zero and treat advection 
!     using the high res method.
!     IF LUMP then lump the mass matrix else dont. 
!     BETA controls the conservative discretisation 
!     BETA=0. is the usual advection form, BETA=1. is divergence form. 
!     BETA=0.5 is an average of the divergence and advection form.(recommended).
!     ABSORB =absorbtion of energy etc. 
!     MUPTXX,MUPTYY =components of diffusivities. 
!     ML2MXX,ML2MXY etc contain the length scales for LES. 
!     NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
!     UG,VG,WG are the grid velocities. 
!     NDGLNO=element pter for unknowns. 
!     SONDGL=element pter for materials and sources. 
!     VONDGL=element pter for advection NU,UG.
!     XONDGL=element pter for coordinates(x,y,z)
!     If INMAT then return matricie(s). 
!     R0=the minimum distance from ocean bed to centre of earth.
!     D0=H/r0 where H is the height of ocean above R0 (H not used here)
!     For dimensional run D0=0.(R0 is not used used here)
!     NB This sub is only used for field equations in ocean modelling. 
!     ------------------------------------------------------------------------
!     If we are solving for another variable like temperature 
!     or chemical species then NU,NV,NW will be the velocities 
!     and U=TEMPERATURE OR CHEMICAL SPECIES. 
!     ------------------------------------------------------------------------
!     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
       REAL RWIND,INFINY
      LOGICAL ADVECT_MEAN,SUF_CONSERV,CV_PRES
       LOGICAL VERTICAL_BOY
      PARAMETER(VERTICAL_BOY=.true.)
!      PARAMETER(RWIND=0.5,INFINY=1.0E+20,ADVECT_MEAN=.FALSE.,CV_PRES=.false.)
      PARAMETER(RWIND=0.5,INFINY=1.0E+20,ADVECT_MEAN=.TRUE.,CV_PRES=.false.)
      INTEGER SNONOD,VNONOD,XNONOD,NCT,NCOLM,NLOC,BIG_NLOC,NGI
      INTEGER MLOC,DISOPT2,DISOPT,DISOPN,DISCRA
      REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
      REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
      REAL DENPT(SNONOD)
      REAL ABSORB(SNONOD)
      REAL  DT,THETA,BETA
      INTEGER TOTELE,NONODS,FREDOP
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
      REAL DGU(TOTELE*NLOC),DGV(TOTELE*NLOC),DGW(TOTELE*NLOC)
      REAL DGUNEW(TOTELE*NLOC),DGVNEW(TOTELE*NLOC),DGWNEW(TOTELE*NLOC)
      REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
      REAL DGNU(TOTELE*NLOC),DGNV(TOTELE*NLOC),DGNW(TOTELE*NLOC)
      
! This is what we add into the SGS model - a discretised source.      
      REAL VECX_SGSADD(TOTELE*NLOC),VECY_SGSADD(TOTELE*NLOC)
      REAL VECZ_SGSADD(TOTELE*NLOC)
       INTEGER NBUOY,NSOGRASUB
! soxgra,soygra,sozgra is the direction of gravity force when geobal.le.-10.or.(HYDROS.NE.0)
       REAL SOXGRA(SNONOD),SOYGRA(SNONOD),SOZGRA(SNONOD)
       REAL SOXGRASUB(TOTELE*NSOGRASUB),SOYGRASUB(TOTELE*NSOGRASUB)
       REAL SOZGRASUB(TOTELE*NSOGRASUB)
       REAL COMAT11(TOTELE,NLOC),COMAT12(TOTELE,NLOC),COMAT13(TOTELE,NLOC)
       REAL COMAT21(TOTELE,NLOC),COMAT22(TOTELE,NLOC),COMAT23(TOTELE,NLOC)
       REAL COMAT31(TOTELE,NLOC),COMAT32(TOTELE,NLOC),COMAT33(TOTELE,NLOC)
       REAL BOYMA11(NBUOY*TOTELE,NLOC),BOYMA12(NBUOY*TOTELE,NLOC),BOYMA13(NBUOY*TOTELE,NLOC)
       REAL BOYMA21(NBUOY*TOTELE,NLOC),BOYMA22(NBUOY*TOTELE,NLOC),BOYMA23(NBUOY*TOTELE,NLOC)
       REAL BOYMA31(NBUOY*TOTELE,NLOC),BOYMA32(NBUOY*TOTELE,NLOC),BOYMA33(NBUOY*TOTELE,NLOC)
      REAL VECX(TOTELE*NLOC),VECY(TOTELE*NLOC)
      REAL VECZ(TOTELE*NLOC)
      INTEGER NBIGMELE
      REAL BIGMELE(NBIGMELE)
      LOGICAL SMALL_DGMAT_MEM
      LOGICAL COGRAX,WEAK_BC
      INTEGER NOBCFSimp
      INTEGER BCT2FSimp(NOBCFSimp)
      REAL MLELE(TOTELE*NLOC),ML(NONODS)
      REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
      REAL SOURCX(SNONOD),SOURCY(SNONOD),SOURCZ(SNONOD)
      INTEGER ISUBSOU
      REAL SUBSOUX(ISUBSOU*TOTELE*NLOC),SUBSOUY(ISUBSOU*TOTELE*NLOC)
      REAL SUBSOUZ(ISUBSOU*TOTELE*NLOC)
      INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
      INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
      INTEGER XONDGL(TOTELE*NLOC)
      INTEGER NCOLMELE
      INTEGER FINDRMELE(TOTELE+1),COLMELE(NCOLMELE)
      INTEGER CENTRMELE(TOTELE)
      INTEGER FINDCT(FREDOP+1),COLCT(NCT)
      REAL MASSDGELE(TOTELE,NLOC,NLOC)
      REAL C1TDGELE(TOTELE,MLOC,NLOC)
      REAL C2TDGELE(TOTELE,MLOC,NLOC),C3TDGELE(TOTELE,MLOC,NLOC)
      REAL P(FREDOP)
      INTEGER IUNST_FREE
      REAL FTHETA,GRAVTY
      REAL G1TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC)
      REAL G2TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC),G3TDGELE(TOTELE*IUNST_FREE,MLOC,NLOC)
      REAL POLD(FREDOP*IUNST_FREE)
      INTEGER NCMC
      REAL DGCMC(NCMC*IUNST_FREE)
      INTEGER FINCMC(FREDOP+1),COLCMC(NCMC),MIDCMC(FREDOP)
      LOGICAL SUF_TOP_OCEAN
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY     

      REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
      REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
      REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
      LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
      REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
      REAL PRES_RHS(FREDOP)
      INTEGER NOBCU,NOBCV,NOBCW
      REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
      INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
      INTEGER NDPSET
      character(len=*), intent(in):: velocity_option_path

      REAL M(MLOC,NGI)
      REAL MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
      REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL WEIGHT(NGI)
      LOGICAL LUMP,UPWIND,OPWIND,XTPETV
      INTEGER PREOPT
! IF PREOPT=2 then treat every boundary like a free surface - need 
! to do this for open boundaries. 
      LOGICAL ABSLUM,SLUMP
      LOGICAL INMAT
      INTEGER VERSIO
      INTEGER nnodp,para,halo_tag
      INTEGER GEOBAL
      REAL SCFACTH0
      
      REAL MASS,MATRIX,SOUMAT,VLMGI,VLKGI,VABS
      REAL::VLM=0.0
      REAL VF1M,VF2M,VF2,VF3M
      REAL THETA1,THETA2,THETA3,AHAT1,ALPHA,ALPHA2,ALPHA3,THTWDT,INDT,TWOTHI
      REAL PE,HOVERQ,PWEIGI
      REAL GALERK,LEASQR
      REAL UD(NGI),VD(NGI),WD(NGI)
      REAL UDDX(NGI),VDDY(NGI),WDDZ(NGI)
      REAL UDDX2(NGI),VDDY2(NGI),WDDZ2(NGI)
      REAL UGDX(NGI),VGDY(NGI),WGDZ(NGI)
      REAL MUGIXX,MUGIXY,MUGIXZ
      REAL MUGIYY,MUGIYZ,MUGIZZ
      REAL L2GIXX,L2GIXY,L2GIXZ
      REAL L2GIYY,L2GIYZ,L2GIZZ
      REAL DENGI(NGI),DENGI2(NGI),L1GI(NGI)
      REAL ABDGI(NGI)
      REAL AGI,BGI,CGI,DGI
      REAL EGI,FGI,GGI,HGI,KGI
      REAL DETJ,DETWEI(NGI)
      REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
      REAL GIVOLX,GIVOLY,GIVOLZ
!     HX,HY-characteristic length scales in x,y directions.
      REAL HXGI,HYGI,HZGI
      REAL AHAT(NGI),GAMMA(NGI),GAMMAG(NGI),MASCON(NGI)
      REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
      REAL BECTWE(NGI)

      INTEGER ELE,ILOC,JLOC,L
      INTEGER GI,COUNT,I,GLOBI,GLOBJ,GLOBIS,GLOBJS
      INTEGER POS,LOWER,UPPER,PGLOBJ

      REAL D0,DD0,INDD0
      LOGICAL NODIM
      LOGICAL GALWIN
      LOGICAL CTWIND,LES
      REAL ABLUMP,RSLUMP
      REAL R,RUD,RVD,RWD,R1,R2,RR,RMXD,RMXMXX,RMXMYY,RMXMZZ,RMXABS
      REAL::RNN=0.0
      REAL RMXSX,RMXSY,RMXSZ,RMID,RMIMXX,RMIMYY,RMIMZZ,RMIABS
      REAL RMISX,RMISY,RMISZ,RLUMP,RSYM,RCONSI,RGAMMA
      INTEGER IGLS,IGLV,IGLX,INUM
      INTEGER ISPHERE,ICENT
      REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
      REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
      REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
      REAL MYZDTW(NGI),MZZDTW(NGI)
      REAL MUGIXX2,MUGIXY2,MUGIXZ2
      REAL MUGIYY2,MUGIYZ2,MUGIZZ2
      REAL MUPTXXGI,MUPTXYGI,MUPTXZGI,MUPTYYGI,MUPTYZGI,MUPTZZGI
      REAL XD(NGI),YD(NGI),ZD(NGI)
      REAL VLND,VLN2,VLND2,VLN3,VLND3
      REAL TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ
      REAL L2GIXX2,L2GIXY2,L2GIXZ2,L2GIYY2,L2GIYZ2,L2GIZZ2
!  for the extra terms in the adjoint model, local variables....
      REAL VLND_UX_ADJ,VLND_UY_ADJ,VLND_UZ_ADJ,    &
     &       VLND_VX_ADJ,VLND_VY_ADJ,VLND_VZ_ADJ, &
     &       VLND_WX_ADJ,VLND_WY_ADJ,VLND_WZ_ADJ  
      REAL HINX,HINY,HINZ,VLKTEN,NLN,VLMAB
      REAL XNDMUJLOC,YNDMUJLOC,ZNDMUJLOC,XNDMUILOC,YNDMUILOC,ZNDMUILOC
      REAL TEN,TENAB,ABTEN
      REAL VLN,ADVILOC,ADVJLOC,RDIFF
      REAL UATTHETA,VATTHETA,WATTHETA
      LOGICAL BALDIF,BALGLS
      REAL L1(NGI),L2(NGI),L3(NGI),L4(NGI),PGWEIG(NGI) 
      REAL A11(NGI),A12(NGI),A13(NGI)
      REAL A21(NGI),A22(NGI),A23(NGI)
      REAL A31(NGI),A32(NGI),A33(NGI) 
      REAL GAMB(NGI),GAMBAL(NGI)
      REAL ABGI(NGI),DDETWE(NGI)
       
      REAL MASSMM(NLOC,NLOC),MASSINVMM(NLOC,NLOC)
      REAL MASSLUMMM(NLOC,NLOC)
       
      REAL C1TDGELELOC(MLOC,NLOC)
      REAL C2TDGELELOC(MLOC,NLOC),C3TDGELELOC(MLOC,NLOC)
      REAL G1TDGELELOC(MLOC,NLOC)
      REAL G2TDGELELOC(MLOC,NLOC),G3TDGELELOC(MLOC,NLOC)
      REAL TOTSOUX(NLOC),TOTSOUY(NLOC),TOTSOUZ(NLOC)
      INTEGER LESNOD,JNODSUB
      REAL RINMAT,RHSMAT,SRHSMAT
      REAL VLNXN,VLNYN,VLNZN
      REAL VLNNX,VLNNY,VLNNZ
      REAL VLMXN,VLMYN,VLMZN
      REAL C1TCONU,C2TCONV,C3TCONW
! Local variables...
      INTEGER SNLOC,SNGI,SNCLOC
      PARAMETER(SNLOC=3,SNCLOC=6,SNGI=7)
      REAL MX(MLOC,NGI),MY(MLOC,NGI),MZ(MLOC,NGI)
      REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      REAL SNX(SNLOC,SNGI),SNY(SNLOC,SNGI),SNZ(SNLOC,SNGI)
      REAL SWEIGH(SNGI)
      REAL SNC(SNCLOC,SNGI),SNCLX(SNCLOC,SNGI),SNCLY(SNCLOC,SNGI)
      REAL SNCX(SNCLOC,SNGI),SNCY(SNCLOC,SNGI),SNCZ(SNCLOC,SNGI)
      REAL XSL(SNCLOC),YSL(SNCLOC),ZSL(SNCLOC)
      REAL SUD(SNGI),SVD(SNGI),SWD(SNGI)
      REAL IN_SUD(SNGI),IN_SVD(SNGI),IN_SWD(SNGI)
      REAL OUT_SUD(SNGI),OUT_SVD(SNGI),OUT_SWD(SNGI)
      REAL SUDGLOB(SNGI),SVDGLOB(SNGI),SWDGLOB(SNGI)
      REAL NDOTQ(SNGI),NDOTQGLOB(SNGI)
      REAL IN_NDOTQ(SNGI),OUT_NDOTQ(SNGI)
      REAL SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI)
      REAL SDETWE(SNGI)
      REAL SDIRX(SNGI),SDIRY(SNGI),SDIRZ(SNGI)
      REAL SGS_UPWIND_FAC(NGI)
      REAL MASSDGELE5(NLOC+1,NLOC+1)
      REAL MATDIFX5(NLOC+1,NLOC+NLOC),MATDIFY5(NLOC+1,NLOC+NLOC),MATDIFZ5(NLOC+1,NLOC+NLOC)
      REAL DMATDIFX5(NLOC+1,NLOC+NLOC),DMATDIFY5(NLOC+1,NLOC+NLOC),DMATDIFZ5(NLOC+1,NLOC+NLOC)
      REAL MMATDIFX5(NLOC+1,NLOC+NLOC),MMATDIFY5(NLOC+1,NLOC+NLOC),MMATDIFZ5(NLOC+1,NLOC+NLOC)
      REAL MASSMM5(NLOC+1,NLOC+1),MASSINVMM5(NLOC+1,NLOC+1)
      REAL SUFMATX(NLOC,NLOC+1),SUFMATY(NLOC,NLOC+1),SUFMATZ(NLOC,NLOC+1)
      REAL SMMATDIFX5(NLOC,NLOC+NLOC),SMMATDIFY5(NLOC,NLOC+NLOC),SMMATDIFZ5(NLOC,NLOC+NLOC)
      REAL SMMATDIF5(NLOC,NLOC+NLOC)
      REAL BIG_AMAT(BIG_NLOC,BIG_NLOC),FAMAT(NLOC,5*NLOC),AMAT(NLOC,NLOC)
      REAL MASSLOC(NLOC,NLOC),MASSLOCINV(NLOC,NLOC)
      REAL LOCDGNUDX2(NLOC),LOCDGNUDY2(NLOC),LOCDGNUDZ2(NLOC)
      REAL LOCDGNVDX2(NLOC),LOCDGNVDY2(NLOC),LOCDGNVDZ2(NLOC)
      REAL LOCDGNWDX2(NLOC),LOCDGNWDY2(NLOC),LOCDGNWDZ2(NLOC)
         
      REAL LOCDGNUDX(NLOC),LOCDGNUDY(NLOC),LOCDGNUDZ(NLOC)
      REAL LOCDGNVDX(NLOC),LOCDGNVDY(NLOC),LOCDGNVDZ(NLOC)
      REAL LOCDGNWDX(NLOC),LOCDGNWDY(NLOC),LOCDGNWDZ(NLOC)
      INTEGER IFACE,NFACE
      PARAMETER(NFACE=4)
      INTEGER LOCLIST(NFACE,3),SILOC2ILOC(SNLOC),ILOC_OTHER_SIDE(SNLOC)
      INTEGER ELE2_LOC_NODS(NLOC)
      INTEGER SILOC,SJLOC,SKLOC,IGL,COL,INOD,JNOD,SGI,KLOC
      INTEGER BIG_ILOC,BIG_JLOC,ELE2,IJDISP,ID,JD,IXNOD,NODDG
      INTEGER IDGNOD,NOD,II,ICOUNT,PGLOBI,ILOC2,COUNTELE,JJLOC5,IILOC5
      INTEGER::IDG,ILOCOTH,ILOCELE,INOD2,JDG2,JDG,JELE,JLOC2,JLOCELE=-1,IDG2
      INTEGER ICOUNTU,ICOUNTV,ICOUNTW,ICOUNTUVW,JGL
      LOGICAL OUTLET_ADV,INLET_ADV,DIFF_BC,XDOM,YDOM,ZDOM
      REAL RADI,RADJ,RADP,RAD,RNU,RNV,RNW,SAREA,VOL
      REAL NORMX,NORMY,NORMZ,VLNDOTQ,allVLNDOTQ,STABINEL,UPWIN,TDIVU
      REAL VLOMEGA,OMEGAGI,STABOME
      REAL A11MAT,A12MAT,A21MAT,A22MAT,O11MAT,O12MAT,O22MAT,O21MAT
      REAL RO11MAT,RO12MAT,RO22MAT,RO21MAT
      REAL RNUSUBGI,RNVSUBGI,RNWSUBGI
      REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
      REAL AVE_SNORMXN,AVE_SNORMYN,AVE_SNORMZN,ONE_OR_ZERO,ROCEAN
      REAL LOCPNORMX,LOCPNORMY,LOCPNORMZ
      REAL XC,YC,ZC,RN,sum,stVLKTEN,stTEN
      REAL ALPHA11,ALPHA12,ALPHA21,ALPHA22,SUMP,SUMU
      REAL RLES,RHAVEVIS,VLNDOTQ_IN,VLNDOTQ_OUT,VLNDOTQ_ALL
      REAL SMEANXX,SMEANXY,SMEANXZ,SMEANYY,SMEANYZ,SMEANZZ
      REAL THETA_FAMAT,VLM_NORX,VLM_NORY,VLM_NORZ
      REAL MNNX,MNNY,MNNZ,NN,FORU,FORV,FORW,RNXN,RNYN,RNZN
      REAL RPRES
      LOGICAL DMATLUM,SUF_TOP,STREAM_UPWIND,BNDSUF
           REAL BOY_ABSXX(NGI),BOY_ABSXY(NGI),BOY_ABSXZ(NGI)
           REAL BOY_ABSYX(NGI),BOY_ABSYY(NGI),BOY_ABSYZ(NGI)
           REAL BOY_ABSZX(NGI),BOY_ABSZY(NGI),BOY_ABSZZ(NGI)
           REAL RSOXGRAGI,RSOYGRAGI,RSOZGRAGI
           REAL GRAV_DEN_NOD_NX,GRAV_DEN_NOD_NY,GRAV_DEN_NOD_NZ
           REAL GRAV_DEN_NOD,GRA_SIG_CHAN
           REAL OVER_RELAX_BOY
! Local coords...
      INTEGER NCLOC,SMLOC,NBLOCK
      REAL, ALLOCATABLE, DIMENSION(:,:)::NC
      REAL, ALLOCATABLE, DIMENSION(:,:)::NCLX
      REAL, ALLOCATABLE, DIMENSION(:,:)::NCLY
      REAL, ALLOCATABLE, DIMENSION(:,:)::NCLZ
      REAL, ALLOCATABLE, DIMENSION(:,:)::NCX
      REAL, ALLOCATABLE, DIMENSION(:,:)::NCY
      REAL, ALLOCATABLE, DIMENSION(:,:)::NCZ
      REAL, ALLOCATABLE, DIMENSION(:)::VEDUDX
      REAL, ALLOCATABLE, DIMENSION(:)::VEDUDY
      REAL, ALLOCATABLE, DIMENSION(:)::VEDUDZ
      REAL, ALLOCATABLE, DIMENSION(:)::VEDVDX
      REAL, ALLOCATABLE, DIMENSION(:)::VEDVDY
      REAL, ALLOCATABLE, DIMENSION(:)::VEDVDZ
      REAL, ALLOCATABLE, DIMENSION(:)::VEDWDX
      REAL, ALLOCATABLE, DIMENSION(:)::VEDWDY
      REAL, ALLOCATABLE, DIMENSION(:)::VEDWDZ

      REAL, ALLOCATABLE, DIMENSION(:)::HIMATX
      REAL, ALLOCATABLE, DIMENSION(:)::HIMATY
      REAL, ALLOCATABLE, DIMENSION(:)::HIMATZ
      REAL, ALLOCATABLE, DIMENSION(:)::HIMASS
 
      REAL, ALLOCATABLE, DIMENSION(:)::VECU
      REAL, ALLOCATABLE, DIMENSION(:)::VECV
      REAL, ALLOCATABLE, DIMENSION(:)::VECW

      REAL, ALLOCATABLE, DIMENSION(:)::VECYT

      REAL, ALLOCATABLE, DIMENSION(:,:,:)::EV
      INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
      INTEGER, ALLOCATABLE, DIMENSION(:)::COLELE4
         
      REAL, ALLOCATABLE, DIMENSION(:)::MASSBIGM
      REAL, ALLOCATABLE, DIMENSION(:)::MASSVEC
         
      REAL, ALLOCATABLE, DIMENSION(:,:)::SNSUB
      INTEGER, ALLOCATABLE, DIMENSION(:)::MLOCSILOC2ILOC 
      INTEGER, ALLOCATABLE, DIMENSION(:)::P_MLOCSILOC2ILOC
      REAL, ALLOCATABLE, DIMENSION(:,:)::SM
      REAL, ALLOCATABLE, DIMENSION(:,:)::ctemp
      INTEGER, ALLOCATABLE, DIMENSION(:)::BNDCON
      REAL, ALLOCATABLE, DIMENSION(:)::VOLELE
      REAL, ALLOCATABLE, DIMENSION(:)::MEANXX,MEANXY,MEANXZ,MEANYY,MEANYZ,MEANZZ
      REAL, ALLOCATABLE, DIMENSION(:,:,:)::MATDIFX,MATDIFY,MATDIFZ,MATVLK
      INTEGER, ALLOCATABLE, DIMENSION(:)::DGMARKU,DGMARKV,DGMARKW
      REAL, ALLOCATABLE, DIMENSION(:)::DIRI_BC_X,DIRI_BC_Y,DIRI_BC_Z
      REAL, ALLOCATABLE, DIMENSION(:)::TESTX,TESTY,TESTZ
! BNDCON marks the boundary condition nodes...
      ALLOCATE(BNDCON(NONODS))

      ewrite(1, *) "in subroutine diff3dsgs_dg"
      ewrite(2, *) "isphere=",isphere

      BNDCON(1:NONODS)=0
      do II=1,NOBCFSimp
         GLOBI=BCT2FSimp(II)
         BNDCON(GLOBI)=1
      END DO
         
      ALLOCATE(DGMARKU(NLOC*TOTELE))
      ALLOCATE(DGMARKV(NLOC*TOTELE))
      ALLOCATE(DGMARKW(NLOC*TOTELE))
      ALLOCATE(DIRI_BC_X(NLOC*TOTELE))
      ALLOCATE(DIRI_BC_Y(NLOC*TOTELE))
      ALLOCATE(DIRI_BC_Z(NLOC*TOTELE))
      CALL DGDM123_DG_WEAK(&
     &                     NLOC,TOTELE,NDGLNO,NONODS,&
     &                     DGU,DGV,DGW,DT,&
     &                     NOBCU,NOBCV,NOBCW,&
     &                     BCU1,BCV1,BCW1,&
     &                     BCU2,BCV2,BCW2,&
     &                     DGMARKU,DGMARKV,DGMARKW,&
     &                     DIRI_BC_X,DIRI_BC_Y,DIRI_BC_Z)
 

      SUF_CONSERV=(ABS(BETA).GT.1.E-3)
      
      DISOPT=DISOPT2
      IF(DISOPT2.EQ.47) DISOPT=147
      IF(DISOPT2.EQ.48) DISOPT=148

      NCLOC=4
      IF(ISPHERE.EQ.2) NCLOC=10

      ALLOCATE(NC(NCLOC,NGI))
      ALLOCATE(NCLX(NCLOC,NGI))
      ALLOCATE(NCLY(NCLOC,NGI))
      ALLOCATE(NCLZ(NCLOC,NGI))
      ALLOCATE(NCX(NCLOC,NGI))
      ALLOCATE(NCY(NCLOC,NGI))
      ALLOCATE(NCZ(NCLOC,NGI))
       
      IF(ISPHERE.NE.0) THEN
         ASSERT(ncloc==10)
          
         call QuadTetShapes(nc, nclx, ncly, nclz, ncloc, ngi)
         
      ELSE
         NC = N
         NCLX = NLX
         NCLY = NLY
         NCLZ = NLZ
      ENDIF
       
      COMAT11=0.0
      COMAT12=0.0
      COMAT13=0.0
      COMAT21=0.0
      COMAT22=0.0
      COMAT23=0.0
      COMAT31=0.0
      COMAT32=0.0
      COMAT33=0.0
      
      BOYMA11=0.0
      BOYMA12=0.0
      BOYMA13=0.0
      BOYMA21=0.0
      BOYMA22=0.0
      BOYMA23=0.0
      BOYMA31=0.0
      BOYMA32=0.0
      BOYMA33=0.0
      
       
      ML=0.0
      MLELE=0.0
      do I=1,TOTELE*NLOC
         VECX(I)=0.
         VECY(I)=0.
         VECZ(I)=0.
      end do
      BIGMELE(1:NBIGMELE) = 0.0
       
      PRES_RHS(1:FREDOP)=0.0

      RLUMP=1.
      IF(LUMP) RLUMP=1.
      RCONSI=1.-RLUMP

      ABLUMP=1.
      RSLUMP=1.
      IF(ABSLUM) ABLUMP=1.
      IF(SLUMP)  RSLUMP=1.
! ALPHA contains the inner element stabilisation that ensures 
! a non-singular system...

!     ALPHA=0.01 
     ALPHA=0.0


     ewrite(3,*) "InnerElement: Filter section"
     if (have_option(trim(velocity_option_path) // "/prognostic/&
          &spatial_discretisation/inner_element")) then
       ewrite(3,*) "InnerElement: Filter section - detected model"
       if (have_option(trim(velocity_option_path) // "/prognostic/&
            &spatial_discretisation/inner_element/&
            &use_filter")) then
          ewrite(3,*) "InnerElement: SGS filtering turned on"
          call get_option(trim(velocity_option_path) // "/prognostic/&
            &spatial_discretisation/inner_element/&
            &use_filter/strength", alpha, default=0.0) 
       else
          ewrite(3,*) "InnerElement: SGS filtering turned off"
          alpha=0.0
       end if
     end if
     
    ewrite(1,*) "InnerElement: Filter section, alpha[",alpha,']'


! alpha=1.0 is good or DMATLUM=.TRUE. this way it works without 
! transport
      DMATLUM=.false.

! alpha2 control works for void but others do not.
      ALPHA2=0.0
      
      ALPHA3=0.0
      ROCEAN=0.0
      IF(SUF_TOP_OCEAN) ROCEAN=1.0
! DISOPT=151 is LES applied to the global system. 
! DISOPT=150 is LES applied to the local Residual system (involving pertation soln) only RECOMMENDED
! DISOPT=149 is LES applied to the local system only RECOMMENDED
! 150 and 149 produce the same results.
!      
      RLES=1.0
      IF((DISOPT.EQ.152).OR.(DISOPT.EQ.153).OR.(DISOPT.EQ.155)) RLES=0.1
      RHAVEVIS=1.0
      IF((DISOPT.EQ.153).OR.(DISOPT.EQ.154)) RHAVEVIS=0.0

      NBLOCK=(NLOC*3)**2

! Calculate surface shape functions SN, SM etc...
      CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
      CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
     &             SNLOC,SNGI, &
     &             SN,SNLX,SNLY,SNLX)
    
      CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
     &             SNCLOC,SNGI, &
     &             SNC,SNCLX,SNCLY,SNCLX)
     
      ALLOCATE(MLOCSILOC2ILOC(SNLOC))
       
      IF(MLOC.EQ.10) THEN
! Quadratic...
         SMLOC=6
         ALLOCATE(SM(SMLOC,SNGI))
         SM=SNC
      ELSE
! linear...
         SMLOC=3
         ALLOCATE(SM(SMLOC,SNGI))
         SM=SN
      ENDIF
      ALLOCATE(P_MLOCSILOC2ILOC(SMLOC))

      ALLOCATE(FINELE(TOTELE+1))
      ALLOCATE(COLELE4(TOTELE*5))
      ewrite(1,*) 'in subroutine DIFF3DSGS_DG going into GETFINELE4'
      CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE4,NONODS)
      ewrite(1,*) 'finished GETFINELE4'
      
      ROT_THIS_DGNOD=.FALSE.
      DGDM1=0.0
      DGDM2=0.0
      DGDM3=0.0
        
      DGT1X=1.0
      DGT1Y=0.0
      DGT1Z=0.0
      DGT2X=0.0
      DGT2Y=1.0
      DGT2Z=0.0 
      DGNX =0.0
      DGNY =0.0
      DGNZ =1.0

! The local cords are in anti-clockwise order
      IFACE=1
      LOCLIST(IFACE,1)=1
      LOCLIST(IFACE,2)=2
      LOCLIST(IFACE,3)=3
      IFACE=2
      LOCLIST(IFACE,1)=1
      LOCLIST(IFACE,2)=2
      LOCLIST(IFACE,3)=4
      IFACE=3
      LOCLIST(IFACE,1)=3
      LOCLIST(IFACE,2)=2
      LOCLIST(IFACE,3)=4
      IFACE=4
      LOCLIST(IFACE,1)=1
      LOCLIST(IFACE,2)=3
      LOCLIST(IFACE,3)=4
! calculate approximate volume of eache element...

      IF(IUNST_FREE.EQ.1) THEN
        DGCMC=0.0
      ENDIF

! Form critical matrices for DG diffusion:
      ALLOCATE(MATDIFX(TOTELE,NLOC,NLOC))
      ALLOCATE(MATDIFY(TOTELE,NLOC,NLOC))
      ALLOCATE(MATDIFZ(TOTELE,NLOC,NLOC))
      ALLOCATE(MATVLK(TOTELE,NLOC,NLOC))
      ALLOCATE(VOLELE(TOTELE))
      ALLOCATE(MEANXX(TOTELE))
      ALLOCATE(MEANXY(TOTELE))
      ALLOCATE(MEANXZ(TOTELE))
      ALLOCATE(MEANYY(TOTELE))
      ALLOCATE(MEANYZ(TOTELE))
      ALLOCATE(MEANZZ(TOTELE))
       
      MASSDGELE=0.0
       
      VOLELE=0.0
      MEANXX=0.0
      MEANXY=0.0
      MEANXZ=0.0
      MEANYY=0.0
      MEANYZ=0.0
      MEANZZ=0.0
      ewrite(1,*) 'in subroutine DIFF3DSGS_DG going into 1st element loop'
       
      element_loop: do  ELE=1,TOTELE

! Get determinant and derivatives...
         CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
     &               N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
     &               NX,NY,NZ, NCX,NCY,NCZ,&
     &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
     &               XD,YD,ZD,&
     &               ISPHERE) 
     
          DO GI=1,NGI
             IF(COGRAX) THEN
               DIRX(GI)=0.0
               DIRY(GI)=0.0
               DIRZ(GI)=1.0
             ELSE
               DIRX(GI)=XD(GI)
               DIRY(GI)=YD(GI)
               DIRZ(GI)=ZD(GI)
               RN=SQRT(DIRX(GI)**2+DIRY(GI)**2+DIRZ(GI)**2)
               DIRX(GI)=DIRX(GI)/RN
               DIRY(GI)=DIRY(GI)/RN
               DIRZ(GI)=DIRZ(GI)/RN
             ENDIF
          END DO
     
! Calculate the derivatives...
         MASSLOC=0.0
           
         LOCDGNUDX2=0.0
         LOCDGNUDY2=0.0
         LOCDGNUDZ2=0.0
            
         LOCDGNVDX2=0.0
         LOCDGNVDY2=0.0
         LOCDGNVDZ2=0.0
            
         LOCDGNWDX2=0.0
         LOCDGNWDY2=0.0
         LOCDGNWDZ2=0.0
           
         do JLOC=1,NLOC
            JDG=(ELE-1)*NLOC+JLOC
            JGL=NDGLNO((ELE-1)*NLOC+JLOC)
            do ILOC=1,NLOC
               IDG=(ELE-1)*NLOC+ILOC
               do GI=1,NGI
                  MASSLOC(ILOC,JLOC)=MASSLOC(ILOC,JLOC)+N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                  RNXN=NX(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                  RNYN=NY(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                  RNZN=NZ(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                  LOCDGNUDX2(ILOC)=LOCDGNUDX2(ILOC)-RNXN*(DGNU(JDG)-UG(JGL))
                  LOCDGNUDY2(ILOC)=LOCDGNUDY2(ILOC)-RNYN*(DGNU(JDG)-UG(JGL))
                  LOCDGNUDZ2(ILOC)=LOCDGNUDZ2(ILOC)-RNZN*(DGNU(JDG)-UG(JGL))
            
                  LOCDGNVDX2(ILOC)=LOCDGNVDX2(ILOC)-RNXN*(DGNV(JDG)-VG(JGL))
                  LOCDGNVDY2(ILOC)=LOCDGNVDY2(ILOC)-RNYN*(DGNV(JDG)-VG(JGL))
                  LOCDGNVDZ2(ILOC)=LOCDGNVDZ2(ILOC)-RNZN*(DGNV(JDG)-VG(JGL))
            
                  LOCDGNWDX2(ILOC)=LOCDGNWDX2(ILOC)-RNXN*(DGNW(JDG)-WG(JGL))
                  LOCDGNWDY2(ILOC)=LOCDGNWDY2(ILOC)-RNYN*(DGNW(JDG)-WG(JGL))
                  LOCDGNWDZ2(ILOC)=LOCDGNWDZ2(ILOC)-RNZN*(DGNW(JDG)-WG(JGL))
               END DO
            END DO
         END DO

! calculate SNX,SNY,SNZ***
! Surface element of domain...
         face_loop: do  IFACE=1,NFACE
         
! NB colele4 is ordered in terms of faces.
            ELE2=COLELE4((ELE-1)*5+IFACE)
         
            do SILOC=1,SNLOC
               SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
            END DO
            IF(ELE2.NE.0) THEN
! Calculate the nodes on the other side of the face:
               do SILOC=1,SNLOC
                  ILOC=SILOC2ILOC(SILOC)
                  INOD=XONDGL((ELE-1)*NLOC+ILOC)
                  do ILOC2=1,NLOC
                     INOD2=XONDGL((ELE2-1)*NLOC+ILOC2)
                     IF(INOD2.EQ.INOD) ILOC_OTHER_SIDE(SILOC)=ILOC2
                  END DO
               END DO
            ENDIF

            do SILOC=1,SNLOC
               ILOC=SILOC2ILOC(SILOC)
               INOD=XONDGL((ELE-1)*NLOC+ILOC)
               do SJLOC=SILOC,SNLOC,1
                  JLOC=SILOC2ILOC(SJLOC)
                  JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                  SKLOC=SILINK2(SILOC,SJLOC)
                  XSL(SKLOC)=0.5*(X(INOD)+X(JNOD))
                  YSL(SKLOC)=0.5*(Y(INOD)+Y(JNOD))
                  ZSL(SKLOC)=0.5*(Z(INOD)+Z(JNOD))
! Assume the mide side nodes are on the sphere
                  IF(ISPHERE.GE.2) THEN
                     IF(ILOC.NE.JLOC) THEN
                        RADI=SQRT(X(INOD)**2+Y(INOD)**2+Z(INOD)**2)
                        RADJ=SQRT(X(JNOD)**2+Y(JNOD)**2+Z(JNOD)**2)
                        RADP=0.5*(RADI+RADJ)
                        RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                        XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                        YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                        ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
                     ENDIF
                  ENDIF
               END DO
            END DO
       
! Form approximate surface normal (NORMX,NORMY,NORMZ)
            CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
     &                     X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)
! Recalculate the normal...
            CALL DGSDETNXLOC2(SNCLOC,SNGI, &
     &        XSL,YSL,ZSL, &
     &        SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
     &        SNORMXN,SNORMYN,SNORMZN, &
     &        NORMX,NORMY,NORMZ)

            do SJLOC=1,SNLOC
               JLOC=SILOC2ILOC(SJLOC)
               JGL =NDGLNO((ELE-1)*NLOC+JLOC)
               JDG=(ELE-1)*NLOC+JLOC
               IF(ELE2.EQ.0) THEN
                  JDG2=JDG
               ELSE
                  JLOC2=ILOC_OTHER_SIDE(SJLOC)
                  JDG2=(ELE2-1)*NLOC+JLOC2
               ENDIF
               FORU=(0.5*(DGNU(JDG)+DGNU(JDG2))-UG(JGL))
               FORV=(0.5*(DGNV(JDG)+DGNV(JDG2))-VG(JGL))
               FORW=(0.5*(DGNW(JDG)+DGNW(JDG2))-WG(JGL))
               do SILOC=1,SNLOC
                  ILOC=SILOC2ILOC(SILOC)
                  IGL =NDGLNO((ELE-1)*NLOC+ILOC)
                  IDG=(ELE-1)*NLOC+ILOC
                  do SGI=1,SNGI
                     RN=SN(SILOC,SGI)*SN(SJLOC,SGI)*SDETWE(SGI)
                     LOCDGNUDX2(ILOC)=LOCDGNUDX2(ILOC) + RN*SNORMXN(SGI)*FORU
                     LOCDGNUDY2(ILOC)=LOCDGNUDY2(ILOC) + RN*SNORMYN(SGI)*FORU
                     LOCDGNUDZ2(ILOC)=LOCDGNUDZ2(ILOC) + RN*SNORMZN(SGI)*FORU
                 
                     LOCDGNVDX2(ILOC)=LOCDGNVDX2(ILOC) + RN*SNORMXN(SGI)*FORV
                     LOCDGNVDY2(ILOC)=LOCDGNVDY2(ILOC) + RN*SNORMYN(SGI)*FORV
                     LOCDGNVDZ2(ILOC)=LOCDGNVDZ2(ILOC) + RN*SNORMZN(SGI)*FORV
                 
                     LOCDGNWDX2(ILOC)=LOCDGNWDX2(ILOC) + RN*SNORMXN(SGI)*FORW
                     LOCDGNWDY2(ILOC)=LOCDGNWDY2(ILOC) + RN*SNORMYN(SGI)*FORW
                     LOCDGNWDZ2(ILOC)=LOCDGNWDZ2(ILOC) + RN*SNORMZN(SGI)*FORW
                  END DO
               END DO
            END DO

! calculate MASSLOC...
            CALL MATDMATINV(MASSLOC,MASSLOCINV,NLOC)
            LOCDGNUDX=0.0
            LOCDGNUDY=0.0
            LOCDGNUDZ=0.0
            
            LOCDGNVDX=0.0
            LOCDGNVDY=0.0
            LOCDGNVDZ=0.0
            
            LOCDGNWDX=0.0
            LOCDGNWDY=0.0
            LOCDGNWDZ=0.0
            do ILOC=1,NLOC
               do JLOC=1,NLOC
                  LOCDGNUDX(ILOC)=LOCDGNUDX(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNUDX2(JLOC)
                  LOCDGNUDY(ILOC)=LOCDGNUDY(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNUDY2(JLOC)
                  LOCDGNUDZ(ILOC)=LOCDGNUDZ(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNUDZ2(JLOC)
                 
                  LOCDGNVDX(ILOC)=LOCDGNVDX(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNVDX2(JLOC)
                  LOCDGNVDY(ILOC)=LOCDGNVDY(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNVDY2(JLOC)
                  LOCDGNVDZ(ILOC)=LOCDGNVDZ(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNVDZ2(JLOC)
                 
                  LOCDGNWDX(ILOC)=LOCDGNWDX(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNWDX2(JLOC)
                  LOCDGNWDY(ILOC)=LOCDGNWDY(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNWDY2(JLOC)
                  LOCDGNWDZ(ILOC)=LOCDGNWDZ(ILOC) + MASSLOCINV(ILOC,JLOC)*LOCDGNWDZ2(JLOC)
               END DO
            END DO

         END DO face_loop

         gauss_loop: do  GI=1,NGI
            VOLELE(ELE)=VOLELE(ELE)+DETWEI(GI)
             
            MUPTXXGI=0.0
            MUPTXYGI=0.0
            MUPTXZGI=0.0
            MUPTYYGI=0.0
            MUPTYZGI=0.0
            MUPTZZGI=0.0
            do  L=1,NLOC
               IGLV=VONDGL((ELE-1)*NLOC+L)
      
               MUPTXXGI=MUPTXXGI+N(L,GI)*MUPTXX(IGLV)
               MUPTXYGI=MUPTXYGI+N(L,GI)*MUPTXY(IGLV)
               MUPTXZGI=MUPTXZGI+N(L,GI)*MUPTXZ(IGLV)
               MUPTYYGI=MUPTYYGI+N(L,GI)*MUPTYY(IGLV)
               MUPTYZGI=MUPTYZGI+N(L,GI)*MUPTYZ(IGLV)
               MUPTZZGI=MUPTZZGI+N(L,GI)*MUPTZZ(IGLV)
            end do

! calculate tensor for oceans:
! this has a QUADRATURE POINT varying element size and shape
            CALL SIZEGIELETENS(NX,NY,NZ,NLOC,NGI,GI,&
     &     L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ )

            IF((DISOPT.GE.147).AND.(DISOPT.LE.156)) THEN

               MXXDTW(GI)=L2GIXX
               MXYDTW(GI)=L2GIXY
               MXZDTW(GI)=L2GIXZ
               MYYDTW(GI)=L2GIYY
               MYZDTW(GI)=L2GIYZ
               MZZDTW(GI)=L2GIZZ
               CALL DG_DG_LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
     &             N,NX,NY,NZ, &
     &             DGNU,DGNV,DGNW, &
     &             MXXDTW,MXYDTW,MXZDTW,&
     &             MYYDTW,MYZDTW,MZZDTW,&
     &             LOCDGNUDX,LOCDGNUDY,LOCDGNUDZ,&
     &             LOCDGNVDX,LOCDGNVDY,LOCDGNVDZ,&
     &             LOCDGNWDX,LOCDGNWDY,LOCDGNWDZ, &
     &             (DISOPT.EQ.155).OR.(DISOPT.EQ.156),DIRX,DIRY,DIRZ)

               MUGIXX=MXXDTW(GI)
               MUGIXY=MXYDTW(GI)
               MUGIXZ=MXZDTW(GI)
               MUGIYY=MYYDTW(GI)
               MUGIYZ=MYZDTW(GI)
               MUGIZZ=MZZDTW(GI)
            ELSE

               MUGIXX=0.0
               MUGIXY=0.0
               MUGIXZ=0.0
               MUGIYY=0.0
               MUGIYZ=0.0
               MUGIZZ=0.0
            ENDIF
       
! 0.25 is used to scale down LES dissipation from ocean dissipation
! factor of 4 is used to reduce the length scale. 
            MXXDTW(GI)=(RLES*MUGIXX+RHAVEVIS*MUPTXXGI)*DETWEI(GI)
            MXYDTW(GI)=(RLES*MUGIXY+RHAVEVIS*MUPTXYGI)*DETWEI(GI)
            MXZDTW(GI)=(RLES*MUGIXZ+RHAVEVIS*MUPTXZGI)*DETWEI(GI)

            MYYDTW(GI)=(RLES*MUGIYY+RHAVEVIS*MUPTYYGI)*DETWEI(GI)
            MYZDTW(GI)=(RLES*MUGIYZ+RHAVEVIS*MUPTYZGI)*DETWEI(GI)

            MZZDTW(GI)=(RLES*MUGIZZ+RHAVEVIS*MUPTZZGI)*DETWEI(GI)

            MEANXX(ELE)=MEANXX(ELE)+MXXDTW(GI)
            MEANXY(ELE)=MEANXY(ELE)+MXYDTW(GI)
            MEANXZ(ELE)=MEANXZ(ELE)+MXZDTW(GI)
            MEANYY(ELE)=MEANYY(ELE)+MYYDTW(GI)
            MEANYZ(ELE)=MEANYZ(ELE)+MYZDTW(GI)
            MEANZZ(ELE)=MEANZZ(ELE)+MZZDTW(GI)
         end do gauss_loop

         MEANXX(ELE)=MEANXX(ELE)/VOLELE(ELE)
         MEANXY(ELE)=MEANXY(ELE)/VOLELE(ELE)
         MEANXZ(ELE)=MEANXZ(ELE)/VOLELE(ELE)
         MEANYY(ELE)=MEANYY(ELE)/VOLELE(ELE)
         MEANYZ(ELE)=MEANYZ(ELE)/VOLELE(ELE)
         MEANZZ(ELE)=MEANZZ(ELE)/VOLELE(ELE)

! now form matrices: MATDIFX,MATDIFY,MATDIFZ,MASSDGELE

         do ILOC=1,NLOC
            do JLOC=1,NLOC
               NN=0.0
               MNNX=0.0
               MNNY=0.0
               MNNZ=0.0
               VLKTEN=0.0
               do GI=1,NGI
                  NN=NN+N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                  MNNX=MNNX+N(ILOC,GI)*(MXXDTW(GI)*NX(JLOC,GI)+MXYDTW(GI)*NY(JLOC,GI)&
     &                            +MXZDTW(GI)*NZ(JLOC,GI))
                  MNNY=MNNY+N(ILOC,GI)*(MXYDTW(GI)*NX(JLOC,GI)+MYYDTW(GI)*NY(JLOC,GI)&
     &                            +MYZDTW(GI)*NZ(JLOC,GI))
                  MNNZ=MNNZ+N(ILOC,GI)*(MXZDTW(GI)*NX(JLOC,GI)+MYZDTW(GI)*NY(JLOC,GI)&
     &                            +MZZDTW(GI)*NZ(JLOC,GI))
                  XNDMUJLOC=(NX(JLOC,GI)*MXXDTW(GI)&
     &                 +NY(JLOC,GI)*MXYDTW(GI)+NZ(JLOC,GI)*MXZDTW(GI))
                  YNDMUJLOC=(NX(JLOC,GI)*MXYDTW(GI)&
     &                 +NY(JLOC,GI)*MYYDTW(GI)+NZ(JLOC,GI)*MYZDTW(GI))
                  ZNDMUJLOC=(NX(JLOC,GI)*MXZDTW(GI)&
     &                 +NY(JLOC,GI)*MYZDTW(GI)+NZ(JLOC,GI)*MZZDTW(GI))
     
                  TEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC+NZ(ILOC,GI)*ZNDMUJLOC
                  VLKTEN=VLKTEN+TEN
               END DO
               MASSDGELE(ELE,ILOC,JLOC)=NN
               MATDIFX(ELE,ILOC,JLOC)=MNNX
               MATDIFY(ELE,ILOC,JLOC)=MNNY
               MATDIFZ(ELE,ILOC,JLOC)=MNNZ
               MATVLK(ELE,ILOC,JLOC)=VLKTEN
            END DO
         END DO
             
      END DO element_loop
      
      ewrite(1,*) 'in subroutine DIFF3DSGS_DG going into 2nd element loop'

      second_element_loop: do ELE=1,TOTELE

! Get determinant and derivatives...
         CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
     &               N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
     &               NX,NY,NZ, NCX,NCY,NCZ,&
     &               A11,A12,A13, A21,A22,A23, A31,A32,A33,&
     &               XD,YD,ZD,&
     &               ISPHERE) 
     
         do ILOC=1,MLOC
            do GI=1,NGI
               MX(ILOC,GI)= A11(GI)*MLX(ILOC,GI)+A12(GI)*MLY(ILOC,GI)+A13(GI)*MLZ(ILOC,GI)
               MY(ILOC,GI)= A21(GI)*MLX(ILOC,GI)+A22(GI)*MLY(ILOC,GI)+A23(GI)*MLZ(ILOC,GI)
               MZ(ILOC,GI)= A31(GI)*MLX(ILOC,GI)+A32(GI)*MLY(ILOC,GI)+A33(GI)*MLZ(ILOC,GI)
            END DO
         END DO
! ***PUT ALL SOURCES IN TOTSOU*****
         do ILOC=1,NLOC
            IGL=SONDGL((ELE-1)*NLOC+ILOC)
            IF(ISUBSOU.EQ.0) THEN
               TOTSOUX(ILOC)=SOURCX(IGL)
               TOTSOUY(ILOC)=SOURCY(IGL)
               TOTSOUZ(ILOC)=SOURCZ(IGL)
            ELSE
               NODDG=(ELE-1)*NLOC+ILOC
               TOTSOUX(ILOC)=SOURCX(IGL)+SUBSOUX(NODDG)
               TOTSOUY(ILOC)=SOURCY(IGL)+SUBSOUY(NODDG)
               TOTSOUZ(ILOC)=SOURCZ(IGL)+SUBSOUZ(NODDG)
            ENDIF
         END DO
         
         second_gauss_loop: do GI=1,NGI

            UD(GI)=0.    
            VD(GI)=0.
            WD(GI)=0.
            UDDX(GI)=0.    
            VDDY(GI)=0.
            WDDZ(GI)=0.
            UDDX2(GI)=0.    
            VDDY2(GI)=0.
            WDDZ2(GI)=0.
            
            UGDX(GI)=0.0
            VGDY(GI)=0.0
            WGDZ(GI)=0.0
             
            ABGI(GI)=0.0
             
            do ILOC=1,NLOC
               UD(GI)=UD(GI)+N(ILOC,GI)*DGNU((ELE-1)*NLOC+ILOC)
               VD(GI)=VD(GI)+N(ILOC,GI)*DGNV((ELE-1)*NLOC+ILOC)
               WD(GI)=WD(GI)+N(ILOC,GI)*DGNW((ELE-1)*NLOC+ILOC)
               UDDX(GI)=UDDX(GI) + NX(ILOC,GI)*DGNU((ELE-1)*NLOC+ILOC)
               VDDY(GI)=VDDY(GI) + NY(ILOC,GI)*DGNV((ELE-1)*NLOC+ILOC)
               WDDZ(GI)=WDDZ(GI) + NZ(ILOC,GI)*DGNW((ELE-1)*NLOC+ILOC)
            END DO
            IF(COGRAX) THEN
               DIRX(GI)=0.0
               DIRY(GI)=0.0
               DIRZ(GI)=1.0
            ELSE
               DIRX(GI)=XD(GI)
               DIRY(GI)=YD(GI)
               DIRZ(GI)=ZD(GI)
               RN=SQRT(DIRX(GI)**2+DIRY(GI)**2+DIRZ(GI)**2)
               DIRX(GI)=DIRX(GI)/RN
               DIRY(GI)=DIRY(GI)/RN
               DIRZ(GI)=DIRZ(GI)/RN
            ENDIF
            
               RSOXGRAGI=0.0
               RSOYGRAGI=0.0
               RSOZGRAGI=0.0
               
               GRAV_DEN_NOD_NX=0.0
               GRAV_DEN_NOD_NY=0.0
               GRAV_DEN_NOD_NZ=0.0

            do  L=1,NLOC
               IGLV=VONDGL((ELE-1)*NLOC+L)
        
               IF(DISOPN.EQ.0) THEN
                  UD(GI)=UD(GI) + N(L,GI)*(-UG(IGLV))
                  VD(GI)=VD(GI) + N(L,GI)*(-VG(IGLV))
                  WD(GI)=WD(GI) + N(L,GI)*(-WG(IGLV))
                  
                  UGDX(GI)=UGDX(GI) + NX(L,GI)*UG(IGLV)
                  VGDY(GI)=VGDY(GI) + NY(L,GI)*VG(IGLV)
                  WGDZ(GI)=WGDZ(GI) + NZ(L,GI)*WG(IGLV)
               ENDIF
           
               ABGI(GI)=ABGI(GI)+N(L,GI)*ABSORB(IGLV)
              IF(NBUOY.NE.0) THEN
               RSOXGRAGI=RSOXGRAGI+N(L,GI)*SOXGRA(IGLV)
               RSOYGRAGI=RSOYGRAGI+N(L,GI)*SOYGRA(IGLV)
               RSOZGRAGI=RSOZGRAGI+N(L,GI)*SOZGRA(IGLV)
               
               GRAV_DEN_NOD=SQRT(SOXGRA(IGLV)**2+SOYGRA(IGLV)**2+SOZGRA(IGLV)**2)
               
               GRAV_DEN_NOD_NX=GRAV_DEN_NOD_NX + NX(L,GI)*GRAV_DEN_NOD
               GRAV_DEN_NOD_NY=GRAV_DEN_NOD_NY + NY(L,GI)*GRAV_DEN_NOD
               GRAV_DEN_NOD_NZ=GRAV_DEN_NOD_NZ + NZ(L,GI)*GRAV_DEN_NOD
              ENDIF
            end do
            
              IF(NBUOY.NE.0) THEN
! Calculate the absorbtion BOY_ABSXX etc associated with vertical buoyancy...               
               GRA_SIG_CHAN= &
     &                    SIGN(1.0,DIRX(GI)*RSOXGRAGI+DIRY(GI)*RSOYGRAGI &
     &                     +DIRZ(GI)*RSOZGRAGI)   
               
               OVER_RELAX_BOY=1.0
               GRAV_DEN_NOD_NX=GRA_SIG_CHAN*GRAV_DEN_NOD_NX
               GRAV_DEN_NOD_NY=GRA_SIG_CHAN*GRAV_DEN_NOD_NY
               GRAV_DEN_NOD_NZ=GRA_SIG_CHAN*GRAV_DEN_NOD_NZ 
               RMID=DT* MAX(0.0, DIRX(GI)*GRAV_DEN_NOD_NX+DIRY(GI)*GRAV_DEN_NOD_NY &
     &                     +DIRZ(GI)*GRAV_DEN_NOD_NZ )*OVER_RELAX_BOY
              ELSE
               RMID=0.0
              ENDIF
           IF(NBUOY.NE.0) THEN
             IF(VERTICAL_BOY) THEN
! only add in a vertical bouyancy term.. which acts like absorption. 
               BOY_ABSXX(GI)=RMID*DIRX(GI)*DIRX(GI)
               BOY_ABSXY(GI)=RMID*DIRX(GI)*DIRY(GI)
               BOY_ABSXZ(GI)=RMID*DIRX(GI)*DIRZ(GI)
               
               BOY_ABSYX(GI)=RMID*DIRY(GI)*DIRX(GI)
               BOY_ABSYY(GI)=RMID*DIRY(GI)*DIRY(GI)
               BOY_ABSYZ(GI)=RMID*DIRY(GI)*DIRZ(GI)
               
               BOY_ABSZX(GI)=RMID*DIRZ(GI)*DIRX(GI)
               BOY_ABSZY(GI)=RMID*DIRZ(GI)*DIRY(GI)
               BOY_ABSZZ(GI)=RMID*DIRZ(GI)*DIRZ(GI)
             ELSE
               RMID=DT
               BOY_ABSXX(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NX
               BOY_ABSXY(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NY
               BOY_ABSXZ(GI)=RMID*DIRX(GI)*GRAV_DEN_NOD_NZ
               
               BOY_ABSYX(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NX
               BOY_ABSYY(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NY
               BOY_ABSYZ(GI)=RMID*DIRY(GI)*GRAV_DEN_NOD_NZ
               
               BOY_ABSZX(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NX
               BOY_ABSZY(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NY
               BOY_ABSZZ(GI)=RMID*DIRZ(GI)*GRAV_DEN_NOD_NZ
             ENDIF
           ELSE
               BOY_ABSXX(GI)=0.0
               BOY_ABSXY(GI)=0.0
               BOY_ABSXZ(GI)=0.0
               
               BOY_ABSYX(GI)=0.0
               BOY_ABSYY(GI)=0.0
               BOY_ABSYZ(GI)=0.0
               
               BOY_ABSZX(GI)=0.0
               BOY_ABSZY(GI)=0.0
               BOY_ABSZZ(GI)=0.0
           ENDIF

         end do second_gauss_loop
         
         MASSMM=0.0
         MASSLUMMM=0.0
         FAMAT    =0.0
         ILOCELE=CENTRMELE(ELE)-FINDRMELE(ELE)+1

! AMAT**********************
         do  ILOC=1,NLOC! Was loop 3501
            do  JLOC=1,NLOC! Was loop 3601

               VLM=0. 
      
               VLND  =0.0
               VLND2 =0.0
    
               VLMAB=0.0
               VLOMEGA=0.0

               do  GI=1,NGI
                  OMEGAGI=2.*FUNOME(XD(GI),YD(GI),ZD(GI))
     
                  VLN =-(NX(ILOC,GI)*UD(GI) &
     &                 +NY(ILOC,GI)*VD(GI) &
     &                 +NZ(ILOC,GI)*WD(GI) &
     &           -N(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*N(JLOC,GI)*DETWEI(GI) 
     
                  VLN2 =N(ILOC,GI)*(NX(jLOC,GI)*UD(GI) &
     &                 +NY(jLOC,GI)*VD(GI) &
     &                 +NZ(jLOC,GI)*WD(GI))*DETWEI(GI) 
             
                  VLND=VLND+VLN
                  VLND2=VLND2+VLN2
             
                  RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
! Coriolis terms lumped...
                  COMAT11(ELE,ILOC)=COMAT11(ELE,ILOC)+RNN*(0.0 - DIRX(GI)*DIRY(GI))*OMEGAGI
                  COMAT12(ELE,ILOC)=COMAT12(ELE,ILOC)+RNN*(-1.0+ DIRX(GI)*DIRX(GI))*OMEGAGI
                  COMAT21(ELE,ILOC)=COMAT21(ELE,ILOC)+RNN*(1.0 - DIRY(GI)*DIRY(GI))*OMEGAGI
                  COMAT22(ELE,ILOC)=COMAT22(ELE,ILOC)+RNN*(0.0 + DIRY(GI)*DIRX(GI))*OMEGAGI
                  COMAT31(ELE,ILOC)=COMAT31(ELE,ILOC)+RNN*(0.0 - DIRZ(GI)*DIRY(GI))*OMEGAGI
                  COMAT32(ELE,ILOC)=COMAT32(ELE,ILOC)+RNN*(0.0 + DIRZ(GI)*DIRX(GI))*OMEGAGI
          
           IF(NBUOY.NE.0) THEN
                 COMAT11(ELE,ILOC)=COMAT11(ELE,ILOC)+RNN*BOY_ABSXX(GI)
                 COMAT12(ELE,ILOC)=COMAT12(ELE,ILOC)+RNN*BOY_ABSXY(GI)
                 COMAT13(ELE,ILOC)=COMAT13(ELE,ILOC)+RNN*BOY_ABSXZ(GI)
                 
                 COMAT21(ELE,ILOC)=COMAT21(ELE,ILOC)+RNN*BOY_ABSYX(GI)
                 COMAT22(ELE,ILOC)=COMAT22(ELE,ILOC)+RNN*BOY_ABSYY(GI)
                 COMAT23(ELE,ILOC)=COMAT23(ELE,ILOC)+RNN*BOY_ABSYZ(GI)
                 
                 COMAT31(ELE,ILOC)=COMAT31(ELE,ILOC)+RNN*BOY_ABSZX(GI)
                 COMAT32(ELE,ILOC)=COMAT32(ELE,ILOC)+RNN*BOY_ABSZY(GI)
                 COMAT33(ELE,ILOC)=COMAT33(ELE,ILOC)+RNN*BOY_ABSZZ(GI)
                 
                 BOYMA11(ELE,ILOC)=BOYMA11(ELE,ILOC)+RNN*BOY_ABSXX(GI)
                 BOYMA12(ELE,ILOC)=BOYMA12(ELE,ILOC)+RNN*BOY_ABSXY(GI)
                 BOYMA13(ELE,ILOC)=BOYMA13(ELE,ILOC)+RNN*BOY_ABSXZ(GI)
                 
                 BOYMA21(ELE,ILOC)=BOYMA21(ELE,ILOC)+RNN*BOY_ABSYX(GI)
                 BOYMA22(ELE,ILOC)=BOYMA22(ELE,ILOC)+RNN*BOY_ABSYY(GI)
                 BOYMA23(ELE,ILOC)=BOYMA23(ELE,ILOC)+RNN*BOY_ABSYZ(GI)
                 
                 BOYMA31(ELE,ILOC)=BOYMA31(ELE,ILOC)+RNN*BOY_ABSZX(GI)
                 BOYMA32(ELE,ILOC)=BOYMA32(ELE,ILOC)+RNN*BOY_ABSZY(GI)
                 BOYMA33(ELE,ILOC)=BOYMA33(ELE,ILOC)+RNN*BOY_ABSZZ(GI)
           ENDIF
          
                  VLM=VLM+RNN 
                  VLOMEGA=VLOMEGA+RNN*OMEGAGI
             
                  VLMAB=VLMAB+ABGI(GI)*RNN 
                 
               end do
         
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
! conservative form BETA=1...
! non-conservative form BETA=0...
               VLKTEN=MATVLK(ELE,ILOC,JLOC)
               AMAT(ILOC,JLOC)=BETA*VLND      +(1.-BETA)*(VLND2) +VLMAB + VLKTEN
             
               MASSMM(ILOC,JLOC)=VLM
               MASSLUMMM(ILOC,ILOC)=MASSLUMMM(ILOC,ILOC)+VLM
             
            end do ! Was loop 3601
         end do ! Was loop 3501

! C1TDGELELOC**********************
         do  ILOC=1,MLOC
            do  JLOC=1,NLOC

               VLMXN=0.0
               VLMYN=0.0
               VLMZN=0.0
               do GI=1,NGI
                  VLMXN=VLMXN+MX(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)   
                  VLMYN=VLMYN+MY(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) 
                  VLMZN=VLMZN+MZ(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) 

               end do
               IF(IUNST_FREE.EQ.1) THEN
                 C1TDGELELOC(ILOC,JLOC)=-FTHETA*VLMXN
                 C2TDGELELOC(ILOC,JLOC)=-FTHETA*VLMYN
                 C3TDGELELOC(ILOC,JLOC)=-FTHETA*VLMZN
                 G1TDGELELOC(ILOC,JLOC)=0.0
                 G2TDGELELOC(ILOC,JLOC)=0.0
                 G3TDGELELOC(ILOC,JLOC)=0.0
               ELSE
                 C1TDGELELOC(ILOC,JLOC)=-VLMXN
                 C2TDGELELOC(ILOC,JLOC)=-VLMYN
                 C3TDGELELOC(ILOC,JLOC)=-VLMZN
               ENDIF
             
            end do
         end do
 
! calculate SNX,SNY,SNZ***
! Surface element of domain...
         do  IFACE=1,NFACE! Was loop 3344
         
! NB colele4 is ordered in terms of faces.
            ELE2=COLELE4((ELE-1)*5+IFACE)
            do COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
               JELE=COLMELE(COUNT)
               IF(JELE.EQ.ELE2) THEN
                  COUNTELE=COUNT
                  JLOCELE=COUNT-FINDRMELE(ELE)+1
               END IF
            END DO
         
            do SILOC=1,SNLOC
               SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
            END DO
            IF(MLOC.EQ.10) THEN
! Calculate P_MLOCSILOC2ILOC
               do SILOC=1,SNLOC
                  ILOC=SILOC2ILOC(SILOC)
                  do SJLOC=1,SNLOC
                     JLOC=SILOC2ILOC(SJLOC)
                     SKLOC=SILINK2(SILOC,SJLOC)
                     KLOC=ILINK2(ILOC,JLOC)
                     P_MLOCSILOC2ILOC(SKLOC)=KLOC
                  END DO
               END DO
            ELSE
               P_MLOCSILOC2ILOC=SILOC2ILOC
            ENDIF

! Put the volume elements togoether for element pairing...
            IF(ELE2.EQ.0) THEN
               SMEANXX=MEANXX(ELE)
               SMEANXY=MEANXY(ELE)
               SMEANXZ=MEANXZ(ELE)
               SMEANYY=MEANYY(ELE)
               SMEANYZ=MEANYZ(ELE)
               SMEANZZ=MEANZZ(ELE)
               ELE2_LOC_NODS=0
               ILOC_OTHER_SIDE=0
            ELSE
               SMEANXX=(MEANXX(ELE)*VOLELE(ELE)+MEANXX(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
               SMEANXY=(MEANXY(ELE)*VOLELE(ELE)+MEANXY(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
               SMEANXZ=(MEANXZ(ELE)*VOLELE(ELE)+MEANXZ(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
               SMEANYY=(MEANYY(ELE)*VOLELE(ELE)+MEANYY(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
               SMEANYZ=(MEANYZ(ELE)*VOLELE(ELE)+MEANYZ(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
               SMEANZZ=(MEANZZ(ELE)*VOLELE(ELE)+MEANZZ(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))

! Calculate the nodes on the other side of the face:
               do SILOC=1,SNLOC
                  ILOC=SILOC2ILOC(SILOC)
                  INOD=XONDGL((ELE-1)*NLOC+ILOC)
                  do ILOC2=1,NLOC
                     INOD2=XONDGL((ELE2-1)*NLOC+ILOC2)
                     IF(INOD2.EQ.INOD) ILOC_OTHER_SIDE(SILOC)=ILOC2
                  END DO
               END DO

! Calculate ELE2_LOC_NODS:
               ELE2_LOC_NODS=0
               do SILOC=1,SNLOC
                  ILOC=SILOC2ILOC(SILOC)
                  ILOC2=ILOC_OTHER_SIDE(SILOC)
                  ELE2_LOC_NODS(ILOC2)=ILOC
               END DO
               do ILOC2=1,NLOC
                  IF(ELE2_LOC_NODS(ILOC2).EQ.0) ELE2_LOC_NODS(ILOC2)=NLOC+1
               END DO
            ENDIF

            MATDIFX5=0.0
            MATDIFY5=0.0
            MATDIFZ5=0.0
            DMATDIFX5=0.0
            DMATDIFY5=0.0
            DMATDIFZ5=0.0
            MASSMM5=0.0
            SUFMATX=0.0
            SUFMATY=0.0
            SUFMATZ=0.0

            do SILOC=1,SNLOC
               ILOC=SILOC2ILOC(SILOC)
               INOD=XONDGL((ELE-1)*NLOC+ILOC)
               do SJLOC=SILOC,SNLOC,1
                  JLOC=SILOC2ILOC(SJLOC)
                  JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                  SKLOC=SILINK2(SILOC,SJLOC)
                  XSL(SKLOC)=0.5*(X(INOD)+X(JNOD))
                  YSL(SKLOC)=0.5*(Y(INOD)+Y(JNOD))
                  ZSL(SKLOC)=0.5*(Z(INOD)+Z(JNOD))
! Assume the mide side nodes are on the sphere
                  IF(ISPHERE.GE.2) THEN
                     IF(ILOC.NE.JLOC) THEN
                        RADI=SQRT(X(INOD)**2+Y(INOD)**2+Z(INOD)**2)
                        RADJ=SQRT(X(JNOD)**2+Y(JNOD)**2+Z(JNOD)**2)
                        RADP=0.5*(RADI+RADJ)
                        RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                        XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                        YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                        ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
                     ENDIF
                  ENDIF
               END DO
            END DO
       
! Form approximate surface normal (NORMX,NORMY,NORMZ)
            CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
     &                     X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)
! Recalculate the normal...
            CALL DGSDETNXLOC2(SNCLOC,SNGI, &
     &        XSL,YSL,ZSL, &
     &        SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
     &        SNORMXN,SNORMYN,SNORMZN, &
     &        NORMX,NORMY,NORMZ)

            SUD=0.0
            SVD=0.0
            SWD=0.0
            IN_SUD=0.0
            IN_SVD=0.0
            IN_SWD=0.0
            OUT_SUD=0.0
            OUT_SVD=0.0
            OUT_SWD=0.0
            do SILOC=1,SNLOC
               ILOC=SILOC2ILOC(SILOC)
               IGL =NDGLNO((ELE-1)*NLOC+ILOC)
               IDG=(ELE-1)*NLOC+ILOC
               IF(.NOT.SUF_CONSERV) THEN
                  IDG2=IDG
               ELSE 
                  IF(ELE2.NE.0) THEN

                     ILOC2=ILOC_OTHER_SIDE(SILOC)
                     IDG2=(ELE2-1)*NLOC+ILOC2
                  ELSE
                     IDG2=IDG
                  ENDIF
               ENDIF

               do  SGI=1,SNGI! Was loop 3313
                  SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*(0.5*(DGNU(IDG)+DGNU(IDG2))-UG(IGL))
                  SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*(0.5*(DGNV(IDG)+DGNV(IDG2))-VG(IGL))
                  SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*(0.5*(DGNW(IDG)+DGNW(IDG2))-WG(IGL))
                 
                  IN_SUD(SGI)=IN_SUD(SGI) + SN(SILOC,SGI)*(DGNU(IDG2)-UG(IGL))
                  IN_SVD(SGI)=IN_SVD(SGI) + SN(SILOC,SGI)*(DGNV(IDG2)-VG(IGL))
                  IN_SWD(SGI)=IN_SWD(SGI) + SN(SILOC,SGI)*(DGNW(IDG2)-WG(IGL))
                 
                  OUT_SUD(SGI)=OUT_SUD(SGI) + SN(SILOC,SGI)*(DGNU(IDG)-UG(IGL))
                  OUT_SVD(SGI)=OUT_SVD(SGI) + SN(SILOC,SGI)*(DGNV(IDG)-VG(IGL))
                  OUT_SWD(SGI)=OUT_SWD(SGI) + SN(SILOC,SGI)*(DGNW(IDG)-WG(IGL))
               END DO
            END DO
           
            IF(ADVECT_MEAN) THEN
! Original method...
               IN_SUD=SUD
               IN_SVD=SVD
               IN_SWD=SWD

               OUT_SUD=SUD
               OUT_SVD=SVD 
               OUT_SWD=SWD 
            ENDIF

            do SGI=1,SNGI
               NDOTQ(SGI)=SUD(SGI)*SNORMXN(SGI)+SVD(SGI)*SNORMYN(SGI) &
     &                 +SWD(SGI)*SNORMZN(SGI)
               IF(NDOTQ(SGI).LT.0.0) THEN
                  IN_NDOTQ(SGI)=IN_SUD(SGI)*SNORMXN(SGI)+IN_SVD(SGI)*SNORMYN(SGI) &
     &                      +IN_SWD(SGI)*SNORMZN(SGI)
                  OUT_NDOTQ(SGI)=0.0
               ELSE
                  IN_NDOTQ(SGI)=0.0
                  OUT_NDOTQ(SGI)=OUT_SUD(SGI)*SNORMXN(SGI)+OUT_SVD(SGI)*SNORMYN(SGI) &
     &                       +OUT_SWD(SGI)*SNORMZN(SGI)
               ENDIF
            END DO
           
            AVE_SNORMXN=0.0
            AVE_SNORMYN=0.0
            AVE_SNORMZN=0.0
            do SGI=1,SNGI
               AVE_SNORMXN=AVE_SNORMXN+SNORMXN(SGI)/SNGI
               AVE_SNORMYN=AVE_SNORMYN+SNORMYN(SGI)/SNGI
               AVE_SNORMZN=AVE_SNORMZN+SNORMZN(SGI)/SNGI
            END DO
! ************************
! Perform surface integration...
            ONE_OR_ZERO=0.0
            BNDSUF=.false.
            SUF_TOP=.FALSE.
            INLET_ADV=.FALSE.
            OUTLET_ADV=.FALSE.
            DIFF_BC=.FALSE.
          
            IF(ELE2.EQ.0) THEN
          
               IF(SUF_TOP_OCEAN) THEN
                  INLET_ADV=.FALSE.
                  OUTLET_ADV=.FALSE.
                  DIFF_BC=.FALSE.
               ELSE
                  INLET_ADV=.FALSE.
                  OUTLET_ADV=.FALSE.
                  DIFF_BC=.FALSE.
                  ICOUNTU=0
                  ICOUNTV=0
                  ICOUNTW=0
                  do SILOC=1,SNLOC
                     ILOC=SILOC2ILOC(SILOC)
                     IDG=(ELE-1)*NLOC+ILOC
                     ICOUNTU=ICOUNTU+DGMARKU(IDG)
                     ICOUNTV=ICOUNTV+DGMARKV(IDG)
                     ICOUNTW=ICOUNTW+DGMARKW(IDG)

                  END DO

                  ICOUNTU=ICOUNTU/SNLOC
                  ICOUNTV=ICOUNTV/SNLOC
                  ICOUNTW=ICOUNTW/SNLOC
                  ICOUNTUVW=ICOUNTU+ICOUNTV+ICOUNTW
! OUTLET_ADV=SUF_TOP=true if there is no velocity component defined normal 
! to the boundary...
                  IF(ICOUNTUVW.EQ.3) THEN
                     DIFF_BC=.TRUE.
                  ELSE
                     XDOM=((ABS(AVE_SNORMXN).GT.ABS(AVE_SNORMYN)).AND.(ABS(AVE_SNORMXN).GT.ABS(AVE_SNORMZN))) 
                     YDOM=((ABS(AVE_SNORMYN).GT.ABS(AVE_SNORMXN)).AND.(ABS(AVE_SNORMYN).GT.ABS(AVE_SNORMZN))) 
                     ZDOM=(.NOT.XDOM).AND.(.NOT.YDOM)
                     IF(XDOM.AND.(ICOUNTU.EQ.0)) SUF_TOP=.TRUE.
                     IF(YDOM.AND.(ICOUNTV.EQ.0)) SUF_TOP=.TRUE.
                     IF(ZDOM.AND.(ICOUNTW.EQ.0)) SUF_TOP=.TRUE.
                     IF(XDOM.AND.(ICOUNTU.EQ.1)) INLET_ADV=.TRUE.
                     IF(YDOM.AND.(ICOUNTV.EQ.1)) INLET_ADV=.TRUE.
                     IF(ZDOM.AND.(ICOUNTW.EQ.1)) INLET_ADV=.TRUE.
                     OUTLET_ADV=SUF_TOP

                  ENDIF

! endof IF(SUF_TOP_OCEAN) THEN else...
               ENDIF
          
               XC=0.0
               YC=0.0
               ZC=0.0
               do ILOC=1,NLOC
                  IXNOD=XONDGL((ELE-1)*NLOC+ILOC)
                  XC=XC+X(IXNOD)/NLOC
                  YC=YC+Y(IXNOD)/NLOC
                  ZC=ZC+Z(IXNOD)/NLOC
               END DO
               IF(INLET_ADV.OR.DIFF_BC) THEN
                  do SILOC=1,SMLOC
                     ILOC=P_MLOCSILOC2ILOC(SILOC)
                     PGLOBI=PNDGLN((ELE-1)*MLOC+ILOC)
                     do SJLOC=1,SNLOC
                        JLOC=SILOC2ILOC(SJLOC)
                        JDG=(ELE-1)*NLOC+JLOC
                        C1TCONU=0.0
                        C2TCONV=0.0
                        C3TCONW=0.0
                        do SGI=1,SNGI
                           RNN=SDETWE(SGI)*SM(SILOC,SGI)*SN(SJLOC,SGI)
                           C1TCONU=C1TCONU+SNORMXN(SGI)*RNN
                           C2TCONV=C2TCONV+SNORMYN(SGI)*RNN
                           C3TCONW=C3TCONW+SNORMZN(SGI)*RNN
                        END DO
                        IF(NDPSET.EQ.0) THEN
                           PRES_RHS(PGLOBI)=PRES_RHS(PGLOBI)&
     &            +C1TCONU*DIRI_BC_X(JDG)+C2TCONV*DIRI_BC_Y(JDG)+C3TCONW*DIRI_BC_Z(JDG)
                        ENDIF

                     END DO
                  END DO
               ENDIF
! BNDSUF=.true. if we are on a boundary condition surface i.e. inlet-outlet
! for the free surface.
               ICOUNT=0
               do SILOC=1,SNLOC
                  ILOC=SILOC2ILOC(SILOC)
                  GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                  ICOUNT=ICOUNT+BNDCON(GLOBI)
               ENDDO
               BNDSUF=((ICOUNT.EQ.SNLOC).or.OUTLET_ADV).and.(ele2.eq.0)
             
               IF(SUF_TOP_OCEAN.OR.(PREOPT.EQ.2).or.BNDSUF) THEN
! Only put surface term in if at top of ocean...
                  IF(COGRAX) THEN
                     LOCPNORMX=0.0
                     LOCPNORMY=0.0
                     LOCPNORMZ=1.0
                  ELSE
                     XC=0.0
                     YC=0.0
                     ZC=0.0
                     do ILOC=1,NLOC
                        IXNOD=XONDGL((ELE-1)*NLOC+ILOC)
                        XC=XC+X(IXNOD)/NLOC
                        YC=YC+Y(IXNOD)/NLOC
                        ZC=ZC+Z(IXNOD)/NLOC
                     END DO
                     RN=SQRT(XC**2+YC**2+ZC**2)
                     LOCPNORMX=XC/RN
                     LOCPNORMY=YC/RN
                     LOCPNORMZ=ZC/RN
                  ENDIF
                  IF(AVE_SNORMXN*LOCPNORMX&
                       &       + AVE_SNORMYN*LOCPNORMY&
                       &       + AVE_SNORMZN*LOCPNORMZ.GT.0.8) THEN
                     ONE_OR_ZERO=1.0
                     SUF_TOP=.TRUE.
                  ENDIF
                  IF((PREOPT.EQ.2).or.BNDSUF) ONE_OR_ZERO=1.0
                  IF(ONE_OR_ZERO.NE.0.0) THEN
            
                     do SILOC=1,SMLOC
                        ILOC=P_MLOCSILOC2ILOC(SILOC)

                        do SJLOC=1,SNLOC
                           JLOC=SILOC2ILOC(SJLOC)
! Have a surface integral on outlet boundary... 

                           C1TCONU=0.0
                           C2TCONV=0.0
                           C3TCONW=0.0
                           do SGI=1,SNGI
                              RNN=SDETWE(SGI)*SM(SILOC,SGI)*SN(SJLOC,SGI)

                              C1TCONU=C1TCONU+SNORMXN(SGI)*RNN
                              C2TCONV=C2TCONV+SNORMYN(SGI)*RNN
                              C3TCONW=C3TCONW+SNORMZN(SGI)*RNN
                           END DO
                           IF(IUNST_FREE.EQ.1) THEN
                             G1TDGELELOC(ILOC,JLOC)=G1TDGELELOC(ILOC,JLOC)+FTHETA*ONE_OR_ZERO*C1TCONU
                             G2TDGELELOC(ILOC,JLOC)=G2TDGELELOC(ILOC,JLOC)+FTHETA*ONE_OR_ZERO*C2TCONV
                             G3TDGELELOC(ILOC,JLOC)=G3TDGELELOC(ILOC,JLOC)+FTHETA*ONE_OR_ZERO*C3TCONW
                           ELSE
                             C1TDGELELOC(ILOC,JLOC)=C1TDGELELOC(ILOC,JLOC)+ONE_OR_ZERO*C1TCONU
                             C2TDGELELOC(ILOC,JLOC)=C2TDGELELOC(ILOC,JLOC)+ONE_OR_ZERO*C2TCONV
                             C3TDGELELOC(ILOC,JLOC)=C3TDGELELOC(ILOC,JLOC)+ONE_OR_ZERO*C3TCONW
                           ENDIF
                
                        END DO
                     END DO
                     
                   if(IUNST_FREE.eq.1) then
! Put in the surface mass matrix into DGCMC which is a compressive term.
                     IF(COGRAX) THEN
                       SDIRX=0.0
                       SDIRY=0.0
                       SDIRZ=1.0
                     ELSE
                       SDIRX=0.0
                       SDIRY=0.0
                       SDIRZ=0.0
                       DO SGI=1,SNGI
                         DO SILOC=1,SNCLOC
                           SDIRX(SGI)=SDIRX(SGI)+SNC(SILOC,SGI)*XSL(SILOC)
                           SDIRY(SGI)=SDIRY(SGI)+SNC(SILOC,SGI)*YSL(SILOC)
                           SDIRZ(SGI)=SDIRZ(SGI)+SNC(SILOC,SGI)*ZSL(SILOC)
                         END DO
                         RNN=SQRT(SDIRX(SGI)**2+SDIRY(SGI)**2+SDIRZ(SGI)**2)
                         SDIRX(SGI)=SDIRX(SGI)/RNN
                         SDIRY(SGI)=SDIRY(SGI)/RNN
                         SDIRZ(SGI)=SDIRZ(SGI)/RNN
                       END DO
                     ENDIF
                     do SILOC=1,SMLOC
                        ILOC=P_MLOCSILOC2ILOC(SILOC)
                        PGLOBI=PNDGLN((ELE-1)*MLOC+ILOC)
                        do SJLOC=1,SMLOC
                           JLOC=P_MLOCSILOC2ILOC(SJLOC)
                           PGLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
! Have a surface integral on outlet boundary... 
                           RNN=0.0
                           do SGI=1,SNGI
                              RNN=RNN+SDETWE(SGI)*SM(SILOC,SGI)*SM(SJLOC,SGI) &
           *(SNORMXN(SGI)*SDIRX(SGI) +SNORMYN(SGI)*SDIRY(SGI)  +SNORMZN(SGI)*SDIRZ(SGI)   )
                           END DO
                           CALL POSINMAT(POS,PGLOBI,PGLOBJ,FREDOP,FINCMC,COLCMC,NCMC)
                           DGCMC(POS)=DGCMC(POS)+RNN/(DT*DT*GRAVTY)
                        END DO
                     END DO
                   endif

                  ENDIF
               ENDIF
            ENDIF

            suf_top=suf_top.or.(PREOPT.EQ.2).or.BNDSUF   
            do  SILOC=1,SNLOC
               ILOC   =SILOC2ILOC(SILOC)
               IDG=(ele-1)*nloc+iloc
               do  SJLOC=1,SNLOC
                  JLOC =SILOC2ILOC(SJLOC)
                  JDG=(ELE-1)*NLOC+JLOC
                  JLOC2=ILOC_OTHER_SIDE(SJLOC)
                  VLNDOTQ_IN =0.0
                  VLNDOTQ_OUT=0.0
                  VLNDOTQ_ALL=0.0
                  VLM_NORX=0.0
                  VLM_NORY=0.0
                  VLM_NORZ=0.0
! Have a surface integral on outlet boundary...  
                  do  SGI=1,SNGI
                     RNN=SDETWE(SGI)*SN(SILOC,SGI)*SN(SJLOC,SGI)

                     VLNDOTQ_IN =VLNDOTQ_IN +IN_NDOTQ(SGI)*RNN
                     VLNDOTQ_OUT=VLNDOTQ_OUT+OUT_NDOTQ(SGI)*RNN
                     VLNDOTQ_ALL=VLNDOTQ_ALL+(IN_NDOTQ(SGI)+OUT_NDOTQ(SGI))*RNN
                     VLM_NORX=VLM_NORX+SNORMXN(SGI)*RNN
                     VLM_NORY=VLM_NORY+SNORMYN(SGI)*RNN
                     VLM_NORZ=VLM_NORZ+SNORMZN(SGI)*RNN
                  END DO
                  IF(SUF_TOP.and.(ELE2.EQ.0)) THEN 
! avoid applying any sort of b.c at the free surface
                     IF(BNDSUF.and.(.not.OUTLET_ADV)) THEN
                        FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_ALL*BETA
                     ELSE
                        FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_OUT*BETA
                     ENDIF
                  ELSE IF(ELE2.EQ.0) THEN
! Bottom bathymetry (ZERO BC)...
                     if(INLET_ADV.OR.DIFF_BC) then
                        FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_OUT*BETA
                        FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) -VLNDOTQ_IN*(1.0-BETA)
                        VECX(IDG)=VECX(IDG) -VLNDOTQ_in*DIRI_BC_X(JDG)
                        VECY(IDG)=VECY(IDG) -VLNDOTQ_in*DIRI_BC_Y(JDG)
                        VECZ(IDG)=VECZ(IDG) -VLNDOTQ_in*DIRI_BC_Z(JDG)
                     else IF(SUF_TOP_OCEAN) THEN
                        FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_OUT*BETA
                     ELSE

                        FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_OUT*BETA
                     endif
                     if(DIFF_BC) then
                        MATDIFX5(ILOC,JLOC)=MATDIFX5(ILOC,JLOC)-(SMEANXX*VLM_NORX+SMEANXY*VLM_NORY+SMEANXZ*VLM_NORZ)
                        MATDIFY5(ILOC,JLOC)=MATDIFY5(ILOC,JLOC)-(SMEANXY*VLM_NORX+SMEANYY*VLM_NORY+SMEANYZ*VLM_NORZ)
                        MATDIFZ5(ILOC,JLOC)=MATDIFZ5(ILOC,JLOC)-(SMEANXZ*VLM_NORX+SMEANYZ*VLM_NORY+SMEANZZ*VLM_NORZ)
                    
                        SUFMATX(ILOC,JLOC)=SUFMATX(ILOC,JLOC)+VLM_NORX
                        SUFMATY(ILOC,JLOC)=SUFMATY(ILOC,JLOC)+VLM_NORY
                        SUFMATZ(ILOC,JLOC)=SUFMATZ(ILOC,JLOC)+VLM_NORZ
                     endif
                  ELSE
                     assert(JLOCELE>=0)
                     FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_OUT*BETA&
     &-VLNDOTQ_IN*(1.-BETA)
                    FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC2)=FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC2)+VLNDOTQ_IN*BETA&
     &+VLNDOTQ_IN*(1.-BETA)
! add diffusion term...

                    MATDIFX5(ILOC,JLOC)=MATDIFX5(ILOC,JLOC)-(SMEANXX*VLM_NORX+SMEANXY*VLM_NORY+SMEANXZ*VLM_NORZ)
                    MATDIFY5(ILOC,JLOC)=MATDIFY5(ILOC,JLOC)-(SMEANXY*VLM_NORX+SMEANYY*VLM_NORY+SMEANYZ*VLM_NORZ)
                    MATDIFZ5(ILOC,JLOC)=MATDIFZ5(ILOC,JLOC)-(SMEANXZ*VLM_NORX+SMEANYZ*VLM_NORY+SMEANZZ*VLM_NORZ)
                    MATDIFX5(ILOC,JLOC2+NLOC)=MATDIFX5(ILOC,JLOC2+NLOC)+(SMEANXX*VLM_NORX+SMEANXY*VLM_NORY+SMEANXZ*VLM_NORZ)
                    MATDIFY5(ILOC,JLOC2+NLOC)=MATDIFY5(ILOC,JLOC2+NLOC)+(SMEANXY*VLM_NORX+SMEANYY*VLM_NORY+SMEANYZ*VLM_NORZ)
                    MATDIFZ5(ILOC,JLOC2+NLOC)=MATDIFZ5(ILOC,JLOC2+NLOC)+(SMEANXZ*VLM_NORX+SMEANYZ*VLM_NORY+SMEANZZ*VLM_NORZ)
                    SUFMATX(ILOC,JLOC)=SUFMATX(ILOC,JLOC)+VLM_NORX
                    SUFMATY(ILOC,JLOC)=SUFMATY(ILOC,JLOC)+VLM_NORY
                    SUFMATZ(ILOC,JLOC)=SUFMATZ(ILOC,JLOC)+VLM_NORZ
                 ENDIF
              END DO
           END DO

! Rotate if on sides or bottom of ocean only...
           IF(ROTATDG) THEN
              IF(SUF_TOP_OCEAN) THEN
                 IF((ELE2.EQ.0).AND.(.NOT.SUF_TOP).AND.(.NOT.BNDSUF)) THEN
! Place into rotations...
                    CALL PLAC_IN_SUF_ROT_DG(ELE,NLOC,TOTELE,SNLOC,SNGI,&
     &               SDETWE,SN, SNORMXN,SNORMYN,SNORMZN, SILOC2ILOC,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ,&
     &               DGDM1,DGDM2,DGDM3)
                 ENDIF
              ENDIF
           ENDIF

! Calculate the diffusion terms*****
           IF(DIFF_BC) THEN
              DMATDIFX5=MATDIFX5
              DMATDIFY5=MATDIFY5
              DMATDIFZ5=MATDIFZ5
           ENDIF
           MATDIFX5(1:NLOC,1:NLOC)=MATDIFX5(1:NLOC,1:NLOC)+MATDIFX(ELE,1:NLOC,1:NLOC)
           MATDIFY5(1:NLOC,1:NLOC)=MATDIFY5(1:NLOC,1:NLOC)+MATDIFY(ELE,1:NLOC,1:NLOC)
           MATDIFZ5(1:NLOC,1:NLOC)=MATDIFZ5(1:NLOC,1:NLOC)+MATDIFZ(ELE,1:NLOC,1:NLOC)
           MASSMM5(1:NLOC,1:NLOC) =MASSMM5(1:NLOC,1:NLOC) +MASSDGELE(ELE,1:NLOC,1:NLOC)

! Put contributions from ELE2 into MATDIFX5
           IF(ELE2.NE.0) THEN
              do ILOC2=1,NLOC
                 IILOC5=ELE2_LOC_NODS(ILOC2)
                 do JLOC2=1,NLOC
                    JJLOC5=ELE2_LOC_NODS(JLOC2)
                    MATDIFX5(IILOC5,JLOC2+NLOC)=MATDIFX5(IILOC5,JLOC2+NLOC)+MATDIFX(ELE2,ILOC2,JLOC2)
                    MATDIFY5(IILOC5,JLOC2+NLOC)=MATDIFY5(IILOC5,JLOC2+NLOC)+MATDIFY(ELE2,ILOC2,JLOC2)
                    MATDIFZ5(IILOC5,JLOC2+NLOC)=MATDIFZ5(IILOC5,JLOC2+NLOC)+MATDIFZ(ELE2,ILOC2,JLOC2)
                    MASSMM5(IILOC5,JJLOC5)=MASSMM5(IILOC5,JJLOC5)+MASSDGELE(ELE2,ILOC2,JLOC2)
                 END DO
              END DO
           ELSE
              MASSMM5(1+NLOC,1+NLOC) =1.0

           ENDIF
! calculate MASS2INV...
           CALL MATDMATINV(MASSMM5,MASSINVMM5,NLOC+1)
! M^{-1}*MATDIFX5
           CALL ABMATRIXMUL(MMATDIFX5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              MATDIFX5,  NLOC+1,NLOC+NLOC)
! M^{-1}*MATDIFY5
           CALL ABMATRIXMUL(MMATDIFY5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              MATDIFY5,  NLOC+1,NLOC+NLOC)
! M^{-1}*MATDIFZ5
           CALL ABMATRIXMUL(MMATDIFZ5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              MATDIFZ5,  NLOC+1,NLOC+NLOC)    

! SUFMATX*MMATDIFX5
           CALL ABMATRIXMUL(SMMATDIFX5,SUFMATX,   NLOC,NLOC+1,&
     &                               MMATDIFX5, NLOC+1,NLOC+NLOC)
! SUFMATX*MMATDIFX5
           CALL ABMATRIXMUL(SMMATDIFY5,SUFMATY,   NLOC,NLOC+1,&
     &                               MMATDIFY5, NLOC+1,NLOC+NLOC)
! SUFMATX*MMATDIFX5
           CALL ABMATRIXMUL(SMMATDIFZ5,SUFMATZ,   NLOC,NLOC+1,&
     &                               MMATDIFZ5, NLOC+1,NLOC+NLOC)
           SMMATDIF5=SMMATDIFX5+SMMATDIFY5+SMMATDIFZ5

           do SILOC=1,SNLOC
              ILOC=SILOC2ILOC(SILOC)
              do JLOC=1,NLOC
! add diffusion term (-VE because diffusion has - sign on lhs of eqns)...
                 FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC)   &
     &                                           -SMMATDIF5(ILOC,JLOC)
                 IF(ELE2.NE.0) THEN
                    assert(JLOCELE>=0)
                    FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC)&
     &                                             -SMMATDIF5(ILOC,NLOC+JLOC)
                 ENDIF
              END DO
           END DO
           
           IF(DIFF_BC) THEN
! M^{-1}*DMATDIFX5
              CALL ABMATRIXMUL(MMATDIFX5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              DMATDIFX5,  NLOC+1,NLOC+NLOC)
! M^{-1}*DMATDIFY5
              CALL ABMATRIXMUL(MMATDIFY5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              DMATDIFY5,  NLOC+1,NLOC+NLOC)
! M^{-1}*DMATDIFZ5
              CALL ABMATRIXMUL(MMATDIFZ5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              DMATDIFZ5,  NLOC+1,NLOC+NLOC)    

! SUFMATX*MMATDIFX5
              CALL ABMATRIXMUL(SMMATDIFX5,SUFMATX,   NLOC,NLOC+1,&
     &                               MMATDIFX5, NLOC+1,NLOC+NLOC)
! SUFMATX*MMATDIFX5
              CALL ABMATRIXMUL(SMMATDIFY5,SUFMATY,   NLOC,NLOC+1,&
     &                               MMATDIFY5, NLOC+1,NLOC+NLOC)
! SUFMATX*MMATDIFX5
              CALL ABMATRIXMUL(SMMATDIFZ5,SUFMATZ,   NLOC,NLOC+1,&
     &                               MMATDIFZ5, NLOC+1,NLOC+NLOC)
              SMMATDIF5=SMMATDIFX5+SMMATDIFY5+SMMATDIFZ5

              do SILOC=1,SNLOC
                 ILOC=SILOC2ILOC(SILOC)
                 IDG=(ELE-1)*NLOC+ILOC
                 do JLOC=1,NLOC
                    JDG=(ELE-1)*NLOC+JLOC
! add diffusion term (-VE because diffusion has - sign on lhs of eqns)...
                    VECX(IDG)=VECX(IDG)-SMMATDIF5(ILOC,JLOC)*DIRI_BC_X(JDG)
                    VECY(IDG)=VECY(IDG)-SMMATDIF5(ILOC,JLOC)*DIRI_BC_Y(JDG)
                    VECZ(IDG)=VECZ(IDG)-SMMATDIF5(ILOC,JLOC)*DIRI_BC_Z(JDG)
                 END DO
              END DO
         
           ENDIF
! Calculate the diffusion terms***** 

        end do ! Was loop 3344
! END OF Put surface integrals *********************

! AMAT:
        do ILOC=1,NLOC   
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           IDG=(ELE-1)*NLOC+ILOC 
           do JLOC=1,NLOC
              JDG=(ELE-1)*NLOC+JLOC 
              A11MAT=MASSMM(ILOC,JLOC)/DT -(1.-THETA)*AMAT(ILOC,JLOC)
              VECX(IDG)=VECX(IDG) +A11MAT*DGU(JDG) +MASSMM(ILOC,JLOC)*TOTSOUX(JLOC)
              VECY(IDG)=VECY(IDG) +A11MAT*DGV(JDG) +MASSMM(ILOC,JLOC)*TOTSOUY(JLOC)
              VECZ(IDG)=VECZ(IDG) +A11MAT*DGW(JDG) +MASSMM(ILOC,JLOC)*TOTSOUZ(JLOC)
! Lump...
              MLELE(IDG)=MLELE(IDG) + MASSMM(ILOC,JLOC)
              ML(GLOBI) =ML(GLOBI)  + MASSMM(ILOC,JLOC)
           END DO
           VECX(IDG)=VECX(IDG) +COMAT11(ELE,ILOC)*DGUNEW(IDG)+COMAT12(ELE,ILOC)*DGVNEW(IDG) &
     &                         +COMAT13(ELE,ILOC)*DGWNEW(IDG)
           VECY(IDG)=VECY(IDG) +COMAT21(ELE,ILOC)*DGUNEW(IDG)+COMAT22(ELE,ILOC)*DGVNEW(IDG) &
     &                         +COMAT23(ELE,ILOC)*DGWNEW(IDG)
           VECZ(IDG)=VECZ(IDG) +COMAT31(ELE,ILOC)*DGUNEW(IDG)+COMAT32(ELE,ILOC)*DGVNEW(IDG) &
     &                         +COMAT33(ELE,ILOC)*DGWNEW(IDG)
        END DO
        
! now the surface integral for advection...       
        assert(JLOCELE>=0)
        do COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
           ELE2=COLMELE(COUNT)
           JLOCELE=COUNT-FINDRMELE(ELE)+1
           do ILOC=1,NLOC   
              IDG=(ELE-1)*NLOC+ILOC 
              do JLOC=1,NLOC
                 JDG=(ELE2-1)*NLOC+JLOC 
                 A11MAT= -(1.-THETA)*FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC)
                 VECX(IDG)=VECX(IDG) +A11MAT*DGU(JDG) 
                 VECY(IDG)=VECY(IDG) +A11MAT*DGV(JDG) 
                 VECZ(IDG)=VECZ(IDG) +A11MAT*DGW(JDG) 
! ************PUT INTO MATRIX...
                 THETA_FAMAT=THETA*FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC)
                 IF(SMALL_DGMAT_MEM) THEN
                    IJDISP=(ILOC-1)*NLOC + JLOC
                    BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)=BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)&
     &                                              +THETA_FAMAT
                 ELSE
                    do ID=1,3
                       BIG_ILOC=(ID-1)*NLOC+ILOC
                       JD=ID
                       BIG_JLOC=(JD-1)*NLOC+JLOC
                       IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
                       BIGMELE((COUNT-1)*NBLOCK+IJDISP)=BIGMELE((COUNT-1)*NBLOCK+IJDISP)&
     &                                             +THETA_FAMAT
                    END DO
                 ENDIF
! ************PUT INTO MATRIX...
              END DO
           END DO
        END DO

        do ILOC=1,NLOC  
           IDG=(ELE-1)*NLOC+ILOC 
! ******add in the pre-discretised source*****
           VECX(IDG)=VECX(IDG) + VECX_SGSADD(IDG)
           VECY(IDG)=VECY(IDG) + VECY_SGSADD(IDG)
           VECZ(IDG)=VECZ(IDG) + VECZ_SGSADD(IDG)
           IF(.NOT.CV_PRES) THEN
              do JLOC=1,MLOC
! Put pressure terms in...  
                 PGLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
                 IF(IUNST_FREE.EQ.1) THEN
                   RPRES=P(PGLOBJ)+((1.-FTHETA)/FTHETA)*POLD(PGLOBJ)
                   VECX(IDG)=VECX(IDG) +C1TDGELELOC(JLOC,ILOC)*RPRES
                   VECY(IDG)=VECY(IDG) +C2TDGELELOC(JLOC,ILOC)*RPRES
                   VECZ(IDG)=VECZ(IDG) +C3TDGELELOC(JLOC,ILOC)*RPRES
                 ELSE
                   VECX(IDG)=VECX(IDG) +C1TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                   VECY(IDG)=VECY(IDG) +C2TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                   VECZ(IDG)=VECZ(IDG) +C3TDGELELOC(JLOC,ILOC)*P(PGLOBJ)
                 ENDIF
              END DO
           END IF
        END DO

! Amend AMAT,BMAT,CMAT,DMAT to take into account time dep.
        AMAT=(1./DT)*MASSMM+THETA*AMAT
        
        IF(.NOT.SMALL_DGMAT_MEM) THEN
           BIG_AMAT=0.0
           BIG_AMAT(1:NLOC,1:NLOC)                  =AMAT(1:NLOC,1:NLOC)
           BIG_AMAT(1+NLOC:2*NLOC,1+NLOC:2*NLOC)    =AMAT(1:NLOC,1:NLOC)
           BIG_AMAT(1+2*NLOC:3*NLOC,1+2*NLOC:3*NLOC)=AMAT(1:NLOC,1:NLOC)
! This is implicit Coriolis lumped...
           do ILOC=1,NLOC
            BIG_AMAT(ILOC,       ILOC)          =BIG_AMAT(ILOC,         ILOC)      + COMAT11(ELE,ILOC)
            BIG_AMAT(ILOC+NLOC,  ILOC+NLOC)     =BIG_AMAT(ILOC+NLOC,ILOC+NLOC)     + COMAT22(ELE,ILOC)
            BIG_AMAT(ILOC+2*NLOC,  ILOC+2*NLOC) =BIG_AMAT(ILOC+2*NLOC,ILOC+2*NLOC) + COMAT33(ELE,ILOC)
            BIG_AMAT(ILOC,       ILOC+NLOC)   =COMAT12(ELE,ILOC)
            BIG_AMAT(ILOC,       ILOC+2*NLOC) =COMAT13(ELE,ILOC)
            BIG_AMAT(ILOC+NLOC,  ILOC)        =COMAT21(ELE,ILOC)
            BIG_AMAT(ILOC+NLOC,  ILOC+2*NLOC) =COMAT23(ELE,ILOC)
            BIG_AMAT(ILOC+2*NLOC,ILOC)      =COMAT31(ELE,ILOC)
            BIG_AMAT(ILOC+2*NLOC,ILOC+NLOC) =COMAT32(ELE,ILOC)
           END DO
        ENDIF

! Now make bigger matrix for coupled vels
          
! Just this element...    
        COUNT=CENTRMELE(ELE)
        do ILOC=1,NLOC
           do JLOC=1,NLOC
! Put into matrix...
! Find count. ***************************************
              IF(SMALL_DGMAT_MEM) THEN
                 IJDISP=(ILOC-1)*NLOC+JLOC
                 BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)=BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)&
     &             +AMAT(ILOC,JLOC)
              ELSE
                 do ID=1,3
                    BIG_ILOC=(ID-1)*NLOC+ILOC
                    do JD=1,3
                       BIG_JLOC=(JD-1)*NLOC+JLOC
                       IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
                       BIGMELE((COUNT-1)*NBLOCK+IJDISP)=BIGMELE((COUNT-1)*NBLOCK+IJDISP)&
     &             +BIG_AMAT(BIG_ILOC,BIG_JLOC)
                    END DO
                 END DO
              ENDIF
           END DO
        END DO

! Store pressure matrix...     
      IF(IUNST_FREE.EQ.1) THEN   
        G1TDGELE(ELE,1:MLOC,1:NLOC)=G1TDGELELOC(1:MLOC,1:NLOC)
        G2TDGELE(ELE,1:MLOC,1:NLOC)=G2TDGELELOC(1:MLOC,1:NLOC)
        G3TDGELE(ELE,1:MLOC,1:NLOC)=G3TDGELELOC(1:MLOC,1:NLOC)
      ENDIF
        C1TDGELE(ELE,1:MLOC,1:NLOC)=C1TDGELELOC(1:MLOC,1:NLOC)
        C2TDGELE(ELE,1:MLOC,1:NLOC)=C2TDGELELOC(1:MLOC,1:NLOC)
        C3TDGELE(ELE,1:MLOC,1:NLOC)=C3TDGELELOC(1:MLOC,1:NLOC)
         
     end do second_element_loop

! a CV pressure treatment...
     IF(CV_PRES) THEN
        CALL CV_CTY_SGS( NONODS, XNONOD, TOTELE, FREDOP, &
     &                 NDGLNO, XONDGL, C1TDGELE,C2TDGELE,C3TDGELE,&
     &                        X, Y, Z, .TRUE., .FALSE., ISPHERE, &
     &                        NLOC, MLOC, SNLOC,  &
     &                        COLELE4,LOCLIST,NFACE,BNDCON,COGRAX, &
     &                        SUF_TOP_OCEAN )
        do ELE=1,TOTELE
           do ILOC=1,NLOC  
              IDG=(ELE-1)*NLOC+ILOC 
              do JLOC=1,MLOC
! Put pressure terms in...  
                 PGLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
                 VECX(IDG)=VECX(IDG) +C1TDGELE(ELE,JLOC,ILOC)*P(PGLOBJ)
                 VECY(IDG)=VECY(IDG) +C2TDGELE(ELE,JLOC,ILOC)*P(PGLOBJ)
                 VECZ(IDG)=VECZ(IDG) +C3TDGELE(ELE,JLOC,ILOC)*P(PGLOBJ)
              END DO
           END DO
        END DO
           
        ALLOCATE(TESTX(NONODS))
        ALLOCATE(TESTY(NONODS))
        ALLOCATE(TESTZ(NONODS))
        TESTX=0.0
        TESTY=0.0
        TESTZ=0.0
        do ELE=1,TOTELE
           do ILOC=1,NLOC  
              IDG=NDGLNO((ELE-1)*NLOC+ILOC)
              do JLOC=1,MLOC
! Put pressure terms in...  
                 PGLOBJ=PNDGLN((ELE-1)*MLOC+JLOC)
                 TESTX(IDG)=TESTX(IDG) +C1TDGELE(ELE,JLOC,ILOC)*z(PGLOBJ)/ML(IDG)
                 TESTY(IDG)=TESTY(IDG) +C2TDGELE(ELE,JLOC,ILOC)*z(PGLOBJ)/ML(IDG)
                 TESTZ(IDG)=TESTZ(IDG) +C3TDGELE(ELE,JLOC,ILOC)*z(PGLOBJ)/ML(IDG)

                 if(idg.eq.1) then
                    ewrite(2,*) 'ele,iloc,jloc,c1t,z:',ele,iloc,jloc,C1TDGELE(ELE,JLOC,ILOC),z(PGLOBJ)
                 endif
              END DO
           END DO
        END DO

        CALL PMINMX(TESTX,NONODS,'TESTX*******  ')
        CALL PMINMX(TESTY,NONODS,'TESTY*******  ')
        CALL PMINMX(TESTZ,NONODS,'TESTZ*******  ')
        STOP 3831
           
     ENDIF
         
     ewrite(2,*) 'r2norm(vecx,NONODS,0):',r2norm(vecx,NONODS,0)
     ewrite(2,*) 'r2norm(vecy,NONODS,0):',r2norm(vecy,NONODS,0)
     ewrite(2,*) 'r2norm(vecz,NONODS,0):',r2norm(vecz,NONODS,0)
         
! now the surface integral for advection...  
     do ele=1,-10     
        do COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
           ELE2=COLMELE(COUNT)

           JLOCELE=COUNT-FINDRMELE(ELE)+1
           do ILOC=1,NLOC   
              GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
              IDG=(ELE-1)*NLOC+ILOC 
              do JLOC=1,NLOC
                 JDG=(ELE2-1)*NLOC+JLOC 
! ************PUT INTO MATRIX...
                 do ID=1,3
                    do jD=1,3
                       BIG_ILOC=(ID-1)*NLOC+ILOC
                       BIG_JLOC=(JD-1)*NLOC+JLOC
                       IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
                    END DO
                 END DO
! ************PUT INTO MATRIX...
              END DO
           END DO
        END DO
     end do

     CALL PMINMX(MUPTXX,NONODS,'MUPTXX*******  ')
     CALL PMINMX(MUPTXY,NONODS,'MUPTXY*******  ')
     CALL PMINMX(MUPTXZ,NONODS,'MUPTXZ*******  ')
     CALL PMINMX(MUPTYY,NONODS,'MUPTYY*******  ') 
     CALL PMINMX(MUPTYZ,NONODS,'MUPTYZ*******  ')
     CALL PMINMX(MUPTZZ,NONODS,'MUPTZZ*******  ')
     CALL PMINMX(absorb,NONODS,'absorb*******  ')
     CALL PMINMX(PRES_RHS,fredop,'PRES_RHS*******  ')

     CALL PMINMX(SOURCX,nonods,'SOURCX*******  ')
     CALL PMINMX(SOURCy,nonods,'SOURCy*******  ')
     CALL PMINMX(SOURCz,nonods,'SOURCz*******  ')
         
     CALL PMINMX(VECX_SGSADD,nloc*totele,'VECX_SGSADD*******  ')
     CALL PMINMX(VECy_SGSADD,nloc*totele,'VECy_SGSADD*******  ')
     CALL PMINMX(VECz_SGSADD,nloc*totele,'VECz_SGSADD*******  ')
         
   END SUBROUTINE DIFF3DSGS_DG

   SUBROUTINE GET_QUAD_LOC(ELE,TQUAD,TOLDQUAD,NQLOC,T,TOLD, &
     &            DIFXT,DIFYT,DIFZT, DIFXTOLD,DIFYTOLD,DIFZTOLD, &
     &            NONODS,XNONOD,NLOC,TOTELE,NDGLNO,XNDGLN,X,Y,Z)
! This sub calculates the quadratic variation inside an element
! defined by a quadratic tet. but recovered from a linear tet 
! nodal values and the nodal gradients.
     INTEGER ELE,NONODS,XNONOD,NLOC,TOTELE,NQLOC
     INTEGER NDGLNO(NLOC*TOTELE),XNDGLN(NLOC*TOTELE)
     REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
     REAL TQUAD(NQLOC),TOLDQUAD(NQLOC)
     REAL T(NONODS),TOLD(NONODS)
     REAL DIFXT(NONODS),DIFYT(NONODS),DIFZT(NONODS)
     REAL DIFXTOLD(NONODS),DIFYTOLD(NONODS),DIFZTOLD(NONODS)
! Local variables...
     INTEGER ILOC,JLOC,KLOC,INOD,JNOD,INODX,JNODX
     REAL DIRX,DIRY,DIRZ,DIST,STAN_GRA_I,STAN_GRA_J
! Variables T:
     REAL INGRADT,JNGRADT,GRAD1ST
     REAL ALT_INGRADT,ALT_JNGRADT
     REAL MAG_T_I,MAG_T_J,MAG_T
! Variables TOLD:
     REAL INGRADTOLD,JNGRADTOLD,GRAD1STOLD
     REAL ALT_INGRADTOLD,ALT_JNGRADTOLD
     REAL MAG_TOLD_I,MAG_TOLD_J,MAG_TOLD
         
     do ILOC=1,NLOC
        INOD=NDGLNO((ELE-1)*NLOC+ILOC)
        KLOC=ILINK2(ILOC,ILOC)
        TQUAD(KLOC)=T(INOD)
        TOLDQUAD(KLOC)=TOLD(INOD)
     END DO

     do ILOC=1,NLOC
        INOD =NDGLNO((ELE-1)*NLOC+ILOC)
        INODX=XNDGLN((ELE-1)*NLOC+ILOC)
        do JLOC=ILOC+1,NLOC
           JNOD =NDGLNO((ELE-1)*NLOC+JLOC)
           JNODX=XNDGLN((ELE-1)*NLOC+JLOC)
           KLOC=ILINK2(ILOC,JLOC)
           DIRX=X(JNODX)-X(INODX)
           DIRY=Y(JNODX)-Y(INODX)
           DIRZ=Z(JNODX)-Z(INODX)
           DIST=SQRT(DIRX**2+DIRY**2+DIRZ**2)
           DIRX=DIRX/DIST
           DIRY=DIRY/DIST
           DIRZ=DIRZ/DIST
           STAN_GRA_I=2.0/DIST
           STAN_GRA_J=-2.0/DIST
! FOR T:
           INGRADT=DIRX*DIFXT(INOD)+DIRY*DIFYT(INOD)+DIRZ*DIFZT(INOD)
           JNGRADT=DIRX*DIFXT(JNOD)+DIRY*DIFYT(JNOD)+DIRZ*DIFZT(JNOD)
           GRAD1ST=(T(JNOD)-T(INOD))/DIST
           ALT_INGRADT=INGRADT-GRAD1ST
           ALT_JNGRADT=JNGRADT-GRAD1ST
           MAG_T_I=ALT_INGRADT/STAN_GRA_I
           MAG_T_J=ALT_JNGRADT/STAN_GRA_J
           MAG_T=0.5*(MAG_T_I+MAG_T_J)
           TQUAD(KLOC)=0.5*(T(JNOD)+T(INOD)) + MAG_T
! FOR TOLD:
           INGRADTOLD=DIRX*DIFXTOLD(INOD)+DIRY*DIFYTOLD(INOD)+DIRZ*DIFZTOLD(INOD)
           JNGRADTOLD=DIRX*DIFXTOLD(JNOD)+DIRY*DIFYTOLD(JNOD)+DIRZ*DIFZTOLD(JNOD)
           GRAD1STOLD=(TOLD(JNOD)-TOLD(INOD))/DIST
           ALT_INGRADTOLD=INGRADTOLD-GRAD1STOLD
           ALT_JNGRADTOLD=JNGRADTOLD-GRAD1STOLD
           MAG_TOLD_I=ALT_INGRADTOLD/STAN_GRA_I
           MAG_TOLD_J=ALT_JNGRADTOLD/STAN_GRA_J
           MAG_TOLD=0.5*(MAG_TOLD_I+MAG_TOLD_J)
           TOLDQUAD(KLOC)=0.5*(TOLD(JNOD)+TOLD(INOD)) + MAG_TOLD
        END DO
     END DO

   END SUBROUTINE GET_QUAD_LOC

   SUBROUTINE CV_CTY_SGS( NONODS, XNONOD, TOTELE, FREDOP, &
     &                 NDGLNO, XONDGL, C1TDGELE,C2TDGELE,C3TDGELE,&
     &                        X, Y, Z, D3, DCYL, ISPHERE, &
     &                        NLOC, MLOC, SNLOC,  &
     &                        COLELE4,LOCLIST,NFACE,BNDCON,COGRAX, &
     &                        SUF_TOP_OCEAN )
! =============================================================
! Subroutine to construct the matrices C1TDGELE,C2TDGELE,C3TDGELE
! =============================================================
    
! =============================================================
! inputs 
! =============================================================
     INTEGER, INTENT(IN) :: NONODS, XNONOD, TOTELE
! nonods = total number of nodes
! xnonod = length of co-ordinate vectors
! totele = total number of elements
    
     INTEGER, INTENT(IN) :: FREDOP
! fredop = number of pressure degrees of freedom
      
     INTEGER, INTENT(IN) :: NLOC, MLOC, SNLOC
! nloc  = number of displacement nodes per element
! mloc  = number of pressure nodes per element
! ngi   = number of gauss points (for current element)
! snloc = number of surface nodes per element
      
     LOGICAL, INTENT(IN) :: D3, DCYL
! d3   = flag for three dimensional simulation
! dcyl = flag for cylindircal simulation
      
     INTEGER, INTENT(IN) :: NDGLNO(TOTELE*NLOC)
! ndglno(totele*nloc) = element pointers for unknowns
      
     INTEGER, INTENT(IN) :: XONDGL(TOTELE*NLOC)
! XONDGL(totele*nloc) = element pointer for co-ordinates
      
     REAL, INTENT(IN) :: X(XNONOD), Y(XNONOD), Z(XNONOD)
! x(nonod) = x co-ordinates of the nodes
! y(nonod) = y co-ordinates of the nodes
! z(nonod) = z co-ordinates of the nodes

     INTEGER :: ISPHERE
      
     INTEGER NFACE,COLELE4(5*TOTELE),BNDCON(NONODS)
     INTEGER LOCLIST(NFACE,SNLOC)
     LOGICAL COGRAX,SUF_TOP_OCEAN
    
! =============================================================
! outputs
! =============================================================
     REAL C1TDGELE(TOTELE,MLOC,NLOC)
     REAL C2TDGELE(TOTELE,MLOC,NLOC),C3TDGELE(TOTELE,MLOC,NLOC)
      
! =============================================================
! local variables
! =============================================================
      
     INTEGER :: OPTELM, ELETYP
      ! optelm = defines (to be set) element type based on nloc
      !           optelm = 1 -> linear        triangles
      !           optelm = 2 -> quadratic     triangles
      !           optelm = 3 -> bi-linear     quadrilaterals
      !           optelm = 4 -> bi-quadratic  quadrilaterals
      !           optelm = 5 -> linear        tetrahedra
      !           optelm = 6 -> quadratic     tetrahedra
      !           optelm = 7 -> tri-linear    hexahedra
      !           optelm = 8 -> tri-quadratic hexahedra
      ! eletyp = defines element type as output from setelm
      
     INTEGER :: CVNGI, CVNLOC, SELETYP
      ! cvngi    = number of body gauss points in control volume
      ! cvnloc   = number of control volume nodes (should equal nloc)
      ! seletyp  = surface element type
      
     INTEGER :: SNGI, SNLOC2, SVNGI
      ! sngi  = number of surface gauss points
      ! snloc2  = number of surface gauss points
      ! svngi  = number of control volume face gauss points
      
     LOGICAL :: CVD3, CVDCYL, REDQUAD
      ! cvd3 = control volume three dimensional flag
      ! cvdcyl = control volume cylindrical flag
      ! redquad = flag indicating whether reduced quadrature is to be used
      !           redquad = .true.  -> use reduced quadrature
      !           redquad = .false. -> no reduced quadrature
     PARAMETER( REDQUAD=.FALSE.) 
      
     INTEGER, ALLOCATABLE :: NEILOC(:,:), SNEILOC(:,:)
      ! neiloc(nloc,svngi)     = neighbour node for a given node/gauss-point pair 
      ! sneiloc(snloc,sngi)    = control volume flag for a given surface node/gauss-point pair 
      
     REAL, ALLOCATABLE :: SVN(:,:), SVNLX(:,:), SVNLY(:,:), SVWEIGH(:)
      ! svn(nloc,svngi)       = the shape function evaluated 
      !                         for each node at each surface gauss point
      ! svnl[x/y](nloc,svngi) = the surface derivatives of the shape 
      !                         function for each node at those same points, and
      ! svweigh(svngi)        = the Gauss weights to use when integrating around 
      !                         the control volume surface
      
     REAL, ALLOCATABLE :: SN(:,:), SNLX(:,:), SNLY(:,:), SWEIGH(:)
      ! sn(snloc,sngi)       = the surface shape function evaluated 
      !                        for each node at each surface gauss point
      ! snl[x/y](snloc,sngi) = the surface derivatives of the shape 
      !                         function for each node at those same points, and
      ! sweigh(sngi)         = the Gauss weights to use when integrating around 
      !                        the element surface
      
     REAL, ALLOCATABLE :: CVM(:,:)
      ! cvm(nloc,svngi) = control volume centred shape function
      !                   (evaluated at each body gauss point around a node)
    
     INTEGER :: NCOLGPTS
      ! ncolgpts = length of colgpts(:) - to be determined
      
     INTEGER, ALLOCATABLE :: COLGPTS(:), FINDGPTS(:), TEMPCOLGPTS(:)
      ! colgpts(ncolgpts)       = returns gauss point pointer from index of findgpts
      ! findgpts(nloc+1)        = gives index for gauss points
      ! tempcolgpts(nloc*svngi) = temporary storage while size of colgpts(:) is determined

     REAL, ALLOCATABLE :: CVDETWEI(:)
     REAL, ALLOCATABLE :: CVNORMX(:), CVNORMY(:), CVNORMZ(:)
      ! cvdetwei(svngi) = determinant multiplied by the weight at the 
      !                   control volume face gauss points
      ! cvnormx(svngi)  = normal to control volume face gauss points x-component
      ! cvnormy(svngi)  = normal to control volume face gauss points x-component
      ! cvnormz(svngi)  = normal to control volume face gauss points x-component
      
     REAL :: SNORMX, SNORMY, SNORMZ, SAREA
      ! snorm[x/y/z] = approximation of normal over entire surface element face
      ! sarea        = area of surface element face
      
     REAL, ALLOCATABLE :: SDETWEI(:)
     REAL, ALLOCATABLE :: SNORMXN(:), SNORMYN(:), SNORMZN(:)
      ! sdetwei(sngi) = determinant times weight of surface element at gauss points
      ! snorm[x/y/z]n(sngi) = precise components of normal to surface element face at each gauss point
        
     INTEGER SILOC2ILOC(SNLOC),TEMPSNDGLN(SNLOC)
      
     LOGICAL :: GOTPNTS
      ! gotpnts = .true.  -> use pntx/y/z to determine direction of the face normal
      ! gotpnts = .false. -> use the current node to determine the direction of the face normal
     PARAMETER( GOTPNTS = .FALSE. )
      
     REAL :: PNTX, PNTY, PNTZ ! dummies
    !   PARAMETER( PNTX = 0., PNTY = 0., PNTZ = 0.)
      
      ! integer counters and pointers:
     INTEGER :: ELE, COUNT, GCOUNT, GI, JCOUNT
     INTEGER :: SELE, SGI, SILOC, SJLOC
     INTEGER :: ILOC, INOD
     INTEGER :: JLOC, JNOD
     INTEGER :: OLOC, ONOD
     INTEGER :: KLOC, KGLS, KGLV
     INTEGER :: LOWER, UPPER, POSMAT
     INTEGER :: IPOS, JPOS, MAT, NOD, ICOUNT, IFACE, IXNOD
     INTEGER ELE2,GLOBI
      
      ! end of integer counters
      
     REAL :: DSNDET
      ! dsndet   = surface shape function n times sdetwei(sgi) time nvddengi times volfrac
      
     REAL :: SVNDET, DSVNDET, SNDET
      ! svndet  = control volume n shape function times detwei(gi)
      ! dsvndet = control volume n shape function times detwei(gi) times nvddengi times volfrac
      ! sndet   = surface n shape function times sdetwe(sgi)
      
     REAL AVE_SNORMXN,AVE_SNORMYN,AVE_SNORMZN
     REAL LOCPNORMX,LOCPNORMY,LOCPNORMZ,XC,YC,ZC,RN
     LOGICAL BNDSUF,SUF_TOP
      
! =============================================================
! 
! Description                                   Programmer      Date
! ==================================================================
! Original version  .............................. CP      21/03/08
! =============================================================
       
     C1TDGELE=0.0
     C2TDGELE=0.0
     IF(D3) C3TDGELE=0.0

! Define the element type
     IF(D3) THEN                     
        IF(NLOC.EQ.4)  OPTELM=5       
        IF(NLOC.EQ.8)  OPTELM=7       
        IF(NLOC.EQ.10) OPTELM=6       
        IF(NLOC.EQ.27) OPTELM=8       
     ELSE                            
        IF(NLOC.EQ.3)  OPTELM=1       
        IF(NLOC.EQ.4)  OPTELM=3       
        IF(NLOC.EQ.6)  OPTELM=2       
        IF(NLOC.EQ.9)  OPTELM=4       
     ENDIF
            
! Use the SETELM routine to define the number of gauss points 
! on the control volume surface (SVNGI), based on the element type
! SETELM takes the inputs OPTELM and REDQUAD, all other
! terms are outputs and are not used
     CALL SETELM( ELETYP, CVNGI,   CVNLOC, &
     &               OPTELM, SELETYP, &
     &               SNGI,   SNLOC2,  SVNGI, &
     &               CVD3,   CVDCYL,  REDQUAD ) 
      
! make sure we are using the correct shape function...
     IF(CVNLOC.NE.NLOC) THEN
        FLAbort('CVNLOC.NE.NLOC IN CV_CTY_SGS')
     ENDIF
     IF(SNLOC2.NE.SNLOC) THEN
        FLAbort('SNLOC2.NE.SNLOC IN CV_CTY_SGS')
     ENDIF
         
     ALLOCATE(   SVN(NLOC,SVNGI) )
     ALLOCATE( SVNLX(NLOC,SVNGI) )
     ALLOCATE( SVNLY(NLOC,SVNGI) )
     ALLOCATE(    SVWEIGH(SVNGI) )
        
     ALLOCATE(   CVM(NLOC,CVNGI) )
      
     ALLOCATE(   CVDETWEI(SVNGI) )
     ALLOCATE(    CVNORMX(SVNGI) )
     ALLOCATE(    CVNORMY(SVNGI) )
     ALLOCATE(    CVNORMZ(SVNGI) )
      
     ALLOCATE(  NEILOC(NLOC,SVNGI) )
     ALLOCATE( TEMPCOLGPTS(NLOC*SVNGI) ) 
! The size of this vector is over-estimated
     ALLOCATE(    FINDGPTS(NLOC+1) )
        
! define the sub-control volume shape functions, etc
     CALL SHAPESV( ELETYP,  NEILOC,   &
     &                CVNGI,     NLOC,    &
     &                SVNGI, &
     &                SVN, &
     &                SVNLX,   SVNLY,  &
     &                SVWEIGH, &
     &                CVM )
                      
! Define the gauss points that lie on the surface of the
! control volume surrounding a given local node (iloc)
     CALL GAUSSILOC( FINDGPTS, TEMPCOLGPTS, NCOLGPTS, &
     &                  NEILOC,   NLOC,    SVNGI     )
        
     IF(NCOLGPTS.GT.NLOC*SVNGI) THEN
        EWRITE(-1,*) 'WARNING: NCOLGPTS EXCEEDED NLOC*SVNGI', NCOLGPTS, NLOC*SVNGI
        FLAbort('Eek!')
     ENDIF
    
     ALLOCATE( COLGPTS(NCOLGPTS) )
        
     do GI = 1, NCOLGPTS
        COLGPTS(GI) = TEMPCOLGPTS(GI)
     ENDDO
        
     DEALLOCATE( TEMPCOLGPTS )
        
! allocate memory for surface integration
     ALLOCATE( SN(SNLOC,SNGI)   )
     ALLOCATE( SNLX(SNLOC,SNGI) )
     ALLOCATE( SNLY(SNLOC,SNGI) )
     ALLOCATE( SWEIGH(SNGI)     )
        
     ALLOCATE( SNEILOC(SNLOC,SNGI) )
        
     ALLOCATE(  SDETWEI(SNGI) )
     ALLOCATE(  SNORMXN(SNGI) )
     ALLOCATE(  SNORMYN(SNGI) )
     ALLOCATE(  SNORMZN(SNGI) )   
        
! get surface shape functions
     CALL SHAPESE( SELETYP,  SNEILOC,   &
     &                SNGI,     SNLOC, &
     &                SN, &
     &                SNLX,   SNLY, &
     &                SWEIGH )
                    
! Loop over elements 
     do ELE = 1,TOTELE
          
        do ILOC=1,NLOC 
                
           INOD = NDGLNO((ELE-1)*NLOC+ILOC)
! won't work when MLOC.NE.NLOC
                
! Loop over quadrature points
           do GCOUNT = FINDGPTS(ILOC),FINDGPTS(ILOC+1)-1
                  
              GI = COLGPTS(GCOUNT)
              
! get face normals
              CALL SCVDETNX( ELE,      GI,      INOD, &
     &                       NLOC,     SVNGI,   TOTELE, &
     &                       XONDGL,   XNONOD, &
     &                       CVDETWEI, CVNORMX, CVNORMY, &
     &                       CVNORMZ,  SVN,     SVNLX, &
     &                       SVNLY,    SVWEIGH, PNTX,   &
     &                       PNTY,     PNTZ,    X,  &
     &                       Y,        Z, &
     &                       D3,       DCYL,    GOTPNTS )
                            
! Get the opposite node for node ILOC and Gauss point GI
              do JLOC=1,NLOC
! won't work when MLOC.NE.NLOC  NEILOC(iloc,gi)
                
! Put results into the cts
                 SVNDET=SVN(JLOC,GI)*CVDETWEI(GI)
                
                 C1TDGELE(ELE,ILOC,JLOC)=C1TDGELE(ELE,ILOC,JLOC) + CVNORMX(GI)*SVNDET
                 C2TDGELE(ELE,ILOC,JLOC)=C2TDGELE(ELE,ILOC,JLOC) + CVNORMY(GI)*SVNDET
                 IF(D3) C3TDGELE(ELE,ILOC,JLOC)=C3TDGELE(ELE,ILOC,JLOC) + CVNORMZ(GI)*SVNDET
                
! END DO JLOC=1,NLOC...
              END DO
! END DO gcount = findgpts(iloc),findgpts(iloc+1)-1
           ENDDO
! END DO iloc = 1, nloc    
        ENDDO
          
! do we have a surface element on top or bottom of the domain...
!******************************************

        do IFACE=1,NFACE
          
! NB colele4 is ordered in terms of faces.
           ELE2=COLELE4((ELE-1)*5+IFACE)

! END DO DO IFACE=1,NFACE...
        END DO
! END DO  ele = 1, totele...
     ENDDO
        
   END SUBROUTINE CV_CTY_SGS

   SUBROUTINE DGDM123_DG_WEAK(&
     &                     NLOC,TOTELE,NDGLNO,NONODS,&
     &                     DGU,DGV,DGW,DT,&
     &                     NOBCU,NOBCV,NOBCW,&
     &                     BCU1,BCV1,BCW1,&
     &                     BCU2,BCV2,BCW2,&
     &                     DGMARKU,DGMARKV,DGMARKW,&
     &                     DIRI_BC_X,DIRI_BC_Y,DIRI_BC_Z)
! Amend DGDM1,DGDM2,DGDM3 to take into account b.c's
! Must be called before  APPLY_DG_ROT_BC_MAT
     REAL INFINY
     PARAMETER(INFINY=1.E+20)
     INTEGER NLOC,TOTELE,NDGLNO(TOTELE*NLOC),NONODS
     REAL DGU(TOTELE*NLOC),DGV(TOTELE*NLOC),DGW(TOTELE*NLOC),DT
     INTEGER NOBCU,NOBCV,NOBCW
     REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
     INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
     INTEGER DGMARKU(NLOC*TOTELE),DGMARKV(NLOC*TOTELE),DGMARKW(NLOC*TOTELE)
     REAL DIRI_BC_X(NLOC*TOTELE),DIRI_BC_Y(NLOC*TOTELE),DIRI_BC_Z(NLOC*TOTELE)
! Local variables
     INTEGER II,ILOC,NOD,GLOBI,IDG,ELE
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKU
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKV
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKW
     REAL, ALLOCATABLE, DIMENSION(:)::VALUE_U
     REAL, ALLOCATABLE, DIMENSION(:)::VALUE_V
     REAL, ALLOCATABLE, DIMENSION(:)::VALUE_W
         
     ALLOCATE(MARKU(NONODS))
     ALLOCATE(MARKV(NONODS))
     ALLOCATE(MARKW(NONODS))
     ALLOCATE(VALUE_U(NONODS))
     ALLOCATE(VALUE_V(NONODS))
     ALLOCATE(VALUE_W(NONODS))
         
     DGMARKU=0
     DGMARKV=0
     DGMARKW=0
         
     MARKU=.FALSE.
     MARKV=.FALSE.
     MARKW=.FALSE.
! Amend DGDM1 to include bc's:
     do II=1,NOBCU
        NOD=BCU2(II)
        MARKU(NOD)=.TRUE.
        VALUE_U(NOD)=BCU1(II)
     END DO
     do II=1,NOBCV
        NOD=BCV2(II)
        MARKV(NOD)=.TRUE.
        VALUE_V(NOD)=BCV1(II)
     END DO
     do II=1,NOBCW
        NOD=BCW2(II)
        MARKW(NOD)=.TRUE.
        VALUE_W(NOD)=BCW1(II)
     END DO

     do ELE=1,TOTELE
        do ILOC=1,NLOC
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           IDG=(ELE-1)*NLOC+ILOC
           IF(MARKU(GLOBI)) THEN
              DIRI_BC_X(IDG)=VALUE_U(GLOBI)
              DGMARKU(IDG)=1
           ENDIF
           IF(MARKV(GLOBI)) THEN
              DIRI_BC_Y(IDG)=VALUE_V(GLOBI)
              DGMARKV(IDG)=1
           ENDIF
           IF(MARKW(GLOBI)) THEN
              DIRI_BC_Z(IDG)=VALUE_W(GLOBI)
              DGMARKW(IDG)=1
           ENDIF
        END DO
     END DO
        
   END SUBROUTINE DGDM123_DG_WEAK

   SUBROUTINE UVW_ROT_DG(NLOC,BIG_NLOC,TOTELE,&
     &               TRANSP,DGU,DGV,DGW,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ)
! apply U=R U or if TRANSP then U=R^T U
     INTEGER NLOC,TOTELE
     INTEGER BIG_NLOC
     LOGICAL TRANSP
     REAL DGU(TOTELE*NLOC),DGV(TOTELE*NLOC),DGW(TOTELE*NLOC)
! place into rotations matrices:
! DGT1X,DGT1Y,DGT1Z, 
! DGT2X,DGT2Y,DGT2Z, 
! DGNX,DGNY,DGNZ,
! and b.c's  DGDM1,DGDM2,DGDM3
     LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
! Local variables...
     INTEGER SILOC,ILOC,SGI,DGNODI,ILOC3,JLOC3,NODDG,ELE
     REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
     REAL MATROT(BIG_NLOC,BIG_NLOC)
     REAL SOUSUBU(NLOC),SOUSUBV(NLOC),SOUSUBW(NLOC)
         
     IF(.NOT.ROTATDG) RETURN
         
     do ELE=1,TOTELE
        RLOC=0.0
        do ILOC=1,NLOC
           NODDG=(ELE-1)*NLOC+ILOC
           IF(ROT_THIS_DGNOD(NODDG)) THEN
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
              RLOC(ILOC,ILOC)       =DGT1X(NODDG)
              RLOC(ILOC,ILOC+NLOC)  =DGT1Y(NODDG)
              RLOC(ILOC,ILOC+2*NLOC)=DGT1Z(NODDG)
              RLOC(ILOC+NLOC,ILOC)       =DGT2X(NODDG)
              RLOC(ILOC+NLOC,ILOC+NLOC)  =DGT2Y(NODDG)
              RLOC(ILOC+NLOC,ILOC+2*NLOC)=DGT2Z(NODDG)
              RLOC(ILOC+2*NLOC,ILOC)       =DGNX(NODDG)
              RLOC(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(NODDG)
              RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(NODDG)
           ELSE
              RLOC(ILOC,ILOC)       =1.0
              RLOC(ILOC+NLOC,ILOC+NLOC)  =1.0
              RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=1.0
           ENDIF
        END DO
         
        IF(.NOT.TRANSP) THEN
           do ILOC3=1,BIG_NLOC
              do JLOC3=1,BIG_NLOC
                 RTLOC(ILOC3,JLOC3)=RLOC(JLOC3,ILOC3)
              END DO
           END DO
           MATROT=RTLOC
        ELSE
           MATROT=RLOC
        ENDIF
                
        do ILOC=1,NLOC
           NODDG=(ELE-1)*NLOC+ILOC
           IF(ROT_THIS_DGNOD(NODDG)) THEN
              SOUSUBU(ILOC)=DGU(NODDG)
              SOUSUBV(ILOC)=DGV(NODDG)
              SOUSUBW(ILOC)=DGW(NODDG)
                 
! Rotate rhs vecs and set to zero (for zero bc)...
              DGU(NODDG)=MATROT(ILOC,ILOC)*SOUSUBU(ILOC)&
     &                    +MATROT(ILOC,ILOC+NLOC)*SOUSUBV(ILOC)&
     &                    +MATROT(ILOC,ILOC+2*NLOC)*SOUSUBW(ILOC)   
              DGV(NODDG)=MATROT(ILOC+NLOC,ILOC)*SOUSUBU(ILOC)&
     &                    +MATROT(ILOC+NLOC,ILOC+NLOC)*SOUSUBV(ILOC)&
     &                    +MATROT(ILOC+NLOC,ILOC+2*NLOC)*SOUSUBW(ILOC)
              DGW(NODDG)=MATROT(ILOC+2*NLOC,ILOC)*SOUSUBU(ILOC)&
     &                    +MATROT(ILOC+2*NLOC,ILOC+NLOC)*SOUSUBV(ILOC)&
     &                    +MATROT(ILOC+2*NLOC,ILOC+2*NLOC)*SOUSUBW(ILOC)
           ENDIF
        END DO
     END DO

   END SUBROUTINE UVW_ROT_DG

   SUBROUTINE PACK_BIGMELE_VECS(TOTELE,NLOC,NBIGMELE,BIGMELE,&
     &                     COMAT11,COMAT12,COMAT13,COMAT21, &
     &                     COMAT22,COMAT23,COMAT31,COMAT32,COMAT33, &
     &                     DGT1X,DGT1Y,DGT1Z, &
     &                     DGT2X,DGT2Y,DGT2Z, &
     &                     DGNX,DGNY,DGNZ,&
     &                     DGDM1,DGDM2,DGDM3)
! Pack key vectors into BIGMELE for use in the matrix solver. 
     INTEGER NLOC,TOTELE,NBIGMELE
     REAL BIGMELE(NBIGMELE)
! COMAT11 is normally COMAT11(TOTELE,NLOC) but here is COMAT11(TOTELE*NLOC) 
! for memory managment.
         REAL COMAT11(TOTELE*NLOC),COMAT12(TOTELE*NLOC),COMAT13(TOTELE*NLOC), &
     &        COMAT21(TOTELE*NLOC),COMAT22(TOTELE*NLOC),COMAT23(TOTELE*NLOC), &
     &        COMAT31(TOTELE*NLOC),COMAT32(TOTELE*NLOC),COMAT33(TOTELE*NLOC), &
     &        DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC), &
     &        DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC), &
     &        DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC),&
     &        DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Coriolis...
! place into rotations matrices:
! DGT1X,DGT1Y,DGT1Z, 
! DGT2X,DGT2Y,DGT2Z, 
! DGNX,DGNY,DGNZ,
! and b.c's               DGDM1,DGDM2,DGDM3
! Local variables...
! Integer pointers...
         INTEGER ICOMAT11,ICOMAT12,ICOMAT13, &
     &           ICOMAT21,ICOMAT22,ICOMAT23, &
     &           ICOMAT31,ICOMAT32,ICOMAT33, &
     &           IDGT1X,IDGT1Y,IDGT1Z, &
     &           IDGT2X,IDGT2Y,IDGT2Z, &
     &           IDGNX,IDGNY,IDGNZ,&
     &           IDGDM1,IDGDM2,IDGDM3
     INTEGER RPT
     
     RPT=NBIGMELE-21*TOTELE*NLOC+1
         CALL ALOMEM(ICOMAT11,   RPT, TOTELE*NLOC, NBIGMELE)
         CALL ALOMEM(ICOMAT12,   RPT, TOTELE*NLOC, NBIGMELE)
         CALL ALOMEM(ICOMAT13,   RPT, TOTELE*NLOC, NBIGMELE)
         
         CALL ALOMEM(ICOMAT21,   RPT, TOTELE*NLOC, NBIGMELE)
         CALL ALOMEM(ICOMAT22,   RPT, TOTELE*NLOC, NBIGMELE)
         CALL ALOMEM(ICOMAT23,   RPT, TOTELE*NLOC, NBIGMELE)
         
         CALL ALOMEM(ICOMAT31,   RPT, TOTELE*NLOC, NBIGMELE)
         CALL ALOMEM(ICOMAT32,   RPT, TOTELE*NLOC, NBIGMELE)
         CALL ALOMEM(ICOMAT33,   RPT, TOTELE*NLOC, NBIGMELE)
    
     CALL ALOMEM(IDGT1X,   RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGT1Y,   RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGT1Z,   RPT, TOTELE*NLOC, NBIGMELE)
     
     CALL ALOMEM(IDGT2X,   RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGT2Y,   RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGT2Z,   RPT, TOTELE*NLOC, NBIGMELE)
         
     CALL ALOMEM(IDGNX,    RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGNY,    RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGNZ,    RPT, TOTELE*NLOC, NBIGMELE)
         
     CALL ALOMEM(IDGDM1,   RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGDM2,   RPT, TOTELE*NLOC, NBIGMELE)
     CALL ALOMEM(IDGDM3,   RPT, TOTELE*NLOC, NBIGMELE)
         
         BIGMELE(ICOMAT11:ICOMAT11+TOTELE*NLOC-1)=COMAT11(1:TOTELE*NLOC)
         BIGMELE(ICOMAT12:ICOMAT12+TOTELE*NLOC-1)=COMAT12(1:TOTELE*NLOC)
         BIGMELE(ICOMAT13:ICOMAT13+TOTELE*NLOC-1)=COMAT13(1:TOTELE*NLOC)
         
         BIGMELE(ICOMAT21:ICOMAT21+TOTELE*NLOC-1)=COMAT21(1:TOTELE*NLOC)
         BIGMELE(ICOMAT22:ICOMAT22+TOTELE*NLOC-1)=COMAT22(1:TOTELE*NLOC)
         BIGMELE(ICOMAT23:ICOMAT23+TOTELE*NLOC-1)=COMAT23(1:TOTELE*NLOC)
         
         BIGMELE(ICOMAT31:ICOMAT31+TOTELE*NLOC-1)=COMAT31(1:TOTELE*NLOC)
         BIGMELE(ICOMAT32:ICOMAT32+TOTELE*NLOC-1)=COMAT32(1:TOTELE*NLOC)
         BIGMELE(ICOMAT33:ICOMAT33+TOTELE*NLOC-1)=COMAT33(1:TOTELE*NLOC)
         
     BIGMELE(IDGT1X:IDGT1X+TOTELE*NLOC-1)=DGT1X(1:TOTELE*NLOC)
     BIGMELE(IDGT1Y:IDGT1Y+TOTELE*NLOC-1)=DGT1Y(1:TOTELE*NLOC)
     BIGMELE(IDGT1Z:IDGT1Z+TOTELE*NLOC-1)=DGT1Z(1:TOTELE*NLOC)
         
     BIGMELE(IDGT2X:IDGT2X+TOTELE*NLOC-1)=DGT2X(1:TOTELE*NLOC)
     BIGMELE(IDGT2Y:IDGT2Y+TOTELE*NLOC-1)=DGT2Y(1:TOTELE*NLOC)
     BIGMELE(IDGT2Z:IDGT2Z+TOTELE*NLOC-1)=DGT2Z(1:TOTELE*NLOC)
         
     BIGMELE(IDGNX:IDGNX+TOTELE*NLOC-1)=DGNX(1:TOTELE*NLOC)
     BIGMELE(IDGNY:IDGNY+TOTELE*NLOC-1)=DGNY(1:TOTELE*NLOC)
     BIGMELE(IDGNZ:IDGNZ+TOTELE*NLOC-1)=DGNZ(1:TOTELE*NLOC)
         
     BIGMELE(IDGDM1:IDGDM1+TOTELE*NLOC-1)=DGDM1(1:TOTELE*NLOC)
     BIGMELE(IDGDM2:IDGDM2+TOTELE*NLOC-1)=DGDM2(1:TOTELE*NLOC)
     BIGMELE(IDGDM3:IDGDM3+TOTELE*NLOC-1)=DGDM3(1:TOTELE*NLOC)

     ewrite(1,*) 'allocated all of bigmele'
         
   END SUBROUTINE PACK_BIGMELE_VECS

   SUBROUTINE APPLY_DG_ROT_BC_MAT(BIGMELE,CENTRMELE,FINDRMELE,&
     &                     COLMELE,NCOLMELE,NBIGMELE,&
     &                     NLOC,BIG_NLOC,TOTELE,&
     &                     ROTATDG,ROT_THIS_DGNOD,&
     &                     DGT1X,DGT1Y,DGT1Z, &
     &                     DGT2X,DGT2Y,DGT2Z, &
     &                     DGNX,DGNY,DGNZ,&
     &                     DGDM1,DGDM2,DGDM3)
! Apply DG b.c.'s also rotate R A R^T
! and apply VEC=R VEC
! Must be called after DGDM123_DG_ROT_BC
     REAL INFINY
     PARAMETER(INFINY=1.E+20)
     INTEGER NLOC,TOTELE
     INTEGER BIG_NLOC
     INTEGER CENTRMELE(TOTELE),FINDRMELE(TOTELE+1)
     INTEGER NBIGMELE,NCOLMELE,COLMELE(NCOLMELE)
     REAL BIGMELE(NBIGMELE)
! place into rotations matrices:
! DGT1X,DGT1Y,DGT1Z, 
! DGT2X,DGT2Y,DGT2Z, 
! DGNX,DGNY,DGNZ,
! and b.c's               DGDM1,DGDM2,DGDM3
     LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Local variables...
     INTEGER ILOC,JLOC,ILOC3,JLOC3,ELE,ELE2
     INTEGER COUNT,IDG,JDG,GLOBI,NOD,ID,JD,II
     INTEGER BIG_ILOC,BIG_JLOC,IJDISP,NBLOCK
     REAL RLOC(BIG_NLOC,BIG_NLOC),RTLOC(BIG_NLOC,BIG_NLOC)
     REAL RLOC2(BIG_NLOC,BIG_NLOC)
     REAL RMATBLOKRT(BIG_NLOC,BIG_NLOC),MATBLOKRT(BIG_NLOC,BIG_NLOC)
     REAL MATBLOK(BIG_NLOC,BIG_NLOC)

     NBLOCK=(3*NLOC)**2
         
     IF(ROTATDG) THEN
         
        do ELE=1,TOTELE
           do COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
              ELE2=COLMELE(COUNT)
           
              RLOC=0.0
              do ILOC=1,NLOC
                 IDG=(ELE-1)*NLOC+ILOC
                 IF(ROT_THIS_DGNOD(IDG)) THEN
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
                    RLOC(ILOC,ILOC)       =DGT1X(IDG)
                    RLOC(ILOC,ILOC+NLOC)  =DGT1Y(IDG)
                    RLOC(ILOC,ILOC+2*NLOC)=DGT1Z(IDG)
                    RLOC(ILOC+NLOC,ILOC)       =DGT2X(IDG)
                    RLOC(ILOC+NLOC,ILOC+NLOC)  =DGT2Y(IDG)
                    RLOC(ILOC+NLOC,ILOC+2*NLOC)=DGT2Z(IDG)
                    RLOC(ILOC+2*NLOC,ILOC)       =DGNX(IDG)
                    RLOC(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(IDG)
                    RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(IDG)
                 ELSE
                    RLOC(ILOC,ILOC)              =1.0
                    RLOC(ILOC+NLOC,ILOC+NLOC)    =1.0
                    RLOC(ILOC+2*NLOC,ILOC+2*NLOC)=1.0
                 ENDIF
              END DO
            
              RLOC2=0.0
              do JLOC=1,NLOC
                 JDG=(ELE2-1)*NLOC+JLOC
                 IF(ROT_THIS_DGNOD(JDG)) THEN
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
                    RLOC2(ILOC,ILOC)       =DGT1X(JDG)
                    RLOC2(ILOC,ILOC+NLOC)  =DGT1Y(JDG)
                    RLOC2(ILOC,ILOC+2*NLOC)=DGT1Z(JDG)
                    RLOC2(ILOC+NLOC,ILOC)       =DGT2X(JDG)
                    RLOC2(ILOC+NLOC,ILOC+NLOC)  =DGT2Y(JDG)
                    RLOC2(ILOC+NLOC,ILOC+2*NLOC)=DGT2Z(JDG)
                    RLOC2(ILOC+2*NLOC,ILOC)       =DGNX(JDG)
                    RLOC2(ILOC+2*NLOC,ILOC+NLOC)  =DGNY(JDG)
                    RLOC2(ILOC+2*NLOC,ILOC+2*NLOC)=DGNZ(JDG)
                 ELSE
                    RLOC2(ILOC,ILOC)              =1.0
                    RLOC2(ILOC+NLOC,ILOC+NLOC)    =1.0
                    RLOC2(ILOC+2*NLOC,ILOC+2*NLOC)=1.0
                 ENDIF
              END DO
              do ILOC3=1,BIG_NLOC
                 do JLOC3=1,BIG_NLOC
                    RTLOC(ILOC3,JLOC3)=RLOC2(JLOC3,ILOC3)
                 END DO
              END DO
! Multiply block by R block R^T
              do ILOC3=1,NLOC
                 do JLOC3=1,NLOC
                    do ID=1,3
                       BIG_ILOC=(ID-1)*NLOC+ILOC3
                       do JD=1,3
                          BIG_JLOC=(JD-1)*NLOC+JLOC3
                          IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
                          MATBLOK(BIG_ILOC,BIG_JLOC)=BIGMELE((COUNT-1)*NBLOCK+IJDISP)
                       END DO
                    END DO
                 END DO
              END DO
! R MATBLOK R^T:
! Calculate rotated b.c's...
              CALL ABMATRIXMUL(MATBLOKRT, MATBLOK,  BIG_NLOC,BIG_NLOC, &
     &                                  RTLOC,    BIG_NLOC,BIG_NLOC)
              CALL ABMATRIXMUL(RMATBLOKRT,RLOC,     BIG_NLOC,BIG_NLOC, &
     &                                  MATBLOKRT,BIG_NLOC,BIG_NLOC)
! Copy back into BIGMELE:
              do ILOC3=1,NLOC
                 do JLOC3=1,NLOC
                    do ID=1,3
                       BIG_ILOC=(ID-1)*NLOC+ILOC3
                       do JD=1,3
                          BIG_JLOC=(JD-1)*NLOC+JLOC3
                          IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
                          BIGMELE((COUNT-1)*NBLOCK+IJDISP)=RMATBLOKRT(BIG_ILOC,BIG_JLOC)
                       END DO
                    END DO
                 END DO
              END DO
! NEXT COUNT...
           END DO
! NEXT ELE...
        END DO
! ENDOF IF(ROTATDG) THEN
     ENDIF

! Place DGDM1, on diagonal...
     do ELE=1,TOTELE
        COUNT=CENTRMELE(ELE)
        do ILOC3=1,NLOC
           JLOC3=ILOC3
           ID=1
           IDG=(ELE-1)*NLOC+ILOC3
           BIG_ILOC=(ID-1)*NLOC+ILOC3
           JD=ID
           BIG_JLOC=(JD-1)*NLOC+JLOC3
           IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
           BIGMELE((COUNT-1)*NBLOCK+IJDISP)=BIGMELE((COUNT-1)*NBLOCK+IJDISP)+DGDM1(IDG)

           ID=2
           IDG=(ELE-1)*NLOC+ILOC3
           BIG_ILOC=(ID-1)*NLOC+ILOC3
           JD=ID
           BIG_JLOC=(JD-1)*NLOC+JLOC3
           IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
           BIGMELE((COUNT-1)*NBLOCK+IJDISP)=BIGMELE((COUNT-1)*NBLOCK+IJDISP)+DGDM2(IDG)

           ID=3
           IDG=(ELE-1)*NLOC+ILOC3
           BIG_ILOC=(ID-1)*NLOC+ILOC3
           JD=ID
           BIG_JLOC=(JD-1)*NLOC+JLOC3
           IJDISP=(BIG_ILOC-1)*3*NLOC + BIG_JLOC
           BIGMELE((COUNT-1)*NBLOCK+IJDISP)=BIGMELE((COUNT-1)*NBLOCK+IJDISP)+DGDM3(IDG)
        END DO
     END DO

   END SUBROUTINE APPLY_DG_ROT_BC_MAT

   SUBROUTINE DGDM123_DG_ROT_BC(&
     &                     NLOC,BIG_NLOC,TOTELE,NDGLNO,NONODS,&
     &                     VECX,VECY,VECZ,&
     &                     DGU,DGV,DGW,DT,&
     &                     NOBCU,NOBCV,NOBCW,&
     &                     BCU1,BCV1,BCW1,&
     &                     BCU2,BCV2,BCW2, &
     &                     ROTATDG,ROT_THIS_DGNOD,&
     &                     DGT1X,DGT1Y,DGT1Z, &
     &                     DGT2X,DGT2Y,DGT2Z, &
     &                     DGNX,DGNY,DGNZ,&
     &                     DGDM1,DGDM2,DGDM3)
! Amend DGDM1,DGDM2,DGDM3 to take into account b.c's
! Must be called before  APPLY_DG_ROT_BC_MAT
     REAL INFINY
     PARAMETER(INFINY=1.E+20)
     INTEGER NLOC,BIG_NLOC,TOTELE,NDGLNO(TOTELE*NLOC),NONODS
     REAL VECX(TOTELE*NLOC),VECY(TOTELE*NLOC),VECZ(TOTELE*NLOC)
     REAL DGU(TOTELE*NLOC),DGV(TOTELE*NLOC),DGW(TOTELE*NLOC),DT
     INTEGER NOBCU,NOBCV,NOBCW
     REAL BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
     INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
! place into rotations matrices:
! DGT1X,DGT1Y,DGT1Z, 
! DGT2X,DGT2Y,DGT2Z, 
! DGNX,DGNY,DGNZ,
! and b.c's               DGDM1,DGDM2,DGDM3
     LOGICAL ROTATDG,ROT_THIS_DGNOD(TOTELE*NLOC)
     REAL DGT1X(TOTELE*NLOC),DGT1Y(TOTELE*NLOC),DGT1Z(TOTELE*NLOC)
     REAL DGT2X(TOTELE*NLOC),DGT2Y(TOTELE*NLOC),DGT2Z(TOTELE*NLOC)
     REAL DGNX(TOTELE*NLOC),DGNY(TOTELE*NLOC),DGNZ(TOTELE*NLOC)
     REAL DGDM1(TOTELE*NLOC),DGDM2(TOTELE*NLOC),DGDM3(TOTELE*NLOC)
! Local variables
     INTEGER II,ILOC,NOD,GLOBI,IDG,ELE
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKU
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKV
     LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKW
     REAL, ALLOCATABLE, DIMENSION(:)::VALUE_U
     REAL, ALLOCATABLE, DIMENSION(:)::VALUE_V
     REAL, ALLOCATABLE, DIMENSION(:)::VALUE_W
         
     ALLOCATE(MARKU(NONODS))
     ALLOCATE(MARKV(NONODS))
     ALLOCATE(MARKW(NONODS))
     ALLOCATE(VALUE_U(NONODS))
     ALLOCATE(VALUE_V(NONODS))
     ALLOCATE(VALUE_W(NONODS))
         
     IF(ROTATDG) THEN
! Rotate VEC=R VEC...
        CALL UVW_ROT_DG(NLOC,BIG_NLOC,TOTELE,&
     &               .FALSE.,VECX,VECY,VECZ,&
     &               ROTATDG,ROT_THIS_DGNOD,&
     &               DGT1X,DGT1Y,DGT1Z, &
     &               DGT2X,DGT2Y,DGT2Z, &
     &               DGNX,DGNY,DGNZ)
! ENDOF IF(ROTATDG) THEN
     ENDIF
         
     MARKU=.FALSE.
     MARKV=.FALSE.
     MARKW=.FALSE.
! Amend DGDM1 to include bc's:
     do II=1,NOBCU
        NOD=BCU2(II)
        MARKU(NOD)=.TRUE.
        VALUE_U(NOD)=BCU1(II)
     END DO
     do II=1,NOBCV
        NOD=BCV2(II)
        MARKV(NOD)=.TRUE.
        VALUE_V(NOD)=BCV1(II)
     END DO
     do II=1,NOBCW
        NOD=BCW2(II)
        MARKW(NOD)=.TRUE.
        VALUE_W(NOD)=BCW1(II)
     END DO

     do ELE=1,TOTELE
        do ILOC=1,NLOC
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           IDG=(ELE-1)*NLOC+ILOC
           IF(MARKU(GLOBI)) THEN
              DGDM1(IDG)=INFINY
              VECX(IDG)=INFINY*(DGU(IDG)+DT*VALUE_U(GLOBI))
           ENDIF
           IF(MARKV(GLOBI)) THEN
              DGDM2(IDG)=INFINY
              VECY(IDG)=INFINY*(DGV(IDG)+DT*VALUE_V(GLOBI))
           ENDIF
           IF(MARKW(GLOBI)) THEN
              DGDM3(IDG)=INFINY
              VECZ(IDG)=INFINY*(DGW(IDG)+DT*VALUE_W(GLOBI))
           ENDIF
        END DO
     END DO
        
   END SUBROUTINE DGDM123_DG_ROT_BC

   SUBROUTINE SOLV_FIELD_SGSDG(T,TTNEW,&
     &     NU,NV,NW,UG,VG,WG,&
     &     SOURCT,X,Y,Z,D0,&
     &     SUBSOUT,ISUBSOU,&
     &     VECT_SGSADD,&
     &     BIGM, ML,&
     &     TOTELE, &
     &     FINDRM,COLM,NCOLM,NONODS,CENTRM, &
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT2,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     NNODP,&
! for the SGS...
     &     NSUBVLOC,&
     &     NUSUB,NVSUB,NWSUB,INUSUB,&
     &     TSUB,TNEWSUB,ISUB,&
     &     NBIGM,NOBCFSimp,BCT2FSimp,&
     &     SUF_TOP_OCEAN,&
     &     NOBCT,BCT1,BCT2,&
! GMRES solver...
     &     PARA,halo_tag,&
! form DGCMC
     &     D3, velocity_option_path) 
! This code solves and SGS momentum and cty eqns...
     INTEGER MATSTR
     PARAMETER(MATSTR=0)
! if LIMIT_DGT BOUND THE VALUE OF DGT
     INTEGER SNONOD,VNONOD,XNONOD,NCOLM,NLOC,NGI
     INTEGER DISOPT2,DISOPT,DISOPN,DISCRA
     REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
     REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
     REAL DENPT(SNONOD)
     REAL ABSORB(SNONOD)
     REAL  DT,THETA,BETA
     INTEGER TOTELE,NONODS,FREDOP,NBIGM,VERSIO,ISPHERE
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
     REAL T(NONODS)
     REAL TTNEW(NONODS)
     INTEGER ISUB,INUSUB,NSUBVLOC
     REAL NUSUB(TOTELE*NSUBVLOC),NVSUB(TOTELE*NSUBVLOC)
     REAL NWSUB(TOTELE*NSUBVLOC)
     REAL TSUB(TOTELE*NLOC)
     REAL TNEWSUB(TOTELE*NLOC)
     REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
     REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
! This is what we add into the SGS model - a discretised source.      
     REAL VECT_SGSADD(TOTELE*NLOC)
     REAL BIGM(NBIGM)
     LOGICAL SUF_TOP_OCEAN
     INTEGER NOBCFSimp
     INTEGER BCT2FSimp(NOBCFSimp)
     INTEGER NOBCT
     INTEGER BCT2(NOBCT)
     REAL BCT1(NOBCT)
     REAL ML(NONODS)
     REAL X(XNONOD),Y(XNONOD),Z(XNONOD),D0
     REAL SOURCT(SNONOD)
     INTEGER ISUBSOU
     REAL SUBSOUT(ISUBSOU*TOTELE*NLOC)
     INTEGER NDGLNO(TOTELE*NLOC)
     INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
     INTEGER XONDGL(TOTELE*NLOC)
      
     INTEGER FINDRM(NONODS+1),COLM(NCOLM)
     INTEGER CENTRM(NONODS)

     REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI)
     REAL, target:: NLY(NLOC,NGI),NLZ(NLOC,NGI)
     REAL WEIGHT(NGI)
     LOGICAL LUMP,UPWIND,OPWIND,XTPETV
     INTEGER PREOPT2
     LOGICAL ABSLUM,SLUMP 
     INTEGER NNODP
! GMRES solver...
     INTEGER PARA,halo_tag
     REAL ERROR1,ERROR2
     INTEGER NOITS1,NOITS2
! Pressure matrix soln:
     INTEGER PREPRE
     INTEGER PMISOU,halo_tag_p
     REAL PREERR
     INTEGER PRENOI,PINNOI,NNODPP
     LOGICAL D3
     LOGICAL LIMIT_DGT
     character(len=*), intent(in):: velocity_option_path
! Local vartiables...
       INTEGER IDUM(1),KITS,NR2,NVN,ILENG,PREOPT
       INTEGER ELE,ELE2,ILOC,JLOC,NODDG,NOD,I,count,ILE,GLOBI,IDG,IDG1,MATSTRVEL
       INTEGER NCOLMELE,NBIGMELE,ITS
       INTEGER IDISP,IJBIGM,ISTAB
       REAL RDUM(1),RMIN,MAXDGTOLD,MINDGTOLD,minmin,maxmax,minminnew,maxmaxnew
       LOGICAL ROTATDG,NO_IN_FLO
! Function to decode PREPRE for free surface
      INTEGER DECOPREPRE,NITS
      
     REAL, ALLOCATABLE, DIMENSION(:)::VECT
     REAL, ALLOCATABLE, DIMENSION(:,:,:)::MASSDGELE,MASSDGELE2
     INTEGER, ALLOCATABLE, DIMENSION(:)::FINDRMELE
     INTEGER, ALLOCATABLE, DIMENSION(:)::COLMELE
     INTEGER, ALLOCATABLE, DIMENSION(:)::CENTRMELE
     REAL, ALLOCATABLE, DIMENSION(:)::DGT,DGTNEW,DGNU,DGNV,DGNW
     REAL, ALLOCATABLE, DIMENSION(:)::UG2,VG2,WG2,DGNU2,DGNV2,DGNW2
       REAL, ALLOCATABLE, DIMENSION(:,:,:)::STAB_DIF_MAT,STORE_BDIA_MAT
       REAL, ALLOCATABLE, DIMENSION(:)::STORE_RHS_VEC
       REAL, ALLOCATABLE, DIMENSION(:)::BIGMELE,MLELE,MLELE2
       REAL, ALLOCATABLE, DIMENSION(:)::TZER,WORKV,WORKW,WORKVSUB,WORKWSUB
       REAL, ALLOCATABLE, DIMENSION(:)::DGTOLD,DGTNEWNEW
       
     IF((BETA.LT.0.0).OR.(BETA.GT.1.001)) THEN
        ewrite(-1,*) 'BETA IS NOT IN CORRECT RANG BETA=',BETA
        ewrite(-1,*) 'NEED TO ADD 10.0 TO ITS NORMAL VALUE IN GEM FILE'
        FLAbort('Bleugh!')
     ENDIF
     ILENG=TOTELE*NLOC
     ALLOCATE(VECT(TOTELE*NLOC))
     ALLOCATE(MASSDGELE(TOTELE,NLOC,NLOC))
     ALLOCATE(MASSDGELE2(TOTELE,NLOC,NLOC))
! Form matrix structure...
     ALLOCATE(FINDRMELE(TOTELE+1))
     ALLOCATE(COLMELE(TOTELE*5))
     ALLOCATE(CENTRMELE(TOTELE))
     CALL SIMP_DG_MAT_STRU(FINDRMELE,COLMELE,NCOLMELE,TOTELE,NLOC,&
     &                       NDGLNO,CENTRMELE,NONODS)
! Map variables to DG variles...
     ALLOCATE(DGT(TOTELE*NLOC))
     ALLOCATE(DGTNEW(TOTELE*NLOC))
     ALLOCATE(DGTOLD(TOTELE*NLOC))
     ALLOCATE(DGTNEWNEW(TOTELE*NLOC))
     ALLOCATE(DGNU(TOTELE*NLOC))
     ALLOCATE(DGNV(TOTELE*NLOC))
     ALLOCATE(DGNW(TOTELE*NLOC))
     ALLOCATE(MLELE(TOTELE*NLOC))
     ALLOCATE(MLELE2(TOTELE*NLOC))
     minmin=1.e+20
     maxmax=-1.e+20
     minminnew=1.e+20
     maxmaxnew=-1.e+20
     do ELE=1,TOTELE
        do ILOC=1,NLOC
           GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
           IDG=(ELE-1)*NLOC+ILOC
           DGT(IDG)=T(GLOBI)+TSUB(IDG)
           DGTNEW(IDG)=TTNEW(GLOBI)+TNEWSUB(IDG)
           minmin=min(minmin,dgt(idg))
           maxmax=max(maxmax,dgt(idg))
           minminnew=min(minminnew,dgtnew(idg))
           maxmaxnew=max(maxmaxnew,dgtnew(idg))
           IF(NSUBVLOC*INUSUB.NE.0) THEN
              DGNU(IDG)=NU(GLOBI)+NUSUB(IDG)
              DGNV(IDG)=NV(GLOBI)+NVSUB(IDG)
              DGNW(IDG)=NW(GLOBI)+NWSUB(IDG)
           ELSE
              DGNU(IDG)=NU(GLOBI)
              DGNV(IDG)=NV(GLOBI)
              DGNW(IDG)=NW(GLOBI)
           ENDIF
        END DO
     END DO
     ewrite(3,*) 'just in advection sub'
     ewrite(3,*) 'minmin,maxmax,minminnew,maxmaxnew:',minmin,maxmax,minminnew,maxmaxnew
     ALLOCATE(DGNU2(TOTELE*NLOC))
     ALLOCATE(DGNV2(TOTELE*NLOC))
     ALLOCATE(DGNW2(TOTELE*NLOC))
     ALLOCATE(UG2(NONODS))
     ALLOCATE(VG2(NONODS))
     ALLOCATE(WG2(NONODS))
     NO_IN_FLO=.true.
     IF(NO_IN_FLO) THEN
        UG2=0.0
        VG2=0.0
        WG2=0.0
        do ELE=1,TOTELE
           do ILOC=1,NLOC
              GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
              IDG=(ELE-1)*NLOC+ILOC
              DGNU2(IDG)=DGNU(IDG)-UG(GLOBI)
              DGNV2(IDG)=DGNV(IDG)-VG(GLOBI)
              DGNW2(IDG)=DGNW(IDG)-WG(GLOBI)
           END DO
        END DO
        CALL MAKE_NO_INFLOW(TOTELE,NLOC,XONDGL,3,XNONOD,&
     &                       DGNU2,DGNV2,DGNW2,X,Y,Z)
     ELSE
        UG2=UG
        VG2=VG
        WG2=WG
        DGNU2=DGNU
        DGNV2=DGNV
        DGNW2=DGNW
     ENDIF
     NBIGMELE=NLOC*NLOC*NCOLMELE
       
       DGTOLD=DGT
       DGTNEWNEW=DGTNEW
!       DGTNEWNEW=DGT

       
! if LIMIT_DGT BOUND THE VALUE OF DGT
!     LIMIT_DGT=DISOPT2.EQ.8
     LIMIT_DGT=(DISOPT2.EQ.8).or.(DISOPT2.EQ.158)
     
       NITS=1
       ISTAB=0
       IF((DISOPT2.EQ.157).OR.(DISOPT2.EQ.158)) THEN
         NITS=4 
!         NITS=20
         ISTAB=1
!         ISTAB=0
       ENDIF
       ALLOCATE(STAB_DIF_MAT(TOTELE*ISTAB,NLOC,NLOC))
       ALLOCATE(STORE_BDIA_MAT(TOTELE*ISTAB,NLOC,NLOC))
       ALLOCATE(STORE_RHS_VEC(TOTELE*ISTAB*NLOC))
       ALLOCATE(BIGMELE(NBIGMELE))
       
       
       DO ITS=1,NITS

! form momentum eqns and put pressure into rhs,
! also put b.c's into matrix...
       DGT=DGTOLD
     CALL HART3DSGS_DG(DGT,DGTNEWNEW,&
     &     DGNU2,DGNV2,DGNW2,UG2,VG2,WG2,&
     &     SOURCT,X,Y,Z,D0,&
     &     SUBSOUT,ISUBSOU,&
     &     VECT_SGSADD, &
     &     VECT,&
     &     BIGMELE,MLELE,ML,TOTELE, &
     &     FINDRMELE,COLMELE,NCOLMELE,NONODS,CENTRMELE, &
     &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
     &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT,&
     &     ABSLUM,SLUMP,&
     &     NDGLNO,&
     &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
     &     DENPT,&
     &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
     &     ABSORB,VERSIO,ISPHERE,&
     &     nnodp,para,halo_tag,&
! for the DG...
     &     NBIGMELE,NOBCFSimp,BCT2FSimp,&
     &     MASSDGELE,&
     &     SUF_TOP_OCEAN,(ITS.EQ.1),MIN(ITS-1,ISTAB),STAB_DIF_MAT, velocity_option_path)
     
        IF(ITS.EQ.1) THEN
          MLELE2=MLELE
          MASSDGELE2=MASSDGELE
        ENDIF
        
        IF(ISTAB.EQ.1) THEN
          IF(ITS.EQ.1) THEN
! Store the block diagonal of the matrix and the rhs vector...
            DO ELE=1,TOTELE
              IDISP=(CENTRMELE(ELE)-1)*NLOC*NLOC
              DO ILOC=1,NLOC
                STORE_RHS_VEC((ELE-1)*NLOC+ILOC)=VECT((ELE-1)*NLOC+ILOC)
                DO JLOC=1,NLOC
                  IJBIGM=IDISP+(ILOC-1)*NLOC+JLOC
                  STORE_BDIA_MAT(ELE,ILOC,JLOC)=BIGMELE(IJBIGM)
                END DO
              END DO
            END DO
          ELSE 
! place the stabilisation matrix STORE_RHS_VEC into matrix BIGMELE...
            DO ELE=1,TOTELE
              IDISP=(CENTRMELE(ELE)-1)*NLOC*NLOC
              DO ILOC=1,NLOC
                VECT((ELE-1)*NLOC+ILOC)=STORE_RHS_VEC((ELE-1)*NLOC+ILOC)
                DO JLOC=1,NLOC
                  IJBIGM=IDISP+(ILOC-1)*NLOC+JLOC
                  BIGMELE(IJBIGM)=STORE_BDIA_MAT(ELE,ILOC,JLOC)+THETA*STAB_DIF_MAT(ELE,ILOC,JLOC)
                  VECT((ELE-1)*NLOC+ILOC)=VECT((ELE-1)*NLOC+ILOC) &
     &                 -(1.-THETA)*STAB_DIF_MAT(ELE,ILOC,JLOC)*DGT((ELE-1)*NLOC+JLOC)
                END DO
              END DO
            END DO
          ENDIF
        ENDIF


     ewrite(1,*) 'going into APPLY_FIELD_DG_BC'
! Apply DG BC'S to field eqn...
     CALL APPLY_FIELD_DG_BC(&
     &     NLOC,TOTELE,NDGLNO,NONODS,&
     &     VECT,DGT,DT,NOBCT,BCT1,BCT2,&
     &     FINDRMELE,COLMELE,NCOLMELE,NBIGMELE,CENTRMELE,BIGMELE)
! Solve momentum eqns... 
     ILENG=TOTELE*NLOC
     ILE=TOTELE

     CALL PMINMX(DGT,TOTELE*NLOC,'******DGT  ')
       DGT=DGTNEWNEW
        
     MATSTRVEL=3
     CALL GMRES(DGT,VECT,ILE,ILENG,ILE,.FALSE.,&
     &       BIGMELE,FINDRMELE,COLMELE,NCOLMELE,NBIGMELE,&
     & halo_tag,KITS)
     
! BOUND THE VALUE OF DGT...    
          IF(LIMIT_DGT) THEN
            CALL SIMP_LIMIT_DGT(TOTELE,NLOC,DGT,DGTOLD, &
!            FINDRMELE,NCOLMELE,COLMELE, MLELE,MASSDGELE, &
            FINDRMELE,NCOLMELE,COLMELE, CENTRMELE, MLELE2,MASSDGELE2, &
            SONDGL,SNONOD,SOURCT,ABSORB,SUBSOUT,ISUBSOU,VECT_SGSADD, &
            THETA,DT,halo_tag)
! This sub limits the value of the field DGT based on the 
! only field to satisfy the max principle.
          ENDIF
     
        DGTNEWNEW=DGT

     CALL PMINMX(DGT,TOTELE*NLOC,'******DGT  ')
! next iteration...
       END DO
       
       DEALLOCATE(BIGMELE)
       
! Choose how to distribute dU_DG between global and SGS ********
     ALLOCATE(TZER(TOTELE*NLOC)) 
     TZER=0.0
     ALLOCATE(WORKV(NONODS))
     ALLOCATE(WORKW(NONODS))
     ALLOCATE(WORKVSUB(TOTELE*NLOC))
     ALLOCATE(WORKWSUB(TOTELE*NLOC))   
     T=0.0
     WORKV=0.0
     WORKW=0.0
     TSUB=0.0
     WORKVSUB=0.0
     WORKWSUB=0.0

     CALL GET_USUB_FROM_UDG(T,WORKV,WORKW,TSUB,WORKVSUB,WORKWSUB, &
     &         DGT,TZER,TZER,&
     &         NLOC,TOTELE,NDGLNO,&
     &         NONODS,FINDRM,COLM,CENTRM,NCOLM,&
     &         DT,D3,&
     &         MASSDGELE,&
     &         0,0,BIGM,&
     &         TZER,TZER,TZER,TZER,TZER,TZER,TZER,TZER,TZER, &
     &         NNODP,PARA,halo_tag,&
     &         0,0,0,BCT1,BCT1,BCT1,BCT2,BCT2,BCT2)
     
       TTNEW=T
       TNEWSUB=TSUB
       
       CALL print_DGT_maxmin(TOTELE,NLOC,TSUB,T,NDGLNO,NONODS,0)
     
       DEALLOCATE(TZER) 
       DEALLOCATE(WORKV)
       DEALLOCATE(WORKW)
       DEALLOCATE(WORKVSUB)
       DEALLOCATE(WORKWSUB)

   END SUBROUTINE SOLV_FIELD_SGSDG
   
   
        SUBROUTINE print_DGT_maxmin(TOTELE,NLOC,TSUB,T,NDGLNO,NONODS,POS)
        INTEGER TOTELE,NLOC,NONODS,POS
        INTEGER NDGLNO(TOTELE*NLOC)
        REAL T(NONODS),TSUB(TOTELE*NLOC)
! Local variables...
        INTEGER ELE,ILOC,IDG,NOD
        REAL RDGT,MINMIN,MAXMAX
          MINMIN=1.E+20
          MAXMAX=-1.E+20
            DO ELE=1,TOTELE
                DO ILOC=1,NLOC
                  IDG=(ELE-1)*NLOC+ILOC
                  NOD=ndglno(IDG)
                  RDGT=T(NOD)+TSUB(IDG)
                  MINMIN=MIN(MINMIN,RDGT)
                  MAXMAX=MAX(MAXMAX,RDGT)
                end do
            end do
         ewrite(3,*)'print_DGT_maxmin AT POS=',POS,MINMIN,MAXMAX
        END SUBROUTINE print_DGT_maxmin
        
   

     
        SUBROUTINE SIMP_LIMIT_DGT(TOTELE,NLOC,DGT,DGTOLD, &
          FINDRMELE,NCOLMELE,COLMELE, CENTRMELE,MLELE,MASSDGELE, &
          SONDGL,SNONOD,SOURCT,ABSORB,SUBSOUT,ISUBSOU,VECT_SGSADD, &
          THETA,DT,halo_tag)
! This sub limits the value of the field DGT based on the 
! only field to satisfy the max principle. 
! It tilts the contents of an element to achieve this 
! property.
! halo_tag is the velocity halo tag if this is -ve then it becomes 
! an element halo tag
        LOGICAL ORIG_LIMIT,BOUND_MEANS,NO_IN_LIMIT_REG
        INTEGER NIT_LIMIT_MEAN,NIT_LIMIT,halo_tag
        PARAMETER(ORIG_LIMIT=.false.,BOUND_MEANS=.TRUE.,NO_IN_LIMIT_REG=.TRUE.)
        PARAMETER(NIT_LIMIT_MEAN=50,NIT_LIMIT=50)
! If BOUND_MEANS Adjust the mean values to get boundedness in the means... 
! if ORIG_LIMIT then use original non-conservative limiter
! NIT_LIMIT=no. of limiting iterations. 
! If NO_IN_LIMIT_REG then do not smooth means into regions where boundedness 
! is violated or about to be violated.
        INTEGER TOTELE,NLOC,SNONOD,ISUBSOU
        REAL DGT(TOTELE*NLOC),DGTOLD(TOTELE*NLOC)
        INTEGER FINDRMELE(TOTELE+1),NCOLMELE,COLMELE(NCOLMELE),CENTRMELE(TOTELE)
        REAL MLELE(TOTELE*NLOC)
        REAL MASSDGELE(TOTELE,NLOC,NLOC)
        INTEGER SONDGL(TOTELE*NLOC)
        REAL SOURCT(SNONOD),ABSORB(SNONOD)
        REAL SUBSOUT(ISUBSOU*TOTELE*NLOC)
! This is what we add into the SGS model - a discretised source.      
        REAL VECT_SGSADD(TOTELE*NLOC)
        REAL THETA,DT
! Local variables...
        INTEGER ELE,ELE3,COUNT,COUNT2,ILOC,JLOC,ELE2,IDG,JDG,IT_LIMIT,NOD
        REAL MAXDGTOLD,MINDGTOLD,maxmax,minmin,SOU,ABS,RELAX
        REAL, ALLOCATABLE, DIMENSION(:)::DGT_PRED,DGT_MAX,DGT_MIN
        REAL, ALLOCATABLE, DIMENSION(:)::MEAN_DGT_MAX,MEAN_DGT_MIN
        REAL, ALLOCATABLE, DIMENSION(:)::MEAN_DGT_MAX2,MEAN_DGT_MIN2
       
     
           IF(.NOT.ORIG_LIMIT) THEN
          
            ALLOCATE(DGT_PRED(TOTELE*NLOC))
            
            ALLOCATE(DGT_MAX(TOTELE*NLOC))
            ALLOCATE(DGT_MIN(TOTELE*NLOC))
            
            ALLOCATE(MEAN_DGT_MAX(TOTELE))
            ALLOCATE(MEAN_DGT_MIN(TOTELE))
            ALLOCATE(MEAN_DGT_MAX2(TOTELE))
            ALLOCATE(MEAN_DGT_MIN2(TOTELE))
! Calculated a predicted value of DGT ie DGT_PRED based on no advection
! just for the boundedness calculation. 
     
            if(.TRUE.) then
            DO ELE=1,TOTELE
                DO ILOC=1,NLOC
                  IDG=(ELE-1)*NLOC+ILOC
                  NOD=SONDGL(IDG)
                  SOU=SOURCT(NOD)+VECT_SGSADD(IDG)/MLELE(IDG)
                  IF(ISUBSOU.NE.0) THEN
                    SOU=SOU+SUBSOUT(IDG)
                  ENDIF
                  ABS=ABSORB(NOD)
                  DGT_PRED(IDG)=((1./DT -(1.-THETA)*ABS)*DGTOLD(IDG) + SOU) &
                               /(1./DT + THETA*ABS)
                END DO
            END DO
            else
              DGT_PRED=DGTOLD
            endif
     
            DO ELE=1,TOTELE
              MAXDGTOLD=-1.E+20
              MINDGTOLD=+1.E+20
              DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
                  ELE2=COLMELE(COUNT)
                  DO ILOC=1,NLOC
                    MAXDGTOLD=MAX(MAXDGTOLD,DGT_PRED((ELE2-1)*NLOC+ILOC))
                    MINDGTOLD=MIN(MINDGTOLD,DGT_PRED((ELE2-1)*NLOC+ILOC))
                  END DO
              END DO
! CALCULATE max and min values for DGT...
              MEAN_DGT_MAX(ELE)=MAXDGTOLD
              MEAN_DGT_MIN(ELE)=MINDGTOLD
            END DO
!
! Gethalo of DGT_MAX_ELE,DGT_MIN_ELE
! get halo's -ver halo_tag for elements...
            call flcomms_update(-halo_tag, MEAN_DGT_MAX, 1, 1, 0)
            call flcomms_update(-halo_tag, MEAN_DGT_MIN, 1, 1, 0)
            
            DO ELE=1,TOTELE
              MAXDGTOLD=-1.E+20
              MINDGTOLD=+1.E+20
              DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
                  ELE2=COLMELE(COUNT)
                  MAXDGTOLD=MAX(MAXDGTOLD,MEAN_DGT_MAX(ELE2))
                  MINDGTOLD=MIN(MINDGTOLD,MEAN_DGT_MIN(ELE2))
              END DO
! CALCULATE max and min values for DGT...
              MEAN_DGT_MAX2(ELE)=MAXDGTOLD
              MEAN_DGT_MIN2(ELE)=MINDGTOLD
            END DO
! Redefine DGT_MAX_ELE,DGT_MIN_ELE... 
            MEAN_DGT_MAX=MEAN_DGT_MAX2
            MEAN_DGT_MIN=MEAN_DGT_MIN2 
! get halo's -ver halo_tag for elements...
            call flcomms_update(-halo_tag, MEAN_DGT_MAX, 1, 1, 0)
            call flcomms_update(-halo_tag, MEAN_DGT_MIN, 1, 1, 0)
            
            DO ELE=1,TOTELE
              DO ILOC=1,NLOC
                DGT_MAX((ELE-1)*NLOC+ILOC)=MEAN_DGT_MAX(ELE)
                DGT_MIN((ELE-1)*NLOC+ILOC)=MEAN_DGT_MIN(ELE)
              END DO
            END DO
            
! after adjustment work out max and min values...
            MAXMAX=-1.E+20
            MINMIN= 1.E+20
            DO ELE=1,TOTELE
              DO ILOC=1,NLOC
                MAXMAX=MAX(MAXMAX,DGT((ELE-1)*NLOC+ILOC))
                MINMIN=MIN(MINMIN,DGT((ELE-1)*NLOC+ILOC))
              END DO
            END DO
            ewrite(2,*) 'before ADJUSTMENT MINMIN,MAXMAX:',MINMIN,MAXMAX

! Adjust the mean values to get boundedness in the means...
            IF(BOUND_MEANS) THEN
              CALL BOUND_DG_MEANS(TOTELE,NLOC, DGT,MLELE, &
                  MEAN_DGT_MIN,MEAN_DGT_MAX, &
                  NCOLMELE,FINDRMELE,COLMELE,CENTRMELE, &
                  NO_IN_LIMIT_REG,NIT_LIMIT_MEAN,halo_tag)
            ENDIF
            
! Tilt soln only within an element to achieve boundedness...
            CALL BOUN_IN_ELE(TOTELE,NLOC,MASSDGELE,MLELE, &
                             DGT,DGT_MIN,DGT_MAX,NIT_LIMIT) 
            
            
! After adjustment work out max and min values...
            MAXMAX=-1.E+20
            MINMIN= 1.E+20
            DO ELE=1,TOTELE
              DO ILOC=1,NLOC
                MAXMAX=MAX(MAXMAX,DGT((ELE-1)*NLOC+ILOC))
                MINMIN=MIN(MINMIN,DGT((ELE-1)*NLOC+ILOC))
              END DO
            END DO
            PRINT *,'AFTER ADJUSTMENT MINMIN,MAXMAX:',MINMIN,MAXMAX
            
            DEALLOCATE(DGT_PRED)
            
            DEALLOCATE(DGT_MAX)
            DEALLOCATE(DGT_MIN)
            DEALLOCATE(MEAN_DGT_MAX)
            DEALLOCATE(MEAN_DGT_MIN)
            
           ELSE
! None conservative - use as test...          
            DO ELE=1,TOTELE
              MAXDGTOLD=-1.E+20
              MINDGTOLD=+1.E+20
              DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
                ELE2=COLMELE(COUNT)
                DO ILOC=1,NLOC
                  MAXDGTOLD=MAX(MAXDGTOLD,DGTOLD((ELE2-1)*NLOC+ILOC))
                  MINDGTOLD=MIN(MINDGTOLD,DGTOLD((ELE2-1)*NLOC+ILOC))
                END DO
              END DO
! Limit the value of DGT...
              DO ILOC=1,NLOC
                DGT((ELE-1)*NLOC+ILOC)=MIN(MAXDGTOLD,DGT((ELE-1)*NLOC+ILOC))
                DGT((ELE-1)*NLOC+ILOC)=MAX(MINDGTOLD,DGT((ELE-1)*NLOC+ILOC))
              END DO
            END DO
          ENDIF
          RETURN
          END SUBROUTINE SIMP_LIMIT_DGT
          
      

        SUBROUTINE BOUN_IN_ELE(TOTELE,NLOC,MASSDGELE,MLELE, &
                               DGT,DGT_MIN,DGT_MAX,NIT_LIMIT) 
! tilt soln only within an element to achieve boundedness of DGT...
        INTEGER TOTELE,NLOC,NIT_LIMIT
        REAL DGT(TOTELE*NLOC),DGT_MIN(TOTELE*NLOC),DGT_MAX(TOTELE*NLOC)
        REAL MASSDGELE(TOTELE,NLOC,NLOC),MLELE(TOTELE*NLOC)
! Local variables...
        INTEGER IT_LIMIT,ELE,ILOC,IDG,JLOC,JDG
        REAL, ALLOCATABLE, DIMENSION(:)::DGT_DEV,DGT_ALT
        
            ALLOCATE(DGT_DEV(TOTELE*NLOC))
            ALLOCATE(DGT_ALT(TOTELE*NLOC))
            
            DO IT_LIMIT=1,NIT_LIMIT
              DO ELE=1,TOTELE
! Limit the value of DGT...
! First calulate a field that defines the deviation from boundedness dgT_dev.
                DO ILOC=1,NLOC
                  IDG=(ELE-1)*NLOC+ILOC
                  if (DGT(IDG) .GT. DGT_max(IDG)) then
                    DGT_dev(IDG) = DGT(IDG) - DGT_max(IDG)
                  else if (DGT(IDG) .LT. DGT_min(IDG)) then
                    DGT_dev(IDG) = DGT(IDG) - DGT_min(IDG)
                  else
                    DGT_dev(IDG) = 0.0
                  end if 
                END DO
!Now spread this deviation T_dev around the surrounding nodes in an element
! M_L T_alt= M T_dev:
                DO ILOC=1,NLOC
                  IDG=(ELE-1)*NLOC+ILOC
                  DGT_alt(IDG)=0.0
                  DO JLOC=1,NLOC
                    JDG=(ELE-1)*NLOC+JLOC 
                    DGT_alt(IDG)=DGT_alt(IDG)+MASSDGELE(ELE,ILOC,JLOC)*DGT_DEV(JDG)
                  END DO
                  DGT_alt(IDG)=DGT_alt(IDG)/MLELE(IDG)
                END DO
! T_new=T_new-T_dev+T_alt:
                DO ILOC=1,NLOC
                  IDG=(ELE-1)*NLOC+ILOC
                  DGT(IDG)=DGT(IDG)-DGT_DEV(IDG)+DGT_alt(IDG)
                END DO
              END DO
              
! Next limiter iteration...
            END DO
        END SUBROUTINE BOUN_IN_ELE
        
                
        
        SUBROUTINE BOUND_DG_MEANS(TOTELE,NLOC, DGT,MLELE, &
          MEAN_DGT_MIN,MEAN_DGT_MAX, &
          NCOLMELE,FINDRMELE,COLMELE,CENTRMELE, &
          NO_IN_LIMIT_REG,NIT_LIMIT_MEAN,halo_tag)
! Adjust the mean values of DGT to achieve boundedness of the means...
        INTEGER TOTELE,NLOC
        REAL DGT(TOTELE*NLOC),MLELE(TOTELE*NLOC)
        REAL MEAN_DGT_MIN(TOTELE),MEAN_DGT_MAX(TOTELE)
        INTEGER NCOLMELE
        INTEGER FINDRMELE(TOTELE+1),COLMELE(NCOLMELE),CENTRMELE(TOTELE)
        LOGICAL NO_IN_LIMIT_REG
        INTEGER NIT_LIMIT_MEAN,halo_tag
! Local variables...
        INTEGER ELE,ILOC,IDG
        REAL SUM_MAS,MAXMAX,MINMIN
        REAL, ALLOCATABLE, DIMENSION(:)::MEAN_DGT,PREV_MEAN_DGT
        REAL, ALLOCATABLE, DIMENSION(:)::MATRIX
          
        ALLOCATE(MEAN_DGT(TOTELE))
        ALLOCATE(PREV_MEAN_DGT(TOTELE))
        ALLOCATE(MATRIX(NCOLMELE))

! Calculate the means...
        DO ELE=1,TOTELE
          MEAN_DGT(ELE)=0.0
          SUM_MAS=0.0
          DO ILOC=1,NLOC
            IDG=(ELE-1)*NLOC+ILOC
            MEAN_DGT(ELE)=MEAN_DGT(ELE)+DGT(IDG)*MLELE(IDG)
            SUM_MAS=SUM_MAS+MLELE(IDG)
          END DO
          MEAN_DGT(ELE)=MEAN_DGT(ELE)/SUM_MAS
        END DO
        
        MAXMAX=-1.0E+20
        MINMIN=+1.0E+20
        DO ELE=1,TOTELE
          MAXMAX=MAX(MAXMAX,MEAN_DGT(ELE))
          MINMIN=MIN(MINMIN,MEAN_DGT(ELE))
        END DO
        PRINT *,'BEFORE ADJUSTING THE MEANS MEAN MIN,MAX:',MINMIN,MAXMAX
          

! Calculate element centre distributed mass matrix MATRIX...
        CALL MAK_MASS_FROM_LUMPED(MLELE,MATRIX,TOTELE,NLOC, &
          NCOLMELE,FINDRMELE,COLMELE,CENTRMELE)
          
        PREV_MEAN_DGT=MEAN_DGT
! Obtain bounded soln fo MEAN_DGT...
        CALL BOUND_FIX(MATRIX,TOTELE,NLOC, &
          MEAN_DGT,MEAN_DGT_MIN,MEAN_DGT_MAX, &
          NCOLMELE,FINDRMELE,COLMELE,CENTRMELE, &
          NO_IN_LIMIT_REG,NIT_LIMIT_MEAN,halo_tag)
! now update DGT
        DO ELE=1,TOTELE
          DO ILOC=1,NLOC
            IDG=(ELE-1)*NLOC+ILOC
            DGT(IDG)=DGT(IDG)+(MEAN_DGT(ELE)-PREV_MEAN_DGT(ELE))
          END DO
        END DO
        
        MAXMAX=-1.0E+20
        MINMIN=+1.0E+20
        DO ELE=1,TOTELE
          MAXMAX=MAX(MAXMAX,MEAN_DGT(ELE))
          MINMIN=MIN(MINMIN,MEAN_DGT(ELE))
        END DO
        PRINT *,'AFTER ADJUSTING THE MEANS MEAN MIN,MAX:',MINMIN,MAXMAX
        
        END SUBROUTINE BOUND_DG_MEANS
        
        

        SUBROUTINE MAK_MASS_FROM_LUMPED(MLELE,MATRIX,TOTELE,NLOC, &
          NCOLMELE,FINDRMELE,COLMELE,CENTRMELE)
! Calculate element centre distributed mass matrix...
        INTEGER TOTELE,NLOC,NCOLMELE
        INTEGER FINDRMELE(TOTELE+1),COLMELE(NCOLMELE),CENTRMELE(TOTELE)
        REAL MLELE(TOTELE*NLOC),MATRIX(NCOLMELE)
! Local variables...
        REAL RSUM,MASS_OFF_D
        INTEGER ELE,ILOC,COUNT,ELE2
        REAL, ALLOCATABLE, DIMENSION(:)::MASS_LUMPED
          
        ALLOCATE(MASS_LUMPED(TOTELE))
        
        MASS_LUMPED=0.0
        DO ELE=1,TOTELE
          DO ILOC=1,NLOC
            MASS_LUMPED(ELE)=MASS_LUMPED(ELE)+MLELE((ELE-1)*NLOC+ILOC)
          END DO
        END DO   
               
        MATRIX=0.0
        DO ELE=1,TOTELE
          RSUM=0.0
          DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
             ELE2=COLMELE(COUNT)
             MASS_OFF_D=min(MASS_LUMPED(ELE),MASS_LUMPED(ELE2))
             MATRIX(COUNT)=MASS_OFF_D
             RSUM=RSUM+MASS_OFF_D
          END DO
          MATRIX(CENTRMELE(ELE))=(RSUM-MASS_LUMPED(ELE))
        END DO
        
        END SUBROUTINE MAK_MASS_FROM_LUMPED




         SUBROUTINE BOUND_FIX(MATRIX2,TOTELE,NLOC, &
          MEAN_DGT,MEAN_DGT_MIN,MEAN_DGT_MAX, &
          NCOLMELE,FINDRMELE,COLMELE,CENTRMELE, &
          NO_IN_LIMIT_REG,NIT_LIMIT_MEAN,halo_tag)
! Obtain bounded soln fo MEAN_DGT.
! If NO_IN_LIMIT_REG Recalculate matrix to to stop adjustment from going into 
!regions with bounds exceeded or about to be...
! halo_tag is halo for velocity which is converted to elements. 
! NIT_LIMIT_MEAN=no of iterations...
        INTEGER TOTELE,NLOC,NCOLMELE
        INTEGER FINDRMELE(TOTELE+1),COLMELE(NCOLMELE),CENTRMELE(TOTELE)
        REAL MEAN_DGT(TOTELE),MEAN_DGT_MIN(TOTELE),MEAN_DGT_MAX(TOTELE)
        REAL MATRIX2(NCOLMELE)
        LOGICAL NO_IN_LIMIT_REG
        INTEGER NIT_LIMIT_MEAN,halo_tag
! Local variables...
        INTEGER ELE,ELE2,COUNT,IT_LIMIT,ICENT,JCENT
        REAL RMAT
        REAL, ALLOCATABLE, DIMENSION(:)::MATRIX,MATRIX_LUMPED
        REAL, ALLOCATABLE, DIMENSION(:)::MEAN_DGT_dev,MEAN_DGT_alt,PREV_IT_MEAN_DGT
        LOGICAL, ALLOCATABLE, DIMENSION(:)::LUMP_ROW
          
        ALLOCATE(MATRIX(NCOLMELE))
        ALLOCATE(MATRIX_LUMPED(TOTELE))
        
        ALLOCATE(MEAN_DGT_dev(TOTELE))
        ALLOCATE(MEAN_DGT_alt(TOTELE))
        ALLOCATE(PREV_IT_MEAN_DGT(TOTELE))
        
        ALLOCATE(LUMP_ROW(TOTELE))
         
           MATRIX=MATRIX2
! Calculate lumped matrix...           
           DO ELE=1,TOTELE
             MATRIX_LUMPED(ELE)=0.0
             DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
                MATRIX_LUMPED(ELE)=MATRIX_LUMPED(ELE)+MATRIX(COUNT)
             END DO
           END DO
         
            DO IT_LIMIT=1,NIT_LIMIT_MEAN
            
              DO ELE=1,TOTELE
                  if (MEAN_DGT(ELE) .GT. MEAN_DGT_max(ELE)) then
                    MEAN_DGT_dev(ELE) = MEAN_DGT(ELE) - MEAN_DGT_max(ELE)
                  else if (MEAN_DGT(ELE) .LT. MEAN_DGT_min(ELE)) then
                    MEAN_DGT_dev(ELE) = MEAN_DGT(ELE) - MEAN_DGT_min(ELE)
                  else
                    MEAN_DGT_dev(ELE) = 0.0
                  end if 
              END DO
              
! Limit the value of MEAN_DGT...
!Now spread this deviation T_dev around the surrounding nodes in an element
! M_L T_alt= M T_dev:
              MATRIX=MATRIX2
              PREV_IT_MEAN_DGT=MEAN_DGT
              DO ELE=1,TOTELE
                MEAN_DGT_alt(ELE)=0.0
                DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
                  ELE2=COLMELE(COUNT)
                  MEAN_DGT_alt(ELE)=MEAN_DGT_alt(ELE)+MATRIX(COUNT)*MEAN_DGT_DEV(ELE2)
                END DO
                MEAN_DGT_alt(ELE)=MEAN_DGT_alt(ELE)/MATRIX_LUMPED(ELE)
              END DO 
! T_new=T_new-T_dev+T_alt:
              MEAN_DGT=MEAN_DGT-MEAN_DGT_DEV+MEAN_DGT_alt
! get halo's -ver halo_tag for elements...
              call flcomms_update(-halo_tag, MEAN_DGT, 1, 1, 0)
              
              
         IF(NO_IN_LIMIT_REG) THEN
! See if there are any MEAN_DGT that are worse 
! if so make the matrix diagonal for that row/coln.
! If MEAN_DGT_DEV has the same sign between neigbouring elements dont exchange mass.
              DO ELE=1,TOTELE
                LUMP_ROW(ELE)=.FALSE.
                IF(MEAN_DGT(ELE) .GT. MEAN_DGT_max(ELE)) THEN
                  IF(MEAN_DGT(ELE) .GT. PREV_IT_MEAN_DGT(ELE)) LUMP_ROW(ELE)=.TRUE.
                ELSE IF(MEAN_DGT(ELE) .LT. MEAN_DGT_min(ELE)) THEN
                  IF(MEAN_DGT(ELE) .LT. PREV_IT_MEAN_DGT(ELE)) LUMP_ROW(ELE)=.TRUE.
                ENDIF
              END DO 
              
              DO ELE=1,TOTELE
                ICENT=CENTRMELE(ELE)
                DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
                  ELE2=COLMELE(COUNT)
                  IF(ELE.NE.ELE2) THEN
                    IF(LUMP_ROW(ELE)) THEN
                      RMAT=MATRIX(COUNT)
                      MATRIX(COUNT)=0.0
                      MATRIX(ICENT)=MATRIX(ICENT)+RMAT
                    ENDIF
                    IF(LUMP_ROW(ELE2)) THEN
                      JCENT=CENTRMELE(ELE2)
                      RMAT=MATRIX(COUNT)
                      MATRIX(COUNT)=0.0
                      MATRIX(JCENT)=MATRIX(JCENT)+RMAT
                    ENDIF
                  ENDIF
                END DO
              END DO 
! Now spread this deviation T_dev around the surrounding nodes in an element
! with the modified matrix...
! M_L T_alt= M T_dev:
              MEAN_DGT=PREV_IT_MEAN_DGT
              DO ELE=1,TOTELE
                MEAN_DGT_alt(ELE)=0.0
                DO COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
                  ELE2=COLMELE(COUNT)
                  MEAN_DGT_alt(ELE)=MEAN_DGT_alt(ELE)+MATRIX(COUNT)*MEAN_DGT_DEV(ELE2)
                END DO
                MEAN_DGT_alt(ELE)=MEAN_DGT_alt(ELE)/MATRIX_LUMPED(ELE)
              END DO 
! T_new=T_new-T_dev+T_alt:
              MEAN_DGT=MEAN_DGT-MEAN_DGT_DEV+MEAN_DGT_alt
! get halo's -ver halo_tag for elements...
              call flcomms_update(-halo_tag, MEAN_DGT, 1, 1, 0)
           
         ENDIF
              
! Next limiter iteration...
            END DO
            END SUBROUTINE BOUND_FIX

      



   SUBROUTINE MAKE_NO_INFLOW(TOTELE,NLOC,XONDGL,SNLOC,XNONOD,&
     &                DGNU,DGNV,DGNW,X,Y,Z)
! This sub adjusts the DG velocities so 
! it has no component noraml to the boudaries...
     INTEGER TOTELE,NLOC,XONDGL(TOTELE*NLOC),SNLOC,XNONOD
     REAL DGNU(TOTELE*NLOC),DGNV(TOTELE*NLOC),DGNW(TOTELE*NLOC)
     REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
     INTEGER, ALLOCATABLE, DIMENSION(:)::COLELE4
       
     ALLOCATE(COLELE4(TOTELE*5))
     CALL GETFINELE4(TOTELE,NLOC,SNLOC,XONDGL, COLELE4,XNONOD)

     CALL MAKE_NO_INFLOW2(TOTELE,NLOC,XONDGL,SNLOC,XNONOD,&
     &                COLELE4,&
     &                DGNU,DGNV,DGNW,X,Y,Z)

   END SUBROUTINE MAKE_NO_INFLOW

   SUBROUTINE MAKE_NO_INFLOW2(TOTELE,NLOC,XONDGL,SNLOC,XNONOD,&
     &                COLELE4,&
     &                DGNU,DGNV,DGNW,X,Y,Z)
! This sub adjusts the DG velocities so 
! it has no component noraml to the boudaries...
     INTEGER TOTELE,NLOC,XONDGL(TOTELE*NLOC),SNLOC,XNONOD
     INTEGER COLELE4(5*TOTELE)
     REAL DGNU(TOTELE*NLOC),DGNV(TOTELE*NLOC),DGNW(TOTELE*NLOC)
     REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
! local variables...
     INTEGER IFACE,NFACE
     PARAMETER(NFACE=4)
     INTEGER LOCLIST(NFACE,3)
     INTEGER ELE,ELE2,SILOC,ILOC,INOD,NOD1,NOD2,NOD3,IDG,globi
     REAL VU,VV,VW,RN
     REAL DIRX,DIRY,DIRZ
     REAL RTREST11,RTREST12,RTREST13
     REAL RTREST21,RTREST22,RTREST23
     REAL RTREST31,RTREST32,RTREST33
! The local cords are in anti-clockwise order
     IFACE=1
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=3
     IFACE=2
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=4
     IFACE=3
     LOCLIST(IFACE,1)=3
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=4
     IFACE=4
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=3
     LOCLIST(IFACE,3)=4
       
     do  ELE=1,TOTELE! Was loop 10
       
        do  IFACE=1,NFACE! Was loop 441
         
! NB colele4 is ordered in terms of faces.
           ELE2=COLELE4((ELE-1)*5+IFACE)
           IF(ELE2.EQ.0) THEN
              SILOC=1
              ILOC=LOCLIST(IFACE,SILOC)
              NOD1=XONDGL((ELE-1)*NLOC+ILOC)
              SILOC=2
              ILOC=LOCLIST(IFACE,SILOC)
              NOD2=XONDGL((ELE-1)*NLOC+ILOC)
              SILOC=3
              ILOC=LOCLIST(IFACE,SILOC)
              NOD3=XONDGL((ELE-1)*NLOC+ILOC)

              CALL XPROD(DIRX,DIRY,DIRZ, &
     &              X(NOD1)-X(NOD2),Y(NOD1)-Y(NOD2),Z(NOD1)-Z(NOD2),&
     &              X(NOD1)-X(NOD3),Y(NOD1)-Y(NOD3),Z(NOD1)-Z(NOD3))
              RN=SQRT(DIRX**2+DIRY**2+DIRZ**2)
              DIRX=DIRX/RN
              DIRY=DIRY/RN
              DIRZ=DIRZ/RN
  
              RTREST11=1.0 - DIRX*DIRX
              RTREST12=    - DIRX*DIRY
              RTREST13=    - DIRX*DIRZ
              
              RTREST21=    - DIRY*DIRX
              RTREST22=1.0 - DIRY*DIRY
              RTREST23=    - DIRY*DIRZ
              
              RTREST31=    - DIRZ*DIRX
              RTREST32=    - DIRZ*DIRY
              RTREST33=1.0 - DIRZ*DIRZ
              do SILOC=1,SNLOC
                 ILOC=LOCLIST(IFACE,SILOC)
                 globi=XONDGL((ELE-1)*NLOC+ILOC)
                 IDG=(ELE-1)*NLOC+ILOC

                 VU=DGNU(IDG)
                 VV=DGNV(IDG)
                 VW=DGNW(IDG)
! Redfine the velocity without a normal component...
                 DGNU(IDG)=RTREST11*VU+RTREST12*VV+RTREST13*VW
                 DGNV(IDG)=RTREST21*VU+RTREST22*VV+RTREST23*VW
                 DGNW(IDG)=RTREST31*VU+RTREST32*VV+RTREST33*VW

              END DO
           ENDIF
        END DO
     END DO

   END SUBROUTINE MAKE_NO_INFLOW2




   SUBROUTINE HART3DSGS_DG(DGT,DGTNEW,&
       &     DGNU,DGNV,DGNW,UG,VG,WG,&
       &     SOURCT,X,Y,Z,D0,&
       &     SUBSOUT,ISUBSOU,&
       &     VECT_SGSADD, &
       &     VECT,&
       &     BIGMELE,MLELE,ML,TOTELE, &
       &     FINDRMELE,COLMELE,NCOLMELE,NONODS,CENTRMELE, &
       &     N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
       &     DISOPT2,DISOPN,DT,THETA,BETA,LUMP,PREOPT,&
       &     ABSLUM,SLUMP,&
       &     NDGLNO,&
       &     SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
       &     DENPT,&
       &     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
       &     ABSORB,VERSIO,ISPHERE,&
       &     nnodp,para,halo_tag,&
! for the DG...
       &     NBIGMELE,NOBCFSimp,BCT2FSimp,&
       &     MASSDGELE,&
       &     SUF_TOP_OCEAN,NO_NON_LIN,ISTAB,STAB_DIF_MAT, velocity_option_path)
      
!     This subroutine discretises a field equation using Galerkin Least squares.
!     The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
! If ISTAB=1 construct stabilisation matrix STAB_DIF_MAT else =0
!     
! All methods use a balancing diffusion term and THETA time stepping. 
! The following are absolute values of DISOPT...
! Galerkin (146 is defaulted to if non of the below are used): 
! 8 then Galerkin by use a non-linear boundedness correction after solving. 
! 158 non-linear DG method but diffusion only in the horizontal for oceans. 
! 157 isotropic non-linear DG method
! 156 same as 155 but reduced magnitude by a factor of 8.
! 155 same as 147 but diffusion only in the horizontal.

! 154 same as 147 and no added viscocity.
! 153 same as 152 and no added viscocity.
! 152 same as 147 but reduced magnitude by a factor of 8.
! 151 is LES applied to the global system. 
! 150 is LES applied to the local Residual system (involving pertation soln) only RECOMMENDED
! 149 is LES applied to the local system only RECOMMENDED
! 148 is LES applied to the global system only 
! 147 is LES applied to every thing. 
! 145 PG applied to just the global part and using global velocity.
! 144 PG applied to just the global part and using the full velocity.
      
!     if DISOPN.ne.0 then set advection to zero and treat advection 
!     using the high res method.
!     IF LUMP then lump the mass matrix else dont. 
!     BETA controls the conservative discretisation 
!     BETA=0. is the usual advection form, BETA=1. is divergence form. 
!     BETA=0.5 is an average of the divergence and advection form.(recommended).
!     ABSORB =absorbtion of energy etc. 
!     MUPTXX,MUPTYY =components of diffusivities. 
!     ML2MXX,ML2MXY etc contain the length scales for LES. 
!     NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
!     UG,VG,WG are the grid velocities. 
!     NDGLNO=element pter for unknowns. 
!     SONDGL=element pter for materials and sources. 
!     VONDGL=element pter for advection NU,UG.
!     XONDGL=element pter for coordinates(x,y,z)
!     If INMAT then return matricie(s). 
!     R0=the minimum distance from ocean bed to centre of earth.
!     D0=H/r0 where H is the height of ocean above R0 (H not used here)
!     For dimensional run D0=0.(R0 is not used used here)
!     NB This sub is only used for field equations in ocean modelling. 
!     ----------------------------------------------------------------------
!     If we are solving for another variable like temperature 
!     or chemical species then NU,NV,NW will be the velocities 
!     and U=TEMPERATURE OR CHEMICAL SPECIES. 
!     ----------------------------------------------------------------------
!     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
     REAL RWIND
     LOGICAL ADVECT_MEAN,SUF_CONSERV
!     PARAMETER(RWIND=0.5,ADVECT_MEAN=.FALSE.)
     PARAMETER(RWIND=0.5,ADVECT_MEAN=.TRUE.)
    ! SUF_CONSERV controls conservation of between element fluxes.
     INTEGER SNONOD,VNONOD,XNONOD,NCT,NCOLM,NLOC,NGI
     INTEGER MLOC,DISOPT2,DISOPT,DISOPN,DISCRA
     REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
     REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
     REAL DENPT(SNONOD)
     REAL ABSORB(SNONOD)
     REAL  DT,THETA,BETA
     INTEGER TOTELE,NONODS,FREDOP
!     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
!     If PHI=1.0 then bakward Euler is used in NS equns.
!     Similarly for the temperature equation except width variable THETA.
!     have gravity in -ve y direction.
     REAL DGT(TOTELE*NLOC)
     REAL DGTNEW(TOTELE*NLOC)
     REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
     REAL DGNU(TOTELE*NLOC),DGNV(TOTELE*NLOC),DGNW(TOTELE*NLOC)
      
! This is what we add into the SGS model - a discretised source.      
     REAL VECT_SGSADD(TOTELE*NLOC)
     REAL VECT(TOTELE*NLOC)
     INTEGER NBIGMELE
     REAL BIGMELE(NBIGMELE)
     REAL MASSDGELE(TOTELE,NLOC,NLOC)
     LOGICAL COGRAX
     INTEGER NOBCFSimp
     INTEGER BCT2FSimp(NOBCFSimp)
     REAL MLELE(TOTELE*NLOC),ML(NONODS)
     REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
     REAL SOURCT(SNONOD)
     INTEGER ISUBSOU
     REAL SUBSOUT(ISUBSOU*TOTELE*NLOC)
     INTEGER NDGLNO(TOTELE*NLOC)
     INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
     INTEGER XONDGL(TOTELE*NLOC)
     INTEGER NCOLMELE
     INTEGER FINDRMELE(TOTELE+1),COLMELE(NCOLMELE)
     INTEGER CENTRMELE(TOTELE)
     LOGICAL SUF_TOP_OCEAN,NO_NON_LIN

     REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
     REAL WEIGHT(NGI)
     LOGICAL LUMP,UPWIND,OPWIND,XTPETV
     INTEGER PREOPT
! IF PREOPT=2 then treat every boundary like a free surface - need 
! to do this for open boundaries. 
     LOGICAL ABSLUM,SLUMP
     LOGICAL INMAT
     INTEGER VERSIO
     INTEGER nnodp,para,halo_tag
     
     REAL MASS,MATRIX,SOUMAT,VLMGI,VLKGI,VABS
     REAL VLM
     REAL VF1M,VF2M,VF2,VF3M
     REAL THETA1,THETA2,THETA3,AHAT1,ALPHA,ALPHA2,ALPHA3,THTWDT,INDT,TWOTHI
     REAL PE,HOVERQ,PWEIGI
     REAL GALERK,LEASQR
     REAL UD(NGI),VD(NGI),WD(NGI)
     REAL UDDX(NGI),VDDY(NGI),WDDZ(NGI)
     REAL UDDX2(NGI),VDDY2(NGI),WDDZ2(NGI)
     REAL UGDX(NGI),VGDY(NGI),WGDZ(NGI)
     REAL MUGIXX,MUGIXY,MUGIXZ
     REAL MUGIYY,MUGIYZ,MUGIZZ
     REAL L2GIXX,L2GIXY,L2GIXZ
     REAL L2GIYY,L2GIYZ,L2GIZZ
     REAL DENGI(NGI),DENGI2(NGI),L1GI(NGI)
     REAL ABDGI(NGI)
     REAL AGI,BGI,CGI,DGI
     REAL EGI,FGI,GGI,HGI,KGI
     REAL DETJ,DETWEI(NGI)
     REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
     REAL GIVOLX,GIVOLY,GIVOLZ
!     HX,HY-characteristic length scales in x,y directions.
     REAL HXGI,HYGI,HZGI
     REAL AHAT(NGI),GAMMA(NGI),GAMMAG(NGI),MASCON(NGI)
     REAL DIRX(NGI),DIRY(NGI),DIRZ(NGI)
     REAL BECTWE(NGI)

     INTEGER ELE,ILOC,JLOC,L
     INTEGER GI,COUNT,I,GLOBI,GLOBJ,GLOBIS,GLOBJS
     INTEGER POS,LOWER,UPPER,PGLOBJ

     REAL D0,DD0,INDD0
     LOGICAL NODIM
     LOGICAL GALWIN
     LOGICAL CTWIND,LES
     REAL ABLUMP,RSLUMP
     REAL R,RUD,RVD,RWD,R1,R2,RR,RMXD,RMXMXX,RMXMYY,RMXMZZ,RMXABS
     REAL::RNN=0.0
     REAL RMXSX,RMXSY,RMXSZ,RMID,RMIMXX,RMIMYY,RMIMZZ,RMIABS
     REAL RMISX,RMISY,RMISZ,RLUMP,RSYM,RCONSI,RGAMMA
     INTEGER IGLS,IGLV,IGLX,JGLX,INUM
     INTEGER ISPHERE,ICENT,ISTAB
     REAL STAB_DIF_MAT(TOTELE*ISTAB,NLOC,NLOC)
     character(len=*), intent(in):: velocity_option_path
     REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
     REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
     REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
     REAL MYZDTW(NGI),MZZDTW(NGI)
     REAL MUGIXX2,MUGIXY2,MUGIXZ2
     REAL MUGIYY2,MUGIYZ2,MUGIZZ2
     REAL MUPTXXGI,MUPTXYGI,MUPTXZGI,MUPTYYGI,MUPTYZGI,MUPTZZGI
     REAL XD(NGI),YD(NGI),ZD(NGI)
     REAL VLND,VLN2,VLND2,VLN3,VLND3
     REAL TENSXX,TENSXY,TENSXZ,TENSYY,TENSYZ,TENSZZ
     REAL L2GIXX2,L2GIXY2,L2GIXZ2,L2GIYY2,L2GIYZ2,L2GIZZ2
!  for the extra terms in the adjoint model, local variables....
     REAL VLND_UX_ADJ,VLND_UY_ADJ,VLND_UZ_ADJ,    &
         &       VLND_VX_ADJ,VLND_VY_ADJ,VLND_VZ_ADJ, &
         &       VLND_WX_ADJ,VLND_WY_ADJ,VLND_WZ_ADJ  
     REAL HINX,HINY,HINZ,VLKTEN,NLN,VLMAB
     REAL XNDMUJLOC,YNDMUJLOC,ZNDMUJLOC,XNDMUILOC,YNDMUILOC,ZNDMUILOC
     REAL TEN,TENAB,ABTEN
     REAL VLN,ADVILOC,ADVJLOC,RDIFF
     REAL UATTHETA,VATTHETA,WATTHETA
     LOGICAL BALDIF,BALGLS
     REAL L1(NGI),L2(NGI),L3(NGI),L4(NGI),PGWEIG(NGI) 
     REAL A11(NGI),A12(NGI),A13(NGI)
     REAL A21(NGI),A22(NGI),A23(NGI)
     REAL A31(NGI),A32(NGI),A33(NGI) 
     REAL GAMB(NGI),GAMBAL(NGI)
     REAL ABGI(NGI),SOUGI(NGI),DDETWE(NGI)
    
     REAL MASSMM(NLOC,NLOC),MASSINVMM(NLOC,NLOC)
     REAL MASSLUMMM(NLOC,NLOC)
       
     REAL TOTSOUT(NLOC)
     INTEGER LESNOD,JNODSUB
     REAL RINMAT,RHSMAT,SRHSMAT
     REAL VLNXN,VLNYN,VLNZN
     REAL VLNNX,VLNNY,VLNNZ
     REAL VLMXN,VLMYN,VLMZN
     REAL C1TCONU,C2TCONV,C3TCONW
! Local variables...
     INTEGER SNLOC,SNGI,SNCLOC
     PARAMETER(SNLOC=3,SNCLOC=6)
     REAL XSL(SNCLOC),YSL(SNCLOC),ZSL(SNCLOC)
     REAL SGS_UPWIND_FAC(NLOC,NGI),SGS_UPWIND_FAC_DIF(NLOC,NGI)
     REAL TOLDGI(NGI),TOLDGIX(NGI),TOLDGIY(NGI),TOLDGIZ(NGI)
     REAL TNEWGI(NGI),TNEWGIX(NGI),TNEWGIY(NGI),TNEWGIZ(NGI)
     REAL TTHETAGI(NGI),TTHETAGIX(NGI),TTHETAGIY(NGI),TTHETAGIZ(NGI)
     REAL USTAR(NGI),VSTAR(NGI),WSTAR(NGI)
       
     REAL MASSDGELE5(NLOC+1,NLOC+1)
     REAL MATDIFX5(NLOC+1,NLOC+NLOC),MATDIFY5(NLOC+1,NLOC+NLOC),MATDIFZ5(NLOC+1,NLOC+NLOC)
     REAL MMATDIFX5(NLOC+1,NLOC+NLOC),MMATDIFY5(NLOC+1,NLOC+NLOC),MMATDIFZ5(NLOC+1,NLOC+NLOC)
     REAL MASSMM5(NLOC+1,NLOC+1),MASSINVMM5(NLOC+1,NLOC+1)
     REAL SUFMATX(NLOC,NLOC+1),SUFMATY(NLOC,NLOC+1),SUFMATZ(NLOC,NLOC+1)
     REAL SMMATDIFX5(NLOC,NLOC+NLOC),SMMATDIFY5(NLOC,NLOC+NLOC),SMMATDIFZ5(NLOC,NLOC+NLOC)
     REAL SMMATDIF5(NLOC,NLOC+NLOC)
     REAL FAMAT(NLOC,5*NLOC),AMAT(NLOC,NLOC)
     INTEGER IFACE,NFACE
     PARAMETER(NFACE=4)
     INTEGER LOCLIST(NFACE,3),SILOC2ILOC(SNLOC),ILOC_OTHER_SIDE(SNLOC)
     INTEGER ELE2_LOC_NODS(NLOC)
     INTEGER SILOC,SJLOC,SKLOC,IGL,COL,INOD,JNOD,SGI,KLOC
     INTEGER BIG_ILOC,BIG_JLOC,ELE2,IJDISP,ID,JD,IXNOD,NODDG
     INTEGER IDGNOD,NOD,II,ICOUNT,PGLOBI,ILOC2,COUNTELE,JJLOC5,IILOC5
     INTEGER::IDG,IDG2,ILOCOTH,ILOCELE=-1,INOD2,JDG2,JDG,JELE,JLOC2,JLOCELE=-1
     REAL RADI,RADJ,RADP,RAD,RNU,RNV,RNW,SAREA,VOL
     REAL NORMX,NORMY,NORMZ,VLNDOTQ,allVLNDOTQ,STABINEL,UPWIN,TDIVU
     REAL VLOMEGA,OMEGAGI,STABOME
     REAL A11MAT,A12MAT,A21MAT,A22MAT,O11MAT,O12MAT,O22MAT,O21MAT
     REAL RO11MAT,RO12MAT,RO22MAT,RO21MAT
     REAL RNUSUBGI,RNVSUBGI,RNWSUBGI
     REAL AVE_SNORMXN,AVE_SNORMYN,AVE_SNORMZN,ONE_OR_ZERO,ROCEAN
     REAL LOCPNORMX,LOCPNORMY,LOCPNORMZ
     REAL XC,YC,ZC,RN,sum,stVLKTEN,stTEN
     REAL ALPHA11,ALPHA12,ALPHA21,ALPHA22,SUMP,SUMU
     REAL RLES,RHAVEVIS,VLNDOTQ_IN,VLNDOTQ_OUT,VLNDOTQ_ALL
     REAL SMEANXX,SMEANXY,SMEANXZ,SMEANYY,SMEANYZ,SMEANZZ
     REAL THETA_FAMAT,VLM_NORX,VLM_NORY,VLM_NORZ
     REAL MNNX,MNNY,MNNZ,NN
     LOGICAL DMATLUM,SUF_TOP,STREAM_UPWIND,BNDSUF
     LOGICAL SGSUPWIND,OPTWIND,ISOTROPIC
     REAL ALPHA1,SIG_THETA,HAT_SIG,EXP_SIG,ALPHA_OPT
     REAL A_DOT_GRADT,TTHETADIVU,RESID,SIGCHAN,COEF,DIFSTAB
     REAL::RANISO,DIST2,HHORIZ=-1,HDIST
     ! Real function...
     REAL EXPTOL
! Local coords...
     INTEGER NCLOC,SMLOC,NBLOCK
     REAL, ALLOCATABLE, DIMENSION(:,:)::NC
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCLX
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCLY
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCLZ
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCX
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCY
     REAL, ALLOCATABLE, DIMENSION(:,:)::NCZ
     REAL, ALLOCATABLE, DIMENSION(:)::VEDUDX
     REAL, ALLOCATABLE, DIMENSION(:)::VEDUDY
     REAL, ALLOCATABLE, DIMENSION(:)::VEDUDZ
     REAL, ALLOCATABLE, DIMENSION(:)::VEDVDX
     REAL, ALLOCATABLE, DIMENSION(:)::VEDVDY
     REAL, ALLOCATABLE, DIMENSION(:)::VEDVDZ
     REAL, ALLOCATABLE, DIMENSION(:)::VEDWDX
     REAL, ALLOCATABLE, DIMENSION(:)::VEDWDY
     REAL, ALLOCATABLE, DIMENSION(:)::VEDWDZ
    
     REAL, ALLOCATABLE, DIMENSION(:)::HIMATX
     REAL, ALLOCATABLE, DIMENSION(:)::HIMATY
     REAL, ALLOCATABLE, DIMENSION(:)::HIMATZ
     REAL, ALLOCATABLE, DIMENSION(:)::HIMASS
 
     REAL, ALLOCATABLE, DIMENSION(:)::VECU
     REAL, ALLOCATABLE, DIMENSION(:)::VECV
     REAL, ALLOCATABLE, DIMENSION(:)::VECW
     
     REAL, ALLOCATABLE, DIMENSION(:)::VECYT

     REAL, ALLOCATABLE, DIMENSION(:,:,:)::EV
     INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
     INTEGER, ALLOCATABLE, DIMENSION(:)::COLELE4
         
     REAL, ALLOCATABLE, DIMENSION(:)::MASSBIGM
     REAL, ALLOCATABLE, DIMENSION(:)::MASSVEC
         
     REAL, ALLOCATABLE, DIMENSION(:,:)::SNSUB
     INTEGER, ALLOCATABLE, DIMENSION(:)::BNDCON
     REAL, ALLOCATABLE, DIMENSION(:)::VOLELE
     REAL, ALLOCATABLE, DIMENSION(:)::MEANXX,MEANXY,MEANXZ,MEANYY,MEANYZ,MEANZZ
     REAL, ALLOCATABLE, DIMENSION(:,:,:)::MATDIFX,MATDIFY,MATDIFZ,MATVLK
! surface memory...
         REAL, ALLOCATABLE, DIMENSION(:,:)::SN,SNLX,SNLY
         REAL, ALLOCATABLE, DIMENSION(:,:)::SNX,SNY,SNZ
         REAL, ALLOCATABLE, DIMENSION(:)::SWEIGH
         REAL, ALLOCATABLE, DIMENSION(:,:)::SNC,SNCLX,SNCLY
         REAL, ALLOCATABLE, DIMENSION(:,:)::SNCX,SNCY,SNCZ 
         REAL, ALLOCATABLE, DIMENSION(:)::SUD,SVD,SWD
         REAL, ALLOCATABLE, DIMENSION(:)::IN_SUD,IN_SVD,IN_SWD
         REAL, ALLOCATABLE, DIMENSION(:)::OUT_SUD,OUT_SVD,OUT_SWD 
         REAL, ALLOCATABLE, DIMENSION(:)::SUDGLOB,SVDGLOB,SWDGLOB 
         REAL, ALLOCATABLE, DIMENSION(:)::NDOTQ,NDOTQGLOB
         REAL, ALLOCATABLE, DIMENSION(:)::IN_NDOTQ,OUT_NDOTQ
         REAL, ALLOCATABLE, DIMENSION(:)::SL1, SL2, SL3, SL4
         REAL, ALLOCATABLE, DIMENSION(:)::SDETWE 
         REAL, ALLOCATABLE, DIMENSION(:)::SNORMXN,SNORMYN,SNORMZN

     ewrite(1, *) "in subroutine hart3dsgs_dg"
     ewrite(2, *) "isphere=",isphere

! BNDCON marks the boundary condition nodes...
     ALLOCATE(BNDCON(NONODS))

     BNDCON(1:NONODS)=0
     do II=1,NOBCFSimp
        GLOBI=BCT2FSimp(II)
        BNDCON(GLOBI)=1
     END DO
 
     COGRAX=ISPHERE.NE.2

     SUF_CONSERV=(ABS(BETA).GT.1.E-3)
      
     DISOPT=DISOPT2
     IF(DISOPT2.EQ.47) DISOPT=147
     IF(DISOPT2.EQ.48) DISOPT=148
     
     NCLOC=4
     IF(ISPHERE.EQ.2) NCLOC=10

     ALLOCATE(NC(NCLOC,NGI))
     ALLOCATE(NCLX(NCLOC,NGI))
     ALLOCATE(NCLY(NCLOC,NGI))
     ALLOCATE(NCLZ(NCLOC,NGI))
     ALLOCATE(NCX(NCLOC,NGI))
     ALLOCATE(NCY(NCLOC,NGI))
     ALLOCATE(NCZ(NCLOC,NGI))
       
     CALL PMINMX(MUPTXX,NONODS,'MUPTXX*******  ')
     CALL PMINMX(MUPTXY,NONODS,'MUPTXY*******  ')
     CALL PMINMX(MUPTXZ,NONODS,'MUPTXZ*******  ')
     CALL PMINMX(MUPTYY,NONODS,'MUPTYY*******  ') 
     CALL PMINMX(MUPTYZ,NONODS,'MUPTYZ*******  ')
     CALL PMINMX(MUPTZZ,NONODS,'MUPTZZ*******  ')
       
     IF(ISPHERE.NE.0) THEN
        ASSERT(ncloc==10)
        call QuadTetShapes(nc, nclx, ncly, nclz, ncloc, ngi)
     ELSE
        NC = N
        NCLX = NLX
        NCLY = NLY
        NCLZ = NLZ
     ENDIF

     ML=0.0
     MLELE=0.0
     do  I=1,TOTELE*NLOC
        VECT(I)=0.
     end do

     RLUMP=1.
     IF(LUMP) RLUMP=1.
     RCONSI=1.-RLUMP

     ABLUMP=1.
     RSLUMP=1.
     IF(ABSLUM) ABLUMP=1.
     IF(SLUMP)  RSLUMP=1.
! ALPHA contains the inner element stabilisation that ensures 
! a non-singular system...

     ALPHA=0.01 


     ewrite(3,*) "InnerElement: Filter section"
     if (have_option(trim(velocity_option_path) // "/prognostic/&
          &spatial_discretisation/inner_element")) then
       ewrite(3,*) "InnerElement: Filter section - detected model"
       if (have_option(trim(velocity_option_path) // "/prognostic/&
            &spatial_discretisation/inner_element/&
            &use_filter")) then
          ewrite(3,*) "InnerElement: SGS filtering turned on"
          call get_option(trim(velocity_option_path) // "/prognostic/&
            &spatial_discretisation/inner_element/&
            &use_filter/strength", alpha, default=0.0) 
       else
          ewrite(3,*) "InnerElement: SGS filtering turned off"
          alpha=0.0
       end if
     end if
     
    ewrite(1,*) "InnerElement: Filter section, alpha[",alpha,']'


! alpha=1.0 is good or DMATLUM=.TRUE. this way it works without 
! transport
     DMATLUM=.false.

! alpha2 control works for void but others do not.
     ALPHA2=0.0
      
     ALPHA3=0.0
     ROCEAN=0.0
     IF(SUF_TOP_OCEAN) ROCEAN=1.0
! DISOPT=151 is LES applied to the global system. 
! DISOPT=150 is LES applied to the local Residual system (involving pertation soln) only RECOMMENDED
! DISOPT=149 is LES applied to the local system only RECOMMENDED
! 150 and 149 produce the same results.

     RLES=1.0
     IF((DISOPT.EQ.152).OR.(DISOPT.EQ.153)) RLES=0.125
     RHAVEVIS=1.0
     IF((DISOPT.EQ.153).OR.(DISOPT.EQ.154)) RHAVEVIS=0.0

     NBLOCK=(NLOC*3)**2
     
     SGSUPWIND=(DISOPT.EQ.157).OR.(DISOPT.EQ.158)
     OPTWIND=.false.
     IF(NO_NON_LIN) SGSUPWIND=.FALSE.
     ISOTROPIC=(DISOPT.EQ.157)
     RANISO=1.0
     IF(ISOTROPIC) RANISO=0.0

        IF(NGI.EQ.4) THEN
          SNGI=3
        ELSE
          SNGI=7
        ENDIF

       ALLOCATE( SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI) )
       ALLOCATE( SNX(SNLOC,SNGI),SNY(SNLOC,SNGI),SNZ(SNLOC,SNGI) )
       ALLOCATE( SWEIGH(SNGI) )
       ALLOCATE( SNC(SNCLOC,SNGI),SNCLX(SNCLOC,SNGI),SNCLY(SNCLOC,SNGI) )
       ALLOCATE( SNCX(SNCLOC,SNGI),SNCY(SNCLOC,SNGI),SNCZ(SNCLOC,SNGI) )
       ALLOCATE( SUD(SNGI),SVD(SNGI),SWD(SNGI) )
       ALLOCATE( IN_SUD(SNGI),IN_SVD(SNGI),IN_SWD(SNGI) )
       ALLOCATE( OUT_SUD(SNGI),OUT_SVD(SNGI),OUT_SWD(SNGI) )
       ALLOCATE( SUDGLOB(SNGI),SVDGLOB(SNGI),SWDGLOB(SNGI) )
       ALLOCATE( NDOTQ(SNGI),NDOTQGLOB(SNGI) )  
       ALLOCATE( IN_NDOTQ(SNGI),OUT_NDOTQ(SNGI) )
       ALLOCATE( SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI) )
       ALLOCATE( SDETWE(SNGI) )
       ALLOCATE( SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI) )

! Calculate surface shape functions SN, SM etc...
     CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
     CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
     &             SNLOC,SNGI, &
     &             SN,SNLX,SNLY,SNLX)
    
     CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
     &             SNCLOC,SNGI, &
     &             SNC,SNCLX,SNCLY,SNCLX)
     
     ALLOCATE(FINELE(TOTELE+1))
     ALLOCATE(COLELE4(TOTELE*5))

     ewrite(1,*) 'in subroutine DIFF3DSGS_DG going into GETFINELE4'
     CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE4,NONODS)
     ewrite(1,*) 'finished GETFINELE4'
      
! The local cords are in anti-clockwise order
     IFACE=1
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=3
     IFACE=2
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=4
     IFACE=3
     LOCLIST(IFACE,1)=3
     LOCLIST(IFACE,2)=2
     LOCLIST(IFACE,3)=4
     IFACE=4
     LOCLIST(IFACE,1)=1
     LOCLIST(IFACE,2)=3
     LOCLIST(IFACE,3)=4
! calculate approximate volume of eache element...

! Form critical matrices for DG diffusion:
     ALLOCATE(MATDIFX(TOTELE,NLOC,NLOC))
     ALLOCATE(MATDIFY(TOTELE,NLOC,NLOC))
     ALLOCATE(MATDIFZ(TOTELE,NLOC,NLOC))
     ALLOCATE(MATVLK(TOTELE,NLOC,NLOC))
     ALLOCATE(VOLELE(TOTELE))
     ALLOCATE(MEANXX(TOTELE))
     ALLOCATE(MEANXY(TOTELE))
     ALLOCATE(MEANXZ(TOTELE))
     ALLOCATE(MEANYY(TOTELE))
     ALLOCATE(MEANYZ(TOTELE))
     ALLOCATE(MEANZZ(TOTELE))
       
     VOLELE=0.0
     MEANXX=0.0
     MEANXY=0.0
     MEANXZ=0.0
     MEANYY=0.0
     MEANYZ=0.0
     MEANZZ=0.0
     

     ewrite(1,*) 'in subroutine DIFF3DSGS_DG going into 1st element loop'
      IF(ISTAB.EQ.1) THEN
        STAB_DIF_MAT=0.0
      ELSE
      BIGMELE(1:NBIGMELE) = 0.0 
     first_element_loop: do  ELE=1,TOTELE
    
! Get determinant and derivatives...
        CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
            & N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
            & NX,NY,NZ, NCX,NCY,NCZ,&
            & A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            & XD,YD,ZD,&
            & ISPHERE) 
     
        gauss_point_loop: do  GI=1,NGI
        
             IF(COGRAX) THEN
               DIRX(GI)=0.0
               DIRY(GI)=0.0
               DIRZ(GI)=1.0
             ELSE
               DIRX(GI)=XD(GI)
               DIRY(GI)=YD(GI)
               DIRZ(GI)=ZD(GI)
               RN=SQRT(DIRX(GI)**2+DIRY(GI)**2+DIRZ(GI)**2)
               DIRX(GI)=DIRX(GI)/RN
               DIRY(GI)=DIRY(GI)/RN
               DIRZ(GI)=DIRZ(GI)/RN
             ENDIF
             
          VOLELE(ELE)=VOLELE(ELE)+DETWEI(GI)
             
          MUPTXXGI=0.0
          MUPTXYGI=0.0
          MUPTXZGI=0.0
          MUPTYYGI=0.0
          MUPTYZGI=0.0
          MUPTZZGI=0.0

          do  L=1,NLOC
             IGLV=VONDGL((ELE-1)*NLOC+L)
      
             MUPTXXGI=MUPTXXGI+N(L,GI)*MUPTXX(IGLV)
             MUPTXYGI=MUPTXYGI+N(L,GI)*MUPTXY(IGLV)
             MUPTXZGI=MUPTXZGI+N(L,GI)*MUPTXZ(IGLV)
             MUPTYYGI=MUPTYYGI+N(L,GI)*MUPTYY(IGLV)
             MUPTYZGI=MUPTYZGI+N(L,GI)*MUPTYZ(IGLV)
             MUPTZZGI=MUPTZZGI+N(L,GI)*MUPTZZ(IGLV)
          end do

! calculate tensor for oceans:
! this has a QUADRATURE POINT varying element size and shape
          CALL SIZEGIELETENS(NX,NY,NZ,NLOC,NGI,GI,&
               & L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ )

!          IF((DISOPT.GE.147).AND.(DISOPT.LE.154)) THEN
          IF((DISOPT.GE.147).AND.(DISOPT.LE.156)) THEN

             MXXDTW(GI)=L2GIXX
             MXYDTW(GI)=L2GIXY
             MXZDTW(GI)=L2GIXZ
             MYYDTW(GI)=L2GIYY
             MYZDTW(GI)=L2GIYZ
             MZZDTW(GI)=L2GIZZ
             CALL DG_LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
                  &             NX,NY,NZ, &
                  &             DGNU,DGNV,DGNW, UG,VG,WG, &
                  &             MXXDTW,MXYDTW,MXZDTW,&
                  &             MYYDTW,MYZDTW,MZZDTW, &
     &             (DISOPT.EQ.155).OR.(DISOPT.EQ.156),DIRX,DIRY,DIRZ)

             MUGIXX=MXXDTW(GI)
             MUGIXY=MXYDTW(GI)
             MUGIXZ=MXZDTW(GI)
             MUGIYY=MYYDTW(GI)
             MUGIYZ=MYZDTW(GI)
             MUGIZZ=MZZDTW(GI)
          ELSE
             MUGIXX=0.0
             MUGIXY=0.0
             MUGIXZ=0.0
             MUGIYY=0.0
             MUGIYZ=0.0
             MUGIZZ=0.0
          ENDIF
       
! 0.25 is used to scale down LES dissipation from ocean dissipation
! factor of 4 is used to reduce the length scale.  
          IF(DISOPT.EQ.148) THEN
! Put LES only in the matrix AMAT.
             MXXDTW(GI)=(MUPTXXGI)*DETWEI(GI)
             MXYDTW(GI)=(MUPTXYGI)*DETWEI(GI)
             MXZDTW(GI)=(MUPTXZGI)*DETWEI(GI)

             MYYDTW(GI)=(MUPTYYGI)*DETWEI(GI)
             MYZDTW(GI)=(MUPTYZGI)*DETWEI(GI)

             MZZDTW(GI)=(MUPTZZGI)*DETWEI(GI)
          ELSE
             MXXDTW(GI)=(RLES*MUGIXX+RHAVEVIS*MUPTXXGI)*DETWEI(GI)
             MXYDTW(GI)=(RLES*MUGIXY+RHAVEVIS*MUPTXYGI)*DETWEI(GI)
             MXZDTW(GI)=(RLES*MUGIXZ+RHAVEVIS*MUPTXZGI)*DETWEI(GI)

             MYYDTW(GI)=(RLES*MUGIYY+RHAVEVIS*MUPTYYGI)*DETWEI(GI)
             MYZDTW(GI)=(RLES*MUGIYZ+RHAVEVIS*MUPTYZGI)*DETWEI(GI)

             MZZDTW(GI)=(RLES*MUGIZZ+RHAVEVIS*MUPTZZGI)*DETWEI(GI)
          ENDIF

          MEANXX(ELE)=MEANXX(ELE)+MXXDTW(GI)
          MEANXY(ELE)=MEANXY(ELE)+MXYDTW(GI)
          MEANXZ(ELE)=MEANXZ(ELE)+MXZDTW(GI)
          MEANYY(ELE)=MEANYY(ELE)+MYYDTW(GI)
          MEANYZ(ELE)=MEANYZ(ELE)+MYZDTW(GI)
          MEANZZ(ELE)=MEANZZ(ELE)+MZZDTW(GI)
       end do gauss_point_loop

       MEANXX(ELE)=MEANXX(ELE)/VOLELE(ELE)
       MEANXY(ELE)=MEANXY(ELE)/VOLELE(ELE)
       MEANXZ(ELE)=MEANXZ(ELE)/VOLELE(ELE)
       MEANYY(ELE)=MEANYY(ELE)/VOLELE(ELE)
       MEANYZ(ELE)=MEANYZ(ELE)/VOLELE(ELE)
       MEANZZ(ELE)=MEANZZ(ELE)/VOLELE(ELE)

! now form matrices: MATDIFX,MATDIFY,MATDIFZ,MASSDGELE

       do ILOC=1,NLOC
          do JLOC=1,NLOC
             NN=0.0
             MNNX=0.0
             MNNY=0.0
             MNNZ=0.0
             VLKTEN=0.0
             do GI=1,NGI
                NN=NN+N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                MNNX=MNNX+N(ILOC,GI)*(MXXDTW(GI)*NX(JLOC,GI) &
                     & +MXYDTW(GI)*NY(JLOC,GI)&
                     & +MXZDTW(GI)*NZ(JLOC,GI))
                MNNY=MNNY+N(ILOC,GI)*(MXYDTW(GI)*NX(JLOC,GI) &
                     & +MYYDTW(GI)*NY(JLOC,GI)&
                     & +MYZDTW(GI)*NZ(JLOC,GI))
                MNNZ=MNNZ+N(ILOC,GI)*(MXZDTW(GI)*NX(JLOC,GI) &
                     & +MYZDTW(GI)*NY(JLOC,GI)&
                     & +MZZDTW(GI)*NZ(JLOC,GI))
                XNDMUJLOC=(NX(JLOC,GI)*MXXDTW(GI)&
                     & +NY(JLOC,GI)*MXYDTW(GI)+NZ(JLOC,GI)*MXZDTW(GI))
                YNDMUJLOC=(NX(JLOC,GI)*MXYDTW(GI)&
                     & +NY(JLOC,GI)*MYYDTW(GI)+NZ(JLOC,GI)*MYZDTW(GI))
                ZNDMUJLOC=(NX(JLOC,GI)*MXZDTW(GI)&
                     & +NY(JLOC,GI)*MYZDTW(GI)+NZ(JLOC,GI)*MZZDTW(GI))
     
                TEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC &
                     & +NZ(ILOC,GI)*ZNDMUJLOC
                VLKTEN=VLKTEN+TEN
             END DO
             MASSDGELE(ELE,ILOC,JLOC)=NN
             MATDIFX(ELE,ILOC,JLOC)=MNNX
             MATDIFY(ELE,ILOC,JLOC)=MNNY
             MATDIFZ(ELE,ILOC,JLOC)=MNNZ
             MATVLK(ELE,ILOC,JLOC)=VLKTEN
          END DO
       END DO
             
    END DO first_element_loop
    ENDIF
      
    ewrite(1,*) 'in subroutine DIFF3DSGS_DG going into 2nd element loop'

    second_element_loop: do  ELE=1,TOTELE

! Get determinant and derivatives...
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
            & N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
            & NX,NY,NZ, NCX,NCY,NCZ,&
            & A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            & XD,YD,ZD,&
            & ISPHERE) 
     
! ***PUT ALL SOURCES IN TOTSOU*****
       do ILOC=1,NLOC
          IGL=SONDGL((ELE-1)*NLOC+ILOC)
          IF(ISUBSOU.EQ.0) THEN
             TOTSOUT(ILOC)=SOURCT(IGL)
          ELSE
             NODDG=(ELE-1)*NLOC+ILOC
             TOTSOUT(ILOC)=SOURCT(IGL)+SUBSOUT(NODDG)
          ENDIF
       END DO
     
         IF(SGSUPWIND.AND.(.NOT.ISOTROPIC)) THEN
           HHORIZ=0.0
           DO ILOC=1,NLOC
             IGLX=XONDGL((ELE-1)*NLOC+ILOC)
             DIST2=0.0
             DO JLOC=1+ILOC,NLOC
               JGLX=XONDGL((ELE-1)*NLOC+JLOC)
               DIST2=(X(IGLX)-X(JGLX))**2+(Y(IGLX)-Y(JGLX))**2 &
     &              +(Z(IGLX)-Z(JGLX))**2
             END DO
             HHORIZ=MAX(HHORIZ,DIST2)
           END DO
           HHORIZ=SQRT(HHORIZ)
         ENDIF
         
       do  GI=1,NGI! Was loop 331
       
             IF(COGRAX) THEN
               DIRX(GI)=0.0
               DIRY(GI)=0.0
               DIRZ(GI)=1.0
             ELSE
               DIRX(GI)=XD(GI)
               DIRY(GI)=YD(GI)
               DIRZ(GI)=ZD(GI)
               RN=SQRT(DIRX(GI)**2+DIRY(GI)**2+DIRZ(GI)**2)
               DIRX(GI)=DIRX(GI)/RN
               DIRY(GI)=DIRY(GI)/RN
               DIRZ(GI)=DIRZ(GI)/RN
             ENDIF

          UD(GI)=0.    
          VD(GI)=0.
          WD(GI)=0.
          UDDX(GI)=0.    
          VDDY(GI)=0.
          WDDZ(GI)=0.
          UDDX2(GI)=0.    
          VDDY2(GI)=0.
          WDDZ2(GI)=0.
             
          UGDX(GI)=0.0
          VGDY(GI)=0.0
          WGDZ(GI)=0.0
             
          TOLDGI(GI)=0.0
          TOLDGIX(GI)=0.0
          TOLDGIY(GI)=0.0
          TOLDGIZ(GI)=0.0
          TNEWGI(GI)=0.0
          TNEWGIX(GI)=0.0
          TNEWGIY(GI)=0.0
          TNEWGIZ(GI)=0.0
             
          ABGI(GI)=0.0
          SOUGI(GI)=0.0
             
             DO ILOC=1,NLOC
                IDG=(ELE-1)*NLOC+ILOC
                UD(GI)=UD(GI)+N(ILOC,GI)*DGNU(IDG)
                VD(GI)=VD(GI)+N(ILOC,GI)*DGNV(IDG)
                WD(GI)=WD(GI)+N(ILOC,GI)*DGNW(IDG)
                UDDX(GI)=UDDX(GI) + NX(ILOC,GI)*DGNU(IDG)
                VDDY(GI)=VDDY(GI) + NY(ILOC,GI)*DGNV(IDG)
                WDDZ(GI)=WDDZ(GI) + NZ(ILOC,GI)*DGNW(IDG)
                TOLDGI(GI)=TOLDGI(GI) + N(ILOC,GI)*DGT(IDG)
                TOLDGIX(GI)=TOLDGIX(GI) + NX(ILOC,GI)*DGT(IDG)
                TOLDGIY(GI)=TOLDGIY(GI) + NY(ILOC,GI)*DGT(IDG)
                TOLDGIZ(GI)=TOLDGIZ(GI) + NZ(ILOC,GI)*DGT(IDG)
                TNEWGI(GI)=TNEWGI(GI) + N(ILOC,GI)*DGTNEW(IDG)
                TNEWGIX(GI)=TNEWGIX(GI) + NX(ILOC,GI)*DGTNEW(IDG)
                TNEWGIY(GI)=TNEWGIY(GI) + NY(ILOC,GI)*DGTNEW(IDG)
                TNEWGIZ(GI)=TNEWGIZ(GI) + NZ(ILOC,GI)*DGTNEW(IDG)
             END DO
             TTHETAGI(GI)=THETA*TNEWGI(GI)+(1.-THETA)*TOLDGI(GI)
             TTHETAGIX(GI)=THETA*TNEWGIX(GI)+(1.-THETA)*TOLDGIX(GI)
             TTHETAGIY(GI)=THETA*TNEWGIY(GI)+(1.-THETA)*TOLDGIY(GI)
             TTHETAGIZ(GI)=THETA*TNEWGIZ(GI)+(1.-THETA)*TOLDGIZ(GI)
             
          IF(COGRAX) THEN
             DIRX(GI)=0.0
             DIRY(GI)=0.0
             DIRZ(GI)=1.0
          ELSE
             DIRX(GI)=XD(GI)
             DIRY(GI)=YD(GI)
             DIRZ(GI)=ZD(GI)
             RN=SQRT(DIRX(GI)**2+DIRY(GI)**2+DIRZ(GI)**2)
             DIRX(GI)=DIRX(GI)/RN
             DIRY(GI)=DIRY(GI)/RN
             DIRZ(GI)=DIRZ(GI)/RN
          ENDIF

          do  L=1,NLOC
             IGLV=VONDGL((ELE-1)*NLOC+L)
        
             IF(DISOPN.EQ.0) THEN
                UD(GI)=UD(GI) + N(L,GI)*(-UG(IGLV))
                VD(GI)=VD(GI) + N(L,GI)*(-VG(IGLV))
                WD(GI)=WD(GI) + N(L,GI)*(-WG(IGLV))
     
                UGDX(GI)=UGDX(GI) + NX(L,GI)*UG(IGLV)
                VGDY(GI)=VGDY(GI) + NY(L,GI)*VG(IGLV)
                WGDZ(GI)=WGDZ(GI) + NZ(L,GI)*WG(IGLV)
             ENDIF
           
             ABGI(GI)=ABGI(GI)+N(L,GI)*ABSORB(IGLV)
             SOUGI(GI)=SOUGI(GI)+N(L,GI)*TOTSOUT(L)
          end do

          ALPHA1=1.0
          USTAR(GI)=0.0
          VSTAR(GI)=0.0
          WSTAR(GI)=0.0
          DO L=1,NLOC
            SGS_UPWIND_FAC(L,GI)=0.0
            SGS_UPWIND_FAC_DIF(L,GI)=0.0
          END DO
          IF(SGSUPWIND) THEN
            A_DOT_GRADT=UD(GI)*TTHETAGIX(GI)+VD(GI)*TTHETAGIY(GI)+WD(GI)*TTHETAGIZ(GI)
            TTHETADIVU=TTHETAGI(GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI))
            RESID=(TNEWGI(GI)-TOLDGI(GI))/DT + A_DOT_GRADT + BETA*TTHETADIVU &
     &           +ABGI(GI)*TTHETAGI(GI) - SOUGI(GI)
! Johnson's method...
            IF(ISOTROPIC) THEN
! Use distance in the direction of
                 HXGI=ABS(A11(GI)*TTHETAGIX(GI)+A12(GI)*TTHETAGIY(GI)+A13(GI)*TTHETAGIZ(GI))
                 HYGI=ABS(A21(GI)*TTHETAGIX(GI)+A22(GI)*TTHETAGIY(GI)+A23(GI)*TTHETAGIZ(GI))
                 HZGI=ABS(A31(GI)*TTHETAGIX(GI)+A32(GI)*TTHETAGIY(GI)+A33(GI)*TTHETAGIZ(GI))
                 HOVERQ=1.0/MAX(HXGI,HYGI,HZGI,1.0E-9)
                 HDIST=HOVERQ
                 COEF=1.0
             ELSE
                 HDIST=HHORIZ
                 COEF=1.0/MAX(1.E-9, &
     &             abs(TTHETAGIX(GI)*(1.-ABS(DIRX(GI))))+abs(TTHETAGIY(GI)*(1.-ABS(DIRY(GI)))) &
     &            +abs(TTHETAGIZ(GI)*(1.-ABS(DIRZ(GI)))) )
             ENDIF
             DO L=1,NLOC
               SGS_UPWIND_FAC_DIF(L,GI)=0.2*HDIST*ABS(RESID)*COEF
!               SGS_UPWIND_FAC_DIF(L,GI)=2.0*HDIST*ABS(RESID)*COEF
             END DO
          ENDIF

       end do ! Was loop 331
       
        IF(ISTAB.EQ.1) THEN
           DO ILOC=1,NLOC
            DO JLOC=1,NLOC
               DIFSTAB =0.0
!     
               DO GI=1,NGI
                 DIFSTAB =DIFSTAB &
     &  +(NX(iLOC,GI)*(1.-RANISO*ABS(DIRX(GI)))+NY(iLOC,GI)*(1.-RANISO*ABS(DIRY(GI))) &
     &   +NZ(iLOC,GI)*(1.-RANISO*ABS(DIRZ(GI)))) &
     &  *(NX(jLOC,GI)*(1.-RANISO*ABS(DIRX(GI)))+NY(jLOC,GI)*(1.-RANISO*ABS(DIRY(GI))) &
     &   +NZ(jLOC,GI)*(1.-RANISO*ABS(DIRZ(GI)))) &
     &   *SGS_UPWIND_FAC_DIF(ILOC,GI)*DETWEI(GI)
               END DO
               STAB_DIF_MAT(ELE,ILOC,JLOC)=DIFSTAB
!               STAB_DIF_MAT(ELE,ILOC,JLOC)=0.0
            END DO
          END DO
        
! ENDOF IF(ISTAB.EQ.1) THEN...
        ELSE
         
         
       MASSMM=0.0
       MASSLUMMM=0.0
       FAMAT    =0.0
       ILOCELE=CENTRMELE(ELE)-FINDRMELE(ELE)+1

! AMAT**********************
       do  ILOC=1,NLOC! Was loop 3501
          do  JLOC=1,NLOC! Was loop 3601

             VLM=0. 
             VLND  =0.0
             VLND2 =0.0
             VLMAB=0.0
             DIFSTAB =0.0

             do  GI=1,NGI
     
                 VLN =-(NX(ILOC,GI)*UD(GI) &
     &                 +NY(ILOC,GI)*VD(GI) &
     &                 +NZ(ILOC,GI)*WD(GI) &
     &           -N(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*N(JLOC,GI)*DETWEI(GI) 
     
!                 VLN2 =-N(ILOC,GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI))*N(JLOC,GI)*DETWEI(GI) 
                 
                 VLN2 =N(ILOC,GI)*(NX(jLOC,GI)*UD(GI) &
     &                 +NY(jLOC,GI)*VD(GI) &
     &                 +NZ(jLOC,GI)*WD(GI))*DETWEI(GI) &
     &  +(NX(iLOC,GI)*USTAR(GI)+NY(iLOC,GI)*VSTAR(GI)+NZ(iLOC,GI)*WSTAR(GI)) &
     &  *(NX(jLOC,GI)*UD(GI)+NY(jLOC,GI)*VD(GI)+NZ(jLOC,GI)*WD(GI)) &
     &   *SGS_UPWIND_FAC(ILOC,GI)*DETWEI(GI) 
             
                 VLND=VLND+VLN
                 VLND2=VLND2+VLN2
             
                 RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI) &
     &  +(NX(iLOC,GI)*USTAR(GI)+NY(iLOC,GI)*VSTAR(GI)+NZ(iLOC,GI)*WSTAR(GI)) &
     &  *N(jLOC,GI)*SGS_UPWIND_FAC(ILOC,GI)*DETWEI(GI) 
          
                 VLM=VLM+RNN 
             
                 VLMAB=VLMAB+ABGI(GI)*RNN 
                 
                 DIFSTAB =DIFSTAB &
     &  +(NX(iLOC,GI)*(1.-RANISO*ABS(DIRX(GI)))+NY(iLOC,GI)*(1.-RANISO*ABS(DIRY(GI))) &
     &   +NZ(iLOC,GI)*(1.-RANISO*ABS(DIRZ(GI)))) &
     &  *(NX(jLOC,GI)*(1.-RANISO*ABS(DIRX(GI)))+NY(jLOC,GI)*(1.-RANISO*ABS(DIRY(GI))) &
     &   +NZ(jLOC,GI)*(1.-RANISO*ABS(DIRZ(GI)))) &
     &   *SGS_UPWIND_FAC_DIF(ILOC,GI)*DETWEI(GI)

             end do
         
! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
! conservative form BETA=1...
! non-conservative form BETA=0...
             VLKTEN=MATVLK(ELE,ILOC,JLOC)
!             AMAT(ILOC,JLOC)=BETA*VLND+(1.-BETA)*(VLND2) +VLMAB + VLKTEN
             AMAT(ILOC,JLOC)=BETA*VLND      +(1.-BETA)*(VLND2) +VLMAB + VLKTEN +DIFSTAB
             
             IF(LUMP) THEN
               MASSMM(ILOC,ILOC)=MASSMM(ILOC,ILOC)+VLM 
             ELSE
               MASSMM(ILOC,JLOC)=VLM
             ENDIF

             MASSLUMMM(ILOC,ILOC)=MASSLUMMM(ILOC,ILOC)+VLM
             
          end do ! Was loop 3601
       end do ! Was loop 3501
 
! calculate SNX,SNY,SNZ***
! Surface element of domain...
       do  IFACE=1,NFACE! Was loop 3344
         
! NB colele4 is ordered in terms of faces.
          ELE2=COLELE4((ELE-1)*5+IFACE)
          do COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
             JELE=COLMELE(COUNT)
             IF(JELE.EQ.ELE2) THEN
                COUNTELE=COUNT
                JLOCELE=COUNT-FINDRMELE(ELE)+1
             END IF
          END DO
         
          do SILOC=1,SNLOC
             SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
          END DO

! Put the volume elements togoether for element pairing...
          IF(ELE2.EQ.0) THEN
             SMEANXX=MEANXX(ELE)
             SMEANXY=MEANXY(ELE)
             SMEANXZ=MEANXZ(ELE)
             SMEANYY=MEANYY(ELE)
             SMEANYZ=MEANYZ(ELE)
             SMEANZZ=MEANZZ(ELE)
             ELE2_LOC_NODS=0
             ILOC_OTHER_SIDE=0
          ELSE
             SMEANXX=(MEANXX(ELE)*VOLELE(ELE)+MEANXX(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
             SMEANXY=(MEANXY(ELE)*VOLELE(ELE)+MEANXY(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
             SMEANXZ=(MEANXZ(ELE)*VOLELE(ELE)+MEANXZ(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
             SMEANYY=(MEANYY(ELE)*VOLELE(ELE)+MEANYY(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
             SMEANYZ=(MEANYZ(ELE)*VOLELE(ELE)+MEANYZ(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))
             SMEANZZ=(MEANZZ(ELE)*VOLELE(ELE)+MEANZZ(ELE2)*VOLELE(ELE2))/(VOLELE(ELE)+VOLELE(ELE2))

! Calculate the nodes on the other side of the face:
             do SILOC=1,SNLOC
                ILOC=SILOC2ILOC(SILOC)
                INOD=XONDGL((ELE-1)*NLOC+ILOC)
                do ILOC2=1,NLOC
                   INOD2=XONDGL((ELE2-1)*NLOC+ILOC2)
                   IF(INOD2.EQ.INOD) ILOC_OTHER_SIDE(SILOC)=ILOC2
                END DO
             END DO

! Calculate ELE2_LOC_NODS:
             ELE2_LOC_NODS=0
             do SILOC=1,SNLOC
                ILOC=SILOC2ILOC(SILOC)
                ILOC2=ILOC_OTHER_SIDE(SILOC)
                ELE2_LOC_NODS(ILOC2)=ILOC
             END DO
             do ILOC2=1,NLOC
                IF(ELE2_LOC_NODS(ILOC2).EQ.0) ELE2_LOC_NODS(ILOC2)=NLOC+1
             END DO
          ENDIF

          MATDIFX5=0.0
          MATDIFY5=0.0
          MATDIFZ5=0.0
          MASSMM5=0.0
          SUFMATX=0.0
          SUFMATY=0.0
          SUFMATZ=0.0

          do SILOC=1,SNLOC
             ILOC=SILOC2ILOC(SILOC)
             INOD=XONDGL((ELE-1)*NLOC+ILOC)
             do SJLOC=SILOC,SNLOC,1
                JLOC=SILOC2ILOC(SJLOC)
                JNOD=XONDGL((ELE-1)*NLOC+JLOC)
                SKLOC=SILINK2(SILOC,SJLOC)
                XSL(SKLOC)=0.5*(X(INOD)+X(JNOD))
                YSL(SKLOC)=0.5*(Y(INOD)+Y(JNOD))
                ZSL(SKLOC)=0.5*(Z(INOD)+Z(JNOD))
! Assume the mide side nodes are on the sphere
                IF(ISPHERE.GE.2) THEN
                   IF(ILOC.NE.JLOC) THEN
                      RADI=SQRT(X(INOD)**2+Y(INOD)**2+Z(INOD)**2)
                      RADJ=SQRT(X(JNOD)**2+Y(JNOD)**2+Z(JNOD)**2)
                      RADP=0.5*(RADI+RADJ)
                      RAD=SQRT(XSL(SKLOC)**2+YSL(SKLOC)**2+ZSL(SKLOC)**2)
                      XSL(SKLOC)=XSL(SKLOC)*RADP/RAD
                      YSL(SKLOC)=YSL(SKLOC)*RADP/RAD
                      ZSL(SKLOC)=ZSL(SKLOC)*RADP/RAD
                   ENDIF
                ENDIF
             END DO
          END DO
          
! Form approximate surface normal (NORMX,NORMY,NORMZ)
          CALL DGSIMPLNORM(ELE,SILOC2ILOC,TOTELE,NLOC,SNLOC,XONDGL, &
     &                     X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)
! Recalculate the normal...
          CALL DGSDETNXLOC2(SNCLOC,SNGI, &
     &        XSL,YSL,ZSL, &
     &        SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
     &        SNORMXN,SNORMYN,SNORMZN, &
     &        NORMX,NORMY,NORMZ)

          SUD=0.0
          SVD=0.0
          SWD=0.0
          IN_SUD=0.0
          IN_SVD=0.0
          IN_SWD=0.0
          OUT_SUD=0.0
          OUT_SVD=0.0
          OUT_SWD=0.0
          do SILOC=1,SNLOC
             ILOC=SILOC2ILOC(SILOC)
             IGL =NDGLNO((ELE-1)*NLOC+ILOC)
             IDG=(ELE-1)*NLOC+ILOC
             IF(.NOT.SUF_CONSERV) THEN
                IDG2=IDG
             ELSE 
                IF(ELE2.NE.0) THEN
                   ILOC2=ILOC_OTHER_SIDE(SILOC)
                   IDG2=(ELE2-1)*NLOC+ILOC2
                ELSE
                   IDG2=IDG
                ENDIF
             ENDIF
             do  SGI=1,SNGI! Was loop 3313
                SUD(SGI)=SUD(SGI) + SN(SILOC,SGI)*(0.5*(DGNU(IDG)+DGNU(IDG2))-UG(IGL))
                SVD(SGI)=SVD(SGI) + SN(SILOC,SGI)*(0.5*(DGNV(IDG)+DGNV(IDG2))-VG(IGL))
                SWD(SGI)=SWD(SGI) + SN(SILOC,SGI)*(0.5*(DGNW(IDG)+DGNW(IDG2))-WG(IGL))
                 
                IN_SUD(SGI)=IN_SUD(SGI) + SN(SILOC,SGI)*(DGNU(IDG2)-UG(IGL))
                IN_SVD(SGI)=IN_SVD(SGI) + SN(SILOC,SGI)*(DGNV(IDG2)-VG(IGL))
                IN_SWD(SGI)=IN_SWD(SGI) + SN(SILOC,SGI)*(DGNW(IDG2)-WG(IGL))
                 
                OUT_SUD(SGI)=OUT_SUD(SGI) + SN(SILOC,SGI)*(DGNU(IDG)-UG(IGL))
                OUT_SVD(SGI)=OUT_SVD(SGI) + SN(SILOC,SGI)*(DGNV(IDG)-VG(IGL))
                OUT_SWD(SGI)=OUT_SWD(SGI) + SN(SILOC,SGI)*(DGNW(IDG)-WG(IGL))
             END DO
          END DO
           
          IF(ADVECT_MEAN) THEN
! Original method...
             IN_SUD=SUD
             IN_SVD=SVD
             IN_SWD=SWD

             OUT_SUD=SUD
             OUT_SVD=SVD 
             OUT_SWD=SWD 
          ENDIF

          do SGI=1,SNGI
             NDOTQ(SGI)=SUD(SGI)*SNORMXN(SGI)+SVD(SGI)*SNORMYN(SGI) &
     &                 +SWD(SGI)*SNORMZN(SGI)
             IF(NDOTQ(SGI).LT.0.0) THEN
                IN_NDOTQ(SGI)=IN_SUD(SGI)*SNORMXN(SGI)+IN_SVD(SGI)*SNORMYN(SGI) &
     &                      +IN_SWD(SGI)*SNORMZN(SGI)
                OUT_NDOTQ(SGI)=0.0
             ELSE
                IN_NDOTQ(SGI)=0.0
                OUT_NDOTQ(SGI)=OUT_SUD(SGI)*SNORMXN(SGI) &
                    & +OUT_SVD(SGI)*SNORMYN(SGI) &
                    & +OUT_SWD(SGI)*SNORMZN(SGI)
             ENDIF
          END DO
           
          AVE_SNORMXN=0.0
          AVE_SNORMYN=0.0
          AVE_SNORMZN=0.0
          do SGI=1,SNGI
             AVE_SNORMXN=AVE_SNORMXN+SNORMXN(SGI)/SNGI
             AVE_SNORMYN=AVE_SNORMYN+SNORMYN(SGI)/SNGI
             AVE_SNORMZN=AVE_SNORMZN+SNORMZN(SGI)/SNGI
          END DO
! ************************
! Perform surface integration...
          ONE_OR_ZERO=0.0
          BNDSUF=.false.
          SUF_TOP=.FALSE.
          IF(ELE2.EQ.0) THEN
! BNDSUF=.true. if we are on a boundary condition surface i.e. inlet-outlet
! for the free surface.
             ICOUNT=0
             do SILOC=1,SNLOC
                ILOC=SILOC2ILOC(SILOC)
                GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
                ICOUNT=ICOUNT+BNDCON(GLOBI)
             ENDDO
             BNDSUF=(ICOUNT.EQ.SNLOC)
             IF(SUF_TOP_OCEAN.OR.(PREOPT.EQ.2).or.BNDSUF) THEN
! Only put surface term in if at top of ocean...
                IF(COGRAX) THEN
                   LOCPNORMX=0.0
                   LOCPNORMY=0.0
                   LOCPNORMZ=1.0
                ELSE
                   XC=0.0
                   YC=0.0
                   ZC=0.0
                   do ILOC=1,NLOC
                      IXNOD=XONDGL((ELE-1)*NLOC+ILOC)
                      XC=XC+X(IXNOD)/NLOC
                      YC=YC+Y(IXNOD)/NLOC
                      ZC=ZC+Z(IXNOD)/NLOC
                   END DO
                   RN=SQRT(XC**2+YC**2+ZC**2)
                   LOCPNORMX=XC/RN
                   LOCPNORMY=YC/RN
                   LOCPNORMZ=ZC/RN
                ENDIF
                IF(AVE_SNORMXN*LOCPNORMX&
                    &       + AVE_SNORMYN*LOCPNORMY&
                    &       + AVE_SNORMZN*LOCPNORMZ.GT.0.8) THEN
                   ONE_OR_ZERO=1.0
                   SUF_TOP=.TRUE.
                ENDIF
                IF((PREOPT.EQ.2).or.BNDSUF) ONE_OR_ZERO=1.0
             ENDIF
          ENDIF

          suf_top=suf_top.or.(PREOPT.EQ.2).or.BNDSUF   
          do  SILOC=1,SNLOC! Was loop 35012
             ILOC   =SILOC2ILOC(SILOC)
             do  SJLOC=1,SNLOC! Was loop 35072
                JLOC =SILOC2ILOC(SJLOC)
                JLOC2=ILOC_OTHER_SIDE(SJLOC)
                VLNDOTQ_IN =0.0
                VLNDOTQ_OUT=0.0
                VLNDOTQ_ALL=0.0
                VLM_NORX=0.0
                VLM_NORY=0.0
                VLM_NORZ=0.0
! Have a surface integral on outlet boundary...  
                do SGI=1,SNGI
                   RNN=SDETWE(SGI)*SN(SILOC,SGI)*SN(SJLOC,SGI)
                   VLNDOTQ_IN =VLNDOTQ_IN +IN_NDOTQ(SGI)*RNN
                   VLNDOTQ_OUT=VLNDOTQ_OUT+OUT_NDOTQ(SGI)*RNN
                   VLNDOTQ_ALL=VLNDOTQ_ALL+(IN_NDOTQ(SGI)+OUT_NDOTQ(SGI))*RNN
                   VLM_NORX=VLM_NORX+SNORMXN(SGI)*RNN
                   VLM_NORY=VLM_NORY+SNORMYN(SGI)*RNN
                   VLM_NORZ=VLM_NORZ+SNORMZN(SGI)*RNN
                END DO

                assert(ILOCELE/=-1)
                assert(JLOCELE/=-1)
                IF(SUF_TOP) THEN
                   ! avoid applying any sort of b.c at the free surface
                   FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_ALL*BETA
                ELSE IF(ELE2.EQ.0) THEN
                   FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_ALL*BETA
                ELSE
                   FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) +VLNDOTQ_OUT*BETA&
     &                                                                                   -VLNDOTQ_IN*(1.-BETA)
                   FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC2)=FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC2)+VLNDOTQ_IN*BETA&
     &                                                                                   +VLNDOTQ_IN*(1.-BETA)
! add diffusion term...

                   MATDIFX5(ILOC,JLOC)=MATDIFX5(ILOC,JLOC)-(SMEANXX*VLM_NORX+SMEANXY*VLM_NORY+SMEANXZ*VLM_NORZ)
                   MATDIFY5(ILOC,JLOC)=MATDIFY5(ILOC,JLOC)-(SMEANXY*VLM_NORX+SMEANYY*VLM_NORY+SMEANYZ*VLM_NORZ)
                   MATDIFZ5(ILOC,JLOC)=MATDIFZ5(ILOC,JLOC)-(SMEANXZ*VLM_NORX+SMEANYZ*VLM_NORY+SMEANZZ*VLM_NORZ)
                   MATDIFX5(ILOC,JLOC2+NLOC)=MATDIFX5(ILOC,JLOC2+NLOC)+(SMEANXX*VLM_NORX+SMEANXY*VLM_NORY+SMEANXZ*VLM_NORZ)
                   MATDIFY5(ILOC,JLOC2+NLOC)=MATDIFY5(ILOC,JLOC2+NLOC)+(SMEANXY*VLM_NORX+SMEANYY*VLM_NORY+SMEANYZ*VLM_NORZ)
                   MATDIFZ5(ILOC,JLOC2+NLOC)=MATDIFZ5(ILOC,JLOC2+NLOC)+(SMEANXZ*VLM_NORX+SMEANYZ*VLM_NORY+SMEANZZ*VLM_NORZ)
                   SUFMATX(ILOC,JLOC)=SUFMATX(ILOC,JLOC)+VLM_NORX
                   SUFMATY(ILOC,JLOC)=SUFMATY(ILOC,JLOC)+VLM_NORY
                   SUFMATZ(ILOC,JLOC)=SUFMATZ(ILOC,JLOC)+VLM_NORZ
                ENDIF
             END DO
          END DO

          IF(ELE2.NE.0) THEN

! Calculate the diffusion terms*****
             MATDIFX5(1:NLOC,1:NLOC)=MATDIFX5(1:NLOC,1:NLOC)+MATDIFX(ELE,1:NLOC,1:NLOC)
             MATDIFY5(1:NLOC,1:NLOC)=MATDIFY5(1:NLOC,1:NLOC)+MATDIFY(ELE,1:NLOC,1:NLOC)
             MATDIFZ5(1:NLOC,1:NLOC)=MATDIFZ5(1:NLOC,1:NLOC)+MATDIFZ(ELE,1:NLOC,1:NLOC)
             MASSMM5(1:NLOC,1:NLOC) =MASSMM5(1:NLOC,1:NLOC) +MASSDGELE(ELE,1:NLOC,1:NLOC)

! Put contributions from ELE2 into MATDIFX5
             IF(ELE2.NE.0) THEN
                do ILOC2=1,NLOC
                   IILOC5=ELE2_LOC_NODS(ILOC2)
                   do JLOC2=1,NLOC
                      JJLOC5=ELE2_LOC_NODS(JLOC2)
                      MATDIFX5(IILOC5,JLOC2+NLOC)=MATDIFX5(IILOC5,JLOC2+NLOC)+MATDIFX(ELE2,ILOC2,JLOC2)
                      MATDIFY5(IILOC5,JLOC2+NLOC)=MATDIFY5(IILOC5,JLOC2+NLOC)+MATDIFY(ELE2,ILOC2,JLOC2)
                      MATDIFZ5(IILOC5,JLOC2+NLOC)=MATDIFZ5(IILOC5,JLOC2+NLOC)+MATDIFZ(ELE2,ILOC2,JLOC2)
                      MASSMM5(IILOC5,JJLOC5)=MASSMM5(IILOC5,JJLOC5)+MASSDGELE(ELE2,ILOC2,JLOC2)
                   END DO
                END DO
             ELSE
                MASSMM5(1+NLOC,1+NLOC) =1.0
             ENDIF
! calculate MASS2INV...
             CALL MATDMATINV(MASSMM5,MASSINVMM5,NLOC+1)
! M^{-1}*MATDIFX5
             CALL ABMATRIXMUL(MMATDIFX5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              MATDIFX5,  NLOC+1,NLOC+NLOC)
! M^{-1}*MATDIFY5
             CALL ABMATRIXMUL(MMATDIFY5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              MATDIFY5,  NLOC+1,NLOC+NLOC)
! M^{-1}*MATDIFZ5
             CALL ABMATRIXMUL(MMATDIFZ5,MASSINVMM5,NLOC+1,NLOC+1,&
     &                              MATDIFZ5,  NLOC+1,NLOC+NLOC)    
! SUFMATX*MMATDIFX5
             CALL ABMATRIXMUL(SMMATDIFX5,SUFMATX,   NLOC,NLOC+1,&
     &                               MMATDIFX5, NLOC+1,NLOC+NLOC)
! SUFMATX*MMATDIFX5
             CALL ABMATRIXMUL(SMMATDIFY5,SUFMATY,   NLOC,NLOC+1,&
     &                               MMATDIFY5, NLOC+1,NLOC+NLOC)
! SUFMATX*MMATDIFX5
             CALL ABMATRIXMUL(SMMATDIFZ5,SUFMATZ,   NLOC,NLOC+1,&
     &                               MMATDIFZ5, NLOC+1,NLOC+NLOC)
             SMMATDIF5=SMMATDIFX5+SMMATDIFY5+SMMATDIFZ5

             do SILOC=1,SNLOC
                ILOC=SILOC2ILOC(SILOC)
                do JLOC=1,NLOC
! add diffusion term (-VE because diffusion has - sign on lhs of eqns)...
                   FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(ILOCELE-1)*NLOC+JLOC)   &
     &                                           -SMMATDIF5(ILOC,JLOC)
                   IF(ELE2.NE.0) THEN
                      assert(JLOCELE/=-1)
                      FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC) =FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC)&
     &                                             -SMMATDIF5(ILOC,NLOC+JLOC)
                   ENDIF
                END DO
             END DO
! Calculate the diffusion terms***** 
          ENDIF

       end do ! Was loop 3344

! END OF Put surface integrals *********************

! AMAT:
       do ILOC=1,NLOC   
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          IDG=(ELE-1)*NLOC+ILOC 
          do JLOC=1,NLOC
             JDG=(ELE-1)*NLOC+JLOC 
             A11MAT=MASSMM(ILOC,JLOC)/DT -(1.-THETA)*AMAT(ILOC,JLOC)
             VECT(IDG)=VECT(IDG) +A11MAT*DGT(JDG) +MASSMM(ILOC,JLOC)*TOTSOUT(JLOC)
! Lump...
             MLELE(IDG)=MLELE(IDG) + MASSLUMMM(ILOC,JLOC)
             ML(GLOBI) =ML(GLOBI)  + MASSLUMMM(ILOC,JLOC)
          END DO
       END DO
        
! now the surface integral for advection...       
       do COUNT=FINDRMELE(ELE),FINDRMELE(ELE+1)-1
          ELE2=COLMELE(COUNT)
          JLOCELE=COUNT-FINDRMELE(ELE)+1
          do ILOC=1,NLOC   
             IDG=(ELE-1)*NLOC+ILOC 
             do JLOC=1,NLOC
                JDG=(ELE2-1)*NLOC+JLOC 
                A11MAT= -(1.-THETA)*FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC)
                VECT(IDG)=VECT(IDG) +A11MAT*DGT(JDG) 
! ************PUT INTO MATRIX...
                THETA_FAMAT=THETA*FAMAT(ILOC,(JLOCELE-1)*NLOC+JLOC)
                IJDISP=(ILOC-1)*NLOC + JLOC
                BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)=BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)&
     &                                            +THETA_FAMAT
! ************PUT INTO MATRIX...
             END DO
          END DO
       END DO

       do ILOC=1,NLOC  
          IDG=(ELE-1)*NLOC+ILOC 
! ******add in the pre-discretised source*****
          VECT(IDG)=VECT(IDG) + VECT_SGSADD(IDG)
       END DO

! Amend AMAT,BMAT,CMAT,DMAT to take into account time dep.
       AMAT=(1./DT)*MASSMM+THETA*AMAT

! Now make bigger matrix for coupled vels
          
! Just this element...    
       COUNT=CENTRMELE(ELE)
       do ILOC=1,NLOC
          do JLOC=1,NLOC
! Put into matrix...
! Find count. ***************************************
             IJDISP=(ILOC-1)*NLOC+JLOC
             BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)=BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)&
     &             +AMAT(ILOC,JLOC)
          END DO
       END DO
         
! ENDOF IF(ISTAB.EQ.1) THEN ELSE...
         ENDIF
    end do second_element_loop

    CALL PMINMX(MUPTXX,NONODS,'MUPTXX*******  ')
    CALL PMINMX(MUPTXY,NONODS,'MUPTXY*******  ')
    CALL PMINMX(MUPTXZ,NONODS,'MUPTXZ*******  ')
    CALL PMINMX(MUPTYY,NONODS,'MUPTYY*******  ') 
    CALL PMINMX(MUPTYZ,NONODS,'MUPTYZ*******  ')
    CALL PMINMX(MUPTZZ,NONODS,'MUPTZZ*******  ')

  END SUBROUTINE HART3DSGS_DG
 
  SUBROUTINE APPLY_FIELD_DG_BC(&
       &     NLOC,TOTELE,NDGLNO,NONODS,&
       &     VECT,DGT,DT,NOBCT,BCT1,BCT2,&
       &     FINDRMELE,COLMELE,NCOLMELE,NBIGMELE,CENTRMELE,BIGMELE)
! Apply DG BC'S to field eqn...
    REAL INFINY
    PARAMETER(INFINY=1.0E+20)
    INTEGER NLOC,TOTELE,NDGLNO(TOTELE*NLOC),NONODS
    REAL VECT(TOTELE*NLOC),DGT(TOTELE*NLOC),DT
    INTEGER NOBCT
    INTEGER BCT2(NOBCT)
    REAL BCT1(NOBCT)
    INTEGER NCOLMELE,NBIGMELE
    INTEGER FINDRMELE(TOTELE+1),COLMELE(NCOLMELE),CENTRMELE(TOTELE)
    REAL BIGMELE(NBIGMELE)
! Local variables
    INTEGER II,ILOC,JLOC,NOD,GLOBI,IDG,ELE,COUNT,IJDISP
    LOGICAL, ALLOCATABLE, DIMENSION(:)::MARKT
    REAL, ALLOCATABLE, DIMENSION(:)::VALUE_T
    REAL, ALLOCATABLE, DIMENSION(:)::DGDMI
         
    ALLOCATE(MARKT(NONODS))
    ALLOCATE(VALUE_T(NONODS))
    ALLOCATE(DGDMI(TOTELE*NLOC))
    MARKT=.FALSE.
    VALUE_T=0.0
    DGDMI=0.0
! Amend DGDM1 to include bc's:
    do II=1,NOBCT
       NOD=BCT2(II)
       MARKT(NOD)=.TRUE.
       VALUE_T(NOD)=BCT1(II)
    END DO

    DGDMI=0.0
    do ELE=1,TOTELE
       do ILOC=1,NLOC
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          IDG=(ELE-1)*NLOC+ILOC
          IF(MARKT(GLOBI)) THEN
             DGDMI(IDG)=INFINY
             VECT(IDG)=INFINY*(DGT(IDG)+DT*VALUE_T(GLOBI))
          ENDIF
       END DO
    END DO
          
    do ELE=1,TOTELE
       COUNT=CENTRMELE(ELE)
       do ILOC=1,NLOC
          IDG=(ELE-1)*NLOC+ILOC
          JLOC=ILOC
          IJDISP=(ILOC-1)*NLOC + JLOC
          BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)= &
               & BIGMELE((COUNT-1)*NLOC*NLOC+IJDISP)+DGDMI(IDG)
       END DO
    END DO
        
  END SUBROUTINE APPLY_FIELD_DG_BC

end module diff3d_module
