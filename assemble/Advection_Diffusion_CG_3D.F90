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

module advection_diffusion_cg_3d

  use allsorts
  use Diffusivity
  use Element_Numbering
  use FETools
  use fields
  use FLDebug
  use hart3d_allsorts
  use mesh_connections
  use position_in_matrix
  use spud
  use Shape_Functions
  use shape_transformations
  use state_module
  use Transform_Elements
  use tr2d_module

#ifdef ADJOINT
  use diff3d_adjoint_module
#else
  use diff3d_module
#endif

  implicit none

  private

  public :: hart3d, hart3dsim, hart3dsgs

contains

  SUBROUTINE HART3D(U,NU,NV,NW,UG,VG,WG,&
       SOURCX,X,Y,Z,D0,&
       VECX, &
       C1T,C2T,C3T, BIGM, ML,&
       FINDCT,COLCT,NCT,FREDOP,TOTELE, &
       FINDRM,COLM,NCOLM,NONODS,CENTRM, &
       M,N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,MLOC,&
       DISOPT,NDISOT,DT,THETA,BETA,LUMP,MAKSYM,&
       ABSLUM,SLUMP,&
       GETC12,INMAT,MLCENT,&
       NDGLNO,PNDGLN,&
       SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
       DENPT,&
       MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
       ML2MXX,ML2MXY,ML2MXZ,ML2MYY,ML2MYZ,ML2MZZ, &
       ABSORB,&
       VERSIO, state, field)

    !     This subroutine forms the discretised momentum or a field equation. 
    !     The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
    !     DISOPT .LT. 0 then optimaly upwinded. 
    !     The following are absolute values of DISOPT...
    !     DISOPT=1 - Petrof Galerkin x-t.
    !     DISOPT=2 - Petrof Galerkin x-t(BUT weighted to get correct diffusion at steady state).
    !     DISOPT=3 - Petrof Galerkin (x,y).
    !     DISOPT=4 - Least squares (x,y,t). (symmetric)
    !     DISOPT=5 - Least squares (x,y).
    !     DISOPT=6 - Least squares THETA-method. (symmetric)
    !     DISOPT=7 - Galerkin THETA-method.(symmetric when MAKSYM=.true.)
    !     DISOPT=8 - Space weighting. 
    !     For DISOPT=7 ...
    !     If THETA=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    !     If THETA=1.0 then bakward Euler is used in NS equns. THETA=2/3 Galerkin. 
    !     Similarly for the temperature equation.
    !     DISOPT=(42-45)- LES option which uses no balancing diffusion and Galerkin THETA method.
    !     
    !     
    !     DISOPT=42- LES option using constant length scale.
    !     DISOPT=43- LES option using isotropic length scale.
    !     DISOPT=44- LES option which uses no balancing diffusion.
    !     DISOPT=45- LES option which uses no balancing diffusion.
    !     if NDISOT.GE.10 then set advection to zero and treat advection 
    !     using the high res method.
    !     IF MAKSYM then make BIGM symmetrix I.E take advection terms out of matrix.
    !     IF LUMP then lump the mass matrix else dont. 
    !     IF MLCENT then put ML to the diagonal of BIGM matrix. 
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
    !     -------------------------------------------------------------------------------------------
    !     If we are solving for another variable like temperature 
    !     or chemical species then NU,NV,NW will be the velocities 
    !     and U=TEMPERATURE OR CHEMICAL SPECIES. 
    !     -------------------------------------------------------------------------------------------
    !     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
    REAL RWIND
    PARAMETER(RWIND=0.5)
    INTEGER SNONOD,VNONOD,XNONOD,NCT,NCOLM,NLOC,NGI
    INTEGER MLOC,DISOPT,NDISOT,DISCRA
    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
    REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
    REAL ML2MXX(SNONOD),ML2MXY(SNONOD),ML2MXZ(SNONOD)
    REAL ML2MYY(SNONOD),ML2MYZ(SNONOD),ML2MZZ(SNONOD)
    REAL DENPT(SNONOD)
    REAL ABSORB(SNONOD)
    REAL  DT,THETA,BETA
    INTEGER TOTELE,NONODS,FREDOP
    !     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    !     If PHI=1.0 then bakward Euler is used in NS equns.
    !     Similarly for the temperature equation except width variable THETA.
    !     have gravity in -ve y direction.
    REAL U(NONODS)
    REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
    REAL MASS,MATRIX,RHSMAT,SOUMAT,VLMGI,VLKGI,VABS
    REAL::VLM=0.0
    REAL VF1M,VF2M,VF2,VF3M
    REAL THETA1,THETA2,THETA3,AHAT1,L1,ALPHA,THTWDT,INDT,TWOTHI
    REAL PE,HOVERQ,PWEIGI
    REAL GALERK,LEASQR
    REAL UD(NGI),VD(NGI),WD(NGI)
    REAL L2GIXX(NGI),L2GIXY(NGI),L2GIXZ(NGI)
    REAL L2GIYY(NGI),L2GIYZ(NGI),L2GIZZ(NGI)
    REAL DENGI(NGI),DENGI2(NGI),L1GI(NGI)
    REAL ABDGI(NGI)
    REAL AGI,BGI,CGI,DGI
    REAL EGI,FGI,GGI,HGI,KGI
    REAL A11,A12,A13, A21,A22,A23, A31,A32,A33
    REAL DETJ,DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    INTEGER VERSIO
    REAL GIVOLX,GIVOLY,GIVOLZ
    !     HX,HY-characteristic length scales in x,y directions.
    REAL HXGI,HYGI,HZGI
    REAL AHAT(NGI),GAMMA(NGI),GAMMAG(NGI),MASCON(NGI)
    REAL BECTWE(NGI)

    REAL VECX(NONODS)
    REAL BIGM(NCOLM)
    REAL ML(NONODS)
    REAL C1T(NCT),C2T(NCT),C3T(NCT)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL SOURCX(SNONOD)
    INTEGER NDGLNO(TOTELE*NLOC),PNDGLN(TOTELE*MLOC)
    INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER ELE,ILOC,JLOC,L
    INTEGER GI,COUNTER,I,GLOBI,GLOBJ,GLOBIS,GLOBJS
    INTEGER POS,LOWER,UPPER,PGLOBJ
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER CENTRM(NONODS)
    INTEGER FINDCT(FREDOP+1),COLCT(NCT)

    REAL M(MLOC,NGI)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL GETC12,MLCENT
    LOGICAL LUMP,MAKSYM,UPWIND,OPWIND,XTPETV
    LOGICAL ABSLUM,SLUMP
    LOGICAL INMAT
    REAL D0,DD0,INDD0
    LOGICAL NODIM
    LOGICAL GALWIN
    LOGICAL LES
    REAL GAMMC1(27),GAMMC2(27),GAMMC3(27)
    REAL ABLUMP,RSLUMP
    REAL R,RUD,RVD,RWD,R1,R2,RR
    REAL::RNN=0.0
    REAL RLUMP,RSYM,RCONSI,RGAMMA
    INTEGER IGLS,IGLV,IGLX,INUM

    !     Modern variables connected with diffusivity calculations.
    !     Diffusivity tensor (MU) at each Gauss point (GI). 
    real, dimension(3,3,ngi) :: MUGI
    !     density is the perturbation density field
    type(scalar_field), pointer :: density=>null()
    type(scalar_field), pointer :: topdis=>null()
    type(scalar_field), pointer :: botdis=>null()
    type(scalar_field), pointer :: d=>null()
    !     Shape function types.
    type(quadrature_type), save :: quad
    type(element_type), save :: tracer_shape, density_shape
    type(element_type), save :: topdis_shape, botdis_shape, d_shape
    real :: tracer_dshape(nloc, ngi, 3), density_dshape(nloc, ngi, 3) 
    !     Node positions in the current element
    real :: X_ele(3,nloc)
    !     Flag to ensure shape functions are only initialised once.
    logical, save :: initialised=.false.
    logical :: modern_diffusivity
    !     Degree of tracer elements.
    integer :: n_degree
    !     Encapsulated state information containing fields. -dham
    type(state_type), intent(inout) :: state
    !     The field being solved for.
    type(scalar_field), pointer :: field

    ewrite(1, *) "just in SUBROUTINE HART3D()"

    if (associated(field)) then
       if(have_option(trim(field%option_path)//&
            'prognostic/subgridscale_parameterisation::Gent_McWilliams')) then
          modern_diffusivity=.true.
       else
          modern_diffusivity=.false.
       end if
    else
       modern_diffusivity=.false.
    end if


    ewrite(2,*) 'cjc disopt = ',disopt

    if ((.not.initialised).and.modern_diffusivity) then
       initialised=.true.
       !     Form modern shape functions corresponding to the old-fashioned
       !         functions passed in.

       !     This is the inverse of the tetrahedral number algorithm n(n+1)(n+2)/6
       n_degree=inv_te(nloc)-1

       if (n_degree<0) then
          ewrite(-1,'(i0,a)') nloc, " is not a possible number of nodes in a tet"
          FLAbort("Dieing")
       end if

       !     Tets
       quad=make_quadrature(vertices=4, dim=3, ngi=ngi)

       tracer_shape=make_element_shape(vertices=4, dim=3, degree=n_degree, quad=quad)

       density_shape=tracer_shape
       topdis_shape=tracer_shape
       botdis_shape=tracer_shape
       d_shape=tracer_shape

    end if

    if (modern_diffusivity) then

       !     This loop really only applies to Gent-McWilliams - goodness knows what 
       !     other modern diffusivity schemes might need. But for now it'll have to
       !     do... (Jemma)

       density=>extract_scalar_field(state,"PerturbationDensity")

       if(.not. associated(density)) then
          FLAbort("Perturbation density is required for modern diffusivity")
       endif

       topdis=>extract_scalar_field(state,"DistanceToTop")

       if(.not. associated(topdis)) then
          FLAbort("DistanceToTop is required for modern diffusivity")
       endif

       botdis=>extract_scalar_field(state,"DistanceToBottom")

       if(.not. associated(botdis)) then
          FLAbort("DistanceToBottom is required for modern diffusivity")
       endif

       d=>extract_scalar_field(state,"DistanceToSideBoundaries")

       if(.not. associated(d)) then
          FLAbort("DistanceToSideBoundaries is required for modern diffusivity")
       endif

    endif

    do I=1,NONODS! Was loop 
       ML(I)=0.
       VECX(I)=0.
    END DO
    IF(INMAT) THEN
       do I=1,NCOLM! Was loop 
          BIGM(I)=0.
       END DO
    ENDIF
    IF(GETC12) THEN
       do I=1,NCT! Was loop 
          C1T(I)=0.
          C2T(I)=0.
          C3T(I)=0.
       END DO
    ENDIF

    RLUMP=0.
    RSYM=1.
    IF(LUMP) RLUMP=1.
    !     If MAKSYM then force the matrix to be symmetric 
    !     by treating advection explicitly. 
    IF(MAKSYM) RSYM=0.
    RCONSI=1.-RLUMP

    ABLUMP=0.
    RSLUMP=0.
    IF(ABSLUM) ABLUMP=1.
    IF(SLUMP)  RSLUMP=1.

    !     Ocean scaling.
    NODIM=.TRUE.
    IF(D0.EQ.0) NODIM=.FALSE.
    INDD0=1.
    DD0=1.
    IF(NODIM) THEN
       INDD0=1./D0
       DD0=D0
    ENDIF

    !     *********************************************
    INDT=1./DT
    TWOTHI=2./3.
    !     The defaults...
    GALERK=0.0
    LEASQR=0.0
    ALPHA =1.0
    L1    =1.0
    AHAT1 =1.0
    THETA1=1.0
    THETA2=1.0
    THETA3=1.0
    LES=.FALSE.
    IF((DISOPT.GE.42).AND.(DISOPT.LE.50)) LES=.TRUE.
    UPWIND=.FALSE.
    OPWIND=.FALSE.
    !     OPWIND=optimal upwinding method. 
    XTPETV=.FALSE.
    RGAMMA=1.0
    !     DISOPT=discretization option 
    DISCRA=ABS(DISOPT)
    IF(DISOPT.LT.0) OPWIND=.TRUE.
    IF(DISCRA.EQ.1) THEN
       !     Petrof Galerkin x-t
       XTPETV=.TRUE.
       UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.2) THEN
       !     Petrof Galerkin x-t(BUT weighted to get correct diffusion at steady state)
       UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.3) THEN
       !     Petrof Galerkin (x,y)
       AHAT1=0.
       UPWIND=.TRUE.
    ENDIF
    IF(DISCRA.EQ.4) THEN
       !     Least squares (x,y,t)
       THETA1=2./3.
       L1    =6./DT
       AHAT1 =-3.0
       LEASQR=1.0
       RGAMMA=1.0
    ENDIF
    IF(DISCRA.EQ.5) THEN
       !     Least squares (x,y)
       AHAT1=0.
       L1   =0.
    ENDIF
    IF(DISCRA.EQ.6) THEN
       !     Least squares THETA-method
       THETA3=THETA*3./2.
       AHAT1 =0.
       L1    =1.0/(DT*THETA)
       LEASQR=1.0
       RGAMMA=1.0
    ENDIF
    IF((DISCRA.EQ.7).OR.(DISCRA.EQ.125)) THEN
       !     Galerkin THETA-method. 
       GALERK=1.0
       AHAT1 =-DT*(3.*THETA-2.0)
       L1    =-3.+6.*THETA
       RGAMMA=0.
    ENDIF
    IF(DISCRA.EQ.8) THEN
       !     Space weight method. 
       THETA1=0.
       THETA2=1.
       THETA3=DT*THETA*3./2.
       AHAT1=0.
       UPWIND=.TRUE.
    ENDIF

    IF((DISCRA.GE.42).AND.(DISCRA.LE.45)) THEN
       !     LES METHOD wich uses Galerkin THETA-method. 
       GALERK=1.0
       AHAT1 =-DT*(3.*THETA-2.0)
       L1    =-3.+6.*THETA
       RGAMMA=0.
    ENDIF

    do GI=1,NGI! Was loop 
       AHAT(GI) =AHAT1
       GAMMA(GI)=RGAMMA
       GAMMAG(GI)=0.0
    END DO

    IF(DISOPT.EQ.-7) THEN
       !     Upwind galerkin...
       UPWIND=.TRUE.
       GALWIN=.TRUE.
       OPWIND=.FALSE.
    ELSE
       GALWIN=.FALSE.
    ENDIF
    !     *********************************************

    THTWDT=THETA3*DT*2./3.

    do  ELE=1,TOTELE! Was loop 340

       !     Zero diffusivity (it is calculated in loop 331)
       MUGI=0.

       !     Saner, nicer, happier coordinate transform calculation.
       if (modern_diffusivity) then
          X_ele(X_,:)=X(xondgl((ele-1)*nloc + 1:ele*nloc))         
          X_ele(Y_,:)=Y(xondgl((ele-1)*nloc + 1:ele*nloc))
          X_ele(Z_,:)=Z(xondgl((ele-1)*nloc + 1:ele*nloc))

          ! all shape functions were assumed equal above  
!!$          call transform_to_physical(X=X_ele,&
!!$               x_shape=tracer_shape,   &
!!$               m=density_shape,   &
!!$               dm_t=density_dshape)
          FLAbort("Gent McWilliams not currently supported")
          tracer_dshape=density_dshape
       end if


       do  GI=1,NGI! Was loop 331

          AGI=0.
          BGI=0.
          CGI=0.

          DGI=0.
          EGI=0.
          FGI=0.

          GGI=0.
          HGI=0.
          KGI=0.

          UD(GI)=0.
          VD(GI)=0.
          WD(GI)=0.

          DENGI(GI)=0.

          L2GIXX(GI)=0.
          L2GIXY(GI)=0.
          L2GIXZ(GI)=0.
          L2GIYY(GI)=0.
          L2GIYZ(GI)=0.
          L2GIZZ(GI)=0.

          do  L=1,NLOC! Was loop 79
             IGLS=SONDGL((ELE-1)*NLOC+L)
             IGLV=VONDGL((ELE-1)*NLOC+L)
             IGLX=XONDGL((ELE-1)*NLOC+L)
             AGI=AGI+NLX(L,GI)*X(IGLX) 
             BGI=BGI+NLX(L,GI)*Y(IGLX) 
             CGI=CGI+NLX(L,GI)*Z(IGLX) 

             DGI=DGI+NLY(L,GI)*X(IGLX) 
             EGI=EGI+NLY(L,GI)*Y(IGLX) 
             FGI=FGI+NLY(L,GI)*Z(IGLX) 

             GGI=GGI+NLZ(L,GI)*X(IGLX) 
             HGI=HGI+NLZ(L,GI)*Y(IGLX) 
             KGI=KGI+NLZ(L,GI)*Z(IGLX) 

             !     Introduce grid vels in non-linear terms. 
             IF((NDISOT.LT.10).AND.(DISCRA.NE.125)) THEN
                UD(GI)=UD(GI) + N(L,GI)*(NU(IGLV)-UG(IGLV))
                VD(GI)=VD(GI) + N(L,GI)*(NV(IGLV)-VG(IGLV))
                WD(GI)=WD(GI) + N(L,GI)*(NW(IGLV)-WG(IGLV))
             ENDIF
             DENGI(GI)=DENGI(GI)+ N(L,GI)*DENPT(IGLS)
             !     The diffusivities at Gauss pts. 
             IF((.NOT.LES).OR.(VERSIO.GT.20)) THEN
                !     Note that MUPT is a symmetric tensor while MUGI is generally asymmetric
                MUGI(X_,X_,GI) = MUGI(X_,X_,GI) + N(L,GI)*MUPTXX(IGLS)
                MUGI(X_,Y_,GI) = MUGI(X_,Y_,GI) + N(L,GI)*MUPTXY(IGLS)
                MUGI(X_,Z_,GI) = MUGI(X_,Z_,GI) + N(L,GI)*MUPTXZ(IGLS)
                MUGI(Y_,X_,GI) = MUGI(Y_,X_,GI) + N(L,GI)*MUPTXY(IGLS)
                MUGI(Y_,Y_,GI) = MUGI(Y_,Y_,GI) + N(L,GI)*MUPTYY(IGLS)
                MUGI(Y_,Z_,GI) = MUGI(Y_,Z_,GI) + N(L,GI)*MUPTYZ(IGLS)
                MUGI(Z_,X_,GI) = MUGI(Z_,X_,GI) + N(L,GI)*MUPTXZ(IGLS)
                MUGI(Z_,Y_,GI) = MUGI(Z_,Y_,GI) + N(L,GI)*MUPTYZ(IGLS)
                MUGI(Z_,Z_,GI) = MUGI(Z_,Z_,GI) + N(L,GI)*MUPTZZ(IGLS)
             ENDIF

             !     The LES length scales at Gauss pts. 
             IF(LES) THEN
                L2GIXX(GI) = L2GIXX(GI) + N(L,GI)*ML2MXX(IGLS)
                L2GIXY(GI) = L2GIXY(GI) + N(L,GI)*ML2MXY(IGLS)
                L2GIXZ(GI) = L2GIXZ(GI) + N(L,GI)*ML2MXZ(IGLS)
                L2GIYY(GI) = L2GIYY(GI) + N(L,GI)*ML2MYY(IGLS)
                L2GIYZ(GI) = L2GIYZ(GI) + N(L,GI)*ML2MYZ(IGLS)
                L2GIZZ(GI) = L2GIZZ(GI) + N(L,GI)*ML2MZZ(IGLS)
             ENDIF
          END DO

          DETJ=AGI*(EGI*KGI-FGI*HGI)&
               -BGI*(DGI*KGI-FGI*GGI)&
               +CGI*(DGI*HGI-EGI*GGI)

          DETWEI(GI)=abs(DETJ)*WEIGHT(GI)

          !     For coefficient in the inverse mat of the jacobian. 
          A11= (EGI*KGI-FGI*HGI) /DETJ
          A21=-(DGI*KGI-FGI*GGI) /DETJ
          A31= (DGI*HGI-EGI*GGI) /DETJ

          A12=-(BGI*KGI-CGI*HGI) /DETJ
          A22= (AGI*KGI-CGI*GGI) /DETJ
          A32=-(AGI*HGI-BGI*GGI) /DETJ

          A13= (BGI*FGI-CGI*EGI) /DETJ
          A23=-(AGI*FGI-CGI*DGI) /DETJ
          A33= (AGI*EGI-BGI*DGI) /DETJ

          R=0.
          ABDGI(GI)=0.
          HXGI=0.
          HYGI=0.
          HZGI=0.
          GIVOLX=0.
          GIVOLY=0.
          GIVOLZ=0.
          do L=1,NLOC! Was loop 
             NX(L,GI)= A11*NLX(L,GI)+A12*NLY(L,GI)+A13*NLZ(L,GI)
             NY(L,GI)= A21*NLX(L,GI)+A22*NLY(L,GI)+A23*NLZ(L,GI)
             NZ(L,GI)=(A31*NLX(L,GI)+A32*NLY(L,GI)+A33*NLZ(L,GI))&
                  *INDD0
             HXGI=HXGI+ABS(NX(L,GI))
             HYGI=HYGI+ABS(NY(L,GI))
             HZGI=HZGI+ABS(NZ(L,GI))
             !     N.B.   BECTWE(GI) = BETA*CTY(GI)*DETWEI(GI) . 
             IGLS=SONDGL((ELE-1)*NLOC+L)
             IGLV=VONDGL((ELE-1)*NLOC+L) 
             IF(NDISOT.LT.10) THEN
                R=R+BETA*(NX(L,GI)*NU(IGLV)+NY(L,GI)*NV(IGLV)&
                     + NZ(L,GI)*NW(IGLV))*DENGI(GI)
             ENDIF

             R=R+N(L,GI)*ABSORB(IGLS)*(1.-ABLUMP)

             ! THIS IS A FIX SO THAT ABSORPTION LUMPING DONE CONSISTENTLY WITH SOURCE TERM 
             ABDGI(GI)=ABDGI(GI) + N(L,GI)*ABSORB(IGLS)*ABLUMP
          END DO

          BECTWE(GI)=R

          DENGI2(GI)=GALERK + (1.0-GALERK)*DENGI(GI)
          L1GI(GI)=L1+LEASQR*BECTWE(GI)

          IF(LES) THEN
             !     Amend diffusivity for LES
             CALL LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
                  NX,NY,NZ, &
                  NU,NV,NW, UG,VG,WG, &
                  L2GIXX,L2GIXY,L2GIXZ,&
                  L2GIYY,L2GIYZ,L2GIZZ)
             !     Add in the normal viscosity...
             !     Note that L2GI is a symmetric tensor while MUGI is generally asymmetric
             MUGI(X_,X_,GI)=MUGI(X_,X_,GI)+L2GIXX(GI)
             MUGI(X_,Y_,GI)=MUGI(X_,Y_,GI)+L2GIXY(GI)
             MUGI(X_,Z_,GI)=MUGI(X_,Z_,GI)+L2GIXZ(GI)
             MUGI(Y_,X_,GI)=MUGI(Y_,X_,GI)+L2GIXY(GI)
             MUGI(Y_,Y_,GI)=MUGI(Y_,Y_,GI)+L2GIYY(GI)
             MUGI(Y_,Z_,GI)=MUGI(Y_,Z_,GI)+L2GIYZ(GI)
             MUGI(Z_,X_,GI)=MUGI(Z_,X_,GI)+L2GIXZ(GI)
             MUGI(Z_,Y_,GI)=MUGI(Z_,Y_,GI)+L2GIYZ(GI)
             MUGI(Z_,Z_,GI)=MUGI(Z_,Z_,GI)+L2GIZZ(GI)
             !           if(les)
          ENDIF

          !     Insert diffusivity from diffusivity module. Note that this overrides
          !     all other diffusivities.
          if (modern_diffusivity) then

             MUGI(:,:,gi)=get_diffusivity( &
                  tracer_shape=tracer_shape,&
                  density=ele_val(density,ele),&
                  density_shape=density_shape, &
                  density_dshape=density_dshape,&
                  topdis=ele_val_at_quad(topdis,ele),&
                  botdis=ele_val_at_quad(botdis,ele),&
                  d=ele_val_at_quad(d,ele),&
                  gi=gi)

          end if

          !     Work out GAMMA at each Gauss pt(OR centre of element)
          IF(UPWIND) THEN 
             HXGI=ABS(A11*UD(GI)+A12*VD(GI)+A13*WD(GI))
             HYGI=ABS(A21*UD(GI)+A22*VD(GI)+A23*WD(GI))
             HZGI=ABS(A31*UD(GI)+A32*VD(GI)+A33*WD(GI))

             !     HX,HY are the characteristic length scales in x,y directions. 
             RUD=MAX(ABS(UD(GI)),1E-7)
             RVD=MAX(ABS(VD(GI)),1E-7)
             RWD=MAX(ABS(WD(GI)),1E-7)
             !     If  XTPETV  then (x,y,t) Petrov Galerkin Method  else (x,y).
             IF(XTPETV) THEN
                HOVERQ=MIN(2./MAX(HXGI,HYGI,HZGI,1.E-7), DT)
                !              if(xtpetv)
             ELSE
                HOVERQ=2./MAX(HXGI,HYGI,HZGI,1.E-7)
                !              if(xtpetv)
                !              else
             ENDIF
             IF(GALWIN) THEN
                GAMMAG(GI)=RWIND*ALPHA*HOVERQ 
                !              if(galwin)
             ELSE
                IF(OPWIND) THEN 
                   !     (Near optimal choice of upwind parameter) PE=grid Peclet no. 
                   !     The PE number here does not take into account anisotropic MU.
                   PE=0.5*DENGI(GI)*(RUD**2+RVD**2+RWD**2)*HOVERQ&
                        / MAX(MUGI(X_,X_,GI),1.E-7) 
                   ALPHA=1.0-1.0/MAX(1.0,PE)
                   !                 if(opwind)
                ENDIF
                GAMMA(GI)=RWIND*ALPHA*HOVERQ 
                AHAT(GI) =AHAT1*GAMMA(GI)
                !              if(galwin)
                !              else
             ENDIF
             !           if(upwind)
          ENDIF
          !     NB.  THTWDT=THETA3*DT*2./3.     
          MASCON(GI) = THETA1*2.*AHAT(GI)*INDT*DENGI(GI)&
               +THETA2*(AHAT(GI)*BECTWE(GI)+L1GI(GI)*DENGI(GI))&
               +THTWDT*L1GI(GI)*BECTWE(GI)     
       END DO

       do  ILOC=1,NLOC! Was loop 350
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          GLOBIS=SONDGL((ELE-1)*NLOC+ILOC)
          do  JLOC=1,NLOC! Was loop 360
             GLOBJ =NDGLNO((ELE-1)*NLOC+JLOC)
             GLOBJS=SONDGL((ELE-1)*NLOC+JLOC)

             MASS=0.
             MATRIX=0.
             RHSMAT=0.
             SOUMAT=0.  
             VLM = 0.   
             VABS = 0.

             do  GI=1,NGI! Was loop 458
                R1=UD(GI)*NX(ILOC,GI)&
                     +VD(GI)*NY(ILOC,GI)&
                     +WD(GI)*NZ(ILOC,GI)
                R2=(UD(GI)*NX(JLOC,GI)+VD(GI)*NY(JLOC,GI) &
                     + WD(GI)*NZ(JLOC,GI)  )*DENGI(GI)

                !     N.B.   BECTWE(GI) = BETA*CTY(GI)*DENGI(GI) .   

                RNN  =N(ILOC,GI)*N(JLOC,GI)
                VLMGI=RNN*DETWEI(GI)
                VLKGI= NX(ILOC,GI)*(NX(JLOC,GI)*MUGI(X_,X_,GI)&
                     +NY(JLOC,GI)*MUGI(X_,Y_,GI)+NZ(JLOC,GI)*MUGI(X_,Z_,GI))&
                     +NY(ILOC,GI)*(NX(JLOC,GI)*MUGI(Y_,X_,GI)&
                     +NY(JLOC,GI)*MUGI(Y_,Y_,GI)+NZ(JLOC,GI)*MUGI(Y_,Z_,GI))&
                     +NZ(ILOC,GI)*(NX(JLOC,GI)*MUGI(Z_,X_,GI)&
                     +NY(JLOC,GI)*MUGI(Z_,Y_,GI)+NZ(JLOC,GI)*MUGI(Z_,Z_,GI))

                VF1M =AHAT(GI)*(N(ILOC,GI)*R2*RSYM + VLKGI)

                VF2  =GAMMA(GI)*R1*N(JLOC,GI)  
                VF2M =VF2*DENGI(GI)

                RR =GAMMA(GI)*R1*(R2+BECTWE(GI)*N(JLOC,GI))&
                     +L1GI(GI)*N(ILOC,GI)*R2
                VF3M=RR*RSYM + L1GI(GI)*VLKGI
                !     VF3 =RR      + L1GI(GI)*VLKGI
                !     All the mass matrix contributions...       
                MASS = MASS + MASCON(GI)*VLMGI
                !     Matrix contributions (Excluding mass matrix)...
                MATRIX=MATRIX +  (THETA2*(VF1M + VF2M*RSYM)&
                     +THTWDT*VF3M)*DETWEI(GI) &
                                !     simple upwinding for use with Galerkin method...
                     +DT*THETA*GAMMAG(GI)*R1*R2*DETWEI(GI)
                !     NB.  THTWDT=THETA3*DT*2./3.
                !     RHSMAT contributions.
                RHSMAT=RHSMAT +( (2.*AHAT(GI)*INDT + L1GI(GI))&
                     *(BECTWE(GI)*RNN + N(ILOC,GI)*R2 + VLKGI)&
                     +GAMMA(GI)*R1*(BECTWE(GI)*N(JLOC,GI)+R2) )*DETWEI(GI)&
                                !     simple upwinding for use with Galerkin method...
                     +GAMMAG(GI)*R1*R2*DETWEI(GI)
                !     Source Material. 
                SOUMAT=SOUMAT + (2.*INDT*AHAT(GI) + L1GI(GI))*VLMGI &
                     + VF2*DETWEI(GI)

                VLM  = VLM + RNN*DETWEI(GI)
             END DO

             VABS = VLM*ABSORB(GLOBIS)*ABLUMP

             ML(GLOBI)=ML(GLOBI)+MASS

             VECX(GLOBI)=VECX(GLOBI) - RHSMAT*U(GLOBJ) &
                  +SOUMAT*SOURCX(GLOBJS)*(1.-RSLUMP)&
                  +SOUMAT*SOURCX(GLOBIS)*RSLUMP &
                  -VABS*U(GLOBI)*ABLUMP    

             !     we will now place contributions into the matrices BIGM and BIGT
             !     these are can be symm matrices.
             !     Find count. ***************************************
             IF(INMAT) THEN
                !     lump the mass matrix.(LUMP=1. if we are to LUMP the mass matrix)
                BIGM(CENTRM(GLOBI))=BIGM(CENTRM(GLOBI))+MASS*RLUMP &
                     + DT*VABS*ABLUMP ! Fully implicit
                LOWER=FINDRM(GLOBI) 
                UPPER=FINDRM(GLOBI+1)-1
7000            CONTINUE
                INUM=LOWER+(UPPER-LOWER+1)/2 
                IF(GLOBJ.GE.COLM(INUM) )  THEN 
                   LOWER=INUM
                ELSE
                   UPPER=INUM
                ENDIF
                IF(UPPER-LOWER.LE.1) THEN
                   IF(GLOBJ.EQ.COLM(LOWER)) THEN
                      COUNTER=LOWER
                   ELSE
                      COUNTER=UPPER
                   ENDIF
                   GOTO 9000
                ENDIF
                GOTO 7000
9000            CONTINUE
                !     Find count. ***************************************

                BIGM(COUNTER)=BIGM(COUNTER)+MATRIX + MASS*RCONSI

                !     end of IF(INMAT) THEN
             ENDIF

          end do ! Was loop 360

          !     *** THE FOLLOWING FORMS THE MATRICES C1T AND C2T ***
          IF(GETC12) THEN

             do  JLOC=1,MLOC! Was loop 123

                PGLOBJ=PNDGLN( (ELE-1)*MLOC +JLOC)
                !     Find count. ***************************************
                LOWER=FINDCT(PGLOBJ) 
                UPPER=FINDCT(PGLOBJ+1)-1
7001            CONTINUE
                INUM=LOWER+(UPPER-LOWER+1)/2 
                IF(GLOBI.GE.COLCT(INUM) )  THEN 
                   LOWER=INUM
                ELSE
                   UPPER=INUM
                ENDIF
                IF(UPPER-LOWER.LE.1) THEN
                   IF(GLOBI.EQ.COLCT(LOWER)) THEN
                      POS=LOWER
                   ELSE
                      POS=UPPER
                   ENDIF
                   GOTO 9001
                ENDIF
                GOTO 7001
9001            CONTINUE
                !     Find count. ***************************************
                 do GI=1,NGI! Was loop 
                    PWEIGI=2.*AHAT(GI)*INDT+L1GI(GI)
                    C1T(POS)=C1T(POS)+PWEIGI*M(JLOC,GI)*NX(ILOC,GI)*DETWEI(GI)
                    C2T(POS)=C2T(POS)+PWEIGI*M(JLOC,GI)*NY(ILOC,GI)*DETWEI(GI)
                    C3T(POS)=C3T(POS)+PWEIGI*M(JLOC,GI)*NZ(ILOC,GI)*DETWEI(GI)
                 END DO
             end do ! Was loop 123

          ENDIF
       end do ! Was loop 350
    end do ! Was loop 340
    CALL PMINMX(MUPTXX,NONODS,'******MUPTXX  ')
    CALL PMINMX(MUPTXY,NONODS,'******MUPTXY  ')
    CALL PMINMX(MUPTXZ,NONODS,'******MUPTXZ  ')
    CALL PMINMX(MUPTYY,NONODS,'******MUPTYY  ')
    CALL PMINMX(MUPTYZ,NONODS,'******MUPTYZ  ')
    CALL PMINMX(MUPTZZ,NONODS,'******MUPTZZ  ')
    CALL PMINMX(DENPT,NONODS,'******DENPT  ')

    IF(MLCENT) THEN
       do I=1,NONODS! Was loop 
          ML(I)=BIGM(CENTRM(I))
       END DO
    ENDIF
    EWRITE(1,*) 'EXITING HART3D()'
  END SUBROUTINE HART3D

  SUBROUTINE HART3DSIM(U,&
       NU,NV,NW,UG,VG,WG,&
       SOURCX,X,Y,Z,&
       VECX, &
       BIGM, ML,&
       TOTELE, &
       FINDRM,COLM,NCOLM,NONODS,CENTRM, &
       N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
       DISOPT2,DISOPN,DT,THETA,LUMP,&
       ABSLUM,SLUMP,&
       NDGLNO, &
       SONDGL,SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
       MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
       ABSORB,ISPHERE,&
       nnodp,para,halo_tag)

    !     This subroutine discretises a field equation using Galerkin Least squares. 
    !     The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
    !     
    ! All methods use a balancing diffusion term and THETA time stepping. 
    ! The following are absolute values of DISOPT...
    ! Galerkin (146 is defaulted to if non of the below are used): 
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

    !     if DISOPN.ne.0 then set advection to zero and treat advection 
    !     using the high res method.
    !     IF LUMP then lump the mass matrix else dont. 
    !     ABSORB =absorbtion of energy etc. 
    !     MUPTXX,MUPTYY =components of diffusivities. 
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
    !     --------------------------------------------------------------------
    !     If we are solving for another variable like temperature 
    !     or chemical species then NU,NV,NW will be the velocities 
    !     and U=TEMPERATURE OR CHEMICAL SPECIES. 
    !     --------------------------------------------------------------------
    !     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
    REAL RWIND
    PARAMETER(RWIND=0.5)
    INTEGER SNONOD,VNONOD,XNONOD,NCOLM,NLOC,NGI
    INTEGER DISOPT2,DISOPT,DISOPN
    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
    REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
    REAL ABSORB(SNONOD)
    REAL  DT,THETA
    INTEGER TOTELE,NONODS
    !     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    !     If PHI=1.0 then bakward Euler is used in NS equns.
    !     Similarly for the temperature equation except width variable THETA.
    !     have gravity in -ve y direction.
    REAL U(NONODS)
    REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
    REAL::VLM=0.0
    REAL UD(NGI),VD(NGI),WD(NGI)
    REAL UDDX(NGI),VDDY(NGI),WDDZ(NGI)
    REAL UDDX2(NGI),VDDY2(NGI),WDDZ2(NGI)
    REAL MUGIXX,MUGIXY,MUGIXZ
    REAL MUGIYY,MUGIYZ,MUGIZZ
    REAL L2GIXX,L2GIXY,L2GIXZ
    REAL L2GIYY,L2GIYZ,L2GIZZ
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    INTEGER nnodp,para,halo_tag
    ! HX,HY-characteristic length scales in x,y directions.
    REAL HXGI,HYGI,HZGI

    REAL VECX(NONODS)
    REAL BIGM(NCOLM)
    REAL ML(NONODS)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL SOURCX(SNONOD)
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER SONDGL(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER ELE,ILOC,JLOC,L
    INTEGER GI,COUNT,I,GLOBI,GLOBJ,GLOBIS,GLOBJS
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    INTEGER CENTRM(NONODS)

    REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL LUMP
    LOGICAL ABSLUM,SLUMP
    REAL ABLUMP,RSLUMP
    REAL::RNN=0.0
    REAL RLUMP,RCONSI
    INTEGER IGLV
    INTEGER ISPHERE,ICENT
    REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
    REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
    REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
    REAL MYZDTW(NGI),MZZDTW(NGI)
    REAL MUPTXXGI,MUPTXYGI,MUPTXZGI,MUPTYYGI,MUPTYZGI,MUPTZZGI
    REAL XD(NGI),YD(NGI),ZD(NGI)
    REAL VLND
    REAL HINX,HINY,HINZ,VLKTEN,LHNLN,LHNN,NLN,VLMAB
    REAL BETBAL,ABGITOT
    REAL XNDMUJLOC,YNDMUJLOC,ZNDMUJLOC,XNDMUILOC,YNDMUILOC,ZNDMUILOC
    REAL TEN,TENAB,ABTEN
    REAL VLN,RDIFF
    real::ADVILOC=0.0, ADVJLOC=0.0
    REAL UATTHETA,VATTHETA,WATTHETA
    LOGICAL BALDIF,BALGLS
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI) 
    REAL GAMB(NGI),GAMBAL(NGI)
    REAL ABGI(NGI),ABSDX(NGI),ABSDY(NGI),ABSDZ(NGI),DDETWE(NGI)
    LOGICAL GETHIGH
    INTEGER LESNOD
    REAL RINMAT,RHSMAT,SRHSMAT
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

    REAL, ALLOCATABLE, DIMENSION(:)::HIMATX
    REAL, ALLOCATABLE, DIMENSION(:)::HIMATY
    REAL, ALLOCATABLE, DIMENSION(:)::HIMATZ
    REAL, ALLOCATABLE, DIMENSION(:)::HIMASS

    REAL, ALLOCATABLE, DIMENSION(:)::VECU
    REAL, ALLOCATABLE, DIMENSION(:)::VECV
    REAL, ALLOCATABLE, DIMENSION(:)::VECW

    REAL, ALLOCATABLE, DIMENSION(:)::VECYT

    REAL, ALLOCATABLE, DIMENSION(:,:,:)::EV

    ewrite(1, *) "just in SUBROUTINE HART3DSIM()"
    ewrite(2, *) "isphere=",isphere

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

    do  I=1,NONODS
       ML(I)=0.
       VECX(I)=0.
    end do
    BIGM(1:NCOLM) = 0.0

    RLUMP=1.
    IF(LUMP) RLUMP=1.
    RCONSI=1.-RLUMP

    ABLUMP=1.
    RSLUMP=1.
    IF(ABSLUM) ABLUMP=1.
    IF(SLUMP)  RSLUMP=1.


    element_loop: do  ELE=1,TOTELE
       !     
       ! Get determinant and derivatives...
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NCLOC,NGI, &
            N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
            NX,NY,NZ, NCX,NCY,NCZ,&
            A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            XD,YD,ZD,&
            ISPHERE) 


       gauss_point_loop1: do GI=1, NGI

          UD(GI)=0.    
          VD(GI)=0.
          WD(GI)=0.
          UDDX(GI)=0.    
          VDDY(GI)=0.
          WDDZ(GI)=0.
          UDDX2(GI)=0.    
          VDDY2(GI)=0.
          WDDZ2(GI)=0.
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

          do L=1,NLOC
             IGLV=VONDGL((ELE-1)*NLOC+L)

             IF(DISOPN.EQ.0) THEN
                ! if not using control volumes:
                UD(GI)=UD(GI) + N(L,GI)*(NU(IGLV)-UG(IGLV))
                VD(GI)=VD(GI) + N(L,GI)*(NV(IGLV)-VG(IGLV))
                WD(GI)=WD(GI) + N(L,GI)*(NW(IGLV)-WG(IGLV))

                UDDX(GI)=UDDX(GI) + NX(L,GI)*(NU(IGLV)-UG(IGLV))
                VDDY(GI)=VDDY(GI) + NY(L,GI)*(NV(IGLV)-VG(IGLV))
                WDDZ(GI)=WDDZ(GI) + NZ(L,GI)*(NW(IGLV)-WG(IGLV))

                UDDX2(GI)=UDDX2(GI) + NX(L,GI)*NU(IGLV)
                VDDY2(GI)=VDDY2(GI) + NY(L,GI)*NV(IGLV)
                WDDZ2(GI)=WDDZ2(GI) + NZ(L,GI)*NW(IGLV)
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
          end do
          !
          ! calculate diffusivity tensor:

          IF((DISOPT.EQ.147).or.(DISOPT.EQ.150).or.(DISOPT.EQ.153)&
               .or.(DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154)) THEN
             ! if applying LES:

              ! based on anisotropic length scales...
              CALL ANO_OCEAN_SCAL(&
                   HIMXXDTW(GI),HIMXYDTW(GI),HIMXZDTW(GI),&
                   HIMYYDTW(GI),HIMYZDTW(GI),HIMZZDTW(GI),&
                   MXXDTW(GI),MXYDTW(GI),MXZDTW(GI),&
                   MYYDTW(GI),MYZDTW(GI),MZZDTW(GI),&
                   XD(GI),YD(GI),ZD(GI),&
                   ELE,NLOC,TOTELE,XONDGL,XNONOD,X,Y,Z,ISPHERE) 

              CALL HALF_LESVIS_HART(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
                   NX,NY,NZ, &
                   NU,NV,NW, UG,VG,WG, &
                   HIMXXDTW,HIMXYDTW,HIMXZDTW,&
                   HIMYYDTW,HIMYZDTW,HIMZZDTW,&
                   MXXDTW,MXYDTW,MXZDTW,&
                   MYYDTW,MYZDTW,MZZDTW)

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
             ! no LES:
             MUGIXX=0.0
             MUGIXY=0.0
             MUGIXZ=0.0
             MUGIYY=0.0
             MUGIYZ=0.0
             MUGIZZ=0.0
          ENDIF
          MXXDTW(GI)=(MUGIXX+MUPTXXGI)*DETWEI(GI)
          MXYDTW(GI)=(MUGIXY+MUPTXYGI)*DETWEI(GI)
          MXZDTW(GI)=(MUGIXZ+MUPTXZGI)*DETWEI(GI)

          MYYDTW(GI)=(MUGIYY+MUPTYYGI)*DETWEI(GI)
          MYZDTW(GI)=(MUGIYZ+MUPTYZGI)*DETWEI(GI)

          MZZDTW(GI)=(MUGIZZ+MUPTZZGI)*DETWEI(GI)

       end do gauss_point_loop1

       iloc_loop: do ILOC=1,NLOC
          GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
          GLOBIS=SONDGL((ELE-1)*NLOC+ILOC)
          jloc_loop: do  JLOC=1,NLOC
             GLOBJ =NDGLNO((ELE-1)*NLOC+JLOC)
             GLOBJS=SONDGL((ELE-1)*NLOC+JLOC)

             VLM=0. 

             HINX=0.0
             HINY=0.0
             HINZ=0.0

             VLKTEN=0.0

             LHNLN=0.0

             LHNN=0.0

             NLN=0.0

             VLND=0.0

             VLMAB=0.0


             RINMAT=0.0
             RHSMAT =0.0
             SRHSMAT=0.0

             gauss_point_loop2: do  GI=1,NGI

                BALDIF=.TRUE.
                BALGLS=.TRUE.

                IF((DISOPT.GE.149).AND.(DISOPT.LE.151)) THEN
                   ! Galerkin Least Squares BETBAL=0.0....
                   HXGI=ABS(A11(GI)*UD(GI)+A12(GI)*VD(GI)+A13(GI)*WD(GI))
                   HYGI=ABS(A21(GI)*UD(GI)+A22(GI)*VD(GI)+A23(GI)*WD(GI))
                   HZGI=ABS(A31(GI)*UD(GI)+A32(GI)*VD(GI)+A33(GI)*WD(GI))
                   BETBAL=0.0
                   GAMB(GI)=MIN(2./MAX(HXGI,HYGI,HZGI,1.E-7), DT)
                ELSE IF((DISOPT.GE.152).AND.(DISOPT.LE.154)) THEN
                   !     Galerkin Least Squares 
                   HXGI=ABS(A11(GI)*UD(GI)+A12(GI)*VD(GI)+A13(GI)*WD(GI))
                   HYGI=ABS(A21(GI)*UD(GI)+A22(GI)*VD(GI)+A23(GI)*WD(GI))
                   HZGI=ABS(A31(GI)*UD(GI)+A32(GI)*VD(GI)+A33(GI)*WD(GI))
                   BETBAL=1.0
                   GAMB(GI)=MIN(2./MAX(HXGI,HYGI,HZGI,1.E-7), DT)
                ELSE
                   ! Galerkin....
                   BETBAL=0.0
                   GAMB(GI)=0.0
                   BALDIF=.FALSE.
                   BALGLS=.FALSE.
                ENDIF

                GAMBAL(GI)=GAMB(GI)  / (1.0+GAMB(GI)*BETBAL/DT)

                ABGITOT=ABGI(GI) 

                VLN=N(ILOC,GI)*(UD(GI)*NX(JLOC,GI)+VD(GI)*NY(JLOC,GI)&
                     +WD(GI)*NZ(JLOC,GI))*DETWEI(GI) 
                VLND=VLND+VLN

                RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)

                VLM=VLM+RNN 

                XNDMUJLOC=(NX(JLOC,GI)*MXXDTW(GI)&
                     +NY(JLOC,GI)*MXYDTW(GI)+NZ(JLOC,GI)*MXZDTW(GI))
                YNDMUJLOC=(NX(JLOC,GI)*MXYDTW(GI)&
                     +NY(JLOC,GI)*MYYDTW(GI)+NZ(JLOC,GI)*MYZDTW(GI))
                ZNDMUJLOC=(NX(JLOC,GI)*MXZDTW(GI)&
                     +NY(JLOC,GI)*MYZDTW(GI)+NZ(JLOC,GI)*MZZDTW(GI))

                XNDMUILOC=(NX(ILOC,GI)*MXXDTW(GI)&
                     +NY(ILOC,GI)*MXYDTW(GI)+NZ(ILOC,GI)*MXZDTW(GI))
                YNDMUILOC=(NX(JLOC,GI)*MXYDTW(GI)&
                     +NY(ILOC,GI)*MYYDTW(GI)+NZ(ILOC,GI)*MYZDTW(GI))
                ZNDMUILOC=(NX(ILOC,GI)*MXZDTW(GI)&
                     +NY(ILOC,GI)*MYZDTW(GI)+NZ(ILOC,GI)*MZZDTW(GI))

                TEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC+NZ(ILOC,GI)*ZNDMUJLOC

                TENAB=     +N(ILOC,GI)*(ABSDX(GI)*XNDMUJLOC+ABSDY(GI)*YNDMUJLOC&
                     +ABSDZ(GI)*ZNDMUJLOC)  +ABGITOT*TEN

                ABTEN=     +(ABSDX(GI)*XNDMUILOC+ABSDY(GI)*YNDMUILOC&
                     +ABSDZ(GI)*ZNDMUILOC)*N(JLOC,GI)  +ABGITOT*TEN


                VLKTEN=VLKTEN+TEN 

                IF((DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154)) THEN
                   ! Subtract out extended stencil 2nd order operator to get a 4th order operator...
                   VECX(GLOBI)=VECX(GLOBI) + NX(ILOC,GI)*N(JLOC,GI)*( &
                        MXXDTW(GI)*VEDUDX(GLOBJ) +MXYDTW(GI)*VEDUDY(GLOBJ) +MXZDTW(GI)*VEDUDZ(GLOBJ) )  &
                        + NY(ILOC,GI)*N(JLOC,GI)*( &
                        MXYDTW(GI)*VEDUDX(GLOBJ) +MYYDTW(GI)*VEDUDY(GLOBJ) +MYZDTW(GI)*VEDUDZ(GLOBJ) ) &
                        + NZ(ILOC,GI)*N(JLOC,GI)*( &
                        MXZDTW(GI)*VEDUDX(GLOBJ) +MYZDTW(GI)*VEDUDY(GLOBJ) +MZZDTW(GI)*VEDUDZ(GLOBJ) ) 


                ENDIF

                IF(BALDIF) THEN
                   ADVILOC=UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI)+WD(GI)*NZ(ILOC,GI)
                   ADVJLOC=UD(GI)*NX(JLOC,GI)+VD(GI)*NY(JLOC,GI)+WD(GI)*NZ(JLOC,GI)
                   RDIFF=ADVILOC*ADVJLOC*DETWEI(GI) 

                   LHNLN=LHNLN+GAMBAL(GI)*RDIFF  &
                        +GAMBAL(GI)*ABGITOT*( ABGITOT*RNN + N(ILOC,GI)*ADVJLOC*DETWEI(GI) &
                        + ADVILOC*N(JLOC,GI)*DETWEI(GI)  )&
                        + GAMBAL(GI)*(ABTEN +TENAB)
                ENDIF


                IF(BALGLS) THEN
                   if(.not.BALDIF) then
                      FLAbort("bug - ADVILOC has not been set")
                   else
                      LHNN=LHNN+GAMBAL(GI)*(ADVILOC*N(JLOC,GI)*DETWEI(GI) +TEN +ABGITOT*RNN)
                   end if
                ENDIF


                NLN=NLN + (VLN+TEN +ABGI(GI)*RNN) 

                VLMAB=VLMAB+ABGI(GI)*RNN

             end do gauss_point_loop2


             ! Find count. ***************************************
             CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
             ! lump the mass matrix.
             ICENT=CENTRM(GLOBI)


             ! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...
              VECX(GLOBI)=VECX(GLOBI) -(NLN+THETA*LHNLN)*U(GLOBJ)

              IF(SLUMP) THEN
                 VECX(GLOBI)=VECX(GLOBI) +(VLM+THETA*LHNN)*SOURCX(GLOBI) 
              ELSE
                 VECX(GLOBI)=VECX(GLOBI) +(VLM+THETA*LHNN)*SOURCX(GLOBJ) 
              ENDIF


              ! lump the mass matrix.(LUMP=1. if we are to LUMP the mass matrix)
              IF(LUMP) THEN
                 BIGM(ICENT)=BIGM(ICENT)  +VLM  
                 BIGM(COUNT)=BIGM(COUNT)       +THETA*(DT*NLN+LHNN+DT*THETA*LHNLN)
              ELSE
                 ! Balancing tensor form...
                 BIGM(COUNT)=BIGM(COUNT)  +VLM +THETA*(DT*NLN+LHNN+DT*THETA*LHNLN)
              ENDIF

             ! Lump the absorption to put in to pressure. 
             ML(GLOBI)=ML(GLOBI) + VLM 

          end do jloc_loop
       end do iloc_loop
    end do element_loop

  END SUBROUTINE HART3DSIM

  SUBROUTINE HART3DSGS(U,&
       NU,NV,NW,UG,VG,WG,&
       SOURCX,X,Y,Z,&
       VECX, &
       BIGM, ML,&
       TOTELE, &
       FINDRM,COLM,NCOLM,NONODS, &
       N,NLX,NLY,NLZ, WEIGHT, NLOC,NGI,&
       DISOPT2,DISOPN,DT,THETA,BETA,LUMP,&
       ABSLUM,SLUMP,&
       NDGLNO,&
       SNONOD, VONDGL,VNONOD, XONDGL,XNONOD,&
       MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,&
       ABSORB,ISPHERE,&
       nnodp,para,halo_tag,&
       NSUBTLOC,NSUBVLOC,&
       NUSUB,NVSUB,NWSUB,&
       USUB, &
       DMATINVCSTORE,&
       DINVSOUSUBUSTORE)

    !     This subroutine discretises a field equation using Galerkin Least squares. 
    !     The discretization scheme is controlled through DISOPT(DISCRETIZATION OPTION) 
    !     
    ! All methods use a balancing diffusion term and THETA time stepping. 
    ! The following are absolute values of DISOPT...
    ! Galerkin (146 is defaulted to if non of the below are used): 
    ! 148 is LES applied to the global system only RECOMMENDED
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
    !     NU,NV,NW are for the non-linear terms ordeneraly NU=U,NV=V,NW=W.
    !     UG,VG,WG are the grid velocities. 
    !     NDGLNO=element pter for unknowns. 
    !     VONDGL=element pter for advection NU,UG.
    !     XONDGL=element pter for coordinates(x,y,z)
    !     R0=the minimum distance from ocean bed to centre of earth.
    !     D0=H/r0 where H is the height of ocean above R0 (H not used here)
    !     For dimensional run D0=0.(R0 is not used used here)
    !     NB This sub is only used for field equations in ocean modelling. 
    !     -------------------------------------------------------------------------------------------
    !     If we are solving for another variable like temperature 
    !     or chemical species then NU,NV,NW will be the velocities 
    !     and U=TEMPERATURE OR CHEMICAL SPECIES. 
    !     -------------------------------------------------------------------------------------------
    !     NB SOURCE, DENPT, viscosity are SNONOD in length and have ele-pt SONDGL.
    REAL RWIND
    PARAMETER(RWIND=0.5)
    INTEGER SNONOD,VNONOD,XNONOD,NCOLM,NLOC,NGI
    INTEGER DISOPT2,DISOPT,DISOPN
    REAL MUPTXX(SNONOD),MUPTXY(SNONOD),MUPTXZ(SNONOD)
    REAL MUPTYY(SNONOD),MUPTYZ(SNONOD),MUPTZZ(SNONOD)
    REAL ABSORB(SNONOD)
    REAL  DT,THETA,BETA
    INTEGER TOTELE,NONODS
    !     If PHI=0.5 then Crank-Nicolson is used in Navier Stokes equ's.
    !     If PHI=1.0 then bakward Euler is used in NS equns.
    !     Similarly for the temperature equation except width variable THETA.
    !     have gravity in -ve y direction.
    REAL U(NONODS)
    INTEGER NSUBTLOC,NSUBVLOC
    REAL NUSUB(TOTELE*NSUBVLOC),NVSUB(TOTELE*NSUBVLOC)
    REAL NWSUB(TOTELE*NSUBVLOC)
    REAL USUB(TOTELE*NSUBTLOC)
    REAL DMATINVCSTORE(TOTELE,NSUBTLOC,NLOC)
    REAL DINVSOUSUBUSTORE(TOTELE,NSUBTLOC)
    REAL UG(VNONOD),VG(VNONOD),WG(VNONOD)
    REAL NU(VNONOD),NV(VNONOD),NW(VNONOD)
    REAL::VLM=0.0
    REAL ALPHA,ALPHA2,ALPHA3
    REAL UD(NGI),VD(NGI),WD(NGI)
    REAL UDDX(NGI),VDDY(NGI),WDDZ(NGI)
    REAL UDDX2(NGI),VDDY2(NGI),WDDZ2(NGI)
    REAL UGDX(NGI),VGDY(NGI),WGDZ(NGI)
    REAL MUGIXX,MUGIXY,MUGIXZ
    REAL MUGIYY,MUGIYZ,MUGIZZ
    REAL L2GIXX,L2GIXY,L2GIXZ
    REAL L2GIYY,L2GIYZ,L2GIZZ
    REAL DETWEI(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    INTEGER nnodp,para,halo_tag

    REAL VECX(NONODS)
    REAL BIGM(NCOLM)
    REAL ML(NONODS)
    REAL X(XNONOD),Y(XNONOD),Z(XNONOD)
    REAL SOURCX(SNONOD)
    INTEGER NDGLNO(TOTELE*NLOC)
    INTEGER VONDGL(TOTELE*NLOC)
    INTEGER XONDGL(TOTELE*NLOC)
    INTEGER ELE,ILOC,JLOC,L
    INTEGER GI,COUNT,I,GLOBI,GLOBJ
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)

    REAL, target:: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    LOGICAL LUMP
    LOGICAL ABSLUM,SLUMP
    REAL ABLUMP,RSLUMP
    REAL::RNN=0.0
    REAL RLUMP,RCONSI
    INTEGER IGLV
    INTEGER ISPHERE
    REAL HIMXXDTW(NGI),HIMXYDTW(NGI),HIMXZDTW(NGI),HIMYYDTW(NGI)
    REAL HIMYZDTW(NGI),HIMZZDTW(NGI)
    REAL MXXDTW(NGI),MXYDTW(NGI),MXZDTW(NGI),MYYDTW(NGI)
    REAL MYZDTW(NGI),MZZDTW(NGI)
    REAL AMXXDTW(NGI),AMXYDTW(NGI),AMXZDTW(NGI),AMYYDTW(NGI)
    REAL AMYZDTW(NGI),AMZZDTW(NGI)
    REAL MUPTXXGI,MUPTXYGI,MUPTXZGI,MUPTYYGI,MUPTYZGI,MUPTZZGI
    REAL XD(NGI),YD(NGI),ZD(NGI)
    REAL VLND,VLN2,VLND2,VLN3,VLND3
    REAL VLKTEN,stVLKTEN,VLMAB
    REAL XNDMUJLOC,YNDMUJLOC,ZNDMUJLOC,XNDMUILOC,YNDMUILOC,ZNDMUILOC
    REAL TEN
    REAL VLN
    REAL A11(NGI),A12(NGI),A13(NGI)
    REAL A21(NGI),A22(NGI),A23(NGI)
    REAL A31(NGI),A32(NGI),A33(NGI)
    REAL ABGI(NGI),ABSDX(NGI),ABSDY(NGI),ABSDZ(NGI),DDETWE(NGI)
    REAL AMAT(NLOC,NLOC),BMAT(NLOC,NSUBTLOC)
    REAL CMAT(NSUBTLOC,NLOC),DMAT(NSUBTLOC,NSUBTLOC)
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
    REAL REALMASSM2M2(NSUBTLOC,NSUBTLOC)
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
    REAL DIFMATXC(NSUBTLOC,NLOC),DIFMATYC(NSUBTLOC,NLOC)
    REAL DIFMATZC(NSUBTLOC,NLOC)
    REAL STREAM_MATX(NLOC,NSUBTLOC),STREAM_MATY(NLOC,NSUBTLOC)
    REAL STREAM_MATZ(NLOC,NSUBTLOC)
    REAL VLSTREAM_MATX,VLSTREAM_MATY,VLSTREAM_MATZ
    INTEGER LESNOD,JNODSUB
    REAL VLNXN,VLNYN,VLNZN
    REAL VLNNX,VLNNY,VLNNZ
    ! Local variables...
    INTEGER SNLOC,SNGI,SNCLOC
    PARAMETER(SNLOC=3,SNCLOC=6,SNGI=7)
    REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
    REAL SWEIGH(SNGI)
    REAL SNC(SNCLOC,SNGI),SNCLX(SNCLOC,SNGI),SNCLY(SNCLOC,SNGI)
    REAL XSL(SNCLOC),YSL(SNCLOC),ZSL(SNCLOC)
    REAL SUD(SNGI),SVD(SNGI),SWD(SNGI),SABS(SNGI)
    REAL SUDGLOB(SNGI),SVDGLOB(SNGI),SWDGLOB(SNGI)
    REAL NDOTQ(SNGI),NDOTQGLOB(SNGI)
    REAL SL1(SNGI), SL2(SNGI), SL3(SNGI), SL4(SNGI)
    REAL SDETWE(SNGI)
    INTEGER IFACE,NFACE
    PARAMETER(NFACE=4)
    INTEGER LOCLIST(NFACE,3),SILOC2ILOC(SNLOC)
    INTEGER SILOC,SJLOC,SKLOC,IGL,COL,INOD,JNOD,SGI,KLOC
    INTEGER SNSUBTLOC,ISUBNOD
    REAL RADI,RADJ,RADP,RAD,SAREA,VOL
    REAL NORMX,NORMY,NORMZ,VLNDOTQ,massVLNDOTQ,absVLNDOTQ,STABINEL,TDIVU
    REAL A11MAT,A12MAT,A21MAT,A22MAT
    REAL SNORMXN(SNGI),SNORMYN(SNGI),SNORMZN(SNGI)
    REAL NSUB(NSUBTLOC,NGI)
    REAL NSUBLX(NSUBTLOC,NGI),NSUBLY(NSUBTLOC,NGI),NSUBLZ(NSUBTLOC,NGI)
    REAL NSUBX(NSUBTLOC,NGI),NSUBY(NSUBTLOC,NGI),NSUBZ(NSUBTLOC,NGI)
    LOGICAL DMATLUM,SGSUPWIND,OPTWIND
    REAL HXGI,HYGI,HZGI
    REAL HOVERQ,SIG_THETA,HAT_SIG,EXP_SIG,ALPHA1
    REAL VLM2
    REAL ALPHA_OPT
    REAL SGS_UPWIND_FAC(NGI),SUF_SGS_UPWIND_FAC
    REAL sidWEI_NSUB,WEI_NSUB,sidVLN2,sidVLND2,sidTDIVU,RNN2
    LOGICAL SGSCROS_WIND
    REAL TUD(NGI),TVD(NGI),TWD(NGI),TDERI,RN
    INTEGER IDGNOD,NOD
    LOGICAL STREAM_UPWIND
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

    REAL, ALLOCATABLE, DIMENSION(:)::VECYT

    REAL, ALLOCATABLE, DIMENSION(:,:,:)::EV
    INTEGER, ALLOCATABLE, DIMENSION(:)::FINELE
    INTEGER, ALLOCATABLE, DIMENSION(:)::COLELE4

    REAL, ALLOCATABLE, DIMENSION(:,:)::SNSUB
    INTEGER, ALLOCATABLE, DIMENSION(:)::MLOCSILOC2ILOC

    ewrite(1, *) "just in SUBROUTINE HART3DSGS()"
    ewrite(2, *) "isphere=",isphere
    !         if(isphere.eq.0) pause

    DISOPT=DISOPT2
    IF(DISOPT2.EQ.47) DISOPT=147
    IF(DISOPT2.EQ.48) DISOPT=148

    !       ISPHERE=2
    !
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

       !         CALL TRIQUA(L1, L2, L3, L4, PGWEIG, .true.,NGI)
       ! Work out the shape functions and there derivatives...
       !         CALL SHATRI(L1, L2, L3, L4, PGWEIG, .true., 
       !     &            NCLOC,NGI,                          
       !     &            NC,NCLX,NCLY,NCLZ) 
    ELSE
       NC = N
       NCLX = NLX
       NCLY = NLY
       NCLZ = NLZ
    ENDIF

    !     
    do  I=1,NONODS! Was loop 1370
       ML(I)=0.
       VECX(I)=0.
    end do ! Was loop 1370
    BIGM(1:NCOLM) = 0.0
    !     
    RLUMP=1.
    IF(LUMP) RLUMP=1.
    RCONSI=1.-RLUMP
    !     
    !      ABLUMP=0.
    !      RSLUMP=0.
    ABLUMP=1.
    RSLUMP=1.
    IF(ABSLUM) ABLUMP=1.
    IF(SLUMP)  RSLUMP=1.
    ! ALPHA contains the inner element stabilisation that ensures 
    ! a non-singular system...
    !      ALPHA=0.0001.
    !      ALPHA=0.0
    !      ALPHA=1000000.
    !      ALPHA=0.001
    ALPHA=0.01
    !      ALPHA=1.0
    ! alpha=1.0 is good or DMATLUM=.TRUE. this way it works without 
    ! transport
    DMATLUM=.false.
    !       DMATLUM=.TRUE.
    !      ALPHA=1.0
    !      ALPHA=10.
    ! SGS upwinding only on SGS eqns.
    !       SGSUPWIND=.TRUE.
    SGSUPWIND=.FALSE.
    OPTWIND=.FALSE.
    !       OPTWIND=.true.
    !       SGSCROS_WIND=.TRUE.
    SGSCROS_WIND=.false.
    ! This is for streamline upwinding only on the global part
    ! and not a PG method.       
    !      STREAM_UPWIND=.false. 
    STREAM_UPWIND=(DISOPT.EQ.144).OR.(DISOPT.EQ.145)

    ! alpha2 control works for void but others do not.
    ALPHA2=0.0
    !      ALPHA2=0.01
    !      ALPHA2=0.0
    !      ALPHA2=10.0
    !      ALPHA2=1000000.0
    !      ALPHA2=10000000000.0

    !      ALPHA3=1.0
    !      ALPHA3=1000.0
    ALPHA3=0.0



    ! Calculate surface shape functions SN, SM etc...
    CALL TRIQUA(SL1, SL2, SL3, SL3, SWEIGH, .FALSE.,SNGI)
    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNLOC,SNGI, &
         SN,SNLX,SNLY,SNLX)

    CALL SHATRI(SL1, SL2, SL3, SL4, SWEIGH, .FALSE., &
         SNCLOC,SNGI, &
         SNC,SNCLX,SNCLY,SNCLX)

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
    !       print *,'NSUBTLOC,SNSUBTLOC,beta:',NSUBTLOC,SNSUBTLOC,beta
    !       stop 3931


    ALLOCATE(FINELE(TOTELE+1))
    ALLOCATE(COLELE4(TOTELE*5))
    ewrite(1,*) 'going into GETFINELE4'
    CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE4,NONODS)

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

    !     
    do  ELE=1,TOTELE! Was loop 340
       !     
       ! Get determinant and derivatives...
       CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,xnonod,NLOC,NCLOC,NGI, &
            N,NLX,NLY,NLZ, NC,NCLX,NCLY,NCLZ, WEIGHT, DETWEI,.TRUE.,.FALSE., &
            NX,NY,NZ, NCX,NCY,NCZ,&
            A11,A12,A13, A21,A22,A23, A31,A32,A33,&
            XD,YD,ZD,&
            ISPHERE) 

       do ILOC=1,NSUBTLOC! Was loop 
          do GI=1,NGI! Was loop 
             NSUBX(ILOC,GI)= A11(GI)*NSUBLX(ILOC,GI)+A12(GI)*NSUBLY(ILOC,GI)+A13(GI)*NSUBLZ(ILOC,GI)
             NSUBY(ILOC,GI)= A21(GI)*NSUBLX(ILOC,GI)+A22(GI)*NSUBLY(ILOC,GI)+A23(GI)*NSUBLZ(ILOC,GI)
             NSUBZ(ILOC,GI)= A31(GI)*NSUBLX(ILOC,GI)+A32(GI)*NSUBLY(ILOC,GI)+A33(GI)*NSUBLZ(ILOC,GI)
          END DO
       END DO

       ! ***START CALCULATING DERIVATIVES OF SUBU AT NODES**********
       IF(NSUBVLOC.NE.0) THEN
          do ILOC=1,NSUBVLOC! Was loop 
             do JLOC=1,NSUBVLOC! Was loop 
                VLM=0.0
                VLNXN=0.0
                VLNYN=0.0
                VLNZN=0.0
                do GI=1,NGI! Was loop 
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
          do ILOC=1,NSUBVLOC! Was loop 
             SOUSUBNU(ILOC)=0.0
             SOUSUBNV(ILOC)=0.0
             SOUSUBNW(ILOC)=0.0
             do JLOC=1,NSUBVLOC! Was loop 
                JNODSUB=(ELE-1)*NSUBVLOC+JLOC
                SOUSUBNU(ILOC)=SOUSUBNU(ILOC)-MATN1XN1(ILOC,JLOC)*NUSUB(JNODSUB)
                SOUSUBNV(ILOC)=SOUSUBNV(ILOC)-MATN1YN1(ILOC,JLOC)*NVSUB(JNODSUB)
                SOUSUBNW(ILOC)=SOUSUBNW(ILOC)-MATN1ZN1(ILOC,JLOC)*NWSUB(JNODSUB)
             END DO
          END DO
          do ILOC=1,NSUBVLOC! Was loop 
             DERSUBNUX(ILOC)=0.0
             DERSUBNVY(ILOC)=0.0
             DERSUBNWZ(ILOC)=0.0
             do JLOC=1,NSUBVLOC! Was loop 
                DERSUBNUX(ILOC)=DERSUBNUX(ILOC)+MASSINVM1M1(ILOC,JLOC)*SOUSUBNU(JLOC)
                DERSUBNVY(ILOC)=DERSUBNVY(ILOC)+MASSINVM1M1(ILOC,JLOC)*SOUSUBNV(JLOC)
                DERSUBNWZ(ILOC)=DERSUBNWZ(ILOC)+MASSINVM1M1(ILOC,JLOC)*SOUSUBNW(JLOC)
             END DO
          END DO
       ENDIF
       ! *****END CALCULATING DERIVATIVES OF SUBU AT NODES**********

       do  GI=1,NGI! Was loop 331
          !
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

          ABSDX(GI)=0.0
          ABSDY(GI)=0.0
          ABSDZ(GI)=0.0

          MUPTXXGI=0.0
          MUPTXYGI=0.0
          MUPTXZGI=0.0
          MUPTYYGI=0.0
          MUPTYZGI=0.0
          MUPTZZGI=0.0

          do ILOC=1,NSUBVLOC! Was loop 
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

          do  L=1,NLOC! Was loop 79
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
          end do ! Was loop 79
          !
          TUD(GI)=UD(GI)
          TVD(GI)=VD(GI)
          TWD(GI)=WD(GI)
          IF(SGSCROS_WIND) THEN
             TUD(GI)=0.
             TVD(GI)=0.
             TWD(GI)=0.
             TDERI=0.
             do ILOC=1,NLOC! Was loop 
                NOD=NDGLNO((ELE-1)*NLOC+ILOC)
                IDGNOD=(ELE-1)*NLOC+ILOC
                TDERI=TDERI+(UD(GI)*NX(ILOC,GI)+VD(GI)*NY(ILOC,GI)&
                     +WD(GI)*NZ(ILOC,GI))*(U(NOD)+USUB(IDGNOD))
                TUD(GI)=TUD(GI)+NX(ILOC,GI)*(U(NOD)+USUB(IDGNOD))
                TVD(GI)=TVD(GI)+NY(ILOC,GI)*(U(NOD)+USUB(IDGNOD))
                TWD(GI)=TWD(GI)+NZ(ILOC,GI)*(U(NOD)+USUB(IDGNOD))
             END DO
             RN=MAX(TUD(GI)**2+TVD(GI)**2+TWD(GI)**2,1.0E-9)
             TUD(GI)=TDERI*TUD(GI)/RN
             TVD(GI)=TDERI*TVD(GI)/RN
             TWD(GI)=TDERI*TWD(GI)/RN
             !
             HXGI=ABS(A11(GI)*TUD(GI)+A12(GI)*TVD(GI)+A13(GI)*TWD(GI))
             HYGI=ABS(A21(GI)*TUD(GI)+A22(GI)*TVD(GI)+A23(GI)*TWD(GI))
             HZGI=ABS(A31(GI)*TUD(GI)+A32(GI)*TVD(GI)+A33(GI)*TWD(GI))
             HOVERQ=2./MAX(HXGI,HYGI,HZGI,1.0E-9)
             SIG_THETA=THETA/DT
             HAT_SIG=SIG_THETA*HOVERQ
             SGS_UPWIND_FAC(GI)=0.5*HOVERQ*MIN(1.0,1./MAX(HAT_SIG,1.E-15))
          ENDIF
          !
          ALPHA1=1.0
          IF(SGSUPWIND) THEN
             HXGI=ABS(A11(GI)*UD(GI)+A12(GI)*VD(GI)+A13(GI)*WD(GI))
             HYGI=ABS(A21(GI)*UD(GI)+A22(GI)*VD(GI)+A23(GI)*WD(GI))
             HZGI=ABS(A31(GI)*UD(GI)+A32(GI)*VD(GI)+A33(GI)*WD(GI))
             !     HOVERQ=MIN(HXGI/RUD, HYGI/RVD,HZGI/RWD)
             !            HOVERQ=2./MAX(HXGI,HYGI,HZGI,1.0E-9)
             HOVERQ=0.5/MAX(HXGI,HYGI,HZGI,1.0E-9)
             !            print *,'ele,gi,hoverq:',ele,gi,hoverq
             !            HOVERQ=1./50.
             IF(OPTWIND) THEN
                SIG_THETA=THETA/DT
                HAT_SIG=SIG_THETA*HOVERQ
                EXP_SIG=EXPTOL(HAT_SIG)

                ALPHA_OPT=ABS( (6.0+4.*HAT_SIG+HAT_SIG**2) + (2.0*HAT_SIG-6.0)*EXP_SIG )&
                     /MAX(ABS( (12.0-6.0*HAT_SIG)*EXP_SIG -(12.0+6.0*HAT_SIG)),1.E-15)
                ALPHA_OPT=MIN(ALPHA_OPT, 1./MAX(HAT_SIG,1.E-15))
                SGS_UPWIND_FAC(GI)=HOVERQ*ALPHA_OPT
             ELSE
                !            SGS_UPWIND_FAC(GI)=4.*0.25/50.
                SIG_THETA=THETA/DT
                HAT_SIG=SIG_THETA*HOVERQ
                SGS_UPWIND_FAC(GI)=0.5*HOVERQ*MIN(1.0,1./MAX(HAT_SIG,1.E-15))
             ENDIF
          ELSE
             SGS_UPWIND_FAC(GI)=0.0
          ENDIF
          !
          ! calculate tensor for oceans:
          ! this has a QUADRATURE POINT varying element size and shape
          CALL SIZEGIELETENS(NX,NY,NZ,NLOC,NGI,GI,&
               L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ )

          IF((DISOPT.EQ.147).or.(DISOPT.EQ.150).or.(DISOPT.EQ.153)&
               .or.(DISOPT.EQ.148).or.(DISOPT.EQ.151).or.(DISOPT.EQ.154)) THEN
             !        IF(GETHIGH) THEN
             !        IF(.true.) THEN
             ! based on constant length scales...
             !         CALL TENFOROCEANS(L2GIXX,L2GIXY,L2GIXZ,L2GIYY,L2GIYZ,L2GIZZ,
             !     &                        XD(GI),YD(GI),ZD(GI)) 
             MXXDTW(GI)=L2GIXX
             MXYDTW(GI)=L2GIXY
             MXZDTW(GI)=L2GIXZ
             MYYDTW(GI)=L2GIYY
             MYZDTW(GI)=L2GIYZ
             MZZDTW(GI)=L2GIZZ
             CALL LESVIS(NONODS,TOTELE,NLOC,VONDGL, NGI, ELE,GI, &
                     NX,NY,NZ, &
                     NU,NV,NW, UG,VG,WG, &
                     MXXDTW,MXYDTW,MXZDTW,&
                     MYYDTW,MYZDTW,MZZDTW)


             MUGIXX=MXXDTW(GI)
             MUGIXY=MXYDTW(GI)
             MUGIXZ=MXZDTW(GI)
             MUGIYY=MYYDTW(GI)
             MUGIYZ=MYZDTW(GI)
             MUGIZZ=MZZDTW(GI)
          ELSE
             !         CALL TENFOROCEANSDIFF(MUGIXX,MUGIXY,MUGIXZ,MUGIYY,MUGIYZ,MUGIZZ,
             !     &                        XD(GI),YD(GI),ZD(GI),DISOPT) 
             MUGIXX=0.0
             MUGIXY=0.0
             MUGIXZ=0.0
             MUGIYY=0.0
             MUGIYZ=0.0
             MUGIZZ=0.0
          ENDIF

          ! 0.25 is used to scale down LES dissipation from ocean dissipation
          ! factor of 4 is used to reduce the length scale. 
          AMXXDTW(GI)=(0.125*MUGIXX+MUPTXXGI)*DETWEI(GI)
          AMXYDTW(GI)=(0.125*MUGIXY+MUPTXYGI)*DETWEI(GI)
          AMXZDTW(GI)=(0.125*MUGIXZ+MUPTXZGI)*DETWEI(GI)
          !
          AMYYDTW(GI)=(0.125*MUGIYY+MUPTYYGI)*DETWEI(GI)
          AMYZDTW(GI)=(0.125*MUGIYZ+MUPTYZGI)*DETWEI(GI)
          !
          AMZZDTW(GI)=(0.125*MUGIZZ+MUPTZZGI)*DETWEI(GI)
          IF(DISOPT.EQ.148) THEN
             ! Put LES only in the matrix AMAT.
             MXXDTW(GI)=(MUPTXXGI)*DETWEI(GI)
             MXYDTW(GI)=(MUPTXYGI)*DETWEI(GI)
             MXZDTW(GI)=(MUPTXZGI)*DETWEI(GI)
             !
             MYYDTW(GI)=(MUPTYYGI)*DETWEI(GI)
             MYZDTW(GI)=(MUPTYZGI)*DETWEI(GI)
             !
             MZZDTW(GI)=(MUPTZZGI)*DETWEI(GI)
          ELSE
             MXXDTW(GI)=(0.125*MUGIXX+MUPTXXGI)*DETWEI(GI)
             MXYDTW(GI)=(0.125*MUGIXY+MUPTXYGI)*DETWEI(GI)
             MXZDTW(GI)=(0.125*MUGIXZ+MUPTXZGI)*DETWEI(GI)
             !
             MYYDTW(GI)=(0.125*MUGIYY+MUPTYYGI)*DETWEI(GI)
             MYZDTW(GI)=(0.125*MUGIYZ+MUPTYZGI)*DETWEI(GI)
             !
             MZZDTW(GI)=(0.125*MUGIZZ+MUPTZZGI)*DETWEI(GI)
          ENDIF
          ! for LES...
       end do ! Was loop 331

       MASSLUMM1M1=0.0
       MASSLUMM2M2=0.0
       !     
       ! AMAT**********************
       do  ILOC=1,NLOC! Was loop 3501
          do  JLOC=1,NLOC! Was loop 3601
             !     
             VLM=0. 

             VLKTEN=0.0

             VLND  =0.0
             VLND2 =0.0

             VLMAB=0.0
             stVLKTEN=0.0
             !     
             do  GI=1,NGI! Was loop 4581

                VLN =-(NX(ILOC,GI)*UD(GI) &
                     +NY(ILOC,GI)*VD(GI) &
                     +NZ(ILOC,GI)*WD(GI) &
                     -N(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*N(JLOC,GI)*DETWEI(GI) 

                VLN2 =N(ILOC,GI)*(NX(JLOC,GI)*UD(GI) &
                     +NY(JLOC,GI)*VD(GI) &
                     +NZ(JLOC,GI)*WD(GI) )*DETWEI(GI)

                VLND=VLND+VLN
                VLND2=VLND2+VLN2

                !                 RNN=N(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                !                 RNN=NSUB(ILOC,GI)*( N(JLOC,GI) - (NX(JLOC,GI)*TUD(GI) 
                !     &                 +NY(JLOC,GI)*TVD(GI) 
                !     &                 +NZ(JLOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI)*ALPHA1*1.0 )*DETWEI(GI)
                RNN=( N(ILOC,GI) + (NX(ILOC,GI)*TUD(GI) &
                     +NY(ILOC,GI)*TVD(GI) &
                     +NZ(ILOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI)*ALPHA1*1.0 )*NSUB(JLOC,GI)*DETWEI(GI)


                VLM=VLM+RNN 

                XNDMUJLOC=(NX(JLOC,GI)*AMXXDTW(GI)&
                     +NY(JLOC,GI)*AMXYDTW(GI)+NZ(JLOC,GI)*AMXZDTW(GI))
                YNDMUJLOC=(NX(JLOC,GI)*AMXYDTW(GI)&
                     +NY(JLOC,GI)*AMYYDTW(GI)+NZ(JLOC,GI)*AMYZDTW(GI))
                ZNDMUJLOC=(NX(JLOC,GI)*AMXZDTW(GI)&
                     +NY(JLOC,GI)*AMYZDTW(GI)+NZ(JLOC,GI)*AMZZDTW(GI))

                XNDMUILOC=(NX(ILOC,GI)*AMXXDTW(GI)&
                     +NY(ILOC,GI)*AMXYDTW(GI)+NZ(ILOC,GI)*AMXZDTW(GI))
                YNDMUILOC=(NX(JLOC,GI)*AMXYDTW(GI)&
                     +NY(ILOC,GI)*AMYYDTW(GI)+NZ(ILOC,GI)*AMYZDTW(GI))
                ZNDMUILOC=(NX(ILOC,GI)*AMXZDTW(GI)&
                     +NY(ILOC,GI)*AMYZDTW(GI)+NZ(ILOC,GI)*AMZZDTW(GI))

                TEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC+NZ(ILOC,GI)*ZNDMUJLOC
                VLKTEN=VLKTEN+TEN

                VLMAB=VLMAB+ABGI(GI)*RNN   

                !                 UPWIN=UPWIN+(0.05/(2.0*sqrt(ud(gi)**2+vd(gi)**2+wd(gi)**2)))
                !     &            *(NX(ILOC,GI)*UD(GI)+NY(ILOC,GI)*VD(GI)+NZ(ILOC,GI)*WD(GI))
                !     &            *(NX(JLOC,GI)*UD(GI)+NY(JLOC,GI)*VD(GI)+NZ(JLOC,GI)*WD(GI))
                !     &            *DETWEI(GI)
                XNDMUJLOC=NX(JLOC,GI)*TUD(GI)*UD(GI)&
                      +NY(JLOC,GI)*TUD(GI)*VD(GI)&
                      +NZ(JLOC,GI)*TUD(GI)*WD(GI)
                YNDMUJLOC=NX(JLOC,GI)*TVD(GI)*UD(GI)&
                      +NY(JLOC,GI)*TVD(GI)*VD(GI)&
                      +NZ(JLOC,GI)*TVD(GI)*WD(GI)
                ZNDMUJLOC=NX(JLOC,GI)*TWD(GI)*UD(GI)&
                      +NY(JLOC,GI)*TWD(GI)*WD(GI)&
                      +NZ(JLOC,GI)*TWD(GI)*WD(GI)
                      
                TEN=NX(ILOC,GI)*XNDMUJLOC+NY(ILOC,GI)*YNDMUJLOC+NZ(ILOC,GI)*ZNDMUJLOC
                stVLKTEN=stVLKTEN+TEN*DETWEI(GI)*SGS_UPWIND_FAC(GI)*ALPHA1*1.0
                
             end do ! Was loop 4581

             ! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
             ! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
             ! conservative form BETA=1...
             ! non-conservative form BETA=0...
             AMAT(ILOC,JLOC)=BETA*VLND          +(1.-BETA)*VLND2 +VLMAB +VLKTEN&
                  +stVLKTEN

             MASSM1M1(ILOC,JLOC)=VLM
             MASSLUMM1M1(ILOC,ILOC)=MASSLUMM1M1(ILOC,ILOC)+VLM

          end do ! Was loop 3601
       end do ! Was loop 3501



       ! BMAT**********************
       do  ILOC=1,NLOC! Was loop 3502
          do  JLOC=1,NSUBTLOC! Was loop 3602
             !     
             VLM=0. 

             VLND  =0.0
             VLND3 =0.0

             VLMAB=0.0
             VLNXN=0.0
             VLNYN=0.0
             VLNZN=0.0

             stVLKTEN=0.0
             !     
             do  GI=1,NGI! Was loop 4582

                VLN =-(NX(ILOC,GI)*UD(GI) &
                     +NY(ILOC,GI)*VD(GI) &
                     +NZ(ILOC,GI)*WD(GI) &
                     -N(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI) 


                VLN3 =-(NX(ILOC,GI)*UD(GI) &
                     +NY(ILOC,GI)*VD(GI) &
                     +NZ(ILOC,GI)*WD(GI) &
                     +N(ILOC,GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI)

                VLND=VLND+VLN
                VLND3=VLND3+VLN3

                !                 RNN=N(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)
                !                 RNN=NSUB(ILOC,GI)*( N(JLOC,GI) - (NX(JLOC,GI)*TUD(GI) 
                !     &                 +NY(JLOC,GI)*TVD(GI) 
                !     &                 +NZ(JLOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI)*ALPHA1 )*DETWEI(GI)
                RNN=( N(ILOC,GI) + (NX(ILOC,GI)*TUD(GI) &
                     +NY(ILOC,GI)*TVD(GI) &
                     +NZ(ILOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI)*ALPHA1 )*NSUB(JLOC,GI)*DETWEI(GI)


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
             BMAT(ILOC,JLOC)=BETA*VLND          +(1.-BETA)*VLND3 +VLMAB&
                  +stVLKTEN

             MATN1XN2(ILOC,JLOC)=VLNXN
             MATN1YN2(ILOC,JLOC)=VLNYN
             MATN1ZN2(ILOC,JLOC)=VLNZN

             MASSM1M2(ILOC,JLOC)=VLM

          end do ! Was loop 3602
       end do ! Was loop 3502

       !           SGS_UPWIND_FAC=0.0

       ! CMAT**********************
       do  ILOC=1,NSUBTLOC! Was loop 3503
          do  JLOC=1,NLOC! Was loop 3603
             !     
             VLM=0. 
             VLM2=0. 

             VLND2 =0.0

             VLMAB=0.0

             TDIVU=0.0

             VLSTREAM_MATX=0.0
             VLSTREAM_MATY=0.0
             VLSTREAM_MATZ=0.0
             stVLKTEN=0.0
             !               print *,'SGS_UPWIND_FAC:',SGS_UPWIND_FAC
             !               print *,'ALPHA_OPT,HOVERQ:',ALPHA_OPT,HOVERQ
             !               stop 3893
             !     
             do  GI=1,NGI! Was loop 4583
                VLN2 =NSUB(ILOC,GI)*(NX(JLOC,GI)*UD(GI) &
                     +NY(JLOC,GI)*VD(GI) &
                     +NZ(JLOC,GI)*WD(GI) )*DETWEI(GI)

                VLND2=VLND2+VLN2
                TDIVU=TDIVU+NSUB(ILOC,GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI))*N(JLOC,GI)*DETWEI(GI) 

                RNN=NSUB(ILOC,GI)*( N(JLOC,GI) - (NX(JLOC,GI)*TUD(GI) &
                     +NY(JLOC,GI)*TVD(GI) &
                     +NZ(JLOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI) )*DETWEI(GI)
                !                 RNN=( N(ILOC,GI) + (NX(ILOC,GI)*TUD(GI) 
                !     &                 +NY(ILOC,GI)*TVD(GI) 
                !     &                 +NZ(ILOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI) )*NSUB(JLOC,GI)*DETWEI(GI)

                !                 RNN=NSUB(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)
                RNN2=NSUB(ILOC,GI)*N(JLOC,GI)*DETWEI(GI)

                VLM=VLM+RNN

                VLM2=VLM2+RNN2 

                VLMAB=VLMAB+ABGI(GI)*RNN  

                VLSTREAM_MATX=VLSTREAM_MATX+( NSUBX(ILOC,GI)*TUD(GI)*UD(GI)&
                     +NSUBY(ILOC,GI)*TVD(GI)*UD(GI)&
                     +NSUBZ(ILOC,GI)*TWD(GI)*UD(GI))&
                     *SGS_UPWIND_FAC(GI)*N(JLOC,GI)*DETWEI(GI)
                VLSTREAM_MATY=VLSTREAM_MATY+( NSUBX(ILOC,GI)*TUD(GI)*VD(GI)&
                     +NSUBY(ILOC,GI)*TVD(GI)*VD(GI)&
                     +NSUBZ(ILOC,GI)*TWD(GI)*VD(GI))&
                     *SGS_UPWIND_FAC(GI)*N(JLOC,GI)*DETWEI(GI)
                VLSTREAM_MATZ=VLSTREAM_MATZ+( NSUBX(ILOC,GI)*TUD(GI)*WD(GI)&
                     +NSUBY(ILOC,GI)*TVD(GI)*WD(GI)&
                     +NSUBZ(ILOC,GI)*TWD(GI)*WD(GI))&
                     *SGS_UPWIND_FAC(GI)*N(JLOC,GI)*DETWEI(GI)
             end do ! Was loop 4583

             ! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
             ! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
             ! conservative form BETA=1...
             ! non-conservative form BETA=0...
             CMAT(ILOC,JLOC)=BETA*(VLND2 +TDIVU)+(1.-BETA)*VLND2 +VLMAB&
                  +stVLKTEN
             !        CMAT(ILOC,JLOC)=BETA*(VLND2 )+(1.-BETA)*VLND2 +VLMAB

             MASSM2M1(ILOC,JLOC)=VLM

             STREAM_MATX(ILOC,JLOC)=VLSTREAM_MATX
             STREAM_MATY(ILOC,JLOC)=VLSTREAM_MATY
             STREAM_MATZ(ILOC,JLOC)=VLSTREAM_MATZ

          end do ! Was loop 3603
       end do ! Was loop 3503


       ! DMAT**********************
       do  ILOC=1,NSUBTLOC! Was loop 3504
          do  JLOC=1,NSUBTLOC! Was loop 3604
             !     
             VLM=0.
             VLM2=0. 

             VLND  =0.0
             VLND3 =0.0

             sidVLND2 =0.0
             sidTDIVU=0.0

             VLMAB=0.0
             VLNXN=0.0
             VLNYN=0.0
             VLNZN=0.0

             VLNNX=0.0  
             VLNNY=0.0
             VLNNZ=0.0
             STABINEL=0.0
             !     
             do  GI=1,NGI! Was loop 4584
                sidWEI_NSUB=(NSUBX(ILOC,GI)*TUD(GI) &
                     +NSUBY(ILOC,GI)*TVD(GI) &
                     +NSUBZ(ILOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI)
                WEI_NSUB=NSUB(ILOC,GI)+(NSUBX(ILOC,GI)*TUD(GI) &
                     +NSUBY(ILOC,GI)*TVD(GI) &
                     +NSUBZ(ILOC,GI)*TWD(GI))*SGS_UPWIND_FAC(GI)

                sidVLN2 =sidWEI_NSUB*(NSUBX(JLOC,GI)*TUD(GI) &
                     +NSUBY(JLOC,GI)*TVD(GI) &
                     +NSUBZ(JLOC,GI)*TWD(GI) )*DETWEI(GI)

                sidVLND2=sidVLND2+sidVLN2
                sidTDIVU=sidTDIVU+sidWEI_NSUB*(UDDX(GI)+VDDY(GI)+WDDZ(GI))*N(JLOC,GI)*DETWEI(GI) 

                RNN=WEI_NSUB*NSUB(JLOC,GI)*DETWEI(GI)

                RNN2=NSUB(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)

                VLM=VLM+RNN 
                VLM2=VLM2+RNN2 

                VLN =-(NSUBX(ILOC,GI)*UD(GI) &
                     +NSUBY(ILOC,GI)*VD(GI) &
                     +NSUBZ(ILOC,GI)*WD(GI) &
                     -NSUB(ILOC,GI)*(UGDX(GI)+VGDY(GI)+WGDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI) 


                VLN3 =-(NSUBX(ILOC,GI)*UD(GI) &
                     +NSUBY(ILOC,GI)*VD(GI) &
                     +NSUBZ(ILOC,GI)*WD(GI) &
                     +NSUB(ILOC,GI)*(UDDX(GI)+VDDY(GI)+WDDZ(GI)) )*NSUB(JLOC,GI)*DETWEI(GI)

                VLND=VLND+VLN
                VLND3=VLND3+VLN3


                VLMAB=VLMAB+ABGI(GI)*RNN   
                VLNXN=VLNXN+NSUBX(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)   
                VLNYN=VLNYN+NSUBY(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI) 
                VLNZN=VLNZN+NSUBZ(ILOC,GI)*NSUB(JLOC,GI)*DETWEI(GI)  

                VLNNX=VLNNX+NSUB(ILOC,GI)*NSUBX(JLOC,GI)*DETWEI(GI)   
                VLNNY=VLNNY+NSUB(ILOC,GI)*NSUBY(JLOC,GI)*DETWEI(GI) 
                VLNNZ=VLNNZ+NSUB(ILOC,GI)*NSUBZ(JLOC,GI)*DETWEI(GI)  

                STABINEL=STABINEL+(NSUBX(ILOC,GI)*(NSUBX(JLOC,GI)*L2GIXX&
                     +NSUBY(JLOC,GI)*L2GIXY+NSUBZ(JLOC,GI)*L2GIXZ)&
                     +NSUBY(ILOC,GI)*(NSUBX(JLOC,GI)*L2GIXY&
                     +NSUBY(JLOC,GI)*L2GIYY+NSUBZ(JLOC,GI)*L2GIYZ)&
                     +NSUBZ(ILOC,GI)*(NSUBX(JLOC,GI)*L2GIXZ&
                     +NSUBY(JLOC,GI)*L2GIYZ+NSUBZ(JLOC,GI)*L2GIZZ)&
                     )*DETWEI(GI)
                !
             end do ! Was loop 4584

             ! TENSOR form of viscosity (USED BECAUSE OF COMPLEXITY OF STRESS FORM...      
             ! assume a steady state to calculate AMAT,BMAT,CMAT,DMAT initially.
             ! conservative form BETA=1...
             ! non-conservative form BETA=0...
             DMAT(ILOC,JLOC)=BETA*VLND          +(1.-BETA)*VLND3 +VLMAB &
                  +(ALPHA/DT)*STABINEL
             !     &                 +BETA*(sidVLND2 +sidTDIVU)+(1.-BETA)*sidVLND2

             MATN2XN2(ILOC,JLOC)=VLNXN
             MATN2YN2(ILOC,JLOC)=VLNYN
             MATN2ZN2(ILOC,JLOC)=VLNZN

             MATN2N2X(ILOC,JLOC)=VLNNX
             MATN2N2Y(ILOC,JLOC)=VLNNY
             MATN2N2Z(ILOC,JLOC)=VLNNZ

             MASSM2M2(ILOC,JLOC)=VLM
             REALMASSM2M2(ILOC,JLOC)=VLM2
             MASSLUMM2M2(ILOC,ILOC)=MASSLUMM2M2(ILOC,ILOC)+VLM2

          end do ! Was loop 3604
       end do ! Was loop 3504

       ! calculate ELEMENT average quantities... 
       MUGIXX=0.0
       MUGIXY=0.0
       MUGIXZ=0.0
       MUGIYY=0.0
       MUGIYZ=0.0
       MUGIZZ=0.0
       VOL=0.0
       SUF_SGS_UPWIND_FAC=0.0
       do GI=1,NGI! Was loop 
          !           MUGIXX=MUGIXX+MXXDTW(GI)
          !           MUGIXY=MUGIXY+MXYDTW(GI)
          !           MUGIXZ=MUGIXZ+MXZDTW(GI)
          !           MUGIYY=MUGIYY+MYYDTW(GI)
          !           MUGIYZ=MUGIYZ+MYZDTW(GI)
          !           MUGIZZ=MUGIZZ+MZZDTW(GI)

          MUGIXX=MUGIXX+MXXDTW(GI) +TUD(GI)*UD(GI)*SGS_UPWIND_FAC(GI)*DETWEI(GI)
          MUGIXY=MUGIXY+MXYDTW(GI) +TUD(GI)*VD(GI)*SGS_UPWIND_FAC(GI)*DETWEI(GI)
          MUGIXZ=MUGIXZ+MXZDTW(GI) +TUD(GI)*WD(GI)*SGS_UPWIND_FAC(GI)*DETWEI(GI)
          MUGIYY=MUGIYY+MYYDTW(GI) +TVD(GI)*VD(GI)*SGS_UPWIND_FAC(GI)*DETWEI(GI)
          MUGIYZ=MUGIYZ+MYZDTW(GI) +TVD(GI)*WD(GI)*SGS_UPWIND_FAC(GI)*DETWEI(GI)
          MUGIZZ=MUGIZZ+MZZDTW(GI) +TWD(GI)*WD(GI)*SGS_UPWIND_FAC(GI)*DETWEI(GI)

          SUF_SGS_UPWIND_FAC=SUF_SGS_UPWIND_FAC+SGS_UPWIND_FAC(GI)*DETWEI(GI)
          VOL=VOL+DETWEI(GI)
       END DO
       MUGIXX=MUGIXX/VOL
       MUGIXY=MUGIXY/VOL
       MUGIXZ=MUGIXZ/VOL
       MUGIYY=MUGIYY/VOL
       MUGIYZ=MUGIYZ/VOL
       MUGIZZ=MUGIZZ/VOL
       SUF_SGS_UPWIND_FAC=SUF_SGS_UPWIND_FAC/VOL

       ! Put surface integrals into DMAT *********************
       ! Surface element of domain...
       do  IFACE=1,NFACE! Was loop 3344

          do SILOC=1,SNLOC! Was loop 
             SILOC2ILOC(SILOC)=LOCLIST(IFACE,SILOC)
          END DO

          do SILOC=1,SNLOC! Was loop 
             ILOC=SILOC2ILOC(SILOC)
             INOD=XONDGL((ELE-1)*NLOC+ILOC)
             do SJLOC=SILOC,SNLOC,1! Was loop 
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
               X,Y,Z,XNONOD,NORMX,NORMY,NORMZ)
          ! Recalculate the normal...
          CALL DGSDETNXLOC2(SNCLOC,SNGI, &
               XSL,YSL,ZSL, &
               SNC,SNCLX,SNCLY, SWEIGH, SDETWE,SAREA,.TRUE.,.FALSE., &
               SNORMXN,SNORMYN,SNORMZN, &
               NORMX,NORMY,NORMZ)

          SUD=0.0
          SVD=0.0
          SWD=0.0
          SUDGLOB=0.0
          SVDGLOB=0.0
          SWDGLOB=0.0
          SABS=0.0
          do SILOC=1,SNLOC! Was loop 
             ILOC=SILOC2ILOC(SILOC)
             IGL =NDGLNO((ELE-1)*NLOC+ILOC)
             ISUBNOD=(ELE-1)*NSUBVLOC+ILOC
             do  SGI=1,SNGI! Was loop 3313
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
                SABS(SGI)=SABS(SGI)+SN(SILOC,SGI)*ABSORB(IGL)
             END DO !Was end of loop 3313
          END DO
          !           print *,'sud:',sud
          !           print *,'svd:',svd
          !           print *,'swd:',swd
          !           print *,'SNORMXN=',SNORMXN
          !           print *,'SNORMYN=',SNORMyN
          !           print *,'SNORMZN=',SNORMzN
          !           print *,'NORMX=',NORMX
          !           print *,'NORMY=',NORMy
          !           print *,'NORMZ=',NORMz
          !           print *,'SWEIGH=',SWEIGH
          !           print *,'SDETWE=',SDETWE
          !           print *,'SAREA=',SAREA
          !           stop 3832
          do SGI=1,SNGI! Was loop 
             NDOTQ(SGI)=SUD(SGI)*SNORMXN(SGI)+SVD(SGI)*SNORMYN(SGI) &
                  +SWD(SGI)*SNORMZN(SGI)
             NDOTQGLOB(SGI)=SUDGLOB(SGI)*SNORMXN(SGI)+SVDGLOB(SGI)*SNORMYN(SGI) &
                  +SWDGLOB(SGI)*SNORMZN(SGI)
          END DO
          !           if(iface.eq.2) then
          !           print *,'iface,ndotq:',iface,ndotq
          !                print *,'SNORMzN(1):',SNORMxN(1),SNORMyN(1),SNORMzN(1)
          !                print *,'SUD(1):',SUD(1),SvD(1),SwD(1)
          !           endif
          ! ************************
          ! Perform surface integration...
          IF(NSUBTLOC.EQ.10) THEN
             !           print *,'before MLOCSILOC2ILOC calc'
             !           stop 663
             ! Calculate MLOCSILOC2ILOC
             do SILOC=1,SNLOC! Was loop 
                ILOC=SILOC2ILOC(SILOC)
                do SJLOC=1,SNLOC! Was loop 
                   JLOC=SILOC2ILOC(SJLOC)
                   SKLOC=SILINK2(SILOC,SJLOC)
                   KLOC=ILINK2(ILOC,JLOC)
                   MLOCSILOC2ILOC(SKLOC)=KLOC
                END DO
             END DO
             !           print *,'after MLOCSILOC2ILOC calc',MLOCSILOC2ILOC
             !           stop 684

             do  SILOC=1,SNSUBTLOC! Was loop 3508
                ILOC=MLOCSILOC2ILOC(SILOC)
                do  SJLOC=1,SNSUBTLOC! Was loop 3507
                   JLOC=MLOCSILOC2ILOC(SJLOC)
                   ! Have a surface integral on outlet boundary...                
                   VLNDOTQ=0.0        
                   absVLNDOTQ=0.0         
                   massVLNDOTQ=0.0
                   do  SGI=1,SNGI! Was loop 4585
                      RNN=SDETWE(SGI)*SNSUB(SILOC,SGI)*SNSUB(SJLOC,SGI)
                      VLNDOTQ=VLNDOTQ+MAX(NDOTQ(SGI),0.0)*RNN
                      absVLNDOTQ=absVLNDOTQ+MIN(NDOTQ(SGI),0.0)*RNN*SABS(SGI)
                      massVLNDOTQ=massVLNDOTQ+MIN(NDOTQ(SGI),0.0)*RNN
                      !                  VLNDOTQ=VLNDOTQ+MAX(NDOTQGLOB(SGI),0.0)*RNN
                      !                  absVLNDOTQ=absVLNDOTQ+MIN(NDOTQGLOB(SGI),0.0)*RNN*SABS(SGI)
                      !                  massVLNDOTQ=massVLNDOTQ+MIN(NDOTQGLOB(SGI),0.0)*RNN
                   END DO ! Was end of loop 4585
                   DMAT(ILOC,JLOC)=DMAT(ILOC,JLOC)+VLNDOTQ-absVLNDOTQ*SUF_SGS_UPWIND_FAC
                   MASSM2M2(ILOC,JLOC)=MASSM2M2(ILOC,JLOC)-massVLNDOTQ*SUF_SGS_UPWIND_FAC
                END DO ! Was end of loop 3507
             END DO ! Was end of loop 3508
          ELSE
             do  SILOC=1,SNSUBTLOC! Was loop 35011
                ILOC=SILOC2ILOC(SILOC)
                do  SJLOC=1,SNSUBTLOC! Was loop 35071
                   JLOC=SILOC2ILOC(SJLOC)
                   ! Have a surface integral on outlet boundary...                
                   VLNDOTQ=0.0              
                   absVLNDOTQ=0.0         
                   massVLNDOTQ=0.0
                   do  SGI=1,SNGI! Was loop 45851
                      RNN=SDETWE(SGI)*SNSUB(SILOC,SGI)*SNSUB(SJLOC,SGI)
                      VLNDOTQ=VLNDOTQ+MAX(NDOTQ(SGI),0.0)*RNN
                      absVLNDOTQ=absVLNDOTQ+MIN(NDOTQ(SGI),0.0)*RNN*SABS(SGI)
                      massVLNDOTQ=massVLNDOTQ+MIN(NDOTQ(SGI),0.0)*RNN
                      !                  VLNDOTQ=VLNDOTQ+MAX(NDOTQGLOB(SGI),0.0)*RNN
                      !                  absVLNDOTQ=absVLNDOTQ+MIN(NDOTQGLOB(SGI),0.0)*RNN*SABS(SGI)
                      !                  massVLNDOTQ=massVLNDOTQ+MIN(NDOTQGLOB(SGI),0.0)*RNN
                   END DO ! Was end of loop 45851
                   DMAT(ILOC,JLOC)=DMAT(ILOC,JLOC)+VLNDOTQ-absVLNDOTQ*SUF_SGS_UPWIND_FAC
                   MASSM2M2(ILOC,JLOC)=MASSM2M2(ILOC,JLOC)-massVLNDOTQ*SUF_SGS_UPWIND_FAC
                END DO ! Was end of loop 35071
             END DO ! Was end of loop 35071
          ENDIF
          
      end do ! Was loop 3344
      !           print *,'just after 3344'

      ! Calculate the diffusion terms*****
      ! calculate MASS2INV...
      CALL MATDMATINV(REALMASSM2M2,MASSINVM2M2,NSUBTLOC)
      CALL ABMATRIXMULTI(KMAT22X, MASSINVM2M2, MATN2XN2,NSUBTLOC)
      KMAT22X=-KMAT22X
      CALL ABMATRIXMULTI(KMAT22Y, MASSINVM2M2, MATN2YN2,NSUBTLOC)
      KMAT22Y=-KMAT22Y
      CALL ABMATRIXMULTI(KMAT22Z, MASSINVM2M2, MATN2ZN2,NSUBTLOC)
      KMAT22Z=-KMAT22Z


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
           DIFMATX,NSUBTLOC,NSUBTLOC)
      CALL ABMATRIXMUL(DIFMATY2B, MATN1YN2,NLOC,NSUBTLOC, &
           DIFMATY,NSUBTLOC,NSUBTLOC)
      CALL ABMATRIXMUL(DIFMATZ2B, MATN1ZN2,NLOC,NSUBTLOC, &
           DIFMATZ,NSUBTLOC,NSUBTLOC)
      DIFMATB=DIFMATX2B+DIFMATY2B+DIFMATZ2B

      BMAT=BMAT+DIFMATB*1.
      !         BMAT=BMAT+DIFMAT*1.
      DMAT=DMAT-DIFMAT*1.
      ! The diffusion terms associated with advection stabilisation:
      CALL ABMATRIXMUL(DIFMATXC, STREAM_MATX,NSUBTLOC,NSUBTLOC, &
           KMAT22X,    NSUBTLOC,NSUBTLOC)
      CALL ABMATRIXMUL(DIFMATYC, STREAM_MATY,NSUBTLOC,NSUBTLOC, &
           KMAT22Y,    NSUBTLOC,NSUBTLOC)
      CALL ABMATRIXMUL(DIFMATZC, STREAM_MATZ,NSUBTLOC,NSUBTLOC, &
           KMAT22Z,    NSUBTLOC,NSUBTLOC)
      ! C matrix...         
      do ILOC=1,NSUBTLOC! Was loop 
         do JLOC=1,NLOC ! Was loop 
            CMAT(ILOC,JLOC)=CMAT(ILOC,JLOC)+DIFMATB(JLOC,ILOC)*1.
            !           CMAT(ILOC,JLOC)=CMAT(ILOC,JLOC)+DIFMAT(JLOC,ILOC)*1.
            !     &   -0.5*(DIFMATXC(ILOC,JLOC)+DIFMATYC(ILOC,JLOC)+DIFMATZC(ILOC,JLOC))
            !     &   -0.5*(DIFMATXC(JLOC,ILOC)+DIFMATYC(JLOC,ILOC)+DIFMATZC(JLOC,ILOC))
            !     &   -0.0*(DIFMATXC(ILOC,JLOC)+DIFMATYC(ILOC,JLOC)+DIFMATZC(ILOC,JLOC))
            !     &   -1.0*(DIFMATXC(JLOC,ILOC)+DIFMATYC(JLOC,ILOC)+DIFMATZC(JLOC,ILOC))
         END DO
      END DO


      ! Calculate the diffusion terms***** 
      !
      !          print *,'after diffusion'
      ! 
      ! Calculate the discretised sources...
      SOUSUBU=0.0
      SOUSUBV=0.0
      SOUSUBW=0.0
      ! AMAT:
      do ILOC=1,NLOC   ! Was loop 
         GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
         do JLOC=1,NLOC! Was loop 
            GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
            A11MAT=MASSM1M1(ILOC,JLOC)/DT -(1.-THETA)*AMAT(ILOC,JLOC)
            VECX(GLOBI)=VECX(GLOBI) &
                 +A11MAT*U(GLOBJ) +MASSM1M1(ILOC,JLOC)*SOURCX(GLOBJ)
            ! Lump the absorption to put in to pressure. 
            ML(GLOBI)=ML(GLOBI) + MASSM1M1(ILOC,JLOC)
         END DO
         ! BMAT:
         do JLOC=1,NSUBTLOC! Was loop 
            A12MAT=MASSM1M2(ILOC,JLOC)/DT -(1.-THETA)*BMAT(ILOC,JLOC)
            VECX(GLOBI)=VECX(GLOBI) &
                 +A12MAT*USUB((ELE-1)*NSUBTLOC+JLOC)
         END DO
      END DO

      ! CMAT:
      do ILOC=1,NSUBTLOC   ! Was loop 
         do JLOC=1,NLOC! Was loop 
            GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
            A21MAT=MASSM2M1(ILOC,JLOC)/DT -(1.-THETA)*CMAT(ILOC,JLOC)
            ! Put the SGS source in: 
            SOUSUBU(ILOC)=SOUSUBU(ILOC)&
                 +A21MAT*U(GLOBJ) +MASSM2M1(ILOC,JLOC)*SOURCX(GLOBJ)
         END DO

         ! DMAT:
         do JLOC=1,NSUBTLOC! Was loop 
            IF(DMATLUM) THEN
               A22MAT=MASSlumM2M2(ILOC,JLOC)/DT -(1.-THETA)*DMAT(ILOC,JLOC)
            ELSE
               A22MAT=MASSM2M2(ILOC,JLOC)/DT -(1.-THETA)*DMAT(ILOC,JLOC)
            ENDIF
            ! Put the SGS source in: 
            SOUSUBU(ILOC)=SOUSUBU(ILOC)&
                 +A22MAT*USUB((ELE-1)*NSUBTLOC+JLOC)
         END DO
      END DO
      !
      !
      ! Amend AMAT,BMAT,CMAT,DMAT to take into account time dep.
      AMAT=(1./DT)*MASSM1M1+THETA*AMAT
      BMAT=(1./DT)*MASSM1M2+THETA*BMAT
      CMAT=(1./DT)*MASSM2M1+THETA*CMAT
      IF(DMATLUM) THEN
         DMAT=(1./DT)*MASSLUMM2M2+THETA*DMAT
      ELSE
         DMAT=(1./DT)*MASSM2M2+THETA*DMAT
      ENDIF

      !         print *,'got matrices'

      ! calculate DMATINV...
      CALL MATDMATINV(DMAT,DMATINV,NSUBTLOC)
      ! calculate D^{-1} C
      CALL ABMATRIXMUL(DMATINVC,DMATINV,NSUBTLOC,NSUBTLOC,&
           CMAT,NSUBTLOC,NLOC)
      ! BDINVC=B D^{-1} C
      CALL ABMATRIXMUL(BDINVC,BMAT,NLOC,NSUBTLOC,&
           DMATINVC,NSUBTLOC,NLOC)
      ! BDINVMAT=B D^{-1} 
      CALL ABMATRIXMUL(BDINVMAT,BMAT,NLOC,NSUBTLOC,&
           DMATINV,NSUBTLOC,NSUBTLOC)
      !          print *,'got mats'

      ! store critical matrices to calculate new SGS soln...
      do ILOC=1,NSUBTLOC! Was loop 
         DINVSOUSUBUSTORE(ELE,ILOC)=0.0
         do JLOC=1,NSUBTLOC! Was loop 
            DINVSOUSUBUSTORE(ELE,ILOC)=DINVSOUSUBUSTORE(ELE,ILOC)&
                 +DMATINV(ILOC,JLOC)*SOUSUBU(JLOC)
         END DO
      END DO
      DMATINVCSTORE(ELE,1:NSUBTLOC,1:NLOC)=DMATINVC(1:NSUBTLOC,1:NLOC)

      !         print *,'got to here'
      !
      ! amend VECX to take into account -B D^{-1}*SOUSUBU:
      do ILOC=1,NLOC  ! Was loop 
         GLOBI=NDGLNO((ELE-1)*NLOC+ILOC) 
         do JLOC=1,NSUBTLOC! Was loop 
            VECX(GLOBI)=VECX(GLOBI)-BDINVMAT(ILOC,JLOC)*SOUSUBU(JLOC)
         END DO
      END DO

      do ILOC=1,NLOC! Was loop 
         GLOBI=NDGLNO((ELE-1)*NLOC+ILOC)
         do JLOC=1,NLOC! Was loop 
            GLOBJ=NDGLNO((ELE-1)*NLOC+JLOC)
            ! Put into matrix...
            ! Find count. ***************************************
            CALL POSINMAT(COUNT,GLOBI,GLOBJ,NONODS,FINDRM,COLM,NCOLM)
            BIGM(COUNT)=BIGM(COUNT)+AMAT(ILOC,JLOC)-BDINVC(ILOC,JLOC)
         END DO
      END DO

    end do ! Was loop 340
    !
    ! Now convert the matrix eqn into rate of change form... 
    do GLOBI=1,NONODS! Was loop 
        do COUNT=FINDRM(GLOBI),FINDRM(GLOBI+1)-1! Was loop 
             COL=COLM(COUNT)
             VECX(GLOBI)=VECX(GLOBI)-BIGM(COUNT)*U(COL)
        END DO
    END DO
    BIGM=BIGM*DT

  END SUBROUTINE HART3DSGS

end module advection_diffusion_cg_3d
