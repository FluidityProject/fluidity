#include "fdebug.h"

module reduced_model

  use allSorts
  use coordinates
  use fldebug
  use shape_transformations
  use solvers
  use tr2d_module
  use coriolis_module

  implicit none

  private

  public :: in_output_pg, reduced_model_diffcoef1, &
    & reduced_model_diffcoef1_pg1, vtkoutput1p

contains

SUBROUTINE REDUCED_MODEL_DIFFCOEF1(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,X,Y,Z,                       &
     LEFTSVD,SMEAN,COEF,M,MLX,MLY,MLZ,                                    &
     N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
     NTIME1,DT1,D3,DCYL,GEOBAL,SCFACTH0,IFWIND)

  INTEGER, INTENT(IN)  ::NSVD_TOTAL,NSVD,NVAR,NSVD_U,NSVD_PHI
  INTEGER, INTENT(IN)  ::TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,MXNODS
  INTEGER, INTENT(IN) ::NTIME1
  REAL, INTENT(INOUT),DIMENSION(NONODS) ::X,Y,Z
  REAL, INTENT(INOUT),DIMENSION(NONODS) :: MUPTXX,MUPTXY,MUPTXZ
  REAL, INTENT(INOUT),DIMENSION(NONODS) :: MUPTYY,MUPTYZ,MUPTZZ
  double precision, INTENT(INOUT),DIMENSION(NVAR*MXNODS,NSVD) :: LEFTSVD
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS) ::SMEAN
  REAL, INTENT(INOUT),DIMENSION(0:ntime1,NSVD_TOTAL) :: COEF
  REAL, INTENT(IN)    ::DT1
  LOGICAL, INTENT(IN) ::D3,DCYL
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: N,NLX,NLY,NLZ
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: M,MLX,MLY,MLZ
  REAL, INTENT(IN),DIMENSION(NGI)  ::WEIGHT
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC)::PNDGLN
  ! coriolos force
  INTEGER, INTENT(IN)  :: GEOBAL
  REAL, INTENT(IN)  :: SCFACTH0
  ! Wind stress
  !------------
  LOGICAL IFWIND

  !       Local variables
  REAL    ::THETA,ONEMINTH
  INTEGER ::ISVD,JSVD,II,JJ,KK,NSTEP1,ELE,GI,ILOC,JLOC,I,J,K
  INTEGER ::IGL,JGL,IGLP,COUNT,COUNT1,JGLP
  REAL    ::RNN,RN1,RN2,RN3,FU,FV,U_x,V_y,W_z,UU_x,VU_y,WU_z,                        &
            UV_x,VV_y,WV_z,UW_x,VW_y,WW_z,P_x,P_y,P_z,                               &
            UD_Ux,VD_Uy,WD_Uz,UD_Vx,VD_Vy,WD_Vz,UD_Wx,VD_Wy,WD_Wz,                   &
            UDmean_Ux,VDmean_Uy,WDmean_Uz,UDmean_Vx,VDmean_Vy,WDmean_Vz,             &
            UDmean_Wx,VDmean_Wy,WDmean_Wz,                                           &
            Umeanx_UD,Umeany_VD,Umeanz_WD,Vmeanx_UD,Vmeany_VD,Vmeanz_WD,             &
            Wmeanx_UD,Wmeany_VD,Wmeanz_WD,                                           &
            U_xx,U_yy,U_zz,V_xx,V_yy,V_zz,W_xx,W_yy,W_zz
  REAL    ::RN1_mean,RN2_mean,RN3_mean,FU_mean,FV_mean,U_x_mean,V_y_mean,W_z_mean,         &
            UU_x_mean,VU_y_mean,WU_z_mean,UV_x_mean,VV_y_mean,WV_z_mean,                   &
            UW_x_mean,VW_y_mean,WW_z_mean,P_x_mean,P_y_mean,P_z_mean,                      &
            U_xx_mean,U_yy_mean,U_zz_mean,V_xx_mean,V_yy_mean,V_zz_mean,W_xx_mean,W_yy_mean,W_zz_mean
  INTEGER ::NCOLM,NCOLM1
  INTEGER ::IBL11,IBL12,IBL13,IBL14
  INTEGER ::IBL21,IBL22,IBL23,IBL24
  INTEGER ::IBL31,IBL32,IBL33,IBL34
  INTEGER ::IBL41,IBL42,IBL43,IBL44
  real    ::VOL,ftest,VOLUME
  REAL    ::DT
  INTEGER ::NSTEP,NTIME,NT

  !       Working arrays/local variables            
  !       ********* DEFINE OPTIONS FOR GMRES **************
  REAL     ::ERROR1,ERROR2,RELAX
  INTEGER  ::NOITS1,NOITS2,LEFTP1,LEFTP2,IMM
  INTEGER  ::MINTE1,MINTE2,NR2,IP,MISOUT
  !       ** options for MULMAT 
  LOGICAL  ::SYM,IGUESS
  INTEGER  ::MULPA,UZAWA,CGSOLQ,GMRESQ,CGSQ
  !???? wrong
  REAL     ::DM1(NSVD),DM2(NSVD),R(7*NSVD)
  real     ::rub(2)
  integer  ::iidum(2),kits
  INTEGER  ::MATSTR,NR,NNV
  PARAMETER(MATSTR=0)
  REAL     ::RRR
  INTEGER  ::PARA,halo_tag, halo_tag_p
  PARAMETER(PARA=0,halo_tag=1, halo_tag_p=2)

! all matraces for the calculations:
  !???? wrong
  INTEGER  FINDRM(NSVD+1),COLM(NSVD*NSVD), CENTRM(NSVD)

  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM,VEC,VEC1
  REAL,DIMENSION(:,:),ALLOCATABLE   :: ReMat,ReMat1
  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM0,BMat,BMat1,BMat0,BMat01,VEC0
  REAL,DIMENSION(:),ALLOCATABLE     :: DETWEI,XD,YD,ZD
  REAL,DIMENSION(:),ALLOCATABLE     :: UD_mean,VD_mean,WD_mean
  REAL,DIMENSION(:,:),ALLOCATABLE   :: UD,VD,WD
  REAL,DIMENSION(:),ALLOCATABLE     :: Umeanx,Umeany,Umeanz,Vmeanx,Vmeany,Vmeanz,Wmeanx,Wmeany,Wmeanz
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_u,WEI_v,WEI_w,WEI_p
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_ux,WEI_uy,WEI_uz,WEI_vx,WEI_vy,WEI_vz,WEI_wx,WEI_wy,WEI_wz
  REAL,DIMENSION(:,:),ALLOCATABLE   :: NX,NY,NZ,MX,MY,MZ
  REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF,NEWCOEF1,NEWCOEF2
  REAL,DIMENSION(:),ALLOCATABLE     :: CORIOLIS

  REAL,DIMENSION(:,:),ALLOCATABLE   :: BIGM1
  REAL,DIMENSION(:),ALLOCATABLE     :: UDMat
!  REAL,DIMENSION(:),ALLOCATABLE     :: FINDRM,COLM,CENTRM
  REAL,DIMENSION(:),ALLOCATABLE     :: ML

  REAL,DIMENSION(:),ALLOCATABLE     :: SOURCX,SOURCY,SOURCZ

  NT=1
  DT=DT1/NT
  NTIME=NT*NTIME1

  ALLOCATE(BIGM(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(VEC(NSVD_TOTAL))
  ALLOCATE(ReMat( NSVD_TOTAL-NSVD,NSVD_TOTAL-NSVD))
  ALLOCATE(BIGM0(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BIGM1(NSVD_TOTAL*NSVD_TOTAL,NSVD_u))
  ALLOCATE(VEC0(NSVD_TOTAL))
  ALLOCATE(VEC1(NSVD_TOTAL-NSVD))
  ALLOCATE(BMat(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat0(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat01(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat1(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(DETWEI(NGI))
  ALLOCATE(UD(NGI,NSVD))
  ALLOCATE(VD(NGI,NSVD))
  ALLOCATE(WD(NGI,NSVD))
  ALLOCATE(UD_mean(NGI))
  ALLOCATE(VD_mean(NGI))
  ALLOCATE(WD_mean(NGI))
  ALLOCATE(XD(NGI))
  ALLOCATE(YD(NGI))
  ALLOCATE(ZD(NGI))
  ALLOCATE(NX(NLOC,NGI))
  ALLOCATE(NY(NLOC,NGI))
  ALLOCATE(NZ(NLOC,NGI))
  ALLOCATE(MX(NLOC,NGI))
  ALLOCATE(MY(NLOC,NGI))
  ALLOCATE(MZ(NLOC,NGI))
  ALLOCATE(WEI_u(NGI))
  ALLOCATE(WEI_v(NGI))
  ALLOCATE(WEI_w(NGI))
  ALLOCATE(WEI_p(NGI))
  ALLOCATE(WEI_ux(NGI))
  ALLOCATE(WEI_uy(NGI))
  ALLOCATE(WEI_uz(NGI))
  ALLOCATE(WEI_vx(NGI))
  ALLOCATE(WEI_vy(NGI))
  ALLOCATE(WEI_vz(NGI))
  ALLOCATE(WEI_wx(NGI))
  ALLOCATE(WEI_wy(NGI))
  ALLOCATE(WEI_wz(NGI))
  ALLOCATE(NEWCOEF(NSVD_TOTAL))
  ALLOCATE(NEWCOEF1(NSVD_TOTAL-NSVD))
  ALLOCATE(CORIOLIS(NGI))

  ALLOCATE(SOURCX(NONODS))
  ALLOCATE(SOURCY(NONODS))


!**************************************************!
!       Initialise all matrices                    !
!**************************************************!
  BIGM =0.0
  VEC =0.0
  ReMat=0.0
  BIGM0=0.0
  BIGM1=0.0
  VEC0=0.0
  VEC1=0.0
  BMat=0.0
  BMat0=0.0
  BMat01=0.0
  BMat1=0.0
  DETWEI=0.0
  UD=0.0
  VD=0.0
  WD=0.0
  UD_mean=0.0
  VD_mean=0.0
  WD_mean=0.0
  XD=0.0
  YD=0.0
  ZD=0.0
  WEI_u=0.0
  WEI_v=0.0
  WEI_w=0.0
  WEI_p=0.0
  WEI_ux=0.0
  WEI_uy=0.0
  WEI_uz=0.0
  WEI_vx=0.0
  WEI_vy=0.0
  WEI_vz=0.0
  WEI_wx=0.0
  WEI_wy=0.0
  WEI_wz=0.0
  NEWCOEF=0.0
  NEWCOEF1=0.0
  CORIOLIS=0.0


  NX=0.0
  NY=0.0
  NZ=0.0
  MX=0.0
  MY=0.0
  MZ=0.0

  SOURCX=0.0
  SOURCY=0.0

           if(IFWIND) then
           ! WIND STRESS
           !-------------
           CALL WINDYSOURCE(SOURCX,SOURCY,X,Y,NONODS,NSVD)
           endif


  R(1:7*NSVD)=0.0


! Calculate the ML:

!         CALL RCLEAR(ML,NONODS)
!         DO  ELE=1,TOTELE
! ! Calcultae DETWEI,VOLUME...
!           CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,     &
!                    N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,     &
!                    NX,NY,NZ) 
!          DO ILOC=1,NLOC
!           IGL = NDGLNO((ELE-1)*NLOC + ILOC)
!             DO JLOC=1,NLOC
!               JGL=NDGLNO((ELE-1)*NLOC+JLOC) 
!               RNN=0.
!               DO GI=1,NGI
!                 RNN =RNN +N(ILOC,GI)*N(JLOC,GI) *DETWEI(GI)
!               END DO
!               ML(IGL)=ML(IGL)+RNN
!             ENDDO
!           ENDDO
!         ENDDO
!
!         VOL=0.0
!         DO II=1,NONODS
!            VOL=VOL+ML(II)
!         ENDDO
!         ewrite(3,*) 'vol=',vol

! give the initial COEF..

! for 2D
  DO II=2*NSVD_U+1,3*NSVD_U
     COEF(0,II)=0.0
  ENDDO

!     LEFTSVD(2*nonods+1:3*nonods,:)=0.0
     SMEAN(2*mxnods+1:3*mxnods)=0.0

     NEWCOEF(:)=COEF(0,:)

  THETA=0.5
  ONEMINTH=1.0-THETA

  ewrite(3,*) 'nsvd,nvar,nsvd',  nsvd,nvar,nsvd

  ! clean the matrix BIM and left vector VEC
  !-----------------------------------------

  ! define the row and colomn
  !---------------------------
!  DO II=1,NSVD
!     CENTRM(II)=(II-1)*NSVD+II
!     DO JJ=1,NSVD
!        COLM( (II-1)*NSVD+JJ )=JJ
!     ENDDO
!  ENDDO
!
!  DO II=1,NSVD+1
!     FINDRM(II)=(II-1)*NSVD+1
!  ENDDO

!******************************************************!
! Define the pointer for storing data in the matrices  !
!******************************************************!

  NCOLM=NSVD*NSVD
  NCOLM1=NSVD_U*NSVD_PHI
  IBL11=0
  IBL12=1*NCOLM
  IBL13=2*NCOLM
  IBL14=3*NCOLM

  IBL21=3*NCOLM+NCOLM1
  IBL22=IBL21+NCOLM
  IBL23=IBL21+2*NCOLM
  IBL24=IBL21+3*NCOLM

  IBL31=IBL24+NCOLM1
  IBL32=IBL31+NCOLM
  IBL33=IBL31+2*NCOLM
  IBL34=IBL31+3*NCOLM

  IBL41=IBL34+NCOLM1
  IBL42=IBL41+NCOLM1
  IBL43=IBL41+2*NCOLM1
  IBL44=IBL41+3*NCOLM1

!GO TO 400
!  DO 100 NSTEP=1,NTIME
     ewrite(3,*) '11'
     ! setup the reduced model for BETA_ii(t)
     ! the original momentum equations are multiplied by the POD basic vectors 
     ! PHI_ii^u(x),PHI_ii^v(x),PHI_ii^w(x), and continuity equation is multiplied by PHI_ii^p(x)
     ! and then integrate them over the whole domain
     !  PHI_ii^p*{U_x+V_y+W_z}
     ! =PHI_ii{[dU_mean(X)/dx+COEF(t)*(d/dx)*SUM(PHI_ii^u)]+....}
     !------------------------------------------------------------------------------------------
     !
     DO 30 ELE=1,TOTELE
           !                  print *,'ele=',ele
           !  ewrite(3,*)  TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC
           !  ewrite(3,*)  D3,DCYL
           ! Calculate DETWEI,NX,NY,NZ for element ELE
           CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,      &    
                N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
                NX,NY,NZ) 

           ! Calculate DETWEI,MX,MY,MZ for element ELE (for pressure)
           CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,MLOC,NGI,      &    
                M,MLX,MLY,MLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
                MX,MY,MZ) 
           DO 20 ISVD=1,NSVD

           ! CALCULATE THE WEIGHTING FUNCTION(In POD, that is the POD basic function) at the Gaussian points
           ! e.g., in FEM. 
           !WEI_u =[N(GI)*PHI(X)_ii^u], 
           !here , N(GI)=the basic function, PHI(X)_ii^u is the POD basic vectors
           !WEI_ux=[NX(GI)*PHI(X)_ii^u],
           !WEI_uy=[NY(GI)*PHI(X)_ii^u],
           !WEI_uz=[NZ(GI)*PHI(X)_ii^u],
           !WEI_v =[N(GI)*PHI(X)_ii^v], 
           !WEI_vx=[NX(GI)*PHI(X)_ii^v],
           !WEI_vy=[NY(GI)*PHI(X)_ii^v],
           !WEI_vz=[NZ(GI)*PHI(X)_ii^v],
           !WEI_w =[N(GI)*PHI(X)_ii^w], 
           !WEI_wx=[NX(GI)*PHI(X)_ii^w],
           !WEI_wy=[NY(GI)*PHI(X)_ii^w],
           !WEI_wz=[NZ(GI)*PHI(X)_ii^w],
           !-------------------------------
          DO GI=1,NGI
              WEI_u(GI)=0.0
              WEI_v(GI)=0.0
              WEI_w(GI)=0.0
              WEI_ux(GI)=0.0
              WEI_uy(GI)=0.0
              WEI_uz(GI)=0.0
              WEI_vx(GI)=0.0
              WEI_vy(GI)=0.0
              WEI_vz(GI)=0.0
              WEI_wx(GI)=0.0
              WEI_wy(GI)=0.0
              WEI_wz(GI)=0.0

              DO ILOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+ILOC)
                 
                 WEI_u(GI)=WEI_u(GI)+LEFTSVD(JGL,ISVD)*N(ILOC,GI)
                 WEI_v(GI)=WEI_v(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*N(ILOC,GI)
                 WEI_w(GI)=WEI_w(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*N(ILOC,GI)
                 WEI_ux(GI)=WEI_ux(GI)+LEFTSVD(JGL,ISVD)*NX(ILOC,GI)
                 WEI_uy(GI)=WEI_uy(GI)+LEFTSVD(JGL,ISVD)*NY(ILOC,GI)
                 WEI_uz(GI)=WEI_uz(GI)+LEFTSVD(JGL,ISVD)*NZ(ILOC,GI)
                 WEI_vx(GI)=WEI_vx(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*NX(ILOC,GI)
                 WEI_vy(GI)=WEI_vy(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*NY(ILOC,GI)
                 WEI_vz(GI)=WEI_vz(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*NZ(ILOC,GI)
                 WEI_wx(GI)=WEI_wx(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*NX(ILOC,GI)
                 WEI_wy(GI)=WEI_wy(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*NY(ILOC,GI)
                 WEI_wz(GI)=WEI_wz(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*NZ(ILOC,GI)
              ENDDO
           ENDDO

           ! Weight for pressure
           !WEI_p(GI) =N(GI)*PHI(X)_ii^p
           !---------------
!           IF(ISVD.LE.NSVD_PHI) THEN
              DO GI=1,NGI
                 WEI_p(GI)=0.0
                 DO ILOC=1,MLOC
                    IGLP=PNDGLN((ELE-1)*MLOC+ILOC)
                    WEI_p(GI)=WEI_p(GI)+LEFTSVD(3*MXNODS+IGLP,ISVD)*M(ILOC,GI)
                 ENDDO
              ENDDO
!           ENDIF


           !PHI_jj*U=PHI_jj*[U_mean+SUM_jj=1^nsvd(COEF*PHI_ii)]=PHI_jj*U_mean+COEF(t,jj)
           ! UD(GI)=COEF(t,jj)+phi_jj*[N(GI)*U(X)_mean]  
           ! VD(GI)=COEF(t,jj)+phi_jj*[N(GI)*V(X)_mean]  
           ! WD(GI)=COEF(t,jj)+phi_jj*[N(GI)*W(X)_mean]  
           !-----------------------------------------------
           DO GI=1,NGI
              UD(GI,:)=0.0
              VD(GI,:)=0.0
              WD(GI,:)=0.0
              UD_mean(GI)=0.0
              VD_mean(GI)=0.0
              WD_mean(GI)=0.0
              XD(GI)=0.0
              YD(GI)=0.0
              ZD(GI)=0.0
              CORIOLIS(GI)=0.0
              

              DO JLOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+JLOC)

                 UD_mean(GI)=UD_mean(GI)+SMEAN(JGL)*N(JLOC,GI)
                 VD_mean(GI)=VD_mean(GI)+SMEAN(1*MXNODS+JGL)*N(JLOC,GI)
                 WD_mean(GI)=WD_mean(GI)+SMEAN(2*MXNODS+JGL)*N(JLOC,GI)
                 XD(GI)=XD(GI) + N(JLOC,GI)*X(JGL)
                 YD(GI)=YD(GI) + N(JLOC,GI)*Y(JGL)
                 ZD(GI)=ZD(GI) + N(JLOC,GI)*Z(JGL)
                  DO JSVD = 1,NSVD
                     UD(GI,JSVD)=UD(GI,JSVD)+LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                     VD(GI,JSVD)=VD(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
                     WD(GI,JSVD)=WD(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
                  ENDDO

              ENDDO
              CORIOLIS(GI) = 2.0*FUNOME(XD(GI),YD(GI),ZD(GI))
           ENDDO


           !Calculate
           !------
           !PHI_jj^p*[U_x_mean+V_y_mean+W_z_mean]
           !------------------
           !PHI_jj^u*[U_x_mean+UD*U_x_mean+VD*U_y_mean+WD*U_z_mean-FV_mean+P_x_mean]
           !+D(PHI_jj^u)/DX*DU_mean/DX+D(PHI_jj^u)/DY*DU_mean/DY+D(PHI_jj^u)/DZ*DU_mean/DZ
           !-----------------
           !PHI_jj^v*[V_x_mean+UD*V_x_mean+VD*V_y_mean+WD*V_z_mean+FU_mean+P_y_mean]
           !+D(PHI_jj^v)/DX*DV_mean/DX+D(PHI_jj^v)/DY*DV_mean/DY+D(PHI_jj^v)/DZ*DV_mean/DZ
           !-----------------
           !PHI_jj^w*[W_x_mean+UD*W_x_mean+VD*W_y_mean+WD*W_z_mean+P_z_mean]
           !+D(PHI_jj^w)/DX*DW_mean/DX+D(PHI_jj^w)/DY*DW_mean/DY+D(PHI_jj^w)/DZ*DW_mean/DZ
           !---------------------------------

              FU_mean =0.0
              FV_mean =0.0
              U_x_mean =0.0
              V_y_mean =0.0
              W_Z_mean =0.0
              UU_x_mean=0.0
              VU_y_mean=0.0
              WU_z_mean=0.0
              UV_x_mean=0.0
              VV_y_mean=0.0
              WV_z_mean=0.0
              UW_x_mean=0.0
              VW_y_mean=0.0
              WW_z_mean=0.0
              U_xx_mean=0.0
              U_yy_mean=0.0
              U_zz_mean=0.0
              V_xx_mean=0.0
              V_yy_mean=0.0
              V_zz_mean=0.0
              W_xx_mean=0.0
              W_yy_mean=0.0
              W_zz_mean=0.0
              P_x_mean =0.0
              P_y_mean =0.0
              P_z_mean =0.0

           DO JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
              DO GI=1,NGI
                 
                 FU_mean =FU_mean+WEI_v(GI)*CORIOLIS(GI)*DETWEI(GI)*SMEAN(JGL)*N(JLOC,GI)
                 FV_mean =FV_mean-WEI_u(GI)*CORIOLIS(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*N(JLOC,GI)
                 U_x_mean =U_x_mean+WEI_p(GI)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                 V_y_mean =V_y_mean+WEI_p(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 W_Z_mean =W_Z_mean+WEI_p(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)
                 UU_x_mean=UU_x_mean+WEI_u(GI)*UD_mean(GI)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                 VU_y_mean=VU_y_mean+WEI_u(GI)*VD_mean(GI)*DETWEI(GI)*SMEAN(JGL)*NY(JLOC,GI)
                 WU_z_mean=WU_z_mean+WEI_u(GI)*WD_mean(GI)*DETWEI(GI)*SMEAN(JGL)*NZ(JLOC,GI)
                 UV_x_mean=UV_x_mean+WEI_v(GI)*UD_mean(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                 VV_y_mean=VV_y_mean+WEI_v(GI)*VD_mean(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 WV_z_mean=WV_z_mean+WEI_v(GI)*WD_mean(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)
                 UW_x_mean=UW_x_mean+WEI_w(GI)*UD_mean(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                 VW_y_mean=VW_y_mean+WEI_w(GI)*VD_mean(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                 WW_z_mean=WW_z_mean+WEI_w(GI)*WD_mean(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)

                 U_xx_mean=U_xx_mean+WEI_ux(GI)*MUPTXX(JGL)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                 U_yy_mean=U_yy_mean+WEI_uy(GI)*MUPTYY(JGL)*DETWEI(GI)*SMEAN(JGL)*NY(JLOC,GI)
                 U_zz_mean=U_zz_mean+WEI_uz(GI)*MUPTZZ(JGL)*DETWEI(GI)*SMEAN(JGL)*NZ(JLOC,GI)
                 V_xx_mean=V_xx_mean+WEI_vx(GI)*MUPTXX(JGL)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                 V_yy_mean=V_yy_mean+WEI_vy(GI)*MUPTYY(JGL)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 V_zz_mean=V_zz_mean+WEI_vz(GI)*MUPTZZ(JGL)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)
                 W_xx_mean=W_xx_mean+WEI_wx(GI)*MUPTXX(JGL)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                 W_yy_mean=W_yy_mean+WEI_wy(GI)*MUPTYY(JGL)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                 W_zz_mean=W_zz_mean+WEI_wz(GI)*MUPTZZ(JGL)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)
!                 U_zz_mean=0.0
!                 V_zz_mean=0.0
!                 W_zz_mean=0.0
              ENDDO
           ENDDO

              ! Calculate the average pressure p_x_mean, P_y_mean, p_z_mean

           DO JLOC=1,MLOC
              JGLP=PNDGLN((ELE-1)*MLOC+JLOC)
              DO GI=1,NGI
                 P_x_mean = P_x_mean+WEI_u(GI)*DETWEI(GI)*SMEAN(3*MXNODS+JGLP)*MX(JLOC,GI)
                 P_y_mean = P_y_mean+WEI_v(GI)*DETWEI(GI)*SMEAN(3*MXNODS+JGLP)*MY(JLOC,GI)
                 P_z_mean = P_z_mean+WEI_w(GI)*DETWEI(GI)*SMEAN(3*MXNODS+JGLP)*MZ(JLOC,GI)
              ENDDO
           ENDDO

!           FV_mean=0.0
!           FU_mean=0.0
           VEC0(ISVD)= VEC0(ISVD)-DT*(UU_x_mean+VU_y_mean+WU_z_mean+FV_mean+P_x_mean+                 &
                U_xx_mean+U_yy_mean+U_zz_mean)
           VEC0(1*NSVD+ISVD)= VEC0(1*NSVD+ISVD)-DT*(UV_x_mean+VV_y_mean+WV_z_mean+FU_mean+P_y_mean+   &
                V_xx_mean+V_yy_mean+V_zz_mean)
           VEC0(2*NSVD+ISVD)= VEC0(2*NSVD+ISVD)-DT*(UW_x_mean+VW_y_mean+WW_z_mean+P_z_mean+           &
                W_xx_mean+W_yy_mean+W_zz_mean)

           IF(ISVD.LE.NSVD_PHI) THEN
              VEC(3*NSVD+ISVD)= VEC(3*NSVD+ISVD)-U_x_mean-V_y_mean-W_z_mean
! for 2D
!              VEC0(3*NSVD+ISVD)= VEC0(3*NSVD+ISVD)-U_x_mean-V_y_mean
           ENDIF

           if(IFWIND) then
           ! WIND STRESS
           !-------------
!           CALL WINDYSOURCE(SOURCX,SOURCY,X,Y,NONODS,NSVD)
           DO ILOC=1,NLOC
              IGL = NDGLNO((ELE-1)*NLOC + ILOC)
              DO GI=1,NGI
                 VEC0(ISVD) = VEC0(ISVD) +DT*WEI_u(GI)*DETWEI(GI)*SOURCX(IGL)*N(ILOC,GI)
              END DO
           ENDDO
           endif

!           ewrite(3,*) 'VEC(ISVD)',VEC(ISVD)

           !Clculate
           !PHI_jj^u*U; PHI_jj^v*V; PHI_jj^w*W
           !------
           !PHI_jj^p*[U_x+V_y+W_z]
           !------------------
           !PHI_jj^u*[U_x+UD*U_x+VD*U_y+WD*U_z+FV+P_x]
           !+D(PHI_jj^u)/DX*DU/DX+D(PHI_jj^u)/DY*DU/DY+D(PHI_jj^u)/DZ*DU/DZ
           !-----------------
           !PHI_jj^v*[V_x+UD*V_x+VD*V_y+WD*V_z-FU+P_y]
           !+D(PHI_jj^v)/DX*DV/DX+D(PHI_jj^v)/DY*DV/DY+D(PHI_jj^v)/DZ*DV/DZ
           !-----------------
           !PHI_jj^w*[W_x+UD*W_x+VD*W_y+WD*W_z+P_z]
           !+D(PHI_jj^w)/DX*DW/DX+D(PHI_jj^w)/DY*DW/DY+D(PHI_jj^w)/DZ*DW/DZ
           !---------------------------------
           !   ewrite(3,*) 'vec1',VEC(ISVD)
           DO 40 JSVD=1,NSVD

              COUNT=(ISVD-1)*NSVD+JSVD
              COUNT1=(ISVD-1)*NSVD_PHI+JSVD

                 RN1=0.0
                 RN2=0.0
                 RN3=0.0
                 U_x=0.0
                 V_y=0.0
                 W_Z=0.0

                 UDmean_Ux=0.0
                 VDmean_Uy=0.0
                 WDmean_Uz=0.0
                 UDmean_Vx=0.0
                 VDmean_Vy=0.0
                 WDmean_Vz=0.0
                 UDmean_Wx=0.0
                 VDmean_Wy=0.0
                 WDmean_Wz=0.0
                 UD_Ux=0.0
                 VD_Uy=0.0
                 WD_Uz=0.0
                 UD_Vx=0.0
                 VD_Vy=0.0
                 WD_Vz=0.0
                 UD_Wx=0.0
                 VD_Wy=0.0
                 WD_Wz=0.0
                 Umeanx_UD=0.0
                 Umeany_VD=0.0
                 Umeanz_WD=0.0
                 Vmeanx_UD=0.0
                 Vmeany_VD=0.0
                 Vmeanz_WD=0.0
                 Wmeanx_UD=0.0
                 Wmeany_VD=0.0
                 Wmeanz_WD=0.0

                 U_xx=0.0
                 U_yy=0.0
                 U_zz=0.0
                 V_xx=0.0
                 V_yy=0.0
                 V_zz=0.0
                 W_xx=0.0
                 W_yy=0.0
                 W_zz=0.0
                 FU=0.0
                 FV=0.0
                 P_x=0.0
                 P_y=0.0
                 P_z=0.0

              DO JLOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+JLOC)
                 DO GI=1,NGI
                    RN1=RN1+WEI_u(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    RN2=RN2+WEI_v(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    RN3=RN3+WEI_w(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    IF(ISVD.LE.NSVD_PHI) THEN
                       U_x=U_x+WEI_p(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                       V_y=V_y+WEI_p(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                       W_Z=W_Z+WEI_p(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    ENDIF
                    UDmean_Ux=UDmean_Ux+WEI_u(GI)*UD_mean(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    VDmean_Uy=VDmean_Uy+WEI_u(GI)*VD_mean(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    WDmean_Uz=WDmean_Uz+WEI_u(GI)*WD_mean(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                    UDmean_Vx=UDmean_Vx+WEI_v(GI)*UD_mean(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VDmean_Vy=VDmean_Vy+WEI_v(GI)*VD_mean(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WDmean_Vz=WDmean_Vz+WEI_v(GI)*WD_mean(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    UDmean_Wx=UDmean_Wx+WEI_w(GI)*UD_mean(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VDmean_Wy=VDmean_Wy+WEI_w(GI)*VD_mean(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WDmean_Wz=WDmean_Wz+WEI_w(GI)*WD_mean(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)



                    Umeanx_UD=Umeanx_UD+WEI_u(GI)*UD(GI,JSVD)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                    Umeany_VD=Umeany_VD+WEI_u(GI)*VD(GI,JSVD)*DETWEI(GI)*SMEAN(JGL)*NY(JLOC,GI)
                    Umeanz_WD=Umeanz_WD+WEI_u(GI)*WD(GI,JSVD)*DETWEI(GI)*SMEAN(JGL)*NZ(JLOC,GI)
                    Vmeanx_UD=Vmeanx_UD+WEI_v(GI)*UD(GI,JSVD)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                    Vmeany_VD=Vmeany_VD+WEI_v(GI)*VD(GI,JSVD)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                    Vmeanz_WD=Vmeanz_WD+WEI_v(GI)*WD(GI,JSVD)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)
                    Wmeanx_UD=Wmeanx_UD+WEI_w(GI)*UD(GI,JSVD)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                    Wmeany_VD=Wmeany_VD+WEI_w(GI)*VD(GI,JSVD)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                    Wmeanz_WD=Wmeanz_WD+WEI_w(GI)*WD(GI,JSVD)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)

!
                    UD_Ux=UD_Ux+WEI_u(GI)*UD(GI,JSVD)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    VD_Uy=VD_Uy+WEI_u(GI)*VD(GI,JSVD)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    WD_Uz=WD_Uz+WEI_u(GI)*WD(GI,JSVD)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                    UD_Vx=UD_Vx+WEI_v(GI)*UD(GI,JSVD)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VD_Vy=VD_Vy+WEI_v(GI)*VD(GI,JSVD)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WD_Vz=WD_Vz+WEI_v(GI)*WD(GI,JSVD)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    UD_Wx=UD_Wx+WEI_w(GI)*UD(GI,JSVD)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VD_Wy=VD_Wy+WEI_w(GI)*VD(GI,JSVD)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WD_Wz=WD_Wz+WEI_w(GI)*WD(GI,JSVD)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)


                    U_xx=U_xx+WEI_ux(GI)*MUPTXX(JGL)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    U_yy=U_yy+WEI_uy(GI)*MUPTYY(JGL)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    U_zz=U_zz+WEI_uz(GI)*MUPTZZ(JGL)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                    V_xx=V_xx+WEI_vx(GI)*MUPTXX(JGL)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    V_yy=V_yy+WEI_vy(GI)*MUPTYY(JGL)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    V_zz=V_zz+WEI_vz(GI)*MUPTZZ(JGL)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    W_xx=W_xx+WEI_wx(GI)*MUPTXX(JGL)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    W_yy=W_yy+WEI_wy(GI)*MUPTYY(JGL)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    W_zz=W_zz+WEI_wz(GI)*MUPTZZ(JGL)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    FU = FU+WEI_v(GI)*CORIOLIS(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    FV = FV-WEI_u(GI)*CORIOLIS(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

!                 U_zz=0.0
!                 V_zz=0.0
!                 W_zz=0.0
                    
                    DO JJ=1,NSVD
                       BIGM1(IBL11+COUNT,JJ)=BIGM1(IBL11+COUNT,JJ)+WEI_u(GI)*DETWEI(GI)*UD(GI,JJ)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                       BIGM1(IBL12+COUNT,JJ)=BIGM1(IBL12+COUNT,JJ)+WEI_u(GI)*DETWEI(GI)*VD(GI,JJ)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                       BIGM1(IBL13+COUNT,JJ)=BIGM1(IBL13+COUNT,JJ)+WEI_u(GI)*DETWEI(GI)*WD(GI,JJ)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                       
                       BIGM1(IBL21+COUNT,JJ)=BIGM1(IBL21+COUNT,JJ)+WEI_v(GI)*DETWEI(GI)*UD(GI,JJ)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                       BIGM1(IBL22+COUNT,JJ)=BIGM1(IBL22+COUNT,JJ)+WEI_v(GI)*DETWEI(GI)*VD(GI,JJ)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                       BIGM1(IBL23+COUNT,JJ)=BIGM1(IBL23+COUNT,JJ)+WEI_v(GI)*DETWEI(GI)*WD(GI,JJ)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                       
                       BIGM1(IBL31+COUNT,JJ)=BIGM1(IBL31+COUNT,JJ)+WEI_w(GI)*DETWEI(GI)*UD(GI,JJ)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                       BIGM1(IBL32+COUNT,JJ)=BIGM1(IBL32+COUNT,JJ)+WEI_w(GI)*DETWEI(GI)*VD(GI,JJ)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                       BIGM1(IBL33+COUNT,JJ)=BIGM1(IBL33+COUNT,JJ)+WEI_w(GI)*DETWEI(GI)*WD(GI,JJ)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    ENDDO


                 ENDDO
              ENDDO

              IF(JSVD.LE.NSVD_PHI) THEN
                 DO JLOC=1,MLOC
                    JGLP=PNDGLN((ELE-1)*MLOC+JLOC)
                    DO GI=1,NGI
                       P_x=P_x+WEI_u(GI)*DETWEI(GI)*LEFTSVD(3*MXNODS+JGLP,JSVD)*MX(JLOC,GI)
                       P_y=P_y+WEI_v(GI)*DETWEI(GI)*LEFTSVD(3*MXNODS+JGLP,JSVD)*MY(JLOC,GI)
                       P_z=P_z+WEI_w(GI)*DETWEI(GI)*LEFTSVD(3*MXNODS+JGLP,JSVD)*MZ(JLOC,GI)
                    ENDDO
                 ENDDO
              ENDIF
              ! Find count. ***************************************
              !                 LOWER=FINDRM(ISVD) 
              !                 UPPER=FINDRM(ISVD+1)-1
              ! 7000                 CONTINUE
              !                 INUM=LOWER+(UPPER-LOWER+1)/2 
              !                 IF(JSVD.GE.COLM(INUM) )  THEN 
              !                    LOWER=INUM
              !                 ELSE
              !                    UPPER=INUM
              !                 ENDIF
              !                 IF(UPPER-LOWER.LE.1) THEN
              !                    IF(JSVD.EQ.COLM(LOWER)) THEN
              !                       COUNT=LOWER
              !                    ELSE
              !                       COUNT=UPPER
              !                    ENDIF
              !                    GOTO 9000
              !                 ENDIF
              !                 GOTO 7000
              ! 9000                 CONTINUE

              !  Add into matrix and rhs...
              !--------------------------
!              FV=0.0
!              FU=0.0
 
              BIGM0(IBL11+COUNT)=BIGM0(IBL11+COUNT)+RN1+DT*THETA*(UDmean_Ux+VDmean_Uy+WDmean_Uz+U_xx+U_yy+U_zz)
              BIGM0(IBL12+COUNT)=BIGM0(IBL12+COUNT)+DT*THETA*FV
              BIGM0(IBL13+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                 BIGM0(IBL14+COUNT1)=BIGM0(IBL14+COUNT1)+DT*THETA*P_x
              ENDIF

              BIGM0(IBL21+COUNT)=BIGM0(IBL21+COUNT)+DT*THETA*FU
              BIGM0(IBL22+COUNT)=BIGM0(IBL22+COUNT)+RN2+DT*THETA*(UDmean_Vx+VDmean_Vy+WDmean_Vz+V_xx+V_yy+V_zz)
              BIGM0(IBL23+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                 BIGM0(IBL24+COUNT1)=BIGM0(IBL24+COUNT1)+DT*THETA*P_y
              ENDIF

              BIGM0(IBL31+COUNT)=0.0
              BIGM0(IBL32+COUNT)=0.0
              BIGM0(IBL33+COUNT)=BIGM0(IBL33+COUNT)+RN3+DT*THETA*(UDmean_Wx+VDmean_Wy+WDmean_Wz+W_xx+W_yy+W_zz)
              IF(JSVD.LE.NSVD_PHI) THEN
                 BIGM0(IBL34+COUNT1)=BIGM0(IBL34+COUNT1)+DT*THETA*P_z
              ENDIF

              IF(ISVD.LE.NSVD_PHI) THEN
                 BIGM0(IBL41+COUNT)=BIGM0(IBL41+COUNT)+U_x
                 BIGM0(IBL42+COUNT)=BIGM0(IBL42+COUNT)+V_y
! for 2D
                 BIGM0(IBL43+COUNT)=BIGM(IBL43+COUNT)+W_z
                 IF(JSVD.LE.NSVD_PHI) THEN
                    BIGM0(IBL44+COUNT1)=0.0
                 ENDIF
              ENDIF

              BMat0(IBL11+COUNT)=BMat0(IBL11+COUNT)+RN1-DT*ONEMINTH*(UDmean_Ux+VDmean_Uy+WDmean_Uz+U_xx+U_yy+U_zz)
              BMat0(IBL12+COUNT)=BMat0(IBL12+COUNT)-DT*ONEMINTH*FV
              BMat0(IBL13+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                 BMat0(IBL14+COUNT1)=BMat0(IBL14+COUNT1)-DT*ONEMINTH*P_x
              ENDIF

 
              BMat0(IBL21+COUNT)=BMat0(IBL21+COUNT)-DT*ONEMINTH*FU
              BMat0(IBL22+COUNT)=BMat0(IBL22+COUNT)+RN2-DT*ONEMINTH*(UDmean_Vx+VDmean_Vy+WDmean_Vz+V_xx+V_yy+V_zz)
              BMat0(IBL23+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                 BMat0(IBL24+COUNT1)=BMat0(IBL24+COUNT1)-DT*ONEMINTH*P_y
              ENDIF
              BMat0(IBL31+COUNT)=0.0
              BMat0(IBL32+COUNT)=0.0
              BMat0(IBL33+COUNT)=BMat0(IBL33+COUNT)+RN3-DT*ONEMINTH*(UDmean_Wx+VDmean_Wy+WDmean_Wz+W_xx+W_yy+W_zz)
              IF(JSVD.LE.NSVD_PHI) THEN
                 BMat0(IBL34+COUNT1)=BMat0(IBL34+COUNT1)-DT*ONEMINTH*P_z
              ENDIF

              IF(ISVD.LE.NSVD_PHI) THEN
                 BMat0(IBL41+COUNT)=0.0
                 BMat0(IBL42+COUNT)=0.0
! for 2D
                 BMat0(IBL43+COUNT)=0.0
                 IF(JSVD.LE.NSVD_PHI) THEN
                    BMat0(IBL44+COUNT1)=0.0
                 ENDIF
              ENDIF
!

              BMat0(IBL11+COUNT)=BMat0(IBL11+COUNT)-DT*Umeanx_UD
              BMat0(IBL12+COUNT)=BMat0(IBL12+COUNT)-DT*Umeany_VD
              BMat0(IBL13+COUNT)=BMat0(IBL13+COUNT)-DT*Umeanz_WD
!!!
              BMat0(IBL21+COUNT)=BMat0(IBL21+COUNT)-DT*Vmeanx_UD
              BMat0(IBL22+COUNT)=BMat0(IBL22+COUNT)-DT*Vmeany_VD
              BMat0(IBL23+COUNT)=BMat0(IBL23+COUNT)-DT*Vmeanz_WD


              BMat0(IBL31+COUNT)=BMat0(IBL31+COUNT)-DT*Wmeanx_UD
              BMat0(IBL32+COUNT)=BMat0(IBL32+COUNT)-DT*Wmeany_VD
              BMat0(IBL33+COUNT)=BMat0(IBL33+COUNT)-DT*Wmeanz_WD

40            CONTINUE
!              ewrite(3,*) 'vec2',VEC(ISVD)
!              ewrite(3,*) 'BIGM(COUNT)',BIGM(COUNT)
!              if(ISVD.EQ.2) stop

20            CONTINUE
!              print *,'1',BIGM(COUNT)
30            CONTINUE


400 CONTINUE
open(1,file='PODBIGM0.dat')
  read(1,*) (BIGM0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
close(1)
open(1,file='PODBMat0.dat')
  read(1,*) (BMat0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
close(1)
open(1,file='PODVEC0.dat')
  read(1,*) (VEC0(II),II=1,NSVD_TOTAL)
close(1)
open(1,file='PODBIGM1.dat')
 DO JJ=1,NSVD_U
  read(1,*) (BIGM1(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
 ENDDO
close(1)
!stop 10



   DO 100 NSTEP=1,ntime

      BIGM(:)=BIGM0(:)
      BMat(:)=BMat0(:)
      VEC(:)=VEC0(:)
      open(1,file='bigm44a.dat')
        write(1,*) (BIGM(IBL44+ii),ii=1,NSVD_PHI*NSVD_PHI)
      close(1)

      DO ISVD=1,NSVD

         DO JSVD=1,NSVD
            COUNT=(ISVD-1)*NSVD+JSVD
            COUNT1=(ISVD-1)*NSVD_PHI+JSVD
            DO JJ=1,NSVD
               BIGM(IBL11+COUNT)=BIGM(IBL11+COUNT)+DT*THETA*( BIGM1(IBL11+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )

               BIGM(IBL22+COUNT)=BIGM(IBL22+COUNT)+DT*THETA*( BIGM1(IBL21+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )

               BIGM(IBL33+COUNT)=BIGM(IBL33+COUNT)+DT*THETA*( BIGM1(IBL31+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )


               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)-DT*ONEMINTH*( BIGM1(IBL11+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )


               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)-DT*ONEMINTH*( BIGM1(IBL21+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )

               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)-DT*ONEMINTH*( BIGM1(IBL31+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )
            ENDDO


               VEC(ISVD)=VEC(ISVD)+BMat(IBL11+COUNT)*NEWCOEF(JSVD)+                       &
                    BMat(IBL12+COUNT)*NEWCOEF(1*NSVD+JSVD)+BMat(IBL13+COUNT)*NEWCOEF(2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(ISVD)=VEC(ISVD)+BMat0(IBL14+COUNT1)*NEWCOEF(3*NSVD+JSVD)
               ENDIF

               VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat(IBL21+COUNT)*NEWCOEF(JSVD)+         &
                    BMat(IBL22+COUNT)*NEWCOEF(1*NSVD+JSVD)+BMat(IBL23+COUNT)*NEWCOEF(2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat0(IBL24+COUNT1)*NEWCOEF(3*NSVD+JSVD)
               ENDIF

               VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat(IBL31+COUNT)*NEWCOEF(JSVD)+         &
                    BMat(IBL32+COUNT)*NEWCOEF(1*NSVD+JSVD)+BMat(IBL33+COUNT)*NEWCOEF(2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat0(IBL34+COUNT1)*NEWCOEF(3*NSVD+JSVD)
               ENDIF
         ENDDO

      ENDDO

              if(.false.) then
              DO kk=1,NVAR
                 DO ii=1,NSVD
                    ftest=0.0
                    DO jj=1,NSVD
                       ReMat( (kk-1)*NSVD+ii,jj)=BIGM( (kk-1)*NVAR*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ReMat( (kk-1)*NSVD+ii,1*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+1*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ReMat( (kk-1)*NSVD+ii,2*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+2*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ReMat( (kk-1)*NSVD+ii,3*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+3*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ftest=ftest+ReMat((kk-1)*NVAR*NSVD+ii,jj)
                    ENDDO
                    
                    if(abs(ftest).lt.1.0e-20) then
!                       ewrite(3,*) '(kk-1)*NVAR*NSVD+ii',kk,ii
!                       stop
                    endif
                 ENDDO
              ENDDO
              
              else

      open(1,file='bigm44.dat')
        write(1,*) IBL44
        write(1,*) (BIGM(IBL44+ii),ii=1,NSVD_PHI*NSVD_PHI)
      close(1)

              ReMat=0.0
              VEC1=0.0
              DO kk=1,NVAR
                 DO ii=1,NSVD
                    DO jj=1,NSVD
                       IF(kk.le.2) then
                          ReMat( (kk-1)*NSVD+ii,jj)=BIGM( (kk-1)*NSVD*NSVD_TOTAL+(ii-1)*NSVD+jj)
                          ReMat( (kk-1)*NSVD+ii,1*NSVD+jj)= BIGM( (kk-1)*NSVD*NSVD_TOTAL+1*NSVD*NSVD+(ii-1)*NSVD+jj)
                          !     ReMat( (kk-1)*NVAR*NSVD+ii,2*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+2*NSVD*NSVD+(ii-1)*NSVD+jj)
                          IF(jj.LE.NSVD_PHI) THEN
                             ReMat( (kk-1)*NSVD+ii,2*NSVD+jj)= BIGM( (kk-1)*NSVD*NSVD_TOTAL+3*NSVD*NSVD+(ii-1)*NSVD_PHI+jj)
                          ENDIF
                       ftest=ftest+ReMat((kk-1)*NVAR*NSVD+ii,jj)
                       ELSEIF((kk.eq.4).AND.(ii.LE.NSVD_PHI)) then
                          ReMat( (kk-2)*NSVD+ii,jj)=BIGM( (kk-1)*NSVD*NSVD_TOTAL+(ii-1)*NSVD+jj)
                          ReMat( (kk-2)*NSVD+ii,1*NSVD+jj)= BIGM( (kk-1)*NSVD*NSVD_TOTAL+1*NSVD_PHI*NSVD+(ii-1)*NSVD+jj)
                          !     ReMat( (kk-1)*NVAR*NSVD+ii,2*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+2*NSVD*NSVD+(ii-1)*NSVD+jj)
                          IF(jj.LE.NSVD_PHI) THEN
                             ReMat( (kk-2)*NSVD+ii,2*NSVD+jj)= BIGM(IBL44+(ii-1)*NSVD_PHI+jj)
                          ENDIF
                       ENDIF
                    ENDDO
                    IF(kk.le.2) then
                       VEC1( (kk-1)*NSVD+ii)=VEC( (kk-1)*NSVD+ii)
                    ELSEIF((kk.eq.4).AND.(ii.LE.NSVD_PHI)) then
                       VEC1( (kk-2)*NSVD+ii)=VEC( (kk-1)*NSVD+ii)
                    ENDIF

                 ENDDO
              ENDDO

              endif

   if(.false.) then
              do ii=1,nvar*nsvd
                 ReMat(ii,ii)=ReMat(ii,ii)+vol
                 vec(ii)=vec(ii)+vol
              enddo
              close(1)
              close(2)
!              stop

              open(1,file='Remat2.dat')
              open(2,file='vec2.dat')
              write(1,*) Remat
              write(2,*) vec
              close(1)
              close(2)
              stop

   endif
!              IF(NSTEP.EQ.1) THEN
!              open(1,file='Remat1.dat')
!              open(2,file='vec1.dat')
!              write(1,*) Remat
!              write(2,*) vec
!              ELSEIF(NSTEP.EQ.3) THEN
!              open(1,file='Remat3.dat')
!              open(2,file='vec3.dat')
!              write(1,*) Remat
!              write(2,*) vec
!              STOP 156
!              ENDIF

              open(1,file='ReMatVec5.dat')
!              PRINT *,'BIGM:'
!              DO ISVD=1,NSVD_phi
!                PRINT *,(BIGM(3*NSVD*NSVD_TOTAL+3*NSVD_phi*NSVD+(ISVD-1)*NSVD_phi+JSVD),JSVD=1,NSVD_phi)
!              END DO
 
               Write(1,*) 'Remat:'
              DO ISVD=1,NSVD_TOTAL-NSVD
                write(1,*) (ReMat(ISVD,JSVD),JSVD=1,NSVD_TOTAL-NSVD)
              END DO
              write(1,*) 'matfor pressure'
              DO ISVD=2*nsvd+1,NSVD_TOTAL-NSVD
                write(1,*) (ReMat(ISVD,JSVD),JSVD=2*nsvd+1,NSVD_TOTAL-NSVD)
              END DO
              
               Write(1,*) 'Vec:'
                write(1,*) Vec
                close(1)
!              ewrite(3,*) 'vec',vec
              ! SOLVE matrix equations for gradients of hydrostic pressure
              ! *********DEFINE OPTIONS FOR GMRES**************

              ewrite(3,*) 'VEC'
              ewrite(3,*) 'BIGM'

              ERROR1 = 1.0E-14
              ERROR2 = ERROR1
              NOITS1 = 201
              NOITS2 = NOITS1
              IGUESS = .FALSE.
              MISOUT = 0
              IMM    = 1
              RELAX  = 1.0                
              LEFTP1 = 0
              LEFTP2 = 0
              MINTE1 = 20
              MINTE2 = 20       
              MULPA  = 0
              UZAWA  = 0
              CGSOLQ = 0
              GMRESQ = 1
              CGSQ   = 0 
              !       
              NR2=0
              IP=1
              MISOUT=1 
              NNV=7*NSVD-2*NSVD
 
              ! *********ENDOF OPTIONS FOR GMRES**************
              !       Solev MATRIX X=VEC

!              PRINT *,'BIGM:'
!              DO ISVD=1,NSVD_phi
!                PRINT *,(BIGM(3*NSVD*NSVD_TOTAL+3*NSVD_phi*NSVD+(ISVD-1)*NSVD_phi+JSVD),JSVD=1,NSVD_phi)
!              END DO
! 
!               PRINT *,'Remat:'
!              DO ISVD=1,NSVD
!                PRINT *,(ReMat(ISVD,JSVD),JSVD=1,NSVD_TOTAL-NSVD)
!              END DO
!              ewrite(3,*) 'vec',vec


        open(10, file ='NEWCOEF.dat')
           write(10, *) (NEWCOEF1(ii), ii=1,(nvar-1)*nsvd)
        close(10) 

              open(1,file='ReMatVec6.dat')
!              PRINT *,'BIGM:'
!              DO ISVD=1,NSVD_phi
!                PRINT *,(BIGM(3*NSVD*NSVD_TOTAL+3*NSVD_phi*NSVD+(ISVD-1)*NSVD_phi+JSVD),JSVD=1,NSVD_phi)
!              END DO
 
               Write(1,*) 'Remat:'
              DO ISVD=1,NSVD_TOTAL-NSVD
                write(1,*) (ReMat(ISVD,JSVD),JSVD=1,NSVD_TOTAL-NSVD)
              END DO
              write(1,*) 'matfor pressure'
              DO ISVD=2*nsvd+1,NSVD_TOTAL-NSVD
                write(1,*) (ReMat(ISVD,JSVD),JSVD=2*nsvd+1,NSVD_TOTAL-NSVD)
              END DO
              
               Write(1,*) 'Vec:'
                write(1,*) Vec
                close(1)

!           IF(2*int(nstep/2) .eq. nstep) THEN
              DO II=1,NSVD_TOTAL
                 IF(II.LE.2*NSVD) THEN
                    NEWCOEF(II)=NEWCOEF1(II)
                    IF(NT*int(nstep/NT) .eq. nstep) THEN
                       COEF(int(nstep/NT),II)=NEWCOEF1(II)
                    ENDIF
                 ELSEIF(II.GT.3*NSVD) THEN
                    NEWCOEF(II)=NEWCOEF1(II-NSVD)
                    IF(NT*int(nstep/NT) .eq. nstep) THEN
                       COEF(int(nstep/NT),II)=NEWCOEF1(II-NSVD)
                    ENDIF
                 ELSE
                    NEWCOEF(II)=0.0
                    IF(NT*int(nstep/NT) .eq. nstep) THEN
                       COEF(int(nstep/NT),II)=0.0
                    ENDIF
                 ENDIF
                 IF(NT*int(nstep/NT) .eq. nstep) THEN
                    ewrite(3,*) 'COEF(int(nstep/NT),II)',COEF(int(nstep/NT),II)
                 ENDIF
               ENDDO
!           ENDIF
              
              ewrite(3,*) 'after NSTEP',NSTEP
100           CONTINUE


              open(10, file ='COEF.dat')
              do ii=0,ntime1
                 write(10, *) (COEF(ii,jj), jj=1,nvar*nsvd)
              enddo
              close(10)
              ewrite(3,*) 'before deallocate'
!-------------------------------!
! DEALLOCATE THE LOCAL MATRICES !
!-------------------------------!

  DEALLOCATE(BIGM)
  DEALLOCATE(VEC)
  DEALLOCATE(ReMat)
  DEALLOCATE(BIGM0)
  DEALLOCATE(BIGM1)
  DEALLOCATE(VEC0)
  DEALLOCATE(VEC1)
  DEALLOCATE(BMat)
  DEALLOCATE(BMat0)
  DEALLOCATE(BMat01)
  DEALLOCATE(BMat1)
  DEALLOCATE(DETWEI)
  DEALLOCATE(UD)
  DEALLOCATE(VD)
  DEALLOCATE(WD)
  DEALLOCATE(UD_mean)
  DEALLOCATE(VD_mean)
  DEALLOCATE(WD_mean)
  DEALLOCATE(XD)
  DEALLOCATE(YD)
  DEALLOCATE(ZD)
  DEALLOCATE(NX)
  DEALLOCATE(NY)
  DEALLOCATE(NZ)
  DEALLOCATE(MX)
  DEALLOCATE(MY)
  DEALLOCATE(MZ)
  DEALLOCATE(WEI_u)
  DEALLOCATE(WEI_v)
  DEALLOCATE(WEI_w)
  DEALLOCATE(WEI_p)
  DEALLOCATE(WEI_ux)
  DEALLOCATE(WEI_uy)
  DEALLOCATE(WEI_uz)
  DEALLOCATE(WEI_vx)
  DEALLOCATE(WEI_vy)
  DEALLOCATE(WEI_vz)
  DEALLOCATE(WEI_wx)
  DEALLOCATE(WEI_wy)
  DEALLOCATE(WEI_wz)
  DEALLOCATE(NEWCOEF)
  DEALLOCATE(NEWCOEF1)
  DEALLOCATE(CORIOLIS)
  DEALLOCATE(SOURCX)
  DEALLOCATE(SOURCY)
              ewrite(3,*) 'exit REDUCED_MODEL_DIFFCOEF'
              RETURN
END SUBROUTINE REDUCED_MODEL_DIFFCOEF1


SUBROUTINE REDUCED_MODEL_DIFFCOEF1_PG1(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,X,Y,Z,                       &
     LEFTSVD,SMEAN,COEF,M,MLX,MLY,MLZ,                                    &
     N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
     NTIME1,DT1,D3,DCYL,GEOBAL,SCFACTH0,IFWIND, &
     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO,SMEAN_PG,leftsvd_PG_x,leftsvd_PG_y,  &
     NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM,IFWRITE,GLOITS,LINITS)

  INTEGER, INTENT(IN)  ::NSVD_TOTAL,NSVD,NVAR,NSVD_U,NSVD_PHI
  INTEGER, INTENT(IN)  ::TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,MXNODS
  INTEGER, INTENT(IN) ::NTIME1
  REAL, INTENT(INOUT),DIMENSION(NONODS) ::X,Y,Z
  REAL, INTENT(INOUT),DIMENSION(NONODS) :: MUPTXX,MUPTXY,MUPTXZ
  REAL, INTENT(INOUT),DIMENSION(NONODS) :: MUPTYY,MUPTYZ,MUPTZZ
  double precision, INTENT(INOUT),DIMENSION(NVAR*MXNODS,NSVD) :: LEFTSVD
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS) ::SMEAN
  REAL, INTENT(INOUT),DIMENSION(0:ntime1,NSVD_TOTAL) :: COEF
  REAL, INTENT(IN)    ::DT1
  LOGICAL, INTENT(IN) ::D3,DCYL
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: N,NLX,NLY,NLZ
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: M,MLX,MLY,MLZ
  REAL, INTENT(IN),DIMENSION(NGI)  ::WEIGHT
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC)::PNDGLN
  ! coriolos force
  INTEGER, INTENT(IN)  :: GEOBAL
  REAL, INTENT(IN)  :: SCFACTH0
  ! Wind stress
  !------------
  LOGICAL IFWIND

  ! GEOGRAPHIC PRESSURE
  !--------------------

  INTEGER MPGLOC,PGNODS
  REAL, INTENT(INOUT),DIMENSION(MPGLOC,NGI)  :: MPG,MPGLX,MPGLY,MPGLZ
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MPGLOC):: PGNDGLNO
!  REAL, INTENT(INOUT),DIMENSION(PGNODS,NSVD_PHI) :: LEFTSVD_PG
  REAL  ::SMEAN_PG(PGNODS)
!  REAL, INTENT(INOUT),DIMENSION(0:ntime1,NSVD_PHI) :: COEF_PG
  REAL,DIMENSION(:,:),ALLOCATABLE   ::MPGX,MPGY,MPGZ
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_PGx,WEI_PGy,WEI_PGz
  REAL L1(NGI),L2(NGI),L3(NGI),L4(NGI),PGWEIG(NGI) 
  REAL A11(NGI),A12(NGI),A13(NGI)
  REAL A21(NGI),A22(NGI),A23(NGI)
  REAL A31(NGI),A32(NGI),A33(NGI) 
!  REAL,DIMENSION(:,:),ALLOCATABLE   ::LEFTSVD_PG_x,LEFTSVD_PG_y
  REAL   ::LEFTSVD_PG_x(PGNODS,NSVD_PHI),LEFTSVD_PG_y(PGNODS,NSVD_PHI)
  REAL,DIMENSION(:),ALLOCATABLE   ::PGMATRIX_x,PGMATRIX_y,PGMATRIX_mean
  REAL,DIMENSION(:),ALLOCATABLE            ::PGVEC_x,PGVEC_y,PGVEC_mean
!  integer, dimension(:), allocatable:: PGFINDRM, PGCENTRM, PGCOLM
       INTEGER NPGCOLM
       INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
       INTEGER PGCENTRM(PGNODS) 
 
  !local for geographic pressure
  REAL PG_x_mean,PG_y_mean,PG_z_mean
  REAL PG_x,PG_y,PG_z
  REAL PG1_x,PG1_y,PG1_z
  REAL PG2_x,PG2_y,PG2_z
  REAL Mx_PGx_mean,My_PGy_mean,Mz_PGz_mean
  REAL Mx_PGx,My_PGy,Mz_PGz
  REAL Mx_FV_mean,My_FU_mean,Mx_FV,My_FU
  INTEGER IGLPG,JGLPG

  !       Local variables
  REAL    ::THETA,ONEMINTH
  INTEGER ::ISVD,JSVD,II,JJ,KK,NSTEP1,ELE,GI,ILOC,JLOC,I,J,K
  INTEGER ::IGL,JGL,IGLP,COUNT,COUNT1,JGLP
  REAL    ::RNN,RN1,RN2,RN3,FU,FV,U_x,V_y,W_z,UU_x,VU_y,WU_z,                        &
            UV_x,VV_y,WV_z,UW_x,VW_y,WW_z,P_x,P_y,P_z,                               &
            UD_Ux,VD_Uy,WD_Uz,UD_Vx,VD_Vy,WD_Vz,UD_Wx,VD_Wy,WD_Wz,                   &
            UDmean_Ux,VDmean_Uy,WDmean_Uz,UDmean_Vx,VDmean_Vy,WDmean_Vz,             &
            UDmean_Wx,VDmean_Wy,WDmean_Wz,                                           &
            Umeanx_UD,Umeany_VD,Umeanz_WD,Vmeanx_UD,Vmeany_VD,Vmeanz_WD,             &
            Wmeanx_UD,Wmeany_VD,Wmeanz_WD,                                           &
            U_xx,U_yy,U_zz,V_xx,V_yy,V_zz,W_xx,W_yy,W_zz
  REAL    ::RN1_mean,RN2_mean,RN3_mean,FU_mean,FV_mean,U_x_mean,V_y_mean,W_z_mean,         &
            UU_x_mean,VU_y_mean,WU_z_mean,UV_x_mean,VV_y_mean,WV_z_mean,                   &
            UW_x_mean,VW_y_mean,WW_z_mean,P_x_mean,P_y_mean,P_z_mean,                      &
            U_xx_mean,U_yy_mean,U_zz_mean,V_xx_mean,V_yy_mean,V_zz_mean,W_xx_mean,W_yy_mean,W_zz_mean
  INTEGER ::NCOLM,NCOLM1,NCOLM2
  INTEGER ::IBL11,IBL12,IBL13,IBL14,IBL15
  INTEGER ::IBL21,IBL22,IBL23,IBL24,IBL25
  INTEGER ::IBL31,IBL32,IBL33,IBL34,IBL35
  INTEGER ::IBL41,IBL42,IBL43,IBL44,IBL45
  INTEGER ::IBL51,IBL52,IBL53,IBL54,IBL55
  real    ::VOL,ftest,VOLUME
  REAL    ::DT
  INTEGER ::NSTEP,NTIME,NT,ROW,ROWCNT
  REAL    ::INFIN1,INFIN2,NEAR
      PARAMETER(INFIN1=10.0E+25)
      PARAMETER(INFIN2=10.0E+25)
      PARAMETER(NEAR=1./INFIN1)


  !       Working arrays/local variables            
  !       ********* DEFINE OPTIONS FOR GMRES **************
  REAL     ::ERROR1,ERROR2,RELAX
  INTEGER  ::NOITS1,NOITS2,LEFTP1,LEFTP2,IMM
  INTEGER  ::MINTE1,MINTE2,NR2,IP,MISOUT
  !       ** options for MULMAT 
  LOGICAL  ::SYM,IGUESS
  INTEGER  ::MULPA,UZAWA,CGSOLQ,GMRESQ,CGSQ
  !???? wrong
  REAL     ::DM1(NSVD),DM2(NSVD),R(7*NSVD)
  real     ::rub(2)
  integer  ::iidum(2),kits
  INTEGER  ::MATSTR,NR,NNV
  PARAMETER(MATSTR=0)
  REAL     ::RRR
  REAL ::PREERR
  INTEGER, PARAMETER::halotagT10=3
  INTEGER  ::PARA,halo_tag, halo_tag_p
  PARAMETER(PARA=0,halo_tag=1, halo_tag_p=2)

  INTEGER IMATST,TIMM,TMISOU
  REAL TRELAX
  LOGICAL MOMSYM
  REAL TEROR1
  INTEGER TNOIT1,TNOIT3D1
  !       INTEGER VISUAL
  REAL RD(1)
  INTEGER LOWER,UPPER,INUM
         REAL, ALLOCATABLE, DIMENSION(:)::WORK1
         REAL, ALLOCATABLE, DIMENSION(:)::WORK2
         REAL, ALLOCATABLE, DIMENSION(:)::WORK3
         REAL, ALLOCATABLE, DIMENSION(:)::WORK4
         REAL, ALLOCATABLE, DIMENSION(:)::WORK5
         REAL, ALLOCATABLE, DIMENSION(:)::WORKP


! all matraces for the calculations:
  !???? wrong
  INTEGER  FINDRM(NSVD+1),COLM(NSVD*NSVD), CENTRM(NSVD)

  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM,VEC,VEC1
  REAL,DIMENSION(:,:),ALLOCATABLE   :: ReMat,ReMat1
  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM0,BMat,BMat1,BMat0,BMat01,VEC0
  REAL,DIMENSION(:),ALLOCATABLE     :: DETWEI,XD,YD,ZD
  REAL,DIMENSION(:),ALLOCATABLE     :: UD_mean,VD_mean,WD_mean
  REAL,DIMENSION(:,:),ALLOCATABLE   :: UD,VD,WD
  REAL,DIMENSION(:),ALLOCATABLE     :: Umeanx,Umeany,Umeanz,Vmeanx,Vmeany,Vmeanz,Wmeanx,Wmeany,Wmeanz
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_u,WEI_v,WEI_w,WEI_p
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_ux,WEI_uy,WEI_uz,WEI_vx,WEI_vy,WEI_vz,WEI_wx,WEI_wy,WEI_wz
  REAL,DIMENSION(:,:),ALLOCATABLE   :: NX,NY,NZ,MX,MY,MZ
  REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF,NEWCOEF1,NEWCOEF2
  REAL,DIMENSION(:),ALLOCATABLE     :: CORIOLIS

  REAL,DIMENSION(:,:),ALLOCATABLE   :: BIGM1
  REAL,DIMENSION(:),ALLOCATABLE     :: UDMat
!  REAL,DIMENSION(:),ALLOCATABLE     :: FINDRM,COLM,CENTRM
  REAL,DIMENSION(:),ALLOCATABLE     :: ML

  REAL,DIMENSION(:),ALLOCATABLE     :: SOURCX,SOURCY,SOURCZ
  INTEGER IOS
  LOGICAL IFWRITE 
  INTEGER GLOITS,LINITS
  CHARACTER*240 NAME,NAME1,NAME2

  ewrite(3,*) 'start to run the reduced model_pg'
  NT=1
  DT=DT1/NT
  NTIME=NT*NTIME1

  ALLOCATE(BIGM(NSVD_TOTAL*NSVD_TOTAL) )
  ALLOCATE(VEC(NSVD_TOTAL))
  ALLOCATE(ReMat( NSVD_TOTAL-NSVD,NSVD_TOTAL-NSVD))
  ALLOCATE(BIGM0( (NSVD_TOTAL)*(NSVD_TOTAL) ))
  ALLOCATE(BIGM1((NSVD_TOTAL)*(NSVD_TOTAL),NSVD_u))
  ALLOCATE(VEC0(NSVD_TOTAL))
  ALLOCATE(VEC1(NSVD_TOTAL-NSVD))
  ALLOCATE(BMat(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat0( NSVD_TOTAL*NSVD_TOTAL ))
  ALLOCATE(BMat01( NSVD_TOTAL*NSVD_TOTAL ))
  ALLOCATE(BMat1( NSVD_TOTAL*NSVD_TOTAL ))
  ALLOCATE(DETWEI(NGI))
  ALLOCATE(UD(NGI,NSVD))
  ALLOCATE(VD(NGI,NSVD))
  ALLOCATE(WD(NGI,NSVD))
  ALLOCATE(UD_mean(NGI))
  ALLOCATE(VD_mean(NGI))
  ALLOCATE(WD_mean(NGI))
  ALLOCATE(XD(NGI))
  ALLOCATE(YD(NGI))
  ALLOCATE(ZD(NGI))
  ALLOCATE(NX(NLOC,NGI))
  ALLOCATE(NY(NLOC,NGI))
  ALLOCATE(NZ(NLOC,NGI))
  ALLOCATE(MX(NLOC,NGI))
  ALLOCATE(MY(NLOC,NGI))
  ALLOCATE(MZ(NLOC,NGI))
  ALLOCATE(WEI_u(NGI))
  ALLOCATE(WEI_v(NGI))
  ALLOCATE(WEI_w(NGI))
  ALLOCATE(WEI_p(NGI))
  ALLOCATE(WEI_ux(NGI))
  ALLOCATE(WEI_uy(NGI))
  ALLOCATE(WEI_uz(NGI))
  ALLOCATE(WEI_vx(NGI))
  ALLOCATE(WEI_vy(NGI))
  ALLOCATE(WEI_vz(NGI))
  ALLOCATE(WEI_wx(NGI))
  ALLOCATE(WEI_wy(NGI))
  ALLOCATE(WEI_wz(NGI))
  ALLOCATE(NEWCOEF(NSVD_TOTAL))
  ALLOCATE(NEWCOEF1(NSVD_TOTAL-NSVD))
  ALLOCATE(CORIOLIS(NGI))

  ALLOCATE(SOURCX(NONODS))
  ALLOCATE(SOURCY(NONODS))

  ! for geographic pressure
  ALLOCATE(MPGX(MPGLOC,NGI))
  ALLOCATE(MPGY(MPGLOC,NGI))
  ALLOCATE(MPGZ(MPGLOC,NGI))
  ALLOCATE(WEI_PGx(NGI))
  ALLOCATE(WEI_PGy(NGI))
  ALLOCATE(WEI_PGz(NGI))
!  ALLOCATE(LEFTSVD_PG_x(PGNODS,NSVD_PHI))
!  ALLOCATE(LEFTSVD_PG_y(PGNODS,NSVD_PHI))

   ewrite(3,*) 'NPGCOLM,PGNODS',NPGCOLM,PGNODS
  ALLOCATE(PGMATRIX_x(NPGCOLM))
  ALLOCATE(PGMATRIX_y(NPGCOLM))
  ALLOCATE(PGMATRIX_mean(NPGCOLM))
  ALLOCATE(PGVEC_x(PGNODS))
  ALLOCATE(PGVEC_y(PGNODS))
  ALLOCATE(PGVEC_mean(PGNODS))


!**************************************************!
!       Initialise all matrices                    !
!**************************************************!
  BIGM =0.0
  VEC =0.0
  ReMat=0.0
  BIGM0=0.0
  BIGM1=0.0
  VEC0=0.0
  VEC1=0.0
  BMat=0.0
  BMat0=0.0
  BMat01=0.0
  BMat1=0.0
  DETWEI=0.0
  UD=0.0
  VD=0.0
  WD=0.0
  UD_mean=0.0
  VD_mean=0.0
  WD_mean=0.0
  XD=0.0
  YD=0.0
  ZD=0.0
  WEI_u=0.0
  WEI_v=0.0
  WEI_w=0.0
  WEI_p=0.0
  WEI_ux=0.0
  WEI_uy=0.0
  WEI_uz=0.0
  WEI_vx=0.0
  WEI_vy=0.0
  WEI_vz=0.0
  WEI_wx=0.0
  WEI_wy=0.0
  WEI_wz=0.0
  NEWCOEF=0.0
  NEWCOEF1=0.0
  CORIOLIS=0.0


  NX=0.0
  NY=0.0
  NZ=0.0
  MX=0.0
  MY=0.0
  MZ=0.0

  SOURCX=0.0
  SOURCY=0.0
           if(IFWIND) then
           ! WIND STRESS
           !-------------
           CALL WINDYSOURCE(SOURCX,SOURCY,X,Y,NONODS,NSVD)
          endif

!******************************************************!
! Define the pointer for storing data in the matrices  !
!******************************************************!
  NCOLM=NSVD*NSVD
  NCOLM1=NSVD_U*NSVD_PHI
  IBL11=0
  IBL12=1*NCOLM
  IBL13=2*NCOLM
  IBL14=3*NCOLM

  IBL21=3*NCOLM+NCOLM1
  IBL22=IBL21+NCOLM
  IBL23=IBL21+2*NCOLM
  IBL24=IBL21+3*NCOLM

  IBL31=IBL24+NCOLM1
  IBL32=IBL31+NCOLM
  IBL33=IBL31+2*NCOLM
  IBL34=IBL31+3*NCOLM

  IBL41=IBL34+NCOLM1
  IBL42=IBL41+NCOLM1
  IBL43=IBL41+2*NCOLM1
  IBL44=IBL41+3*NCOLM1

  R(1:7*NSVD)=0.0

! give the initial COEF..

! for 2D
  DO II=2*NSVD_U+1,3*NSVD_U
     COEF(0,II)=0.0
  ENDDO

  open(1,file='leftsvd0c.dat')
  do jj=1,NVAR*MXNODS
  write(1,*) (leftsvd(jj,ii),ii=1,nsvd)
  enddo
  close(1)
!     LEFTSVD(2*nonods+1:3*nonods,:)=0.0
     SMEAN(2*mxnods+1:3*mxnods)=0.0

    NEWCOEF(1:NSVD_TOTAL)=COEF(0,:)

  THETA=0.5
  ONEMINTH=1.0-THETA


  ! for geostraphic pressure

  MPGLX=0.0
  MPGLY=0.0
  MPGLZ=0.0
  MPGX=0.0
  MPGY=0.0
  MPGZ=0.0
  WEI_PGx=0.0
  WEI_PGy=0.0
  WEI_PGz=0.0
  LEFTSVD_PG_x=0.0
  LEFTSVD_PG_y=0.0
  SMEAN_PG=0.0


      MPGLOC=10
         CALL TRIQUA(L1, L2, L3, L4, PGWEIG, D3,NGI)
 ! Work out the shape functions and there derivatives...
        CALL SHATRI(L1, L2, L3, L4, PGWEIG, D3, &
                 MPGLOC,NGI,                     &     
                 MPG,MPGLX,MPGLY,MPGLZ) 
!  ! Calculate the matrix sparcity...
!  allocate(pgfindrm(1:pgnods+1))
!  allocate(pgcentrm(1:pgnods))
!  CALL POSINM(TOTELE,PGNODS,MPGLOC, &
!       R(RPT),NPGCOLM,NRMEM-RPT+1, &
!       PGFINDRM, PGCENTRM, &
!       PGNODS,MPGLOC,PGNDGLNO, &
!       PGNDGLNO)
!   allocate(pgcolm(NPGCOLM))
!   call icopy(pgcolm,R(RPT),NPGCOLM)
 ! 
! Calculate the basis function for geostrophic pressure
! Solve for geostrophic pressure options... ************

IF(IFWRITE.AND.LINITS.EQ.1) then

  DO 300 ISVD=1,NSVD_PHI
       PGMATRIX_x=0.0
       PGMATRIX_y=0.0
       PGMATRIX_mean=0.0
       PGVEC_x=0.0
       PGVEC_y=0.0
       PGVEC_mean=0.0

     DO 200 ELE=1,TOTELE
         CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,MLOC,NGI,      &    
                M,MLX,MLY,MLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
                MX,MY,MZ) 

         CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MPGLOC,NGI, &
                    N,NLX,NLY,NLZ, MPG,MPGLX,MPGLY,MPGLZ,WEIGHT, DETWEI,.TRUE.,.FALSE., &
                    NX,NY,NZ,MPGX,MPGY,MPGZ,                   &
                    A11,A12,A13, A21,A22,A23, A31,A32,A33,     &
                    XD,YD,ZD,                                  &
                    2) 

           DO GI=1,NGI
              CORIOLIS(GI) = 2.0*FUNOME(XD(GI),YD(GI),ZD(GI))
           ENDDO

        DO 160 ILOC=1,MPGLOC
         IGLPG=PGNDGLNO((ELE-1)*MPGLOC+ILOC)
          DO 170 JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
             IF(ISVD.EQ.1) THEN         
              Mx_FV_mean=0.0
              My_FU_mean=0.0
              DO GI=1,NGI
                 My_FU_mean =My_FU_mean+MPGY(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*SMEAN(JGL)*N(JLOC,GI)
                 Mx_FV_mean =Mx_FV_mean-MPGX(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*N(JLOC,GI)
              ENDDO
              PGVEC_mean(IGLPG)=PGVEC_mean(IGLPG)-(Mx_FV_mean+My_FU_mean)
             ENDIF
              Mx_FV=0.0
              My_FU=0.0
              DO GI=1,NGI
                 My_FU = My_FU+MPGY(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*LEFTSVD(JGL,ISVD)*N(JLOC,GI)
                 Mx_FV = Mx_FV-MPGX(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,ISVD)*N(JLOC,GI)
              ENDDO
              PGVEC_x(IGLPG)=PGVEC_x(IGLPG)-Mx_FV
              PGVEC_y(IGLPG)=PGVEC_y(IGLPG)-My_FU
170          CONTINUE

          DO 190 JLOC=1,MPGLOC
              JGLPG=PGNDGLNO((ELE-1)*MPGLOC+JLOC)
!              ewrite(3,*) 'JGLPG=',JGLPG
!               Find count. ***************************************
                               LOWER=PGFINDRM(IGLPG) 
                               UPPER=PGFINDRM(IGLPG+1)-1
               7000                 CONTINUE
                               INUM=LOWER+(UPPER-LOWER+1)/2 
                               IF(JGLPG.GE.PGCOLM(INUM) )  THEN 
                                  LOWER=INUM
                               ELSE
                                  UPPER=INUM
                               ENDIF
                               IF(UPPER-LOWER.LE.1) THEN
                                  IF(JGLPG.EQ.PGCOLM(LOWER)) THEN
                                     COUNT=LOWER
                                  ELSE
                                     COUNT=UPPER
                                  ENDIF
                                  GOTO 9000
                               ENDIF
                               GOTO 7000
               9000                 CONTINUE

              Mx_PGx_mean=0.0
              My_PGy_mean=0.0
              Mz_PGz_mean=0.0
             IF(ISVD.EQ.1) THEN         
              DO GI=1,NGI
                 Mx_PGx_mean=Mx_PGx_mean+MPGX(ILOC,GI)*DETWEI(GI)*MPGX(JLOC,GI)
                 My_PGy_mean=My_PGy_mean+MPGY(ILOC,GI)*DETWEI(GI)*MPGY(JLOC,GI)
                 Mz_PGz_mean=Mz_PGz_mean+MPGZ(ILOC,GI)*DETWEI(GI)*MPGZ(JLOC,GI)
              ENDDO
              PGMATRIX_mean(COUNT)=PGMATRIX_mean(COUNT)+(Mx_PGx_mean+My_PGy_mean+Mz_PGz_mean)
             ENDIF
                 Mx_PGx=0.0
                 My_PGy=0.0
                 Mz_PGz=0.0
                 DO GI=1,NGI
                    Mx_PGx=Mx_PGx+MPGX(ILOC,GI)*DETWEI(GI)*MPGX(JLOC,GI)
                    My_PGy=My_PGy+MPGY(ILOC,GI)*DETWEI(GI)*MPGY(JLOC,GI)
                    Mz_PGz=Mz_PGz+MPGZ(ILOC,GI)*DETWEI(GI)*MPGZ(JLOC,GI)
                 ENDDO
                   PGMATRIX_x(COUNT)=PGMATRIX_x(COUNT) +Mx_PGx+My_PGy+Mz_PGz
                   PGMATRIX_y(COUNT)=PGMATRIX_y(COUNT) +Mx_PGx+My_PGy+Mz_PGz

190              CONTINUE
160              CONTINUE


200              CONTINUE

 
   ALLOCATE(WORK1(PGNODS))
   ALLOCATE(WORK2(PGNODS))
   ALLOCATE(WORK3(PGNODS))
   ALLOCATE(WORK4(PGNODS))
   ALLOCATE(WORK5(PGNODS)) 
   ALLOCATE(WORKP(4*PGNODS))
   WORK1=0.0
   WORK2=0.0
   WORK3=0.0
   WORK4=0.0
   WORK5=0.0
   WORKP=0.0

  IMATST=0
!  TIMM=-1
!  TMISOU=1
  TIMM=1
  TMISOU=1
  TRELAX=1.0
  MOMSYM=.FALSE.
!  TEROR1=PREERR/100000.
   PREERR=1.E-10
  TEROR1=PREERR/10.
!  TNOIT1=PRENOI*20
  TNOIT1=3000
  IIDUM=0
  TNOIT3D1=TNOIT1

   ROW = 1
   ROWCNT = PGCENTRM(ROW)
   PGMATRIX_x(ROWCNT) = INFIN1
   PGMATRIX_y(ROWCNT) = INFIN1
   PGMATRIX_mean(ROWCNT) =  INFIN1

  CALL SOLCG(LEFTSVD_PG_x,PGVEC_x,PGNODS,PGNODS,PGNODS,.TRUE.,  &
       PGMATRIX_x,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM,&
       halotagT10,KITS)
 

  CALL SOLCG(LEFTSVD_PG_y,PGVEC_y,PGNODS,PGNODS,PGNODS,.TRUE.,  &
       PGMATRIX_y,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM,&
       halotagT10, KITS)



2000 format (20(1X,E24.16)) 


IF(ISVD.EQ.1) THEN         
  CALL SOLCG(SMEAN_PG,PGVEC_mean,PGNODS,PGNODS,PGNODS,.TRUE.,  &
       PGMATRIX_mean,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM,&
       halotagT10,&
       KITS)
ENDIF

    DEALLOCATE(WORK1)
    DEALLOCATE(WORK2)
    DEALLOCATE(WORK3)
    DEALLOCATE(WORK4)
    DEALLOCATE(WORK5)
    DEALLOCATE(WORKP)

300              CONTINUE
ENDIF !end (IFWRITE.AND.LINITS.EQ.1) if re-calculating the geostrophic pressure or passing from forward model

      IF(IFWRITE.AND.LINITS.EQ.1) THEN
        write(name1,'(a,i0)') 'LEFTSVD_PG_x.dat.',GLOITS
        write(name2,'(a,i0)') 'LEFTSVD_PG_y.dat.',GLOITS
        OPEN(1,FILE=trim(NAME1),STATUS='UNKNOWN',IOSTAT=IOS)
        OPEN(2,FILE=trim(NAME2),STATUS='UNKNOWN',IOSTAT=IOS)
              DO j=1,NSVD
                write(1,*) (LEFTSVD_PG_x(I,j),I=1,PGNODS)
                write(2,*) (LEFTSVD_PG_y(I,j),I=1,PGNODS)
              ENDDO
        close(1)
        close(2)

        write(name,'(a,i0)') 'SMEAN_PG.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (SMEAN_PG(i),I=1,PGNODS)
        close(1)
      ELSE
        write(name1,'(a,i0)') 'LEFTSVD_PG_x.dat.',GLOITS
        write(name2,'(a,i0)') 'LEFTSVD_PG_y.dat.',GLOITS
        OPEN(1,FILE=trim(NAME1),STATUS='UNKNOWN',IOSTAT=IOS)
        OPEN(2,FILE=trim(NAME2),STATUS='UNKNOWN',IOSTAT=IOS)
              DO j=1,NSVD
                read(1,*) (LEFTSVD_PG_x(I,j),I=1,PGNODS)
                read(2,*) (LEFTSVD_PG_y(I,j),I=1,PGNODS)
              ENDDO
        close(1)
        close(2)

        write(name,'(a,i0)') 'SMEAN_PG.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (SMEAN_PG(i),I=1,PGNODS)
        close(1)

      ENDIF


  ewrite(3,*) 'nsvd,nvar,nsvd',  nsvd,nvar,nsvd

  ! clean the matrix BIM and left vector VEC
  !-----------------------------------------

  ! define the row and colomn
  !---------------------------
!  DO II=1,NSVD
!     CENTRM(II)=(II-1)*NSVD+II
!     DO JJ=1,NSVD
!        COLM( (II-1)*NSVD+JJ )=JJ
!     ENDDO
!  ENDDO
!
!  DO II=1,NSVD+1
!     FINDRM(II)=(II-1)*NSVD+1
!  ENDDO

! go to 400
IF(IFWRITE.AND.LINITS.EQ.1) then
!  DO 100 NSTEP=1,NTIME
     ewrite(3,*) '11'
     ! setup the reduced model for BETA_ii(t)
     ! the original momentum equations are multiplied by the POD basic vectors 
     ! PHI_ii^u(x),PHI_ii^v(x),PHI_ii^w(x), and continuity equation is multiplied by PHI_ii^p(x)
     ! and then integrate them over the whole domain
     !  PHI_ii^p*{U_x+V_y+W_z}
     ! =PHI_ii{[dU_mean(X)/dx+COEF(t)*(d/dx)*SUM(PHI_ii^u)]+....}
     !------------------------------------------------------------------------------------------
     !
     DO 30 ELE=1,TOTELE
           !  ewrite(3,*)  TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC
           !  ewrite(3,*)  D3,DCYL
           ! Calculate DETWEI,NX,NY,NZ for element ELE
!           CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,      &    
!                N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
!                NX,NY,NZ) 

           ! Calculate DETWEI,MX,MY,MZ for element ELE (for pressure)
           CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,MLOC,NGI,      &    
                M,MLX,MLY,MLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
                MX,MY,MZ) 

         CALL CALHIGHNODAI(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,MPGLOC,NGI, &
                    N,NLX,NLY,NLZ, MPG,MPGLX,MPGLY,MPGLZ,WEIGHT, DETWEI,.TRUE.,.FALSE., &
                    NX,NY,NZ,MPGX,MPGY,MPGZ,                   &
                    A11,A12,A13, A21,A22,A23, A31,A32,A33,     &
                    XD,YD,ZD,                                  &
                    2) 

!                    ewrite(3,*) 'MPGLX',MPGLX
!                    ewrite(3,*) 'MPGLY',MPGLY
!                    ewrite(3,*) 'MPGLZ',MPGLZ
!                    ewrite(3,*) 'MPGX',MPGX
!                    ewrite(3,*) 'MPGY',MPGY
!                    ewrite(3,*) 'MPGZ',MPGZ
!                    stop

!     ewrite(3,*) 'ELE=',ELE

           DO 20 ISVD=1,NSVD

           ! CALCULATE THE WEIGHTING FUNCTION(In POD, that is the POD basic function) at the Gaussian points
           ! e.g., in FEM. 
           !WEI_u =[N(GI)*PHI(X)_ii^u], 
           !here , N(GI)=the basic function, PHI(X)_ii^u is the POD basic vectors
           !WEI_ux=[NX(GI)*PHI(X)_ii^u],
           !WEI_uy=[NY(GI)*PHI(X)_ii^u],
           !WEI_uz=[NZ(GI)*PHI(X)_ii^u],
           !WEI_v =[N(GI)*PHI(X)_ii^v], 
           !WEI_vx=[NX(GI)*PHI(X)_ii^v],
           !WEI_vy=[NY(GI)*PHI(X)_ii^v],
           !WEI_vz=[NZ(GI)*PHI(X)_ii^v],
           !WEI_w =[N(GI)*PHI(X)_ii^w], 
           !WEI_wx=[NX(GI)*PHI(X)_ii^w],
           !WEI_wy=[NY(GI)*PHI(X)_ii^w],
           !WEI_wz=[NZ(GI)*PHI(X)_ii^w],
           !-------------------------------
          DO GI=1,NGI
              WEI_u(GI)=0.0
              WEI_v(GI)=0.0
              WEI_w(GI)=0.0
              WEI_ux(GI)=0.0
              WEI_uy(GI)=0.0
              WEI_uz(GI)=0.0
              WEI_vx(GI)=0.0
              WEI_vy(GI)=0.0
              WEI_vz(GI)=0.0
              WEI_wx(GI)=0.0
              WEI_wy(GI)=0.0
              WEI_wz(GI)=0.0

              DO ILOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+ILOC)
                 
                 WEI_u(GI)=WEI_u(GI)+LEFTSVD(JGL,ISVD)*N(ILOC,GI)
                 WEI_v(GI)=WEI_v(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*N(ILOC,GI)
                 WEI_w(GI)=WEI_w(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*N(ILOC,GI)
                 WEI_ux(GI)=WEI_ux(GI)+LEFTSVD(JGL,ISVD)*NX(ILOC,GI)
                 WEI_uy(GI)=WEI_uy(GI)+LEFTSVD(JGL,ISVD)*NY(ILOC,GI)
                 WEI_uz(GI)=WEI_uz(GI)+LEFTSVD(JGL,ISVD)*NZ(ILOC,GI)
                 WEI_vx(GI)=WEI_vx(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*NX(ILOC,GI)
                 WEI_vy(GI)=WEI_vy(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*NY(ILOC,GI)
                 WEI_vz(GI)=WEI_vz(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*NZ(ILOC,GI)
                 WEI_wx(GI)=WEI_wx(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*NX(ILOC,GI)
                 WEI_wy(GI)=WEI_wy(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*NY(ILOC,GI)
                 WEI_wz(GI)=WEI_wz(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*NZ(ILOC,GI)
              ENDDO
           ENDDO


           ! Weight for pressure
           !WEI_p(GI) =N(GI)*PHI(X)_ii^p
           !---------------
!           IF(ISVD.LE.NSVD_PHI) THEN
              DO GI=1,NGI
                 WEI_p(GI)=0.0
                 DO ILOC=1,MLOC
                    IGLP=PNDGLN((ELE-1)*MLOC+ILOC)
                    WEI_p(GI)=WEI_p(GI)+LEFTSVD(3*MXNODS+IGLP,ISVD)*M(ILOC,GI)
                 ENDDO
              ENDDO

!           ENDIF


           !PHI_jj*U=PHI_jj*[U_mean+SUM_jj=1^nsvd(COEF*PHI_ii)]=PHI_jj*U_mean+COEF(t,jj)
           ! UD(GI)=COEF(t,jj)+phi_jj*[N(GI)*U(X)_mean]  
           ! VD(GI)=COEF(t,jj)+phi_jj*[N(GI)*V(X)_mean]  
           ! WD(GI)=COEF(t,jj)+phi_jj*[N(GI)*W(X)_mean]  
           !-----------------------------------------------
           DO GI=1,NGI
              UD(GI,:)=0.0
              VD(GI,:)=0.0
              WD(GI,:)=0.0
              UD_mean(GI)=0.0
              VD_mean(GI)=0.0
              WD_mean(GI)=0.0
!              XD(GI)=0.0
!              YD(GI)=0.0
!              ZD(GI)=0.0
              CORIOLIS(GI)=0.0
              

              DO JLOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+JLOC)

                 UD_mean(GI)=UD_mean(GI)+SMEAN(JGL)*N(JLOC,GI)
                 VD_mean(GI)=VD_mean(GI)+SMEAN(1*MXNODS+JGL)*N(JLOC,GI)
                 WD_mean(GI)=WD_mean(GI)+SMEAN(2*MXNODS+JGL)*N(JLOC,GI)
!                 XD(GI)=XD(GI) + N(JLOC,GI)*X(JGL)
!                 YD(GI)=YD(GI) + N(JLOC,GI)*Y(JGL)
!                 ZD(GI)=ZD(GI) + N(JLOC,GI)*Z(JGL)

                  DO JSVD = 1,NSVD
                     UD(GI,JSVD)=UD(GI,JSVD)+LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                     VD(GI,JSVD)=VD(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
                     WD(GI,JSVD)=WD(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
                  ENDDO

              ENDDO
              CORIOLIS(GI) = 2.0*FUNOME(XD(GI),YD(GI),ZD(GI))
              ! testing
!              CORIOLIS(GI) =0.0
           ENDDO


           !Calculate
           !------
           !PHI_jj^p*[U_x_mean+V_y_mean+W_z_mean]
           !------------------
           !PHI_jj^u*[U_x_mean+UD*U_x_mean+VD*U_y_mean+WD*U_z_mean-FV_mean+P_x_mean]
           !+D(PHI_jj^u)/DX*DU_mean/DX+D(PHI_jj^u)/DY*DU_mean/DY+D(PHI_jj^u)/DZ*DU_mean/DZ
           !-----------------
           !PHI_jj^v*[V_x_mean+UD*V_x_mean+VD*V_y_mean+WD*V_z_mean+FU_mean+P_y_mean]
           !+D(PHI_jj^v)/DX*DV_mean/DX+D(PHI_jj^v)/DY*DV_mean/DY+D(PHI_jj^v)/DZ*DV_mean/DZ
           !-----------------
           !PHI_jj^w*[W_x_mean+UD*W_x_mean+VD*W_y_mean+WD*W_z_mean+P_z_mean]
           !+D(PHI_jj^w)/DX*DW_mean/DX+D(PHI_jj^w)/DY*DW_mean/DY+D(PHI_jj^w)/DZ*DW_mean/DZ
           !---------------------------------

              FU_mean =0.0
              FV_mean =0.0
              U_x_mean =0.0
              V_y_mean =0.0
              W_Z_mean =0.0
              UU_x_mean=0.0
              VU_y_mean=0.0
              WU_z_mean=0.0
              UV_x_mean=0.0
              VV_y_mean=0.0
              WV_z_mean=0.0
              UW_x_mean=0.0
              VW_y_mean=0.0
              WW_z_mean=0.0
              U_xx_mean=0.0
              U_yy_mean=0.0
              U_zz_mean=0.0
              V_xx_mean=0.0
              V_yy_mean=0.0
              V_zz_mean=0.0
              W_xx_mean=0.0
              W_yy_mean=0.0
              W_zz_mean=0.0
              P_x_mean =0.0
              P_y_mean =0.0
              P_z_mean =0.0
              PG_x_mean =0.0
              PG_y_mean =0.0
              PG_z_mean =0.0

           DO JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
              DO GI=1,NGI
                 
                 FU_mean =FU_mean+WEI_v(GI)*CORIOLIS(GI)*DETWEI(GI)*SMEAN(JGL)*N(JLOC,GI)
                 FV_mean =FV_mean-WEI_u(GI)*CORIOLIS(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*N(JLOC,GI)
                 U_x_mean =U_x_mean+WEI_p(GI)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                 V_y_mean =V_y_mean+WEI_p(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 W_Z_mean =W_Z_mean+WEI_p(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)
                 UU_x_mean=UU_x_mean+WEI_u(GI)*UD_mean(GI)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                 VU_y_mean=VU_y_mean+WEI_u(GI)*VD_mean(GI)*DETWEI(GI)*SMEAN(JGL)*NY(JLOC,GI)
                 WU_z_mean=WU_z_mean+WEI_u(GI)*WD_mean(GI)*DETWEI(GI)*SMEAN(JGL)*NZ(JLOC,GI)
                 UV_x_mean=UV_x_mean+WEI_v(GI)*UD_mean(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                 VV_y_mean=VV_y_mean+WEI_v(GI)*VD_mean(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 WV_z_mean=WV_z_mean+WEI_v(GI)*WD_mean(GI)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)
                 UW_x_mean=UW_x_mean+WEI_w(GI)*UD_mean(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                 VW_y_mean=VW_y_mean+WEI_w(GI)*VD_mean(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                 WW_z_mean=WW_z_mean+WEI_w(GI)*WD_mean(GI)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)

                 U_xx_mean=U_xx_mean+WEI_ux(GI)*MUPTXX(JGL)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                 U_yy_mean=U_yy_mean+WEI_uy(GI)*MUPTYY(JGL)*DETWEI(GI)*SMEAN(JGL)*NY(JLOC,GI)
                 U_zz_mean=U_zz_mean+WEI_uz(GI)*MUPTZZ(JGL)*DETWEI(GI)*SMEAN(JGL)*NZ(JLOC,GI)
                 V_xx_mean=V_xx_mean+WEI_vx(GI)*MUPTXX(JGL)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                 V_yy_mean=V_yy_mean+WEI_vy(GI)*MUPTYY(JGL)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 V_zz_mean=V_zz_mean+WEI_vz(GI)*MUPTZZ(JGL)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)
                 W_xx_mean=W_xx_mean+WEI_wx(GI)*MUPTXX(JGL)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                 W_yy_mean=W_yy_mean+WEI_wy(GI)*MUPTYY(JGL)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                 W_zz_mean=W_zz_mean+WEI_wz(GI)*MUPTZZ(JGL)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)
!                 U_zz_mean=0.0
!                 V_zz_mean=0.0
!                 W_zz_mean=0.0


              ENDDO
           ENDDO

              ! Calculate the average pressure p_x_mean, P_y_mean, p_z_mean

          DO JLOC=1,MLOC
              JGLP=PNDGLN((ELE-1)*MLOC+JLOC)
              DO GI=1,NGI
                 P_x_mean = P_x_mean+WEI_u(GI)*DETWEI(GI)*SMEAN(3*MXNODS+JGLP)*MX(JLOC,GI)
                 P_y_mean = P_y_mean+WEI_v(GI)*DETWEI(GI)*SMEAN(3*MXNODS+JGLP)*MY(JLOC,GI)
                 P_z_mean = P_z_mean+WEI_w(GI)*DETWEI(GI)*SMEAN(3*MXNODS+JGLP)*MZ(JLOC,GI)
              ENDDO
           ENDDO

           DO JLOC=1,MPGLOC
              JGLPG=PGNDGLNO((ELE-1)*MPGLOC+JLOC)
              DO GI=1,NGI
                 PG_x_mean = PG_x_mean+WEI_u(GI)*DETWEI(GI)*SMEAN_PG(JGLPG)*MPGX(JLOC,GI)
                 PG_y_mean = PG_y_mean+WEI_v(GI)*DETWEI(GI)*SMEAN_PG(JGLPG)*MPGY(JLOC,GI)
                 PG_z_mean = PG_z_mean+WEI_w(GI)*DETWEI(GI)*SMEAN_PG(JGLPG)*MPGZ(JLOC,GI)
              ENDDO
           ENDDO

           VEC0(ISVD)= VEC0(ISVD)-DT*(UU_x_mean+VU_y_mean+WU_z_mean+FV_mean+P_x_mean+PG_x_mean+               &
                U_xx_mean+U_yy_mean+U_zz_mean)
           VEC0(1*NSVD+ISVD)= VEC0(1*NSVD+ISVD)-DT*(UV_x_mean+VV_y_mean+WV_z_mean+FU_mean+P_y_mean+PG_y_mean+  &
                V_xx_mean+V_yy_mean+V_zz_mean)
           VEC0(2*NSVD+ISVD)= VEC0(2*NSVD+ISVD)-DT*(UW_x_mean+VW_y_mean+WW_z_mean+P_z_mean+PG_z_mean+         &
                W_xx_mean+W_yy_mean+W_zz_mean)

           IF(ISVD.LE.NSVD_PHI) THEN
              VEC(3*NSVD+ISVD)= VEC(3*NSVD+ISVD)-U_x_mean-V_y_mean-W_z_mean
! for 2D
!              VEC0(3*NSVD+ISVD)= VEC0(3*NSVD+ISVD)-U_x_mean-V_y_mean

           ENDIF

           if(IFWIND) then
           ! WIND STRESS
           !-------------
!           CALL WINDYSOURCE(SOURCX,SOURCY,X,Y,NONODS,NSVD)
           DO ILOC=1,NLOC
              IGL = NDGLNO((ELE-1)*NLOC + ILOC)
              DO GI=1,NGI
                 VEC0(ISVD) = VEC0(ISVD) +DT*WEI_u(GI)*DETWEI(GI)*SOURCX(IGL)*N(ILOC,GI)
              END DO
           ENDDO
           endif



           !Clculate
           !PHI_jj^u*U; PHI_jj^v*V; PHI_jj^w*W
           !------
           !PHI_jj^p*[U_x+V_y+W_z]
           !------------------
           !PHI_jj^u*[U_x+UD*U_x+VD*U_y+WD*U_z+FV+P_x]
           !+D(PHI_jj^u)/DX*DU/DX+D(PHI_jj^u)/DY*DU/DY+D(PHI_jj^u)/DZ*DU/DZ
           !-----------------
           !PHI_jj^v*[V_x+UD*V_x+VD*V_y+WD*V_z-FU+P_y]
           !+D(PHI_jj^v)/DX*DV/DX+D(PHI_jj^v)/DY*DV/DY+D(PHI_jj^v)/DZ*DV/DZ
           !-----------------
           !PHI_jj^w*[W_x+UD*W_x+VD*W_y+WD*W_z+P_z]
           !+D(PHI_jj^w)/DX*DW/DX+D(PHI_jj^w)/DY*DW/DY+D(PHI_jj^w)/DZ*DW/DZ
           !---------------------------------
           !   ewrite(3,*) 'vec1',VEC(ISVD)

           DO 40 JSVD=1,NSVD

              COUNT=(ISVD-1)*NSVD+JSVD
              COUNT1=(ISVD-1)*NSVD_PHI+JSVD

                 RN1=0.0
                 RN2=0.0
                 RN3=0.0
                 U_x=0.0
                 V_y=0.0
                 W_Z=0.0

                 UDmean_Ux=0.0
                 VDmean_Uy=0.0
                 WDmean_Uz=0.0
                 UDmean_Vx=0.0
                 VDmean_Vy=0.0
                 WDmean_Vz=0.0
                 UDmean_Wx=0.0
                 VDmean_Wy=0.0
                 WDmean_Wz=0.0

                 UD_Ux=0.0
                 VD_Uy=0.0
                 WD_Uz=0.0
                 UD_Vx=0.0
                 VD_Vy=0.0
                 WD_Vz=0.0
                 UD_Wx=0.0
                 VD_Wy=0.0
                 WD_Wz=0.0

                 Umeanx_UD=0.0
                 Umeany_VD=0.0
                 Umeanz_WD=0.0
                 Vmeanx_UD=0.0
                 Vmeany_VD=0.0
                 Vmeanz_WD=0.0
                 Wmeanx_UD=0.0
                 Wmeany_VD=0.0
                 Wmeanz_WD=0.0

                 U_xx=0.0
                 U_yy=0.0
                 U_zz=0.0
                 V_xx=0.0
                 V_yy=0.0
                 V_zz=0.0
                 W_xx=0.0
                 W_yy=0.0
                 W_zz=0.0
                 FU=0.0
                 FV=0.0
                 P_x=0.0
                 P_y=0.0
                 P_z=0.0
                 PG1_x=0.0
                 PG1_y=0.0
                 PG1_z=0.0
                 PG2_x=0.0
                 PG2_y=0.0
                 PG2_z=0.0

              DO JLOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+JLOC)
                 DO GI=1,NGI
                    RN1=RN1+WEI_u(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    RN2=RN2+WEI_v(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    RN3=RN3+WEI_w(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    IF(ISVD.LE.NSVD_PHI) THEN
                       U_x=U_x+WEI_p(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                       V_y=V_y+WEI_p(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                       W_Z=W_Z+WEI_p(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    ENDIF
                    UDmean_Ux=UDmean_Ux+WEI_u(GI)*UD_mean(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    VDmean_Uy=VDmean_Uy+WEI_u(GI)*VD_mean(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    WDmean_Uz=WDmean_Uz+WEI_u(GI)*WD_mean(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                    UDmean_Vx=UDmean_Vx+WEI_v(GI)*UD_mean(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VDmean_Vy=VDmean_Vy+WEI_v(GI)*VD_mean(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WDmean_Vz=WDmean_Vz+WEI_v(GI)*WD_mean(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    UDmean_Wx=UDmean_Wx+WEI_w(GI)*UD_mean(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VDmean_Wy=VDmean_Wy+WEI_w(GI)*VD_mean(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WDmean_Wz=WDmean_Wz+WEI_w(GI)*WD_mean(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)



                    Umeanx_UD=Umeanx_UD+WEI_u(GI)*UD(GI,JSVD)*DETWEI(GI)*SMEAN(JGL)*NX(JLOC,GI)
                    Umeany_VD=Umeany_VD+WEI_u(GI)*VD(GI,JSVD)*DETWEI(GI)*SMEAN(JGL)*NY(JLOC,GI)
                    Umeanz_WD=Umeanz_WD+WEI_u(GI)*WD(GI,JSVD)*DETWEI(GI)*SMEAN(JGL)*NZ(JLOC,GI)
                    Vmeanx_UD=Vmeanx_UD+WEI_v(GI)*UD(GI,JSVD)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                    Vmeany_VD=Vmeany_VD+WEI_v(GI)*VD(GI,JSVD)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                    Vmeanz_WD=Vmeanz_WD+WEI_v(GI)*WD(GI,JSVD)*DETWEI(GI)*SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)
                    Wmeanx_UD=Wmeanx_UD+WEI_w(GI)*UD(GI,JSVD)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                    Wmeany_VD=Wmeany_VD+WEI_w(GI)*VD(GI,JSVD)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                    Wmeanz_WD=Wmeanz_WD+WEI_w(GI)*WD(GI,JSVD)*DETWEI(GI)*SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)

!
                    UD_Ux=UD_Ux+WEI_u(GI)*UD(GI,JSVD)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    VD_Uy=VD_Uy+WEI_u(GI)*VD(GI,JSVD)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    WD_Uz=WD_Uz+WEI_u(GI)*WD(GI,JSVD)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                    UD_Vx=UD_Vx+WEI_v(GI)*UD(GI,JSVD)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VD_Vy=VD_Vy+WEI_v(GI)*VD(GI,JSVD)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WD_Vz=WD_Vz+WEI_v(GI)*WD(GI,JSVD)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    UD_Wx=UD_Wx+WEI_w(GI)*UD(GI,JSVD)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VD_Wy=VD_Wy+WEI_w(GI)*VD(GI,JSVD)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WD_Wz=WD_Wz+WEI_w(GI)*WD(GI,JSVD)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)


                    U_xx=U_xx+WEI_ux(GI)*MUPTXX(JGL)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    U_yy=U_yy+WEI_uy(GI)*MUPTYY(JGL)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    U_zz=U_zz+WEI_uz(GI)*MUPTZZ(JGL)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                    V_xx=V_xx+WEI_vx(GI)*MUPTXX(JGL)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    V_yy=V_yy+WEI_vy(GI)*MUPTYY(JGL)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    V_zz=V_zz+WEI_vz(GI)*MUPTZZ(JGL)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    W_xx=W_xx+WEI_wx(GI)*MUPTXX(JGL)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    W_yy=W_yy+WEI_wy(GI)*MUPTYY(JGL)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    W_zz=W_zz+WEI_wz(GI)*MUPTZZ(JGL)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)

                    FU = FU+WEI_v(GI)*CORIOLIS(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    FV = FV-WEI_u(GI)*CORIOLIS(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

              if(isvd.eq.2) then
!                ewrite(3,*) 'CORIOLIS',CORIOLIS
!                stop 234
              endif
!                 U_zz=0.0
!                 V_zz=0.0
!                 W_zz=0.0
                    
                    DO JJ=1,NSVD
                       BIGM1(IBL11+COUNT,JJ)=BIGM1(IBL11+COUNT,JJ)+WEI_u(GI)*DETWEI(GI)*UD(GI,JJ)*LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                       BIGM1(IBL12+COUNT,JJ)=BIGM1(IBL12+COUNT,JJ)+WEI_u(GI)*DETWEI(GI)*VD(GI,JJ)*LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                       BIGM1(IBL13+COUNT,JJ)=BIGM1(IBL13+COUNT,JJ)+WEI_u(GI)*DETWEI(GI)*WD(GI,JJ)*LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)
                       
                       BIGM1(IBL21+COUNT,JJ)=BIGM1(IBL21+COUNT,JJ)+WEI_v(GI)*DETWEI(GI)*UD(GI,JJ)*LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                       BIGM1(IBL22+COUNT,JJ)=BIGM1(IBL22+COUNT,JJ)+WEI_v(GI)*DETWEI(GI)*VD(GI,JJ)*LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                       BIGM1(IBL23+COUNT,JJ)=BIGM1(IBL23+COUNT,JJ)+WEI_v(GI)*DETWEI(GI)*WD(GI,JJ)*LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                       
                       BIGM1(IBL31+COUNT,JJ)=BIGM1(IBL31+COUNT,JJ)+WEI_w(GI)*DETWEI(GI)*UD(GI,JJ)*LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                       BIGM1(IBL32+COUNT,JJ)=BIGM1(IBL32+COUNT,JJ)+WEI_w(GI)*DETWEI(GI)*VD(GI,JJ)*LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                       BIGM1(IBL33+COUNT,JJ)=BIGM1(IBL33+COUNT,JJ)+WEI_w(GI)*DETWEI(GI)*WD(GI,JJ)*LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                    ENDDO


                 ENDDO
              ENDDO

              IF(JSVD.LE.NSVD_PHI) THEN
                 DO JLOC=1,MLOC
                    JGLP=PNDGLN((ELE-1)*MLOC+JLOC)
                    DO GI=1,NGI
                       P_x=P_x+WEI_u(GI)*DETWEI(GI)*LEFTSVD(3*MXNODS+JGLP,JSVD)*MX(JLOC,GI)
                       P_y=P_y+WEI_v(GI)*DETWEI(GI)*LEFTSVD(3*MXNODS+JGLP,JSVD)*MY(JLOC,GI)
                       P_z=P_z+WEI_w(GI)*DETWEI(GI)*LEFTSVD(3*MXNODS+JGLP,JSVD)*MZ(JLOC,GI)
                    ENDDO
                 ENDDO
              ENDIF
             ! for geographic pressure
              IF(JSVD.LE.NSVD_PHI) THEN
                 DO JLOC=1,MPGLOC
                    JGLPG=PGNDGLNO((ELE-1)*MPGLOC+JLOC)
                    DO GI=1,NGI
                       PG1_x=PG1_x+WEI_u(GI)*DETWEI(GI)*LEFTSVD_PG_x(JGLPG,JSVD)*MPGX(JLOC,GI)
                       PG2_x=PG2_x+WEI_u(GI)*DETWEI(GI)*LEFTSVD_PG_y(JGLPG,JSVD)*MPGX(JLOC,GI)
                       PG1_y=PG1_y+WEI_v(GI)*DETWEI(GI)*LEFTSVD_PG_x(JGLPG,JSVD)*MPGY(JLOC,GI)
                       PG2_y=PG2_y+WEI_v(GI)*DETWEI(GI)*LEFTSVD_PG_y(JGLPG,JSVD)*MPGY(JLOC,GI)
                    ENDDO
                 ENDDO
              ENDIF


              ! Find count. ***************************************
              !                 LOWER=FINDRM(ISVD) 
              !                 UPPER=FINDRM(ISVD+1)-1
              ! 7000                 CONTINUE
              !                 INUM=LOWER+(UPPER-LOWER+1)/2 
              !                 IF(JSVD.GE.COLM(INUM) )  THEN 
              !                    LOWER=INUM
              !                 ELSE
              !                    UPPER=INUM
              !                 ENDIF
              !                 IF(UPPER-LOWER.LE.1) THEN
              !                    IF(JSVD.EQ.COLM(LOWER)) THEN
              !                       COUNT=LOWER
              !                    ELSE
              !                       COUNT=UPPER
              !                    ENDIF
              !                    GOTO 9000
              !                 ENDIF
              !                 GOTO 7000
              ! 9000                 CONTINUE

              !  Add into matrix and rhs...
              !--------------------------
 
              BIGM0(IBL11+COUNT)=BIGM0(IBL11+COUNT)+RN1+DT*THETA*(UDmean_Ux+VDmean_Uy+WDmean_Uz+U_xx+U_yy+U_zz)
              BIGM0(IBL11+COUNT)=BIGM0(IBL11+COUNT)+DT*THETA*PG2_x
              BIGM0(IBL12+COUNT)=BIGM0(IBL12+COUNT)+DT*THETA*FV
              BIGM0(IBL12+COUNT)=BIGM0(IBL12+COUNT)+DT*THETA*PG1_x
              BIGM0(IBL13+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                 BIGM0(IBL14+COUNT1)=BIGM0(IBL14+COUNT1)+DT*P_x
              ENDIF

              BIGM0(IBL21+COUNT)=BIGM0(IBL21+COUNT)+DT*THETA*FU
              BIGM0(IBL21+COUNT)=BIGM0(IBL21+COUNT)+DT*THETA*PG2_y
              BIGM0(IBL22+COUNT)=BIGM0(IBL22+COUNT)+RN2+DT*THETA*(UDmean_Vx+VDmean_Vy+WDmean_Vz+V_xx+V_yy+V_zz)
              BIGM0(IBL22+COUNT)=BIGM0(IBL22+COUNT)+DT*THETA*PG1_y
              BIGM0(IBL23+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                  BIGM0(IBL24+COUNT1)=BIGM0(IBL24+COUNT1)+DT*P_y
              ENDIF
              BIGM0(IBL31+COUNT)=BIGM0(IBL31+COUNT)+DT*THETA*PG2_z
              BIGM0(IBL32+COUNT)=BIGM0(IBL32+COUNT)+DT*THETA*PG1_z
              BIGM0(IBL33+COUNT)=BIGM0(IBL33+COUNT)+RN3+DT*THETA*(UDmean_Wx+VDmean_Wy+WDmean_Wz+W_xx+W_yy+W_zz)
              IF(JSVD.LE.NSVD_PHI) THEN
                  BIGM0(IBL34+COUNT1)=BIGM0(IBL34+COUNT1)+DT*P_z
              ENDIF

              IF(ISVD.LE.NSVD_PHI) THEN
                 BIGM0(IBL41+COUNT)=BIGM0(IBL41+COUNT)+U_x
                 BIGM0(IBL42+COUNT)=BIGM0(IBL42+COUNT)+V_y

! for 2D
                 BIGM0(IBL43+COUNT)=BIGM(IBL43+COUNT)+W_z
                 IF(JSVD.LE.NSVD_PHI) THEN
                    BIGM0(IBL44+COUNT1)=0.0
                 ENDIF

              ENDIF

              BMat0(IBL11+COUNT)=BMat0(IBL11+COUNT)+RN1-DT*ONEMINTH*(UDmean_Ux+VDmean_Uy+WDmean_Uz+U_xx+U_yy+U_zz)
              BMat0(IBL11+COUNT)=BMat0(IBL11+COUNT)-DT*ONEMINTH*PG2_x
              BMat0(IBL12+COUNT)=BMat0(IBL12+COUNT)-DT*ONEMINTH*FV
              BMat0(IBL12+COUNT)=BMat0(IBL12+COUNT)-DT*ONEMINTH*PG1_x
              BMat0(IBL13+COUNT)=0.0

 
              BMat0(IBL21+COUNT)=BMat0(IBL21+COUNT)-DT*ONEMINTH*FU
              BMat0(IBL21+COUNT)=BMat0(IBL21+COUNT)-DT*ONEMINTH*PG2_y

              BMat0(IBL22+COUNT)=BMat0(IBL22+COUNT)+RN2-DT*ONEMINTH*(UDmean_Vx+VDmean_Vy+WDmean_Vz+V_xx+V_yy+V_zz)
              BMat0(IBL22+COUNT)=BMat0(IBL22+COUNT)-DT*ONEMINTH*PG1_y
              BMat0(IBL23+COUNT)=0.0
              BMat0(IBL31+COUNT)=BMat0(IBL31+COUNT)-DT*ONEMINTH*PG2_z
              BMat0(IBL32+COUNT)=BMat0(IBL32+COUNT)-DT*ONEMINTH*PG1_z
              BMat0(IBL33+COUNT)=BMat0(IBL33+COUNT)+RN3-DT*ONEMINTH*(UDmean_Wx+VDmean_Wy+WDmean_Wz+W_xx+W_yy+W_zz)

              IF(ISVD.LE.NSVD_PHI) THEN
                 BMat0(IBL41+COUNT)=0.0
                 BMat0(IBL42+COUNT)=0.0
! for 2D
                 BMat0(IBL43+COUNT)=0.0
                 IF(JSVD.LE.NSVD_PHI) THEN
                    BMat0(IBL44+COUNT1)=0.0
                 ENDIF
              ENDIF
!

              BMat0(IBL11+COUNT)=BMat0(IBL11+COUNT)-DT*Umeanx_UD
              BMat0(IBL12+COUNT)=BMat0(IBL12+COUNT)-DT*Umeany_VD
              BMat0(IBL13+COUNT)=BMat0(IBL13+COUNT)-DT*Umeanz_WD
!
              BMat0(IBL21+COUNT)=BMat0(IBL21+COUNT)-DT*Vmeanx_UD
              BMat0(IBL22+COUNT)=BMat0(IBL22+COUNT)-DT*Vmeany_VD
              BMat0(IBL23+COUNT)=BMat0(IBL23+COUNT)-DT*Vmeanz_WD

              BMat0(IBL31+COUNT)=BMat0(IBL31+COUNT)-DT*Wmeanx_UD
              BMat0(IBL32+COUNT)=BMat0(IBL32+COUNT)-DT*Wmeany_VD
              BMat0(IBL33+COUNT)=BMat0(IBL33+COUNT)-DT*Wmeanz_WD

40            CONTINUE
!              ewrite(3,*) 'vec2',VEC(ISVD)
!              ewrite(3,*) 'BIGM(COUNT)',BIGM(COUNT)
!              if(ISVD.EQ.2) stop

20            CONTINUE
!              print *,'1',BIGM(COUNT)
30            CONTINUE


!400 CONTINUE
        write(name,'(a,i0)') 'PODBIGM0.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BIGM0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBMat0.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BMat0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODVEC0.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (VEC0(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBIGM1.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            write(1,*) (BIGM1(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)

ELSE  !NOT IF(IFWRITE.AND.LINITS.EQ.1)

        write(name,'(a,i0)') 'PODBIGM0.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (BIGM0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBMat0.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (BMat0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODVEC0.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (VEC0(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBIGM1.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            read(1,*) (BIGM1(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)
ENDIF  !END IF(IFWRITE.AND.LINITS.EQ.1) then


      ewrite(3,*) 'before solving'
   DO 100 NSTEP=1,ntime
      ewrite(3,*) 'NSTEP=',NSTEP
      BIGM(:)=BIGM0(:)
      BMat(:)=BMat0(:)
      VEC(:)=VEC0(:)

      DO ISVD=1,NSVD
!      ewrite(3,*) 'ISVD=',ISVD

         DO JSVD=1,NSVD
!         ewrite(3,*) 'JSVD=',JSVD
            COUNT=(ISVD-1)*NSVD+JSVD
            COUNT1=(ISVD-1)*NSVD_PHI+JSVD
            DO JJ=1,NSVD
               BIGM(IBL11+COUNT)=BIGM(IBL11+COUNT)+DT*THETA*( BIGM1(IBL11+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )

               BIGM(IBL22+COUNT)=BIGM(IBL22+COUNT)+DT*THETA*( BIGM1(IBL21+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )

               BIGM(IBL33+COUNT)=BIGM(IBL33+COUNT)+DT*THETA*( BIGM1(IBL31+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )


               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)-DT*ONEMINTH*( BIGM1(IBL11+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )


               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)-DT*ONEMINTH*( BIGM1(IBL21+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )

               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)-DT*ONEMINTH*( BIGM1(IBL31+COUNT,JJ)*NEWCOEF(JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*NEWCOEF(1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*NEWCOEF(2*NSVD+JJ) )
            ENDDO


               VEC(ISVD)=VEC(ISVD)+BMat(IBL11+COUNT)*NEWCOEF(JSVD)+                       &
                    BMat(IBL12+COUNT)*NEWCOEF(1*NSVD+JSVD)+BMat(IBL13+COUNT)*NEWCOEF(2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(ISVD)=VEC(ISVD)+BMat0(IBL14+COUNT1)*NEWCOEF(3*NSVD+JSVD)
               ENDIF

               VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat(IBL21+COUNT)*NEWCOEF(JSVD)+         &
                    BMat(IBL22+COUNT)*NEWCOEF(1*NSVD+JSVD)+BMat(IBL23+COUNT)*NEWCOEF(2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat0(IBL24+COUNT1)*NEWCOEF(3*NSVD+JSVD)
               ENDIF

               VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat(IBL31+COUNT)*NEWCOEF(JSVD)+         &
                    BMat(IBL32+COUNT)*NEWCOEF(1*NSVD+JSVD)+BMat(IBL33+COUNT)*NEWCOEF(2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat0(IBL34+COUNT1)*NEWCOEF(3*NSVD+JSVD)
               ENDIF
         ENDDO

      ENDDO

              if(.false.) then
              DO kk=1,NVAR
                 DO ii=1,NSVD
                    ftest=0.0
                    DO jj=1,NSVD
                       ReMat( (kk-1)*NSVD+ii,jj)=BIGM( (kk-1)*NVAR*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ReMat( (kk-1)*NSVD+ii,1*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+1*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ReMat( (kk-1)*NSVD+ii,2*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+2*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ReMat( (kk-1)*NSVD+ii,3*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+3*NSVD*NSVD+(ii-1)*NSVD+jj)
                       ftest=ftest+ReMat((kk-1)*NVAR*NSVD+ii,jj)
                    ENDDO
                    
                    if(abs(ftest).lt.1.0e-20) then
!                       ewrite(3,*) '(kk-1)*NVAR*NSVD+ii',kk,ii
!                       stop
                    endif
                 ENDDO
              ENDDO
              
              else

      open(1,file='bigm44.dat')
        write(1,*) IBL44
        write(1,*) (BIGM(IBL44+ii),ii=1,NSVD_PHI*NSVD_PHI)
      close(1)

              ReMat=0.0
              VEC1=0.0
              DO kk=1,NVAR
                 DO ii=1,NSVD
                    DO jj=1,NSVD
                       IF(kk.le.2) then
                          ReMat( (kk-1)*NSVD+ii,jj)=BIGM( (kk-1)*NSVD*NSVD_TOTAL+(ii-1)*NSVD+jj)
                          ReMat( (kk-1)*NSVD+ii,1*NSVD+jj)= BIGM( (kk-1)*NSVD*NSVD_TOTAL+1*NSVD*NSVD+(ii-1)*NSVD+jj)
                         IF(jj.LE.NSVD_PHI) THEN
                             ReMat( (kk-1)*NSVD+ii,2*NSVD+jj)= BIGM( (kk-1)*NSVD*NSVD_TOTAL+3*NSVD*NSVD+(ii-1)*NSVD_PHI+jj)
                          ENDIF
!                       ftest=ftest+ReMat((kk-1)*NVAR*NSVD+ii,jj)
                       ELSEIF((kk.eq.4).AND.(ii.LE.NSVD_PHI)) then
                          ReMat( (kk-2)*NSVD+ii,jj)=BIGM( (kk-1)*NSVD*NSVD_TOTAL+(ii-1)*NSVD+jj)
                          ReMat( (kk-2)*NSVD+ii,1*NSVD+jj)= BIGM( (kk-1)*NSVD*NSVD_TOTAL+1*NSVD_PHI*NSVD+(ii-1)*NSVD+jj)
                          IF(jj.LE.NSVD_PHI) THEN
                             ReMat( (kk-2)*NSVD+ii,2*NSVD+jj)= BIGM(IBL44+(ii-1)*NSVD_PHI+jj)
                          ENDIF
                       ENDIF
                    ENDDO
                    IF(kk.le.2) then
                       VEC1( (kk-1)*NSVD+ii)=VEC( (kk-1)*NSVD+ii)
                    ELSEIF((kk.eq.4).AND.(ii.LE.NSVD_PHI)) then
                       VEC1( (kk-2)*NSVD+ii)=VEC( (kk-1)*NSVD+ii)
                    ENDIF

                 ENDDO
              ENDDO

              endif

   if(.false.) then
              do ii=1,nvar*nsvd
                 ReMat(ii,ii)=ReMat(ii,ii)+vol
                 vec(ii)=vec(ii)+vol
              enddo
              close(1)
              close(2)
!              stop

              open(1,file='Remat2.dat')
              open(2,file='vec2.dat')
              write(1,*) Remat
              write(2,*) vec
              close(1)
              close(2)
              stop

   endif
!              IF(NSTEP.EQ.1) THEN
!              open(1,file='Remat1.dat')
!              open(2,file='vec1.dat')
!              write(1,*) Remat
!              write(2,*) vec
!              ELSEIF(NSTEP.EQ.3) THEN
!              open(1,file='Remat3.dat')
!              open(2,file='vec3.dat')
!              write(1,*) Remat
!              write(2,*) vec
!              STOP 156
!              ENDIF

              ! SOLVE matrix equations for gradients of hydrostic pressure
              ! *********DEFINE OPTIONS FOR GMRES**************

              ewrite(3,*) 'VEC'
              ewrite(3,*) 'BIGM'

              ERROR1 = 1.0E-14
              ERROR2 = ERROR1
              NOITS1 = 201
              NOITS2 = NOITS1
              IGUESS = .FALSE.
              MISOUT = 0
              IMM    = 1
              RELAX  = 1.0                
              LEFTP1 = 0
              LEFTP2 = 0
              MINTE1 = 20
              MINTE2 = 20       
              MULPA  = 0
              UZAWA  = 0
              CGSOLQ = 0
              GMRESQ = 1
              CGSQ   = 0 
              !       
              NR2=0
              IP=1
              MISOUT=1 
              NNV=7*NSVD-2*NSVD
 
              ! *********ENDOF OPTIONS FOR GMRES**************
              !       Solev MATRIX X=VEC


        open(10, file ='NEWCOEF.dat')
           write(10, *) (NEWCOEF1(ii), ii=1,(nvar-1)*nsvd)
        close(10) 


!           IF(2*int(nstep/2) .eq. nstep) THEN
              DO II=1,NSVD_TOTAL
                 IF(II.LE.2*NSVD) THEN
                    NEWCOEF(II)=NEWCOEF1(II)
                    IF(NT*int(nstep/NT) .eq. nstep) THEN
                       COEF(int(nstep/NT),II)=NEWCOEF1(II)
                    ENDIF
                 ELSEIF(II.GT.3*NSVD) THEN
                    NEWCOEF(II)=NEWCOEF1(II-NSVD)
                    IF(NT*int(nstep/NT) .eq. nstep) THEN
                       COEF(int(nstep/NT),II)=NEWCOEF1(II-NSVD)
                    ENDIF
                 ELSE
                    NEWCOEF(II)=0.0
                    IF(NT*int(nstep/NT) .eq. nstep) THEN
                       COEF(int(nstep/NT),II)=0.0
                    ENDIF
                 ENDIF
                 IF(NT*int(nstep/NT) .eq. nstep) THEN
                    ewrite(3,*) 'COEF(int(nstep/NT),II)',COEF(int(nstep/NT),II)
                 ENDIF
               ENDDO
!           ENDIF
              
              ewrite(3,*) 'after NSTEP',NSTEP
100           CONTINUE


              open(10, file ='COEF.dat')
              do ii=0,ntime1
                 write(10, *) (COEF(ii,jj), jj=1,nsvd_total)
              enddo
              close(10)
              ewrite(3,*) 'before deallocate'
!-------------------------------!
! DEALLOCATE THE LOCAL MATRICES !
!-------------------------------!

  DEALLOCATE(BIGM)
  DEALLOCATE(VEC)
  DEALLOCATE(ReMat)
  DEALLOCATE(BIGM0)
  DEALLOCATE(BIGM1)
  DEALLOCATE(VEC0)
  DEALLOCATE(VEC1)
  DEALLOCATE(BMat)
  DEALLOCATE(BMat0)
  DEALLOCATE(BMat01)
  DEALLOCATE(BMat1)
  DEALLOCATE(DETWEI)
  DEALLOCATE(UD)
  DEALLOCATE(VD)
  DEALLOCATE(WD)
  DEALLOCATE(UD_mean)
  DEALLOCATE(VD_mean)
  DEALLOCATE(WD_mean)
  DEALLOCATE(XD)
  DEALLOCATE(YD)
  DEALLOCATE(ZD)
  DEALLOCATE(NX)
  DEALLOCATE(NY)
  DEALLOCATE(NZ)
  DEALLOCATE(MX)
  DEALLOCATE(MY)
  DEALLOCATE(MZ)
  DEALLOCATE(WEI_u)
  DEALLOCATE(WEI_v)
  DEALLOCATE(WEI_w)
  DEALLOCATE(WEI_p)
  DEALLOCATE(WEI_ux)
  DEALLOCATE(WEI_uy)
  DEALLOCATE(WEI_uz)
  DEALLOCATE(WEI_vx)
  DEALLOCATE(WEI_vy)
  DEALLOCATE(WEI_vz)
  DEALLOCATE(WEI_wx)
  DEALLOCATE(WEI_wy)
  DEALLOCATE(WEI_wz)
  DEALLOCATE(NEWCOEF)
  DEALLOCATE(NEWCOEF1)
  DEALLOCATE(CORIOLIS)
  DEALLOCATE(SOURCX)
  DEALLOCATE(SOURCY)

  ewrite(3,*) 'in the middle of allocationg'
  ! for geographic pressure
  DEALLOCATE(MPGX)
  DEALLOCATE(MPGY)
  DEALLOCATE(MPGZ)
  DEALLOCATE(WEI_PGx)
  DEALLOCATE(WEI_PGy)
  DEALLOCATE(WEI_PGz)
  ewrite(3,*) 'bbbbbb'
!  DEALLOCATE(LEFTSVD_PG_x)
!  DEALLOCATE(LEFTSVD_PG_y)
  ewrite(3,*) 'cccccc'
  DEALLOCATE(PGMATRIX_x)
  DEALLOCATE(PGMATRIX_y)
  DEALLOCATE(PGMATRIX_mean)
  ewrite(3,*) 'dddddd'
  DEALLOCATE(PGVEC_x)
  DEALLOCATE(PGVEC_y)
  DEALLOCATE(PGVEC_mean)
             ewrite(3,*) 'exit REDUCED_MODEL_DIFFCOEF'
              RETURN
END SUBROUTINE REDUCED_MODEL_DIFFCOEF1_PG1




      SUBROUTINE OUTPUT_POD(NSTEP,NONODS,XNONOD,TOTELE,NLOC,X,Y,Z,NDGLNO)
      
!************************************************************
!| This subroutine is used to output the results from POD   |
!************************************************************
      include 'paramnew.h'

      INTEGER NONODS,XNONOD,TOTELE,NLOC,NSTEP
      REAL X(NONODS),Y(NONODS),Z(NONODS)
      REAL u0(0:ntimemax,mxpoi),v0(0:ntimemax,mxpoi),w0(0:ntimemax,mxpoi)
      REAL phi0(0:ntimemax,mxpoi)
      INTEGER NDGLNO(TOTELE*NLOC)
      CHARACTER*240 FILNAM
      CHARACTER*240 NAME
 
      INTEGER I,J,K
      common /fwdtraj/ phi0,u0,v0,w0   

      FILNAM='/home/fang/Gyre2-steady-POD'
      ewrite(3,*) FILNAM,nstep
      
      NAME=trim(FILNAM)//'.d'
      K = INDEX(NAME,' ') - 1
         write(name,'(a,i0)') trim(name)//'.', nstep
!           NAME=NAME(1:K)//'.'//CHARA(nstep)
         K = INDEX(NAME,' ') - 1
      ewrite(3,*) NAME(1:K)
      OPEN(1,FILE=NAME(1:K))
         WRITE(1,222) (FUNSIN( U0(nstep,i)),I=1,NONODS)
         WRITE(1,222) (FUNSIN( V0(nstep,i)),I=1,NONODS)
         WRITE(1,222) (FUNSIN( W0(nstep,i)),I=1,NONODS)
         WRITE(1,222) (FUNSIN( PHI0(nstep,i)),I=1,NONODS)

         WRITE(1,222) ( FUNSIN(X(I)),I=1,NONODS)
         WRITE(1,222) ( FUNSIN(Y(I)),I=1,NONODS)
         WRITE(1,222) ( FUNSIN(Z(I)),I=1,NONODS)
         WRITE(1,111) ( NDGLNO(I),I=1,TOTELE*NLOC)

         CLOSE(1)
111        FORMAT(10I9)
222     FORMAT(4E15.7)

      RETURN
      END SUBROUTINE OUTPUT_POD

        FUNCTION FUNSIN(RVEC)
        REAL FUNSIN
        REAL RVEC
        IF(ABS(RVEC).LT.1.E-30) THEN
           FUNSIN=0.
        ELSE
           FUNSIN=RVEC
        ENDIF
        RETURN
        END FUNCTION



SUBROUTINE VTKOUTPUT1P(TOTELE,NONODS,NTsol,NLoc, &
     D3, &
     NDGLNO, &
     X,Y,Z, &
     DATA1,Data2,Data3,Data4,DataM,&
     Task, &
     DumNo, filename)
  
  INTEGER, Intent(In)::TOTELE,NONODS,Ntsol,Nloc
  Logical, Intent(In)::D3
  INTEGER, Intent(In)::NDGLNO(TOTELE*NLOC)

  REAL, Intent(In)::X(NONODS),Y(NONODS),Z(NONODS)
  
  REAL, Intent(In)::DATA1(NONODS),Data2(Nonods), Data3(Nonods) 
  Real, Intent(In)::Data4(Nonods), DataM(ntsol*Nonods)

  INTEGER, Intent(InOut)::DUMNO
  character, intent(in) ::  filename*240
!  Character filename*7
  Integer, Intent(In)::Task
  integer lname1, lname2, lname3, lname4,lname5,lname6,lenvec
  character*240 name1, name2, name3, name4,name5,name6,namevec

! Local stuff...
  INTEGER, Allocatable::EleType(:), EleSize(:)
  Real, Allocatable::ZStore(:),data5(:),data6(:)
  CHARACTER*240 outName, vtkTitle

  Integer OutNameLth,vtkTitleLth,k,i,it

  if(ntsol.gt.2) then
     write(3,*) 'ntsol gt 2 not supported yet..'
     stop 88
  end if

  Allocate(data5(Nonods))
  if(ntsol.eq.2) then
     Allocate(data6(Nonods))
  end if

  do i=1,Nonods
     data5(i) = dataM(i)
     if(ntsol.eq.2) then
        data6(i) = dataM(Nonods+i)
     end if
  end do
   write(3,*) 'Inside subroutine VtkOutput1P'
  k=index(filename,' ')-1
     
  Allocate( EleType(TotEle) )
  Allocate( EleSize(TotEle) )
  Allocate( ZStore(Nonods))

  if(D3) then
     do i=1,totele
        eletype(i) = 10
        elesize(i) = 4
     end do
     
     do i=1,nonods
        ZStore(i) = Z(i)   
     end do
     
  else
! adjusted for triangles..
     do i=1,totele
        eletype(i) = 5
        elesize(i) = 3
     end do

     do i=1,nonods
        ZStore(i) = 0.0
     end do

  end if

 
!  CHECK(DUMNO)

! Write the file names and descriptions depending on OutOpt..  

  WRITE(outName, FMT='(A,A,I5.5,A)') filename(1:k),"_",DUMNO,".vtu"
  OutNameLth = k + 10

  ewrite(3,*) 'outName=',outName

  WRITE(vtkTitle, FMT='(A)') "Output field"
  vtkTitleLth = 12

  WRITE(namevec, FMT='(A)') "vectors"
  lenvec = 7

  if(task.eq.1) then
     WRITE(name1, FMT='(A)') "FwdU"
     WRITE(name2, FMT='(A)') "FwdV"
     WRITE(name3, FMT='(A)') "FwdW"
     WRITE(name4, FMT='(A)') "FwdP"
     WRITE(name5, FMT='(A)') "FwdT1"
     WRITE(name6, FMT='(A)') "FwdT2"
     lname1=4
     lname2=4
     lname3=4
     lname4=4
     lname5=5
     lname6=5
  else if(task.eq.2) then
     WRITE(name1, FMT='(A)') "AdjU"
     WRITE(name2, FMT='(A)') "AdjV"
     WRITE(name3, FMT='(A)') "AdjW"
     WRITE(name4, FMT='(A)') "AdjP"
     WRITE(name5, FMT='(A)') "AdjT1"
     WRITE(name6, FMT='(A)') "AdjT2"
     lname1=4
     lname2=4
     lname3=4
     lname4=4
     lname5=5
     lname6=5
  else
     stop 9666
  end if
     
  CALL VTKOPEN(outName, OutNameLth, vtkTitle, vtkTitleLth)
!  CALL VTKOPEN('fred.vtu', len('fred'), 'fred', len('fred'))
 
!#ifdef DOUBLEP


if(.true.) then
  write(3,*) 'in double prec'

  CALL VTKWRITEMESHD(NONODS, TOTELE, X, Y, ZStore, NDGLNO, EleType, EleSize)
    write(3,*) 'written mesh'
  CALL VTKWRITEDSN(data1, name1, lname1)
    write(3,*) 'written data1'
  CALL VTKWRITEDSN(data2, name2, lname2)
    write(3,*) 'written data2'
  CALL VTKWRITEDSN(data3, name3, lname3)
    write(3,*) 'written data3'
  CALL VTKWRITEDSN(data4, name4, lname4)
   write(3,*) 'written data4'
   CALL VTKWRITEDSN(data5, name5, lname5)
    write(3,*) 'written data5'
   if(ntsol.eq.2) then
      CALL VTKWRITEDSN(data6, name6, lname6)
       write(3,*) 'written data6'
   end if
   Call VTKWRITEDVN(data1, data2, data3, namevec, lenvec)

else

  ewrite(3,*) 'in single p'

  CALL VTKWRITEMESH(NONODS, TOTELE, X, Y, ZStore, NDGLNO, EleType, EleSize)
  
  CALL VTKWRITEFSN(data1, name1, lname1)
  CALL VTKWRITEFSN(data2, name2, lname2)
  CALL VTKWRITEFSN(data3, name3, lname3)
  CALL VTKWRITEFSN(data4, name4, lname4)
  CALL VTKWRITEFSN(data5, name5, lname5)
  if(ntsol.eq.2) then
     CALL VTKWRITEFSN(data6, name6, lname6)
  end if
  Call VTKWRITEFVN(data1, data2, data3, namevec, lenvec)
endif

  CALL VTKCLOSE()

  Deallocate(ZStore)
  Deallocate(EleType)
  Deallocate(EleSize, data5)
  
  if(ntsol.eq.2) then
     Deallocate(data6)
  end if

  Return
END SUBROUTINE VTKOUTPUT1P


SUBROUTINE Source_POD(NONODS,LTIME,NSVD,DT,NTIME,ITSTIME,                         &
                    SOURCX,SOURCY,SOURCZ,                                         &
                    SMEAN,LEFTSVD,COEF,                                           &
                    TOTELE,NLOC,NGI,MXNODS,                                       &
                    NDGLNO,XNONOD,xondgl,X,Y,Z,                                   &
                    N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)

  INTEGER, INTENT(IN) ::NONODS,XNONOD,NTIME,ITSTIME,NSVD,MXNODS
  LOGICAL  D3,DCYL
  INTEGER  TOTELE,NLOC,NGI
  REAL     N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
  INTEGER  NDGLNO(TOTELE*NLOC),XONDGL(TOTELE*NLOC)
  REAL     WEIGHT(NGI)
  REAL, INTENT(IN) ::LTIME,DT
  REAL,        INTENT(OUT),DIMENSION(NSVD) :: SOURCX,SOURCY,SOURCZ
  REAL, INTENT(IN),DIMENSION(NONODS) ::X,Y,Z
  REAL, INTENT(IN),DIMENSION(3*MXNODS,NSVD) :: LEFTSVD
  REAL, INTENT(IN),DIMENSION(3*MXNODS) ::SMEAN
  REAL, INTENT(IN),DIMENSION(3*NSVD) ::COEF
  ! Local variables.....
  INTEGER      ::NOSTIM
  REAL            ::RRmean_X,RRmean_Y,RRmean_Z,RRmean_XT,RRmean_YT,RRmean_ZT
  REAL            ::RR_X,RR_Y,RR_Z,RR_XT,RR_YT,RR_ZT
  INTEGER   ::NUMCASE
  REAL,PARAMETER ::PIE = 3.141592654
  REAL      ::ALPHA1,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  INTEGER   ::NOSNOD
  REAL      ::WEI_u(NGI),WEI_v(NGI),WEI_w(NGI),DETWEI(NGI)
  REAL      ::VOLUME
  REAL      ::NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
  REAL,DIMENSION(:,:),ALLOCATABLE ::UEXAC,VEXAC,WEXAC
  REAL,DIMENSION(:),ALLOCATABLE ::SX,SY,SZ
  REAL,DIMENSION(:),ALLOCATABLE ::U,V,W
  REAL,DIMENSION(:),ALLOCATABLE ::STIME
  REAL,DIMENSION(:),ALLOCATABLE::ADWEI
  REAL,DIMENSION(:),ALLOCATABLE ::SOX,SOY,SOZ
  REAL            ::GLOWEI,GLOWEI_X,GLOWEI_Y,GLOWEI_Z,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
  INTEGER   ::I,J,K1,K2
  INTEGER   ::ELE,IGL,JGL,ILOC,JLOC,GI,ISVD,JSVD
  INTEGER   ::ii,jj,kk,k0
  CHARACTER*240 ::FILNAM_SOURCE

  INTEGER   ::SNOD_LIST_NO
  REAL,DIMENSION(:),ALLOCATABLE ::ML
  INTEGER,DIMENSION(:),ALLOCATABLE ::SNOD_LIST
  
  OPEN(3,FILE='GaussianSpread.dat',STATUS='OLD')
  READ(3,*) NOSNOD,NOSTIM,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  write(3,*) 'NOSNOD,NOSTIM,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
       NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  READ(3,*) NUMCASE,ALPHA1
  ewrite(3,*) 'NUMCASE,ALPHA1',NUMCASE,ALPHA1
  CLOSE(3)

  ALLOCATE(ML(NONODS))
  ALLOCATE(UEXAC(NOSNOD,NOSTIM))  
  ALLOCATE(VEXAC(NOSNOD,NOSTIM))  
  ALLOCATE(WEXAC(NOSNOD,NOSTIM))  
  ALLOCATE(U(NOSNOD))  
  ALLOCATE(V(NOSNOD))  
  ALLOCATE(W(NOSNOD))  
  ALLOCATE(SX(NOSNOD))  
  ALLOCATE(SY(NOSNOD))  
  ALLOCATE(SZ(NOSNOD))  
  ALLOCATE(STIME(NOSTIM))  
  ALLOCATE(SOX(NOSNOD))  
  ALLOCATE(SOY(NOSNOD))  
  ALLOCATE(SOZ(NOSNOD))  
  ALLOCATE(SNOD_LIST(NONODS))

  ML=0.0
  UEXAC=0.0
  VEXAC=0.0
  WEXAC=0.0
  U=0.0
  V=0.0
  W=0.0
  SX=0.0
  SY=0.0
  SZ=0.0
  STIME=0.0
  SOX=0.0
  SOY=0.0
  SOZ=0.0

  SOURCX=0.0
  SOURCY=0.0
  SOURCZ=0.0
  SNOD_LIST=0

  CALL CALCULATE_ML(ML,NONODS,XNONOD,xondgl,                           &
                   TOTELE,NLOC,NGI,                                    &
                   NDGLNO,X,Y,Z,                                       &
                   N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)                      

!**********************************
! get the observational data..... !
!**********************************


  FILNAM_SOURCE='AdjSourceData.dat'
  CALL EXASOL_POD(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,                    &
              SX,SY,SZ,STIME,.true.,FILNAM_SOURCE,0.0,1.0,.false.)

  OPEN(1,file='exacxsolution.dat')
   WRITE(1,*) 'UEXAC'
   DO KK=1,NOSTIM
     WRITE(1,*) 'NOSTIM=',NOSTIM
     WRITE(1,*) (UEXAC(II,KK),II=1,NOSNOD)
   ENDDO
   WRITE(1,*) 'VEXAC'
   DO KK=1,NOSTIM
     WRITE(1,*) 'NOSTIM=',NOSTIM
     WRITE(1,*) (VEXAC(II,KK),II=1,NOSNOD)
   ENDDO
  CLOSE(1)

  ewrite(3,*) 'after reading the obseravational data'

!-------
!**********************************
! calculate the pod u,v,w   ..... !
!**********************************

  DO II=1,NONODS
     U(II)=SMEAN(II)
     V(II)=SMEAN(1*MXNODS+II)
     W(II)=SMEAN(2*MXNODS+II)
     DO JJ=1,NSVD
       U(II)=U(II)+ LEFTSVD(II,JJ)*COEF(JJ)
       V(II)=V(II)+ LEFTSVD(1*MXNODS+II,JJ)*COEF(1*NSVD+JJ)
       W(II)=W(II)+ LEFTSVD(2*MXNODS+II,JJ)*COEF(2*NSVD+JJ)
     ENDDO
  ENDDO


!-----------------------------
  SELECT CASE(NUMCASE)
   

! 1D case
!...........

  CASE(1)

  ALLOCATE(ADWEI(NOSNOD))

  GLOWEI_X1 = 0.0
  DO JJ = 1,NONODS
     RR_XT = 0.0
     RR_YT = 0.0
     RR_ZT = 0.0
     ADWEI=0.0
     GLOWEI = 0.0
     RR_X = 0.0
     DO II=1,NOSNOD
        ADWEI(II) = EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

        RR_X =EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*      &
             ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )
        RR_XT =RR_XT+RR_X   !RR_XT =RR_X in this case
        GLOWEI = GLOWEI+ADWEI(II)                                                
     ENDDO   ! end II
     GLOWEI_X1 = GLOWEI_X1+GLOWEI
     SOX(JJ) = RR_XT
  ENDDO    ! end jj

  DEALLOCATE(ADWEI)

  DO JJ=1,NONODS
     IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
        SOX(JJ) = ALPHA1*SOX(JJ)/GLOWEI_X1
     ENDIF
  ENDDO

! 2D case
  CASE(2)
  ALLOCATE(ADWEI(NOSNOD))

  GLOWEI_X1 = 0.0
  SNOD_LIST_NO=0
  SNOD_LIST=0
  DO JJ = 1,NONODS
     RRmean_XT = 0.0
     RRmean_YT = 0.0
     RRmean_ZT = 0.0
     RR_XT = 0.0
     RR_YT = 0.0
     RR_ZT = 0.0
     ADWEI(1:NOSNOD)=0.0
     GLOWEI = 0.0
     RRmean_X = 0.0
     RRmean_Y = 0.0
     RRmean_Z = 0.0
     RR_X = 0.0
     RR_Y = 0.0
     RR_Z = 0.0
     DO II=1,NOSNOD

        ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

        RR_X =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
             ( U(JJ)-UEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
        RR_Y =EXP(-((Y(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
             ( V(JJ)-VEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
        RR_Z =EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*      &
             ( W(JJ)-WEXAC(II,NTIME-ITSTIME+1) )  ! NOTE:NTIME-ITSTIME+1
        RR_XT =RR_XT+RR_X   !RR_T =RR_X in this case
        RR_YT =RR_YT+RR_Y   !RR_T =RR_X in this case
        RR_ZT =RR_ZT+RR_Z   !RR_T =RR_X in this case

        GLOWEI = GLOWEI+ADWEI(II)
! to judge whether the nonods is been taked into account the list of nodes around the detectors
        IF( SQRT((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2).LT.3*ISPAN_X.AND.         &
             SNOD_LIST(SNOD_LIST_NO).NE.JJ) THEN
           SNOD_LIST_NO=SNOD_LIST_NO+1
           SNOD_LIST(SNOD_LIST_NO)=JJ
        ENDIF
             
     ENDDO   ! end II


     GLOWEI_X1 = GLOWEI_X1+GLOWEI

     SOX(JJ) = RR_XT
     SOY(JJ) = RR_YT
     IF(D3) THEN
        SOZ(JJ) = RR_ZT
     ENDIF
  ENDDO    ! end jj


  DEALLOCATE(ADWEI)

  DO JJ=1,NONODS
     IF( (ABS(GLOWEI_X1)-0.0).GT.1.0E-12) THEN
      open(2,file='sourcx1.dat')
      write(2,*) ALPHA1, GLOWEI_X1,sox(JJ)
        SOX(JJ) = ALPHA1*SOX(JJ)/GLOWEI_X1
        SOY(JJ) = ALPHA1*SOY(JJ)/GLOWEI_X1
        SOZ(JJ) = ALPHA1*SOZ(JJ)/GLOWEI_X1
      open(1,file='sourcx.dat')
      write(1,*) sox(JJ)
                
     ENDIF
  ENDDO
     close(1)
     close(2)

! 3D case

  CASE DEFAULT

  END SELECT


  DO 30 ELE=1,TOTELE
     DO ILOC=1,NLOC
        JGL=NDGLNO((ELE-1)*NLOC+ILOC)
        DO II=1,SNOD_LIST_NO
           IF(JGL.EQ.SNOD_LIST(SNOD_LIST_NO)) GOTO 10
        ENDDO
        GOTO 30
     ENDDO
10   CONTINUE
     CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,      &    
          N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,               &
          NX,NY,NZ) 

     DO 20 ISVD=1,NSVD
        ! Calcultae DETWEI,VOLUME...
        DO GI=1,NGI
           WEI_u(GI)=0.0
           WEI_v(GI)=0.0
           WEI_w(GI)=0.0
           DO ILOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+ILOC)
              
              WEI_u(GI)=WEI_u(GI)+LEFTSVD(JGL,ISVD)*N(ILOC,GI)
              WEI_v(GI)=WEI_v(GI)+LEFTSVD(1*MXNODS+JGL,ISVD)*N(ILOC,GI)
              WEI_w(GI)=WEI_w(GI)+LEFTSVD(2*MXNODS+JGL,ISVD)*N(ILOC,GI)
           ENDDO
        ENDDO
        
        SOURCX=0.
        SOURCY=0.
        SOURCZ=0.
        DO ILOC=1,NLOC
           IGL = NDGLNO((ELE-1)*NLOC + ILOC)
           DO GI=1,NGI
              SOURCX(ISVD) =SOURCX(ISVD) +WEI_u(GI)*N(ILOC,GI) *DETWEI(GI)*SOX(IGL)
              SOURCY(ISVD) =SOURCY(ISVD) +WEI_v(GI)*N(ILOC,GI) *DETWEI(GI)*SOY(IGL)
              SOURCZ(ISVD) =SOURCZ(ISVD) +WEI_w(GI)*N(ILOC,GI) *DETWEI(GI)*SOZ(IGL)
           END DO
        ENDDO
20      CONTINUE ! END ELE loop
30      CONTINUE !END ISVD LOOP
     DEALLOCATE(ML)
     DEALLOCATE(UEXAC)  
     DEALLOCATE(VEXAC)  
     DEALLOCATE(WEXAC)  
     DEALLOCATE(U)  
     DEALLOCATE(V)  
     DEALLOCATE(W)  
     DEALLOCATE(SX)  
     DEALLOCATE(SY)  
     DEALLOCATE(SZ)  
     DEALLOCATE(STIME)  
     DEALLOCATE(SOX)  
     DEALLOCATE(SOY)  
     DEALLOCATE(SOZ)  
     DEALLOCATE(SNOD_LIST)

END SUBROUTINE Source_POD


SUBROUTINE CALCULATE_ML(ML,NONODS,XNONOD,xondgl,                         &
                   TOTELE,NLOC,NGI,                                    &
                   NDGLNO,X,Y,Z,                                       &
                   N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)                      

      LOGICAL  D3,DCYL
      INTEGER  TOTELE,NLOC,NGI
      INTEGER  NONODS,XNONOD
      REAL     ML(NONODS)
      INTEGER  NDGLNO(TOTELE*NLOC),XONDGL(TOTELE*NLOC)
      REAL     X(XNONOD),Y(XNONOD),Z(XNONOD)
      REAL     N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL     WEIGHT(NGI)
      REAL     RNN,VOLUME
      INTEGER  ELE,ILOC,JLOC,IGL,JGL,GI
      REAL     NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
      REAL,DIMENSION(:),ALLOCATABLE  :: MMASS,MML,DETWEI

      ALLOCATE( DETWEI(NGI) )
      DETWEI(1:NGI)=0.0

         ML(1:NONODS)=0.0
         DO  ELE=1,TOTELE
 ! Calcultae DETWEI,VOLUME...
           CALL DETNLX(ELE, X,Y,Z, XONDGL, TOTELE,XNONOD,NLOC,NGI,     &
                    N,NLX,NLY,NLZ, WEIGHT, DETWEI,VOLUME, D3,DCYL,     &
                    NX,NY,NZ) 
          DO ILOC=1,NLOC
           IGL = NDGLNO((ELE-1)*NLOC + ILOC)
             DO JLOC=1,NLOC
               JGL=NDGLNO((ELE-1)*NLOC+JLOC) 
               RNN=0.
               DO GI=1,NGI
                 RNN =RNN +N(ILOC,GI)*N(JLOC,GI) *DETWEI(GI)
               END DO
               ML(IGL)=ML(IGL)+RNN
             ENDDO
           ENDDO
         ENDDO

      DEALLOCATE( DETWEI )

END SUBROUTINE CALCULATE_ML


SUBROUTINE EXASOL_POD(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,SX,SY,SZ,STIME,D3,FILNAM_SOURCE,  &
     FSZERO,FSSCALE,INOUT)
  
  !----------------------------------------------------
  ! This subroutine is used to get the exact solutions
  !----------------------------------------------------
  
  INTEGER ::NOSTIM,NOSNOD,ITSTIME
  REAL    ::UEXAC(NOSNOD,NOSTIM),VEXAC(NOSNOD,NOSTIM),WEXAC(NOSNOD,NOSTIM)
  REAL    ::SX(NOSNOD),SY(NOSNOD),SZ(NOSNOD)
  REAL    ::STIME(NOSTIM),ST,DT,FSZERO,FSSCALE
  LOGICAL ::D3
  LOGICAL ::INOUT    !true--writeout; false--read data
  CHARACTER(40) ::UNAME,VNAME,WNAME,SXNAME,SYNAME,SZNAME,STIMENAME
  CHARACTER*240 ::FILNAM_SOURCE
  ! Local variables.......
  INTEGER ::II,KK,I
  REAL    ::AA
  
  
  STIME=0.0
  
  DO KK =2,NOSTIM
     STIME(KK) = STIME(KK-1)+DT
  ENDDO
  
  IF(INOUT) THEN

    OPEN(1,FILE=FILNAM_SOURCE)
    ! Read the obervational data on free surface.... ignor the height change for the time being
    
    DO II=1,NOSNOD
        READ(1,*) UNAME
        write(3,*) UNAME
        WRITE(1,*) SX(II),SY(II)
        DO KK =1,NOSTIM
           WRITE(1,*) ST,UEXAC(II,KK)
           UEXAC(II,KK)=UEXAC(II,KK)*FSSCALE+FSZERO
        ENDDO
    ENDDO
 
    CLOSE(1)

 ELSE

  OPEN(1,FILE=FILNAM_SOURCE,STATUS='OLD')
  ! Read the obervational data on free surface.... ignor the height change for the time being

    DO II=1,NOSNOD
        READ(1,*) UNAME
        READ(1,*) SX(II),SY(II)
        DO KK =1,NOSTIM
           READ(1,*) ST,UEXAC(II,KK)
           UEXAC(II,KK)=UEXAC(II,KK)*FSSCALE+FSZERO
        ENDDO
    ENDDO
 
  CLOSE(1)

 ENDIF

END SUBROUTINE EXASOL_POD


SUBROUTINE WINDYSOURCE(SOURCX,SOURCY,X,Y,NONODS,NSVD)
  !----------------------------------------------------
  ! This subroutine is used to get wind stress
  !----------------------------------------------------
  
  INTEGER ::NONODS,NSVD
  REAL    ::SOURCX(NONODS),SOURCY(NONODS)
  REAL    ::VECX(NSVD),VECY(NSVD)
  REAL    ::X(NONODS),Y(NONODS)
  REAL    ::T0,PI
  INTEGER ::II,JJ,KK,NOD

  PI = 3.141592653
  T0 = 163.2653

  DO NOD=1,NONODS        
     SOURCX(NOD) = -T0*COS(PI*Y(NOD))
  ENDDO

END SUBROUTINE WINDYSOURCE


     !.......................

SUBROUTINE FUNCTION_POD(FUNCT,NONODS,LTIME,DT,NTIME,ITSTIME,ACCTIM,    &
           U,V,W,T,X,Y,Z,                                             &
           TOTELE,NLOC,NGI,                                           &
           NDGLNO,XNONOD,xondgl,                                      &
           N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)


  !******************************************************************************************
  ! This subroutine is used to calculate the objective function
  ! Variables:
  !******************************************************************************************
  !

  REAL, INTENT(OUT) ::FUNCT
  INTEGER, INTENT(IN) ::NONODS,XNONOD,NTIME,ITSTIME
  REAL, INTENT(IN) ::LTIME,DT,ACCTIM
  REAL, INTENT(IN),DIMENSION(NONODS) ::U,V,W,T
  REAL, INTENT(IN),DIMENSION(NONODS) ::X,Y,Z
  LOGICAL  D3,DCYL
  INTEGER  TOTELE,NLOC,NGI
  REAL     N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
  INTEGER  NDGLNO(TOTELE*NLOC),XONDGL(TOTELE*NLOC)
  REAL     WEIGHT(NGI)

  ! Local variables.......
  INTEGER      ::NOSTIM
  INTEGER      ::NOSNOD
  REAL         ::ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  REAL,PARAMETER ::PIE=3.141592654
  REAL         ::ALPHA1
  INTEGER ::NUMCASE
  REAL          ::FUNCT_X,FUNCT_T,FUNCT_XT,FUNCT0
  REAL          ::GLOWEI_X,GLOWEI_T,GLOWEI_XT,GLOWEI_X1
  CHARACTER*240 ::FILNAM_SOURCE
  REAL      ::VOLUME
  REAL      ::NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
  REAL,DIMENSION(:),ALLOCATABLE::ADWEI
  REAL,DIMENSION(:),ALLOCATABLE::UEXAC_WEI
  REAL,DIMENSION(:),ALLOCATABLE::VEXAC_WEI
  REAL,DIMENSION(:),ALLOCATABLE::WEXAC_WEI

  REAL,DIMENSION(:,:),ALLOCATABLE ::UEXAC,VEXAC,WEXAC
  REAL,DIMENSION(:),ALLOCATABLE ::SX,SY,SZ
  REAL,DIMENSION(:),ALLOCATABLE ::STIME
  INTEGER   ::SNOD_LIST_NO
  REAL,DIMENSION(:),ALLOCATABLE ::ML
  INTEGER,DIMENSION(:),ALLOCATABLE ::SNOD_LIST

  INTEGER,PARAMETER ::REGU=1  ! if regulazation of the cost function
  REAL,PARAMETER    ::LAMDA = 0.0  !0.1 
  REAL    ::VARI_WEI,VARI
  INTEGER ::II,KK,JJ,MM
  INTEGER   ::I,J,K,K1,K2
  INTEGER   ::ELE,IGL,JGL,ILOC,JLOC,GI,ISVD,JSVD


  ALLOCATE(ML(NONODS))
  ALLOCATE(UEXAC(NOSNOD,NOSTIM))  
  ALLOCATE(VEXAC(NOSNOD,NOSTIM))  
  ALLOCATE(WEXAC(NOSNOD,NOSTIM))  
  ALLOCATE(SX(NOSNOD))  
  ALLOCATE(SY(NOSNOD))  
  ALLOCATE(SZ(NOSNOD))  
  ALLOCATE(STIME(NOSTIM))  
  ALLOCATE(SNOD_LIST(NONODS))
  ML=0.0
  UEXAC=0.0
  VEXAC=0.0
  WEXAC=0.0
  SX=0.0
  SY=0.0
  SZ=0.0
  STIME=0.0
  SNOD_LIST =0

  OPEN(3,FILE='GaussianSpread.dat',STATUS='OLD')
  READ(3,*) NOSNOD,NOSTIM,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  ewrite(3,*) 'NOSNOD,NOSTIM,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z',                &
       NOSNOD,ISPAN_T,ISPAN_X,ISPAN_Y,ISPAN_Z
  READ(3,*) NUMCASE,ALPHA1
  ewrite(3,*) 'NUMCASE,ALPHA1',NUMCASE,ALPHA1
  CLOSE(3)


  CALL CALCULATE_ML(ML,NONODS,XNONOD,xondgl,                           &
                   TOTELE,NLOC,NGI,                                    &
                   NDGLNO,X,Y,Z,                                       &
                   N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)                      

!**********************************
! get the observational data..... !
!**********************************

  FILNAM_SOURCE='AdjSourceData.dat'
  CALL EXASOL_POD(NOSNOD,NOSTIM,DT,UEXAC,VEXAC,WEXAC,                    &
              SX,SY,SZ,STIME,.true.,FILNAM_SOURCE,0.0,1.0,.false.)

  ewrite(3,*) 'after reading the obseravational data'


! Calculate the cost function FUNCT: which is the functional gauge of data miss fit.
!***********************************************

    SELECT CASE (NUMCASE)


! Method ONE
!-----------
      CASE (1)
      !In this case, One-dimensional

  FUNCT0 =0.0
  GLOWEI_XT = 0.0
  ALLOCATE(ADWEI(NOSNOD))
  DO JJ=1,NONODS
        FUNCT_XT = 0.0
        ADWEI(1:NOSNOD)=0.0
        GLOWEI_X = 0.0
        DO II=1,NOSNOD
! this is for 1D for the time being, we will fix it later on.... 
          IF( ABS(Z(JJ)-SZ(II)).LT.1.0E-6) THEN
            ADWEI(II) = EXP(-((X(JJ)-SX(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

             FUNCT_X = 0.5*( (U(JJ)-UEXAC(II,ITSTIME))**2 )
           FUNCT_XT=FUNCT_XT+ADWEI(II)*FUNCT_X
           GLOWEI_X = GLOWEI_X+ADWEI(II)
         ENDIF !IF( ABS(Z(JJ)-SZ(II)).LT.1.0E-6) THEN
       ENDDO  ! end ii
       GLOWEI_XT=GLOWEI_XT+GLOWEI_X
         FUNCT0 =FUNCT0+FUNCT_XT
  ENDDO     ! end jj

  DEALLOCATE(ADWEI)
  IF( (ABS(GLOWEI_XT)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_XT
  ENDIF

 FUNCT= FUNCT+ALPHA1*FUNCT0 
 ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0

! if regulation of the cost function...
  FUNCT0 =0.0
  GLOWEI_X = 0.0
  DO JJ=1,NONODS
     IF(ABS(X(JJ)).LT.1.0E-6) THEN
             FUNCT_X = 0.5*( U(JJ)**2 )     
           FUNCT0=FUNCT0+ABS(ML(JJ))*FUNCT_X
           GLOWEI_X = GLOWEI_X+ABS(ML(JJ))
     ENDIF
  ENDDO     ! end jj

  IF( (ABS(GLOWEI_X)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_X
  ENDIF

  FUNCT= FUNCT+LAMDA*FUNCT0 
  ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0


! Method TWO
!-----------
      CASE (2)
!In this case, the method isthe same as the method two 
!except for without spreading in time
  FUNCT0 =0.0
  GLOWEI_XT = 0.0
  ALLOCATE(ADWEI(NOSNOD))
  DO JJ=1,NONODS
        FUNCT_XT = 0.0
        ADWEI(1:NOSNOD)=0.0
        GLOWEI_X = 0.0
        DO II=1,NOSNOD
! this is for 2D for the time being, we will fix it later on.... 
          IF( ABS(Z(JJ)-SZ(II)).LT.1.0E-6) THEN
            ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

           IF(D3) THEN
             FUNCT_X = 0.5*( (U(JJ)-UEXAC(II,ITSTIME))**2+(V(JJ)-VEXAC(II,ITSTIME))**2+       &
                             (W(JJ)-WEXAC(II,ITSTIME))**2 )
           ELSE
             FUNCT_X = 0.5*( (U(JJ)- UEXAC(II,ITSTIME))**2+(V(JJ)-VEXAC(II,ITSTIME))**2)
           ENDIF
           FUNCT_XT=FUNCT_XT+ADWEI(II)*FUNCT_X
           GLOWEI_X = GLOWEI_X+ADWEI(II)
         ENDIF !IF( ABS(Z(JJ)-SZ(II)).LT.1.0E-6) THEN
       ENDDO  ! end ii
       GLOWEI_XT=GLOWEI_XT+GLOWEI_X
         FUNCT0 =FUNCT0+FUNCT_XT
  ENDDO     ! end jj

  DEALLOCATE(ADWEI)
  IF( (ABS(GLOWEI_XT)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_XT
  ENDIF

 FUNCT= FUNCT+ALPHA1*FUNCT0 
 ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0

! if regulation of the cost function...
  FUNCT0 =0.0
  GLOWEI_X = 0.0
  DO JJ=1,NONODS
     IF(ABS(X(JJ)).LT.1.0E-6) THEN
          IF(D3) THEN
             FUNCT_X = 0.5*( U(JJ)**2+V(JJ)**2+W(JJ)**2 )     
           ELSE
             FUNCT_X = 0.5*( U(JJ)**2+V(JJ)**2 )     
           ENDIF
           FUNCT0=FUNCT0+ABS(ML(JJ))*FUNCT_X
           GLOWEI_X = GLOWEI_X+ABS(ML(JJ))
     ENDIF
  ENDDO     ! end jj

  IF( (ABS(GLOWEI_X)).GT.1.0E-6) THEN
     FUNCT0=FUNCT0/GLOWEI_X
  ENDIF

  FUNCT= FUNCT+LAMDA*FUNCT0 
  ewrite(3,*) 'FUNCT,FUNCT0',FUNCT,FUNCT0

! Method Three
!-----------
      CASE(3)

  FUNCT0 = 0.0
  GLOWEI_XT = 0.0
  FUNCT_XT = 0.0
  ALLOCATE(ADWEI(NOSNOD))
  DO JJ=1,NONODS
     DO KK = 1,NOSTIM
        ADWEI(1:NOSNOD)=0.0
        GLOWEI_X = 0.0
        DO II=1,NOSNOD
! this is for 2D for the time being, we will fix it later on.... 
           ADWEI(II) = EXP(-((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2)/(2.*ISPAN_X**2))*ABS(ML(JJ))

           GLOWEI_XT = GLOWEI_XT+(EXP(-(ACCTIM-STIME(KK))**2/(2.*ISPAN_T**2)))* ADWEI(II)
           IF(D3) THEN
             FUNCT_X = 0.5*( (U(JJ)-UEXAC(II,KK))**2+(V(JJ)-VEXAC(II,KK))**2+       &
                             (W(JJ)-WEXAC(II,KK))**2 )
           ELSE
             FUNCT_X = 0.5*( (U(JJ)- UEXAC(II,KK))**2+(V(JJ)-VEXAC(II,KK))**2)
           ENDIF
           FUNCT_XT=FUNCT_XT+ GLOWEI_XT*FUNCT_X
        ENDDO  ! end ii
     ENDDO ! end kk
! to judge whether the nonods is been taked into account the list of nodes around the detectors
        IF( SQRT((X(JJ)-SX(II))**2+(Y(JJ)-SY(II))**2).LT.3*ISPAN_X.AND.         &
             SNOD_LIST(SNOD_LIST_NO).NE.JJ) THEN
           SNOD_LIST_NO=SNOD_LIST_NO+1
           SNOD_LIST(SNOD_LIST_NO)=JJ
        ENDIF

  ENDDO     ! end jj

  DEALLOCATE(ADWEI)
  IF( (ABS(GLOWEI_XT)).GT.1.0E-6) THEN
    FUNCT0 = FUNCT0+FUNCT_XT/(GLOWEI_XT)
  ENDIF

  FUNCT= FUNCT+ALPHA1*FUNCT0 

! Method Three
!--------------

      CASE DEFAULT

    END SELECT

  DEALLOCATE(ML)
     DEALLOCATE(UEXAC)  
     DEALLOCATE(VEXAC)  
     DEALLOCATE(WEXAC)  
     DEALLOCATE(SX)  
     DEALLOCATE(SY)  
     DEALLOCATE(SZ)  
     DEALLOCATE(STIME)  
  DEALLOCATE(SNOD_LIST)

END SUBROUTINE FUNCTION_POD


SUBROUTINE IN_OUTPUT_PG(IN,PG,PGNODS,TOTELE,MLOC,NGI,PGNDGLNO,M,MLX,MLY,MLZ,  &
    D3,ITSTIME,  &
    NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM)
  INTEGER ::PGNODS,TOTELE,MLOC,NGI,IN
  INTEGER ::PGNDGLNO(TOTELE*MLOC)
  REAL    ::PG(PGNODS)
  REAL M(MLOC,NGI),MLX(MLOC,NGI),MLY(MLOC,NGI),MLZ(MLOC,NGI)
  LOGICAL ::D3
  REAL    ::ACCTIM,DT
  INTEGER ::ITSTIME
  CHARACTER(40) DIRNAME
  CHARACTER(40) FILENAME_PG,FILENAME_MPG,FILENAME1,FILENAME2
  INTEGER I,K1,K2,K3,GI
  INTEGER NPGCOLM
  INTEGER PGFINDRM(PGNODS+1),PGCOLM(NPGCOLM)
  INTEGER PGCENTRM(PGNODS) 

!############
if(itstime.eq.1) then
  WRITE(FILENAME_PG, 10) ITSTIME
10 FORMAT('PG_data',I5.5)
  WRITE(FILENAME_MPG, 20) ITSTIME
20 FORMAT('MPG_data',I5.5)

!  ewrite(3,*)  FILENAME_PG
!  K1 = index( DIRNAME, ' ' ) - 1
  K2 = index( FILENAME_PG, ' ' ) - 1
  K3 = index( FILENAME_MPG, ' ' ) - 1
!  FILENAME=DIRNAME(1:K1) // FILENAME_PG(1:K2)
  FILENAME2=FILENAME_MPG(1:K3)
!  ewrite(3,*)  FILENAME2
  OPEN(2,FILE = FILENAME2)

  FILENAME1=FILENAME_PG(1:K2)
!  ewrite(3,*)  FILENAME1
  OPEN(1,FILE = FILENAME1)

IF(IN.EQ.1) THEN
  READ(2,*) PGNODS,TOTELE,MLOC,D3
  READ(2,*) NPGCOLM
  read(2,*) (PGFINDRM(I),I=1,PGNODS+1)
  read(2,*) (PGCOLM(I),I=1,NPGCOLM)
  read(2,*) (PGCENTRM(I),I=1,PGNODS)

!  DO GI=1,NGI
!     READ(1,*) (M(I,GI),I=1,MLOC)
!     READ(1,*) (MLX(I,GI),I=1,MLOC)
!     READ(1,*) (MLY(I,GI),I=1,MLOC)
!     IF(D3) THEN
!        READ(1,*) (MLZ(I,GI),I=1,MLOC)
!     ENDIF            
!  ENDDO
!  DO GI=1,NGI
!     READ(1,*) (M(GI,I),I=1,MLOC)
!     READ(1,*) (MLX(GI,I),I=1,MLOC)
!     READ(1,*) (MLY(GI,I),I=1,MLOC)
!     IF(D3) THEN
!        READ(1,*) (MLZ(GI,I),I=1,MLOC)
!     ENDIF            
!  ENDDO

  READ(1,*) (PGNDGLNO(I),I=1,TOTELE*MLOC)
  READ(1,*) (PG(I),I=1,PGNODS)

ELSE
        ewrite(3,*) 'wrting out...'
  WRITE(2,*) PGNODS,TOTELE,MLOC,D3
  WRITE(2,*) NPGCOLM
  WRITE(2,*) (PGFINDRM(I),I=1,PGNODS+1)
  WRITE(2,*) (PGCOLM(I),I=1,NPGCOLM)
  WRITE(2,*) (PGCENTRM(I),I=1,PGNODS)
  ewrite(3,*) '11'
!  DO GI=1,NGI
!     WRITE(1,*) (M(I,GI),I=1,MLOC)
!     WRITE(1,*) (MLX(I,GI),I=1,MLOC)
!     WRITE(1,*) (MLY(I,GI),I=1,MLOC)
!     IF(D3) THEN
!        WRITE(1,*) (MLZ(I,GI),I=1,MLOC)
!     ENDIF            
!  ENDDO

!  DO GI=1,NGI
!     WRITE(1,*) (M(GI,I),I=1,MLOC)
!     WRITE(1,*) (MLX(GI,I),I=1,MLOC)
!     WRITE(1,*) (MLY(GI,I),I=1,MLOC)
!     IF(D3) THEN
!        WRITE(1,*) (MLZ(GI,I),I=1,MLOC)
!     ENDIF            
!  ENDDO

  WRITE(1,*) (PGNDGLNO(I),I=1,TOTELE*MLOC)
  WRITE(1,*) (PG(I),I=1,PGNODS)

ENDIF
  CLOSE(2)
  CLOSE(1)

endif
!#########
END SUBROUTINE IN_OUTPUT_PG

end module reduced_model
