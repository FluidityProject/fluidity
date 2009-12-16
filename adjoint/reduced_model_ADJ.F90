#include "fdebug.h"



SUBROUTINE REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,X,Y,Z,                       &
     LEFTSVD,SMEAN,COEF,M,MLX,MLY,MLZ,                                    &
     N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
     NTIME,DT,D3,DCYL,GEOBAL,SCFACTH0,COEF_ADJ,           &
     obsv,iflagobs,istate,coef_multivar_obs,leftsvd_obs,smean_obs,SnapNDT_obs,IFWIND,  &
     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO,SMEAN_PG,leftsvd_PG_x,leftsvd_PG_y,  &
     NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM,NOCVA,G,IFWRITE,GLOITS,LINITS,IF_REMESH_OBS,IF_REMESH_VORTICITY)

      use FLDebug
      use Solvers
      use AllSorts
      use Shape_Transformations
      use Coordinates
      use coriolis_module
      use tr2d_module
      use detnlxr_module
  IMPLICIT NONE
  INTEGER, INTENT(IN)  ::NSVD_TOTAL,NSVD,NVAR,NSVD_U,NSVD_PHI
  INTEGER, INTENT(IN)  ::TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,MXNODS
  REAL, INTENT(IN),DIMENSION(NONODS) ::X,Y,Z
  REAL, INTENT(IN),DIMENSION(NONODS) :: MUPTXX,MUPTXY,MUPTXZ
  REAL, INTENT(IN),DIMENSION(NONODS) :: MUPTYY,MUPTYZ,MUPTZZ
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS,NSVD) :: LEFTSVD
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS) ::SMEAN
  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF
  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF_ADJ
  REAL, INTENT(IN)    ::DT
  INTEGER, INTENT(IN) ::NTIME
  LOGICAL, INTENT(IN) ::D3,DCYL,IF_REMESH_VORTICITY,IF_REMESH_OBS
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: N,NLX,NLY,NLZ
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: M,MLX,MLY,MLZ
  REAL, INTENT(IN),DIMENSION(NGI)  ::WEIGHT
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC)::PNDGLN
  ! coriolos force
  INTEGER, INTENT(IN)  :: GEOBAL
  REAL, INTENT(IN)  :: SCFACTH0
  INTEGER NOCVA
  REAL G(NOCVA)
  ! Wind stress
  !------------
  LOGICAL IFWIND

  ! observational
  integer istate
  integer SnapNDT_obs
  integer,INTENT(IN),DIMENSION(0:ntime,mxnods) ::iflagobs
  REAL,INTENT(IN),DIMENSION(0:ntime,istate):: obsv
  REAL coef_multivar_obs(0:ntime,NSVD_TOTAL)
  REAL leftsvd_obs(ISTATE,nsvd)
  REAl soux_mean,souy_mean,souz_mean,soux1,souy1,souz1,soux2,souy2,souz2
  REAL smean_obs(istate)
  REAL,DIMENSION(:),ALLOCATABLE   :: Bmat_obs

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
  INTEGER ::ISVD,JSVD,II,JJ,KK,NSTEP,ELE,GI,ILOC,JLOC,I,J,K
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
            U_xx_mean,U_yy_mean,U_zz_mean,V_xx_mean,V_yy_mean,V_zz_mean,W_xx_mean,W_yy_mean,W_zz_mean,  &
            UDmean_UADJx,VDmean_UADJy,WDmean_UADJz,UDmean_VADJx,VDmean_VADJy,WDmean_VADJz, &
            UDmean_WADJx,VDmean_WADJy,WDmean_WADJz,                                        &
            UD_UADJx,VD_UADJy,WD_UADJz,UD_VADJx,VD_VADJy,WD_VADJz,UD_WADJx,VD_WADJy,WD_WADJz, &
            Umeanx_UADJ,Umeany_VADJ,Umeanz_WADJ,Vmeanx_UADJ,Vmeany_VADJ,                   &
            Vmeanz_WADJ,Wmeanx_UADJ,Wmeany_VADJ,Wmeanz_WADJ,                               &         
            Ux_UADJ,Uy_VADJ,Uz_WADJ,Vx_UADJ,Vy_VADJ,Vz_WADJ,Wx_UADJ,Wy_VADJ,Wz_WADJ


  INTEGER ::NCOLM,NCOLM1,NCOLM2
  INTEGER ::IBL11,IBL12,IBL13,IBL14,IBL15
  INTEGER ::IBL21,IBL22,IBL23,IBL24,IBL25
  INTEGER ::IBL31,IBL32,IBL33,IBL34,IBL35
  INTEGER ::IBL41,IBL42,IBL43,IBL44,IBL45
  INTEGER ::IBL51,IBL52,IBL53,IBL54,IBL55
  real    ::VOL,ftest,VOLUME
  INTEGER ::NT,ROW,ROWCNT
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
  INTEGER  FINDRM(NSVD+1),COLM(NSVD*NSVD), CENTRM(NSVD)

  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM,VEC,VEC1,VEC_extr
  REAL,DIMENSION(:,:),ALLOCATABLE   :: ReMat,ReMat1
  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM0,BMat,BMat1,BMat0,BMat01,VEC0,VEC0_source,BMat0_source,BMat_extr
  REAL,DIMENSION(:),ALLOCATABLE     :: DETWEI,XD,YD,ZD
  REAL,DIMENSION(:),ALLOCATABLE     :: UD_mean,VD_mean,WD_mean
  REAL,DIMENSION(:,:),ALLOCATABLE   :: UD,VD,WD
  REAL,DIMENSION(:),ALLOCATABLE     :: UDx_mean,UDy_mean,UDz_mean,VDx_mean,VDy_mean,VDz_mean,WDx_mean,WDy_mean,WDz_mean
  REAL,DIMENSION(:,:),ALLOCATABLE   :: UDx,UDy,UDz,VDx,VDy,VDz,WDx,WDy,WDz
  REAL,DIMENSION(:),ALLOCATABLE     :: Umeanx,Umeany,Umeanz,Vmeanx,Vmeany,Vmeanz,Wmeanx,Wmeany,Wmeanz
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_u,WEI_v,WEI_w,WEI_p
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_ux,WEI_uy,WEI_uz,WEI_vx,WEI_vy,WEI_vz,WEI_wx,WEI_wy,WEI_wz
  REAL,DIMENSION(:,:),ALLOCATABLE   :: NX,NY,NZ,MX,MY,MZ
  REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF,NEWCOEF1
  REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF_ADJ,NEWCOEF1_ADJ
  REAL,DIMENSION(:),ALLOCATABLE     :: CORIOLIS

  REAL,DIMENSION(:,:),ALLOCATABLE   :: BIGM1,BIGM1consi
  REAL,DIMENSION(:),ALLOCATABLE     :: UDMat
!  REAL,DIMENSION(:),ALLOCATABLE     :: FINDRM,COLM,CENTRM
  REAL,DIMENSION(:),ALLOCATABLE     :: ML
  REAL,DIMENSION(:),ALLOCATABLE     :: SOURCX,SOURCY,SOURCZ
  REAL,DIMENSION(:),ALLOCATABLE     :: VORTICITY_U_mean,VORTICITY_V_mean
  REAL,DIMENSION(:,:),ALLOCATABLE   :: VORTICITY_U_leftsvd,VORTICITY_V_leftsvd

  INTEGER IOS
  LOGICAL IFWRITE
  INTEGER GLOITS,LINITS
  CHARACTER*240 NAME,NAME1,NAME2

!  LEFTSVD=leftsvd_obs

  ALLOCATE(BIGM(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(VEC(NSVD_TOTAL))
  ALLOCATE(ReMat( NSVD_TOTAL-NSVD,NSVD_TOTAL-NSVD))
  ALLOCATE(BIGM0(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BIGM1(NSVD_TOTAL*NSVD_TOTAL,NSVD_u))
  ALLOCATE(BIGM1consi(NSVD_TOTAL*NSVD_TOTAL,NSVD_u))
  ALLOCATE(VEC0(NSVD_TOTAL))
  ALLOCATE(VEC_extr(NSVD_TOTAL))
  ALLOCATE(BMat_extr(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(VEC1(NSVD_TOTAL-NSVD))
  ALLOCATE(BMat(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat0(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat01(NSVD_TOTAL*NSVD_TOTAL))
!  ALLOCATE(BMat1(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(DETWEI(NGI))
  ALLOCATE(UD(NGI,NSVD))
  ALLOCATE(VD(NGI,NSVD))
  ALLOCATE(WD(NGI,NSVD))
  ALLOCATE(UD_mean(NGI))
  ALLOCATE(VD_mean(NGI))
  ALLOCATE(WD_mean(NGI))
  ALLOCATE(UDx(NGI,NSVD))
  ALLOCATE(UDy(NGI,NSVD))
  ALLOCATE(UDz(NGI,NSVD))
  ALLOCATE(VDx(NGI,NSVD))
  ALLOCATE(VDy(NGI,NSVD))
  ALLOCATE(VDz(NGI,NSVD))
  ALLOCATE(WDx(NGI,NSVD))
  ALLOCATE(WDy(NGI,NSVD))
  ALLOCATE(WDz(NGI,NSVD))
  ALLOCATE(UDx_mean(NGI))
  ALLOCATE(UDy_mean(NGI))
  ALLOCATE(UDz_mean(NGI))
  ALLOCATE(VDx_mean(NGI))
  ALLOCATE(VDy_mean(NGI))
  ALLOCATE(VDz_mean(NGI))
  ALLOCATE(WDx_mean(NGI))
  ALLOCATE(WDy_mean(NGI))
  ALLOCATE(WDz_mean(NGI))

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
  ALLOCATE(NEWCOEF_ADJ(NSVD_TOTAL))
  ALLOCATE(NEWCOEF1_ADJ(NSVD_TOTAL-NSVD))
  ALLOCATE(CORIOLIS(NGI))

!  ALLOCATE(SOURCX(NSVD))
!  ALLOCATE(SOURCY(NSVD))
!  ALLOCATE(SOURCZ(NSVD))
  ALLOCATE(BMat_OBS(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(VEC0_source(NSVD_TOTAL))
  ALLOCATE(BMat0_source(NSVD_TOTAL*NSVD_TOTAL))

  ! for geographic pressure
  ALLOCATE(MPGX(MPGLOC,NGI))
  ALLOCATE(MPGY(MPGLOC,NGI))
  ALLOCATE(MPGZ(MPGLOC,NGI))
  ALLOCATE(WEI_PGx(NGI))
  ALLOCATE(WEI_PGy(NGI))
  ALLOCATE(WEI_PGz(NGI))
  ewrite(3,*) '111'
!  ALLOCATE(LEFTSVD_PG_x(PGNODS,NSVD_PHI))
!  ALLOCATE(LEFTSVD_PG_y(PGNODS,NSVD_PHI))
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
  BIGM1consi=0.0
  VEC0=0.0
  VEC1=0.0
  VEC_extr=0.0
  BMat_extr=0.0
  BMat=0.0
  BMat0=0.0
  BMat01=0.0
!  BMat1=0.0
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
  UDx=0.0
  UDy=0.0
  UDz=0.0              
  VDx=0.0
  VDy=0.0
  VDz=0.0              
  WDx=0.0
  WDy=0.0
  WDz=0.0
  UDx_mean=0.0
  UDy_mean=0.0
  UDz_mean=0.0              
  VDx_mean=0.0
  VDy_mean=0.0
  VDz_mean=0.0              
  WDx_mean=0.0
  WDy_mean=0.0
  WDz_mean=0.0

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
  NEWCOEF_ADJ=0.0
  NEWCOEF1_ADJ=0.0
  CORIOLIS=0.0


  NX=0.0
  NY=0.0
  NZ=0.0
  MX=0.0
  MY=0.0
  MZ=0.0

!  SOURCX=0.0
!  SOURCY=0.0
!  SOURCZ=0.0

   VEC0_source=0.0
   BMat0_source=0.0
   BMat_OBS =0.0
  R(1:7*NSVD)=0.0
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
!  LEFTSVD_PG_x=0.0
!  LEFTSVD_PG_y=0.0
!  SMEAN_PG=0.0

!****************
! VORTICITY
!****************
IF(IF_REMESH_VORTICITY) THEN
  allocate(VORTICITY_U_mean(nonods))
  allocate(VORTICITY_V_mean(nonods))
  allocate(VORTICITY_U_leftsvd(nonods,nsvd))
  allocate(VORTICITY_V_leftsvd(nonods,nsvd))
  VORTICITY_U_mean=0.0
  VORTICITY_V_mean=0.0
  VORTICITY_U_leftsvd=0.0
  VORTICITY_V_leftsvd=0.0
  call Vorticity2D(X,Y,Z, &
     SMEAN(1:nonods), SMEAN(mxnods+1:mxnods+nonods), &
     ndglno, totele, nonods, nloc, ngi, &
     D3,DCYL, & 
     VORTICITY_U_mean, VORTICITY_V_mean, &
     N,NLX,NLY,NLZ, &
     WEIGHT)
 DO ii=1,nsvd
  call Vorticity2D(X,Y,Z, &
     LEFTSVD(1:nonods,ii), LEFTSVD(mxnods+1:mxnods+nonods,ii), &
     ndglno, totele, nonods, nloc, ngi, &
     D3,DCYL, & 
     VORTICITY_U_leftsvd(:,ii), VORTICITY_V_leftsvd(:,ii), &
     N,NLX,NLY,NLZ, &
     WEIGHT)
 ENDDO
ENDIF
!*******************






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
!              open(1,file='LEFTSVD_PG_x.dat')
!              open(2,file='LEFTSVD_PG_y.dat')
!              DO j=1,NSVD
!                read(1,*) (LEFTSVD_PG_x(I,j),I=1,PGNODS)
!                read(2,*) (LEFTSVD_PG_y(I,j),I=1,PGNODS)
!              ENDDO
!              close(1)
!              close(2)

!              open(1,file='SMEAN_PG.dat')
!                read(1,*) (SMEAN_PG(i),I=1,PGNODS)
!              close(1)

IF(.false.) then

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
!!CORIOLIS=0.0

!           ewrite(3,*) 'MPGX',MPGX
!           ewrite(3,*) 'MPGY',MPGY
!           ewrite(3,*) 'CORIOLIS',CORIOLIS
!           ewrite(3,*) 'DETWEI',DETWEI

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

  CALL SOLCG(LEFTSVD_PG_x(:,ISVD),PGVEC_x,PGNODS,PGNODS,PGNODS,.TRUE.,  &
       PGMATRIX_x,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM, &
       halotagT10,KITS)

  CALL SOLCG(LEFTSVD_PG_y(:,ISVD),PGVEC_y,PGNODS,PGNODS,PGNODS,.TRUE.,  &
       PGMATRIX_y,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM,&
       halotagT10,KITS)

IF(ISVD.EQ.1) THEN         
  CALL SOLCG(SMEAN_PG,PGVEC_mean,PGNODS,PGNODS,PGNODS,.TRUE.,  &
       PGMATRIX_mean,PGFINDRM,PGCOLM,NPGCOLM, NPGCOLM,&
       halotagT10,KITS)
ENDIF

    DEALLOCATE(WORK1)
    DEALLOCATE(WORK2)
    DEALLOCATE(WORK3)
    DEALLOCATE(WORK4)
    DEALLOCATE(WORK5)
    DEALLOCATE(WORKP)

300              CONTINUE

ENDIF !end if re-calculating the geostrophic pressure or passing from forward model

!LEFTSVD_PG_x=0.0
!LEFTSVD_PG_y=0.0


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

! only for 2D!!!!
  DO II=2*NSVD_U+1,3*NSVD_U
     COEF_ADJ(0,II)=0.0
  ENDDO

!     LEFTSVD(2*nonods+1:3*nonods,:)=0.0
     SMEAN(2*mxnods+1:3*mxnods)=0.0

! for adjoint model
     NEWCOEF_ADJ(:)=0.0

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
IF(IFWRITE.AND.LINITS.EQ.1) then
!  DO 100 NSTEP=1,ntime
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
           ewrite(3,*) 'ele=',ele
             ewrite(3,*)  TOTELE,NONODS,XNONOD,NLOC,NGI,MLOC
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
              UD_mean(GI)=0.0
              VD_mean(GI)=0.0
              WD_mean(GI)=0.0
              UD(GI,:)=0.0
              VD(GI,:)=0.0
              WD(GI,:)=0.0
!              XD(GI)=0.0
!              YD(GI)=0.0
!              ZD(GI)=0.0
! for the nonlinear term
              UDx(GI,:)=0.0
              UDy(GI,:)=0.0
              UDz(GI,:)=0.0              
              VDx(GI,:)=0.0
              VDy(GI,:)=0.0
              VDz(GI,:)=0.0              
              WDx(GI,:)=0.0
              WDy(GI,:)=0.0
              WDz(GI,:)=0.0
              UDx_mean(GI)=0.0
              UDy_mean(GI)=0.0
              UDz_mean(GI)=0.0              
              VDx_mean(GI)=0.0
              VDy_mean(GI)=0.0
              VDz_mean(GI)=0.0              
              WDx_mean(GI)=0.0
              WDy_mean(GI)=0.0
              WDz_mean(GI)=0.0

              CORIOLIS(GI)=0.0

              DO JLOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+JLOC)

                 UD_mean(GI)=UD_mean(GI)+SMEAN(JGL)*N(JLOC,GI)
                 VD_mean(GI)=VD_mean(GI)+SMEAN(1*MXNODS+JGL)*N(JLOC,GI)
                 WD_mean(GI)=WD_mean(GI)+SMEAN(2*MXNODS+JGL)*N(JLOC,GI)

                 UDx_mean(GI)=UDx_mean(GI)+SMEAN(JGL)*NX(JLOC,GI)
                 UDy_mean(GI)=UDy_mean(GI)+SMEAN(JGL)*NY(JLOC,GI)
                 UDz_mean(GI)=UDz_mean(GI)+SMEAN(JGL)*NZ(JLOC,GI)

                 VDx_mean(GI)=VDx_mean(GI)+SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                 VDy_mean(GI)=VDy_mean(GI)+SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 VDz_mean(GI)=VDz_mean(GI)+SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)

                 WDx_mean(GI)=WDx_mean(GI)+SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                 WDy_mean(GI)=WDy_mean(GI)+SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                 WDz_mean(GI)=WDz_mean(GI)+SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)

!                 XD(GI)=XD(GI) + N(JLOC,GI)*X(JGL)
!                 YD(GI)=YD(GI) + N(JLOC,GI)*Y(JGL)
!                 ZD(GI)=ZD(GI) + N(JLOC,GI)*Z(JGL)

                 DO JSVD = 1,NSVD
                    UD(GI,JSVD)=UD(GI,JSVD)+LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    VD(GI,JSVD)=VD(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    WD(GI,JSVD)=WD(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)

                    UDx(GI,JSVD)=UDx(GI,JSVD)+LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    UDy(GI,JSVD)=UDy(GI,JSVD)+LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    UDz(GI,JSVD)=UDz(GI,JSVD)+LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)

                    VDx(GI,JSVD)=VDx(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VDy(GI,JSVD)=VDy(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    VDz(GI,JSVD)=VDz(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)

                    WDx(GI,JSVD)=WDx(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    WDy(GI,JSVD)=WDy(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WDz(GI,JSVD)=WDz(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                 ENDDO
              ENDDO
                 CORIOLIS(GI) = 2.0*FUNOME(XD(GI),YD(GI),ZD(GI))
           ENDDO

! Source terms which represet the misfit between the numerical results and observations

           soux_mean=0.0
           souy_mean=0.0
           souz_mean=0.0

        IF(IF_REMESH_VORTICITY) THEN
           DO JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
              DO GI=1,NGI
                 soux_mean=soux_mean + WEI_u(GI)*DETWEI(GI)*VORTICITY_U_mean(JGL)*N(JLOC,GI)
                 souy_mean=souy_mean + WEI_v(GI)*DETWEI(GI)*VORTICITY_V_mean(JGL)*N(JLOC,GI)
              ENDDO
           ENDDO
        ELSE
           DO JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
              DO GI=1,NGI
                 soux_mean=soux_mean + WEI_u(GI)*DETWEI(GI)*(SMEAN(JGL)-SMEAN_OBS(JGL))*N(JLOC,GI)
                 souy_mean=souy_mean + WEI_v(GI)*DETWEI(GI)*(SMEAN(1*MXNODS+JGL)-SMEAN_OBS(1*MXNODS+JGL))*N(JLOC,GI)
                 souz_mean=souz_mean + WEI_w(GI)*DETWEI(GI)*(SMEAN(2*MXNODS+JGL)-SMEAN_OBS(2*MXNODS+JGL))*N(JLOC,GI)
              ENDDO
           ENDDO
        ENDIF

           VEC0_source(ISVD)= VEC0_source(ISVD)+DT*soux_mean
           VEC0_source(1*NSVD+ISVD)= VEC0_source(1*NSVD+ISVD)+DT*souy_mean
           VEC0_source(2*NSVD+ISVD)= VEC0_source(2*NSVD+ISVD)+DT*souz_mean
! end the source terms


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
                 
                 Ux_UADJ=0.0
                 Uy_VADJ=0.0
                 Uz_WADJ=0.0
                 Vx_UADJ=0.0
                 Vy_VADJ=0.0
                 Vz_WADJ=0.0
                 Wx_UADJ=0.0
                 Wy_VADJ=0.0
                 Wz_WADJ=0.0
                 
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
                 ! source terms for misfit
                 SOUX1 =0.0
                 SOUY1 =0.0
                 SOUZ1 =0.0
                 SOUX2 =0.0
                 SOUY2 =0.0
                 SOUZ2 =0.0

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
                    ! source terms for misfit

                  IF(IF_REMESH_VORTICITY) THEN
                    SOUX1 = SOUX1 + WEI_u(GI)*DETWEI(GI)*VORTICITY_U_leftsvd(JGL,JSVD)*N(JLOC,GI)
                    SOUY1 = SOUY1 + WEI_v(GI)*DETWEI(GI)*VORTICITY_V_leftsvd(JGL,JSVD)*N(JLOC,GI)

                  ELSE
                    SOUX1 = SOUX1 + WEI_u(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    SOUY1 = SOUY1 + WEI_v(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
!                    SOUZ1 = SOUZ1 + WEI_w(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    
                    SOUX2 = SOUX2 + WEI_u(GI)*DETWEI(GI)*LEFTSVD_OBS(JGL,JSVD)*N(JLOC,GI)
                    SOUY2 = SOUY2 + WEI_v(GI)*DETWEI(GI)*LEFTSVD_OBS(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
!                    SOUZ2 = SOUZ2 +WEI_w(GI)*DETWEI(GI)*LEFTSVD_OBS(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
                  ENDIF

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


                      !BIGM1consi(IBL11+COUNT,JJ)=BIGM1consi(IBL11+COUNT,JJ)+Ux_UADJ
                       !--------------------------------------------------------------------------------
                       BIGM1consi(IBL11+COUNT,JJ)=BIGM1consi(IBL11+COUNT,JJ)+   &
                                WEI_u(GI)*UDx(GI,JJ)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)                          

                       !BIGM1consi(IBL12+COUNT,JJ)=BIGM1consi(IBL12+COUNT,JJ)+Uy_VADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL12+COUNT,JJ)=BIGM1consi(IBL12+COUNT,JJ)+   &
                          WEI_u(GI)*UDy(GI,JJ)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL13+COUNT,JJ)=BIGM1consi(IBL13+COUNT,JJ)+Uz_WADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL13+COUNT,JJ)=BIGM1consi(IBL13+COUNT,JJ)+   &
                          WEI_u(GI)*UDz(GI,JJ)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
!                      
!
                       !BIGM1consi(IBL21+COUNT,JJ)=BIGM1consi(IBL21+COUNT,JJ)+Vx_UADJ
                       !----------------------------------------------------
                       BIGM1consi(IBL21+COUNT,JJ)=BIGM1consi(IBL21+COUNT,JJ)+   &
                          WEI_v(GI)*VDx(GI,JJ)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL22+COUNT,JJ)=BIGM1consi(IBL22+COUNT,JJ)+Vy_VADJ)
                       !--------------------------------------------------------------------------------
                       BIGM1consi(IBL22+COUNT,JJ)=BIGM1consi(IBL22+COUNT,JJ)+   & 
                          WEI_v(GI)*VDy(GI,JJ)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL23+COUNT,JJ)=BIGM1consi(IBL23+COUNT,JJ)+Vz_WADJ
                       !----------------------------------------------------
                       BIGM1consi(IBL23+COUNT,JJ)=BIGM1consi(IBL23+COUNT,JJ)+  &
                          WEI_v(GI)*VDz(GI,JJ)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
!
!
                       !BIGM1consi(IBL31+COUNT,JJ)=BIGM1consi(IBL31+COUNT,JJ)+Wx_UADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL31+COUNT,JJ)=BIGM1consi(IBL31+COUNT,JJ)+  &
                          WEI_w(GI)*WDx(GI,JJ)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL32+COUNT,JJ)=BIGM1consi(IBL32+COUNT,JJ)+Wy_VADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL32+COUNT,JJ)=BIGM1consi(IBL32+COUNT,JJ)+  &
                          WEI_w(GI)*WDz(GI,JJ)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL33+COUNT,JJ)=BIGM1consi(IBL33+COUNT,JJ)+Wz_WADJ
                       !--------------------------------------------------------------------------------
                       BIGM1consi(IBL33+COUNT,JJ)=BIGM1consi(IBL33+COUNT,JJ)+  &
                          WEI_w(GI)*WDz(GI,JJ)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)

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

              ! source terms for misfit
              BMat0_source(IBL11+COUNT)=BMat0_source(IBL11+COUNT)+DT*SOUX1
              BMat0_source(IBL22+COUNT)=BMat0_source(IBL22+COUNT)+DT*SOUY1
              BMat0_source(IBL33+COUNT)=BMat0_source(IBL33+COUNT)+DT*SOUZ1
!
              BMat_obs(IBL11+COUNT)=BMat_obs(IBL11+COUNT)+DT*SOUX2
              BMat_obs(IBL22+COUNT)=BMat_obs(IBL22+COUNT)+DT*SOUY2
              BMat_obs(IBL33+COUNT)=BMat_obs(IBL33+COUNT)+DT*SOUZ2


40            CONTINUE

20            CONTINUE

30            CONTINUE

!400 CONTINUE
        write(name,'(a,i0)') 'PODBIGM0b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BIGM0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBMat0b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BMat0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODVEC0b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (VEC0(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBIGM1b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            write(1,*) (BIGM1(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)


        write(name,'(a,i0)') 'vec0_source.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (VEC0_source(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'BMat0_source.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BMat0_source(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'BMat_obs.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BMat_obs(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBIGM1consi.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            write(1,*) (BIGM1consi(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)

ELSE  !NOT IF(IFWRITE.AND.LINITS.EQ.1)
      ewrite(3,*) 'read the submatrices for running the POD adjoint model'
        write(name,'(a,i0)') 'PODBIGM0b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (BIGM0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBMat0b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (BMat0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODVEC0b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (VEC0(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBIGM1b.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            read(1,*) (BIGM1(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)

        write(name,'(a,i0)') 'vec0_source.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (VEC0_source(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'BMat0_source.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (BMat0_source(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'BMat_obs.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                read(1,*) (BMat_obs(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBIGM1consi.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            read(1,*) (BIGM1consi(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)

ENDIF  !END IF(IFWRITE.AND.LINITS.EQ.1) then
      ewrite(3,*) 'checkout the submatrices for running the POD adjoint model'
        write(name,'(a,i0)') 'PODBIGM0b.dat.',GLOITS

        write(name,'(a,i0)') 'PODBIGM1consi-.dat.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            write(1,*) (BIGM1consi(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)

        write(name,'(a,i0)') 'PODBIGM0.dat-.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BIGM0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBMat0.dat-.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BMat0(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'PODBIGM1.dat-.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
         DO JJ=1,NSVD_U
            write(1,*) (BIGM1(II,JJ),II=1,NSVD_TOTAL*NSVD_TOTAL)
         ENDDO
        close(1)

        write(name,'(a,i0)') 'PODVEC0b.dat-.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (VEC0(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'vec0_source.dat-.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (VEC0_source(II),II=1,NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'BMat0_source.dat-.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BMat0_source(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)

        write(name,'(a,i0)') 'BMat_obs.dat-.',GLOITS
        OPEN(1,FILE=trim(NAME),STATUS='UNKNOWN',IOSTAT=IOS)
                write(1,*) (BMat_obs(II),II=1,NSVD_TOTAL*NSVD_TOTAL)
        close(1)


  open(1,file='coef_multivar2.dat')
  do jj=0,ntime
  write(1,*) (coef(jj,ii),ii=1,nsvd_total)
  enddo
  close(1)

  open(1,file='coef_multivar_obs2.dat')
  do jj=0,ntime
  write(1,*) (coef_multivar_obs(jj,ii),ii=1,nsvd_total)
  enddo
  close(1)

  open(1,file='leftsvd2.dat')
  do jj=1,istate
  write(1,*) (leftsvd(jj,ii),ii=1,nsvd)
  enddo
  close(1)

  open(1,file='leftsvd_obs2.dat')
  do jj=1,istate
  write(1,*) (leftsvd_obs(jj,ii),ii=1,nsvd)
  enddo
  close(1)


   DO 100 NSTEP=1,ntime
   ewrite(3,*) 'nstep=',nstep
!      CALL Source_POD(NONODS,NTIME*DT,NSVD,DT,NTIME+1,NSTEP+1,     &
!                    SOURCX,SOURCY,SOURCZ,                          &
!                    SMEAN,LEFTSVD,COEF(:,NSTEP),                   &
!                    TOTELE,NLOC,NGI,                               &
!                    NDGLNO,XNONOD,xondgl,X,Y,Z,                    &
!                    N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)

      BIGM(:)=BIGM0(:)
      BMat(:)=BMat0(:)
      VEC(:)=VEC0(:)
!      if (SnapNDT_obs*int(nstep/SnapNDT_obs) .eq. nstep.or.nstep.eq.1) then  
!      if (SnapNDT_obs*int(nstep/SnapNDT_obs)) then  
!if(SnapNDT_obs*int(nstep/SnapNDT_obs).eq. nstep.and.nstep.le.int(ntime/2)) then

    IF(IF_REMESH_OBS) THEN
      if (SnapNDT_obs*int(nstep/SnapNDT_obs).eq. nstep.and.nstep.gt.1) then  
         VEC(:)=VEC(:)+VEC0_source(:)
      endif
    ELSE
         VEC(:)=VEC(:)+VEC0_source(:)
    ENDIF

      DO ISVD=1,NSVD
         ewrite(3,*) 'ISVD=',ISVD
         DO JSVD=1,NSVD
           ewrite(3,*) 'JSVD=',ISVD
            COUNT=(ISVD-1)*NSVD+JSVD
            COUNT1=(ISVD-1)*NSVD_PHI+JSVD
            DO JJ=1,NSVD
               BIGM(IBL11+COUNT)=BIGM(IBL11+COUNT)+DT*THETA*( BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ))

               BIGM(IBL22+COUNT)=BIGM(IBL22+COUNT)+DT*THETA*( BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ))

               BIGM(IBL33+COUNT)=BIGM(IBL33+COUNT)+DT*THETA*( BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ))



               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)-DT*ONEMINTH*( BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)  )


               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)-DT*ONEMINTH*( BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)  )

               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)-DT*ONEMINTH*( BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ) )
!        ewrite(3,*) 'cheking the results--marices 111'

              ! extra term for nonlinear            
              if(.true.) then
               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)+  &
                    DT*THETA*BIGM1consi(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+          &
                    DT*ONEMINTH*BIGM1consi(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL12+COUNT)=BMat(IBL12+COUNT)+DT*THETA*BIGM1consi(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1consi(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL13+COUNT)=BMat(IBL13+COUNT)+DT*THETA*BIGM1consi(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1consi(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)


               BMat(IBL21+COUNT)=BMat(IBL21+COUNT)+DT*THETA*BIGM1consi(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1consi(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)+  &
                    DT*THETA*BIGM1consi(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+  &
                    DT*ONEMINTH*BIGM1consi(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL23+COUNT)=BMat(IBL23+COUNT)+DT*THETA*BIGM1consi(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1consi(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)


               BMat(IBL31+COUNT)=BMat(IBL31+COUNT)+DT*THETA*BIGM1consi(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1consi(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL32+COUNT)=BMat(IBL32+COUNT)+DT*THETA*BIGM1consi(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+               & 
                    DT*ONEMINTH*BIGM1consi(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)+    &
                    DT*THETA*BIGM1consi(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+   &
                    DT*ONEMINTH*BIGM1consi(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
              else
               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)+  &
                    DT*THETA*BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+          &
                    DT*ONEMINTH*BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL12+COUNT)=BMat(IBL12+COUNT)+DT*THETA*BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL13+COUNT)=BMat(IBL13+COUNT)+DT*THETA*BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)


               BMat(IBL21+COUNT)=BMat(IBL21+COUNT)+DT*THETA*BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)+  &
                    DT*THETA*BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+  &
                    DT*ONEMINTH*BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL23+COUNT)=BMat(IBL23+COUNT)+DT*THETA*BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)


               BMat(IBL31+COUNT)=BMat(IBL31+COUNT)+DT*THETA*BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL32+COUNT)=BMat(IBL32+COUNT)+DT*THETA*BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+               & 
                    DT*ONEMINTH*BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)+    &
                    DT*THETA*BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+   &
                    DT*ONEMINTH*BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
              endif

            ENDDO
               VEC(ISVD)=VEC(ISVD)+BMat(IBL11+COUNT)*COEF_ADJ(NSTEP-1,JSVD)+                       &
                    BMat(IBL12+COUNT)*COEF_ADJ(NSTEP-1,1*NSVD+JSVD)+BMat(IBL13+COUNT)*COEF_ADJ(NSTEP-1,2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(ISVD)=VEC(ISVD)+BMat(IBL14+COUNT1)*COEF_ADJ(NSTEP-1,3*NSVD+JSVD)
               ENDIF

               VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat(IBL21+COUNT)*COEF_ADJ(NSTEP-1,JSVD)+         &
                    BMat(IBL22+COUNT)*COEF_ADJ(NSTEP-1,1*NSVD+JSVD)+BMat(IBL23+COUNT)*COEF_ADJ(NSTEP-1,2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat(IBL24+COUNT1)*COEF_ADJ(NSTEP-1,3*NSVD+JSVD)
               ENDIF

               VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat(IBL31+COUNT)*COEF_ADJ(NSTEP-1,JSVD)+         &
                    BMat(IBL32+COUNT)*COEF_ADJ(NSTEP-1,1*NSVD+JSVD)+BMat(IBL33+COUNT1)*COEF_ADJ(NSTEP-1,2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat(IBL34+COUNT1)*COEF_ADJ(NSTEP-1,3*NSVD+JSVD)
               ENDIF
        ewrite(3,*) 'cheking the results--marices 222'

            ! for source terms
!          if (SnapNDT_obs*int(nstep/SnapNDT_obs) .eq. nstep.or.nstep.eq.1) then  
!          if (SnapNDT_obs*int(nstep/SnapNDT_obs) .eq. nstep) then  
!          if (SnapNDT_obs*int(nstep/SnapNDT_obs).eq. nstep.and.nstep.le.int(ntime/2)) then  
    IF(IF_REMESH_OBS) THEN
          if (SnapNDT_obs*int(nstep/SnapNDT_obs).eq. nstep.and.nstep.gt.1) then  
            VEC(ISVD)=VEC(ISVD)+BMat0_source(IBL11+COUNT)*coef(NTIME-NSTEP,JSVD) &
                   -BMat_OBS(IBL11+COUNT)*coef_multivar_obs(NTIME-NSTEP,JSVD)
            VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat0_source(IBL22+COUNT)*coef(NTIME-NSTEP,1*NSVD+JSVD)  &
                   -BMat_OBS(IBL22+COUNT)*coef_multivar_obs(NTIME-NSTEP,1*NSVD+JSVD)
            VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat0_source(IBL33+COUNT)*coef(NTIME-NSTEP,2*NSVD+JSVD)  &
                   -BMat_OBS(IBL33+COUNT)*coef_multivar_obs(NTIME-NSTEP,2*NSVD+JSVD)
          endif
    ENDIF
         ENDDO
      ENDDO

      if(nstep.eq.1) then                
              open(1,file='Bmat41.dat')
              open(2,file='vec41.dat')
              write(1,*) Bmat
              write(2,*) vec
              close(1)
              close(2)
      endif

      if(nstep.eq.2) then                
              open(1,file='Bmat42.dat')
              open(2,file='vec42.dat')
              write(1,*) Bmat
              write(2,*) vec
              close(1)
              close(2)
      endif
      if(nstep.eq.3) then                
              open(1,file='Bmat43.dat')
              open(2,file='vec43.dat')
              write(1,*) Bmat
              write(2,*) vec
              close(1)
              close(2)
      endif
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

              ReMat = 0.0
              VEC1 = 0.0

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
!                       ftest=ftest+ReMat((kk-1)*NVAR*NSVD+ii,jj)
                       ELSEIF((kk.eq.4).AND.(ii.LE.NSVD_PHI)) then
                          ReMat( (kk-2)*NSVD+ii,jj)=BIGM( (kk-1)*NSVD*NSVD_TOTAL+(ii-1)*NSVD+jj)
                          ReMat( (kk-2)*NSVD+ii,1*NSVD+jj)= BIGM( (kk-1)*NSVD*NSVD_TOTAL+1*NSVD_PHI*NSVD+(ii-1)*NSVD+jj)
                          !     ReMat( (kk-1)*NVAR*NSVD+ii,2*NSVD+jj)= BIGM( (kk-1)*NVAR*NSVD*NSVD+2*NSVD*NSVD+(ii-1)*NSVD+jj)
                          IF(jj.LE.NSVD_PHI) THEN
                             ReMat( (kk-2)*NSVD+ii,2*NSVD+jj)= BIGM( IBL44+(ii-1)*NSVD_PHI+jj)
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
!        ewrite(3,*) 'cheking the results--marices 333'

           if(nstep.eq.1) then
              open(1,file='Remat1.dat')
              open(2,file='vec1.dat')
              write(1,*) Remat
              write(2,*) vec
              close(1)
              close(2)
           else if(nstep.eq.2) then
              open(1,file='Remat2.dat')
              open(2,file='vec2.dat')
              write(1,*) Remat
              write(2,*) vec
              close(1)
              close(2)
           endif
              open(1,file='Remat.dat')
              open(2,file='vec.dat')
              write(1,*) Remat
              write(2,*) vec
              close(1)
              close(2)

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
           write(10, *) (NEWCOEF1_ADJ(ii), ii=1,(nvar-1)*nsvd)
        close(10) 


              CALL GMRES(NEWCOEF_ADJ,VEC,NSVD,NSVD,NSVD,.FALSE.,                        & 
                   BIGM,FINDRM,COLM,NSVD*NSVD,NSVD*NSVD,&
                   halo_tag,KITS)



              DO II=1,NSVD_TOTAL
                 IF(II.LE.2*NSVD) THEN
                    COEF_ADJ(NSTEP,II)=NEWCOEF1_ADJ(II)
                 ELSEIF(II.GT.3*NSVD) THEN
                    COEF_ADJ(NSTEP,II)=NEWCOEF1_ADJ(II-NSVD)
                 ELSE
                    COEF_ADJ(NSTEP,II)=0.0
                 ENDIF
                 ewrite(3,*) 'COEF_ADJ(NSTEP,II)',COEF_ADJ(NSTEP,II)
              ENDDO

              ewrite(3,*) 'after NSTEP',NSTEP
100           CONTINUE


              open(10, file ='COEF_ADJ.dat')
              do ii=0,ntime
                 write(10, *) (COEF_ADJ(ii,jj), jj=1,nsvd)
              enddo
              close(10)
              ewrite(3,*) 'exit REDUCED_MODEL_DIFFCOEF'

! Calculate the gradient for the intial inversion
!-------------------------------------------------
      DO ISVD=1,NSVD
         DO JSVD=1,NSVD
           IF(ISVD.EQ.JSVD) THEN
             COUNT=(ISVD-1)*NSVD+JSVD
             COUNT1=(ISVD-1)*NSVD_PHI+JSVD
            DO JJ=1,NSVD
              
              if(.true.) then
               BMat_extr(IBL11+COUNT)=BMat_extr(IBL11+COUNT)+ &
                    DT*ONEMINTH*BIGM1consi(IBL11+COUNT,JJ)*COEF(0,JJ)+BMat(IBL11+COUNT)

               BMat_extr(IBL22+COUNT)=BMat_extr(IBL22+COUNT)+  &
                    DT*ONEMINTH*BIGM1consi(IBL22+COUNT,JJ)*COEF(0,1*NSVD+JJ)+BMat(IBL22+COUNT)

               BMat_extr(IBL33+COUNT)=BMat_extr(IBL33+COUNT)+    &
                    DT*ONEMINTH*BIGM1consi(IBL33+COUNT,JJ)*COEF(0,2*NSVD+JJ)+BMat(IBL33+COUNT)
              else
               BMat_extr(IBL11+COUNT)=BMat_extr(IBL11+COUNT)+  &
                    DT*ONEMINTH*BIGM1(IBL11+COUNT,JJ)*COEF(0,JJ)+BMat(IBL11+COUNT)

               BMat_extr(IBL22+COUNT)=BMat_extr(IBL22+COUNT)+  &
                    DT*ONEMINTH*BIGM1(IBL22+COUNT,JJ)*COEF(0,1*NSVD+JJ)+BMat(IBL22+COUNT)


               BMat_extr(IBL33+COUNT)=BMat_extr(IBL33+COUNT)+    &
                    DT*ONEMINTH*BIGM1(IBL33+COUNT,JJ)*COEF(0,2*NSVD+JJ)+BMat(IBL33+COUNT)
              endif
            ENDDO
            VEC_extr(ISVD)=BMat_extr(IBL11+COUNT)
            VEC_extr(1*NSVD+ISVD)=BMat_extr(IBL22+COUNT)
            VEC_extr(2*NSVD+ISVD)=BMat_extr(IBL33+COUNT)
           ENDIF
         ENDDO
      ENDDO

      DO II=1,3*nsvd
        G(ii) = coef_ADJ(NTIME-1,ii)*VEC_extr(ii)
      ENDDO


!-------------------------------!
! DEALLOCATE THE LOCAL MATRICES !
!-------------------------------!

  DEALLOCATE(BIGM)
  DEALLOCATE(VEC)
  DEALLOCATE(ReMat)
  DEALLOCATE(BIGM0)
  DEALLOCATE(BIGM1)
  DEALLOCATE(BIGM1consi)
  DEALLOCATE(VEC0)
  DEALLOCATE(VEC1)
  DEALLOCATE(VEC_extr)
  DEALLOCATE(BMat_extr)
  DEALLOCATE(BMat)
  DEALLOCATE(BMat0)
  DEALLOCATE(BMat01)
  DEALLOCATE(DETWEI)
  DEALLOCATE(UD)
  DEALLOCATE(VD)
  DEALLOCATE(WD)
  DEALLOCATE(UD_mean)
  DEALLOCATE(VD_mean)
  DEALLOCATE(WD_mean)
  DEALLOCATE(UDx)
  DEALLOCATE(UDy)
  DEALLOCATE(UDz)
  DEALLOCATE(VDx)
  DEALLOCATE(VDy)
  DEALLOCATE(VDz)
  DEALLOCATE(WDx)
  DEALLOCATE(WDy)
  DEALLOCATE(WDz)

  DEALLOCATE(UDx_mean)
  DEALLOCATE(UDy_mean)
  DEALLOCATE(UDz_mean)
  DEALLOCATE(VDx_mean)
  DEALLOCATE(VDy_mean)
  DEALLOCATE(VDz_mean)
  DEALLOCATE(WDx_mean)
  DEALLOCATE(WDy_mean)
  DEALLOCATE(WDz_mean)

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
  DEALLOCATE(NEWCOEF_ADJ)
  DEALLOCATE(NEWCOEF1_ADJ)
  DEALLOCATE(CORIOLIS)

!DEALLOCATE(SOURCX)
!DEALLOCATE(SOURCY)
!DEALLOCATE(SOURCZ)

  DEALLOCATE(VEC0_source)
  DEALLOCATE(BMat_OBS)
  DEALLOCATE(BMat0_source)
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
IF(IF_REMESH_VORTICITY) THEN
  deallocate(VORTICITY_U_mean)
  deallocate(VORTICITY_V_mean)
  deallocate(VORTICITY_U_leftsvd)
  deallocate(VORTICITY_V_leftsvd)
ENDIF

  RETURN
END SUBROUTINE REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ

SUBROUTINE REDUCED_MODEL_DIFFCOEF1consistent_ADJ(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,X,Y,Z,                       &
     LEFTSVD,SMEAN,COEF,M,MLX,MLY,MLZ,                                    &
     N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
     NTIME,DT,D3,DCYL,GEOBAL,SCFACTH0,COEF_ADJ,           &
     obsv,iflagobs,istate,coef_multivar_obs,leftsvd_obs,smean_obs,SnapNDT_obs,IFWIND,NOCVA,G)

      use FLDebug
      use Solvers
      use AllSorts
      use Shape_Transformations
      use Coordinates
      use coriolis_module
  IMPLICIT NONE
  INTEGER, INTENT(IN)  ::NSVD_TOTAL,NSVD,NVAR,NSVD_U,NSVD_PHI
  INTEGER, INTENT(IN)  ::TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,MXNODS
  REAL, INTENT(IN),DIMENSION(NONODS) ::X,Y,Z
  REAL, INTENT(IN),DIMENSION(NONODS) :: MUPTXX,MUPTXY,MUPTXZ
  REAL, INTENT(IN),DIMENSION(NONODS) :: MUPTYY,MUPTYZ,MUPTZZ
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS,NSVD) :: LEFTSVD
  REAL, INTENT(INOUT),DIMENSION(NVAR*MXNODS) ::SMEAN
  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF
  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF_ADJ
  REAL, INTENT(IN)    ::DT
  INTEGER, INTENT(IN) ::NTIME
  LOGICAL, INTENT(IN) ::D3,DCYL
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: N,NLX,NLY,NLZ
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: M,MLX,MLY,MLZ
  REAL, INTENT(IN),DIMENSION(NGI)  ::WEIGHT
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC)::PNDGLN
  ! coriolos force
  INTEGER, INTENT(IN)  :: GEOBAL
  REAL, INTENT(IN)  :: SCFACTH0
  INTEGER NOCVA
  REAL G(NOCVA)
  ! Wind stress
  !------------
  LOGICAL IFWIND

  ! observational
  integer istate
  integer SnapNDT_obs
  integer,INTENT(IN),DIMENSION(0:ntime,mxnods) ::iflagobs
  REAL,INTENT(IN),DIMENSION(0:ntime,istate):: obsv
  REAL coef_multivar_obs(0:ntime,NSVD_TOTAL)
  REAL leftsvd_obs(ISTATE,nsvd)
  REAl soux_mean,souy_mean,souz_mean,soux1,souy1,souz1,soux2,souy2,souz2
  REAL smean_obs(istate)
  REAL,DIMENSION(:),ALLOCATABLE   :: Bmat_obs
  !       Local variables
  REAL    ::THETA,ONEMINTH
  INTEGER ::ISVD,JSVD,II,JJ,KK,NSTEP,ELE,GI,ILOC,JLOC,I,J,K
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
            U_xx_mean,U_yy_mean,U_zz_mean,V_xx_mean,V_yy_mean,V_zz_mean,W_xx_mean,W_yy_mean,W_zz_mean,  &
            UDmean_UADJx,VDmean_UADJy,WDmean_UADJz,UDmean_VADJx,VDmean_VADJy,WDmean_VADJz, &
            UDmean_WADJx,VDmean_WADJy,WDmean_WADJz,                                        &
            UD_UADJx,VD_UADJy,WD_UADJz,UD_VADJx,VD_VADJy,WD_VADJz,UD_WADJx,VD_WADJy,WD_WADJz, &
            Umeanx_UADJ,Umeany_VADJ,Umeanz_WADJ,Vmeanx_UADJ,Vmeany_VADJ,                   &
            Vmeanz_WADJ,Wmeanx_UADJ,Wmeany_VADJ,Wmeanz_WADJ,                               &         
            Ux_UADJ,Uy_VADJ,Uz_WADJ,Vx_UADJ,Vy_VADJ,Vz_WADJ,Wx_UADJ,Wy_VADJ,Wz_WADJ

  INTEGER ::NCOLM,NCOLM1
  INTEGER ::IBL11,IBL12,IBL13,IBL14
  INTEGER ::IBL21,IBL22,IBL23,IBL24
  INTEGER ::IBL31,IBL32,IBL33,IBL34
  INTEGER ::IBL41,IBL42,IBL43,IBL44
  real    ::VOL,ftest,VOLUME


  !       Working arrays/local variables            
  !       ********* DEFINE OPTIONS FOR GMRES **************
  REAL     ::ERROR1,ERROR2,RELAX
  INTEGER  ::NOITS1,NOITS2,LEFTP1,LEFTP2,IMM
  INTEGER  ::MINTE1,MINTE2,NR2,IP,MISOUT
  !       ** options for MULMAT 
  LOGICAL  ::SYM,IGUESS
  INTEGER  ::MULPA,UZAWA,CGSOLQ,GMRESQ,CGSQ
  REAL     ::DM1(NSVD),DM2(NSVD),R(7*NSVD)
  real     ::rub(2)
  integer  ::iidum(2),kits
  INTEGER  ::MATSTR,NR,NNV
  PARAMETER(MATSTR=0)
  REAL     ::RRR
  INTEGER  ::PARA,halo_tag, halo_tag_p
  PARAMETER(PARA=0,halo_tag=1, halo_tag_p=2)

! all matraces for the calculations:
  INTEGER  FINDRM(NSVD+1),COLM(NSVD*NSVD), CENTRM(NSVD)

  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM,VEC,VEC1,VEC_extr
  REAL,DIMENSION(:,:),ALLOCATABLE   :: ReMat,ReMat1
  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM0,BMat,BMat1,BMat0,BMat01,VEC0,VEC0_source,BMat0_source,BMat_extr
  REAL,DIMENSION(:),ALLOCATABLE     :: DETWEI,XD,YD,ZD
  REAL,DIMENSION(:),ALLOCATABLE     :: UD_mean,VD_mean,WD_mean
  REAL,DIMENSION(:,:),ALLOCATABLE   :: UD,VD,WD
  REAL,DIMENSION(:),ALLOCATABLE     :: UDx_mean,UDy_mean,UDz_mean,VDx_mean,VDy_mean,VDz_mean,WDx_mean,WDy_mean,WDz_mean
  REAL,DIMENSION(:,:),ALLOCATABLE   :: UDx,UDy,UDz,VDx,VDy,VDz,WDx,WDy,WDz
  REAL,DIMENSION(:),ALLOCATABLE     :: Umeanx,Umeany,Umeanz,Vmeanx,Vmeany,Vmeanz,Wmeanx,Wmeany,Wmeanz
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_u,WEI_v,WEI_w,WEI_p
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_ux,WEI_uy,WEI_uz,WEI_vx,WEI_vy,WEI_vz,WEI_wx,WEI_wy,WEI_wz
  REAL,DIMENSION(:,:),ALLOCATABLE   :: NX,NY,NZ,MX,MY,MZ
  REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF,NEWCOEF1
  REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF_ADJ,NEWCOEF1_ADJ
  REAL,DIMENSION(:),ALLOCATABLE     :: CORIOLIS

  REAL,DIMENSION(:,:),ALLOCATABLE   :: BIGM1,BIGM1consi
  REAL,DIMENSION(:),ALLOCATABLE     :: UDMat
!  REAL,DIMENSION(:),ALLOCATABLE     :: FINDRM,COLM,CENTRM
  REAL,DIMENSION(:),ALLOCATABLE     :: ML
  REAL,DIMENSION(:),ALLOCATABLE     :: SOURCX,SOURCY,SOURCZ


  ALLOCATE(BIGM(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(VEC(NSVD_TOTAL))
  ALLOCATE(ReMat( NSVD_TOTAL-NSVD,NSVD_TOTAL-NSVD))
  ALLOCATE(BIGM0(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BIGM1(NSVD_TOTAL*NSVD_TOTAL,NSVD_u))
  ALLOCATE(BIGM1consi(NSVD_TOTAL*NSVD_TOTAL,NSVD_u))
  ALLOCATE(VEC0(NSVD_TOTAL))
  ALLOCATE(VEC1(NSVD_TOTAL-NSVD))
  ALLOCATE(VEC_extr(NSVD_TOTAL))
  ALLOCATE(BMat_extr(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat0(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(BMat01(NSVD_TOTAL*NSVD_TOTAL))
!  ALLOCATE(BMat1(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(DETWEI(NGI))
  ALLOCATE(UD(NGI,NSVD))
  ALLOCATE(VD(NGI,NSVD))
  ALLOCATE(WD(NGI,NSVD))
  ALLOCATE(UD_mean(NGI))
  ALLOCATE(VD_mean(NGI))
  ALLOCATE(WD_mean(NGI))
  ALLOCATE(UDx(NGI,NSVD))
  ALLOCATE(UDy(NGI,NSVD))
  ALLOCATE(UDz(NGI,NSVD))
  ALLOCATE(VDx(NGI,NSVD))
  ALLOCATE(VDy(NGI,NSVD))
  ALLOCATE(VDz(NGI,NSVD))
  ALLOCATE(WDx(NGI,NSVD))
  ALLOCATE(WDy(NGI,NSVD))
  ALLOCATE(WDz(NGI,NSVD))
  ALLOCATE(UDx_mean(NGI))
  ALLOCATE(UDy_mean(NGI))
  ALLOCATE(UDz_mean(NGI))
  ALLOCATE(VDx_mean(NGI))
  ALLOCATE(VDy_mean(NGI))
  ALLOCATE(VDz_mean(NGI))
  ALLOCATE(WDx_mean(NGI))
  ALLOCATE(WDy_mean(NGI))
  ALLOCATE(WDz_mean(NGI))

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
  ALLOCATE(NEWCOEF_ADJ(NSVD_TOTAL))
  ALLOCATE(NEWCOEF1_ADJ(NSVD_TOTAL-NSVD))
  ALLOCATE(CORIOLIS(NGI))

!  ALLOCATE(SOURCX(NSVD))
!  ALLOCATE(SOURCY(NSVD))
!  ALLOCATE(SOURCZ(NSVD))
  ALLOCATE(BMat_OBS(NSVD_TOTAL*NSVD_TOTAL))
  ALLOCATE(VEC0_source(NSVD_TOTAL))
  ALLOCATE(BMat0_source(NSVD_TOTAL*NSVD_TOTAL))

!**************************************************!
!       Initialise all matrices                    !
!**************************************************!
  BIGM =0.0
  VEC =0.0
  ReMat=0.0
  BIGM0=0.0
  BIGM1=0.0
  BIGM1consi=0.0
  VEC0=0.0
  VEC1=0.0
  VEC_extr=0.0
  BMat_extr=0.0
  BMat=0.0
  BMat0=0.0
  BMat01=0.0
!  BMat1=0.0
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
  UDx=0.0
  UDy=0.0
  UDz=0.0              
  VDx=0.0
  VDy=0.0
  VDz=0.0              
  WDx=0.0
  WDy=0.0
  WDz=0.0
  UDx_mean=0.0
  UDy_mean=0.0
  UDz_mean=0.0              
  VDx_mean=0.0
  VDy_mean=0.0
  VDz_mean=0.0              
  WDx_mean=0.0
  WDy_mean=0.0
  WDz_mean=0.0

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
  NEWCOEF_ADJ=0.0
  NEWCOEF1_ADJ=0.0
  CORIOLIS=0.0


  NX=0.0
  NY=0.0
  NZ=0.0
  MX=0.0
  MY=0.0
  MZ=0.0

!  SOURCX=0.0
!  SOURCY=0.0
!  SOURCZ=0.0

   VEC0_source=0.0
   BMat0_source=0.0
   BMat_OBS =0.0
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

! only for 2D!!!!
  DO II=2*NSVD_U+1,3*NSVD_U
     COEF_ADJ(0,II)=0.0
  ENDDO

!     LEFTSVD(2*nonods+1:3*nonods,:)=0.0
     SMEAN(2*mxnods+1:3*mxnods)=0.0

! for adjoint model
     NEWCOEF_ADJ(:)=0.0

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


!  DO 100 NSTEP=1,1
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
              UD_mean(GI)=0.0
              VD_mean(GI)=0.0
              WD_mean(GI)=0.0
              UD(GI,:)=0.0
              VD(GI,:)=0.0
              WD(GI,:)=0.0
              XD(GI)=0.0
              YD(GI)=0.0
              ZD(GI)=0.0
              UDx(GI,:)=0.0
              UDy(GI,:)=0.0
              UDz(GI,:)=0.0              
              VDx(GI,:)=0.0
              VDy(GI,:)=0.0
              VDz(GI,:)=0.0              
              WDx(GI,:)=0.0
              WDy(GI,:)=0.0
              WDz(GI,:)=0.0
              UDx_mean(GI)=0.0
              UDy_mean(GI)=0.0
              UDz_mean(GI)=0.0              
              VDx_mean(GI)=0.0
              VDy_mean(GI)=0.0
              VDz_mean(GI)=0.0              
              WDx_mean(GI)=0.0
              WDy_mean(GI)=0.0
              WDz_mean(GI)=0.0

              CORIOLIS(GI)=0.0

              DO JLOC=1,NLOC
                 JGL=NDGLNO((ELE-1)*NLOC+JLOC)

                 UD_mean(GI)=UD_mean(GI)+SMEAN(JGL)*N(JLOC,GI)
                 VD_mean(GI)=VD_mean(GI)+SMEAN(1*MXNODS+JGL)*N(JLOC,GI)
                 WD_mean(GI)=WD_mean(GI)+SMEAN(2*MXNODS+JGL)*N(JLOC,GI)

                 UDx_mean(GI)=UDx_mean(GI)+SMEAN(JGL)*NX(JLOC,GI)
                 UDy_mean(GI)=UDy_mean(GI)+SMEAN(JGL)*NY(JLOC,GI)
                 UDz_mean(GI)=UDz_mean(GI)+SMEAN(JGL)*NZ(JLOC,GI)

                 VDx_mean(GI)=VDx_mean(GI)+SMEAN(1*MXNODS+JGL)*NX(JLOC,GI)
                 VDy_mean(GI)=VDy_mean(GI)+SMEAN(1*MXNODS+JGL)*NY(JLOC,GI)
                 VDz_mean(GI)=VDz_mean(GI)+SMEAN(1*MXNODS+JGL)*NZ(JLOC,GI)

                 WDx_mean(GI)=WDx_mean(GI)+SMEAN(2*MXNODS+JGL)*NX(JLOC,GI)
                 WDy_mean(GI)=WDy_mean(GI)+SMEAN(2*MXNODS+JGL)*NY(JLOC,GI)
                 WDz_mean(GI)=WDz_mean(GI)+SMEAN(2*MXNODS+JGL)*NZ(JLOC,GI)

                 XD(GI)=XD(GI) + N(JLOC,GI)*X(JGL)
                 YD(GI)=YD(GI) + N(JLOC,GI)*Y(JGL)
                 ZD(GI)=ZD(GI) + N(JLOC,GI)*Z(JGL)

                 DO JSVD = 1,NSVD
                    UD(GI,JSVD)=UD(GI,JSVD)+LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    VD(GI,JSVD)=VD(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    WD(GI,JSVD)=WD(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)

                    UDx(GI,JSVD)=UDx(GI,JSVD)+LEFTSVD(JGL,JSVD)*NX(JLOC,GI)
                    UDy(GI,JSVD)=UDy(GI,JSVD)+LEFTSVD(JGL,JSVD)*NY(JLOC,GI)
                    UDz(GI,JSVD)=UDz(GI,JSVD)+LEFTSVD(JGL,JSVD)*NZ(JLOC,GI)

                    VDx(GI,JSVD)=VDx(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    VDy(GI,JSVD)=VDy(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    VDz(GI,JSVD)=VDz(GI,JSVD)+LEFTSVD(1*MXNODS+JGL,JSVD)*NZ(JLOC,GI)

                    WDx(GI,JSVD)=WDx(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*NX(JLOC,GI)
                    WDy(GI,JSVD)=WDy(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*NY(JLOC,GI)
                    WDz(GI,JSVD)=WDz(GI,JSVD)+LEFTSVD(2*MXNODS+JGL,JSVD)*NZ(JLOC,GI)
                 ENDDO
              ENDDO
                 CORIOLIS(GI) = 2.0*FUNOME(XD(GI),YD(GI),ZD(GI))
           ENDDO

! Source terms which represet the misfit between the numerical results and observations

           soux_mean=0.0
           souy_mean=0.0
           souz_mean=0.0

           DO JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
              DO GI=1,NGI
                 soux_mean=soux_mean + WEI_u(GI)*DETWEI(GI)*(SMEAN(JGL)-SMEAN_OBS(JGL))*N(JLOC,GI)
                 souy_mean=souy_mean + WEI_v(GI)*DETWEI(GI)*(SMEAN(1*MXNODS+JGL)-SMEAN_OBS(1*MXNODS+JGL))*N(JLOC,GI)
                 souz_mean=souz_mean + WEI_w(GI)*DETWEI(GI)*(SMEAN(2*MXNODS+JGL)-SMEAN_OBS(2*MXNODS+JGL))*N(JLOC,GI)
              ENDDO
           ENDDO

           VEC0_source(ISVD)= VEC0_source(ISVD)+DT*soux_mean
           VEC0_source(1*NSVD+ISVD)= VEC0_source(1*NSVD+ISVD)+DT*souy_mean
           VEC0_source(2*NSVD+ISVD)= VEC0_source(2*NSVD+ISVD)+DT*souz_mean
! end the source terms


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

                 ! source terms for misfit
                 SOUX1 =0.0
                 SOUY1 =0.0
                 SOUZ1 =0.0
                 SOUX2 =0.0
                 SOUY2 =0.0
                 SOUZ2 =0.0

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


                    ! source terms for misfit
                    SOUX1 = SOUX1 + WEI_u(GI)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)
                    SOUY1 = SOUY1 + WEI_v(GI)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
!                    SOUZ1 = SOUZ1 + WEI_w(GI)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
                    
                    SOUX2 = SOUX2 + WEI_u(GI)*DETWEI(GI)*LEFTSVD_OBS(JGL,JSVD)*N(JLOC,GI)
                    SOUY2 = SOUY2 + WEI_v(GI)*DETWEI(GI)*LEFTSVD_OBS(1*MXNODS+JGL,JSVD)*N(JLOC,GI)
!                    SOUZ2 = SOUZ2 + WEI_w(GI)*DETWEI(GI)*LEFTSVD_OBS(2*MXNODS+JGL,JSVD)*N(JLOC,GI)

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       !BIGM1consi(IBL11+COUNT,JJ)=BIGM1consi(IBL11+COUNT,JJ)+Ux_UADJ
                       !--------------------------------------------------------------------------------
                       BIGM1consi(IBL11+COUNT,JJ)=BIGM1consi(IBL11+COUNT,JJ)+   &
                                WEI_u(GI)*UDx(GI,JJ)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)                          

                       !BIGM1consi(IBL12+COUNT,JJ)=BIGM1consi(IBL12+COUNT,JJ)+Uy_VADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL12+COUNT,JJ)=BIGM1consi(IBL12+COUNT,JJ)+   &
                          WEI_u(GI)*UDy(GI,JJ)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL13+COUNT,JJ)=BIGM1consi(IBL13+COUNT,JJ)+Uz_WADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL13+COUNT,JJ)=BIGM1consi(IBL13+COUNT,JJ)+   &
                          WEI_u(GI)*UDz(GI,JJ)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
!                      
!
                       !BIGM1consi(IBL21+COUNT,JJ)=BIGM1consi(IBL21+COUNT,JJ)+Vx_UADJ
                       !----------------------------------------------------
                       BIGM1consi(IBL21+COUNT,JJ)=BIGM1consi(IBL21+COUNT,JJ)+   &
                          WEI_v(GI)*VDx(GI,JJ)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL22+COUNT,JJ)=BIGM1consi(IBL22+COUNT,JJ)+Vy_VADJ)
                       !--------------------------------------------------------------------------------
                       BIGM1consi(IBL22+COUNT,JJ)=BIGM1consi(IBL22+COUNT,JJ)+   & 
                          WEI_v(GI)*VDy(GI,JJ)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL23+COUNT,JJ)=BIGM1consi(IBL23+COUNT,JJ)+Vz_WADJ
                       !----------------------------------------------------
                       BIGM1consi(IBL23+COUNT,JJ)=BIGM1consi(IBL23+COUNT,JJ)+  &
                          WEI_v(GI)*VDz(GI,JJ)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)
!
!
                       !BIGM1consi(IBL31+COUNT,JJ)=BIGM1consi(IBL31+COUNT,JJ)+Wx_UADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL31+COUNT,JJ)=BIGM1consi(IBL31+COUNT,JJ)+  &
                          WEI_w(GI)*WDx(GI,JJ)*DETWEI(GI)*LEFTSVD(JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL32+COUNT,JJ)=BIGM1consi(IBL32+COUNT,JJ)+Wy_VADJ
                       !---------------------------------------------------
                       BIGM1consi(IBL32+COUNT,JJ)=BIGM1consi(IBL32+COUNT,JJ)+  &
                          WEI_w(GI)*WDz(GI,JJ)*DETWEI(GI)*LEFTSVD(1*MXNODS+JGL,JSVD)*N(JLOC,GI)

                       !BIGM1consi(IBL33+COUNT,JJ)=BIGM1consi(IBL33+COUNT,JJ)+Wz_WADJ
                       !--------------------------------------------------------------------------------
                       BIGM1consi(IBL33+COUNT,JJ)=BIGM1consi(IBL33+COUNT,JJ)+  &
                          WEI_w(GI)*WDz(GI,JJ)*DETWEI(GI)*LEFTSVD(2*MXNODS+JGL,JSVD)*N(JLOC,GI)

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
 
              BIGM0(IBL11+COUNT)=BIGM0(IBL11+COUNT)+RN1+DT*THETA*(UDmean_Ux+VDmean_Uy+WDmean_Uz+U_xx+U_yy+U_zz)
              BIGM0(IBL12+COUNT)=BIGM0(IBL12+COUNT)+DT*THETA*FV
              BIGM0(IBL13+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                 BIGM0(IBL14+COUNT1)=BIGM0(IBL14+COUNT1)+DT*P_x
              ENDIF

              BIGM0(IBL21+COUNT)=BIGM0(IBL21+COUNT)+DT*THETA*FU
              BIGM0(IBL22+COUNT)=BIGM0(IBL22+COUNT)+RN2+DT*THETA*(UDmean_Vx+VDmean_Vy+WDmean_Vz+V_xx+V_yy+V_zz)
              BIGM0(IBL23+COUNT)=0.0
              IF(JSVD.LE.NSVD_PHI) THEN
                  BIGM0(IBL24+COUNT1)=BIGM0(IBL24+COUNT1)+DT*P_y
              ENDIF
              BIGM0(IBL31+COUNT)=0.0
              BIGM0(IBL32+COUNT)=0.0
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
              BMat0(IBL12+COUNT)=BMat0(IBL12+COUNT)-DT*ONEMINTH*FV
              BMat0(IBL13+COUNT)=0.0

 
              BMat0(IBL21+COUNT)=BMat0(IBL21+COUNT)-DT*ONEMINTH*FU
              BMat0(IBL22+COUNT)=BMat0(IBL22+COUNT)+RN2-DT*ONEMINTH*(UDmean_Vx+VDmean_Vy+WDmean_Vz+V_xx+V_yy+V_zz)
              BMat0(IBL23+COUNT)=0.0

              BMat0(IBL31+COUNT)=0.0
              BMat0(IBL32+COUNT)=0.0
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

              ! source terms for misfit
              BMat0_source(IBL11+COUNT)=BMat0_source(IBL11+COUNT)+DT*SOUX1
              BMat0_source(IBL22+COUNT)=BMat0_source(IBL22+COUNT)+DT*SOUY1
              BMat0_source(IBL33+COUNT)=BMat0_source(IBL33+COUNT)+DT*SOUZ1
!
              BMat_obs(IBL11+COUNT)=BMat_obs(IBL11+COUNT)+DT*SOUX2
              BMat_obs(IBL22+COUNT)=BMat_obs(IBL22+COUNT)+DT*SOUY2
              BMat_obs(IBL33+COUNT)=BMat_obs(IBL33+COUNT)+DT*SOUZ2
40            CONTINUE
!              ewrite(3,*) 'vec2',VEC(ISVD)
!              ewrite(3,*) 'BIGM(COUNT)',BIGM(COUNT)
!              if(ISVD.EQ.2) stop

!60            CONTINUE
20            CONTINUE

30            CONTINUE




   DO 100 NSTEP=1,ntime

!      CALL Source_POD(NONODS,NTIME*DT,NSVD,DT,NTIME+1,NSTEP+1,     &
!                    SOURCX,SOURCY,SOURCZ,                          &
!                    SMEAN,LEFTSVD,COEF(:,NSTEP),                   &
!                    TOTELE,NLOC,NGI,                               &
!                    NDGLNO,XNONOD,xondgl,X,Y,Z,                    &
!                    N,NLX,NLY,NLZ,WEIGHT,D3,DCYL)

      BIGM(:)=BIGM0(:)
      BMat(:)=BMat0(:)
      VEC(:)=VEC0(:)
!      if (SnapNDT_obs*int(nstep/SnapNDT_obs) .eq. nstep) then  
         VEC(:)=VEC(:)+VEC0_source(:)
!         BMat(:)=BMat(:)+BMat0_source(:)
!      endif

      DO ISVD=1,NSVD
         DO JSVD=1,NSVD
            COUNT=(ISVD-1)*NSVD+JSVD
            COUNT1=(ISVD-1)*NSVD_PHI+JSVD
            DO JJ=1,NSVD
               BIGM(IBL11+COUNT)=BIGM(IBL11+COUNT)+DT*THETA*( BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ))

               BIGM(IBL22+COUNT)=BIGM(IBL22+COUNT)+DT*THETA*( BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ))

               BIGM(IBL33+COUNT)=BIGM(IBL33+COUNT)+DT*THETA*( BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ))


               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)-DT*ONEMINTH*( BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ) )


               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)-DT*ONEMINTH*( BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ) )

               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)-DT*ONEMINTH*( BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)+       &
                    BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)+BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ) )


                    ! extra term for nonlinear
              if(.true.) then
               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)+  &
                    DT*THETA*BIGM1consi(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+          &
                    DT*ONEMINTH*BIGM1consi(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL12+COUNT)=BMat(IBL12+COUNT)+DT*THETA*BIGM1consi(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1consi(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL13+COUNT)=BMat(IBL13+COUNT)+DT*THETA*BIGM1consi(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1consi(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)


               BMat(IBL21+COUNT)=BMat(IBL21+COUNT)+DT*THETA*BIGM1consi(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1consi(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)+  &
                    DT*THETA*BIGM1consi(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+  &
                    DT*ONEMINTH*BIGM1consi(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL23+COUNT)=BMat(IBL23+COUNT)+DT*THETA*BIGM1consi(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1consi(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)


               BMat(IBL31+COUNT)=BMat(IBL31+COUNT)+DT*THETA*BIGM1consi(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1consi(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL32+COUNT)=BMat(IBL32+COUNT)+DT*THETA*BIGM1consi(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+               & 
                    DT*ONEMINTH*BIGM1consi(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)+    &
                    DT*THETA*BIGM1consi(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+   &
                    DT*ONEMINTH*BIGM1consi(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
              else
               BMat(IBL11+COUNT)=BMat(IBL11+COUNT)+  &
                    DT*THETA*BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+          &
                    DT*ONEMINTH*BIGM1(IBL11+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL12+COUNT)=BMat(IBL12+COUNT)+DT*THETA*BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1(IBL12+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)
               BMat(IBL13+COUNT)=BMat(IBL13+COUNT)+DT*THETA*BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP+1,JJ)+        &
                    DT*ONEMINTH*BIGM1(IBL13+COUNT,JJ)*COEF(NTIME-NSTEP,JJ)


               BMat(IBL21+COUNT)=BMat(IBL21+COUNT)+DT*THETA*BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1(IBL21+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL22+COUNT)=BMat(IBL22+COUNT)+  &
                    DT*THETA*BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+  &
                    DT*ONEMINTH*BIGM1(IBL22+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)
               BMat(IBL23+COUNT)=BMat(IBL23+COUNT)+DT*THETA*BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP+1,1*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1(IBL23+COUNT,JJ)*COEF(NTIME-NSTEP,1*NSVD+JJ)


               BMat(IBL31+COUNT)=BMat(IBL31+COUNT)+DT*THETA*BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+                & 
                    DT*ONEMINTH*BIGM1(IBL31+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL32+COUNT)=BMat(IBL32+COUNT)+DT*THETA*BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+               & 
                    DT*ONEMINTH*BIGM1(IBL32+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
               BMat(IBL33+COUNT)=BMat(IBL33+COUNT)+    &
                    DT*THETA*BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP+1,2*NSVD+JJ)+   &
                    DT*ONEMINTH*BIGM1(IBL33+COUNT,JJ)*COEF(NTIME-NSTEP,2*NSVD+JJ)
              endif

            ENDDO


               VEC(ISVD)=VEC(ISVD)+BMat(IBL11+COUNT)*COEF_ADJ(NSTEP-1,JSVD)+                       &
                    BMat(IBL12+COUNT)*COEF_ADJ(NSTEP-1,1*NSVD+JSVD)+BMat(IBL13+COUNT)*COEF_ADJ(NSTEP-1,2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(ISVD)=VEC(ISVD)+BMat(IBL14+COUNT1)*COEF_ADJ(NSTEP-1,3*NSVD+JSVD)
               ENDIF

               VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat(IBL21+COUNT)*COEF_ADJ(NSTEP-1,JSVD)+         &
                    BMat(IBL22+COUNT)*COEF_ADJ(NSTEP-1,1*NSVD+JSVD)+BMat(IBL23+COUNT)*COEF_ADJ(NSTEP-1,2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat(IBL24+COUNT1)*COEF_ADJ(NSTEP-1,3*NSVD+JSVD)
               ENDIF

               VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat(IBL31+COUNT)*COEF_ADJ(NSTEP-1,JSVD)+         &
                    BMat(IBL32+COUNT)*COEF_ADJ(NSTEP-1,1*NSVD+JSVD)+BMat(IBL33+COUNT1)*COEF_ADJ(NSTEP-1,2*NSVD+JSVD)
               IF(JSVD.LE.NSVD_PHI) THEN
                 VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat(IBL34+COUNT1)*COEF_ADJ(NSTEP-1,3*NSVD+JSVD)
               ENDIF

            ! for source terms
            !       if (SnapNDT_obs*int(nstep/SnapNDT_obs) .eq. nstep) then  
            VEC(ISVD)=VEC(ISVD)+BMat0_source(IBL11+COUNT)*coef(NTIME-NSTEP,JSVD) &
                   -BMat_OBS(IBL11+COUNT)*coef_multivar_obs(NTIME-NSTEP,JSVD)
            VEC(1*NSVD+ISVD)=VEC(1*NSVD+ISVD)+BMat0_source(IBL22+COUNT)*coef(NTIME-NSTEP,1*NSVD+JSVD)  &
                   -BMat_OBS(IBL22+COUNT)*coef_multivar_obs(NTIME-NSTEP,1*NSVD+JSVD)
            VEC(2*NSVD+ISVD)=VEC(2*NSVD+ISVD)+BMat0_source(IBL33+COUNT)*coef(NTIME-NSTEP,2*NSVD+JSVD)  &
                   -BMat_OBS(IBL33+COUNT)*coef_multivar_obs(NTIME-NSTEP,2*NSVD+JSVD)
            !       endif
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

              Remat = 0.0
              VEC1 = 0.0

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
                             ReMat( (kk-2)*NSVD+ii,2*NSVD+jj)= BIGM( IBL44+(ii-1)*NSVD_PHI+jj)
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

              open(1,file='Remat1.dat')
              open(2,file='vec1.dat')
              write(1,*) Remat
              write(2,*) vec
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
           write(10, *) (NEWCOEF1_ADJ(ii), ii=1,(nvar-1)*nsvd)
        close(10) 


              CALL GMRES(NEWCOEF_ADJ,VEC,NSVD,NSVD,NSVD,.FALSE.,                        & 
                   BIGM,BIGM,FINDRM,COLM,NSVD*NSVD,NSVD*NSVD,&
                   halo_tag,KITS)



              DO II=1,NSVD_TOTAL
                 IF(II.LE.2*NSVD) THEN
                    COEF_ADJ(NSTEP,II)=NEWCOEF1_ADJ(II)
                 ELSEIF(II.GT.3*NSVD) THEN
                    COEF_ADJ(NSTEP,II)=NEWCOEF1_ADJ(II-NSVD)
                 ELSE
                    COEF_ADJ(NSTEP,II)=0.0
                 ENDIF
                 ewrite(3,*) 'COEF_ADJ(NSTEP,II)',COEF_ADJ(NSTEP,II)
              ENDDO

              ewrite(3,*) 'after NSTEP',NSTEP
100           CONTINUE


              open(10, file ='COEF_ADJ.dat')
              do ii=0,ntime
                 write(10, *) (COEF_ADJ(ii,jj), jj=1,nsvd)
              enddo
              close(10)
              ewrite(3,*) 'exit REDUCED_MODEL_DIFFCOEF'

! Calculate the gradient for the intial inversion
!-------------------------------------------------
      DO ISVD=1,NSVD
         DO JSVD=1,NSVD
           IF(ISVD.EQ.JSVD) THEN
             COUNT=(ISVD-1)*NSVD+JSVD
             COUNT1=(ISVD-1)*NSVD_PHI+JSVD
            DO JJ=1,NSVD
              
              if(.true.) then
               BMat_extr(IBL11+COUNT)=BMat_extr(IBL11+COUNT)+ &
                    DT*ONEMINTH*BIGM1consi(IBL11+COUNT,JJ)*COEF(0,JJ)+BMat(IBL11+COUNT)

               BMat_extr(IBL22+COUNT)=BMat_extr(IBL22+COUNT)+  &
                    DT*ONEMINTH*BIGM1consi(IBL22+COUNT,JJ)*COEF(0,1*NSVD+JJ)+BMat(IBL22+COUNT)

               BMat_extr(IBL33+COUNT)=BMat_extr(IBL33+COUNT)+    &
                    DT*ONEMINTH*BIGM1consi(IBL33+COUNT,JJ)*COEF(0,2*NSVD+JJ)+BMat(IBL33+COUNT)
              else
               BMat_extr(IBL11+COUNT)=BMat_extr(IBL11+COUNT)+  &
                    DT*ONEMINTH*BIGM1(IBL11+COUNT,JJ)*COEF(0,JJ)+BMat(IBL11+COUNT)

               BMat_extr(IBL22+COUNT)=BMat_extr(IBL22+COUNT)+  &
                    DT*ONEMINTH*BIGM1(IBL22+COUNT,JJ)*COEF(0,1*NSVD+JJ)+BMat(IBL22+COUNT)


               BMat_extr(IBL33+COUNT)=BMat_extr(IBL33+COUNT)+    &
                    DT*ONEMINTH*BIGM1(IBL33+COUNT,JJ)*COEF(0,2*NSVD+JJ)+BMat(IBL33+COUNT)
              endif
            ENDDO
            VEC_extr(ISVD)=BMat_extr(IBL11+COUNT)
            VEC_extr(1*NSVD+ISVD)=BMat_extr(IBL22+COUNT)
            VEC_extr(2*NSVD+ISVD)=BMat_extr(IBL33+COUNT)
           ENDIF
         ENDDO
      ENDDO

      DO II=1,3*nsvd
        G(ii) = coef_ADJ(NTIME-1,ii)*VEC_extr(ii)
      ENDDO


!-------------------------------!
! DEALLOCATE THE LOCAL MATRICES !
!-------------------------------!

  DEALLOCATE(BIGM)
  DEALLOCATE(VEC)
  DEALLOCATE(ReMat)
  DEALLOCATE(BIGM0)
  DEALLOCATE(BIGM1)
  DEALLOCATE(BIGM1consi)
  DEALLOCATE(VEC0)
  DEALLOCATE(VEC1)
  DEALLOCATE(BMat)
  DEALLOCATE(BMat0)
  DEALLOCATE(BMat01)
  DEALLOCATE(DETWEI)
  DEALLOCATE(UD)
  DEALLOCATE(VD)
  DEALLOCATE(WD)
  DEALLOCATE(UD_mean)
  DEALLOCATE(VD_mean)
  DEALLOCATE(WD_mean)
  DEALLOCATE(UDx)
  DEALLOCATE(UDy)
  DEALLOCATE(UDz)
  DEALLOCATE(VDx)
  DEALLOCATE(VDy)
  DEALLOCATE(VDz)
  DEALLOCATE(WDx)
  DEALLOCATE(WDy)
  DEALLOCATE(WDz)

  DEALLOCATE(UDx_mean)
  DEALLOCATE(UDy_mean)
  DEALLOCATE(UDz_mean)
  DEALLOCATE(VDx_mean)
  DEALLOCATE(VDy_mean)
  DEALLOCATE(VDz_mean)
  DEALLOCATE(WDx_mean)
  DEALLOCATE(WDy_mean)
  DEALLOCATE(WDz_mean)

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
  DEALLOCATE(NEWCOEF_ADJ)
  DEALLOCATE(NEWCOEF1_ADJ)
  DEALLOCATE(CORIOLIS)

!DEALLOCATE(SOURCX)
!DEALLOCATE(SOURCY)
!DEALLOCATE(SOURCZ)

  DEALLOCATE(VEC0_source)
  DEALLOCATE(BMat_OBS)
  DEALLOCATE(BMat0_source)
  DEALLOCATE(VEC_extr)
  DEALLOCATE(BMat_extr)

  RETURN
END SUBROUTINE REDUCED_MODEL_DIFFCOEF1consistent_ADJ

Subroutine Vorticity2D(X,Y,Z, &
     U, V, &
     ndglno, totele, nonods, nloc, ngi, &
     D3,DCYL, & 
     Fu, Fv, &
     N,NLX,NLY,NLZ, &
     WEIGHT)
  use detnlxr_module
Implicit None

Real, Intent(In)::X(Nonods),Y(Nonods),Z(Nonods)
Real, Intent(In)::U(Nonods),V(Nonods)
Real, Intent(Out)::FU(Nonods),FV(Nonods)
Integer, intent(in)::totele, nloc, ndglno(nloc*totele), nonods, ngi
Logical, Intent(in)::D3,DCYL

REAL, Intent(In)::N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
REAL, Intent(in):: WEIGHT(NGI) 

Real, Allocatable::DetWei(:), RA(:)
REAL, Allocatable::NX(:,:),NY(:,:),NZ(:,:), newmass(:)

Integer ele, gi, inod, iloc, jNod, jloc
REAL  nXx, nyx, nxy, nyY, VOLUME, nn

Allocate(DetWei(NGI), RA(NGI), newmass(Nonods))
Allocate(NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI))


Fu =0.0
Fv=0.0
newmass=0.0

DO ELE=1,Totele

! Calculate DETWEI,RA,NX,NY,NZ for element ELE
   CALL DETNLXR(ELE, X,Y,Z, NDGLNO, TOTELE,NONODS,NLOC,NGI, &
        N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, D3,DCYL, & 
        NX,NY,NZ) 

! Vorticity

   DO iLoc=1,NLOC
      iNod = NDGLNO((ELE-1)*NLOC+iLoc)
      
      DO jLoc=1,NLOC
         jNod = NDGLNO((ELE-1)*NLOC+jLoc)

         nXx=0.
         nyx=0.
         nxy=0.
         nyy=0.
         nn=0.0

         DO GI=1,NGI
         
            NXX=NXX+NX(iLoc,GI)*NX(JLoc,GI)*DETWEI(GI)
            NXY=NXY+NX(iLoc,GI)*NY(JLoc,GI)*DETWEI(GI)
            NYX=NYX+NY(iLoc,GI)*NX(JLoc,GI)*DETWEI(GI)
            NYY=NYY+NY(iLoc,GI)*NY(JLoc,GI)*DETWEI(GI)
            
            nn = nn + N(iLoc,GI)*N(JLoc,GI)*DETWEI(GI)

      !      FU(iNod) = FU(iNod) +(NX(jloc,gi)*V(JNod) - NY(jloc,gi)*U(JNod))*DetWei(gi)
       !     FV(iNod) = FV(iNod) +(NX(jloc,gi)*V(JNod) - NY(jloc,gi)*U(JNod))*DetWei(gi)

         END DO

         newmass(iNod) = newmass(iNod) + nn

         FU(iNod) = FU(iNod) -NYX*V(JNod) + NYY*U(JNod)
         FV(iNod) = FV(iNod) +NXX*V(JNod) - NXY*U(JNod) 
                 
      End Do
   End do
   
End Do

DeAllocate(DetWei,RA)
DeAllocate(NX,NY,NZ)

do inod=1,nonods
   FU(iNod) = FU(iNod)/newmass(iNod)
   FV(iNod) = FV(iNod)/newmass(iNod)
end do

Deallocate(newmass)

End Subroutine Vorticity2D

     SUBROUTINE REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ_errormeasure(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,fredop,&
          X,Y,Z,            &
          FERROR,COEF_ADJ,U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ,  &
          U0,V0,W0,PHI0,U0_OBS,V0_OBS,W0_OBS,NCOLM,COLM,FINDRM,NCT,COLCT,FINDCT,NCT_PG,COLCT_PG,FINDCT_PG,   &
          M,MLX,MLY,MLZ,  &
          N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       & 
          MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
          NTIME,DT,D3,DCYL,GEOBAL,SCFACTH0,           &
          obsv,iflagobs,istate,coef_multivar_obs,leftsvd_obs,smean_obs,SnapNDT_obs,IFWIND,  &
          PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO,SMEAN_PG,leftsvd_PG_x,leftsvd_PG_y,  &
          NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM,DF,IF_REMESH_OBS,IF_REMESH_VORTICITY)
       
       use AllSorts
       use FLDebug
       use Shape_Transformations
       use Coordinates
       use coriolis_module
       use position_in_matrix
       use tr2d_module
       IMPLICIT NONE
       INTEGER, INTENT(IN)  ::NSVD_TOTAL,NSVD,NVAR,NSVD_U,NSVD_PHI
       INTEGER, INTENT(IN)  ::TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,MXNODS,fredop
       REAL, INTENT(out),DIMENSION(NONODS*9) ::FERROR
       REAL, INTENT(IN),DIMENSION(0:NTIME,NONODS) ::U0,V0,W0,PHI0
       REAL, INTENT(IN),DIMENSION(0:NTIME,NONODS) ::U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ
       REAL, INTENT(IN),DIMENSION(0:NTIME,NONODS) ::U0_OBS,V0_OBS,W0_OBS
       REAL, INTENT(IN),DIMENSION(NONODS) ::X,Y,Z
       REAL, INTENT(IN),DIMENSION(NONODS) :: MUPTXX,MUPTXY,MUPTXZ
       REAL, INTENT(IN),DIMENSION(NONODS) :: MUPTYY,MUPTYZ,MUPTZZ
       REAL,        DIMENSION(NVAR*MXNODS,NSVD) :: LEFTSVD
       REAL,        DIMENSION(NVAR*MXNODS) ::SMEAN
       !  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF
       REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_TOTAL) :: COEF_ADJ
       REAL, INTENT(IN)    ::DT
       INTEGER, INTENT(IN) ::NTIME
       LOGICAL, INTENT(IN) ::D3,DCYL
       REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: N,NLX,NLY,NLZ
       REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: M,MLX,MLY,MLZ
       REAL, INTENT(IN),DIMENSION(NGI)  ::WEIGHT
       INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
       INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC)::PNDGLN
       ! coriolos force
       INTEGER, INTENT(IN)  :: GEOBAL
       REAL, INTENT(IN)  :: SCFACTH0
       REAL tcoef(NSVD)
       REAL DF
       
       INTEGER NCOLM
       INTEGER FINDRM(NONODS+1),COLM(NCOLM)
       INTEGER  :: NCT
       INTEGER  :: FINDCT(FREDOP+1), COLCT(NCT)
       INTEGER  :: NCT_PG
       INTEGER MPGLOC,PGNODS
       INTEGER  :: FINDCT_PG(PGNODS+1), COLCT_PG(NCT_PG)
       ! Wind stress
       !------------
       LOGICAL IFWIND

       LOGICAL IF_REMESH_OBS,IF_REMESH_VORTICITY
       REAL,DIMENSION(:),ALLOCATABLE     :: VORTICITY_U, VORTICITY_V
       REAL,DIMENSION(:,:),ALLOCATABLE   :: VORTICITY_U_leftsvd,VORTICITY_V_leftsvd
       
       ! observational
       integer istate
       integer SnapNDT_obs
       integer,INTENT(IN),DIMENSION(0:ntime,mxnods) ::iflagobs
       REAL,INTENT(IN),DIMENSION(0:ntime,istate):: obsv
       REAL coef_multivar_obs(0:ntime,NSVD_TOTAL)
       REAL leftsvd_obs(ISTATE,nsvd)
       REAl soux_mean,souy_mean,souz_mean,soux1,souy1,souz1,soux2,souy2,souz2
       REAL SOUX,SOUY,SOUZ
       REAL smean_obs(istate)
       REAL,DIMENSION(:),ALLOCATABLE   :: Bmat_obs
       
       ! GEOGRAPHIC PRESSURE
       !--------------------
       
       REAL, INTENT(INOUT),DIMENSION(MPGLOC,NGI)  :: MPG,MPGLX,MPGLY,MPGLZ
       INTEGER, INTENT(IN),DIMENSION(TOTELE*MPGLOC):: PGNDGLNO
       !  REAL, INTENT(INOUT),DIMENSION(PGNODS,NSVD_PHI) :: LEFTSVD_PG
       REAL  ::SMEAN_PG(PGNODS)
       !  REAL, INTENT(INOUT),DIMENSION(0:ntime,NSVD_PHI) :: COEF_PG
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
       INTEGER ::ISVD,JSVD,II,JJ,KK,NSTEP,ELE,GI,ILOC,JLOC,I,J,K
       INTEGER ::IGL,JGL,IGLP,COUNT,COUNT1,JGLP
       REAL    ::RNN,RN1a,RN2a,RN3a,RN1b,RN2b,RN3b,FU,FV,U_x,V_y,W_z,UU_x,VU_y,WU_z,  &
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
            U_xx_mean,U_yy_mean,U_zz_mean,V_xx_mean,V_yy_mean,V_zz_mean,W_xx_mean,W_yy_mean,W_zz_mean,  &
            UDmean_UADJx,VDmean_UADJy,WDmean_UADJz,UDmean_VADJx,VDmean_VADJy,WDmean_VADJz, &
            UDmean_WADJx,VDmean_WADJy,WDmean_WADJz,                                        &
            UD_UADJx,VD_UADJy,WD_UADJz,UD_VADJx,VD_VADJy,WD_VADJz,UD_WADJx,VD_WADJy,WD_WADJz, &
            Umeanx_UADJ,Umeany_VADJ,Umeanz_WADJ,Vmeanx_UADJ,Vmeany_VADJ,                   &
            Vmeanz_WADJ,Wmeanx_UADJ,Wmeany_VADJ,Wmeanz_WADJ,                               &         
            Ux_UADJ,Uy_VADJ,Uz_WADJ,Vx_UADJ,Vy_VADJ,Vz_WADJ,Wx_UADJ,Wy_VADJ,Wz_WADJ
       
       
       INTEGER ::NCOLM1,NCOLM2
       INTEGER ::IBL11,IBL12,IBL13,IBL14,IBL15
       INTEGER ::IBL21,IBL22,IBL23,IBL24,IBL25
       INTEGER ::IBL31,IBL32,IBL33,IBL34,IBL35
       INTEGER ::IBL41,IBL42,IBL43,IBL44,IBL45
       INTEGER ::IBL51,IBL52,IBL53,IBL54,IBL55
       real    ::VOL,ftest,VOLUME
       INTEGER ::NT,ROW,ROWCNT
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
       
       REAL,DIMENSION(:),ALLOCATABLE     :: BIGM,VEC,VEC1,VEC_extr,BIGM_P,BIGM_PG
       REAL,DIMENSION(:,:),ALLOCATABLE   :: ReMat,ReMat1
       REAL,DIMENSION(:),ALLOCATABLE     :: BIGM0,BMat,BMat1,BMat0,BMat01,VEC0,VEC0_source,BMat0_source,BMat_extr
       REAL,DIMENSION(:),ALLOCATABLE     :: DETWEI,XD,YD,ZD
       REAL,DIMENSION(:),ALLOCATABLE     :: UD_mean,VD_mean,WD_mean
       REAL,DIMENSION(:),ALLOCATABLE     :: UD,VD,WD,UD_ADJ,VD_ADJ,WD_ADJ
       REAL,DIMENSION(:),ALLOCATABLE     :: UDx_mean,UDy_mean,UDz_mean,VDx_mean,VDy_mean,VDz_mean,WDx_mean,WDy_mean,WDz_mean
       REAL,DIMENSION(:),ALLOCATABLE     :: UDx,UDy,UDz,VDx,VDy,VDz,WDx,WDy,WDz
       REAL,DIMENSION(:),ALLOCATABLE     :: Umeanx,Umeany,Umeanz,Vmeanx,Vmeany,Vmeanz,Wmeanx,Wmeany,Wmeanz
       REAL,DIMENSION(:),ALLOCATABLE     :: WEI_u,WEI_v,WEI_w,WEI_p
       REAL,DIMENSION(:),ALLOCATABLE     :: WEI_ux,WEI_uy,WEI_uz,WEI_vx,WEI_vy,WEI_vz,WEI_wx,WEI_wy,WEI_wz
       REAL,DIMENSION(:,:),ALLOCATABLE   :: NX,NY,NZ,MX,MY,MZ
       REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF,NEWCOEF1
       REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF_ADJ,NEWCOEF1_ADJ
       REAL,DIMENSION(:),ALLOCATABLE     :: CORIOLIS
       
       REAL,DIMENSION(:,:),ALLOCATABLE   :: BIGM1,BIGM1consi
       REAL,DIMENSION(:),ALLOCATABLE     :: UDMat
       !  REAL,DIMENSION(:),ALLOCATABLE     :: FINDRM,COLM,CENTRM
       REAL,DIMENSION(:),ALLOCATABLE     :: ML,uu
       REAL,DIMENSION(:),ALLOCATABLE     :: SOURCX,SOURCY,SOURCZ
       REAL,DIMENSION(:,:),ALLOCATABLE   ::PHI0_ADJ_PG
       REAL varincr_pg(PGNODS)
       
       ! FOR HESSIAN AND ERROR MEASURE
       REAL              :: VEC_TOTAL(NONODS)
       REAL,DIMENSION(NONODS*9)        ::HESSIAN_U,HESSIAN_V,HESSIAN_W
       REAL U_ADJ_t,V_ADJ_t,W_ADJ_t,U_ADJ_t1,V_ADJ_t1,W_ADJ_t1
       
       write(3,*) '1'
       ALLOCATE(BIGM(3*NCOLM))
       ALLOCATE(BIGM_P(3*NCT))
       ALLOCATE(BIGM_PG(3*NCT_PG))
       write(3,*) '2'
       ALLOCATE(VEC(4*NONODS))
       ALLOCATE(DETWEI(NGI))
       ALLOCATE(UD(NGI))
       ALLOCATE(VD(NGI))
       ALLOCATE(WD(NGI))
       ALLOCATE(UDx(NGI))
       ALLOCATE(UDy(NGI))
       ALLOCATE(UDz(NGI))
       ALLOCATE(VDx(NGI))
       ALLOCATE(VDy(NGI))
       ALLOCATE(VDz(NGI))
       ALLOCATE(WDx(NGI))
       ALLOCATE(WDy(NGI))
       ALLOCATE(WDz(NGI))
       write(3,*) '3'
       ALLOCATE(UDx_mean(NGI))
       ALLOCATE(UDy_mean(NGI))
       ALLOCATE(UDz_mean(NGI))
       ALLOCATE(VDx_mean(NGI))
       ALLOCATE(VDy_mean(NGI))
       ALLOCATE(VDz_mean(NGI))
       ALLOCATE(WDx_mean(NGI))
       ALLOCATE(WDy_mean(NGI))
       ALLOCATE(WDz_mean(NGI))
       
       ALLOCATE(XD(NGI))
       ALLOCATE(YD(NGI))
       ALLOCATE(ZD(NGI))
       ALLOCATE(NX(NLOC,NGI))
       ALLOCATE(NY(NLOC,NGI))
       ALLOCATE(NZ(NLOC,NGI))
       ALLOCATE(MX(NLOC,NGI))
       ALLOCATE(MY(NLOC,NGI))
       ALLOCATE(MZ(NLOC,NGI))
       ALLOCATE(CORIOLIS(NGI))
       
       write(3,*) '4'
       
       ! for geographic pressure
       ALLOCATE(MPGX(MPGLOC,NGI))
       ALLOCATE(MPGY(MPGLOC,NGI))
       ALLOCATE(MPGZ(MPGLOC,NGI))
       write(3,*) ntime,PGNODS
       ALLOCATE(PHI0_ADJ_PG(0:ntime,PGNODS))
       !**************************************************!
       !       Initialise all matrices                    !
       !**************************************************!
       write(3,*) '5'
       VEC_TOTAL=0.0
       FERROR=0.0
       BIGM=0.0
       BIGM_P=0.0
       BIGM_PG=0.0
       VEC =0.0
       DETWEI=0.0
       UD=0.0
       VD=0.0
       WD=0.0
       XD=0.0
       YD=0.0
       ZD=0.0
       UDx=0.0
       UDy=0.0
       UDz=0.0              
       VDx=0.0
       VDy=0.0
       VDz=0.0              
       WDx=0.0
       WDy=0.0
       WDz=0.0
       
       CORIOLIS=0.0
       
       
       NX=0.0
       NY=0.0
       NZ=0.0
       MX=0.0
       MY=0.0
       MZ=0.0
       
       R(1:7*NSVD)=0.0
       
       ! for geostraphic pressure
       
       MPGLX=0.0
       MPGLY=0.0
       MPGLZ=0.0
       MPGX=0.0
       MPGY=0.0
       MPGZ=0.0
       PHI0_ADJ_PG=0.0
       MPGLOC=10
       CALL TRIQUA(L1, L2, L3, L4, PGWEIG, D3,NGI)
       ! Work out the shape functions and there derivatives...
       CALL SHATRI(L1, L2, L3, L4, PGWEIG, D3, &
            MPGLOC,NGI,                     &     
            MPG,MPGLX,MPGLY,MPGLZ) 
       
       
       THETA=0.5
       ONEMINTH=1.0-THETA
       
       IBL11=0
       IBL12=1*NCOLM
       IBL13=2*NCOLM
       IBL21=3*NCOLM
       IBL22=4*NCOLM
       IBL23=5*NCOLM
       IBL31=6*NCOLM
       IBL32=7*NCOLM
       IBL33=8*NCOLM
       

       ! calculate PHI0_PG....
       !-----------------------
       PHI0_ADJ_PG=0.0
       DO NSTEP = 0,ntime
          ewrite(3,*) 'NSTEP =',NSTEP
          !****************project back from rom to full space
          varincr_pg=0.0
          tcoef(:)=coef_adj(NSTEP,1:NSVD)
          call podtofullproj(PGNODS,nsvd_phi,varincr_pg,leftsvd_PG_x,tcoef)
          DO II=1,PGNODS        
             PHI0_ADJ_PG(NSTEP,II)=PHI0_ADJ_PG(NSTEP,II)+varincr_pg(II)+smean_pg(II)
          ENDDO
          varincr_pg=0.0
          tcoef(:)=coef_adj(NSTEP,NSVD+1:2*NSVD)
          call podtofullproj(PGNODS,nsvd_phi,varincr_pg,leftsvd_PG_y,tcoef)
          DO II=1,PGNODS        
             PHI0_ADJ_PG(NSTEP,II)=PHI0_ADJ_PG(NSTEP,II)+varincr_pg(II)
          ENDDO
       ENDDO
       
       
       ! Start to calculate the residual
       !---------------------------------
       !****************
       ! VORTICITY
       !****************
       IF(IF_REMESH_VORTICITY) THEN
          allocate(VORTICITY_U(nonods))
          allocate(VORTICITY_V(nonods))
          VORTICITY_U=0.0
          VORTICITY_V=0.0
       ENDIF
       
       DO 100 NSTEP=1,ntime
          !  DO 100 NSTEP=1,6
          write(3,*) 'NSTEP =',NSTEP
          ! setup the reduced model for BETA_ii(t)
          ! the original momentum equations are multiplied by the POD basic vectors 
          ! PHI_ii^u(x),PHI_ii^v(x),PHI_ii^w(x), and continuity equation is multiplied by PHI_ii^p(x)
          ! and then integrate them over the whole domain
          !  PHI_ii^p*{U_x+V_y+W_z}
          ! =PHI_ii{[dU_mean(X)/dx+COEF(t)*(d/dx)*SUM(PHI_ii^u)]+....}
          !------------------------------------------------------------------------------------------
          
          !****************
          ! VORTICITY
          !****************
          IF(IF_REMESH_VORTICITY) THEN
             VORTICITY_U=0.0
             VORTICITY_V=0.0
             call Vorticity2D(X,Y,Z, &
                  U0(nstep,1:nonods), V0(nstep,1:nonods), &
                  ndglno, totele, nonods, nloc, ngi, &
                  D3,DCYL, & 
                  VORTICITY_U, VORTICITY_V, &
                  N,NLX,NLY,NLZ, &
                  WEIGHT)
          ENDIF
          !*******************
          
          BIGM=0.0
          
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
             
             
             DO 20 ILOC=1,NLOC
                
                IGL=NDGLNO((ELE-1)*NLOC+ILOC)
                
                !PHI_jj*U=PHI_jj*[U_mean+SUM_jj=1^nsvd(COEF*PHI_ii)]=PHI_jj*U_mean+COEF(t,jj)
                ! UD(GI)=COEF(t,jj)+phi_jj*[N(GI)*U(X)_mean]  
                ! VD(GI)=COEF(t,jj)+phi_jj*[N(GI)*V(X)_mean]  
                ! WD(GI)=COEF(t,jj)+phi_jj*[N(GI)*W(X)_mean]  
           !-----------------------------------------------
                
                DO GI=1,NGI
                   UD(GI)=0.0
                   VD(GI)=0.0
                   WD(GI)=0.0
                   
                   UDx(GI)=0.0
                   UDy(GI)=0.0
                   UDz(GI)=0.0
                   
                   VDx(GI)=0.0
                   VDy(GI)=0.0
                   VDz(GI)=0.0
                   
                   WDx(GI)=0.0
                   WDy(GI)=0.0
                   WDz(GI)=0.0
                   CORIOLIS(GI)=0.0
                   
                   DO JLOC=1,NLOC
                      JGL=NDGLNO((ELE-1)*NLOC+JLOC)
                      
                      UD(GI)=UD(GI)+U0(NTIME-NSTEP,JGL)*N(JLOC,GI)
                      VD(GI)=VD(GI)+V0(NTIME-NSTEP,JGL)*N(JLOC,GI)
                      WD(GI)=WD(GI)+W0(NTIME-NSTEP,JGL)*N(JLOC,GI)
                      
                      UDx(GI)=UDx(GI)+U0(NTIME-NSTEP,JGL)*NX(JLOC,GI)
                      UDy(GI)=UDy(GI)+U0(NTIME-NSTEP,JGL)*NY(JLOC,GI)
                      UDz(GI)=UDz(GI)+U0(NTIME-NSTEP,JGL)*NZ(JLOC,GI)
                      
                      VDx(GI)=VDx(GI)+V0(NTIME-NSTEP,JGL)*NX(JLOC,GI)
                      VDy(GI)=VDy(GI)+V0(NTIME-NSTEP,JGL)*NY(JLOC,GI)
                      VDz(GI)=VDz(GI)+V0(NTIME-NSTEP,JGL)*NZ(JLOC,GI)
                      
                      WDx(GI)=WDx(GI)+W0(NTIME-NSTEP,JGL)*NX(JLOC,GI)
                      WDy(GI)=WDy(GI)+W0(NTIME-NSTEP,JGL)*NY(JLOC,GI)
                      WDz(GI)=WDz(GI)+W0(NTIME-NSTEP,JGL)*NZ(JLOC,GI)
                      
                   ENDDO
                   CORIOLIS(GI) = 2.0*FUNOME(XD(GI),YD(GI),ZD(GI))
                ENDDO
                
                
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
                DO JLOC=1,NLOC
                   JGL=NDGLNO((ELE-1)*NLOC+JLOC)
                   IF(IGL.EQ.JGL) THEN
                      U_ADJ_t =0.0
                      V_ADJ_t =0.0
                      W_ADJ_t =0.0
                      U_ADJ_t1 =0.0
                      V_ADJ_t1 =0.0
                      W_ADJ_t1 =0.0 
                      K =0
                      DO II=findrm(igl),findrm(igl+1)-1
                         K=k+1
                         U_ADJ_t = U_ADJ_t+ U0_ADJ(NSTEP,COLM(II))
                         V_ADJ_t = V_ADJ_t+ V0_ADJ(NSTEP,COLM(II))
                         W_ADJ_t = W_ADJ_t+ W0_ADJ(NSTEP,COLM(II))
                         U_ADJ_t1 = U_ADJ_t1+ U0_ADJ(NSTEP-1,COLM(II))
                         V_ADJ_t1 = V_ADJ_t1+ V0_ADJ(NSTEP-1,COLM(II))
                         W_ADJ_t1 = W_ADJ_t1+ W0_ADJ(NSTEP-1,COLM(II))
                      ENDDO
                      U_ADJ_t = U_ADJ_t/k
                      V_ADJ_t = V_ADJ_t/k
                      W_ADJ_t = W_ADJ_t/k
                      U_ADJ_t1 = U_ADJ_t1/k
                      V_ADJ_t1 = V_ADJ_t1/k
                      W_ADJ_t1 = W_ADJ_t1/k
                      
                   ELSE
                      U_ADJ_t=U0_ADJ(NSTEP,JGL)
                      V_ADJ_t=V0_ADJ(NSTEP,JGL)
                      W_ADJ_t=W0_ADJ(NSTEP,JGL)
                      U_ADJ_t1=U0_ADJ(NSTEP-1,JGL)
                      V_ADJ_t1=V0_ADJ(NSTEP-1,JGL)
                      W_ADJ_t1=W0_ADJ(NSTEP-1,JGL)
                   ENDIF
                   
                   
                   RN1a=0.0
                   RN2a=0.0
                   RN3a=0.0
                   RN1b=0.0
                   RN2b=0.0
                   RN3b=0.0
                   U_x=0.0
                   V_y=0.0
                   W_Z=0.0
                   
                   UD_Ux=0.0
                   VD_Uy=0.0
                   WD_Uz=0.0
                   UD_Vx=0.0
                   VD_Vy=0.0
                   WD_Vz=0.0
                   UD_Wx=0.0
                   VD_Wy=0.0
                   WD_Wz=0.0
                   
                   Ux_UADJ=0.0
                   Uy_VADJ=0.0
                   Uz_WADJ=0.0
                   Vx_UADJ=0.0
                   Vy_VADJ=0.0
                   Vz_WADJ=0.0
                   Wx_UADJ=0.0
                   Wy_VADJ=0.0
                   Wz_WADJ=0.0
                   
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
                   PG_x=0.0
                   PG_y=0.0
                   PG_z=0.0
                   soux =0.0
                   souy =0.0
                   souz =0.0
                   
                   DO GI=1,NGI
                      RN1a=RN1a+ N(ILOC,GI)*DETWEI(GI)*U_ADJ_t1*N(JLOC,GI)
                      RN2a=RN2a+ N(ILOC,GI)*DETWEI(GI)*V_ADJ_t1*N(JLOC,GI)
                      RN3a=RN3a+ N(ILOC,GI)*DETWEI(GI)*W_ADJ_t1*N(JLOC,GI)
                      
                      RN1b=RN1b+ N(ILOC,GI)*DETWEI(GI)*U_ADJ_t*N(JLOC,GI)
                      RN2b=RN2b+ N(ILOC,GI)*DETWEI(GI)*V_ADJ_t*N(JLOC,GI)
                      RN3b=RN3b+ N(ILOC,GI)*DETWEI(GI)*W_ADJ_t*N(JLOC,GI)
                      
                      U_x=U_x+ N(ILOC,GI)*DETWEI(GI)*U_ADJ_t*NX(JLOC,GI)
                      V_y=V_y+ N(ILOC,GI)*DETWEI(GI)*V_ADJ_t*NY(JLOC,GI)
                      W_Z=W_Z+ N(ILOC,GI)*DETWEI(GI)*W_ADJ_t*NZ(JLOC,GI)
                      
                      UD_Ux=UD_Ux+ N(ILOC,GI)*UD(GI)*DETWEI(GI)*U_ADJ_t*NX(JLOC,GI)
                      VD_Uy=VD_Uy+ N(ILOC,GI)*VD(GI)*DETWEI(GI)*U_ADJ_t*NY(JLOC,GI)
                      WD_Uz=WD_Uz+ N(ILOC,GI)*WD(GI)*DETWEI(GI)*U_ADJ_t*NZ(JLOC,GI)
                      UD_Vx=UD_Vx+ N(ILOC,GI)*UD(GI)*DETWEI(GI)*V_ADJ_t*NX(JLOC,GI)
                      VD_Vy=VD_Vy+ N(ILOC,GI)*VD(GI)*DETWEI(GI)*V_ADJ_t*NY(JLOC,GI)
                      WD_Vz=WD_Vz+ N(ILOC,GI)*WD(GI)*DETWEI(GI)*V_ADJ_t*NZ(JLOC,GI)
                      UD_Wx=UD_Wx+ N(ILOC,GI)*UD(GI)*DETWEI(GI)*W_ADJ_t*NX(JLOC,GI)
                      VD_Wy=VD_Wy+ N(ILOC,GI)*VD(GI)*DETWEI(GI)*W_ADJ_t*NY(JLOC,GI)
                      WD_Wz=WD_Wz+ N(ILOC,GI)*WD(GI)*DETWEI(GI)*W_ADJ_t*NZ(JLOC,GI)
                      
                      
                    Ux_UADJ=Ux_UADJ-N(ILOC,GI)*UDx(GI)*DETWEI(GI)*U_ADJ_t*N(JLOC,GI)
                    Uy_VADJ=Uy_VADJ-N(ILOC,GI)*UDy(GI)*DETWEI(GI)*V_ADJ_t*N(JLOC,GI)
                    Uz_WADJ=Uz_WADJ-N(ILOC,GI)*UDz(GI)*DETWEI(GI)*W_ADJ_t*N(JLOC,GI)
                    Vx_UADJ=Vx_UADJ-N(ILOC,GI)*VDx(GI)*DETWEI(GI)*U_ADJ_t*N(JLOC,GI)
                    Vy_VADJ=Vy_VADJ-N(ILOC,GI)*VDy(GI)*DETWEI(GI)*V_ADJ_t*N(JLOC,GI)
                    Vz_WADJ=Vz_WADJ-N(ILOC,GI)*VDz(GI)*DETWEI(GI)*W_ADJ_t*N(JLOC,GI)
                    Wx_UADJ=Wx_UADJ-N(ILOC,GI)*WDx(GI)*DETWEI(GI)*U_ADJ_t*N(JLOC,GI)
                    Wy_VADJ=Wy_VADJ-N(ILOC,GI)*WDy(GI)*DETWEI(GI)*V_ADJ_t*N(JLOC,GI)
                    Wz_WADJ=Wz_WADJ-N(ILOC,GI)*WDz(GI)*DETWEI(GI)*W_ADJ_t*N(JLOC,GI)
                    
                    U_xx=U_xx+NX(ILOC,GI)*MUPTXX(JGL)*DETWEI(GI)*U_ADJ_t*NX(JLOC,GI)
                    U_yy=U_yy+NY(ILOC,GI)*MUPTYY(JGL)*DETWEI(GI)*U_ADJ_t*NY(JLOC,GI)
                    U_zz=U_zz+NZ(ILOC,GI)*MUPTZZ(JGL)*DETWEI(GI)*U_ADJ_t*NZ(JLOC,GI)
                    V_xx=V_xx+NX(ILOC,GI)*MUPTXX(JGL)*DETWEI(GI)*V_ADJ_t*NX(JLOC,GI)
                    V_yy=V_yy+NY(ILOC,GI)*MUPTYY(JGL)*DETWEI(GI)*V_ADJ_t*NY(JLOC,GI)
                    V_zz=V_zz+NZ(ILOC,GI)*MUPTZZ(JGL)*DETWEI(GI)*V_ADJ_t*NZ(JLOC,GI)
                    W_xx=W_xx+NX(ILOC,GI)*MUPTXX(JGL)*DETWEI(GI)*W_ADJ_t*NX(JLOC,GI)
                    W_yy=W_yy+NY(ILOC,GI)*MUPTYY(JGL)*DETWEI(GI)*W_ADJ_t*NY(JLOC,GI)
                    W_zz=W_zz+NZ(ILOC,GI)*MUPTZZ(JGL)*DETWEI(GI)*W_ADJ_t*NZ(JLOC,GI)

                    FU = FU+N(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*U_ADJ_t*N(JLOC,GI)
                    FV = FV-N(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*V_ADJ_t*N(JLOC,GI)

                    ! Source terms which represet the misfit between the numerical results and observations
                    
                    
                    !                    IF(IF_REMESH_VORTICITY) THEN
                    !                      SOUX = SOUX +N(ILOC,GI)*DETWEI(GI)*VORTICITY_U(JGL)*N(JLOC,GI)
                    !                      SOUY = SOUY +N(ILOC,GI)*DETWEI(GI)*VORTICITY_V(JGL)*N(JLOC,GI)
                    !                    ELSE IF(IF_REMESH_OBS) THEN
                    !                      soux =soux+N(ILOC,GI)*DETWEI(GI)*(U0(NTIME-NSTEP,JGL)-U0_OBS(NTIME-NSTEP,JGL))*N(JLOC,GI)
                    !                      souy =souy+N(ILOC,GI)*DETWEI(GI)*(V0(NTIME-NSTEP,JGL)-V0_OBS(NTIME-NSTEP,JGL))*N(JLOC,GI)
                    !                      souz =souz+N(ILOC,GI)*DETWEI(GI)*(W0(NTIME-NSTEP,JGL)-W0_OBS(NTIME-NSTEP,JGL))*N(JLOC,GI)
                    !                    ELSE
                    !                      soux =soux+N(ILOC,GI)*DETWEI(GI)*(U0(NTIME-NSTEP,JGL))*N(JLOC,GI)
                    !                      souy =souy+N(ILOC,GI)*DETWEI(GI)*(V0(NTIME-NSTEP,JGL))*N(JLOC,GI)
                    !                      souz =souz+N(ILOC,GI)*DETWEI(GI)*(W0(NTIME-NSTEP,JGL))*N(JLOC,GI)
                    !                   ENDIF
                    
                    
                 ENDDO
                 
                 LOWER=FINDRM(IGL) 
                 UPPER=FINDRM(IGL+1)-1
7000             CONTINUE
                 INUM=LOWER+(UPPER-LOWER+1)/2 
                 IF(JGL.GE.COLM(INUM) )  THEN 
                    LOWER=INUM
                 ELSE
                    UPPER=INUM
                 ENDIF
                 IF(UPPER-LOWER.LE.1) THEN
                    IF(JGL.EQ.COLM(LOWER)) THEN
                       COUNT=LOWER
                    ELSE
                       COUNT=UPPER
                    ENDIF
                    GOTO 9000
                 ENDIF
                 GOTO 7000
9000             CONTINUE
                 
                 BIGM(COUNT)=BIGM(COUNT)+RN1b-RN1a+DT*(UD_Ux+VD_Uy+WD_Uz+U_xx+U_yy+U_zz+FV+  &
                      Ux_UADJ+Uy_VADJ+Uz_WADJ)-DT*soux
                 BIGM(NCOLM+COUNT)=BIGM(NCOLM+COUNT)+RN2b-RN2a+DT*(UD_Vx+VD_Vy+WD_Vz+V_xx+V_yy+V_zz+FU + &
                      Vx_UADJ+Vy_VADJ+Vz_WADJ)-DT*souy
                 BIGM(2*NCOLM+COUNT)=BIGM(2*NCOLM+COUNT)+RN3b-RN3a+DT*(UD_Wx+VD_Wy+WD_Wz+W_xx+W_yy+W_zz+ &
                      Wx_UADJ+Wy_VADJ+Wz_WADJ)-DT*souz
              ENDDO
              
              DO JLOC=1,MLOC
                 JGLP=PNDGLN((ELE-1)*MLOC+JLOC)
                 P_x=0.0
                 P_y=0.0
                 P_z=0.0
                 DO GI=1,NGI
                    P_x=P_x+N(ILOC,GI)*DETWEI(GI)*PHI0_ADJ(NSTEP,JGLP)*MX(JLOC,GI)
                    P_y=P_y+N(ILOC,GI)*DETWEI(GI)*PHI0_ADJ(NSTEP,JGLP)*MY(JLOC,GI)
                    P_z=P_z+N(ILOC,GI)*DETWEI(GI)*PHI0_ADJ(NSTEP,JGLP)*MZ(JLOC,GI)
                 ENDDO
                 CALL POSINMAT(COUNT,JGLP,IGL,FREDOP,FINDCT,COLCT,NCT)
                 BIGM_P(COUNT)=BIGM_P(COUNT)+DT*P_x
                 BIGM_P(NCT+COUNT)=BIGM_P(NCT+COUNT)+DT*P_y
                 BIGM_P(2*NCT+COUNT)=BIGM_P(2*NCT+COUNT)+DT*P_z
              ENDDO
              
              DO JLOC=1,MPGLOC
                 JGLPG=PGNDGLNO((ELE-1)*MPGLOC+JLOC)
                 if(JGLPG.GT.PGNODS.OR.JGLPG.LE.0) THEN
                    write(3,*) PGNDGLNO((ELE-1)*MPGLOC+JLOC),ELE,MPGLOC,JLOC,JGLPG
                    stop 344
                 endif
                 PG_x=0.0
                 PG_y=0.0
                 PG_z=0.0
                 DO GI=1,NGI
                    PG_x =PG_x + N(ILOC,GI)*DETWEI(GI)*PHI0_ADJ_PG(NSTEP,JGLPG)*MPGX(JLOC,GI)
                    PG_y =PG_y + N(ILOC,GI)*DETWEI(GI)*PHI0_ADJ_PG(NSTEP,JGLPG)*MPGY(JLOC,GI)
                    PG_z =PG_z + N(ILOC,GI)*DETWEI(GI)*PHI0_ADJ_PG(NSTEP,JGLPG)*MPGZ(JLOC,GI)
                 ENDDO
                 CALL POSINMAT(COUNT,JGLPG,IGL,PGNODS,FINDCT_PG,COLCT_PG,NCT_PG)
                 BIGM_PG(COUNT)=BIGM_PG(COUNT)+DT*PG_x
                 BIGM_PG(NCT_PG+COUNT)=BIGM_PG(NCT_PG+COUNT)+DT*PG_y
                 BIGM_PG(2*NCT_PG+COUNT)=BIGM_PG(2*NCT_PG+COUNT)+DT*PG_z
              ENDDO
              
              
20            CONTINUE
              !              write(3,*) '1',BIGM(COUNT)
              
              
30            CONTINUE
              

!              BIGM=ABS(BIGM)
              !              BIGM_P=ABS(BIGM_P)
!              BIGM_PG=ABS(BIGM_PG)
              VEC=0.0
              DO II=1,NONODS
                 
                 LOWER=FINDRM(II)
                 UPPER=FINDRM(II+1)-1
                 DO KK=LOWER,UPPER
                    VEC(II)=VEC(II)+BIGM( KK )
                    VEC(NONODS+II)=VEC(NONODS+II)+BIGM( NCOLM+KK )
                    VEC(2*NONODS+II)=VEC(2*NONODS+II)+BIGM( 2*NCOLM+KK )
                 ENDDO
                 
              ENDDO
              DO JJ=1,FREDOP
                 
                 LOWER=FINDCT(JJ)
                 UPPER=FINDCT(JJ+1)-1
                 DO KK=LOWER,UPPER
                    II=COLCT(KK)
                    VEC(II)=VEC(II)+BIGM_P(KK)
                    VEC(NONODS+II)=VEC(NONODS+II)+BIGM_P( NCT+KK )
                    VEC(2*NONODS+II)=VEC(2*NONODS+II)+BIGM_P( 2*NCT+KK )
                 ENDDO
              ENDDO

              DO JJ=1,PGNODS
                 
                 LOWER=FINDCT_PG(JJ)
                 UPPER=FINDCT_PG(JJ+1)-1
                 DO KK=LOWER,UPPER
                    II=COLCT_PG(KK)
                    VEC(II)=VEC(II)+BIGM_PG(KK)
                    VEC(NONODS+II)=VEC(NONODS+II)+BIGM_PG(NCT_PG+KK )
                    VEC(2*NONODS+II)=VEC(2*NONODS+II)+BIGM_PG( 2*NCT_PG+KK )
                 ENDDO
              ENDDO
              
              VEC = ABS(VEC)
              VEC_TOTAL =0.0
              DO II=1,NONODS
                 VEC_TOTAL(II)=VEC_TOTAL(II)+VEC(II)+VEC(NONODS+II)+VEC(2*NONODS+II)
                 !                write(3,*) II,VEC_TOTAL(II),VEC(II),VEC(NONODS+II),VEC(2*NONODS+II)
              ENDDO
              open(1,file='vec.dat')
              write(1,*) (VEC_TOTAL(ii),ii=1,NONODS)
              write(1,*) (VEC(ii),ii=1,3*NONODS)
              close(1)
              HESSIAN_U=0.0
              HESSIAN_V=0.0
              HESSIAN_W=0.0
              allocate(uu(nonods))
              uu(1:nonods)=U0(NSTEP,1:nonods)
              CALL CalHessian(NONODS,XNONOD,NLOC,NGI,TOTELE,3,UU,HESSIAN_U,         &
                   X,Y,Z,XONDGL,NDGLNO,N,NLX,NLY,NLZ,DCYL,D3,WEIGHT)
              uu(1:nonods)=V0(NSTEP,1:nonods)
              CALL CalHessian(NONODS,XNONOD,NLOC,NGI,TOTELE,3,UU,HESSIAN_V,         &
                   X,Y,Z,XONDGL,NDGLNO,N,NLX,NLY,NLZ,DCYL,D3,WEIGHT)
              uu(1:nonods)=W0(NSTEP,1:nonods)
              CALL CalHessian(NONODS,XNONOD,NLOC,NGI,TOTELE,3,UU,HESSIAN_W,         &
                   X,Y,Z,XONDGL,NDGLNO,N,NLX,NLY,NLZ,DCYL,D3,WEIGHT)
              deallocate(uu)
              DO II=1,NONODS
                 DO JJ=1,9
                    FERROR((II-1)*9+JJ)=FERROR((II-1)*9+JJ)+ABS(HESSIAN_U((II-1)*9+JJ))*ABS(VEC(II))
                    FERROR((II-1)*9+JJ)=FERROR((II-1)*9+JJ)+ABS(HESSIAN_V((II-1)*9+JJ))*ABS(VEC(NONODS+II))
                    FERROR((II-1)*9+JJ)=FERROR((II-1)*9+JJ)+ABS(HESSIAN_W((II-1)*9+JJ))*ABS(VEC(2*NONODS+II))
                 ENDDO
              ENDDO
              
100           CONTINUE
              
              DO II=1,NONODS
                 DO JJ=1,9
                    !                    FERROR((II-1)*9+JJ)=FERROR((II-1)*9+JJ)/ABS(VEC_TOTAL(II))
                    FERROR((II-1)*9+JJ)=FERROR((II-1)*9+JJ)/DF
                 ENDDO
              ENDDO
              
              
              
              !-------------------------------!
              ! DEALLOCATE THE LOCAL MATRICES !
              !-------------------------------!
              write(3,*) 'deallocate and end SUBROUTINE REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ_errormeasure'
              DEALLOCATE(BIGM)
              DEALLOCATE(BIGM_P)
              DEALLOCATE(BIGM_PG)
  DEALLOCATE(VEC)
  DEALLOCATE(DETWEI)
  DEALLOCATE(UD)
  DEALLOCATE(VD)
  DEALLOCATE(WD)
  DEALLOCATE(UDx)
  DEALLOCATE(UDy)
  DEALLOCATE(UDz)
  DEALLOCATE(VDx)
  DEALLOCATE(VDy)
  DEALLOCATE(VDz)
  DEALLOCATE(WDx)
  DEALLOCATE(WDy)
  DEALLOCATE(WDz)
  
  DEALLOCATE(XD)
  DEALLOCATE(YD)
  DEALLOCATE(ZD)
  DEALLOCATE(NX)
  DEALLOCATE(NY)
  DEALLOCATE(NZ)
  DEALLOCATE(MX)
  DEALLOCATE(MY)
  DEALLOCATE(MZ)
  DEALLOCATE(CORIOLIS)

  ! for geographic pressure
  DEALLOCATE(MPGX)
  DEALLOCATE(MPGY)
  DEALLOCATE(MPGZ)

  DEALLOCATE(PHI0_ADJ_PG)
  IF(IF_REMESH_VORTICITY) THEN
     deallocate(VORTICITY_U)
     deallocate(VORTICITY_V)
  ENDIF
  
  RETURN
END SUBROUTINE REDUCED_MODEL_DIFFCOEF1consistent_PG_ADJ_errormeasure


