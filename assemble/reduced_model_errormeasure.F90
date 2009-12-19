#include "fdebug.h"

module reduced_model_error_measure

  use allsorts
  use coordinates
  use fldebug
  use position_in_matrix
  use shape_functions
  use tr2d_module
  use coriolis_module
  
  implicit none
  
  private
  
  !public :: 
  
contains

SUBROUTINE REDUCED_MODEL_DIFFCOEF1_PG1_errormeasure(NSVD_TOTAL,NVAR,NSVD,NSVD_U,NSVD_PHI,MXNODS,NONODS,XNONOD,fredop,  &
     X,Y,Z,            &
     FERROR,COEF,U0,V0,W0,PHI0,U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ, &
     NCOLM,COLM,FINDRM,NCT,COLCT,FINDCT,NCT_PG,COLCT_PG,FINDCT_PG,  &
     M,MLX,MLY,MLZ,                                    &
     N,NLX,NLY,NLZ,NDGLNO,XONDGL,PNDGLN,WEIGHT,TOTELE,NLOC,NGI,MLOC,       &
     MUPTXX,MUPTXY,MUPTXZ,MUPTYY,MUPTYZ,MUPTZZ,                            &
     NTIME1,DT1,D3,DCYL,GEOBAL,SCFACTH0,IFWIND, &
     PGNODS,MPGLOC,MPG,MPGLX,MPGLY,MPGLZ,PGNDGLNO,SMEAN_PG,leftsvd_PG_x,leftsvd_PG_y,  &
     NPGCOLM,PGFINDRM,PGCOLM,PGCENTRM,DF)

  INTEGER, INTENT(IN)  ::NSVD_TOTAL,NSVD,NVAR,NSVD_U,NSVD_PHI
  INTEGER, INTENT(IN)  ::TOTELE,NONODS,XNONOD,NLOC,MLOC,NGI,MXNODS,fredop
  INTEGER, INTENT(IN) ::NTIME1
  INTEGER MPGLOC,PGNODS
  REAL, INTENT(OUT)  :: FERROR(NONODS*9)
  REAL, INTENT(IN)  :: u0(0:ntime1,nonods)
  REAL, INTENT(IN)  :: v0(0:ntime1,nonods)
  REAL, INTENT(IN)  :: w0(0:ntime1,nonods)
  REAL, INTENT(IN)  :: phi0(0:ntime1,nonods)
  REAL, INTENT(IN),DIMENSION(0:NTIME1,NONODS) ::U0_ADJ,V0_ADJ,W0_ADJ,PHI0_ADJ
  
  
  REAL, INTENT(INOUT),DIMENSION(NONODS) ::X,Y,Z
  REAL, INTENT(INOUT),DIMENSION(NONODS) :: MUPTXX,MUPTXY,MUPTXZ
  REAL, INTENT(INOUT),DIMENSION(NONODS) :: MUPTYY,MUPTYZ,MUPTZZ
  REAL, DIMENSION(NVAR*MXNODS,NSVD) :: LEFTSVD
  REAL, DIMENSION(NVAR*MXNODS) ::SMEAN
  REAL, DIMENSION(0:ntime1,NSVD_TOTAL) :: COEF
  REAL, INTENT(IN)    ::DT1
  LOGICAL, INTENT(IN) ::D3,DCYL
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: N,NLX,NLY,NLZ
  REAL, INTENT(IN),DIMENSION(NLOC,NGI)  :: M,MLX,MLY,MLZ
  REAL, INTENT(IN),DIMENSION(NGI)  ::WEIGHT
  INTEGER, INTENT(IN),DIMENSION(TOTELE*NLOC):: NDGLNO,XONDGL
  INTEGER, INTENT(IN),DIMENSION(TOTELE*MLOC)::PNDGLN
  INTEGER  :: NCT
  INTEGER  :: FINDCT(FREDOP+1), COLCT(NCT)
  INTEGER  :: NCT_PG
  INTEGER  :: FINDCT_PG(PGNODS+1), COLCT_PG(NCT_PG)
  ! coriolos force
  INTEGER, INTENT(IN)  :: GEOBAL
  REAL, INTENT(IN)  :: SCFACTH0
  REAL tcoef(NSVD)
  REAL DF
  ! Wind stress
  !------------
  LOGICAL IFWIND

  ! GEOGRAPHIC PRESSURE
  !--------------------

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
  REAL,DIMENSION(:,:),ALLOCATABLE   ::PHI0_PG
  REAL varincr_pg(PGNODS)

  !       Local variables
  REAL    ::THETA,ONEMINTH
  INTEGER ::ISVD,JSVD,II,JJ,KK,NSTEP1,ELE,GI,ILOC,JLOC,I,J,K,GLOBI,GLOBJ
  INTEGER ::IGL,JGL,IGLP,COUNT,COUNT1,JGLP
  REAL    ::RNN,RN1a,RN2a,RN3a,RN1b,RN2b,RN3b,FU,FV,U_x,V_y,W_z,UU_x,VU_y,WU_z,      &
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
  REAL          ::SOUX,SOUY,SOUZ
  INTEGER ::NSTEP,NTIME,NT,ROW,ROWCNT
  REAL    ::INFIN1,INFIN2,NEAR
  PARAMETER(INFIN1=10.0E+25)
  PARAMETER(INFIN2=10.0E+25)
  PARAMETER(NEAR=1./INFIN1)
  
  INTEGER        :: FINDRM(NONODS+1),COLM(NCOLM)
  INTEGER        :: CENTRM(NONODS) 
  
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
  
  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM,VEC,VEC1,BIGM_P,BIGM_PG
  REAL,DIMENSION(:,:),ALLOCATABLE   :: ReMat,ReMat1
  REAL,DIMENSION(:),ALLOCATABLE     :: BIGM0,BMat,BMat1,BMat0,BMat01,VEC0
  REAL,DIMENSION(:),ALLOCATABLE     :: DETWEI,XD,YD,ZD
  REAL,DIMENSION(:),ALLOCATABLE     :: UD_mean,VD_mean,WD_mean
  REAL,DIMENSION(:),ALLOCATABLE            :: UD,VD,WD
  REAL,DIMENSION(:),ALLOCATABLE     :: Umeanx,Umeany,Umeanz,Vmeanx,Vmeany,Vmeanz,Wmeanx,Wmeany,Wmeanz
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_u,WEI_v,WEI_w,WEI_p
  REAL,DIMENSION(:),ALLOCATABLE     :: WEI_ux,WEI_uy,WEI_uz,WEI_vx,WEI_vy,WEI_vz,WEI_wx,WEI_wy,WEI_wz
  REAL,DIMENSION(:,:),ALLOCATABLE   :: NX,NY,NZ,MX,MY,MZ
  REAL,DIMENSION(:),ALLOCATABLE     :: NEWCOEF,NEWCOEF1,NEWCOEF2
  REAL,DIMENSION(:),ALLOCATABLE     :: CORIOLIS
  
  REAL,DIMENSION(:,:),ALLOCATABLE   :: BIGM1
  REAL,DIMENSION(:),ALLOCATABLE     :: UDMat
  !  REAL,DIMENSION(:),ALLOCATABLE     :: FINDRM,COLM,CENTRM
  REAL,DIMENSION(:),ALLOCATABLE     :: ML,UU
  
  REAL,DIMENSION(:),ALLOCATABLE     :: SOURCX,SOURCY,SOURCZ
  
  ! FOR HESSIAN AND ERROR MEASURE
  
  ! FOR HESSIAN AND ERROR MEASURE
  REAL              :: VEC_TOTAL(NONODS)
  REAL,DIMENSION(NONODS*9)        ::HESSIAN_U,HESSIAN_V,HESSIAN_W
  
  REAL, DIMENSION(NONODS)        :: T
  REAL,DIMENSION(NONODS)        :: TX,TY,TZ,TXX,TXY,TXZ,TYY,TYZ,TZZ
  REAL U_t,V_t,W_t,U_t1,V_t1,W_t1
  
  
  write(3,*) 'start to run the reduced model_pg'
  NT=1
  DT=DT1/NT
  NTIME=NT*NTIME1
  
  ALLOCATE(BIGM(3*NCOLM))
  ALLOCATE(BIGM_P(3*NCT))
  ALLOCATE(BIGM_PG(3*NCT_PG))
  ALLOCATE(VEC(3*NONODS))
  ALLOCATE(DETWEI(NGI))
  ALLOCATE(UD(NGI))
  ALLOCATE(VD(NGI))
  ALLOCATE(WD(NGI))
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

  ALLOCATE(SOURCX(NONODS))
  ALLOCATE(SOURCY(NONODS))

  ! for geographic pressure
  ALLOCATE(MPGX(MPGLOC,NGI))
  ALLOCATE(MPGY(MPGLOC,NGI))
  ALLOCATE(MPGZ(MPGLOC,NGI))
  ALLOCATE(PHI0_PG(0:ntime1,PGNODS))
!**************************************************!
!       Initialise all matrices                    !
!**************************************************!
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
     CALL WINDYSOURCE(SOURCX,SOURCY,X,Y,NONODS,1)
  endif
  
  THETA=0.5
  ONEMINTH=1.0-THETA
  
  ! for geostraphic pressure

  MPGLX=0.0
  MPGLY=0.0
  MPGLZ=0.0
  MPGX=0.0
  MPGY=0.0
  MPGZ=0.0
  PHI0_PG=0.0
  
  MPGLOC=10
  CALL TRIQUA(L1, L2, L3, L4, PGWEIG, D3,NGI)
  ! Work out the shape functions and there derivatives...
  CALL SHATRI(L1, L2, L3, L4, PGWEIG, D3, &
       MPGLOC,NGI,                     &     
       MPG,MPGLX,MPGLY,MPGLZ) 
  
  ! calculate PHI0_PG....
  !-----------------------
  PHI0_PG=0.0
  DO NSTEP = 0,ntime1
     !     ewrite(3,*) 'NSTEP =',NSTEP
     !****************project back from rom to full space
     varincr_pg=0.0
     tcoef(:)=coef(NSTEP,1:NSVD)
     call podtofullproj(PGNODS,nsvd_phi,varincr_pg,leftsvd_PG_x,tcoef)
     DO II=1,PGNODS        
        PHI0_PG(NSTEP,II)=PHI0_PG(NSTEP,II)+varincr_pg(II)+smean_pg(II)
     ENDDO
     varincr_pg=0.0
     tcoef(:)=coef(NSTEP,NSVD+1:2*NSVD)
     call podtofullproj(PGNODS,nsvd_phi,varincr_pg,leftsvd_PG_y,tcoef)
     DO II=1,PGNODS        
        PHI0_PG(NSTEP,II)=PHI0_PG(NSTEP,II)+varincr_pg(II)
     ENDDO
  ENDDO

  
  ! Start to calculate the residual
  !---------------------------------
  DO 100 NSTEP=1,NTIME1
     write(3,*) 'nstep', nstep
     ewrite(3,*) 'nstep',nstep
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
        !     ewrite(3,*) 'ELE=',ELE
        
        
        
        !PHI_jj*U=PHI_jj*[U_mean+SUM_jj=1^nsvd(COEF*PHI_ii)]=PHI_jj*U_mean+COEF(t,jj)
        ! UD(GI)=COEF(t,jj)+phi_jj*[N(GI)*U(X)_mean]  
        ! VD(GI)=COEF(t,jj)+phi_jj*[N(GI)*V(X)_mean]  
        ! WD(GI)=COEF(t,jj)+phi_jj*[N(GI)*W(X)_mean]  
        !-----------------------------------------------
        DO GI=1,NGI
           UD(GI)=0.0
           VD(GI)=0.0
           WD(GI)=0.0
           CORIOLIS(GI)=0.0
           
           
           DO JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
              
              UD(GI)=UD(GI)+U0(NSTEP,JGL)*N(JLOC,GI)
              VD(GI)=VD(GI)+V0(NSTEP,JGL)*N(JLOC,GI)
              WD(GI)=WD(GI)+W0(NSTEP,JGL)*N(JLOC,GI)
              
           ENDDO
           CORIOLIS(GI) = 2.0*FUNOME(XD(GI),YD(GI),ZD(GI))
           ! testing
           !              CORIOLIS(GI) =0.0
        ENDDO
        
        DO 20 ILOC =1,NLOC
           IGL=NDGLNO((ELE-1)*NLOC+ILOC)
           
           !Calculate
           !------
           !PHI_jj^p*[U_x+V_y+W_z]
           !------------------
           !PHI_jj^u*[U_x+UD*U_x+VD*U_y+WD*U_z-FV+P_x]
           !+D(PHI_jj^u)/DX*DU/DX+D(PHI_jj^u)/DY*DU/DY+D(PHI_jj^u)/DZ*DU/DZ
           !-----------------
           !PHI_jj^v*[V_x+UD*V_x+VD*V_y+WD*V_z+FU+P_y]
           !+D(PHI_jj^v)/DX*DV/DX+D(PHI_jj^v)/DY*DV/DY+D(PHI_jj^v)/DZ*DV/DZ
           !-----------------
           !PHI_jj^w*[W_x+UD*W_x+VD*W_y+WD*W_z+P_z]
           !+D(PHI_jj^w)/DX*DW/DX+D(PHI_jj^w)/DY*DW/DY+D(PHI_jj^w)/DZ*DW/DZ
           !---------------------------------


           DO JLOC=1,NLOC
              JGL=NDGLNO((ELE-1)*NLOC+JLOC)
              IF(IGL.EQ.JGL) THEN
                 U_t =0.0
                 V_t =0.0
                 W_t =0.0
                 U_t1 =0.0
                 V_t1 =0.0
                 W_t1 =0.0 
                 K =0
                 DO II=findrm(igl),findrm(igl+1)-1
                    K=k+1
                    U_t = U_t+ U0(NSTEP,COLM(II))
                    V_t = V_t+ V0(NSTEP,COLM(II))
                    W_t = W_t+ W0(NSTEP,COLM(II))
                    U_t1 = U_t1+ U0(NSTEP-1,COLM(II))
                    V_t1 = V_t1+ V0(NSTEP-1,COLM(II))
                    W_t1 = W_t1+ W0(NSTEP-1,COLM(II))
                 ENDDO
                 U_t = U_t/k
                 V_t = V_t/k
                 W_t = W_t/k
                 U_t1 = U_t1/k
                 V_t1 = V_t1/k
                 W_t1 = W_t1/k
                 
              ELSE
                 U_t=U0(NSTEP,JGL)
                 V_t=V0(NSTEP,JGL)
                 W_t=W0(NSTEP,JGL)
                 U_t1=U0(NSTEP-1,JGL)
                 V_t1=V0(NSTEP-1,JGL)
                 W_t1=W0(NSTEP-1,JGL)
              ENDIF
              
              
              RN1a=0.0
              RN2a=0.0
              RN3a=0.0
              RN1b=0.0
              RN2b=0.0
              RN3b=0.0
              
              FU  =0.0
              FV  =0.0
              U_x  =0.0
              V_y  =0.0
              W_Z  =0.0
              UU_x =0.0
              VU_y =0.0
              WU_z =0.0
              UV_x =0.0
              VV_y =0.0
              WV_z =0.0
              UW_x =0.0
              VW_y =0.0
              WW_z =0.0
              U_xx =0.0
              U_yy =0.0
              U_zz =0.0
              V_xx =0.0
              V_yy =0.0
              V_zz =0.0
              W_xx =0.0
              W_yy =0.0
              W_zz =0.0
              P_x  =0.0
              P_y  =0.0
              P_z  =0.0
              PG_x  =0.0
              PG_y  =0.0
              PG_z  =0.0
              soux=0.0
              souy=0.0
              souz=0.0
              
              DO GI=1,NGI
                 
                 RN1a=RN1a+ N(ILOC,GI)*DETWEI(GI)*U_t1*N(JLOC,GI)
                 RN2a=RN2a+ N(ILOC,GI)*DETWEI(GI)*V_t1*N(JLOC,GI)
                 RN3a=RN3a+ N(ILOC,GI)*DETWEI(GI)*W_t1*N(JLOC,GI)
                 
                 RN1b=RN1b+ N(ILOC,GI)*DETWEI(GI)*U_t*N(JLOC,GI)
                 RN2b=RN2b+ N(ILOC,GI)*DETWEI(GI)*V_t*N(JLOC,GI)
                 RN3b=RN3b+ N(ILOC,GI)*DETWEI(GI)*W_t*N(JLOC,GI)

                 FU  =FU+ N(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*U_t*N(JLOC,GI)
                 FV  =FV -N(ILOC,GI)*CORIOLIS(GI)*DETWEI(GI)*V_t*N(JLOC,GI)
                 U_x  =U_x+ N(ILOC,GI)*DETWEI(GI)*U_t*NX(JLOC,GI)
                 V_y  =V_y+ N(ILOC,GI)*DETWEI(GI)*V_t*NY(JLOC,GI)
                 W_Z  =W_z+ N(ILOC,GI)*DETWEI(GI)*W_t*NZ(JLOC,GI)
                 UU_x =UU_x+ N(ILOC,GI)*UD (GI)*DETWEI(GI)*U_t*NX(JLOC,GI)
                 VU_y =VU_y+ N(ILOC,GI)*VD (GI)*DETWEI(GI)*U_t*NY(JLOC,GI)
                 WU_z =WU_z+ N(ILOC,GI)*WD (GI)*DETWEI(GI)*U_t*NZ(JLOC,GI)
                 UV_x =UV_x+ N(ILOC,GI)*UD (GI)*DETWEI(GI)*V_t*NX(JLOC,GI)
                 VV_y =VV_y+ N(ILOC,GI)*VD (GI)*DETWEI(GI)*V_t*NY(JLOC,GI)
                 WV_z =WV_z+ N(ILOC,GI)*WD (GI)*DETWEI(GI)*V_t*NZ(JLOC,GI)
                 UW_x =UW_x+ N(ILOC,GI)*UD (GI)*DETWEI(GI)*W_t*NX(JLOC,GI)
                 VW_y =VW_y+ N(ILOC,GI)*VD (GI)*DETWEI(GI)*W_t*NY(JLOC,GI)
                 WW_z =WW_z+ N(ILOC,GI)*WD (GI)*DETWEI(GI)*W_t*NZ(JLOC,GI)

                 U_xx =U_xx+ NX(ILOC,GI)*MUPTXX(JGL)*DETWEI(GI)*U_t*NX(JLOC,GI)
                 U_yy =U_yy+ NY(ILOC,GI)*MUPTYY(JGL)*DETWEI(GI)*U_t*NY(JLOC,GI)
                 U_zz =U_zz+ NZ(ILOC,GI)*MUPTZZ(JGL)*DETWEI(GI)*U_t*NZ(JLOC,GI)
                 V_xx =V_xx+ NX(ILOC,GI)*MUPTXX(JGL)*DETWEI(GI)*V_t*NX(JLOC,GI)
                 V_yy =V_yy+ NY(ILOC,GI)*MUPTYY(JGL)*DETWEI(GI)*V_t*NY(JLOC,GI)
                 V_zz =V_zz+ NZ(ILOC,GI)*MUPTZZ(JGL)*DETWEI(GI)*V_t*NZ(JLOC,GI)
                 W_xx =W_xx+ NX(ILOC,GI)*MUPTXX(JGL)*DETWEI(GI)*W_t*NX(JLOC,GI)
                 W_yy =W_yy+ NY(ILOC,GI)*MUPTYY(JGL)*DETWEI(GI)*W_t*NY(JLOC,GI)
                 W_zz =W_zz+ NZ(ILOC,GI)*MUPTZZ(JGL)*DETWEI(GI)*W_t*NZ(JLOC,GI)
!                 U_zz =0.0
!                 V_zz =0.0
!                 W_zz =0.0

                 if(IFWIND) then
                    ! WIND STRESS
                    !-------------
                    !           CALL WINDYSOURCE(SOURCX,SOURCY,X,Y,NONODS,1)
                    soux = soux -DT*N(ILOC,GI)*DETWEI(GI)*SOURCX(IGL)*N(JLOC,GI)
                 endif
                 
                 
              ENDDO
              
              CALL POSINMAT(COUNT,IGL,JGL,NONODS,FINDRM,COLM,NCOLM)
              
              BIGM(COUNT)=BIGM(COUNT)+RN1b-RN1a+DT*(UU_x+VU_y+WU_z+U_xx+U_yy+U_zz+FV)-DT*soux
              BIGM(NCOLM+COUNT)=BIGM(NCOLM+COUNT)+RN2b-RN2a+DT*(UV_X+VV_y+WV_z+V_xx+V_yy+V_zz+FU)-DT*souy
              BIGM(2*NCOLM+COUNT)=BIGM(2*NCOLM+COUNT)+RN3b-RN3a+DT*(UW_x+VW_y+WW_z+W_xx+W_yy+W_zz)-DT*souz
              
           ENDDO
           
           DO JLOC=1,MLOC
              JGLP=PNDGLN((ELE-1)*MLOC+JLOC)
              P_x=0.0
              P_y=0.0
              P_z=0.0
              DO GI=1,NGI
                 P_x=P_x + N(ILOC,GI)*DETWEI(GI)*PHI0(NSTEP,JGLP)*MX(JLOC,GI)
                 P_y=P_y + N(ILOC,GI)*DETWEI(GI)*PHI0(NSTEP,JGLP)*MY(JLOC,GI)
                 P_z=P_z + N(ILOC,GI)*DETWEI(GI)*PHI0(NSTEP,JGLP)*MZ(JLOC,GI)
              ENDDO
              CALL POSINMAT(COUNT,JGLP,IGL,FREDOP,FINDCT,COLCT,NCT)
              BIGM_P(COUNT)=BIGM_P(COUNT)+DT*P_x
              BIGM_P(NCT+COUNT)=BIGM_P(NCT+COUNT)+DT*P_y
              BIGM_P(2*NCT+COUNT)=BIGM_P(2*NCT+COUNT)+DT*P_z
           ENDDO
           
           DO JLOC=1,MPGLOC
              JGLPG=PGNDGLNO((ELE-1)*MPGLOC+JLOC)
              PG_x=0.0
              PG_y=0.0
              PG_z=0.0
              if(JGLPG.GT.PGNODS.OR.JGLPG.LE.0) THEN
                 write(3,*) PGNDGLNO((ELE-1)*MPGLOC+JLOC),ELE,MPGLOC,JLOC,JGLPG
                 stop 344
              endif
              DO GI=1,NGI
                 PG_x = PG_x + N(ILOC,GI)*DETWEI(GI)*PHI0_PG(NSTEP,JGLPG)*MPGX(JLOC,GI)
                 PG_y = PG_y + N(ILOC,GI)*DETWEI(GI)*PHI0_PG(NSTEP,JGLPG)*MPGY(JLOC,GI)
                 PG_z = PG_z + N(ILOC,GI)*DETWEI(GI)*PHI0_PG(NSTEP,JGLPG)*MPGZ(JLOC,GI)
              ENDDO
              CALL POSINMAT(COUNT,JGLPG,IGL,PGNODS,FINDCT_PG,COLCT_PG,NCT_PG)
              BIGM_PG(COUNT)=BIGM_PG(COUNT)+DT*PG_x
              BIGM_PG(NCT_PG+COUNT)=BIGM_PG(NCT_PG+COUNT)+DT*PG_y
              BIGM_PG(2*NCT_PG+COUNT)=BIGM_PG(2*NCT_PG+COUNT)+DT*PG_z
           ENDDO
           
20         CONTINUE
           !              write(3,*) '1',BIGM(COUNT)
           
           
30         CONTINUE
           
           !              BIGM=ABS(BIGM)
           !              BIGM_P=ABS(BIGM_P)
           !              BIGM_PG=ABS(BIGM_PG)
           
           VEC=0.0
           DO II=1,NONODS
              LOWER=FINDRM(II)
              UPPER=FINDRM(II+1)-1
              DO KK=LOWER,UPPER
                 VEC(II)=VEC(II)+BIGM(KK)
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
          VEC_TOTAL(II)=VEC_TOTAL(II)+ABS(VEC(II))+ABS(VEC(NONODS+II))+ABS(VEC(2*NONODS+II))
          !                write(3,*) II,VEC_TOTAL(II),VEC(II),VEC(NONODS+II),VEC(2*NONODS+II)
       ENDDO
       
       HESSIAN_U=0.0
       HESSIAN_V=0.0
       HESSIAN_W=0.0
       allocate(uu(nonods))
       uu(1:nonods)=U0_ADJ(NSTEP,1:nonods)
       CALL CalHessian(NONODS,XNONOD,NLOC,NGI,TOTELE,3,UU,HESSIAN_U,         &
            X,Y,Z,XONDGL,NDGLNO,N,NLX,NLY,NLZ,DCYL,D3,WEIGHT)
       uu(1:nonods)=V0_ADJ(NSTEP,1:nonods)
       CALL CalHessian(NONODS,XNONOD,NLOC,NGI,TOTELE,3,UU,HESSIAN_V,         &
            X,Y,Z,XONDGL,NDGLNO,N,NLX,NLY,NLZ,DCYL,D3,WEIGHT)
       uu(1:nonods)=W0_ADJ(NSTEP,1:nonods)
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
       
100    CONTINUE
       
       DO II=1,NONODS
          DO JJ=1,9
             !                    FERROR((II-1)*9+JJ)=FERROR((II-1)*9+JJ)/ABS(VEC_TOTAL(II))
             FERROR((II-1)*9+JJ)=FERROR((II-1)*9+JJ)/DF
          ENDDO
       ENDDO
       
1100   continue
       !-------------------------------!
       ! DEALLOCATE THE LOCAL MATRICES !
       !-------------------------------!
       write(3,*) '22'
       DEALLOCATE(BIGM)
       DEALLOCATE(BIGM_P)
       DEALLOCATE(BIGM_PG)
       DEALLOCATE(VEC)
       DEALLOCATE(DETWEI)
       DEALLOCATE(UD)
       DEALLOCATE(VD)
       DEALLOCATE(WD)
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
       DEALLOCATE(SOURCX)
       DEALLOCATE(SOURCY)
       
       ewrite(3,*) 'in the middle of allocationg'
       ! for geographic pressure
       DEALLOCATE(MPGX)
       DEALLOCATE(MPGY)
       DEALLOCATE(MPGZ)
       
       DEALLOCATE(PHI0_PG)
       ewrite(3,*) 'exit REDUCED_MODEL_DIFFCOEF'
       RETURN
     END SUBROUTINE REDUCED_MODEL_DIFFCOEF1_PG1_errormeasure
     
end module reduced_model_error_measure
