
PROGRAM MAIN 
! **********************************************************
! THIS PERFORMS A 1D CV-ADVECTION TEST CASE ****************
! **********************************************************
  use spact
  use multiphase_1D_engine
  use multiphase_EOS
  use solvers_module

  IMPLICIT NONE

  INTEGER, PARAMETER :: NPHASE = 1, TOTELE = 20,  NDIM = 1, &
       U_NLOC = 3, XU_NLOC = 3, CV_NLOC = 3, X_NLOC = 3, P_NLOC = 3, &
       NCOEF = 10, &
       CV_SNLOC = 1, U_SNLOC = 1,P_SNLOC = 1, X_SNLOC = 1, STOTEL = 2,  &
       MAT_NLOC = 3, &
       U_ELE_TYPE = 1, CV_ELE_TYPE = 1, CV_SELE_TYPE = 1,  &
       P_ELE_TYPE = 1, MAT_ELE_TYPE = 1, U_SELE_TYPE = 1, &
       NITS = 1, T_DISOPT = 8, T_DG_VEL_INT_OPT = 1, NUABS_COEFS = 1

  INTEGER, PARAMETER :: CV_NONODS = (CV_NLOC - 1) * TOTELE + 1,  &
       U_NONODS = U_NLOC * TOTELE, &
       P_NONODS = P_NLOC * TOTELE, &
       CV_PHA_NONODS = CV_NONODS * NPHASE, &
       U_PHA_NONODS = U_NONODS * NPHASE,   &
       NLENMCY = U_NONODS * NPHASE + CV_NONODS, & 
       MAT_NONODS = MAT_NLOC * TOTELE, &
       X_NONODS = CV_NONODS, &
       NCP_COEFS = NPHASE

  REAL     :: PATMOS, DT, ACCTIM ! Reference pressure, time-step size and accumulated time
  INTEGER  :: ITIME, NTIME
  
  REAL, PARAMETER :: V_BETA=0.0, T_BETA=0.0, T_THETA=-1.0

  ! Working vectors
  REAL, DIMENSION( U_PHA_NONODS )  :: U, V, W, UOLD, VOLD, WOLD
  REAL, DIMENSION( CV_PHA_NONODS ) :: DEN, DENOLD, DERIV, VOLFRA, VOLOLD, T, TOLD, CPDEN, CPDENOLD
  REAL, DIMENSION( CV_PHA_NONODS ) :: CV_ONE
  REAL, DIMENSION( CV_NONODS )     :: P, POLD
  REAL, DIMENSION( MAT_NONODS,NDIM,NDIM,NPHASE )  :: TDIFFUSION

  ! Pointers
  INTEGER, DIMENSION( TOTELE * U_NLOC )    :: U_NDGLN
  INTEGER, DIMENSION( TOTELE * XU_NLOC )   :: XU_NDGLN
  INTEGER, DIMENSION( TOTELE * CV_NLOC )   :: P_NDGLN, CV_NDGLN, X_NDGLN
  INTEGER, DIMENSION( TOTELE * MAT_NLOC )  :: MAT_NDGLN

  ! Surface mapping
  INTEGER, DIMENSION( STOTEL * CV_SNLOC ) :: CV_SNDGLN
  INTEGER, DIMENSION( STOTEL * U_SNLOC )  :: U_SNDGLN, P_SNDGLN, X_SNDGLN
  REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ) :: SUF_ONE_BC

  ! For EOS:
  REAL, DIMENSION( NPHASE, NCOEF )    :: EOS_COEFS
  INTEGER, DIMENSION( NPHASE )        :: EOS_OPTION, CP_OPTION
  REAL, DIMENSION( NPHASE, NCP_COEFS ):: CP_COEFS

  ! For velocities and positions:
  REAL, DIMENSION( X_NONODS ) :: X, Y, Z
  REAL, DIMENSION( U_NONODS*NPHASE ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG
  REAL, DIMENSION( CV_PHA_NONODS ) :: SOURCV
  REAL, DIMENSION( CV_PHA_NONODS,NPHASE,NPHASE ) :: ABSORBV

  ! For matrix sparcity
  INTEGER :: NELE_PHA, NPHA_TOTELE, NCOLDGM_PHA

  INTEGER, DIMENSION( NLENMCY + 1 )       :: FINMCY ! Force balance plus cty multi-phase eqns
  INTEGER, DIMENSION( NLENMCY )           :: MIDMCY

  INTEGER, DIMENSION( CV_PHA_NONODS )     :: MIDACV ! CV multi-phase eqns 
  INTEGER, DIMENSION( CV_PHA_NONODS + 1 ) :: FINACV

  INTEGER, DIMENSION( TOTELE + 1)         :: FINELE ! Element connectivity
  INTEGER, DIMENSION( TOTELE )            :: MIDELE

  INTEGER, DIMENSION( U_NONODS )          :: CENTC ! C sparcity operating on pressure in force balance
  INTEGER, DIMENSION( U_NONODS + 1 )      :: FINDC

  INTEGER, DIMENSION( CV_NONODS )         :: MIDCMC ! pressure matrix for projection method
  INTEGER, DIMENSION( CV_NONODS + 1 )     :: FINDCMC

  INTEGER, DIMENSION( CV_NONODS )         :: CENTCT ! CT sparcity - global cty eqn.
  INTEGER, DIMENSION( CV_NONODS + 1 )     :: FINDCT

  INTEGER, DIMENSION( U_PHA_NONODS + 1 )  :: FINDGM_PHA ! Force balance sparcity
  INTEGER, DIMENSION( U_PHA_NONODS )      :: MIDDGM_PHA

  INTEGER, DIMENSION( CV_NONODS + 1 )     :: FINDM ! Sparcity for the CV-FEM
  INTEGER, DIMENSION( CV_NONODS )         :: MIDM
  
    INTEGER, DIMENSION( STOTEL * NPHASE ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
    
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ) :: SUF_T_BC, SUF_D_BC
    REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2


  INTEGER, DIMENSION( : ), allocatable    :: COLMCY, COLACV, ACV, COLELE, COLCT, COLC, &
                                   COLDGM_PHA, COLCMC, COLM
  REAL, DIMENSION( : ), allocatable       :: DGM_PHA

  ! Absorption term
  INTEGER, DIMENSION( NPHASE ) :: UABS_OPTION
  REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ) :: SIGMA 
  REAL, DIMENSION( NPHASE, NUABS_COEFS ) :: UABS_COEFS
  REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ) :: U_ABSORB
  REAL, DIMENSION( U_NONODS * NPHASE ) :: U_SOURCE

  ! For assembling and solving
  REAL, DIMENSION( CV_NONODS ) :: RHS, RHS_CV, DIAG_PRES
  REAL, DIMENSION( TOTELE ) :: VOLFRA_PORE
  REAL, DIMENSION( TOTELE, NDIM, NDIM ) :: PERM
  REAL, DIMENSION( CV_PHA_NONODS ) :: T_SOURCE, V_SOURCE
  REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ) :: T_ABSORB, V_ABSORB
  LOGICAL :: GETCV_DISC, GETCT

  ! Solver and Assembling options (NITS: no of non-linear iterations)
  INTEGER :: U_DISOPT, V_DISOPT, W_DISOPT, U_DG_VEL_INT_OPT, &
       V_DG_VEL_INT_OPT, W_DG_VEL_INT_OPT
  REAL ::  U_THETA, V_THETA


  ! For the BC's:
  INTEGER                                 :: NOBCT, NOBCV, NOBCU
  INTEGER, DIMENSION( : ), allocatable    :: BCT1, BCV1, BCU1, BCT2, BCV2, BCU2, dummy

  INTEGER :: NOD, II, COL, COUNT, NACV, NCOLCT, NCOLACV, NCT, NCOLC, NCMC, NCOLELE, NCOLMCY, NELE, &
       ITS, PNOD,   &
       BOT_BC_TYPE, TOP_BC_TYPE, ELE, ILOC, JLOC, CV_NOD, U_NOD, XU_NOD, X_NOD, &
       NCOLCMC, NCOLM, NC, MAT_NOD
  INTEGER :: MXNELE, MX_NCT, MX_NC, MX_NCOLCMC, MX_NCOLDGM_PHA,  & 
       MX_NCOLMCY, MX_NCOLACV, MX_NCOLM, MX_NACV

  REAL    :: DOMAIN_LENGTH, DX, UBOT, UTOP, DEN_IN_TOP, DEN_IN_BOT,            &
       P_INI, PERM_VALUE, XACC

  INTEGER, PARAMETER :: IZERO = 0
  REAL, PARAMETER :: RZERO = 0.


  ! Initialise to zero some variables
  CV_SNDGLN( 1 ) = 1
  CV_SNDGLN( 2 ) = CV_NONODS
  U_SNDGLN( 1 ) = 1
  U_SNDGLN( 2 ) = U_NONODS
  V_SOURCE  =0.0
  V_ABSORB =0.0
  
  WIC_T_BC=0
  WIC_D_BC=0 
  WIC_U_BC=0
  WIC_T_BC(1)=1
  SUF_T_BC=0.0 
  SUF_D_BC=0.0
  SUF_U_BC=0.0 
  SUF_V_BC=0.0 
  SUF_W_BC=0.0
  SUF_T_BC_ROB1=0.0
  SUF_T_BC_ROB2=0.0
  SUF_T_BC(1)=1.0 
  
!  NTIME = 10
  NTIME = 200

!  NTIME = 1
  
  DOMAIN_LENGTH = 1.0
  DX = DOMAIN_LENGTH / REAL(TOTELE)
  DT = 0.001

  ! U_NDGLN, XU_NDGLN, P_NDGLN, CV_NDGLN:
  U_NOD  = 0
  XU_NOD  = 0
  CV_NOD = 0
  X_NOD  = 0
  MAT_NOD = 0.0
  Loop_Elements: DO ELE = 1, TOTELE

     DO ILOC = 1, U_NLOC ! storing velocity nodes
        U_NOD = U_NOD + 1
        U_NDGLN( ( ELE - 1 ) * U_NLOC + ILOC ) = U_NOD
     END DO

     DO ILOC = 1, XU_NLOC ! storing velocity nodes
        XU_NOD = XU_NOD + 1
        XU_NDGLN( ( ELE - 1 ) * XU_NLOC + ILOC ) = XU_NOD
     END DO
     XU_NOD = XU_NOD - 1

     DO ILOC = 1, CV_NLOC ! storing scalar nodes
        CV_NOD = CV_NOD + 1
        P_NDGLN( ( ELE - 1 ) * CV_NLOC + ILOC ) = CV_NOD
        CV_NDGLN( ( ELE - 1 ) * CV_NLOC + ILOC ) = CV_NOD
     END DO

     CV_NOD = CV_NOD - 1

     DO ILOC = 1, X_NLOC ! storing scalar nodes
        X_NOD = X_NOD + 1
        X_NDGLN( ( ELE - 1 ) * X_NLOC + ILOC ) = X_NOD
	X(X_NOD)=REAL(ELE-1)*DX + REAL(ILOC-1)*DX/REAL(X_NLOC-1)
     END DO

     X_NOD = X_NOD - 1

     DO ILOC = 1, MAT_NLOC ! storing scalar nodes
        MAT_NOD = MAT_NOD + 1
        MAT_NDGLN( ( ELE - 1 ) * MAT_NLOC + ILOC ) = MAT_NOD
     END DO

  END DO Loop_Elements
  
! Initialise T and TOLD...
  T=0.0
  Loop_Elements2: DO ELE = 1, TOTELE
     DO ILOC = 1, CV_NLOC 
        X_NOD  = X_NDGLN( ( ELE - 1 ) * X_NLOC + ILOC )
        CV_NOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC+ ILOC )
	IF((X(X_NOD).GT.0.199).AND.(X(X_NOD).LT.0.401)) T(CV_NOD)=1.0
     END DO
  END DO Loop_Elements2
  TOLD=T


  ! Initialising working arrays:
  NU=1.0 
  NV=0.0 
  NW=0.0 
  NUOLD=1.0 
  NVOLD=0.0 
  NWOLD=0.0 
  UG=0.0 
  VG=0.0 
  WG=0.0
  Y=0.0
  Z=0.0
  CV_ONE=1.0
  
  ! Allocate boundary conditions...  
  ! Left boundary for phase 1: 
  SUF_ONE_BC( 1 ) = 1.0
  ! Right boundary for phase 1: 
  SUF_ONE_BC( 2 ) = 1.0

  VOLFRA_PORE=1.0
  PERM=0.0
!  TDIFFUSION=0.1
!  TDIFFUSION=0.01
  TDIFFUSION=0.0

  CP_COEFS = 1.0
  EOS_OPTION = 1

  ! Define matrix sparcity...          
  MXNELE = ( 2 * NDIM + 1 ) * TOTELE
  MX_NCT = CV_NONODS * ( 2 * U_NLOC + 1 ) * NDIM * NPHASE
  MX_NC  = MX_NCT
  MX_NCOLCMC = ( 2 * CV_NLOC + 1 ) * CV_NONODS
  MX_NCOLDGM_PHA = MXNELE * ( U_NLOC * NDIM )**2 * NPHASE + TOTELE * ( U_NLOC * NDIM * NPHASE )**2
  MX_NCOLMCY= MX_NCOLDGM_PHA + MX_NCT + MX_NC + CV_NONODS
  MX_NCOLACV = ( 2 * NDIM + 1 ) * CV_NONODS * NPHASE + CV_NONODS * ( NPHASE - 1 ) * NPHASE
  MX_NCOLM = MXNELE * U_NLOC * U_NLOC

  ALLOCATE( COLMCY( MX_NCOLMCY ))
  ALLOCATE( COLACV( MX_NCOLACV )) 
  ALLOCATE( COLELE( MXNELE ))
  ALLOCATE( COLCT(  MX_NCT ))
  ALLOCATE( COLC(   MX_NC  ))
  ALLOCATE( COLDGM_PHA( MX_NCOLDGM_PHA ))
  ALLOCATE( COLCMC( MX_NCOLCMC ))
  ALLOCATE( COLM( MX_NCOLM ))

  CALL GET_SPARS_PATS(   &
       NDIM, U_PHA_NONODS, CV_PHA_NONODS,   &
       U_NONODS, CV_NONODS, &
       U_NLOC, CV_NLOC, NPHASE, TOTELE, U_NDGLN,  &
       MX_NCOLACV, NCOLACV, FINACV, COLACV, MIDACV,  & ! CV multi-phase eqns (e.g. vol frac, temp)
       NLENMCY, MX_NCOLMCY, NCOLMCY, FINMCY, COLMCY, MIDMCY, & ! Force balance plus cty multi-phase eqns
       MXNELE, NCOLELE, MIDELE, FINELE, COLELE, & ! Element connectivity
       MX_NCOLDGM_PHA, NCOLDGM_PHA, COLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, & ! Force balance sparcity
       MX_NCT, NCOLCT, FINDCT, COLCT, &! CT sparcity - global cty eqn
       MX_NC, NCOLC, FINDC, COLC,  & ! C sparcity operating on pressure in force balance
       MX_NCOLCMC, NCOLCMC, FINDCMC, COLCMC, MIDCMC, & ! pressure matrix for projection method
       MX_NCOLM, NCOLM, FINDM, COLM, MIDM ) ! CV-FEM matrix

  ALLOCATE( DGM_PHA( NCOLDGM_PHA ))

  if( .false. ) then
     open( 15, file = 'CheckSparcityMatrix.dat', status = 'unknown' )
     write( 15, * )'########## FINMCY, MIDMCY, COLMCY ##################'
     write(15, * )'NCOLMCY:', NCOLMCY
     call checksparcity( .true., 15, NCOLMCY, NLENMCY, MX_NCOLMCY, FINMCY, MIDMCY, COLMCY )

     write( 15, * )'########## FINACV, COLACV, MIDACV ##################'
     write(15, * )'NCOLACV:', NCOLACV
     call checksparcity( .true., 15, NCOLACV, CV_PHA_NONODS, MX_NCOLACV, FINACV, MIDACV, COLACV  )

     write( 15, * )'########## FINELE, MIDELE, COLELE  ##################'
     write(15, * )'NCOLELE:',NCOLELE 
     call checksparcity( .true., 15, NCOLELE, TOTELE, MXNELE, FINELE, MIDELE, COLELE )

     allocate( dummy( CV_NONODS ))
     write( 15, * )'########## FINDCT, COLCT ##################'
     write(15, * )'NCOLCT:', NCOLCT
     call checksparcity( .false., 15, NCOLCT, CV_NONODS, MX_NCT, FINDCT, dummy, COLCT  )
     deallocate( dummy )

     allocate( dummy( U_NONODS ))
     write( 15, * )'########## FINDC, COLC ##################'
     write(15, * )'NCOLC:', NCOLC
     call checksparcity( .false., 15, NCOLC, U_NONODS, MX_NC, FINDC, dummy, COLC )
     deallocate( dummy )

     write( 15, * )'########## FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA ##################'
     write(15, * )'NCOLDGM_PHA:',NCOLDGM_PHA 
     call checksparcity( .true., 15, NCOLDGM_PHA, U_PHA_NONODS, MX_NCOLDGM_PHA, FINDGM_PHA, MIDDGM_PHA, COLDGM_PHA )

     write( 15, * )'########## FINDCMC, MIDCMC, COLCMC ##################'
     write(15, * )'NCOLCMC:',NCOLCMC 
     call checksparcity( .true., 15, NCOLCMC, CV_NONODS, MX_NCOLCMC, FINDCMC, MIDCMC, COLCMC )

     write( 15, * )'########## FINDM, MIDM, COLM ##################'
     write(15, * )'NCOLM:',NCOLM 
     call checksparcity( .true., 15, NCOLM, CV_NONODS, MX_NCOLM, FINDM, MIDM, COLM )
  end if

  Loop_Time: DO ITIME=1,NTIME

     ACCTIM = ACCTIM + DT
     TOLD=T

     ! Non linear its:
     Loop_ITS: DO ITS = 1, NITS
        print *,'itime,its:',itime,its

        CALL INTENERGE_ASSEM_SOLVE(  &
             NCOLACV, FINACV, COLACV, MIDACV, & ! CV sparcity pattern matrix
             NCOLCT, FINDCT, COLCT, &
             CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
             U_ELE_TYPE, CV_ELE_TYPE, CV_SELE_TYPE,  &
             NPHASE,  &
             CV_NLOC, U_NLOC, X_NLOC,  &
             CV_NDGLN, X_NDGLN, U_NDGLN, &
             CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
             X, Y, Z, &
             NU, NV, NW, NUOLD, NVOLD, NWOLD, UG, VG, WG, T, TOLD, CV_ONE, CV_ONE, &
             MAT_NLOC,MAT_NDGLN,MAT_NONODS,TDIFFUSION, &
             T_DISOPT, T_DG_VEL_INT_OPT, DT, T_THETA, T_BETA, &
             SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
             SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
             WIC_T_BC, WIC_D_BC, WIC_U_BC, &
             DERIV, P, &
             T_SOURCE, T_ABSORB, &
             NDIM,  &
             NCOLM, FINDM, COLM, MIDM, &
             XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE )
      END DO Loop_ITS
	 
  END DO Loop_Time


  ! Output the results...
  PRINT *, 'T:'
  XACC = 0.0
  DO ELE = 1, TOTELE
     DO ILOC = 1, CV_NLOC
        CV_NOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + ILOC )
        PRINT *, XACC, T( CV_NOD )
        XACC = XACC + DX / REAL( CV_NLOC - 1 )
     END DO
     XACC = XACC - DX / REAL( CV_NLOC - 1 )
  END DO

END PROGRAM MAIN


