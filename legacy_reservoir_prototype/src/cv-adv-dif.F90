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
!    amcgsoftware@imperial.ac.uk
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

module cv_advection

  use fldebug

contains

  SUBROUTINE CV_ASSEMB( CV_RHS, &
       NCOLACV, ACV, FINACV, COLACV, MIDACV, &
       NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
       CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
       CV_ELE_TYPE,  &
       NPHASE,  &
       CV_NLOC, U_NLOC, X_NLOC, &
       CV_NDGLN, X_NDGLN, U_NDGLN, &
       CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
       X, Y, Z, &
       NU, NV, NW, NUOLD, NVOLD, NWOLD, &
       T, TOLD, DEN, DENOLD, &
       MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
       CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, &
       SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, &
       SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
       WIC_T_BC, WIC_D_BC, WIC_U_BC, &
       DERIV, CV_P, &
       SOURCT, ABSORBT, VOLFRA_PORE, & 
       NDIM, GETCV_DISC, GETCT, &
       NCOLM, FINDM, COLM, MIDM, &
       XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
       OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, T_FEMT, DEN_FEMT, &
       IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
       THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
       SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
       NOIT_DIM, &
       MASS_ERROR_RELAX2_NOIT, MEAN_PORE_CV )
    !  =====================================================================
    !     In this subroutine the advection terms in the advection-diffusion
    !     equation (in the matrix and RHS) are calculated as ACV and CV_RHS. 
    !
    !     This routine uses a Control Volume (CV) formulation to compute
    !     the advection terms. The general procedure is as follows:
    !
    !        1. For each node-pair, define which node is the donor, which is
    !           the receptor and define an "upwind" value of the field being 
    !           advected (and the accompanying "density"; see note below)
    !        2. Calculate the volume flux across the CV face that separates
    !           these two nodes
    !        3. Estimate the value of the advected variable at the control
    !           volume face.
    !        4. Using information from the donor, receptor and upwind nodes,
    !           limit the field face value (removes oscillations from the
    !           solution)
    !        5. Assemble the fluxes to form the matrix and rhs of the 
    !           advection equation
    !
    !     This procedure is implemented by considering the CV to be made up
    !     of a number of sub-control-volumes, which represent the part of
    !     the control volume within a given element.  The assembly of terms
    !     considers each of these sub-CVs in turn, calculating (and limiting)
    !     the flux across sub-CV faces that are external to the CV...
    !
    !     NOTE: Add in note about what density is in this sub!!!
    !
    !     To define the "upwind" value of the field variable, which is
    !     necessary for the limiting scheme, either:
    !
    !        A. The upwind value of the field variable to be advected is
    !           found by interpolation and stored in a matrix (TUPWIND)
    !        B. The neighbouring nodes are searched for the local maximum 
    !           and minimum
    !     
    !     The subroutine has several options...
    !
    !     Discretisation option
    !     ---------------------
    !      - The estimate of the face value may be determined in one of
    !        several ways.
    !      - The face value may be centered in time by either a specified 
    !        CV_THETA value, or a non-linear CV_THETA value that is determined 
    !        automatically.  
    !      - The face value may be limited using a univeral-limiter-type
    !        scheme, or a limited-downwind scheme that is ideal for INTERFACE
    !        TRACKING.  Alternatively no limiting can be applied.  
    !     
    !     These options are defined by the value of CV_DISOPT, which corresponds 
    !     to the clast digit of the GEM option NDISOT for the field in question.
    !
    !     CV_DISOPT=discretisation option in space and time
    !     ------------------------------------------------------------------
    !     CV_DISOPT   Method for face-value est.    Time-stepping     Limiting   
    !     ------------------------------------------------------------------
    !       =0      1st order in space          Theta=specified    UNIVERSAL
    !       =1      1st order in space          Theta=non-linear   UNIVERSAL
    !       =2      Trapezoidal rule in space   Theta=specified    UNIVERSAL
    !       =2 if isotropic limiter then FEM-quadratic & stratification adjust. Theta=non-linear 
    !       =3      Trapezoidal rule in space   Theta=non-linear   UNIVERSAL
    !       =4      Finite elements in space    Theta=specified    UNIVERSAL
    !       =5      Finite elements in space    Theta=non-linear   UNIVERSAL
    !       =6      Finite elements in space    Theta=specified    NONE
    !       =7      Finite elements in space    Theta=non-linear   NONE
    !       =8      Finite elements in space    Theta=specified    DOWNWIND+
    !       =9      Finite elements in space    Theta=non-linear   DOWNWIND+
    !
    !     CV_DG_VEL_INT_OPT=interface velocity calculation option between elements
    !
    !     Limiting scheme
    !     ---------------
    !     The limiting scheme is defined in the subroutine NVDFUNNEW; 
    !     the limited values are computed in subroutine ANONVDLIM/ONVDLIM.
    !     
    !     ONVDLIM is the original limiting algorithm
    !
    !     ANONVDLIM is a new anisoptropic limiting algorithm, which is 
    !     called if either ALOLIM=1 (where ALOLIM is an option flag set 
    !     in this subroutine), or if the interface tracking limiting option 
    !     is selected (CV_DISOPT=8/9).  ***In general ALOLIM appears to be set to 1 (GSC)
    !     
    !     NOTE: ANONVDLIM only works for TETS; for all other element types 
    !     ONVDLIM is used.
    !
    !
    !     IMPORTANT INPUTS:
    !     ----------------
    !     
    !     ACV   - Matrix for assembling the advection terms (empty on input)
    !     CV_RHS      - Right-hand side vector for advection-diffusion terms
    !     X,Y,Z    - Node co-ordinates
    !     NU       - Nodal velocity component
    !     T,TOLD   - New and old advected field values at nodes
    !     DEN,  - New and old "density" at nodes, which is actually a constant
    !     DENOLD     multiplying the advection diffusion equation for the field
    !     CV_DISOPT   - The discretisation/limiting option (see above)
    !     DT       - The time step
    !     CV_THETA    - The time-stepping discretisation parameter
    !     CV_BETA     - Conservative(1.)/non-conservative(0.) flag
    !     ELE_TYP   - Integer flag definining element type   
    !
    !     IMPORTANT OUTPUTS:
    !     -----------------
    !
    !     ACV   - Matrix updated to include the advection terms
    !     CV_RHS      - Right-hand side vector updated to include advection terms
    !
    !
    !     IMPORTANT LOCAL PARAMETERS:
    !     --------------------------
    !
    !     TIMOPT    - Temporal discretisation option, derived from CV_DISOPT.
    !                (1 for non-linear theta; 0 for theta specified (THETA))
    !
    !
    !***********************************************************************
    use shape_functions
    use matrix_operations
    use printout
    ! Inputs/Outputs
    IMPLICIT NONE
    INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
         TOTELE, &
         CV_ELE_TYPE, &
         NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_DISOPT, CV_DG_VEL_INT_OPT, NDIM, &
         NCOLM, XU_NLOC, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
         IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND

    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: T_FEMT, DEN_FEMT

    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN 
    INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
    INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
    INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
    INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN 
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
    INTEGER, DIMENSION( STOTEL * NPHASE * IGOT_T2 ), intent( in ) ::  WIC_T2_BC
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: CV_RHS
    REAL, DIMENSION( NCOLACV ), intent( inout ) :: ACV
    INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
    INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
    INTEGER, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: MIDACV 
    REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( inout ) :: CT
    ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
    REAL, DIMENSION( CV_NONODS ), intent( inout ) :: DIAG_SCALE_PRES 
    REAL, DIMENSION( CV_NONODS  ), intent( inout ) :: CT_RHS
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
    INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: T, TOLD, DEN, DENOLD
    REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( in ) :: T2, T2OLD
    REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: THETA_GDIFF
    REAL, DIMENSION( TOTELE*IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
         intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
    REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: TDIFFUSION
    REAL, intent( in ) :: DT, CV_THETA, CV_BETA
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC, SUF_D_BC
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ), intent( in ) :: SUF_T2_BC
    REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2 ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
    REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: DERIV
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: CV_P
    REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: SOURCT
    REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: ABSORBT
    REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE 
    LOGICAL, intent( in ) :: GETCV_DISC, GETCT, GET_THETA_FLUX, USE_THETA_FLUX
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
    INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
    INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
    INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
    INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, INTENT( IN ) :: NOIT_DIM
    REAL, DIMENSION( NOIT_DIM ), intent( in ) :: MASS_ERROR_RELAX2_NOIT
    REAL, DIMENSION( CV_NONODS ), intent( inout ) :: MEAN_PORE_CV

    ! Local variables - Allocatable Arrays
    INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, &
         CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
         JCOUNT_KLOC, JCOUNT_KLOC2, COLGPTS, CV_SLOC2LOC, U_SLOC2LOC
    INTEGER, DIMENSION( : , : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, &
         FACE_ELE, CV_NEILOC
    REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
         CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
         UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
         SCVFENLX, SCVFENLY, SCVFENLZ, &
         SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY,SUFENLZ
    REAL, DIMENSION( : , : ), allocatable :: SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
         SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
         SBUFENLX, SBUFENLY, SBUFENLZ
    REAL, DIMENSION( : ), allocatable :: CVWEIGHT,CVWEIGHT_SHORT,SCVFEWEIGH,SBCVFEWEIGH
    REAL, DIMENSION( : ), allocatable :: TMAX, TMIN, TOLDMAX, &
         TOLDMIN, DENMAX, DENMIN, DENOLDMAX, DENOLDMIN, CVNORMX, &
         CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
         FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
         SCVDETWEI, SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
         SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, T2MAX, T2MIN, T2OLDMAX, &
         T2OLDMIN
    INTEGER, DIMENSION( : ), allocatable :: TMAX_NOD, TMIN_NOD, TOLDMAX_NOD, &
         TOLDMIN_NOD, DENMAX_NOD, DENMIN_NOD, DENOLDMAX_NOD, DENOLDMIN_NOD, &
         T2MAX_NOD, T2MIN_NOD, T2OLDMAX_NOD, T2OLDMIN_NOD
    REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE,  &
         DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
    REAL, DIMENSION( : ), allocatable :: UP_WIND_NOD

    LOGICAL, DIMENSION( : ), allocatable :: X_SHARE,LOG_ON_BOUND
    LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE,U_ON_FACE
    LOGICAL, PARAMETER :: INCLUDE_PORE_VOL_IN_DERIV = .FALSE.
    ! Local variables - pointers, indeces, etc
    ! This is for decifering WIC_T_BC,WIC_D_BC,WIC_U_BC
    INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
         WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
         WIC_U_BC_DIRICHLET = 1
    !          ===> INTEGERS <====
    INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
         ELE, ELE2, GI, GCOUNT, SELE,   &
         NCOLGPTS, &
         CV_SILOC, U_KLOC, &
         CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
         CV_NODJ, CV_NODJ_IPHA, &
         CV_NODI, CV_NODI_IPHA, CV_NODI_JPHA, U_NODK, TIMOPT, &
         JCOUNT_IPHA, IMID_IPHA, &
         NFACE, X_NODI,  &
         CV_INOD, MAT_NODI, FACE_ITS, NFACE_ITS
    !        ===>  REALS  <===
    REAL :: NDOTQ, NDOTQOLD,  &
         INCOME, INCOMEOLD, HDC, FVT, FVTOLD, FVT2, FVT2OLD, &
         FVD, FVDOLD, LIMT, LIMTOLD, LIMT2, LIMT2OLD,&
         LIMD, LIMDOLD, FTHETA, VTHETA, &
         LIMDT, LIMDTOLD, LIMDTT2, LIMDTT2OLD, &
         FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
         TMID, TOLDMID, &
         DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, BCZERO, ROBIN1, ROBIN2, &
         SUM, &
         SUM_LIMT, SUM_LIMTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD
    ! Functions...
    !REAL :: R2NORM,FACE_THETA  
    !        ===>  LOGICALS  <===
    LOGICAL :: GETMAT, &
         D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
         NORMALISE, SUM2ONE, GET_GTHETA

    ewrite(3,*) 'In CV_ASSEMB'
    ewrite(3,*)'ENTERING CONSTRUCT_ADVECTION_DIFFUSION_CV()'
    ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA', &
         CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA

    CALL RETRIEVE_NGI( CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI,NFACE, &
         NDIM, CV_ELE_TYPE, CV_NLOC, U_NLOC) 

    GOT_DIFFUS = ( R2NORM( TDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0 )

    ! Allocate memory for the control volume surface shape functions, etc.
    ALLOCATE( JCOUNT_KLOC(  U_NLOC ))
    ALLOCATE( JCOUNT_KLOC2(  U_NLOC ))

    ALLOCATE( CVNORMX( SCVNGI ))
    ALLOCATE( CVNORMY( SCVNGI ))
    ALLOCATE( CVNORMZ( SCVNGI ))
    ALLOCATE( COLGPTS( CV_NLOC * SCVNGI )) !The size of this vector is over-estimated
    ALLOCATE( FINDGPTS( CV_NLOC + 1 ))
    ALLOCATE( SNDOTQ( SCVNGI ))
    ALLOCATE( SNDOTQOLD( SCVNGI ))
    ALLOCATE( CV_ON_FACE( CV_NLOC, SCVNGI ))
    ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))
    ALLOCATE( CV_OTHER_LOC( CV_NLOC ))
    ALLOCATE( U_OTHER_LOC( U_NLOC ))
    ALLOCATE( MAT_OTHER_LOC( MAT_NLOC ))
    ALLOCATE( X_SHARE( X_NONODS ))
    ALLOCATE( CVWEIGHT( CV_NGI ))
    ALLOCATE( CVN( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFEN( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENLX( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENLY( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENLZ( CV_NLOC, CV_NGI ))

    ALLOCATE( CVWEIGHT_SHORT( CV_NGI_SHORT ))
    ALLOCATE( CVN_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFEN_SHORT( CV_NLOC, CV_NGI_SHORT))
    ALLOCATE( CVFENLX_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFENLY_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFENLZ_SHORT( CV_NLOC, CV_NGI_SHORT ))

    ALLOCATE( UFEN( U_NLOC, CV_NGI)) 
    ALLOCATE( UFENLX( U_NLOC, CV_NGI ))
    ALLOCATE( UFENLY( U_NLOC, CV_NGI ))
    ALLOCATE( UFENLZ( U_NLOC, CV_NGI ))

    ALLOCATE( SCVFEN( CV_NLOC, SCVNGI ))
    ALLOCATE( SCVFENSLX( CV_NLOC, SCVNGI ))
    ALLOCATE( SCVFENSLY( CV_NLOC, SCVNGI ))
    ALLOCATE( SCVFENLX( CV_NLOC, SCVNGI ))
    ALLOCATE( SCVFENLY( CV_NLOC, SCVNGI ))
    ALLOCATE( SCVFENLZ( CV_NLOC, SCVNGI ))
    ALLOCATE( SCVFEWEIGH( SCVNGI ))

    ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
    ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
    ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
    ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
    ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
    ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

    ALLOCATE( SCVDETWEI( SCVNGI ))
    ALLOCATE( SRA( SCVNGI ))
    ALLOCATE( LOG_ON_BOUND(CV_NONODS))

    ALLOCATE( SBCVFEN( CV_SNLOC, SBCVNGI ))
    ALLOCATE( SBCVFENSLX( CV_SNLOC, SBCVNGI ))
    ALLOCATE( SBCVFENSLY( CV_SNLOC, SBCVNGI ))
    ALLOCATE( SBCVFEWEIGH( SBCVNGI ))
    ALLOCATE( SBCVFENLX( CV_SNLOC, SBCVNGI ))
    ALLOCATE( SBCVFENLY( CV_SNLOC, SBCVNGI ))
    ALLOCATE( SBCVFENLZ( CV_SNLOC, SBCVNGI ))
    ALLOCATE( SBUFEN( U_SNLOC, SBCVNGI ))
    ALLOCATE( SBUFENSLX( U_SNLOC, SBCVNGI ))
    ALLOCATE( SBUFENSLY( U_SNLOC, SBCVNGI ))
    ALLOCATE( SBUFENLX( U_SNLOC, SBCVNGI ))
    ALLOCATE( SBUFENLY( U_SNLOC, SBCVNGI ))
    ALLOCATE( SBUFENLZ( U_SNLOC, SBCVNGI ))

    ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
    ALLOCATE( U_SLOC2LOC( U_SNLOC )) 
    ALLOCATE( CV_SLOCLIST( NFACE, CV_SNLOC ))
    ALLOCATE( U_SLOCLIST( NFACE, U_SNLOC ))
    ALLOCATE( CV_NEILOC( CV_NLOC, SCVNGI ))

    ALLOCATE( SELE_OVERLAP_SCALE(CV_NLOC) )

    ALLOCATE( UGI_COEF_ELE(U_NLOC),  VGI_COEF_ELE(U_NLOC),  WGI_COEF_ELE(U_NLOC) )
    ALLOCATE( UGI_COEF_ELE2(U_NLOC), VGI_COEF_ELE2(U_NLOC), WGI_COEF_ELE2(U_NLOC) )
    ! The procity mapped to the CV nodes
    ALLOCATE( SUM_CV( CV_NONODS ))
    ALLOCATE( UP_WIND_NOD( CV_NONODS * NPHASE ))
    UP_WIND_NOD = 0.0

    ALLOCATE( ONE_PORE( TOTELE ))
    IF( INCLUDE_PORE_VOL_IN_DERIV ) THEN
       ! solve for actual velocity
       ONE_PORE = VOLFRA_PORE
    ELSE
       ! solve for vel=porcity*actual velocity
       ONE_PORE = 1.0
    ENDIF

    ewrite(3,*)'here1'

    D1 = ( NDIM == 1 )
    D3 = ( NDIM == 3 )
    DCYL= ( NDIM == -2 )

    GETMAT = .TRUE.

    X_SHARE = .FALSE.

    ! If using the original limiting scheme, the first step is to estimate 
    ! the upwind field value from the surrounding nodes

    ! Allocate memory for terms needed by GETGXYZ OR ONVDLIM
    ALLOCATE(      TMIN( CV_NONODS * NPHASE) )
    ALLOCATE(      TMAX( CV_NONODS * NPHASE) )
    ALLOCATE(   TOLDMIN( CV_NONODS * NPHASE) )
    ALLOCATE(   TOLDMAX( CV_NONODS * NPHASE) )
    ALLOCATE(      T2MIN( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(      T2MAX( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(   T2OLDMIN( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(   T2OLDMAX( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(    DENMIN( CV_NONODS * NPHASE) )
    ALLOCATE(    DENMAX( CV_NONODS * NPHASE) )
    ALLOCATE( DENOLDMIN( CV_NONODS * NPHASE) )
    ALLOCATE( DENOLDMAX( CV_NONODS * NPHASE) )

    ALLOCATE(      TMIN_NOD( CV_NONODS * NPHASE) )
    ALLOCATE(      TMAX_NOD( CV_NONODS * NPHASE) )
    ALLOCATE(   TOLDMIN_NOD( CV_NONODS * NPHASE) )
    ALLOCATE(   TOLDMAX_NOD( CV_NONODS * NPHASE) )
    ALLOCATE(      T2MIN_NOD( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(      T2MAX_NOD( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(   T2OLDMIN_NOD( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(   T2OLDMAX_NOD( CV_NONODS * NPHASE* IGOT_T2) )
    ALLOCATE(    DENMIN_NOD( CV_NONODS * NPHASE) )
    ALLOCATE(    DENMAX_NOD( CV_NONODS * NPHASE) )
    ALLOCATE( DENOLDMIN_NOD( CV_NONODS * NPHASE) )
    ALLOCATE( DENOLDMAX_NOD( CV_NONODS * NPHASE) )

    ewrite(3,*)'here1.2'
    NCOLGPTS = 0
    COLGPTS = 0
    FINDGPTS = 0
    ewrite(3,*)'in ASSEMB_FORCE_CTY',NCOLGPTS 
    ewrite(3,*)'in ASSEMB_FORCE_CTY, COLGPTS',size(COLGPTS), COLGPTS
    ewrite(3,*)'in ASSEMB_FORCE_CTY, FINDGPTS',size(FINDGPTS), FINDGPTS


    !     ======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========

    CALL CV_FEM_SHAPE_FUNS( &
                                ! Volume shape functions...
         NDIM, CV_ELE_TYPE,  & 
         CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
         CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
         CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
         UFEN, UFENLX, UFENLY, UFENLZ, &
                                ! Surface of each CV shape functions...
         SCVNGI, CV_NEILOC, CV_ON_FACE,  &  
         SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
         SCVFENLX, SCVFENLY, SCVFENLZ,  &
         SUFEN, SUFENSLX, SUFENSLY,  &
         SUFENLX, SUFENLY, SUFENLZ,  &
                                ! Surface element shape funcs...
         U_ON_FACE, NFACE, & 
         SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
         SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
         CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
         FINDGPTS, COLGPTS, NCOLGPTS, &
         SELE_OVERLAP_SCALE)  

    ! Determine FEMT (finite element wise) etc from T (control volume wise)
    ! Also determine the CV mass matrix MASS_CV and centre of the CV's XC_CV,YC_CV,ZC_CV. 
    ! This is for projecting to finite element basis functions... 
    ALLOCATE( FEMT( CV_NONODS * NPHASE ))
    ALLOCATE( FEMTOLD( CV_NONODS * NPHASE ))
    ALLOCATE( FEMDEN( CV_NONODS * NPHASE ))
    ALLOCATE( FEMDENOLD( CV_NONODS * NPHASE ))
    ALLOCATE( FEMT2( CV_NONODS * NPHASE * IGOT_T2 ))
    ALLOCATE( FEMT2OLD( CV_NONODS * NPHASE * IGOT_T2 ))
    ALLOCATE( MASS_CV( CV_NONODS ))
    ALLOCATE( MASS_ELE( TOTELE ))
    ALLOCATE( XC_CV( CV_NONODS ))
    ALLOCATE( YC_CV( CV_NONODS ))
    ALLOCATE( ZC_CV( CV_NONODS ))
    ALLOCATE( DTX_ELE( CV_NLOC, NPHASE,TOTELE ))
    ALLOCATE( DTY_ELE( CV_NLOC, NPHASE,TOTELE ))
    ALLOCATE( DTZ_ELE( CV_NLOC, NPHASE,TOTELE ))
    ALLOCATE( DTOLDX_ELE( CV_NLOC, NPHASE,TOTELE ))
    ALLOCATE( DTOLDY_ELE( CV_NLOC, NPHASE,TOTELE ))
    ALLOCATE( DTOLDZ_ELE( CV_NLOC, NPHASE,TOTELE ))

    ewrite(3,*)'here2'

    CALL PROJ_CV_TO_FEM_4( FEMT, FEMTOLD, FEMDEN, FEMDENOLD, T, TOLD, DEN, DENOLD, &
         IGOT_T2, T2,T2OLD, FEMT2,FEMT2OLD, &
         XC_CV, YC_CV, ZC_CV, MASS_CV, MASS_ELE,  &
         NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI_SHORT, CV_NLOC, CVN_SHORT, CVWEIGHT_SHORT, &
         CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
         X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
         MASS_ERROR_RELAX2_NOIT ) 

    NORMALISE = .false.
    IF( NORMALISE ) THEN
       ! make sure the FEM representation sums to unity so we dont get surprising results...
       DO CV_INOD = 1, CV_NONODS
          SUM = 0.0
          DO IPHASE = 1, NPHASE
             SUM = SUM + FEMT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS )
          END DO
          DO IPHASE = 1, NPHASE
             FEMT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMT( CV_INOD +( IPHASE - 1 ) &
                  * CV_NONODS ) / SUM
          END DO
       END DO
    ENDIF

    ! Calculate MEAN_PORE_CV   
    SUM_CV = 0.0
    MEAN_PORE_CV = 0.0
    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
          CV_INOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
          SUM_CV( CV_INOD )      = SUM_CV( CV_INOD ) + MASS_ELE( ELE )
          MEAN_PORE_CV( CV_INOD ) = MEAN_PORE_CV( CV_INOD ) + & 
               MASS_ELE( ELE ) * VOLFRA_PORE( ELE )
       END DO
    END DO
    MEAN_PORE_CV = MEAN_PORE_CV / SUM_CV


    ! For each node, find the largest and smallest value of T and 
    ! DENSITY for both the current and previous timestep, out of 
    ! the node value and all its surrounding nodes including Dirichlet b.c's.
    CALL SURRO_CV_MINMAX( TMAX, TMIN, TOLDMAX, TOLDMIN, DENMAX, DENMIN, DENOLDMAX, DENOLDMIN, &
         T2MAX, T2MIN, T2OLDMAX, T2OLDMIN, &
         T, TOLD, T2, T2OLD, DEN, DENOLD, IGOT_T2, NPHASE, CV_NONODS, NCOLACV, FINACV, COLACV, &
         STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC, SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
         WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
         WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
         T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
         DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD )


    ALLOCATE( FACE_ELE( NFACE, TOTELE ))
    ! Calculate FACE_ELE
    CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
         NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
         CV_SLOCLIST, X_NLOC, X_NDGLN )

    IF( GOT_DIFFUS ) THEN

       !       femt(:)=-(x(:)-0.5)*(x(:)-0.5)
       !       femt(:)=x(:)
       CALL DG_DERIVS( FEMT, FEMTOLD, &
            NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
                                !            NGI2, CV_NLOC, N2, WEIGHT2, N2, NLX2, NLY2, NLY2, &
            CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, CVFEN_SHORT, &
            CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
            X_NONODS, X, Y, Z,  &
            DTX_ELE, DTY_ELE, DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
            NFACE, FACE_ELE, CV_SLOCLIST, STOTEL, CV_SNLOC, WIC_T_BC, SUF_T_BC, &
            WIC_T_BC_DIRICHLET, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH )
    ENDIF

    !     =============== DEFINE THETA FOR TIME-STEPPING ===================
    ! Define the type of time integration:
    ! Timopt is 0 if CV_DISOPT is even (theta specified); 
    ! Timopt is 1 if CV_DISOPT is odd (non-linear theta scheme) 
    TIMOPT = MOD( CV_DISOPT, 2 )

    VTHETA = 1.0    
    FTHETA = CV_THETA    

    IF( GETCT ) THEN ! Obtain the CV discretised CT eqns plus RHS       
       CT_RHS = 0.0
       CT = 0.0
    ENDIF

    IF(GETCV_DISC) THEN ! Obtain the CV discretised advection/diffusion eqns
       CV_RHS = 0.0
       IF( GETMAT ) THEN
          ACV = 0.0
       ENDIF
    ENDIF

    GET_GTHETA=.FALSE.
    IF( IGOT_THETA_FLUX == 1 ) THEN
       IF( GET_THETA_FLUX ) THEN
          THETA_FLUX = 0.0
          ONE_M_THETA_FLUX = 0.0
          GET_GTHETA = .TRUE.
          THETA_GDIFF = 0.0
       ENDIF
    ENDIF

    ! Now we begin the loop over elements to assemble the advection terms
    ! into the matrix (ACV) and the RHS

    ewrite(3,*)'x,ftheta:'

    Loop_Elements: DO ELE = 1, TOTELE


       Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC ! Loop over the nodes of the element

          ! Global node number of the local node
          CV_NODI = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )

          ! Loop over quadrature (gauss) points in ELE neighbouring ILOC
          Loop_GCOUNT: DO GCOUNT = FINDGPTS( CV_ILOC ), FINDGPTS( CV_ILOC + 1 ) - 1, 1

             ! COLGPTS stores the local Gauss-point number in the ELE
             GI = COLGPTS( GCOUNT )

             ! Get the neighbouring node for node ILOC and Gauss point GI
             CV_JLOC = CV_NEILOC( CV_ILOC, GI )

             ELE2 = 0
             SELE = 0
             INTEGRAT_AT_GI=.TRUE.

             IF( CV_JLOC == -1 ) THEN

                ! We are on the boundary or next to another element.  Determine CV_OTHER_LOC
                ! CV_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
                ! Look for these nodes on the other elements.
                CALL FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC,  &
                     MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI,  &
                     TOTELE, X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                     CV_SNLOC, CV_ON_FACE, SCVNGI, GI, X_SHARE, X_NONODS, ELE, ELE2,  &
                     FINELE, COLELE, NCOLELE ) 

                IF(INTEGRAT_AT_GI) THEN
                   CV_JLOC = CV_OTHER_LOC( CV_ILOC )
                   SELE=0

                   IF( CV_JLOC == 0 ) THEN ! We are on the boundary of the domain
                      ewrite(3,*)'---here1.1'
                      CV_JLOC = CV_ILOC
                      ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC
                      CALL CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, CV_NGI, U_SLOC2LOC, CV_SLOC2LOC, &
                           FACE_ELE, TOTELE, NFACE, CV_ON_FACE, GI, &
                           CV_NONODS, LOG_ON_BOUND, CV_NLOC, U_NLOC, CV_SNLOC,U_SNLOC, STOTEL, &
                           CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN ) 
                      !           EWRITE(3,*)'*****AFTER CALC_SELE SELE,CV_SILOC,CV_SNLOC:',SELE,CV_SILOC,CV_SNLOC
                   ENDIF
                   INTEGRAT_AT_GI=.NOT.((ELE==ELE2).AND.(SELE==0))
                   !        ewrite(3,*)'---her2'
                ENDIF

             ENDIF


             ! avoid indegrating across the middle of a CV on the boundaries of elements        
             Conditional_integration: IF(INTEGRAT_AT_GI) THEN

                ! if necessary determine the derivatives between elements ELE and ELE2

                ! Calculate the control volume normals at the Gauss pts.
                !     EWRITE(3,*)'************CV_ILOC=',CV_ILOC
                CALL SCVDETNX( ELE,      GI,          &
                     X_NLOC,  SCVNGI,  TOTELE,  &
                     X_NDGLN,  X_NONODS,         &
                     SCVDETWEI, CVNORMX, CVNORMY, &  
                     CVNORMZ,  SCVFEN,     SCVFENLX,   &   
                     SCVFENLY, SCVFEWEIGH, XC_CV(CV_NODI),    & 
                     YC_CV(CV_NODI),     ZC_CV(CV_NODI),    X,       & 
                     Y,        Z,                &
                     D1,       D3,      DCYL )


                ! ================ COMPUTE THE FLUX ACROSS SUB-CV FACE ===============

                ! Find its global node number
                IF(ELE2==0) THEN
                   CV_NODJ = CV_NDGLN(( ELE - 1 )  * CV_NLOC + CV_JLOC )
                ELSE
                   CV_NODJ = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_JLOC )
                ENDIF
                X_NODI  =  X_NDGLN(( ELE - 1 ) * X_NLOC  + CV_ILOC )

                Conditional_GETCT1: IF( GETCT ) THEN ! Obtain the CV discretised CT equations plus RHS

                   DO U_KLOC = 1, U_NLOC
                      U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
                      JCOUNT = 0
                      DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1, 1
                         IF(COLCT( COUNT ) == U_NODK) JCOUNT = COUNT
                      END DO
                      JCOUNT_KLOC( U_KLOC ) = JCOUNT
                   END DO
                   IF((ELE2 /= 0).AND.(ELE2 /= ELE)) THEN
                      ewrite(3,*)'CV_NODI,CV_NODJ,ELE,ELE2:',CV_NODI,CV_NODJ,ELE,ELE2
                      DO U_KLOC = 1, U_NLOC
                         U_NODK = U_NDGLN(( ELE2 - 1 ) * U_NLOC + U_KLOC )
                         !             ewrite(3,*)'U_NODK:',U_NODK
                         JCOUNT = 0
                         DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1, 1
                            !           EWRITE(3,*)'COUNT,COLCT( COUNT ):',COUNT,COLCT( COUNT )
                            IF(COLCT( COUNT ) == U_NODK) JCOUNT = COUNT
                         END DO
                         !        STOP 56
                         JCOUNT_KLOC2( U_KLOC ) = JCOUNT
                      END DO

                   ENDIF

                ENDIF Conditional_GETCT1

                ! Compute the distance HDC between the nodes either side of the CV face 
                ! (this is needed to compute the local courant number and the non-linear
                ! theta)
                IF(SELE == 0) THEN
                   HDC = SQRT( (XC_CV(CV_NODI)-XC_CV(CV_NODJ))**2+(YC_CV(CV_NODI)-YC_CV(CV_NODJ))**2 &
                        +(ZC_CV(CV_NODI)-ZC_CV(CV_NODJ))**2  )
                ELSE
                   HDC = SQRT( (XC_CV(CV_NODI)-X(X_NODI))**2+(YC_CV(CV_NODI)-Y(X_NODI))**2 &
                        +(ZC_CV(CV_NODI)-Z(X_NODI))**2  )
                ENDIF

                ! get the sum of limiting functions correct...************
                SUM2ONE=.false.
                Conditional_SUMLimiting: IF(SUM2ONE) THEN
                   SUM_LIMT   =0.0
                   SUM_LIMTOLD=0.0
                   Loop_IPHASE5: DO IPHASE = 1, NPHASE

                      CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                      CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS

                      ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
                      IF(IGOT_T2==1) THEN 
                         CALL GET_INT_VEL( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                              HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                              T2, T2OLD, FEMT2, FEMT2OLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                              CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                              CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                              SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                              WIC_U_BC_DIRICHLET, &
                              UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                              ONE_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                              MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                              1, LIMT2, LIMT2OLD, FEMDGI, FEMT2GI, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                              T2MIN, T2MAX, T2OLDMIN, T2OLDMAX, T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                              IN_ELE_UPWIND, DG_ELE_UPWIND )
                      ELSE
                         CALL GET_INT_VEL( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                              HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                              T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                              CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                              CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                              SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                              WIC_U_BC_DIRICHLET, &
                              UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                              ONE_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                              MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                              1, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
                              TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
                              IN_ELE_UPWIND, DG_ELE_UPWIND )
                      ENDIF


                      !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                      ! Calculate T and DEN on the CV face at quadrature point GI.
                      CALL GET_INT_T_DEN( FVT, FVTOLD, FVT2, FVT2OLD, FVD, FVDOLD, LIMD, LIMT, LIMT2, &
                           LIMDOLD, LIMTOLD, LIMT2OLD, LIMDT, LIMDTOLD, LIMDTT2, LIMDTT2OLD,&
                           FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
                           CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                           CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, INCOMEOLD, &
                           T, TOLD, T2, T2OLD, DEN, DENOLD, FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, &
                           TMIN, TOLDMIN, T2MIN, T2OLDMIN, DENMIN, DENOLDMIN, &
                           TMAX, TOLDMAX, T2MAX, T2OLDMAX, DENMAX, DENOLDMAX, &
                           SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                           WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                           WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                           MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
                           DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD, &
                           T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2)

                      SUM_LIMT    = SUM_LIMT    + LIMT
                      SUM_LIMTOLD = SUM_LIMTOLD + LIMTOLD

                   END DO Loop_IPHASE5
                ENDIF Conditional_SUMLimiting
                ! get the sum of limiting functions correct...************

                Loop_IPHASE: DO IPHASE = 1, NPHASE

                   CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                   CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS

                   If_GOT_DIFFUS: IF( GOT_DIFFUS ) THEN
                      ! This sub caculates the effective diffusion coefficientd DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
                      CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX,  &
                           CV_NLOC, MAT_NLOC, CV_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
                           SCVFEN,SCVNGI,GI,IPHASE,NDIM,TDIFFUSION,HDC, T,TOLD,CV_NODJ_IPHA,CV_NODI_IPHA,ELE,ELE2, &
                           CVNORMX,CVNORMY,CVNORMZ,  &
                           DTX_ELE,DTY_ELE,DTZ_ELE,DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE, &
                           SELE,STOTEL,WIC_T_BC,WIC_T_BC_DIRICHLET, CV_OTHER_LOC,MAT_OTHER_LOC )
                   ELSE ! IF(GOT_DIFFUS) THEN...
                      DIFF_COEF_DIVDX    = 0.0
                      DIFF_COEFOLD_DIVDX = 0.0
                   END IF If_GOT_DIFFUS

                   NFACE_ITS=1
                   !                   IF((CV_ELE_TYPE==2).AND.(CV_NONODS==TOTELE*CV_NLOC)) NFACE_ITS=2
                   DO FACE_ITS = 1, NFACE_ITS
                      ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
                      IF(IGOT_T2==1) THEN
                         CALL GET_INT_VEL( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                              HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                              T2, T2OLD, FEMT2, FEMT2OLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                              CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                              CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                              SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                              WIC_U_BC_DIRICHLET, &
                              UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                              ONE_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                              MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                              FACE_ITS, LIMT2, LIMT2OLD, FEMDGI, FEMT2GI, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                              T2MIN, T2MAX, T2OLDMIN, T2OLDMAX, T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                              IN_ELE_UPWIND, DG_ELE_UPWIND )
                      ELSE
                         CALL GET_INT_VEL( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                              HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                              T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                              CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                              CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                              SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                              WIC_U_BC_DIRICHLET, &
                              UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                              ONE_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                              MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                              FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
                              TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
                              IN_ELE_UPWIND, DG_ELE_UPWIND )
                      ENDIF

                      !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                      ! Calculate T and DEN on the CV face at quadrature point GI.
                      CALL GET_INT_T_DEN( FVT, FVTOLD, FVT2, FVT2OLD, FVD, FVDOLD, LIMD, LIMT, LIMT2, &
                           LIMDOLD, LIMTOLD, LIMT2OLD, LIMDT, LIMDTOLD, LIMDTT2, LIMDTT2OLD,&
                           FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
                           CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                           CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, INCOMEOLD, &
                           T, TOLD, T2, T2OLD, DEN, DENOLD, FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, &
                           TMIN, TOLDMIN, T2MIN, T2OLDMIN, DENMIN, DENOLDMIN, &
                           TMAX, TOLDMAX, T2MAX, T2OLDMAX, DENMAX, DENOLDMAX, &
                           SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                           WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                           WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                           MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
                           DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD, &
                           T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2)

                   END DO

                   IF( SUM2ONE ) THEN
                      LIMT   =LIMT/SUM_LIMT
                      LIMTOLD=LIMTOLD/SUM_LIMTOLD
                      LIMDT=LIMT*LIMD
                      LIMDTOLD=LIMTOLD*LIMDOLD
                   ENDIF


                   ! Define face value of theta
                   IF(IGOT_T2==1) THEN
                      FTHETA=FACE_THETA( DT, CV_THETA, HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                           T(CV_NODJ_IPHA)*DEN(CV_NODJ_IPHA)*T2(CV_NODJ_IPHA), &
                           T(CV_NODI_IPHA)*DEN(CV_NODI_IPHA)*T2(CV_NODI_IPHA), &
                           NDOTQOLD, LIMDTT2OLD, DIFF_COEFOLD_DIVDX, &
                           TOLD(CV_NODJ_IPHA)*DENOLD(CV_NODJ_IPHA)*T2OLD(CV_NODJ_IPHA), &
                           TOLD(CV_NODI_IPHA)*DENOLD(CV_NODI_IPHA)*T2OLD(CV_NODI_IPHA) )
                   ELSE
                      FTHETA=FACE_THETA( DT, CV_THETA, HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                           T(CV_NODJ_IPHA)*DEN(CV_NODJ_IPHA), T(CV_NODI_IPHA)*DEN(CV_NODI_IPHA), &
                           NDOTQOLD, LIMDTT2OLD, DIFF_COEFOLD_DIVDX, &
                           TOLD(CV_NODJ_IPHA)*DENOLD(CV_NODJ_IPHA), &
                           TOLD(CV_NODI_IPHA)*DENOLD(CV_NODI_IPHA) )
                   ENDIF

                   FTHETA_T2         =FTHETA*LIMT2
                   ONE_M_FTHETA_T2OLD=(1.0-FTHETA)*LIMT2OLD


                   IF(IGOT_THETA_FLUX==1) THEN
                      IF(GET_THETA_FLUX) THEN
                         THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)      =FTHETA*LIMT
                         ONE_M_THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)=(1.0-FTHETA)*LIMTOLD
                      ENDIF
                      IF(USE_THETA_FLUX) THEN
                         FTHETA_T2         =THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)
                         ONE_M_FTHETA_T2OLD=ONE_M_THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)
                      ENDIF
                   ENDIF

                   ROBIN1=0.0
                   ROBIN2=0.0
                   IF(SELE /= 0) THEN
                      IF(WIC_T_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_ROBIN) THEN
                         ROBIN1=SUF_T_BC_ROB1((SELE-1)*CV_SNLOC+CV_SILOC +(IPHASE-1)*STOTEL*CV_SNLOC)
                         ROBIN2=SUF_T_BC_ROB2((SELE-1)*CV_SNLOC+CV_SILOC +(IPHASE-1)*STOTEL*CV_SNLOC)
                      ENDIF
                   ENDIF

                   !====================== ACV AND RHS ASSEMBLY ===================
                   Conditional_GETCT2 : IF( GETCT ) THEN ! Obtain the CV discretised CT eqations plus RHS

                      CALL PUT_IN_CT_RHS(CT, CT_RHS, U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
                           CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE, ELE, ELE2, SELE, &
                           JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC, U_NDGLN, NU, NV, NW,  &
                           SUFEN, SCVDETWEI, CVNORMX, CVNORMY, CVNORMZ, DENOLD, CV_NODI, CV_NODI_IPHA, &
                           UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
                           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                           NDOTQ, NDOTQOLD, LIMDT, LIMDTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD )

                   ENDIF Conditional_GETCT2

                   Conditional_GETCV_DISC: IF(GETCV_DISC) THEN 
                      ! Obtain the CV discretised advection/diffusion equations

                      IF(GETMAT) THEN
                         ! - Calculate the integration of the limited, high-order flux over a face  
                         ! Conservative discretisation. The matrix (PIVOT ON LOW ORDER SOLN)
                         IF( ( CV_NODI_IPHA /= CV_NODJ_IPHA ).AND.( CV_NODJ_IPHA /= 0 ) ) THEN
                            DO COUNT = FINACV( CV_NODI_IPHA ), FINACV( CV_NODI_IPHA + 1 ) - 1, 1
                               IF( COLACV( COUNT ) == CV_NODJ_IPHA )  JCOUNT_IPHA = COUNT
                            END DO

                            ACV( JCOUNT_IPHA ) =  ACV( JCOUNT_IPHA ) &
                                 +  FTHETA_T2 * SCVDETWEI( GI ) * NDOTQ * INCOME * LIMD  &  ! advection
                                 -  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX   ! Diffusion contribution
                            IF(GET_GTHETA) THEN
                               THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                    +  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX *T(CV_NODJ_IPHA) !Diffusion contribution
                            ENDIF
                         ELSE IF(WIC_T_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRICHLET) THEN
                            CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA )  &
                                 + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &
                                 *SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC + (IPHASE-1)*STOTEL*CV_SNLOC)
                            IF(GET_GTHETA) THEN
                               THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                    + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &
                                    *SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC + (IPHASE-1)*STOTEL*CV_SNLOC)
                            ENDIF
                         ENDIF

                         IMID_IPHA = MIDACV( CV_NODI_IPHA )

                         ACV( IMID_IPHA ) =  ACV( IMID_IPHA ) &
                              +  FTHETA_T2 * SCVDETWEI( GI ) * NDOTQ * ( 1. - INCOME ) * LIMD   & ! advection
                              +  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &  ! Diffusion contribution
                              +  SCVDETWEI( GI ) * ROBIN1  ! Robin bc
                         IF(GET_GTHETA) THEN
                            THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                 -  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX*T( CV_NODI_IPHA )  &!Diffusion contribution
                                 -  SCVDETWEI( GI ) * ROBIN1*T( CV_NODI_IPHA )  ! Robin bc
                         ENDIF

                         ! CV_BETA=0 for Non-conservative discretisation (CV_BETA=1 for conservative disc)
                         ACV( IMID_IPHA ) =  ACV( IMID_IPHA )  &
                              - FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQ * LIMD 
                      ENDIF

                      TMID   =     T( CV_NODI_IPHA )
                      TOLDMID = TOLD( CV_NODI_IPHA )

                      ! Make allowances for no matrix stencil operating from outside the boundary.
                      BCZERO=1.0
                      IF( (SELE /= 0) .AND. (INCOME > 0.5) ) BCZERO=0.0

                      ! Put results into the RHS vector
                      CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA )  &
                                ! subtract 1st order adv. soln.
                           + FTHETA_T2 * NDOTQ * SCVDETWEI( GI ) * LIMD * FVT * BCZERO & 
                           - SCVDETWEI( GI ) * (FTHETA_T2 * NDOTQ *LIMDT &
                           + ONE_M_FTHETA_T2OLD * NDOTQOLD *LIMDTOLD ) ! hi order adv

                      !  Subtract out 1st order term non-conservative adv.
                      CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA ) &
                           - FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQ * LIMD * TMID

                      ! High-order non-conservative advection contribution 
                      CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA ) &
                           +  ( 1. - CV_BETA) * SCVDETWEI( GI ) &
                           * ( FTHETA_T2 * NDOTQ *TMID * LIMD  &
                           +  ONE_M_FTHETA_T2OLD * NDOTQOLD *LIMDOLD * TOLDMID )

                      ! Diffusion contribution
                      CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA ) &
                           + (1.-FTHETA)*SCVDETWEI(GI)*DIFF_COEFOLD_DIVDX  &
                           *(TOLD(CV_NODJ_IPHA)-TOLD(CV_NODI_IPHA)) &
                                ! Robin bc
                           + SCVDETWEI(GI)*ROBIN2
                      IF(GET_GTHETA) THEN
                         THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                              + (1.-FTHETA)*SCVDETWEI(GI)*DIFF_COEFOLD_DIVDX  &
                              *(TOLD(CV_NODJ_IPHA)-TOLD(CV_NODI_IPHA)) &
                                ! Robin bc
                              + SCVDETWEI(GI)*ROBIN2
                      ENDIF


                   ENDIF Conditional_GETCV_DISC

                END DO Loop_IPHASE

             END IF Conditional_integration

          END DO Loop_GCOUNT

          !       IF(CV_ILOC.EQ.3) stop 39831
          !       IF((CV_ILOC.EQ.1).and.(ELE.EQ.2)) stop 39831
       END DO Loop_CV_ILOC

    END DO Loop_Elements



    IF(GET_GTHETA) THEN
       DO CV_NODI = 1, CV_NONODS
          DO IPHASE = 1, NPHASE    
             CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
             THETA_GDIFF(CV_NODI_IPHA) = THETA_GDIFF(CV_NODI_IPHA) / MASS_CV(CV_NODI)
          END DO
       END DO
    ENDIF


    Conditional_GETCV_DISC2: IF( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion equations

       ewrite(3,*)'before adding extra bits*****DEN:',DEN
       ewrite(3,*)'before adding extra bits*****DENOLD:',DENOLD
       ewrite(3,*)'before adding extra bits*****TOLD:',TOLD
       ewrite(3,*)'before adding extra bits*****MEAN_PORE_CV:',MEAN_PORE_CV

       !sourct2( 1 : cv_nonods ) = -.0 !-981. * ( 1.05 - .71 )
       !sourct2( cv_nonods + 1 : cv_nonods * nphase ) = -0.!-10.0e-1 !-981. * ( 1.05 - .71 )

       Loop_CVNODI2: DO CV_NODI = 1, CV_NONODS! Put onto the diagonal of the matrix 

          Loop_IPHASE2: DO IPHASE = 1, NPHASE    
             CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
             ! For the gravity term
!!!! SOURCT2( CV_NODI_IPHA ) = SOURCT( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA ) 
!!!! SOURCT2( CV_NODI_IPHA ) = SOURCT2( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA ) 
             IMID_IPHA = MIDACV( CV_NODI_IPHA )

             IF(IGOT_T2==1) THEN
                CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                     + MASS_CV( CV_NODI ) * SOURCT( CV_NODI_IPHA ) &
                     - CV_BETA * MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) & ! conservative time term. 
                     * TOLD(CV_NODI_IPHA)  &
                     *(DEN(CV_NODI_IPHA)*T2(CV_NODI_IPHA)-DENOLD(CV_NODI_IPHA)*T2OLD(CV_NODI_IPHA))/DT

                ACV( IMID_IPHA ) =  ACV( IMID_IPHA ) &
                     + (CV_BETA *DENOLD(CV_NODI_IPHA)*T2OLD(CV_NODI_IPHA) &
                     + (1.-CV_BETA) *DEN(CV_NODI_IPHA)*T2(CV_NODI_IPHA))  &
                     * MEAN_PORE_CV( CV_NODI ) *MASS_CV( CV_NODI )/DT
                ewrite(3,*)'sec time - CV_RHS::',CV_NODI_IPHA, CV_RHS( CV_NODI_IPHA )

                CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                     + (CV_BETA *DENOLD(CV_NODI_IPHA)*T2OLD(CV_NODI_IPHA) &
                     + (1.-CV_BETA) *DEN(CV_NODI_IPHA)*T2(CV_NODI_IPHA))  &
                     * MEAN_PORE_CV( CV_NODI ) *MASS_CV( CV_NODI ) * TOLD(CV_NODI_IPHA)/DT
             ELSE
                CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                     + MASS_CV( CV_NODI ) * SOURCT( CV_NODI_IPHA ) &
                     - CV_BETA * MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) & ! conservative time term. 
                     * TOLD(CV_NODI_IPHA)*(DEN(CV_NODI_IPHA)-DENOLD(CV_NODI_IPHA))/DT

                ACV( IMID_IPHA ) =  ACV( IMID_IPHA ) &
                     + (CV_BETA *DENOLD(CV_NODI_IPHA) + (1.-CV_BETA) *DEN(CV_NODI_IPHA))  &
                     * MEAN_PORE_CV( CV_NODI ) *MASS_CV( CV_NODI )/DT

                CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                     + (CV_BETA *DENOLD(CV_NODI_IPHA) + (1.-CV_BETA) *DEN(CV_NODI_IPHA))  &
                     * MEAN_PORE_CV( CV_NODI ) *MASS_CV( CV_NODI ) * TOLD(CV_NODI_IPHA)/DT
             ENDIF


             Conditional_GETMAT2: IF( GETMAT ) THEN

                DO JPHASE = 1, NPHASE   
                   CV_NODI_JPHA = CV_NODI + ( JPHASE - 1 ) * CV_NONODS
                   DO COUNT = FINACV( CV_NODI_IPHA ), FINACV( CV_NODI_IPHA + 1 ) - 1, 1
                      IF( CV_NODI_JPHA == COLACV( COUNT )) JCOUNT_IPHA = COUNT
                   END DO
                   ACV(JCOUNT_IPHA) = ACV(JCOUNT_IPHA) &
                        +  MASS_CV( CV_NODI ) * ABSORBT( CV_NODI, IPHASE, JPHASE )
                END DO

             ENDIF Conditional_GETMAT2

          END DO Loop_IPHASE2

       END DO Loop_CVNODI2


    ENDIF Conditional_GETCV_DISC2






    IF(GETCT) THEN! Form the rhs of the discretised equations
       ! Determine CV pressure CV_P

       DIAG_SCALE_PRES = 0.0 ! Obtain the CV discretised CT eqations plus RHS

       DO IPHASE = 1, NPHASE    
          DO CV_NODI = 1, CV_NONODS
             CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS

             CT_RHS( CV_NODI ) = CT_RHS( CV_NODI )  -  MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) &
                  * TOLD( CV_NODI_IPHA )*( DEN( CV_NODI_IPHA ) - DENOLD( CV_NODI_IPHA )  &
                  - DERIV( CV_NODI_IPHA ) * CV_P( CV_NODI )) / ( DT * DENOLD( CV_NODI_IPHA ))   

             DIAG_SCALE_PRES( CV_NODI ) = DIAG_SCALE_PRES( CV_NODI ) + &
                  MEAN_PORE_CV( CV_NODI ) * TOLD( CV_NODI_IPHA ) * DERIV( CV_NODI_IPHA )  &
                  / ( DT * DENOLD( CV_NODI_IPHA ))

             CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) + MASS_CV( CV_NODI ) * SOURCT( CV_NODI_IPHA ) 
             !CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) + MASS_CV( CV_NODI ) * SOURCT2( CV_NODI_IPHA ) 

             DO JPHASE = 1, NPHASE   
                CV_NODI_JPHA = CV_NODI + ( JPHASE - 1 ) * CV_NONODS
                CT_RHS( CV_NODI ) = CT_RHS( CV_NODI )  &
                     - MASS_CV( CV_NODI ) * ABSORBT( CV_NODI, IPHASE, JPHASE )* T( CV_NODI_JPHA )
             END DO
          END DO
       END DO

    ENDIF


    ewrite(3,*)'upwind fraction:'
    do iphase=1,nphase
       ewrite(3,*)'for phase iphase=',iphase
       do ele=1,totele-1
          do cv_iloc=1,cv_nloc
             cv_nodi = cv_ndgln((ele-1)*cv_nloc+cv_iloc)
             mat_nodi = mat_ndgln((ele-1)*cv_nloc+cv_iloc)
             cv_nodi_IPHA=cv_nodi +(IPHASE-1)*CV_NONODS

             if(cv_nonods==x_nonods) then
                ewrite(3,*)0.5*(x(cv_nodi)+x(cv_nodi+1)),UP_WIND_NOD(cv_nodi_IPHA)
             else
                if(cv_iloc==cv_nloc) then
                   ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),  &
                        UP_WIND_NOD(cv_nodi_IPHA)
                else
                   ewrite(3,*)0.5*(x(x_ndgln((ele-1)*x_nloc+cv_iloc))+x(x_ndgln((ele-1)*x_nloc+cv_iloc+1))),  &
                        UP_WIND_NOD(cv_nodi_IPHA)
                endif
             endif

          end do
       end do
    end do

    ! for the output
    T_FEMT = femt
    DEN_FEMT = FEMDEN

    ewrite(3,*)'IN cv_adv_dif a CV representation t:'
    CALL PRINT_CV_DIST(CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE, &
         T, X_NDGLN, CV_NDGLN, X) 
    ewrite(3,*) 'just print out - in cv_assemb'
    ! Deallocating temporary working arrays
    DEALLOCATE( JCOUNT_KLOC )
    DEALLOCATE( CVNORMX )
    DEALLOCATE( CVNORMY )
    DEALLOCATE( CVNORMZ )
    DEALLOCATE( COLGPTS ) !The size of this vector is over-estimated
    DEALLOCATE( FINDGPTS )
    DEALLOCATE( SNDOTQ )
    DEALLOCATE( SNDOTQOLD )
    DEALLOCATE( CV_ON_FACE )
    DEALLOCATE( U_ON_FACE )
    DEALLOCATE( CV_OTHER_LOC )
    DEALLOCATE( U_OTHER_LOC )
    DEALLOCATE( MAT_OTHER_LOC )
    DEALLOCATE( X_SHARE )

    DEALLOCATE( CVWEIGHT )
    DEALLOCATE( CVN )
    DEALLOCATE( CVFEN )
    DEALLOCATE( CVFENLX )
    DEALLOCATE( CVFENLY )
    DEALLOCATE( CVFENLZ )

    DEALLOCATE( CVWEIGHT_SHORT )
    DEALLOCATE( CVN_SHORT )
    DEALLOCATE( CVFEN_SHORT )
    DEALLOCATE( CVFENLX_SHORT )
    DEALLOCATE( CVFENLY_SHORT )
    DEALLOCATE( CVFENLZ_SHORT )

    DEALLOCATE( UFEN ) 
    DEALLOCATE( UFENLX )
    DEALLOCATE( UFENLY )
    DEALLOCATE( UFENLZ )

    DEALLOCATE( SCVFEN )
    DEALLOCATE( SCVFENSLX )
    DEALLOCATE( SCVFENSLY )
    DEALLOCATE( SCVFENLX )
    DEALLOCATE( SCVFENLY )
    DEALLOCATE( SCVFENLZ )
    DEALLOCATE( SCVFEWEIGH )

    DEALLOCATE( SUFEN )
    DEALLOCATE( SUFENSLX )
    DEALLOCATE( SUFENSLY )
    DEALLOCATE( SUFENLX )
    DEALLOCATE( SUFENLY )
    DEALLOCATE( SUFENLZ )

    DEALLOCATE( SCVDETWEI )
    DEALLOCATE( SRA )
    DEALLOCATE( LOG_ON_BOUND )

    DEALLOCATE( SBCVFEN )
    DEALLOCATE( SBCVFENSLX )
    DEALLOCATE( SBCVFENSLY )
    DEALLOCATE( SBCVFEWEIGH )
    DEALLOCATE( SBCVFENLX )
    DEALLOCATE( SBCVFENLY )
    DEALLOCATE( SBCVFENLZ )
    DEALLOCATE( SBUFEN )
    DEALLOCATE( SBUFENSLX )
    DEALLOCATE( SBUFENSLY )
    DEALLOCATE( SBUFENLX )
    DEALLOCATE( SBUFENLY )
    DEALLOCATE( SBUFENLZ )

    DEALLOCATE( CV_SLOC2LOC )
    DEALLOCATE( U_SLOC2LOC ) 
    DEALLOCATE( CV_SLOCLIST )
    DEALLOCATE( U_SLOCLIST )
    DEALLOCATE( CV_NEILOC )
    DEALLOCATE( TMIN )
    DEALLOCATE( TMAX )
    DEALLOCATE( TOLDMIN )
    DEALLOCATE( TOLDMAX )
    DEALLOCATE( DENMIN )
    DEALLOCATE( DENMAX )
    DEALLOCATE( DENOLDMIN )
    DEALLOCATE( DENOLDMAX )
    DEALLOCATE( FEMT )
    DEALLOCATE( FEMTOLD )
    DEALLOCATE( FEMDEN )
    DEALLOCATE( FEMDENOLD )
    DEALLOCATE( MASS_CV )
    DEALLOCATE( XC_CV )
    DEALLOCATE( YC_CV )
    DEALLOCATE( ZC_CV )
    DEALLOCATE( DTX_ELE )
    DEALLOCATE( DTY_ELE )
    DEALLOCATE( DTZ_ELE )
    DEALLOCATE( DTOLDX_ELE )
    DEALLOCATE( DTOLDY_ELE )
    DEALLOCATE( DTOLDZ_ELE )
    DEALLOCATE( FACE_ELE )

    ewrite(3,*) 'Leaving CV_ASSEMB'

  END SUBROUTINE CV_ASSEMB




  SUBROUTINE FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC,  &
       MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI, &
       TOTELE, X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
       CV_SNLOC, CV_ON_FACE, SCVNGI, GI, X_SHARE, X_NONODS, ELE, ELE2,  &
       FINELE, COLELE, NCOLELE) 
    ! We are on the boundary or next to another element. Determine CV_OTHER_LOC,
    ! U_OTHER_LOC. 
    ! CV_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
    ! Look for these nodes on the other elements. 
    ! ELE2=0 alo when we are between elements but are trying to integrate across 
    ! the middle of a CV. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NLOC, CV_NODI, U_NLOC, X_NONODS, TOTELE, X_NLOC, XU_NLOC, &
         CV_SNLOC, MAT_NLOC, SCVNGI, GI, ELE, NCOLELE 
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * XU_NLOC ), intent( in ) :: XU_NDGLN
    INTEGER, DIMENSION( CV_NLOC ), intent( inout ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( U_NLOC ), intent( inout ) :: U_OTHER_LOC
    INTEGER, DIMENSION( MAT_NLOC ), intent( inout ) :: MAT_OTHER_LOC
    LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: CV_ON_FACE
    LOGICAL, DIMENSION( X_NONODS ), intent( inout ) :: X_SHARE
    INTEGER, intent( inout ) :: ELE2    
    LOGICAL, intent( inout ) :: INTEGRAT_AT_GI  
    INTEGER, DIMENSION( TOTELE  +  1 ), intent( in ) :: FINELE
    INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
    ! Local variables
    INTEGER :: X_KLOC, X_NODK, X_NODK2, COUNT, ELE3, SUF_COUNT, CV_KLOC, CV_KLOC2, &
         U_KLOC, U_KLOC2, CV_NODK, XU_NODK, XU_NODK2, ILEV, JLEV

    !ewrite(3,*) 'In FIND_OTHER_SIDE'

    DO X_KLOC = 1, X_NLOC
       X_NODK = X_NDGLN(( ELE - 1) * X_NLOC + X_KLOC )
       X_SHARE( X_NODK ) = CV_ON_FACE( X_KLOC, GI )
    END DO

    ELE3 = 0
    DO COUNT = FINELE( ELE ), FINELE( ELE + 1 ) - 1, 1
       ELE2 = COLELE( COUNT )
       !       ewrite(3,*)'ele2',ele2

       ! See if we share the same nodes
       SUF_COUNT = 0
       IF(ELE2.NE.ELE) THEN
          DO X_KLOC = 1, X_NLOC
             X_NODK = X_NDGLN(( ELE2 - 1 ) * X_NLOC + X_KLOC )
             IF( X_SHARE( X_NODK )) SUF_COUNT = SUF_COUNT + 1
          END DO
       ENDIF

       IF( SUF_COUNT == CV_SNLOC ) ELE3 = ELE2

    END DO

    ELE2 = ELE3
    DO X_KLOC = 1, X_NLOC
       X_NODK = X_NDGLN(( ELE - 1 ) * X_NLOC + X_KLOC )
       X_SHARE( X_NODK ) = .FALSE.
    END DO

    IF(ELE2 /= 0) THEN 
       ! Is CV_NODI in element ELE2 if yes set ELE2=0 as we dont want to integrate in 
       ! the middle of a CV. 
       DO CV_KLOC = 1, CV_NLOC
          CV_NODK = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_KLOC )
          IF(CV_NODK.EQ.CV_NODI) INTEGRAT_AT_GI = .FALSE.
       END DO
    ENDIF

    IF((ELE2 /= 0).AND.INTEGRAT_AT_GI) THEN ! Determine CV_OTHER_LOC(CV_KLOC)
       CV_OTHER_LOC = 0
       DO CV_KLOC = 1, CV_NLOC

          IF( CV_ON_FACE( CV_KLOC, GI )) THEN ! Find opposite local node
             X_NODK = X_NDGLN(( ELE - 1 ) * X_NLOC + CV_KLOC )

             DO CV_KLOC2 = 1, CV_NLOC
                X_NODK2 = X_NDGLN(( ELE2 - 1 ) * X_NLOC + CV_KLOC2 )
                IF( X_NODK2 == X_NODK ) CV_OTHER_LOC( CV_KLOC ) = CV_KLOC2
             END DO

          ENDIF

       END DO

       U_OTHER_LOC = 0 ! Determine U_OTHER_LOC(U_KLOC)
       IF(XU_NLOC /= U_NLOC) THEN ! This is for the overlapping approach
          DO U_KLOC = 1, XU_NLOC ! Find opposite local node
             XU_NODK = XU_NDGLN(( ELE - 1 ) * XU_NLOC + U_KLOC )

             DO U_KLOC2 = 1, XU_NLOC
                XU_NODK2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_KLOC2 )
                ! XU_NLOC==1 is a special case...
                IF(( XU_NODK2 ==XU_NODK ).OR.(XU_NLOC==1)) THEN
                   DO ILEV=1,CV_NLOC
                      JLEV=CV_OTHER_LOC( ILEV )
                      IF(JLEV /= 0) THEN 
                         U_OTHER_LOC( U_KLOC +(ILEV-1)*XU_NLOC) = U_KLOC2 +(JLEV-1)*XU_NLOC
                      ENDIF
                   END DO
                END IF
             END DO
          END DO
          !   if((ele==2).and.(ele2==1)) then
          !       ewrite(3,*)'U_OTHER_LOC:',U_OTHER_LOC
          !       stop 32823
          !   endif
       ELSE
          DO U_KLOC = 1, U_NLOC ! Find opposite local node
             XU_NODK = XU_NDGLN(( ELE - 1 ) * XU_NLOC + U_KLOC )

             DO U_KLOC2 = 1, U_NLOC
                XU_NODK2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_KLOC2 )
                IF(( XU_NODK2 ==XU_NODK ).OR.(XU_NLOC==1)) U_OTHER_LOC( U_KLOC ) = U_KLOC2
             END DO
          END DO
       ENDIF

       MAT_OTHER_LOC = CV_OTHER_LOC
    ELSE
       CV_OTHER_LOC =0
       U_OTHER_LOC  =0
       MAT_OTHER_LOC=0
    ENDIF
    !!    ewrite(3,*)'CV_OTHER_LOC:',CV_OTHER_LOC
    !!    ewrite(3,*)'u_OTHER_LOC:',u_OTHER_LOC
    !    stop

    !ewrite(3,*) 'Leaving FIND_OTHER_SIDE'

    RETURN

  END SUBROUTINE FIND_OTHER_SIDE




  SUBROUTINE PROJ_CV_TO_FEM_4( FEMT, FEMTOLD, FEMDEN, FEMDENOLD, T, TOLD, DEN, DENOLD, &
       IGOT_T2,T2,T2OLD, FEMT2,FEMT2OLD, &
       XC_CV,YC_CV,ZC_CV, MASS_CV, MASS_ELE, &
       NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
       CV_NGI, CV_NLOC, CVN, CVWEIGHT, N, NLX, NLY, NLZ, &
       X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
       MASS_ERROR_RELAX2_NOIT )

    ! determine FEMT (finite element wise) etc from T (control volume wise) 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, NPHASE, CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
         X_NONODS, NCOLM, IGOT_T2
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
    REAL, DIMENSION( CV_NONODS*NPHASE ), intent( inout ) :: FEMT, FEMTOLD, FEMDEN, FEMDENOLD
    REAL, DIMENSION( CV_NONODS*NPHASE*IGOT_T2 ), intent( inout ) :: FEMT2, FEMT2OLD
    REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: T, TOLD, DEN, DENOLD
    REAL, DIMENSION( CV_NONODS*NPHASE*IGOT_T2 ), intent( in ) :: T2, T2OLD
    REAL, DIMENSION( CV_NONODS ), intent( inout ) :: MASS_CV, XC_CV,YC_CV,ZC_CV
    REAL, DIMENSION( TOTELE ), intent( inout ) :: MASS_ELE
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in ) :: CVN
    REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in ) :: N, NLX, NLY, NLZ 
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
    INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
    INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
    REAL, DIMENSION( 5 ), intent( in ) :: MASS_ERROR_RELAX2_NOIT
    ! Local variables
    REAL, DIMENSION( : ), allocatable :: PSI, FEMPSI, PSI_AVE, PSI_INT
    INTEGER :: NTSOL,NTSOL_AVE,NTSOL_INT,ELE,CV_ILOC,X_INOD,CV_INOD,NL,NFIELD

    ewrite(3,*) 'In PROJ_CV_TO_FEM_4'

    NFIELD=4 + IGOT_T2*2
    ALLOCATE( PSI( NFIELD * CV_NONODS*NPHASE ))
    ALLOCATE( FEMPSI( NFIELD * CV_NONODS*NPHASE ))
    ALLOCATE( PSI_AVE( 3 * CV_NONODS ))
    ALLOCATE( PSI_INT( 1 * CV_NONODS ))

    NTSOL = NFIELD*NPHASE
    NTSOL_AVE = 3
    NTSOL_INT = 1

    NL=CV_NONODS*NPHASE
    PSI( 1 + 0 * NL : NL + 0 * NL )  =      T( 1 : NL ) 
    PSI( 1 + 1 * NL : NL + 1 * NL )  =   TOLD( 1 : NL ) 
    PSI( 1 + 2 * NL : NL + 2 * NL )  =    DEN( 1 : NL ) 
    PSI( 1 + 3 * NL : NL + 3 * NL )  = DENOLD( 1 : NL ) 
    IF(IGOT_T2==1) THEN
       PSI( 1 + 4 * NL : NL + 4 * NL )  =      T2( 1 : NL ) 
       PSI( 1 + 5 * NL : NL + 5 * NL )  =   T2OLD( 1 : NL ) 
    ENDIF

    DO ELE=1,TOTELE
       DO CV_ILOC=1,CV_NLOC
          X_INOD = X_NDGLN((ELE-1)*X_NLOC +CV_ILOC)
          CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
          PSI_AVE(CV_INOD)            =X(X_INOD)
          PSI_AVE(CV_INOD+CV_NONODS)  =Y(X_INOD)
          PSI_AVE(CV_INOD+2*CV_NONODS)=Z(X_INOD)
       END DO
    END DO
    PSI_INT=1.0

    FEMPSI = PSI

    CALL PROJ_CV_TO_FEM( FEMPSI, PSI, NTSOL, NDIM, &
         PSI_AVE,NTSOL_AVE, PSI_INT,NTSOL_INT, MASS_ELE, &
         CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI, CV_NLOC, CVN, CVWEIGHT, N, NLX, NLY, NLZ, &
         X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
         MASS_ERROR_RELAX2_NOIT )

    NL=CV_NONODS*NPHASE
    FEMT( 1 : NL ) = FEMPSI( 1 + 0 * NL : NL + 0 * NL ) 
    FEMTOLD( 1 : NL ) = FEMPSI( 1 + 1 * NL : NL + 1 * NL) 
    FEMDEN( 1 : NL ) = FEMPSI( 1 + 2 * NL : NL + 2 * NL ) 
    FEMDENOLD( 1 : NL ) = FEMPSI( 1 + 3 * NL : NL + 3*NL) 
    IF(IGOT_T2==1) THEN 
       FEMT2( 1 : NL ) = FEMPSI( 1 + 4 * NL : NL + 4 * NL ) 
       FEMT2OLD( 1 : NL ) = FEMPSI( 1 + 5 * NL : NL + 5*NL) 
    ENDIF

    XC_CV( 1 : CV_NONODS ) = PSI_AVE( 1 : CV_NONODS )
    YC_CV( 1 : CV_NONODS ) = PSI_AVE( 1 +CV_NONODS:   2*CV_NONODS )
    ZC_CV( 1 : CV_NONODS ) = PSI_AVE( 1 +2*CV_NONODS: 3*CV_NONODS )
    MASS_CV( 1 : CV_NONODS ) = PSI_INT( 1 : CV_NONODS )

    DEALLOCATE( PSI )
    DEALLOCATE( FEMPSI )
    DEALLOCATE( PSI_AVE )
    DEALLOCATE( PSI_INT )

    ewrite(3,*) 'Leaving PROJ_CV_TO_FEM_4'

    RETURN

  END SUBROUTINE PROJ_CV_TO_FEM_4




  SUBROUTINE PROJ_CV_TO_FEM( FEMPSI, PSI, NTSOL, NDIM, &
       PSI_AVE, NTSOL_AVE, PSI_INT, NTSOL_INT, MASS_ELE, &
       CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
       CV_NGI, CV_NLOC, CVN, CVWEIGHT, N, NLX, NLY, NLZ, &
       X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
       MASS_ERROR_RELAX2_NOIT )

    ! Determine FEMT (finite element wise) etc from T (control volume wise)
    ! Also integrate PSI_INT over each CV and avergae PSI_AVE over each CV. 
    use shape_functions 
    use solvers_module
    use spact
    IMPLICIT NONE
    INTEGER, intent( in ) :: NTSOL, NTSOL_AVE, NTSOL_INT, NDIM, CV_NONODS, TOTELE, &
         X_NLOC, CV_NGI, CV_NLOC, X_NONODS, NCOLM     
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
    REAL, DIMENSION( CV_NONODS * NTSOL ), intent( inout ) :: FEMPSI
    REAL, DIMENSION( CV_NONODS * NTSOL ), intent( in ) :: PSI
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in ) :: CVN
    REAL, DIMENSION( CV_NONODS * NTSOL_AVE ), intent( inout ) :: PSI_AVE
    REAL, DIMENSION( CV_NONODS * NTSOL_INT ), intent( inout ) :: PSI_INT
    REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in ) :: N, NLX, NLY, NLZ 
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( TOTELE ), intent( inout ) :: MASS_ELE
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
    INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
    INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
    REAL, DIMENSION( 5 ), intent( in ) :: MASS_ERROR_RELAX2_NOIT
    ! Local variables
    LOGICAL :: D1, D3, DCYL
    REAL, DIMENSION( : ), allocatable :: MASS_CV, FEMPSI_RHS, &
         PSI_AVE2, PSI_INT2, MAT, DETWEI, RA
    REAL, DIMENSION( :, : ), allocatable :: NX, NY, NZ
    REAL :: VOLUME, NN, NM, MN, MM
    INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_GI, COUNT, IT

    ewrite(3,*) 'In PROJ_CV_TO_FEM'

    ALLOCATE( MASS_CV( CV_NONODS )) 
    ALLOCATE( FEMPSI_RHS( CV_NONODS * NTSOL )) 
    ALLOCATE( PSI_AVE2( CV_NONODS * NTSOL_AVE )) 
    ALLOCATE( PSI_INT2( CV_NONODS * NTSOL_INT )) 
    ALLOCATE( MAT( NCOLM ))
    ALLOCATE( DETWEI( CV_NGI )) 
    ALLOCATE( RA( CV_NGI ))
    ALLOCATE( NX( CV_NLOC, CV_NGI ))
    ALLOCATE( NY( CV_NLOC, CV_NGI ))
    ALLOCATE( NZ( CV_NLOC, CV_NGI ))

    PSI_AVE2 = PSI_AVE
    PSI_INT2 = PSI_INT
    PSI_AVE = 0.0
    PSI_INT = 0.0

    FEMPSI_RHS = 0.0
    MAT = 0.0
    MASS_CV = 0.0

    D1 = ( NDIM == 1 )
    D3 = ( NDIM == 3 )
    DCYL = .FALSE. 

    Loop_Elements: DO ELE = 1, TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNLXR( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, CV_NLOC, CV_NGI, &
            N, NLX, NLY, NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
            NX, NY, NZ ) 
       MASS_ELE( ELE ) = VOLUME

       Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC

          CV_NODI = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
          !    ewrite(3,*)'ele,CV_NODI,CV_ILOC:',ele,CV_NODI,CV_ILOC, x(CV_NODI)

          Loop_CV_JLOC: DO CV_JLOC = 1, CV_NLOC

             CV_NODJ = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )
             NN = 0.0
             NM = 0.0
             MN = 0.0
             MM = 0.0
             DO CV_GI = 1, CV_NGI
                NN = NN + N( CV_ILOC, CV_GI )   * N(   CV_JLOC, CV_GI ) * DETWEI( CV_GI )
                NM = NM + N( CV_ILOC, CV_GI )   * CVN( CV_JLOC, CV_GI ) * DETWEI( CV_GI )
                MN = MN + CVN( CV_ILOC, CV_GI ) * N( CV_JLOC, CV_GI )   * DETWEI( CV_GI )
                MM = MM + CVN( CV_ILOC, CV_GI ) * CVN( CV_JLOC, CV_GI ) * DETWEI( CV_GI )
             END DO

             CALL POSINMAT( COUNT, CV_NODI, CV_NODJ, CV_NONODS, FINDM, COLM, NCOLM )

             MAT( COUNT ) = MAT( COUNT ) + NN
             MASS_CV( CV_NODI ) = MASS_CV( CV_NODI ) + MM

             DO IT = 1, NTSOL
                FEMPSI_RHS( CV_NODI + ( IT - 1 ) * CV_NONODS ) = FEMPSI_RHS( CV_NODI +( IT - 1 ) &
                     * CV_NONODS ) + NM * PSI( CV_NODJ + ( IT - 1 ) * CV_NONODS )
             END DO
             DO IT = 1, NTSOL_AVE
                PSI_AVE( CV_NODI + ( IT - 1 ) * CV_NONODS ) = PSI_AVE( CV_NODI +( IT - 1 ) &
                     * CV_NONODS ) + MN * PSI_AVE2( CV_NODJ + ( IT - 1) * CV_NONODS )
             END DO
             DO IT = 1, NTSOL_INT
                PSI_INT( CV_NODI + ( IT - 1 ) * CV_NONODS ) = PSI_INT( CV_NODI +( IT - 1 ) &
                     * CV_NONODS ) + MN * PSI_INT2( CV_NODJ + ( IT - 1) * CV_NONODS )
             END DO

          END DO Loop_CV_JLOC

       END DO Loop_CV_ILOC

    END DO Loop_Elements
    !    ewrite(3,*)'before solving for averages' 
    !   ewrite(3,*)'psi_int: ',psi_int

    ! Form average...
    DO CV_NODI=1,CV_NONODS
       DO IT = 1, NTSOL_AVE
          PSI_AVE( CV_NODI + ( IT - 1 ) * CV_NONODS ) = PSI_AVE( CV_NODI +( IT - 1 ) &
               * CV_NONODS ) / MASS_CV( CV_NODI )
       END DO
    END DO
    !    ewrite(3,*)'psi_ave: ',PSI_AVE

    ! Solve...
    DO IT = 1, NTSOL
       CALL SOLVER( MAT,  FEMPSI( 1 + (IT - 1 ) * CV_NONODS ),  FEMPSI_RHS( 1 + ( IT - 1 ) * CV_NONODS ),  &
            NCOLM, CV_NONODS, FINDM, COLM, MIDM,   &
            MASS_ERROR_RELAX2_NOIT(1), MASS_ERROR_RELAX2_NOIT(2), MASS_ERROR_RELAX2_NOIT(3), &
            MASS_ERROR_RELAX2_NOIT(4), INT(MASS_ERROR_RELAX2_NOIT(5)+0.1) )
       !            1.0, 0.0, 1.0, 200 )
    END DO

    DEALLOCATE( MASS_CV )
    DEALLOCATE( FEMPSI_RHS )
    DEALLOCATE( PSI_AVE2 )
    DEALLOCATE( PSI_INT2 )
    DEALLOCATE( MAT )
    DEALLOCATE( DETWEI )
    DEALLOCATE( RA )
    DEALLOCATE( NX )
    DEALLOCATE( NY )
    DEALLOCATE( NZ )

    ewrite(3,*) 'Leaving PROJ_CV_TO_FEM'

    RETURN

  END SUBROUTINE PROJ_CV_TO_FEM





  SUBROUTINE DG_DERIVS_UVW( U, UOLD, V, VOLD, W, WOLD, &
       NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, X_NLOC, X_NDGLN, &
       CV_NGI, U_NLOC, CVWEIGHT, N, NLX, NLY, NLZ, &
       X_NONODS, X, Y, Z, &
       DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
       DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
       DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
       NFACE, FACE_ELE, U_SLOCLIST, STOTEL, U_SNLOC, WIC_U_BC,  &
       SUF_U_BC,SUF_V_BC,SUF_W_BC, &
       WIC_U_BC_DIRICHLET, SBCVNGI,SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH ) 

    ! determine FEMT (finite element wise) etc from T (control volume wise)
    use shape_functions 
    use matrix_operations
    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, NPHASE, U_NONODS, TOTELE, X_NLOC, CV_NGI, U_NLOC, &
         X_NONODS, STOTEL, U_SNLOC, WIC_U_BC_DIRICHLET, SBCVNGI, NFACE
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, UOLD, V, VOLD, W, WOLD
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( STOTEL*NPHASE ), intent( in ) ::  WIC_U_BC
    REAL, DIMENSION( STOTEL*U_SNLOC*NPHASE ), intent( in ) ::  SUF_U_BC,SUF_V_BC,SUF_W_BC
    INTEGER, DIMENSION( NFACE,U_SNLOC ), intent( in ) ::  U_SLOCLIST
    INTEGER, DIMENSION( NFACE,TOTELE ), intent( in ) ::  FACE_ELE
    REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( U_NLOC, CV_NGI ), intent( in ) :: N, NLX, NLY, NLZ 
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( U_NLOC, NPHASE, TOTELE ), intent( inout ) :: DUX_ELE,DUY_ELE,DUZ_ELE, &
         DUOLDX_ELE,DUOLDY_ELE,DUOLDZ_ELE
    REAL, DIMENSION( U_NLOC, NPHASE, TOTELE ), intent( inout ) :: DVX_ELE,DVY_ELE,DVZ_ELE, &
         DVOLDX_ELE,DVOLDY_ELE,DVOLDZ_ELE
    REAL, DIMENSION( U_NLOC, NPHASE, TOTELE ), intent( inout ) :: DWX_ELE,DWY_ELE,DWZ_ELE, &
         DWOLDX_ELE,DWOLDY_ELE,DWOLDZ_ELE
    REAL, DIMENSION( U_SNLOC, SBCVNGI ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
    REAL, DIMENSION( SBCVNGI ), intent( in ) :: SBWEIGH

    CALL DG_DERIVS( U, UOLD, &
         NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI, U_NLOC, CVWEIGHT, N, NLX, NLY, NLZ, &
         X_NONODS, X, Y, Z, &
         DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
         NFACE, FACE_ELE, U_SLOCLIST, STOTEL, U_SNLOC, WIC_U_BC, SUF_U_BC, &
         WIC_U_BC_DIRICHLET, SBCVNGI,SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH ) 

    CALL DG_DERIVS( V, VOLD, &
         NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI, U_NLOC, CVWEIGHT, N, NLX, NLY, NLZ, &
         X_NONODS, X, Y, Z, &
         DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
         NFACE, FACE_ELE, U_SLOCLIST, STOTEL, U_SNLOC, WIC_U_BC, SUF_V_BC, &
         WIC_U_BC_DIRICHLET, SBCVNGI,SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH ) 

    CALL DG_DERIVS( W, WOLD, &
         NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI, U_NLOC, CVWEIGHT, N, NLX, NLY, NLZ, &
         X_NONODS, X, Y, Z, &
         DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
         NFACE, FACE_ELE, U_SLOCLIST, STOTEL, U_SNLOC, WIC_U_BC, SUF_W_BC, &
         WIC_U_BC_DIRICHLET, SBCVNGI,SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH ) 

    RETURN
  END SUBROUTINE DG_DERIVS_UVW




  SUBROUTINE DG_DERIVS( FEMT, FEMTOLD, &
       NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
       CV_NGI, CV_NLOC, CVWEIGHT, N, NLX, NLY, NLZ, &
       X_NONODS, X, Y, Z, &
       DTX_ELE, DTY_ELE, DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
       NFACE, FACE_ELE, CV_SLOCLIST, STOTEL, CV_SNLOC, WIC_T_BC, SUF_T_BC, &
       WIC_T_BC_DIRICHLET, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH ) 

    ! determine FEMT (finite element wise) etc from T (control volume wise)
    use shape_functions 
    use matrix_operations
    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, NPHASE, CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
         X_NONODS, STOTEL, CV_SNLOC, WIC_T_BC_DIRICHLET, SBCVNGI, NFACE
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: FEMT, FEMTOLD
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( STOTEL*NPHASE ), intent( in ) ::  WIC_T_BC
    REAL, DIMENSION( STOTEL*CV_SNLOC*NPHASE ), intent( in ) ::  SUF_T_BC
    INTEGER, DIMENSION( NFACE,CV_SNLOC ), intent( in ) ::  CV_SLOCLIST
    INTEGER, DIMENSION( NFACE,TOTELE ), intent( in ) ::  FACE_ELE
    REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in ) :: N, NLX, NLY, NLZ 
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( CV_NLOC, NPHASE, TOTELE ), intent( inout ) :: DTX_ELE, DTY_ELE, DTZ_ELE, &
         DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE
    REAL, DIMENSION( CV_SNLOC, SBCVNGI ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
    REAL, DIMENSION( SBCVNGI ), intent( in ) :: SBWEIGH
    ! Local variables
    LOGICAL :: D1, D3, DCYL, APPLYBC
    REAL, DIMENSION( : ), allocatable :: DETWEI, RA
    REAL, DIMENSION( :, : ), allocatable :: NX, NY, NZ
    REAL, DIMENSION( :, :, : ), allocatable :: MASELE
    REAL, DIMENSION( :, :, : ), allocatable :: VTX_ELE, VTY_ELE, VTZ_ELE, VTOLDX_ELE, VTOLDY_ELE, VTOLDZ_ELE
    REAL, DIMENSION( :, : ), allocatable :: MASS, INV_MASS
    REAL, DIMENSION( : ), allocatable :: VTX, VTY, VTZ, VTOLDX, VTOLDY, VTOLDZ, DTX, DTY, DTZ, DTOLDX, DTOLDY, DTOLDZ
    REAL, DIMENSION( : ), allocatable :: XSL, YSL, ZSL, SNORMXN, SNORMYN, SNORMZN, SDETWE
    INTEGER, DIMENSION( : ), allocatable :: SLOC2LOC, ILOC_OTHER_SIDE
    REAL :: VOLUME, NN, NNX, NNY, NNZ, NORMX, NORMY, NORMZ, SAREA, NRBC, RNN, RTBC, &
         VLM_NORX, VLM_NORY, VLM_NORZ
    INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_GI, CV_ILOC2, &
         CV_INOD, CV_INOD2, CV_JLOC2, CV_NODJ2, CV_NODJ2_IPHA, CV_NODJ_IPHA, &
         CV_SILOC, CV_SJLOC, ELE2, IFACE, IPHASE, SELE2, SUF_CV_SJ2, SUF_CV_SJ2_IPHA, &
         X_INOD, SGI

    ewrite(3,*)'in DG_DERIVS sbrt'
    ALLOCATE( DETWEI( CV_NGI )) 
    ALLOCATE( RA( CV_NGI ))
    ALLOCATE( NX( CV_NLOC, CV_NGI ))
    ALLOCATE( NY( CV_NLOC, CV_NGI ))
    ALLOCATE( NZ( CV_NLOC, CV_NGI ))
    ALLOCATE( MASELE( CV_NLOC, CV_NLOC, TOTELE ))
    ALLOCATE( VTX_ELE( CV_NLOC, NPHASE, TOTELE ))
    ALLOCATE( VTY_ELE( CV_NLOC, NPHASE, TOTELE ))
    ALLOCATE( VTZ_ELE( CV_NLOC, NPHASE, TOTELE ))
    ALLOCATE( VTOLDX_ELE( CV_NLOC, NPHASE, TOTELE ))
    ALLOCATE( VTOLDY_ELE( CV_NLOC, NPHASE, TOTELE ))
    ALLOCATE( VTOLDZ_ELE( CV_NLOC, NPHASE, TOTELE ))
    ALLOCATE( MASS( CV_NLOC, CV_NLOC) )
    ALLOCATE( INV_MASS( CV_NLOC, CV_NLOC) )
    ALLOCATE( VTX( CV_NLOC ))
    ALLOCATE( VTY( CV_NLOC ))
    ALLOCATE( VTZ( CV_NLOC ))
    ALLOCATE( VTOLDX( CV_NLOC ))
    ALLOCATE( VTOLDY( CV_NLOC ))
    ALLOCATE( VTOLDZ( CV_NLOC ))
    ALLOCATE( SLOC2LOC( CV_SNLOC ))
    ALLOCATE( ILOC_OTHER_SIDE( CV_SNLOC ))
    ALLOCATE( XSL( CV_SNLOC ))
    ALLOCATE( YSL( CV_SNLOC ))
    ALLOCATE( ZSL( CV_SNLOC ))
    ALLOCATE( SNORMXN( SBCVNGI ))
    ALLOCATE( SNORMYN( SBCVNGI ))
    ALLOCATE( SNORMZN( SBCVNGI ))
    ALLOCATE( SDETWE( SBCVNGI ))
    ALLOCATE( DTX( CV_NLOC ))
    ALLOCATE( DTY( CV_NLOC ))
    ALLOCATE( DTZ( CV_NLOC ))
    ALLOCATE( DTOLDX( CV_NLOC ))
    ALLOCATE( DTOLDY( CV_NLOC ))
    ALLOCATE( DTOLDZ( CV_NLOC ))

    MASELE = 0.0
    VTX_ELE = 0.0
    VTY_ELE = 0.0
    VTZ_ELE = 0.0
    VTOLDX_ELE = 0.0
    VTOLDY_ELE = 0.0
    VTOLDZ_ELE = 0.0

    D1 = ( NDIM == 1 )
    D3 = ( NDIM == 3 )
    DCYL = .FALSE. 

    Loop_Elements1: DO ELE = 1, TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNLXR( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, CV_NLOC, CV_NGI, &
            N, NLX, NLY, NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
            NX, NY, NZ ) 

       !       ewrite(3,*)'N',N
       !       ewrite(3,*)'nlx:',nlx
       !       ewrite(3,*)'nx:',nx
       !       ewrite(3,*)'CVWEIGHT:',CVWEIGHT
       !       ewrite(3,*)'DETWEI:',DETWEI
       !       ewrite(3,*)'volume=',volume
       !       stop 12


       Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC

          CV_NODI = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )

          Loop_CV_JLOC: DO CV_JLOC = 1, CV_NLOC

             CV_NODJ = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )

             NN = 0.0
             NNX = 0.0
             NNY = 0.0
             NNZ = 0.0
             DO CV_GI = 1, CV_NGI
                NN = NN +  N( CV_ILOC, CV_GI )  * N(  CV_JLOC, CV_GI ) * DETWEI( CV_GI )
                NNX = NNX + N( CV_ILOC, CV_GI )  * NX( CV_JLOC, CV_GI ) * DETWEI( CV_GI )
                NNY = NNY + N( CV_ILOC, CV_GI )  * NY( CV_JLOC, CV_GI ) * DETWEI( CV_GI )
                NNZ = NNZ + N( CV_ILOC, CV_GI )  * NZ( CV_JLOC, CV_GI ) * DETWEI( CV_GI )
             END DO

             MASELE( CV_ILOC,CV_JLOC, ELE)  = MASELE( CV_ILOC,CV_JLOC, ELE ) + NN

             DO IPHASE = 1, NPHASE
                CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS
                VTX_ELE( CV_ILOC, IPHASE, ELE ) = VTX_ELE( CV_ILOC, IPHASE, ELE ) + NNX * FEMT( CV_NODJ_IPHA )
                VTY_ELE( CV_ILOC, IPHASE, ELE ) = VTY_ELE( CV_ILOC, IPHASE, ELE ) + NNY * FEMT( CV_NODJ_IPHA )
                VTZ_ELE( CV_ILOC, IPHASE, ELE ) = VTZ_ELE( CV_ILOC, IPHASE, ELE ) + NNZ * FEMT( CV_NODJ_IPHA )

                VTOLDX_ELE( CV_ILOC, IPHASE, ELE ) = VTOLDX_ELE( CV_ILOC, IPHASE, ELE ) + NNX * FEMTOLD( CV_NODJ_IPHA )
                VTOLDY_ELE( CV_ILOC, IPHASE, ELE ) = VTOLDY_ELE( CV_ILOC, IPHASE, ELE ) + NNY * FEMTOLD( CV_NODJ_IPHA )
                VTOLDZ_ELE( CV_ILOC, IPHASE, ELE ) = VTOLDZ_ELE( CV_ILOC, IPHASE, ELE ) + NNZ * FEMTOLD( CV_NODJ_IPHA )
             END DO

          END DO Loop_CV_JLOC

       END DO Loop_CV_ILOC



    END DO Loop_Elements1


    ! Example of    CV_SLOCLIST for a tet element.
    !         INTEGER CV_SLOCLIST(NFACE,CV_SNLOC)
    !         ! The local cords are in anti-clockwise order
    !         IFACE=1
    !         CV_SLOCLIST(IFACE,1)=1
    !         CV_SLOCLIST(IFACE,2)=2
    !         CV_SLOCLIST(IFACE,3)=3
    !         IFACE=2
    !         CV_SLOCLIST(IFACE,1)=1
    !         CV_SLOCLIST(IFACE,2)=2
    !         CV_SLOCLIST(IFACE,3)=4
    !         IFACE=3
    !         CV_SLOCLIST(IFACE,1)=3
    !         CV_SLOCLIST(IFACE,2)=2
    !         CV_SLOCLIST(IFACE,3)=4
    !         IFACE=4
    !         CV_SLOCLIST(IFACE,1)=1
    !         CV_SLOCLIST(IFACE,2)=3
    !         CV_SLOCLIST(IFACE,3)=4

    ! Loop over surface elements
    ! EWRITE(3,*)'VTX_ELE(1,1,1 ):',VTX_ELE(1,1,1)   


    Loop_Elements2: DO ELE=1,TOTELE

       Between_Elements_And_Boundary: DO IFACE = 1, NFACE
          ELE2 = FACE_ELE( IFACE, ELE )
          SELE2 = MAX( 0, - ELE2 )
          ELE2 = MAX( 0, + ELE2 )
          !ewrite(3,*)'FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE ):',FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE )

          ! The surface nodes on element face IFACE.  
          SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )

          ! Form approximate surface normal (NORMX,NORMY,NORMZ)
          CALL DGSIMPLNORM( ELE, SLOC2LOC, TOTELE, CV_NLOC, CV_SNLOC, X_NDGLN, &
               X, Y, Z, X_NONODS, NORMX, NORMY, NORMZ )

          ! Recalculate the normal...
          DO CV_SILOC = 1, CV_SNLOC
             CV_ILOC = SLOC2LOC( CV_SILOC )
             X_INOD = X_NDGLN(( ELE - 1 ) * X_NLOC + CV_ILOC )
             XSL( CV_SILOC ) = X( X_INOD )
             YSL( CV_SILOC ) = Y( X_INOD )
             ZSL( CV_SILOC ) = Z( X_INOD )
          END DO

          CALL DGSDETNXLOC2 (CV_SNLOC, SBCVNGI, &
               XSL, YSL, ZSL, &
               SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH, SDETWE, SAREA, &
               (NDIM==1), (NDIM==3), (NDIM==-2), &
               SNORMXN, SNORMYN, SNORMZN, &
               NORMX, NORMY, NORMZ )

          !ewrite(3,*)'*********************'
          !ewrite(3,*)'ele,ele2,sele2:',ele,ele2,sele2
          !ewrite(3,*)'iface=',iface
          IF(SELE2 == 0) THEN
             ! Calculate the nodes on the other side of the face:
             ewrite(3,*)'SLOC2LOC:',SLOC2LOC
             DO CV_SILOC = 1, CV_SNLOC
                CV_ILOC = SLOC2LOC( CV_SILOC )
                CV_INOD = X_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
                DO CV_ILOC2 = 1, CV_NLOC
                   CV_INOD2 = X_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )
                   IF( CV_INOD2 == CV_INOD ) ILOC_OTHER_SIDE( CV_SILOC ) = CV_ILOC2
                END DO
             END DO
             APPLYBC=(ELE /= ELE2).AND.(ELE2 /= 0)
             !ewrite(3,*)'ele,ele2:',ele,ele2
             !ewrite(3,*)'iface=',iface
             !ewrite(3,*)'CV_SLOCLIST:',CV_SLOCLIST
             !ewrite(3,*)'SLOC2LOC:',SLOC2LOC
             !ewrite(3,*)'ILOC_OTHER_SIDE:',ILOC_OTHER_SIDE
          ELSE
             APPLYBC = ( WIC_T_BC(SELE2) == WIC_T_BC_DIRICHLET )
          ENDIF

          IF(APPLYBC) THEN
             DO CV_SILOC=1,CV_SNLOC
                CV_ILOC   =SLOC2LOC(CV_SILOC)
                DO CV_SJLOC=1,CV_SNLOC
                   CV_JLOC =SLOC2LOC(CV_SJLOC)
                   CV_NODJ=CV_NDGLN((ELE-1)*CV_NLOC+CV_JLOC)
                   IF(SELE2 /= 0) THEN
                      CV_JLOC2=CV_JLOC
                      CV_NODJ2=CV_NODJ
                      SUF_CV_SJ2 = CV_SJLOC + CV_SNLOC * ( SELE2 - 1 ) 
                      NRBC=0.0
                   ELSE
                      CV_JLOC2=ILOC_OTHER_SIDE(CV_SJLOC)
                      !ewrite(3,*)'(ELE2-1)*CV_NLOC+CV_JLOC2,ELE2,CV_NLOC,CV_JLOC2:', &
                      !     (ELE2-1)*CV_NLOC+CV_JLOC2,ELE2,CV_NLOC,CV_JLOC2
                      !ewrite(3,*)'ILOC_OTHER_SIDE:',ILOC_OTHER_SIDE
                      !ewrite(3,*)'ELE,ELE2,SELE2,CV_JLOC2=',ELE,ELE2,SELE2,CV_JLOC2
                      CV_NODJ2=CV_NDGLN((ELE2-1)*CV_NLOC+CV_JLOC2)
                      NRBC=1.0
                   ENDIF
                   VLM_NORX=0.0
                   VLM_NORY=0.0
                   VLM_NORZ=0.0
                   ! Have a surface integral on element boundary...  
                   DO SGI=1,SBCVNGI
                      RNN=SDETWE(SGI)*SBCVFEN(CV_SILOC,SGI)*SBCVFEN(CV_SJLOC,SGI)
                      VLM_NORX=VLM_NORX+SNORMXN(SGI)*RNN
                      VLM_NORY=VLM_NORY+SNORMYN(SGI)*RNN
                      VLM_NORZ=VLM_NORZ+SNORMZN(SGI)*RNN
                   END DO
                   !         EWRITE(3,*)'IFACE,CV_SILOC,CV_SJLOC:',IFACE,CV_SILOC,CV_SJLOC
                   !         EWRITE(3,*)'VLM_NORX,VLM_NORY,VLM_NORZ:',VLM_NORX,VLM_NORY,VLM_NORZ
                   !         EWRITE(3,*)'SNORMXN:',SNORMXN
                   !         EWRITE(3,*)'SDETWE:',SDETWE
                   ! add diffusion term...
                   DO IPHASE=1,NPHASE
                      CV_NODJ_IPHA =CV_NODJ  + (IPHASE-1)*CV_NONODS
                      CV_NODJ2_IPHA=CV_NODJ2 + (IPHASE-1)*CV_NONODS
                      VTX_ELE(CV_ILOC, IPHASE,ELE ) = VTX_ELE(CV_ILOC, IPHASE,ELE ) &
                           - VLM_NORX*0.5*(FEMT(CV_NODJ_IPHA)-FEMT(CV_NODJ2_IPHA)*NRBC)
                      VTY_ELE(CV_ILOC, IPHASE,ELE ) = VTY_ELE(CV_ILOC, IPHASE,ELE ) &
                           - VLM_NORY*0.5*(FEMT(CV_NODJ_IPHA)-FEMT(CV_NODJ2_IPHA)*NRBC)
                      VTZ_ELE(CV_ILOC, IPHASE,ELE ) = VTZ_ELE(CV_ILOC, IPHASE,ELE ) &
                           - VLM_NORZ*0.5*(FEMT(CV_NODJ_IPHA)-FEMT(CV_NODJ2_IPHA)*NRBC)

                      VTOLDX_ELE(CV_ILOC, IPHASE,ELE ) = VTOLDX_ELE(CV_ILOC, IPHASE,ELE ) &
                           - VLM_NORX*0.5*(FEMTOLD(CV_NODJ_IPHA)-FEMTOLD(CV_NODJ2_IPHA)*NRBC)
                      VTOLDY_ELE(CV_ILOC, IPHASE,ELE ) = VTOLDY_ELE(CV_ILOC, IPHASE,ELE ) &
                           - VLM_NORY*0.5*(FEMTOLD(CV_NODJ_IPHA)-FEMTOLD(CV_NODJ2_IPHA)*NRBC)
                      VTOLDZ_ELE(CV_ILOC, IPHASE,ELE ) = VTOLDZ_ELE(CV_ILOC, IPHASE,ELE ) &
                           - VLM_NORZ*0.5*(FEMTOLD(CV_NODJ_IPHA)-FEMTOLD(CV_NODJ2_IPHA)*NRBC)
                      IF(SELE2.NE.0) THEN
                         IF(WIC_T_BC(SELE2+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRICHLET) THEN
                            SUF_CV_SJ2_IPHA = SUF_CV_SJ2 + STOTEL * CV_SNLOC * ( IPHASE - 1 )
                            RTBC=SUF_T_BC( SUF_CV_SJ2_IPHA )
                            VTX_ELE(CV_ILOC, IPHASE,ELE ) = VTX_ELE(CV_ILOC, IPHASE,ELE ) &
                                 + VLM_NORX*0.5*RTBC
                            VTY_ELE(CV_ILOC, IPHASE,ELE ) = VTY_ELE(CV_ILOC, IPHASE,ELE ) &
                                 + VLM_NORY*0.5*RTBC
                            VTZ_ELE(CV_ILOC, IPHASE,ELE ) = VTZ_ELE(CV_ILOC, IPHASE,ELE ) &
                                 + VLM_NORZ*0.5*RTBC

                            VTOLDX_ELE(CV_ILOC, IPHASE,ELE ) = VTOLDX_ELE(CV_ILOC, IPHASE,ELE ) &
                                 + VLM_NORX*0.5*RTBC
                            VTOLDY_ELE(CV_ILOC, IPHASE,ELE ) = VTOLDY_ELE(CV_ILOC, IPHASE,ELE ) &
                                 + VLM_NORY*0.5*RTBC
                            VTOLDZ_ELE(CV_ILOC, IPHASE,ELE ) = VTOLDZ_ELE(CV_ILOC, IPHASE,ELE ) &
                                 + VLM_NORZ*0.5*RTBC
                         ENDIF
                      ENDIF
                      !       IF((ELE==1).AND.(CV_ILOC.NE.10)) EWRITE(3,*)'VTX_ELE(CV_ILOC, IPHASE,ELE ):', &
                      !                                                VTX_ELE(CV_ILOC, IPHASE,ELE )
                      !       EWRITE(3,*)'CV_ILOC=',CV_ILOC
                      !       STOP 3331
                   END DO
                END DO
             END DO
          ENDIF


       END DO Between_Elements_And_Boundary


    END DO Loop_Elements2



    ! Solve local system for the gradients DTX_ELE etc:     

    Loop_Elements3: DO ELE=1,TOTELE

       MASS(:,:)=MASELE(:,:,ELE)
       CALL MATDMATINV( MASS, INV_MASS, CV_NLOC)

       Loop_IPHASE: DO IPHASE=1,NPHASE

          VTX(:)=VTX_ELE(:, IPHASE,ELE )
          VTY(:)=VTY_ELE(:, IPHASE,ELE )
          VTZ(:)=VTZ_ELE(:, IPHASE,ELE )
          VTOLDX(:)=VTOLDX_ELE(:, IPHASE,ELE )
          VTOLDY(:)=VTOLDY_ELE(:, IPHASE,ELE )
          VTOLDZ(:)=VTOLDZ_ELE(:, IPHASE,ELE )
          DTX=0.0
          DTY=0.0
          DTZ=0.0
          DTOLDX=0.0
          DTOLDY=0.0
          DTOLDZ=0.0
          DO CV_ILOC=1,CV_NLOC
             DO CV_JLOC=1,CV_NLOC
                DTX(CV_ILOC)=DTX(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*VTX(CV_JLOC)
                DTY(CV_ILOC)=DTY(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*VTY(CV_JLOC)
                DTZ(CV_ILOC)=DTZ(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*VTZ(CV_JLOC)
                DTOLDX(CV_ILOC)=DTOLDX(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*VTOLDX(CV_JLOC)
                DTOLDY(CV_ILOC)=DTOLDY(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*VTOLDY(CV_JLOC)
                DTOLDZ(CV_ILOC)=DTOLDZ(CV_ILOC) +INV_MASS(CV_ILOC,CV_JLOC)*VTOLDZ(CV_JLOC)
             END DO
          END DO
          DTX_ELE(:, IPHASE,ELE )=DTX(:)
          DTY_ELE(:, IPHASE,ELE )=DTY(:)
          DTZ_ELE(:, IPHASE,ELE )=DTZ(:)
          DTOLDX_ELE(:, IPHASE,ELE )=DTOLDX(:)
          DTOLDY_ELE(:, IPHASE,ELE )=DTOLDY(:)
          DTOLDZ_ELE(:, IPHASE,ELE )=DTOLDZ(:)

       END DO Loop_IPHASE

    END DO Loop_Elements3

    ewrite(3,*)'--derivative:'
    do ele=1,-totele
       do cv_iloc=1,cv_nloc
          !        ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),femt(cv_ndgln((ele-1)*cv_nloc+cv_iloc))
          ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),DTX_ELE(cv_iloc,1,ele)
          !        ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),VTX_ELE(cv_iloc,1,ele)
       end do
    end do
    !    stop 331

    DEALLOCATE( DETWEI ) 
    DEALLOCATE( RA )
    DEALLOCATE( NX )
    DEALLOCATE( NY )
    DEALLOCATE( NZ )
    DEALLOCATE( MASELE )
    DEALLOCATE( VTX_ELE )
    DEALLOCATE( VTY_ELE )
    DEALLOCATE( VTZ_ELE )
    DEALLOCATE( VTOLDX_ELE )
    DEALLOCATE( VTOLDY_ELE )
    DEALLOCATE( VTOLDZ_ELE )
    DEALLOCATE( MASS )
    DEALLOCATE( INV_MASS )
    DEALLOCATE( VTX )
    DEALLOCATE( VTY )
    DEALLOCATE( VTZ )
    DEALLOCATE( VTOLDX )
    DEALLOCATE( VTOLDY )
    DEALLOCATE( VTOLDZ )
    DEALLOCATE( SLOC2LOC )
    DEALLOCATE( ILOC_OTHER_SIDE )
    DEALLOCATE( XSL )
    DEALLOCATE( YSL )
    DEALLOCATE( ZSL )
    DEALLOCATE( SNORMXN )
    DEALLOCATE( SNORMYN )
    DEALLOCATE( SNORMZN )
    DEALLOCATE( SDETWE )
    DEALLOCATE( DTX )
    DEALLOCATE( DTY )
    DEALLOCATE( DTZ )
    DEALLOCATE( DTOLDX )
    DEALLOCATE( DTOLDY )
    DEALLOCATE( DTOLDZ )

    RETURN

  END SUBROUTINE DG_DERIVS




  SUBROUTINE ONVDLIM( TOTELE, &
       TDLIM, TDCEN, INCOME, PELE, PELEOT, &
       ETDNEW, TDMIN, TDMAX, FIRORD, NOLIMI )
    implicit none 
    ! This sub calculates the limited face values TDADJ(1...SNGI) from the central 
    ! difference face values TDCEN(1...SNGI) using a NVD shceme.  
    ! INCOME(1...SNGI)=1 for incomming to element ELE  else =0.
    ! LIBETA is the flux limiting parameter. 
    ! TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
    ! TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
    ! PELEOT=element at other side of current face. 
    ! ELEOT2=element at other side of the element ELEOTH. 
    ! ELESID=element next to oposing current face. 
    ! The elements are arranged in this order: ELEOT2,ELE, PELEOT, ELESID. 
    ! This sub finds the neighbouring elements. Suppose that this is the face IFACE. 
    !---------------------------------------------------
    !|   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |  
    !--------------------------------------------------- 
    ! TAIN         THALF       TAOUT
    !--------------------------------------------------- 
    !>TEXTIN
    !TEXOUT<
    !--------------------------------------------------- 
    INTEGER, intent( in ) :: TOTELE 
    REAL, intent( inout ) :: TDLIM  
    REAL, intent( in ) :: TDCEN, INCOME
    INTEGER, intent( in ) :: PELE, PELEOT
    REAL, DIMENSION( TOTELE ), intent( in ) :: ETDNEW, TDMIN, TDMAX
    LOGICAL, intent( in ) :: FIRORD, NOLIMI
    ! Local variables   
    REAL, PARAMETER :: TOLER=1.0E-10
    REAL :: UCIN, UCOU, TUPWIN, TUPWI2, TDELE, DENOIN, CTILIN, DENOOU, &
         CTILOU, FTILIN, FTILOU

    IF( NOLIMI ) THEN
       TDLIM = TDCEN
       RETURN
    ENDIF

    Conditional_PELEOT: IF( PELEOT /= PELE ) THEN

       IF( ETDNEW( PELEOT ) > ETDNEW( PELE )) THEN
          TUPWIN = TDMAX( PELEOT )
          TUPWI2 = TDMIN( PELE )
       ELSE
          TUPWIN = TDMIN( PELEOT )
          TUPWI2 = TDMAX( PELE ) 
       ENDIF

       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW( PELE ) 
       DENOIN = TDELE - TUPWIN 

       IF( ABS( DENOIN ) < TOLER ) DENOIN = SIGN( TOLER, DENOIN )

       UCIN = ETDNEW( PELEOT )
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW( PELEOT ) 
       DENOOU = TDELE - TUPWI2

       IF( ABS( DENOOU ) < TOLER) DENOOU = SIGN( TOLER, DENOOU )
       UCOU = ETDNEW( PELE )
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU

    ELSE

       ! Calculate normalisation parameters for incomming velocities 
       TUPWIN = ETDNEW( PELE )
       UCIN = ETDNEW( PELE )
       DENOIN = 1.
       CTILIN = 0.

       ! Calculate normalisation parameters for out going velocities 
       TUPWI2 = ETDNEW( PELE )
       UCOU = ETDNEW( PELE )
       DENOOU = 1.
       CTILOU = 0.

    ENDIF Conditional_PELEOT

    Conditional_FIRORD: IF( FIRORD ) THEN ! Velocity is pointing into element

       ! Velocity is going out of element
       TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 

    ELSE

       FTILIN = ( TDCEN - TUPWIN ) / DENOIN
       FTILOU = ( TDCEN - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM= INCOME*( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, -1.0 ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, -1.0 ) &
            * DENOOU )

    ENDIF Conditional_FIRORD

    RETURN

  END SUBROUTINE ONVDLIM





  SUBROUTINE ONVDLIMsqrt( TOTELE, &
       TDLIM, TDCEN, INCOME, PELE, PELEOT, &
       ETDNEW, TDMIN, TDMAX, FIRORD, NOLIMI )
    implicit none 
    ! This sub calculates the limited face values TDADJ(1...SNGI) from the central 
    ! difference face values TDCEN(1...SNGI) using a NVD shceme.  
    ! INCOME(1...SNGI)=1 for incomming to element ELE  else =0.
    ! LIBETA is the flux limiting parameter. 
    ! TDMAX(PELE)=maximum of the surrounding 6 element values of element PELE.
    ! TDMIN(PELE)=minimum of the surrounding 6 element values of element PELE.
    ! PELEOT=element at other side of current face. 
    ! ELEOT2=element at other side of the element ELEOTH. 
    ! ELESID=element next to oposing current face. 
    ! The elements are arranged in this order: ELEOT2,ELE, PELEOT, ELESID. 
    ! This sub finds the neighbouring elements. Suppose that this is the face IFACE. 
    !---------------------------------------------------
    !|   ELEOT2   |   ELEOTH   |   ELE     |   ELESID   |  
    !--------------------------------------------------- 
    ! TAIN         THALF       TAOUT
    !--------------------------------------------------- 
    !>TEXTIN
    !TEXOUT<
    !--------------------------------------------------- 
    INTEGER, intent( in ) :: TOTELE 
    REAL, intent( inout ) :: TDLIM  
    REAL, intent( in ) :: TDCEN, INCOME
    INTEGER, intent( in ) :: PELE, PELEOT
    REAL, DIMENSION( TOTELE ), intent( in ) :: ETDNEW, TDMIN, TDMAX
    LOGICAL, intent( in ) :: FIRORD, NOLIMI
    ! Local variables   
    REAL, PARAMETER :: TOLER=1.0E-10
    REAL :: UCIN, UCOU, TUPWIN, TUPWI2, TDELE, DENOIN, CTILIN, DENOOU, &
         CTILOU, FTILIN, FTILOU

    IF( NOLIMI ) THEN
       TDLIM = TDCEN
       RETURN
    ENDIF

    Conditional_PELEOT: IF( PELEOT /= PELE ) THEN

       IF( ETDNEW( PELEOT ) > ETDNEW( PELE )) THEN
          TUPWIN = TDMAX( PELEOT )**2
          TUPWI2 = TDMIN( PELE )**2
       ELSE
          TUPWIN = TDMIN( PELEOT )**2
          TUPWI2 = TDMAX( PELE ) **2
       ENDIF

       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW( PELE )**2 
       DENOIN = TDELE - TUPWIN

       IF( ABS( DENOIN ) < TOLER ) DENOIN = SIGN( TOLER, DENOIN )

       UCIN = ETDNEW( PELEOT )**2
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW( PELEOT )**2 
       DENOOU = TDELE - TUPWI2

       IF( ABS( DENOOU ) < TOLER) DENOOU = SIGN( TOLER, DENOOU )
       UCOU = ETDNEW( PELE )**2
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU

    ELSE

       ! Calculate normalisation parameters for incomming velocities 
       TUPWIN = ETDNEW( PELE )**2
       UCIN = ETDNEW( PELE )**2
       DENOIN = 1.
       CTILIN = 0.

       ! Calculate normalisation parameters for out going velocities 
       TUPWI2 = ETDNEW( PELE )**2
       UCOU = ETDNEW( PELE )**2
       DENOOU = 1.
       CTILOU = 0.

    ENDIF Conditional_PELEOT

    Conditional_FIRORD: IF( FIRORD ) THEN ! Velocity is pointing into element

       ! Velocity is going out of element
       TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 
       TDLIM=sqrt(max(TDLIM,0.0))

    ELSE

       FTILIN = ( max(0.0,TDCEN)**2 - TUPWIN ) / DENOIN
       FTILOU = ( max(0.0,TDCEN)**2 - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM= INCOME*( TUPWIN + NVDFUNNEWsqrt( FTILIN, CTILIN, -1.0 ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEWsqrt( FTILOU, CTILOU, -1.0 ) &
            * DENOOU )

       TDLIM=sqrt(max(TDLIM,0.0))

    ENDIF Conditional_FIRORD

    RETURN

  END SUBROUTINE ONVDLIMsqrt




  REAL FUNCTION NVDFUNNEWsqrt( UF, UC, COURAT )
    implicit none
    ! The function computes NVDFUNNEW, the normalised value of the 
    ! advected variable on the face of the control volume, based on 
    ! the normalised value of the advected variable in the donor CV,
    ! UC, and the high-order estimate of the face value UF. 
    ! NVDFUNNEW is limited so that it is in the non-oscillatory 
    ! region of normalised variable diagram (NVD).
    !
    ! This version can also use an extremely compressive advection 
    ! scheme, appropriate for advecting volume/mass fractions in
    ! multi-material problems.  This scheme is applied if the courant
    ! number (COURAT) is positive.  IF YOU DO NOT WANT COMPRESSIVE
    ! ADVECTION: call this function with COURAT=-1.
    !
    ! # INPUTS:
    ! UC     - Normalised advected field value at donor CV (node)
    ! UF     - High-order estimate of normalised value at CV face
    ! COURAT - Local courant number (flow_vel*dt/dx)
    !      (set to -ve number if not using interface tracking scheme)
    !
    ! # OUTPUTS:
    ! NVDFUNNEW - Limited normalised advected field value at CV face
    !
    ! The function is based on the flux limiting function in the paper: Riemann 
    ! solvers on 3-D unstructured meshes for radiation transport, written by 
    ! Dr C.C. Pain et al. See also the paper: The ULTIMATE conservative 
    ! difference scheme applied to unsteady 1D advection, 
    ! Leonard B. P., 1991, Computer Methods in Applied Mechanics and 
    ! Engineering, 88, 17-74.
    !
    ! Description                                   Programmer      Date
    ! ==================================================================
    ! Original subroutine .............................CCP    2002-01-21
    ! Interface tracking option added..................GSC    2006-03-14
    ! ==================================================================
    !
    ! XI is the parameter in equation 38 of the Riemann paper. If XI is equal
    ! to 2 then this corresponds to a TVD condition in 1-D, a value of XI 
    ! equal to 3 has been recommended elsewhere
    ! 
    REAL :: UC, UF, COURAT
    ! Local variables
    REAL, PARAMETER :: XI = 2. ! XI = 1.
    !    REAL, PARAMETER :: XI = 10. ! XI = 1.
    !    REAL, PARAMETER :: XI = 5. ! XI = 1.
    REAL :: TILDEUF, MAXUF

    ! For the region 0 < UC < 1 on the NVD, define the limiter
    IF( ( UC > 0.0 ) .AND. ( UC < 1.0 ) ) THEN
       ! For the interface tracking scheme, use Hyper-C (NOTE: at high Courant numbers
       !  limit is defined by XI, not the Hyper-C scheme.)

       IF( COURAT > 0.0 ) THEN
          TILDEUF = MIN( 1.0, max( UC / COURAT, XI * UC ))
       ELSE !For the normal limiting
          MAXUF = MAX( 0.0, UF )
          !          MAXUF = MAX( UC, UF )
          TILDEUF = MIN( 1.0, XI * UC, MAXUF )
          !          TILDEUF = MIN( 1.0, 1.5* UC, MAXUF )

          !          TILDEUF = MIN( 1.0, SQRT(UC), MAXUF )
       ENDIF

    ELSE ! Outside the region 0<UC<1 on the NVD, use first-order upwinding
       TILDEUF = UC
    ENDIF

    NVDFUNNEWsqrt = TILDEUF

  end function nvdfunnewsqrt




  REAL FUNCTION NVDFUNNEW( UF, UC, COURAT )
    implicit none
    ! The function computes NVDFUNNEW, the normalised value of the 
    ! advected variable on the face of the control volume, based on 
    ! the normalised value of the advected variable in the donor CV,
    ! UC, and the high-order estimate of the face value UF. 
    ! NVDFUNNEW is limited so that it is in the non-oscillatory 
    ! region of normalised variable diagram (NVD).
    !
    ! This version can also use an extremely compressive advection 
    ! scheme, appropriate for advecting volume/mass fractions in
    ! multi-material problems.  This scheme is applied if the courant
    ! number (COURAT) is positive.  IF YOU DO NOT WANT COMPRESSIVE
    ! ADVECTION: call this function with COURAT=-1.
    !
    ! # INPUTS:
    ! UC     - Normalised advected field value at donor CV (node)
    ! UF     - High-order estimate of normalised value at CV face
    ! COURAT - Local courant number (flow_vel*dt/dx)
    !      (set to -ve number if not using interface tracking scheme)
    !
    ! # OUTPUTS:
    ! NVDFUNNEW - Limited normalised advected field value at CV face
    !
    ! The function is based on the flux limiting function in the paper: Riemann 
    ! solvers on 3-D unstructured meshes for radiation transport, written by 
    ! Dr C.C. Pain et al. See also the paper: The ULTIMATE conservative 
    ! difference scheme applied to unsteady 1D advection, 
    ! Leonard B. P., 1991, Computer Methods in Applied Mechanics and 
    ! Engineering, 88, 17-74.
    !
    ! Description                                   Programmer      Date
    ! ==================================================================
    ! Original subroutine .............................CCP    2002-01-21
    ! Interface tracking option added..................GSC    2006-03-14
    ! ==================================================================
    !
    ! XI is the parameter in equation 38 of the Riemann paper. If XI is equal
    ! to 2 then this corresponds to a TVD condition in 1-D, a value of XI 
    ! equal to 3 has been recommended elsewhere
    ! 
    REAL :: UC, UF, COURAT
    ! Local variables
    REAL, PARAMETER :: XI = 2. ! XI = 1.
    REAL :: TILDEUF, MAXUF

    ! For the region 0 < UC < 1 on the NVD, define the limiter
    IF( ( UC > 0.0 ) .AND. ( UC < 1.0 ) ) THEN
       ! For the interface tracking scheme, use Hyper-C (NOTE: at high Courant numbers
       !  limit is defined by XI, not the Hyper-C scheme.)

       IF( COURAT > 0.0 ) THEN
          TILDEUF = MIN( 1.0, max( UC / COURAT, XI * UC ))
       ELSE !For the normal limiting
          MAXUF = MAX( 0.0, UF )
          TILDEUF = MIN( 1.0, XI * UC, MAXUF )
       ENDIF

    ELSE ! Outside the region 0<UC<1 on the NVD, use first-order upwinding
       TILDEUF = UC
    ENDIF

    NVDFUNNEW = TILDEUF

  end function nvdfunnew


  SUBROUTINE CALC_MASS_CV( MASS_CV, CV_NONODS, X_NONODS, &
       X, X_NDGLN, TOTELE, CV_NLOC )
    implicit none
    INTEGER, intent( in ) :: CV_NONODS, X_NONODS, TOTELE, CV_NLOC
    REAL, DIMENSION( CV_NONODS ), intent( inout ) ::  MASS_CV
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: X_NDGLN
    ! Local variables
    INTEGER ::ELE, NOD1, NOD2, NOD3
    REAL :: DX

    DO ELE = 1, TOTELE

       NOD1 = X_NDGLN(( ELE - 1 ) * CV_NLOC + 1 )
       NOD2 = X_NDGLN(( ELE - 1 ) * CV_NLOC + 2 )
       NOD3 = X_NDGLN(( ELE - 1 ) * CV_NLOC + 3 )

       DX = X( NOD3 ) - X( NOD1 )

       MASS_CV( NOD1 ) = 0.25 * DX
       MASS_CV( NOD2 ) = 0.50 * DX
       MASS_CV( NOD3 ) = 0.25 * DX

    END DO

    RETURN

  END SUBROUTINE CALC_MASS_CV




  SUBROUTINE SCVDETNX( ELE,      GI,        &
                                !     - INTEGERS
       NLOC,     SVNGI,   TOTELE,   &
       XNDGLN,   XNONOD,&
                                !     - REALS
       CVDETWEI, CVNORMX, CVNORMY,  &
       CVNORMZ,  SVN,     SVNLX,    &
       SVNLY,    SVWEIGH, XC,     &
       YC,       ZC,      X,        &
       Y,        Z,  &
                                !     - LOGICALS
       D1,       D3,       DCYL )
    !     --------------------------------------------------
    !     
    !     - this subroutine calculates the control volume (CV) 
    !     - CVNORMX, CVNORMY, CVNORMZ normals at the Gaussian 
    !     - integration points GI. NODI = the current global 
    !     - node number for the co-ordinates. 
    !     - (XC,YC,ZC) is the centre of CV NODI
    !     
    !     -------------------------------
    !     - date last modified : 15/03/2003
    !     -------------------------------
    INTEGER, intent( in ) :: ELE,    GI
    INTEGER, intent( in ) ::  NLOC
    INTEGER , intent( in ) ::  SVNGI,  TOTELE
    INTEGER, intent( in ) ::   XNONOD     
    REAL, intent( in ) ::   XC,YC,ZC
    INTEGER, DIMENSION( TOTELE * NLOC ), intent( in ) :: XNDGLN
    REAL, DIMENSION( SVNGI ), intent( inout ) :: CVDETWEI, CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( NLOC, SVNGI ), intent( in ) :: SVN, SVNLX, SVNLY
    REAL, DIMENSION( SVNGI ), intent( in ) :: SVWEIGH
    REAL, DIMENSION( XNONOD ), intent( in ) :: X, Y, Z
    LOGICAL, intent( in ) ::  D1, D3, DCYL

    !     - Local variables
    INTEGER :: NODJ,  JLOC
    REAL :: A, B, C
    REAL :: DETJ
    REAL :: DXDLX, DXDLY, DYDLX
    REAL :: DYDLY, DZDLX, DZDLY
    REAL :: TWOPI
    REAL, PARAMETER :: PI = 3.14159265
    REAL :: POSVGIX, POSVGIY, POSVGIZ
    REAL :: RGI

    Conditional_Dimension: IF( D3 ) THEN

       DXDLX = 0.0
       DXDLY = 0.0

       DYDLX = 0.0
       DYDLY = 0.0

       DZDLX = 0.0
       DZDLY = 0.0

       POSVGIX = 0.0
       POSVGIY = 0.0
       POSVGIZ = 0.0

       do  JLOC = 1, NLOC

          NODJ = XNDGLN((ELE-1)*NLOC+JLOC)

          DXDLX = DXDLX + SVNLX(JLOC,GI)*X(NODJ)
          DXDLY = DXDLY + SVNLY(JLOC,GI)*X(NODJ) 
          DYDLX = DYDLX + SVNLX(JLOC,GI)*Y(NODJ) 
          DYDLY = DYDLY + SVNLY(JLOC,GI)*Y(NODJ) 
          DZDLX = DZDLX + SVNLX(JLOC,GI)*Z(NODJ) 
          DZDLY = DZDLY + SVNLY(JLOC,GI)*Z(NODJ) 

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X(NODJ)
          POSVGIY = POSVGIY + SVN(JLOC,GI)*Y(NODJ)
          POSVGIZ = POSVGIZ + SVN(JLOC,GI)*Z(NODJ)
       end do

       !     - Note that POSVGIX,POSVGIY and POSVGIZ can be considered as the 
       !     - components of the Gauss pnt GI with the co-ordinate origin 
       !     - positioned at the current control volume NODI.

       POSVGIX = POSVGIX - XC
       POSVGIY = POSVGIY - YC
       POSVGIZ = POSVGIZ - ZC

       A = DYDLX*DZDLY - DYDLY*DZDLX
       B = DXDLX*DZDLY - DXDLY*DZDLX
       C = DXDLX*DYDLY - DXDLY*DYDLX
       !
       !     - Calculate the determinant of the Jacobian at Gauss pnt GI.
       !     
       DETJ = SQRT( A**2 + B**2 + C**2 )
       !     
       !     - Calculate the determinant times the surface weight at Gauss pnt GI.
       !     
       CVDETWEI(GI) = DETJ*SVWEIGH(GI)
       !     
       !     - Calculate the normal at the Gauss pts
       !     - TANX1 = DXDLX, TANY1 = DYDLX, TANZ1 = DZDLX,    
       !     - TANX2 = DXDLY, TANY2 = DYDLY, TANZ2 = DZDLY
       !     - Perform cross-product. N = T1 x T2
       !     
       CALL NORMGI( CVNORMX(GI), CVNORMY(GI), CVNORMZ(GI),&
            DXDLX,       DYDLX,       DZDLX, &
            DXDLY,       DYDLY,       DZDLY,&
            POSVGIX,     POSVGIY,     POSVGIZ ) 

    ELSE IF(.NOT.D1) THEN

       TWOPI = 1.0

       IF( DCYL ) TWOPI = 2.0*PI

       RGI   = 0.0

       DXDLX = 0.0
       DXDLY = 0.0

       DYDLX = 0.0
       DYDLY = 0.0

       DZDLX = 0.0
       !     
       !     - Note that we set the derivative wrt to y of coordinate z to 1.0
       !     
       DZDLY = 1.0

       POSVGIX = 0.0
       POSVGIY = 0.0
       POSVGIZ = 0.0

       do  JLOC = 1, NLOC! Was loop 300

          NODJ = XNDGLN((ELE-1)*NLOC+JLOC)

          DXDLX = DXDLX + SVNLX(JLOC,GI)*X(NODJ) 
          DYDLX = DYDLX + SVNLX(JLOC,GI)*Y(NODJ) 

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X(NODJ)
          POSVGIY = POSVGIY + SVN(JLOC,GI)*Y(NODJ)

          RGI = RGI + SVN(JLOC,GI)*Y(NODJ)

       end do ! Was loop 300
       !     
       !     - Note that POSVGIX and POSVGIY can be considered as the components 
       !     - of the Gauss pnt GI with the co-ordinate origin positioned at the
       !     - current control volume NODI.
       !     

       POSVGIX = POSVGIX - XC
       POSVGIY = POSVGIY - YC

       IF( .NOT. DCYL ) RGI = 1.0 

       DETJ = SQRT( DXDLX**2 + DYDLX**2 )
       CVDETWEI(GI)  = TWOPI*RGI*DETJ*SVWEIGH(GI)
       !
       !     - Calculate the normal at the Gauss pts
       !     - TANX1 = DXDLX, TANY1 = DYDLX, TANZ1 = DZDLX,    
       !     - TANX2 = DXDLY, TANY2 = DYDLY, TANZ2 = DZDLY
       !     - Perform cross-product. N = T1 x T2
       !     
       CALL NORMGI( CVNORMX(GI), CVNORMY(GI), CVNORMZ(GI),&
            DXDLX,       DYDLX,       DZDLX, &
            DXDLY,       DYDLY,       DZDLY,&
            POSVGIX,     POSVGIY,     POSVGIZ )
       !     
       !     - End of GI loop

    ELSE 
       ! For 1D...

       POSVGIX = 0.0

       do  JLOC = 1, NLOC! Was loop 300

          NODJ = XNDGLN((ELE-1)*NLOC+JLOC)

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X(NODJ)

       end do ! Was loop 300
       !     
       !     - Note that POSVGIX and POSVGIY can be considered as the components 
       !     - of the Gauss pnt GI with the co-ordinate origin positioned at the
       !     - current control volume NODI.
       !     
       !          EWRITE(3,*)'POSVGIX, XC,POSVGIX - XC:',POSVGIX, XC,POSVGIX - XC
       POSVGIX = POSVGIX - XC
       IF(POSVGIX > 0 ) THEN
          CVNORMX(GI) = +1.0
       ELSE
          CVNORMX(GI) = -1.0
       ENDIF
       CVNORMY(GI)=0.0 
       CVNORMZ(GI)=0.0

       DETJ = 1.0
       CVDETWEI(GI)  = DETJ*SVWEIGH(GI)

       !       IF(GI.EQ.3) THEN
       !         EWRITE(3,*)'CVNORMX(GI),POSVGIX,XC:',CVNORMX(GI),POSVGIX,XC
       !         STOP 39344
       !       ENDIF

    ENDIF Conditional_Dimension

  END SUBROUTINE SCVDETNX





  SUBROUTINE SIMPNORM( NORMX, NORMY, NORMZ, D3, &
       SNDGLN, STOTEL, SNLOC, X_NONODS, NONODS, ELE, &
       X, Y, Z, &
       NCOLM, FINDRM, COLM)
    ! Calculate an approx normal (normx,normy,normz)
    IMPLICIT NONE

    REAL, intent( inout ) :: NORMX, NORMY, NORMZ
    LOGICAL, intent( in ) :: D3
    INTEGER, intent( in ) :: STOTEL, SNLOC, X_NONODS, NONODS, ELE, NCOLM
    INTEGER, DIMENSION( STOTEL * SNLOC ), intent( in ) :: SNDGLN
    REAL, DIMENSION( X_NONODS), intent( in ) :: X, Y, Z
    INTEGER, DIMENSION( NONODS + 1 ), intent( in ) :: FINDRM
    INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
    ! Local variables
    REAL :: XCS, YCS, ZCS, XCV, YCV, ZCV, NORM
    INTEGER :: ICV, L, IGL, L2, IGL2, COUNT, COUNT2, NODJ, NODJ2
    LOGICAL :: FOUND, SURF, ALLSUF

    XCS = 0.
    YCS = 0.
    ZCS = 0.
    DO L=1,SNLOC
       IGL = SNDGLN(( ELE - 1 ) * SNLOC + L )
       XCS = XCS + X( IGL ) / REAL( SNLOC )
       YCS = YCS + Y( IGL ) / REAL( SNLOC )
       IF(D3) ZCS = ZCS + Z( IGL ) / REAL( SNLOC )
    END DO
    ewrite(3,*)'XCS,YCS,ZCS:',XCS,YCS,ZCS

    XCV = 0.
    YCV = 0.
    ZCV = 0.
    ICV = 0
    L = 1
    IGL = SNDGLN(( ELE - 1 ) * SNLOC + L )

    Loop_Row: DO COUNT = FINDRM( IGL ) , FINDRM( IGL + 1 ) - 1, 1
       NODJ = COLM( COUNT )
       Inside_matrix: IF(NODJ <= NONODS) THEN
          ! Make sure its not a surface node of surface element ELE
          SURF = .FALSE.
          DO L2 = 1, SNLOC
             IGL2 = SNDGLN(( ELE - 1 ) * SNLOC + L2 )
             IF( IGL2 == NODJ ) SURF = .TRUE.
          END DO

          Cond_Surf: IF( .NOT. SURF) THEN
             ! make sure NODJ is connected to all surface nodes and is not a surface node for 
             ! surface element.
             ALLSUF = .TRUE.
             DO L2 = 1, SNLOC
                IGL2 = SNDGLN(( ELE - 1 ) * SNLOC + L2 )
                FOUND = .FALSE.
                DO COUNT2 = FINDRM( NODJ ) , FINDRM( NODJ + 1 ) - 1, 1
                   NODJ2 = COLM( COUNT2 )
                   IF( IGL2 == NODJ2) FOUND = .TRUE.
                END DO
                IF( .NOT. FOUND) ALLSUF = .FALSE.
             END DO

             IF( ALLSUF ) THEN
                XCV = XCV + X( NODJ )
                YCV = YCV + Y( NODJ )
                IF( D3 ) ZCV = ZCV + Z( NODJ )
                ICV = ICV + 1
             ENDIF

          ENDIF Cond_Surf
       ENDIF Inside_matrix

    END DO Loop_Row

    XCV = XCV / REAL( ICV )
    YCV = YCV / REAL( ICV )
    IF( D3 ) ZCV = ZCV / REAL( ICV )
    NORMX = XCS - XCV
    NORMY = YCS - YCV
    NORMZ = ZCS - ZCV

    NORM = NORMX **2 + NORMY **2
    IF( D3 ) NORM = NORM + NORMZ ** 2
    NORM = MAX( 1.E-15, SQRT( NORM ))

    NORMX = NORMX / NORM
    NORMY = NORMY / NORM
    IF( D3 ) NORMZ = NORMZ / NORM

    RETURN

  END SUBROUTINE SIMPNORM





  SUBROUTINE SDETNX2_PLUS( SELE, CV_SNDGLN,&
       STOTEL, X_NONODS, CV_SNLOC, SBNGI, X, Y, Z, NDIM, &
       SBWEIGH, SBDETWEI, &
       SBUFEN, SBUFENSLX, SBUFENSLY,  &
       SAREA, SNORMXN, SNORMYN, SNORMZN, SNORMX, SNORMY, SNORMZ )
    ! use AdvectionDiffusion
    IMPLICIT NONE
    INTEGER, intent( in ) :: SELE, STOTEL, X_NONODS, CV_SNLOC, SBNGI, NDIM
    INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
    REAL, DIMENSION( CV_SNLOC, SBNGI ), intent( in ) :: SBUFEN, SBUFENSLX, SBUFENSLY
    REAL, DIMENSION( SBNGI ), intent( in ) :: SBWEIGH
    REAL, DIMENSION( SBNGI ), intent( inout ) :: SBDETWEI
    REAL, DIMENSION( SBNGI ), intent( inout ) :: SNORMXN, SNORMYN, SNORMZN
    REAL, intent( in ) :: SNORMX, SNORMY, SNORMZ
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X,Y,Z
    ! Local variables
    REAL, PARAMETER :: PIE = 3.141592654
    INTEGER :: GI, SL, IGLX 
    REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY, A, B, C, &
         DETJ, RGI, TWOPIE, SAREA, RSNORMZ
    LOGICAL :: D1, D3, DCYL

    D1 = NDIM == 1
    D3 = NDIM == 3
    DCYL = NDIM == -2

    SAREA=0.

    Conditional_Dimension: IF( D1 ) THEN
       SNORMXN = SNORMX
       SNORMYN = 0.0
       SNORMZN = 0.0
       SBDETWEI = 1.0
       SAREA = 1.0

    ELSE IF( D3 ) THEN
       Loop_Gauss3D: DO GI = 1, SBNGI
          DXDLX = 0.
          DXDLY = 0.
          DYDLX = 0.
          DYDLY = 0.
          DZDLX = 0.
          DZDLY = 0.

          DO SL = 1, CV_SNLOC
             IGLX  = CV_SNDGLN(( SELE - 1 ) * CV_SNLOC + SL )
             DXDLX = DXDLX + SBUFENSLX( SL, GI ) * X( IGLX )
             DXDLY = DXDLY + SBUFENSLY( SL, GI ) * X( IGLX ) 
             DYDLX = DYDLX + SBUFENSLX( SL, GI ) * Y( IGLX ) 
             DYDLY = DYDLY + SBUFENSLY( SL, GI ) * Y( IGLX ) 
             DZDLX = DZDLX + SBUFENSLX( SL, GI ) * Z( IGLX ) 
             DZDLY = DZDLY + SBUFENSLY( SL, GI ) * Z( IGLX ) 
          END DO

          A = DYDLX * DZDLY - DYDLY * DZDLX
          B = DXDLX * DZDLY - DXDLY * DZDLX
          C = DXDLX * DYDLY - DXDLY * DYDLX

          DETJ = SQRT( A**2 + B**2 + C**2 )
          SBDETWEI( GI ) = DETJ * SBWEIGH( GI )
          SAREA = SAREA + SBDETWEI( GI )

          ! Calculate the normal at the Gauss pts:
          !  TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
          ! Perform x-product. N=T1 x T2
          CALL NORMGI( SNORMXN(GI), SNORMYN(GI), SNORMZN(GI), &
               DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
               SNORMX, SNORMY, SNORMZ)  

       END DO Loop_Gauss3D

    ELSE

       TWOPIE = 1.0 
       IF( DCYL ) TWOPIE = 2. * PIE
       !
       DO GI = 1, SBNGI
          RGI = 0.
          DXDLX = 0.
          DXDLY = 0.
          DYDLX = 0.
          DYDLY = 0.
          DZDLX = 0.
          DZDLY=1. ! DZDLY=1 is to calculate the normal. 

          DO SL = 1, CV_SNLOC

             IGLX =  CV_SNDGLN(( SELE - 1 ) * CV_SNLOC + SL )
             DXDLX = DXDLX + SBUFENSLX( SL, GI ) * X( IGLX ) 
             DYDLX = DYDLX + SBUFENSLX( SL, GI ) * Y( IGLX )

             RGI =   RGI + SBUFEN( SL, GI ) * Y( IGLX )

          END DO

          IF( .NOT. DCYL) RGI = 1.0 
          DETJ = SQRT( DXDLX**2 + DYDLX**2 )
          SBDETWEI( GI ) = TWOPIE * RGI * DETJ * SBWEIGH( GI )
          SAREA = SAREA + SBDETWEI( GI ) 

          ! Calculate the normal at the Gauss pts:
          ! TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
          ! Perform x-product. N=T1 x T2

          RSNORMZ = 0.
          CALL NORMGI( SNORMXN( GI ), SNORMYN( GI ), SNORMZN( GI ), &
               DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
               SNORMX, SNORMY, RSNORMZ ) 

       END DO

    ENDIF Conditional_Dimension

    RETURN

  END SUBROUTINE SDETNX2_PLUS





  SUBROUTINE SDETNX2( ELE, SNDGLN, &
       STOTEL, X_NONODS, SNLOC, SNGI, & 
       X, Y, Z, &
       SN, SNLX, SNLY, SWEIGH, DETWEI, SAREA, D1, D3, DCYL,  &
       NORMXN, NORMYN, NORMZN, &
       NORMX, NORMY, NORMZ ) 
    ! use AdvectionDiffusion
    IMPLICIT NONE
    INTEGER, intent( in ) :: ELE, STOTEL, X_NONODS, SNLOC, SNGI
    INTEGER, DIMENSION( STOTEL * SNLOC ), intent( in ) :: SNDGLN
    REAL, DIMENSION( SNLOC, SNGI ), intent( in ) :: SN, SNLX, SNLY
    REAL, DIMENSION( SNGI ), intent( in ) :: SWEIGH
    REAL, DIMENSION( SNGI ), intent( inout ) :: DETWEI
    REAL, intent( inout ) :: SAREA
    LOGICAL, intent( in ) :: D1, D3, DCYL
    REAL, DIMENSION( SNGI ), intent( inout ) :: NORMXN, NORMYN, NORMZN
    REAL, intent( inout ) :: NORMX, NORMY, NORMZ
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X,Y,Z
    ! Local variables
    REAL, PARAMETER :: PIE = 3.141592654
    INTEGER :: GI, SL, IGLX
    REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY, A, B, C, &
         DETJ,RGI,TWOPIE
    !
    SAREA=0.
    !

    Conditional_Dimension: IF( D1 ) THEN
       NORMXN = NORMX
       NORMYN = 0.0
       NORMZN = 0.0
       DETWEI = 1.0
       SAREA = 1.0

    ELSE IF( D3 ) THEN
       Loop_Gauss3D: DO GI = 1, SNGI
          DXDLX = 0.
          DXDLY = 0.
          DYDLX = 0.
          DYDLY = 0.
          DZDLX = 0.
          DZDLY = 0.

          DO SL = 1, SNLOC
             IGLX  = SNDGLN(( ELE - 1 ) * SNLOC + SL )
             DXDLX = DXDLX + SNLX( SL, GI ) * X( IGLX )
             DXDLY = DXDLY + SNLY( SL, GI ) * X( IGLX ) 
             DYDLX = DYDLX + SNLX( SL, GI ) * Y( IGLX ) 
             DYDLY = DYDLY + SNLY( SL, GI ) * Y( IGLX ) 
             DZDLX = DZDLX + SNLX( SL, GI ) * Z( IGLX ) 
             DZDLY = DZDLY + SNLY( SL, GI ) * Z( IGLX ) 
          END DO

          A = DYDLX * DZDLY - DYDLY * DZDLX
          B = DXDLX * DZDLY - DXDLY * DZDLX
          C = DXDLX * DYDLY - DXDLY * DYDLX

          DETJ = SQRT( A**2 + B**2 + C**2 )
          DETWEI( GI ) = DETJ * SWEIGH( GI )
          SAREA = SAREA + DETWEI( GI )

          ! Calculate the normal at the Gauss pts:
          !  TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
          ! Perform x-product. N=T1 x T2
          CALL NORMGI( NORMXN(GI), NORMYN(GI), NORMZN(GI), &
               DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
               NORMX, NORMY, NORMZ)  

       END DO Loop_Gauss3D

    ELSE

       TWOPIE = 1.0 
       IF( DCYL ) TWOPIE = 2. * PIE
       !
       DO GI = 1, SNGI
          RGI = 0.
          DXDLX = 0.
          DXDLY = 0.
          DYDLX = 0.
          DYDLY = 0.
          DZDLX = 0.
          DZDLY=1. ! DZDLY=1 is to calculate the normal. 

          DO SL = 1, SNLOC

             IGLX =  SNDGLN(( ELE - 1 ) * SNLOC + SL )
             DXDLX = DXDLX + SNLX( SL, GI ) * X( IGLX ) 
             DYDLX = DYDLX + SNLX( SL, GI ) * Y( IGLX )

             RGI =   RGI + SN( SL, GI ) * Y( IGLX )

          END DO

          IF( .NOT. DCYL) RGI = 1.0 
          DETJ = SQRT( DXDLX**2 + DYDLX**2 )
          DETWEI( GI ) = TWOPIE *RGI * DETJ * SWEIGH( GI )
          SAREA = SAREA + DETWEI( GI ) 

          ! Calculate the normal at the Gauss pts:
          ! TANX1=DXDLX,TANY1=DYDLX,TANZ1=DZDLX,    TANX2=DXDLY,TANY2=DYDLY,TANZ2=DZDLY
          ! Perform x-product. N=T1 x T2

          NORMZ = 0.
          CALL NORMGI( NORMXN( GI ), NORMYN( GI ), NORMZN( GI ), &
               DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
               NORMX, NORMY, NORMZ ) 

       END DO

    ENDIF Conditional_Dimension

    RETURN

  END SUBROUTINE SDETNX2



  SUBROUTINE NORMGI( NORMXN, NORMYN, NORMZN, &
       DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
       NORMX, NORMY, NORMZ) 
    ! Calculate the normal at the Gauss pts
    ! Perform x-product. N=T1 x T2
    implicit none
    REAL, intent( inout ) :: NORMXN, NORMYN, NORMZN
    REAL, intent( in )    :: DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY
    REAL, intent( in )    :: NORMX, NORMY, NORMZ
    ! Local variables
    REAL :: RN, SIRN

    CALL XPROD( NORMXN, NORMYN, NORMZN, &
         DXDLX, DYDLX, DZDLX, &
         DXDLY, DYDLY, DZDLY )

    RN = SQRT( NORMXN**2 + NORMYN**2 + NORMZN**2 )

    SIRN = SIGN( 1.0 / RN, NORMXN * NORMX + NORMYN * NORMY + NORMZN * NORMZ )

    NORMXN = SIRN * NORMXN
    NORMYN = SIRN * NORMYN
    NORMZN = SIRN * NORMZN

    RETURN

  END SUBROUTINE NORMGI




  SUBROUTINE XPROD( AX, AY, AZ, &
       BX, BY, BZ, &
       CX, CY, CZ )
    implicit none
    REAL, intent( inout ) :: AX, AY, AZ
    REAL, intent( in )    :: BX, BY, BZ, CX, CY, CZ

    ! Perform x-product. a=b x c 
    AX =    BY * CZ - BZ * CY
    AY = -( BX * CZ - BZ * CX )
    AZ =    BX * CY - BY * CX

    RETURN
  END subroutine XPROD




  REAL FUNCTION R2NORM( VEC, NVEC )
    IMPLICIT NONE
    INTEGER :: NVEC
    REAL, DIMENSION( NVEC ) :: VEC
    ! Local variables
    INTEGER :: I
    REAL :: RSUM

    RSUM = 0.0
    DO I = 1, NVEC
       RSUM = RSUM + VEC( I )**2
    END DO
    R2NORM = SQRT( RSUM )

    RETURN
  END FUNCTION R2NORM



  real function tolfun(value)
    ! This function is a tolerance function for a value which is used as a denominator.
    ! If the absolute value of VALUE less than 1E-10, then it returns SIGN(A,B) i.e. 
    ! the absolute value of A times the sign of B where A is TOLERANCE and B is VALUE.

    implicit none
    real :: value
    ! Local
    real, parameter :: tolerance = 1.e-10

    if( abs( value ) < tolerance ) then
       tolfun = sign( tolerance, value )
    else
       tolfun = value
    endif

    return

  end function tolfun




  SUBROUTINE DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX,  &
       CV_NLOC, MAT_NLOC, CV_NONODS, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
       SCVFEN,SCVNGI,GI,IPHASE,NDIM,TDIFFUSION,HDC, T,TOLD,CV_NODJ_IPHA,CV_NODI_IPHA,ELE,ELE2, &
       CVNORMX,CVNORMY,CVNORMZ, &
       DTX_ELE,DTY_ELE,DTZ_ELE,DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE, &
       SELE,STOTEL,WIC_T_BC,WIC_T_BC_DIRICHLET, CV_OTHER_LOC,MAT_OTHER_LOC )
    ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
    ! based on a non-linear method and a non-oscillating scheme.
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NLOC, MAT_NLOC, CV_NONODS,NPHASE, TOTELE, MAT_NONODS, &
         SCVNGI,GI,IPHASE,NDIM,CV_NODJ_IPHA,CV_NODI_IPHA,ELE,ELE2, &
         SELE,STOTEL,WIC_T_BC_DIRICHLET
    REAL, intent( in ) :: HDC
    REAL, intent( inout ) :: DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
    INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) ::MAT_NDGLN
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::WIC_T_BC
    INTEGER, DIMENSION( CV_NLOC ), intent( in ) ::CV_OTHER_LOC
    INTEGER, DIMENSION( MAT_NLOC ), intent( in ) ::MAT_OTHER_LOC
    REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) ::T,TOLD
    REAL, DIMENSION( CV_NLOC,SCVNGI  ), intent( in ) :: SCVFEN
    REAL, DIMENSION( MAT_NONODS,NDIM,NDIM,NPHASE  ), intent( in ) :: TDIFFUSION
    REAL, DIMENSION( CV_NLOC, NPHASE, TOTELE ), intent( in ) :: DTX_ELE,DTY_ELE,DTZ_ELE, &
         DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
    REAL, DIMENSION( SCVNGI ), intent( in ) :: CVNORMX,CVNORMY,CVNORMZ

    ! local variables
    !        ===>  REALS  <===
    ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use 
    ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction. 
    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.05, DIFF_MAX_FRAC = 20.0
    REAL, DIMENSION( : , : ), allocatable :: DIFF_GI,DIFF_GI2
    REAL :: DTDX_GI,DTDY_GI,DTDZ_GI,DTOLDDX_GI,DTOLDDY_GI,DTOLDDZ_GI, &
         DTDX_GI2,DTDY_GI2,DTDZ_GI2,DTOLDDX_GI2,DTOLDDY_GI2,DTOLDDZ_GI2, &
         N_DOT_DKDT,N_DOT_DKDTOLD,  &
         N_DOT_DKDT2,N_DOT_DKDTOLD2, &
         DIFF_STAND_DIVDX,DIFF_STAND_DIVDX2
    INTEGER :: CV_KLOC,CV_KLOC2,MAT_KLOC,MAT_KLOC2,MAT_NODK,MAT_NODK2
    LOGICAL :: ZER_DIFF

    ALLOCATE( DIFF_GI(3,3) )
    ALLOCATE( DIFF_GI2(3,3) )

    ZER_DIFF=.FALSE.
    IF(SELE /= 0) ZER_DIFF= (WIC_T_BC(SELE+(IPHASE-1)*STOTEL) /= WIC_T_BC_DIRICHLET)

    Cond_ZerDiff: IF(ZER_DIFF) THEN

       DIFF_COEF_DIVDX    = 0.0 
       DIFF_COEFOLD_DIVDX = 0.0

    ELSE

       DTDX_GI = 0.0
       DTDY_GI = 0.0
       DTDZ_GI = 0.0
       DTOLDDX_GI = 0.0
       DTOLDDY_GI = 0.0
       DTOLDDZ_GI = 0.0

       DO CV_KLOC = 1, CV_NLOC
          DTDX_GI = DTDX_GI + SCVFEN( CV_KLOC, GI ) * DTX_ELE(CV_KLOC,IPHASE,ELE)
          DTDY_GI = DTDY_GI + SCVFEN( CV_KLOC, GI ) * DTY_ELE(CV_KLOC,IPHASE,ELE)
          DTDZ_GI = DTDZ_GI + SCVFEN( CV_KLOC, GI ) * DTZ_ELE(CV_KLOC,IPHASE,ELE)

          DTOLDDX_GI = DTOLDDX_GI + SCVFEN( CV_KLOC, GI ) * DTOLDX_ELE(CV_KLOC,IPHASE,ELE)
          DTOLDDY_GI = DTOLDDY_GI + SCVFEN( CV_KLOC, GI ) * DTOLDY_ELE(CV_KLOC,IPHASE,ELE)
          DTOLDDZ_GI = DTOLDDZ_GI + SCVFEN( CV_KLOC, GI ) * DTOLDZ_ELE(CV_KLOC,IPHASE,ELE)
       END DO

       DIFF_GI(:,:) = 0.0
       DO MAT_KLOC = 1, MAT_NLOC
          MAT_NODK = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_KLOC )
          DIFF_GI( 1:NDIM , 1:NDIM ) = DIFF_GI( 1:NDIM , 1:NDIM ) &
               + SCVFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK, 1:NDIM , 1:NDIM , IPHASE )
       END DO

       N_DOT_DKDT=CVNORMX(GI)*(DIFF_GI(1,1)*DTDX_GI+DIFF_GI(1,2)*DTDY_GI+DIFF_GI(1,3)*DTDZ_GI) &
            +CVNORMY(GI)*(DIFF_GI(2,1)*DTDX_GI+DIFF_GI(2,2)*DTDY_GI+DIFF_GI(2,3)*DTDZ_GI) &
            +CVNORMZ(GI)*(DIFF_GI(3,1)*DTDX_GI+DIFF_GI(3,2)*DTDY_GI+DIFF_GI(3,3)*DTDZ_GI) 

       N_DOT_DKDTOLD=CVNORMX(GI)*(DIFF_GI(1,1)*DTOLDDX_GI  &
            +DIFF_GI(1,2)*DTOLDDY_GI+DIFF_GI(1,3)*DTOLDDZ_GI) &
            +CVNORMY(GI)*(DIFF_GI(2,1)*DTOLDDX_GI  &
            +DIFF_GI(2,2)*DTOLDDY_GI+DIFF_GI(2,3)*DTOLDDZ_GI) &
            +CVNORMZ(GI)*(DIFF_GI(3,1)*DTOLDDX_GI  &
            +DIFF_GI(3,2)*DTOLDDY_GI+DIFF_GI(3,3)*DTOLDDZ_GI)

       ! This is the minimum diffusion...
       DIFF_STAND_DIVDX=( ABS(CVNORMX(GI))*DIFF_GI(1,1) &
            +ABS(CVNORMY(GI))*DIFF_GI(2,2) + ABS(CVNORMZ(GI))*DIFF_GI(3,3) ) /HDC

       Conditional_MAT_DISOPT_ELE2: IF( ( ELE2 /= 0 ).AND.( ELE2 /= ELE) ) THEN
          DIFF_GI2(:,:) = 0.0
          DO MAT_KLOC = 1, MAT_NLOC
             MAT_KLOC2 = MAT_OTHER_LOC( MAT_KLOC )
             IF(MAT_KLOC2 /= 0 )THEN
                MAT_NODK2 = MAT_NDGLN(( ELE2 - 1 ) * MAT_NLOC + MAT_KLOC2 )
                DIFF_GI2( 1:NDIM, 1:NDIM )=DIFF_GI2( 1:NDIM, 1:NDIM ) +SCVFEN( MAT_KLOC, GI ) &
                     *TDIFFUSION(MAT_NODK2, 1:NDIM, 1:NDIM ,IPHASE)
             ENDIF
          END DO

          DTDX_GI2 = 0.0
          DTDY_GI2 = 0.0
          DTDZ_GI2 = 0.0
          DTOLDDX_GI2 = 0.0
          DTOLDDY_GI2 = 0.0
          DTOLDDZ_GI2 = 0.0
          DO CV_KLOC = 1, CV_NLOC
             CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
             IF(CV_KLOC2 /= 0 )THEN
                DTDX_GI2 = DTDX_GI2 + SCVFEN( CV_KLOC, GI ) * DTX_ELE(CV_KLOC2,IPHASE,ELE2)
                DTDY_GI2 = DTDY_GI2 + SCVFEN( CV_KLOC, GI ) * DTY_ELE(CV_KLOC2,IPHASE,ELE2)
                DTDZ_GI2 = DTDZ_GI2 + SCVFEN( CV_KLOC, GI ) * DTZ_ELE(CV_KLOC2,IPHASE,ELE2)

                DTOLDDX_GI2 = DTOLDDX_GI2 + SCVFEN( CV_KLOC, GI ) *DTOLDX_ELE(CV_KLOC2,IPHASE,ELE2)
                DTOLDDY_GI2 = DTOLDDY_GI2 + SCVFEN( CV_KLOC, GI ) *DTOLDY_ELE(CV_KLOC2,IPHASE,ELE2)
                DTOLDDZ_GI2 = DTOLDDZ_GI2 + SCVFEN( CV_KLOC, GI ) *DTOLDZ_ELE(CV_KLOC2,IPHASE,ELE2)
             ENDIF
          END DO
          N_DOT_DKDT2=CVNORMX(GI)*(DIFF_GI2(1,1)*DTDX_GI2  &
               +DIFF_GI2(1,2)*DTDY_GI2+DIFF_GI2(1,3)*DTDZ_GI2) &
               +CVNORMY(GI)*(DIFF_GI2(2,1)*DTDX_GI2  &
               +DIFF_GI2(2,2)*DTDY_GI2+DIFF_GI2(2,3)*DTDZ_GI2) &
               +CVNORMZ(GI)*(DIFF_GI2(3,1)*DTDX_GI2  &
               +DIFF_GI2(3,2)*DTDY_GI2+DIFF_GI2(3,3)*DTDZ_GI2) 

          N_DOT_DKDTOLD2=CVNORMX(GI)*(DIFF_GI2(1,1)*DTOLDDX_GI2  &
               +DIFF_GI2(1,2)*DTOLDDY_GI2+DIFF_GI2(1,3)*DTOLDDZ_GI2) &
               +CVNORMY(GI)*(DIFF_GI2(2,1)*DTOLDDX_GI2  &
               +DIFF_GI2(2,2)*DTOLDDY_GI2+DIFF_GI2(2,3)*DTOLDDZ_GI2) &
               +CVNORMZ(GI)*(DIFF_GI2(3,1)*DTOLDDX_GI2  &
               +DIFF_GI2(3,2)*DTOLDDY_GI2+DIFF_GI2(3,3)*DTOLDDZ_GI2) 

          ! This is the minimum diffusion...
          DIFF_STAND_DIVDX2 = ( ABS(CVNORMX(GI))*DIFF_GI2(1,1) &
               +ABS(CVNORMY(GI))*DIFF_GI2(2,2) + ABS(CVNORMZ(GI))*DIFF_GI2(3,3) ) /HDC

          N_DOT_DKDT = 0.5*( N_DOT_DKDT + N_DOT_DKDT2 )
          N_DOT_DKDTOLD= 0.5*( N_DOT_DKDTOLD + N_DOT_DKDTOLD2 )

          !   ewrite(3,*)'DToldDX_GI,DToldDX_GI2,N_DOT_DKDTOLD,N_DOT_DKDTOLD2:',  &
          !            DToldDX_GI,DToldDX_GI2,N_DOT_DKDTOLD,N_DOT_DKDTOLD2

          ! This is the minimum diffusion...
          DIFF_STAND_DIVDX = MIN( DIFF_STAND_DIVDX, DIFF_STAND_DIVDX2 ) 

       ENDIF Conditional_MAT_DISOPT_ELE2


       DIFF_COEF_DIVDX    = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX, N_DOT_DKDT / &
            TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA )) )
       DIFF_COEFOLD_DIVDX = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX, N_DOT_DKDTOLD /  &
            TOLFUN( TOLD( CV_NODJ_IPHA ) - TOLD( CV_NODI_IPHA )))

       ! Make sure the diffusion has an upper bound...       
       DIFF_COEF_DIVDX    = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX, DIFF_COEF_DIVDX )
       DIFF_COEFOLD_DIVDX = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX, DIFF_COEFOLD_DIVDX )

    END IF Cond_ZerDiff

    !    ewrite(3,*)'HDC,DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX:', &
    !             HDC,DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
    !    ewrite(3,*)'TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA )):',  &
    !             TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA ))

    DEALLOCATE( DIFF_GI )
    DEALLOCATE( DIFF_GI2 )

    RETURN            

  END SUBROUTINE DIFFUS_CAL_COEFF





  SUBROUTINE GET_INT_VEL( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
       HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
       T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
       CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
       WIC_U_BC_DIRICHLET, &
       UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
       VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
       MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
       FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
       TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
       IN_ELE_UPWIND, DG_ELE_UPWIND )
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, WIC_U_BC_DIRICHLET, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI
    REAL, intent( inout ) :: NDOTQ, INCOME, NDOTQOLD, INCOMEOLD
    INTEGER, DIMENSION( TOTELE*U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( U_NLOC ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( U_SNLOC ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_U_BC
    INTEGER, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
         TOLDMAX_NOD
    REAL, DIMENSION( U_NLOC, SCVNGI  ), intent( in ) :: SUFEN
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: T, TOLD, FEMT, FEMTOLD, DEN, DENOLD
    REAL, DIMENSION( U_NONODS * NPHASE  ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
    REAL, DIMENSION( SCVNGI ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( U_NLOC ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( TOTELE  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( TOTELE*CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( CV_NLOC ), intent( in ) :: CV_OTHER_LOC
    REAL, DIMENSION( CV_NONODS  ), intent( in ) :: MASS_CV
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( CV_NONODS*NPHASE  ), intent( inout ) :: UP_WIND_NOD
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN, TOLDMIN, TMAX, TOLDMAX

    IF(CV_ELE_TYPE == 2) THEN
       ! For overlapping basis function approach.
       CALL GET_INT_VEL_OVERLAP( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
            HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
            WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
            TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
            IN_ELE_UPWIND, DG_ELE_UPWIND )
    ElSE
       CALL GET_INT_VEL_ORIG( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
            HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
            WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
            TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
            IN_ELE_UPWIND, DG_ELE_UPWIND )
    ENDIF
    RETURN
  END SUBROUTINE GET_INT_VEL





  SUBROUTINE GET_INT_VEL_OVERLAP( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
       HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
       T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
       CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
       WIC_U_BC_DIRICHLET, &
       UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
       VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
       MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
       FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
       TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
       IN_ELE_UPWIND, DG_ELE_UPWIND )
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    ! it assumes an overlapping decomposition approach for velocity. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, WIC_U_BC_DIRICHLET, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI
    REAL, intent( inout ) :: NDOTQ, INCOME, NDOTQOLD, INCOMEOLD
    INTEGER, DIMENSION( TOTELE*U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( U_NLOC ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( U_SNLOC ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_U_BC
    INTEGER, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
         TOLDMAX_NOD
    REAL, DIMENSION( U_NLOC, SCVNGI  ), intent( in ) :: SUFEN
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: T, TOLD, FEMT, FEMTOLD, DEN, DENOLD
    REAL, DIMENSION( U_NONODS * NPHASE  ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
    REAL, DIMENSION( SCVNGI ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( U_NLOC ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( TOTELE  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( TOTELE*CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( CV_NLOC ), intent( in ) :: CV_OTHER_LOC
    REAL, DIMENSION( CV_NONODS  ), intent( in ) :: MASS_CV
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( CV_NONODS*NPHASE  ), intent( inout ) :: UP_WIND_NOD
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN, TOLDMIN, TMAX, TOLDMAX

    ! Local variables
    REAL :: UDGI,VDGI,WDGI,UOLDDGI,VOLDDGI,WOLDDGI,  &
         UDGI2,VDGI2,WDGI2,UOLDDGI2,VOLDDGI2,WOLDDGI2, DT_I,DT_J,DTOLD_I,DTOLD_J, &
         UDGI_INT,VDGI_INT,WDGI_INT,UOLDDGI_INT,VOLDDGI_INT,WOLDDGI_INT, &
         NDOTQ_INT,NDOTQOLD_INT,FEMTOLDGI2,NDOTQ2,NDOTQOLD2, &
         FEMTOLDGI_IPHA, OVER_RELAX, ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
         GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA, V_NODI, G_NODI, V_NODJ, G_NODJ, &
         GEOMTOLDGI_IPHA, W_UPWIND, INCOME3, INCOME4, INCOMEOLD3, INCOMEOLD4, &
         LIMT3, LIMTOLD3, FEMTGI_IPHA, GEOMTGI_IPHA, UPWIND_FRAC
    REAL :: NVEC(3)
    INTEGER :: U_KLOC,U_NODK,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_SKLOC, &
         U_SNODK,U_SNODK_IPHA, II,  &
         U_KLOC_LEV, U_NLOC_LEV, U_SKLOC_LEV, U_SNLOC_LEV, CV_KLOC, CV_KNOD, &
         CV_KLOC2, CV_KNOD2, CV_NODI, CV_NODJ, IDIM, JDIM, IJ, MAT_NODI, MAT_NODJ
    LOGICAL :: CONSERV, MAX_OPER, SAT_BASED
    ! IN_ELE_UPWIND=1 switches on upwinding within and element (=3 recommended).
    !    INTEGER, PARAMETER :: IN_ELE_UPWIND = 3
    !    INTEGER, PARAMETER :: IN_ELE_UPWIND = 2
    ! DG_ELE_UPWIND=1 switches on upwinding between elements (=3 recommended).
    !    INTEGER, PARAMETER :: DG_ELE_UPWIND = 3
    !    INTEGER, PARAMETER :: DG_ELE_UPWIND = 2
    LOGICAL, PARAMETER :: LIM_VOL_ADJUST = .true.
    LOGICAL :: RESET_STORE
    REAL :: TMIN_STORE, TMAX_STORE, TOLDMIN_STORE, TOLDMAX_STORE
    ! coefficients for this element ELE
    UGI_COEF_ELE=0.0
    VGI_COEF_ELE=0.0
    WGI_COEF_ELE=0.0 

    ! coefficients for this element ELE2
    UGI_COEF_ELE2=0.0
    VGI_COEF_ELE2=0.0
    WGI_COEF_ELE2=0.0 

    U_NLOC_LEV =U_NLOC /CV_NLOC
    U_SNLOC_LEV=U_SNLOC/CV_NLOC


    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain. 
       IF( WIC_U_BC( SELE + ( IPHASE - 1 ) * STOTEL) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UOLDDGI = 0.0
          VOLDDGI = 0.0
          WOLDDGI = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
             UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
             VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
             WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )
             UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+1.0
             VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+1.0
             WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+1.0
             UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK_IPHA ) 
             VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK_IPHA ) 
             WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK_IPHA ) 
          END DO
       ELSE ! Specified vel bc.
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UOLDDGI = 0.0
          VOLDDGI = 0.0
          WOLDDGI = 0.0
          DO U_SKLOC_LEV = 1, U_SNLOC_LEV
             U_SKLOC = (CV_ILOC-1)*U_SNLOC_LEV + U_SKLOC_LEV
             U_KLOC = U_SLOC2LOC( U_SKLOC )
             U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
             U_NODK_IPHA = U_NODK + ( IPHASE - 1 ) * U_NONODS
             U_SNODK = ( SELE - 1 ) * U_SNLOC + U_SKLOC 
             U_SNODK_IPHA = U_SNODK + ( IPHASE - 1 ) * STOTEL*U_SNLOC
             IF(WIC_U_BC(SELE+(IPHASE-1)*STOTEL) == 10) THEN
                UDGI = UDGI + SUFEN( U_KLOC, GI ) * 0.5*(NU( U_NODK_IPHA )+SUF_U_BC( U_SNODK_IPHA ))
                VDGI = VDGI + SUFEN( U_KLOC, GI ) * 0.5*(NV( U_NODK_IPHA )+SUF_V_BC( U_SNODK_IPHA ))
                WDGI = WDGI + SUFEN( U_KLOC, GI ) * 0.5*(NW( U_NODK_IPHA )+SUF_W_BC( U_SNODK_IPHA ))
                UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+0.5
                VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+0.5
                WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+0.5
                UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * 0.5*(NUOLD( U_NODK_IPHA )+SUF_U_BC( U_SNODK_IPHA ))
                VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * 0.5*(NVOLD( U_NODK_IPHA )+SUF_V_BC( U_SNODK_IPHA ))
                WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * 0.5*(NWOLD( U_NODK_IPHA )+SUF_W_BC( U_SNODK_IPHA ))
             ELSE
                UDGI = UDGI + SUFEN( U_KLOC, GI ) * SUF_U_BC( U_SNODK_IPHA )
                VDGI = VDGI + SUFEN( U_KLOC, GI ) * SUF_V_BC( U_SNODK_IPHA )
                WDGI = WDGI + SUFEN( U_KLOC, GI ) * SUF_W_BC( U_SNODK_IPHA )
                UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * SUF_U_BC( U_SNODK_IPHA ) 
                VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * SUF_V_BC( U_SNODK_IPHA ) 
                WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * SUF_W_BC( U_SNODK_IPHA )
             ENDIF
          END DO
       ENDIF

    ELSE ! Conditional_SELE. Not on the boundary of the domain.
       Conditional_ELE2: IF(( ELE2 == 0 ).OR.( ELE2 == ELE)) THEN
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UOLDDGI = 0.0
          VOLDDGI = 0.0
          WOLDDGI = 0.0
          UDGI2 = 0.0
          VDGI2 = 0.0
          WDGI2 = 0.0
          UOLDDGI2 = 0.0
          VOLDDGI2 = 0.0
          WOLDDGI2 = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC =(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_KLOC2=(CV_JLOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_NODK_IPHA  = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
             U_NODK2_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC2 )+(IPHASE-1)*U_NONODS
             UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
             VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
             WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )
             UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK_IPHA ) 
             VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK_IPHA ) 
             WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK_IPHA ) 

             UDGI2 = UDGI2 + SUFEN( U_KLOC2, GI ) * NU( U_NODK2_IPHA )
             VDGI2 = VDGI2 + SUFEN( U_KLOC2, GI ) * NV( U_NODK2_IPHA )
             WDGI2 = WDGI2 + SUFEN( U_KLOC2, GI ) * NW( U_NODK2_IPHA )
             UOLDDGI2 = UOLDDGI2 + SUFEN( U_KLOC2, GI ) * NUOLD( U_NODK2_IPHA ) 
             VOLDDGI2 = VOLDDGI2 + SUFEN( U_KLOC2, GI ) * NVOLD( U_NODK2_IPHA ) 
             WOLDDGI2 = WOLDDGI2 + SUFEN( U_KLOC2, GI ) * NWOLD( U_NODK2_IPHA ) 
          END DO

          NDOTQ =  0.5*( CVNORMX( GI ) * (UDGI+UDGI2) + CVNORMY( GI ) * (VDGI+VDGI2)  &
               + CVNORMZ(GI) * (WDGI+WDGI2) )
          NDOTQOLD =0.5*(CVNORMX( GI ) * (UOLDDGI+UOLDDGI2) + CVNORMY( GI ) * (VOLDDGI+VOLDDGI2) &
               + CVNORMZ(GI) * (WOLDDGI+WOLDDGI2) )

          CV_NODI=CV_NODI_IPHA -(IPHASE-1)*CV_NONODS
          CV_NODJ=CV_NODJ_IPHA -(IPHASE-1)*CV_NONODS


          IF(IN_ELE_UPWIND==1) THEN
             IF(NDOTQ < 0.0) THEN
                INCOME=1.0
             ELSE
                INCOME=0.0
             ENDIF

             IF(NDOTQOLD < 0.0) THEN
                INCOMEOLD=1.0
             ELSE
                INCOMEOLD=0.0
             ENDIF
          ELSE IF(IN_ELE_UPWIND==2) THEN ! the best
             UPWIND_FRAC=0.8
             IF(NDOTQ < 0.0) THEN
                INCOME=0.8
                !         INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ELSE
                INCOME=0.2
                !         INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ENDIF
             IF(NDOTQOLD < 0.0) THEN
                INCOMEOLD=0.8
                !         INCOMEOLD=1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ELSE
                INCOMEOLD=0.2
                !         INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ENDIF
          ELSE IF(IN_ELE_UPWIND==3) THEN ! the best optimal upwind frac.

             NDOTQ = CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI  &
                  + CVNORMZ(GI) * WDGI 
             NDOTQOLD =CVNORMX( GI ) * UOLDDGI + CVNORMY( GI ) * VOLDDGI &
                  + CVNORMZ(GI) * WOLDDGI 
             NDOTQ2 =  CVNORMX( GI ) * UDGI2 + CVNORMY( GI ) * VDGI2  &
                  + CVNORMZ(GI) * WDGI2 
             NDOTQOLD2 =CVNORMX( GI ) * UOLDDGI2 + CVNORMY( GI ) * VOLDDGI2 &
                  + CVNORMZ(GI) * WOLDDGI2 

             MAT_NODI=MAT_NDGLN((ELE-1)*MAT_NLOC+CV_ILOC)
             MAT_NODJ=MAT_NDGLN((ELE-1)*MAT_NLOC+CV_JLOC)
             FEMTGI_IPHA = 0.0
             FEMTOLDGI_IPHA = 0.0
             DO CV_KLOC=1,CV_NLOC
                CV_KNOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                FEMTGI_IPHA = FEMTGI_IPHA &
                     + SCVFEN(CV_KLOC,GI)*FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
                !!     + SCVFEN(CV_KLOC,GI)*T(CV_KNOD+(IPHASE-1)*CV_NONODS)
                FEMTOLDGI_IPHA = FEMTOLDGI_IPHA &
                     + SCVFEN(CV_KLOC,GI)*FEMTOLD(CV_KNOD+(IPHASE-1)*CV_NONODS)
                !!     + SCVFEN(CV_KLOC,GI)*TOLD(CV_KNOD+(IPHASE-1)*CV_NONODS)
             END DO
             !      FEMTGI_IPHA = (MASS_CV(CV_NODJ)* T(CV_NODI_IPHA) + &
             !                  MASS_CV(CV_NODI)*T(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             !      FEMTOLDGI_IPHA = (MASS_CV(CV_NODJ)* TOLD(CV_NODI_IPHA) + &
             !                  MASS_CV(CV_NODI)*TOLD(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))

             !             GEOMTGI_IPHA = (MASS_CV(CV_NODJ)* T(CV_NODI_IPHA) + &
             !                  MASS_CV(CV_NODI)*T(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             !             GEOMTOLDGI_IPHA = (MASS_CV(CV_NODJ)* TOLD(CV_NODI_IPHA) + &
             !                  MASS_CV(CV_NODI)*TOLD(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ! Central is fine as its within an element with equally spaced nodes. 
             FEMTGI_IPHA = 0.5*(T(CV_NODI_IPHA) + T(CV_NODJ_IPHA) )
             FEMTOLDGI_IPHA = 0.5*(TOLD(CV_NODI_IPHA) + TOLD(CV_NODJ_IPHA) )
             GEOMTGI_IPHA = 0.5*(T(CV_NODI_IPHA) + T(CV_NODJ_IPHA) )
             GEOMTOLDGI_IPHA = 0.5*(TOLD(CV_NODI_IPHA) + TOLD(CV_NODJ_IPHA) )
             !             GEOMTGI_IPHA = FEMTGI_IPHA
             !                          GEOMTOLDGI_IPHA = FEMTOLDGI_IPHA 
             NVEC(1)=CVNORMX(GI)
             NVEC(2)=CVNORMY(GI)
             NVEC(3)=CVNORMZ(GI)
             ABS_CV_NODI_IPHA      = 0.0
             GRAD_ABS_CV_NODI_IPHA = 0.0
             ABS_CV_NODJ_IPHA      = 0.0
             GRAD_ABS_CV_NODJ_IPHA = 0.0
             DO IDIM=1,NDIM
                V_NODI=0.0
                G_NODI=0.0
                V_NODJ=0.0
                G_NODJ=0.0
                DO JDIM=1,NDIM
                   IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODI-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                   V_NODI = V_NODI + OPT_VEL_UPWIND_COEFS(IJ) * NVEC(JDIM)
                   G_NODI = G_NODI + OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM)*NVEC(JDIM)
                   IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODJ-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                   V_NODJ = V_NODJ + OPT_VEL_UPWIND_COEFS(IJ) * NVEC(JDIM)
                   G_NODJ = G_NODJ + OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM)*NVEC(JDIM)
                END DO
                ABS_CV_NODI_IPHA      = ABS_CV_NODI_IPHA + NVEC(IDIM)*V_NODI
                GRAD_ABS_CV_NODI_IPHA = GRAD_ABS_CV_NODI_IPHA + NVEC(IDIM)*G_NODI
                ABS_CV_NODJ_IPHA      = ABS_CV_NODJ_IPHA + NVEC(IDIM)*V_NODJ
                GRAD_ABS_CV_NODJ_IPHA = GRAD_ABS_CV_NODJ_IPHA + NVEC(IDIM)*G_NODJ
             END DO
             ! OPT_VEL_UPWIND_COEFS contains the coefficients
             !             IF(FACE_ITS>1) FEMTOLDGI_IPHA=0.5*(LIMT+FEMTOLDGI_IPHA)
             OVER_RELAX=1.0
             MAX_OPER=.TRUE.
             CONSERV=.true.
             SAT_BASED=.FALSE.
             !             IF(FACE_ITS>1) THEN
             IF(.false.) THEN
                !         FEMTOLDGI_IPHA=(LIMT-FEMTGI) + GEOMTOLDGI_IPHA
                FEMTOLDGI_IPHA=LIMT
                !         FEMTOLDGI_IPHA=sqrt(LIMT-FEMTOLDGI_IPHA) + FEMTOLDGI_IPHA
                !  GEOMTOLDGI_IPHA=FEMTGI

                OVER_RELAX=1.0
                MAX_OPER=.FALSE.
                SAT_BASED=.TRUE.
             ENDIF

             CALL FIND_OPT_INCOME_INTERP(INCOME, W_UPWIND, NDOTQ, NDOTQ2, IPHASE, &
                  T(CV_NODI_IPHA),T(CV_NODJ_IPHA), FEMTGI_IPHA, GEOMTGI_IPHA,OVER_RELAX, &
                  ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                  GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

             CALL FIND_OPT_INCOME_INTERP(INCOMEOLD, W_UPWIND, NDOTQOLD, NDOTQOLD2, IPHASE, &
                  TOLD(CV_NODI_IPHA),TOLD(CV_NODJ_IPHA), FEMTOLDGI_IPHA, GEOMTOLDGI_IPHA,OVER_RELAX, &
                  ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                  GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)



             if(.true.) then

                IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                   INCOME3=1.0
                ELSE
                   INCOME3=0.0
                ENDIF
                IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                   INCOMEOLD3=1.0
                ELSE
                   INCOMEOLD3=0.0
                ENDIF
                IF(LIM_VOL_ADJUST) THEN
                   RESET_STORE=.FALSE. 
                   CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
                   CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )

                   CALL CAL_LIM_VOL_ADJUST(TOLDMIN_STORE,TOLDMIN,TOLD,TOLDMIN_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                   CALL CAL_LIM_VOL_ADJUST(TOLDMAX_STORE,TOLDMAX,TOLD,TOLDMAX_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                ENDIF


                CALL ONVDLIMsqrt( CV_NONODS, &
                     LIMT3, FEMTGI_IPHA, INCOME3, CV_NODI, CV_NODJ, &
                     T( (IPHASE-1)*CV_NONODS+1 ), TMIN( (IPHASE-1)*CV_NONODS+1 ), TMAX( (IPHASE-1)*CV_NONODS+1 ), .FALSE. , .FALSE.)

                CALL ONVDLIMsqrt( CV_NONODS, &
                     LIMTOLD3, FEMTOLDGI_IPHA, INCOMEOLD3, CV_NODI, CV_NODJ, &
                     TOLD( (IPHASE-1)*CV_NONODS+1 ), TOLDMIN( (IPHASE-1)*CV_NONODS+1 ), TOLDMAX( (IPHASE-1)*CV_NONODS+1 ), .FALSE., .FALSE. )


                INCOME4=(LIMT3 - T( CV_NODI_IPHA )) &
                     /TOLFUN(T( CV_NODJ_IPHA )- T( CV_NODI_IPHA ))

                INCOMEOLD4=(LIMTOLD3 - TOLD( CV_NODI_IPHA )) &
                     /TOLFUN(TOLD( CV_NODJ_IPHA )- TOLD( CV_NODI_IPHA ))


                INCOME    = INCOME3*MAX(INCOME4,INCOME) + (1.-INCOME3)*MIN(INCOME4,INCOME) 
                INCOMEOLD = INCOMEOLD3*MAX(INCOMEOLD4,INCOMEOLD) + (1.-INCOMEOLD3)*MIN(INCOMEOLD4,INCOMEOLD)     

                IF(LIM_VOL_ADJUST) THEN
                   RESET_STORE=.TRUE. 
                   CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
                   CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )

                   CALL CAL_LIM_VOL_ADJUST(TOLDMIN_STORE,TOLDMIN,TOLD,TOLDMIN_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                   CALL CAL_LIM_VOL_ADJUST(TOLDMAX_STORE,TOLDMAX,TOLD,TOLDMAX_NOD,RESET_STORE,MASS_CV, &
                        CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                ENDIF

             endif

          ELSE
             UPWIND_FRAC=0.5
             IF(NDOTQ < 0.0) THEN
                !         INCOME=0.8
                INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ELSE
                !         INCOME=0.2
                INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ENDIF
             IF(NDOTQOLD < 0.0) THEN
                !         INCOMEOLD=0.8
                INCOMEOLD=1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ELSE
                !         INCOMEOLD=0.2
                INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ENDIF
          ENDIF
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UOLDDGI = 0.0
          VOLDDGI = 0.0
          WOLDDGI = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC =(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_KLOC2=(CV_JLOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_NODK_IPHA  = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
             U_NODK2_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC2 )+(IPHASE-1)*U_NONODS

             UDGI = UDGI + SUFEN( U_KLOC, GI ) * (NU( U_NODK_IPHA )*(1.0-INCOME) &
                  + NU( U_NODK2_IPHA )*INCOME )
             VDGI = VDGI + SUFEN( U_KLOC, GI ) * (NV( U_NODK_IPHA )*(1.0-INCOME) &
                  + NV( U_NODK2_IPHA )*INCOME )
             WDGI = WDGI + SUFEN( U_KLOC, GI ) * (NW( U_NODK_IPHA )*(1.0-INCOME) &
                  + NW( U_NODK2_IPHA )*INCOME )
             UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * (NUOLD( U_NODK_IPHA )*(1.0-INCOMEOLD) &
                  +NUOLD( U_NODK2_IPHA )*INCOMEOLD  )
             VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * (NVOLD( U_NODK_IPHA )*(1.0-INCOMEOLD) &
                  +NVOLD( U_NODK2_IPHA )*INCOMEOLD  )
             WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * (NWOLD( U_NODK_IPHA )*(1.0-INCOMEOLD) &
                  +NWOLD( U_NODK2_IPHA )*INCOMEOLD  )

             UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+1.0-INCOME
             VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+1.0-INCOME
             WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+1.0-INCOME
             UGI_COEF_ELE(U_KLOC2)=UGI_COEF_ELE(U_KLOC2)+INCOME
             VGI_COEF_ELE(U_KLOC2)=VGI_COEF_ELE(U_KLOC2)+INCOME
             WGI_COEF_ELE(U_KLOC2)=WGI_COEF_ELE(U_KLOC2)+INCOME
          END DO


       ELSE ! Conditional_ELE2: IF( ELE2 /= 0 ) THEN

          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UOLDDGI = 0.0
          VOLDDGI = 0.0
          WOLDDGI = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
             UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
             VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
             WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )
             UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+1.0
             VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+1.0
             WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+1.0
             UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK_IPHA ) 
             VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK_IPHA ) 
             WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK_IPHA ) 
          END DO

          UDGI2 = 0.0
          VDGI2 = 0.0
          WDGI2 = 0.0
          UOLDDGI2 = 0.0
          VOLDDGI2 = 0.0
          WOLDDGI2 = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_KLOC2 = U_OTHER_LOC( U_KLOC )
             IF( U_KLOC2 /= 0 ) THEN
                U_NODK2_IPHA = U_NDGLN((ELE2-1)*U_NLOC+U_KLOC2) +(IPHASE-1)*U_NONODS
                UDGI2  = UDGI2 + SUFEN( U_KLOC, GI ) * NU( U_NODK2_IPHA )
                VDGI2  = VDGI2 + SUFEN( U_KLOC, GI ) * NV( U_NODK2_IPHA )
                WDGI2  = WDGI2 + SUFEN( U_KLOC, GI ) * NW( U_NODK2_IPHA )
                UGI_COEF_ELE2(U_KLOC2)=UGI_COEF_ELE2(U_KLOC2)+1.0
                VGI_COEF_ELE2(U_KLOC2)=VGI_COEF_ELE2(U_KLOC2)+1.0
                WGI_COEF_ELE2(U_KLOC2)=WGI_COEF_ELE2(U_KLOC2)+1.0
                UOLDDGI2  = UOLDDGI2 + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK2_IPHA )
                VOLDDGI2  = VOLDDGI2 + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK2_IPHA )
                WOLDDGI2  = WOLDDGI2 + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK2_IPHA )
             ENDIF
          END DO

          NDOTQ = CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI  &
               + CVNORMZ(GI) * WDGI 
          NDOTQOLD =CVNORMX( GI ) * UOLDDGI + CVNORMY( GI ) * VOLDDGI &
               + CVNORMZ(GI) * WOLDDGI 
          NDOTQ2 =  CVNORMX( GI ) * UDGI2 + CVNORMY( GI ) * VDGI2  &
               + CVNORMZ(GI) * WDGI2 
          NDOTQOLD2 =CVNORMX( GI ) * UOLDDGI2 + CVNORMY( GI ) * VOLDDGI2 &
               + CVNORMZ(GI) * WOLDDGI2 

          CV_NODI=CV_NODI_IPHA -(IPHASE-1)*CV_NONODS
          CV_NODJ=CV_NODJ_IPHA -(IPHASE-1)*CV_NONODS


          IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
             DT_I=1.0
             DT_J=1.0
             DTOLD_I=1.0
             DTOLD_J=1.0
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
             DT_I=MAX(1.E-2,T( CV_NODI_IPHA ))
             DT_J=MAX(1.E-2,T( CV_NODJ_IPHA ))
             DTOLD_I=MAX(1.E-2,TOLD( CV_NODI_IPHA ))
             DTOLD_J=MAX(1.E-2,TOLD( CV_NODJ_IPHA ))
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
             DT_I=DEN( CV_NODI_IPHA )*T( CV_NODI_IPHA )
             DT_J=DEN( CV_NODJ_IPHA )*T( CV_NODJ_IPHA )
             DTOLD_I=DENOLD( CV_NODI_IPHA )*TOLD( CV_NODI_IPHA )
             DTOLD_J=DENOLD( CV_NODJ_IPHA )*TOLD( CV_NODJ_IPHA )
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 4) THEN
             IF(DG_ELE_UPWIND==1) THEN
                IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                   INCOME=1.0
                ELSE
                   INCOME=0.0
                ENDIF
                IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                   INCOMEOLD=1.0
                ELSE
                   INCOMEOLD=0.0
                ENDIF
             ELSE IF(DG_ELE_UPWIND==2) THEN ! the best
                UPWIND_FRAC=0.8
                IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                   !            INCOME=0.8
                   INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ELSE
                   !                 INCOME=0.2
                   INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ENDIF
                IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                   !            INCOMEOLD=0.8
                   INCOMEOLD=1.0 -2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ELSE
                   !            INCOMEOLD=0.2
                   INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ENDIF
             ELSE IF(DG_ELE_UPWIND==3) THEN ! the best optimal upwind frac.

                MAT_NODI=MAT_NDGLN((ELE-1)*MAT_NLOC+CV_ILOC)
                MAT_NODJ=MAT_NDGLN((ELE2-1)*MAT_NLOC+CV_JLOC)
                FEMTGI_IPHA = 0.0
                FEMTOLDGI_IPHA = 0.0
                DO CV_KLOC = 1, CV_NLOC
                   CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                   IF(CV_KLOC2 /= 0 ) THEN
                      if(.false.) then
                         CV_KNOD =CV_NDGLN((ELE-1) *CV_NLOC+CV_KLOC)
                         CV_KNOD2=CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
                         FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) & 
                              *0.5*(FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS) &
                              +FEMT(CV_KNOD2+(IPHASE-1)*CV_NONODS))
                         FEMTOLDGI_IPHA = FEMTOLDGI_IPHA + SCVFEN(CV_KLOC,GI) & 
                              *0.5*(FEMTOLD(CV_KNOD+(IPHASE-1)*CV_NONODS) &
                              +FEMTOLD(CV_KNOD2+(IPHASE-1)*CV_NONODS))
                      endif
                      if(.true.) then
                         if(0.5*(NDOTQ+NDOTQ2)>0.0) then
                            CV_KNOD =CV_NDGLN((ELE-1) *CV_NLOC+CV_KLOC)
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 *FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
                            !!           *T(CV_KNOD+(IPHASE-1)*CV_NONODS)
                         else
                            CV_KNOD2=CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 *FEMT(CV_KNOD2+(IPHASE-1)*CV_NONODS)
                            !!           *T(CV_KNOD2+(IPHASE-1)*CV_NONODS)
                         endif
                         if(0.5*(NDOTQOLD+NDOTQOLD2)>0.0) then
                            CV_KNOD =CV_NDGLN((ELE-1) *CV_NLOC+CV_KLOC)
                            FEMTOLDGI_IPHA = FEMTOLDGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 *FEMTOLD(CV_KNOD+(IPHASE-1)*CV_NONODS)
                            !!           *TOLD(CV_KNOD+(IPHASE-1)*CV_NONODS)
                         else
                            CV_KNOD2=CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
                            FEMTOLDGI_IPHA = FEMTOLDGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 *FEMTOLD(CV_KNOD2+(IPHASE-1)*CV_NONODS)
                            !!           *TOLD(CV_KNOD2+(IPHASE-1)*CV_NONODS)
                         endif
                      endif
                   ENDIF
                END DO
                !      FEMTOLDGI_IPHA = (MASS_CV(CV_NODJ)* T(CV_NODI_IPHA) + &
                !         MASS_CV(CV_NODI)*T(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                !      FEMTOLDGI_IPHA = (MASS_CV(CV_NODi)* T(CV_NODI_IPHA) + &
                !         MASS_CV(CV_NODj)*T(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                !      FEMTOLDGI_IPHA = 0.5*(T(CV_NODI_IPHA) + T(CV_NODJ_IPHA) )
                GEOMTGI_IPHA = (MASS_CV(CV_NODJ)* T(CV_NODI_IPHA) + &
                     MASS_CV(CV_NODI)*T(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                GEOMTOLDGI_IPHA = (MASS_CV(CV_NODJ)* TOLD(CV_NODI_IPHA) + &
                     MASS_CV(CV_NODI)*TOLD(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                !                             GEOMTGI_IPHA =FEMTGI_IPHA
                !                             GEOMTOLDGI_IPHA =FEMTOLDGI_IPHA
                NVEC(1)=CVNORMX(GI)
                NVEC(2)=CVNORMY(GI)
                NVEC(3)=CVNORMZ(GI)
                ABS_CV_NODI_IPHA      = 0.0
                GRAD_ABS_CV_NODI_IPHA = 0.0
                ABS_CV_NODJ_IPHA      = 0.0
                GRAD_ABS_CV_NODJ_IPHA = 0.0
                DO IDIM=1,NDIM
                   V_NODI=0.0
                   G_NODI=0.0
                   V_NODJ=0.0
                   G_NODJ=0.0
                   DO JDIM=1,NDIM
                      IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODI-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                      V_NODI = V_NODI + OPT_VEL_UPWIND_COEFS(IJ) * NVEC(JDIM)
                      G_NODI = G_NODI + OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM)*NVEC(JDIM)
                      IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODJ-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                      V_NODJ = V_NODJ + OPT_VEL_UPWIND_COEFS(IJ) * NVEC(JDIM)
                      G_NODJ = G_NODJ + OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM)*NVEC(JDIM)
                   END DO
                   ABS_CV_NODI_IPHA      = ABS_CV_NODI_IPHA + NVEC(IDIM)*V_NODI
                   GRAD_ABS_CV_NODI_IPHA = GRAD_ABS_CV_NODI_IPHA + NVEC(IDIM)*G_NODI
                   ABS_CV_NODJ_IPHA      = ABS_CV_NODJ_IPHA + NVEC(IDIM)*V_NODJ
                   GRAD_ABS_CV_NODJ_IPHA = GRAD_ABS_CV_NODJ_IPHA + NVEC(IDIM)*G_NODJ
                END DO
                ! OPT_VEL_UPWIND_COEFS contains the coefficients...
                OVER_RELAX=1.
                MAX_OPER=.TRUE.
                CONSERV=.true.
                SAT_BASED=.FALSE.

                CALL FIND_OPT_INCOME_INTERP(INCOME, W_UPWIND, NDOTQ, NDOTQ2, IPHASE, &
                     T(CV_NODI_IPHA),T(CV_NODJ_IPHA), FEMTGI_IPHA, GEOMTGI_IPHA,OVER_RELAX, &
                     ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                     GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

                CALL FIND_OPT_INCOME_INTERP(INCOMEOLD, W_UPWIND, NDOTQOLD, NDOTQOLD2, IPHASE,  &
                     TOLD(CV_NODI_IPHA),TOLD(CV_NODJ_IPHA), FEMTOLDGI_IPHA, GEOMTOLDGI_IPHA,OVER_RELAX, &
                     ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                     GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)


                IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                   INCOME3=1.0
                ELSE
                   INCOME3=0.0
                ENDIF
                IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                   INCOMEOLD3=1.0
                ELSE
                   INCOMEOLD3=0.0
                ENDIF

                if(.true.) then

                   IF(LIM_VOL_ADJUST) THEN
                      RESET_STORE=.FALSE. 
                      CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
                      CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )

                      CALL CAL_LIM_VOL_ADJUST(TOLDMIN_STORE,TOLDMIN,TOLD,TOLDMIN_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                      CALL CAL_LIM_VOL_ADJUST(TOLDMAX_STORE,TOLDMAX,TOLD,TOLDMAX_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                   ENDIF


                   CALL ONVDLIMsqrt( CV_NONODS, &
                        LIMT3, FEMTGI_IPHA, INCOME3, CV_NODI, CV_NODJ, &
                        T( (IPHASE-1)*CV_NONODS+1 ), TMIN( (IPHASE-1)*CV_NONODS+1 ), TMAX( (IPHASE-1)*CV_NONODS+1 ), .FALSE. , .FALSE.)


                   CALL ONVDLIMsqrt( CV_NONODS, &
                        LIMTOLD3, FEMTOLDGI_IPHA, INCOMEOLD3, CV_NODI, CV_NODJ, &
                        TOLD( (IPHASE-1)*CV_NONODS+1 ), TOLDMIN( (IPHASE-1)*CV_NONODS+1 ), TOLDMAX( (IPHASE-1)*CV_NONODS+1 ), .FALSE., .FALSE. )


                   IF(LIM_VOL_ADJUST) THEN
                      RESET_STORE=.TRUE. 
                      CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
                      CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )

                      CALL CAL_LIM_VOL_ADJUST(TOLDMIN_STORE,TOLDMIN,TOLD,TOLDMIN_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                      CALL CAL_LIM_VOL_ADJUST(TOLDMAX_STORE,TOLDMAX,TOLD,TOLDMAX_NOD,RESET_STORE,MASS_CV, &
                           CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD3 )
                   ENDIF

                endif

                OVER_RELAX=1.
                MAX_OPER=.TRUE.
                CONSERV=.true.
                SAT_BASED=.FALSE.
                !             IF(FACE_ITS>1) THEN
                !             IF(.false.) THEN
                IF(.true.) THEN
                   FEMTGI_IPHA=LIMT3
                   FEMTOLDGI_IPHA=LIMTOLD3
                   !  GEOMTOLDGI_IPHA=FEMTGI
                   OVER_RELAX=1.0
                   MAX_OPER=.FALSE.
                   SAT_BASED=.TRUE.
                ENDIF

                CALL FIND_OPT_INCOME_INTERP(INCOME4, W_UPWIND, NDOTQ, NDOTQ2, IPHASE, &
                     T(CV_NODI_IPHA),T(CV_NODJ_IPHA), FEMTGI_IPHA, GEOMTGI_IPHA,OVER_RELAX, &
                     ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                     GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

                CALL FIND_OPT_INCOME_INTERP(INCOMEOLD4, W_UPWIND, NDOTQOLD, NDOTQOLD2, IPHASE,  &
                     TOLD(CV_NODI_IPHA),TOLD(CV_NODJ_IPHA), FEMTOLDGI_IPHA, GEOMTOLDGI_IPHA,OVER_RELAX, &
                     ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                     GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

                INCOME    = INCOME3*MAX(INCOME4,INCOME) + (1.-INCOME3)*MIN(INCOME4,INCOME) 
                INCOMEOLD = INCOMEOLD3*MAX(INCOMEOLD4,INCOMEOLD) + (1.-INCOMEOLD3)*MIN(INCOMEOLD4,INCOMEOLD)
             ELSE
                UPWIND_FRAC=0.5
                IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                   !            INCOME=0.8
                   INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ELSE
                   !                 INCOME=0.2
                   INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ENDIF
                IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                   !            INCOMEOLD=0.8
                   INCOMEOLD=1.0 -2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ELSE
                   !            INCOMEOLD=0.2
                   INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                ENDIF
             ENDIF

             DT_I = 1.0 - INCOME
             DT_J = 1.0 - DT_I

             DTOLD_I = 1.0 - INCOMEOLD
             DTOLD_J = 1.0 - DTOLD_I
          ENDIF
          ! Amend weighting for porosity only across elements...
          IF((ABS(CV_DG_VEL_INT_OPT ) == 2).OR.(ABS(CV_DG_VEL_INT_OPT ) == 3)) THEN 
             IF(ELE /= ELE2) THEN 
                DT_I=VOLFRA_PORE(ELE) *DT_I
                DT_J=VOLFRA_PORE(ELE2)*DT_J
                DTOLD_I=VOLFRA_PORE(ELE) *DTOLD_I
                DTOLD_J=VOLFRA_PORE(ELE2)*DTOLD_J
             ENDIF
          ENDIF

          !   dt_i=0.5
          !   dt_j=0.5
          !   dtold_i=0.5
          !   dtold_j=0.5


          UDGI_INT = ( DT_I * UDGI + DT_J * UDGI2 ) / (DT_I + DT_J)
          VDGI_INT = ( DT_I * VDGI + DT_J * VDGI2 ) / (DT_I + DT_J)
          WDGI_INT = ( DT_I * WDGI + DT_J * WDGI2 ) / (DT_I + DT_J)

          UOLDDGI_INT = ( DTOLD_I * UOLDDGI + DTOLD_J * UOLDDGI2 ) / (DTOLD_I + DTOLD_J)
          VOLDDGI_INT = ( DTOLD_I * VOLDDGI + DTOLD_J * VOLDDGI2 ) / (DTOLD_I + DTOLD_J)
          WOLDDGI_INT = ( DTOLD_I * WOLDDGI + DTOLD_J * WOLDDGI2 ) / (DTOLD_I + DTOLD_J)


          IF( CV_DG_VEL_INT_OPT < 0 ) THEN

             NDOTQ_INT = CVNORMX( GI ) * UDGI_INT + CVNORMY( GI ) * VDGI_INT + &
                  CVNORMZ( GI ) * WDGI_INT
             NDOTQOLD_INT = CVNORMX( GI ) * UOLDDGI_INT + CVNORMY( GI ) * VOLDDGI_INT + &
                  CVNORMZ( GI ) * WOLDDGI_INT

             IF( NDOTQ_INT <= 0.0 ) THEN  !Incoming
                !   DT_I=1.0
                DT_J=DT_I+DT_J
             ELSE
                DT_I=DT_I+DT_J
                !   DT_J=1.0
             ENDIF
             IF( NDOTQOLD_INT <= 0.0 ) THEN  !Incoming
                !   DTOLD_I=1.0
                DTOLD_J=DTOLD_I+DTOLD_J
             ELSE
                DTOLD_I=DTOLD_I+DTOLD_J
                !   DTOLD_J=1.0
             ENDIF

             UDGI_INT = ( DT_I * UDGI + DT_J * UDGI2 ) / (DT_I + DT_J)
             VDGI_INT = ( DT_I * VDGI + DT_J * VDGI2 ) / (DT_I + DT_J)
             WDGI_INT = ( DT_I * WDGI + DT_J * WDGI2 ) / (DT_I + DT_J)

             UOLDDGI_INT = ( DTOLD_I * UOLDDGI + DTOLD_J * UOLDDGI2 ) / (DTOLD_I + DTOLD_J)
             VOLDDGI_INT = ( DTOLD_I * VOLDDGI + DTOLD_J * VOLDDGI2 ) / (DTOLD_I + DTOLD_J)
             WOLDDGI_INT = ( DTOLD_I * WOLDDGI + DTOLD_J * WOLDDGI2 ) / (DTOLD_I + DTOLD_J)
          ENDIF

          UDGI  = UDGI_INT 
          VDGI  = VDGI_INT 
          WDGI  = WDGI_INT

          UGI_COEF_ELE(:)=DT_I * UGI_COEF_ELE(:) / (DT_I + DT_J)
          VGI_COEF_ELE(:)=DT_I * VGI_COEF_ELE(:) / (DT_I + DT_J)
          WGI_COEF_ELE(:)=DT_I * WGI_COEF_ELE(:) / (DT_I + DT_J)

          UGI_COEF_ELE2(:)=DT_J * UGI_COEF_ELE2(:) / (DT_I + DT_J)
          VGI_COEF_ELE2(:)=DT_J * VGI_COEF_ELE2(:) / (DT_I + DT_J)
          WGI_COEF_ELE2(:)=DT_J * WGI_COEF_ELE2(:) / (DT_I + DT_J)

          UOLDDGI  = UOLDDGI_INT 
          VOLDDGI  = VOLDDGI_INT 
          WOLDDGI  = WOLDDGI_INT

       ENDIF Conditional_ELE2

    ENDIF Conditional_SELE

    NDOTQ =  CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI + CVNORMZ(GI) * WDGI
    NDOTQOLD =  CVNORMX( GI ) * UOLDDGI + CVNORMY( GI ) * VOLDDGI + CVNORMZ(GI) * WOLDDGI

    ! Define whether flux is incoming or outgoing, depending on direction of flow
    IF( NDOTQ >  0.0 ) THEN
       INCOME = 0.0  !Outgoing
    ELSE
       INCOME = 1.0  !Incoming
    ENDIF
    IF( NDOTQOLD >  0.0 ) THEN
       INCOMEOLD = 0.0  !Outgoing
    ELSE
       INCOMEOLD = 1.0  !Incoming
    ENDIF

    RETURN  

  END SUBROUTINE GET_INT_VEL_OVERLAP




  SUBROUTINE GET_INT_VEL_ORIG( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
       HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
       T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
       CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
       WIC_U_BC_DIRICHLET, &
       UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
       VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
       MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
       FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
       TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
       IN_ELE_UPWIND, DG_ELE_UPWIND )
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, WIC_U_BC_DIRICHLET, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI
    REAL, intent( inout ) :: NDOTQ, INCOME, NDOTQOLD, INCOMEOLD
    INTEGER, DIMENSION( TOTELE*U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( U_NLOC ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( U_SNLOC ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_U_BC
    INTEGER, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
         TOLDMAX_NOD
    REAL, DIMENSION( U_NLOC, SCVNGI  ), intent( in ) :: SUFEN
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: T, TOLD, FEMT, FEMTOLD, DEN, DENOLD
    REAL, DIMENSION( U_NONODS * NPHASE  ), intent( in ) :: NU, NV, NW, NUOLD, NVOLD, NWOLD
    REAL, DIMENSION( SCVNGI ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( U_NLOC ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( TOTELE  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( TOTELE*CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( CV_NLOC ), intent( in ) :: CV_OTHER_LOC
    REAL, DIMENSION( CV_NONODS  ), intent( in ) :: MASS_CV
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( CV_NONODS*NPHASE  ), intent( inout ) :: UP_WIND_NOD
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN, TOLDMIN, TMAX, TOLDMAX

    ! Local variables
    REAL :: UDGI,VDGI,WDGI,UOLDDGI,VOLDDGI,WOLDDGI,  &
         UDGI2,VDGI2,WDGI2,UOLDDGI2,VOLDDGI2,WOLDDGI2, DT_I,DT_J,DTOLD_I,DTOLD_J, &
         UDGI_INT,VDGI_INT,WDGI_INT,UOLDDGI_INT,VOLDDGI_INT,WOLDDGI_INT, &
         NDOTQ_INT, NDOTQOLD_INT
    INTEGER :: U_KLOC,U_NODK,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_SKLOC, &
         U_SNODK,U_SNODK_IPHA, II
    ! coefficients for this element ELE
    UGI_COEF_ELE=0.0
    VGI_COEF_ELE=0.0
    WGI_COEF_ELE=0.0 

    ! coefficients for this element ELE2
    UGI_COEF_ELE2=0.0
    VGI_COEF_ELE2=0.0
    WGI_COEF_ELE2=0.0 


    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain. 
       IF( WIC_U_BC( SELE + ( IPHASE - 1 ) * STOTEL) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UOLDDGI = 0.0
          VOLDDGI = 0.0
          WOLDDGI = 0.0
          DO U_KLOC = 1, U_NLOC
             U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
             UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
             VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
             WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )
             UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+1.0
             VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+1.0
             WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+1.0
             UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK_IPHA ) 
             VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK_IPHA ) 
             WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK_IPHA ) 
          END DO
       ELSE ! Specified vel bc.
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UOLDDGI = 0.0
          VOLDDGI = 0.0
          WOLDDGI = 0.0
          DO U_SKLOC = 1, U_SNLOC
             U_KLOC = U_SLOC2LOC( U_SKLOC )
             U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
             U_NODK_IPHA = U_NODK + ( IPHASE - 1 ) * U_NONODS
             U_SNODK = ( SELE - 1 ) * U_SNLOC + U_SKLOC 
             U_SNODK_IPHA = U_SNODK + ( IPHASE - 1 ) * STOTEL*U_SNLOC
             IF(WIC_U_BC(SELE+(IPHASE-1)*STOTEL) == 10) THEN
                UDGI = UDGI + SUFEN( U_KLOC, GI ) * 0.5*(NU( U_NODK_IPHA )+SUF_U_BC( U_SNODK_IPHA ))
                VDGI = VDGI + SUFEN( U_KLOC, GI ) * 0.5*(NV( U_NODK_IPHA )+SUF_V_BC( U_SNODK_IPHA ))
                WDGI = WDGI + SUFEN( U_KLOC, GI ) * 0.5*(NW( U_NODK_IPHA )+SUF_W_BC( U_SNODK_IPHA ))
                UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+0.5
                VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+0.5
                WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+0.5
                UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * 0.5*(NUOLD( U_NODK_IPHA )+SUF_U_BC( U_SNODK_IPHA ))
                VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * 0.5*(NVOLD( U_NODK_IPHA )+SUF_V_BC( U_SNODK_IPHA ))
                WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * 0.5*(NWOLD( U_NODK_IPHA )+SUF_W_BC( U_SNODK_IPHA ))
             ELSE
                UDGI = UDGI + SUFEN( U_KLOC, GI ) * SUF_U_BC( U_SNODK_IPHA )
                VDGI = VDGI + SUFEN( U_KLOC, GI ) * SUF_V_BC( U_SNODK_IPHA )
                WDGI = WDGI + SUFEN( U_KLOC, GI ) * SUF_W_BC( U_SNODK_IPHA )
                UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * SUF_U_BC( U_SNODK_IPHA ) 
                VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * SUF_V_BC( U_SNODK_IPHA ) 
                WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * SUF_W_BC( U_SNODK_IPHA )
             ENDIF
          END DO
       ENDIF

    ELSE ! Conditional_SELE. Not on the boundary of the domain.
       UDGI = 0.0
       VDGI = 0.0
       WDGI = 0.0
       UOLDDGI = 0.0
       VOLDDGI = 0.0
       WOLDDGI = 0.0
       DO U_KLOC = 1, U_NLOC
          U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
          UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
          VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
          WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )
          UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+1.0
          VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+1.0
          WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+1.0
          UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK_IPHA ) 
          VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK_IPHA ) 
          WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK_IPHA ) 
       END DO

       Conditional_ELE2: IF( ELE2 /= 0 ) THEN
          UDGI2 = 0.0
          VDGI2 = 0.0
          WDGI2 = 0.0
          UOLDDGI2 = 0.0
          VOLDDGI2 = 0.0
          WOLDDGI2 = 0.0
          DO U_KLOC = 1, U_NLOC
             U_KLOC2 = U_OTHER_LOC( U_KLOC )
             IF( U_KLOC2 /= 0 ) THEN
                U_NODK2_IPHA = U_NDGLN((ELE2-1)*U_NLOC+U_KLOC2) +(IPHASE-1)*U_NONODS
                UDGI2  = UDGI2 + SUFEN( U_KLOC, GI ) * NU( U_NODK2_IPHA )
                VDGI2  = VDGI2 + SUFEN( U_KLOC, GI ) * NV( U_NODK2_IPHA )
                WDGI2  = WDGI2 + SUFEN( U_KLOC, GI ) * NW( U_NODK2_IPHA )
                UGI_COEF_ELE2(U_KLOC2)=UGI_COEF_ELE2(U_KLOC2)+1.0
                VGI_COEF_ELE2(U_KLOC2)=VGI_COEF_ELE2(U_KLOC2)+1.0
                WGI_COEF_ELE2(U_KLOC2)=WGI_COEF_ELE2(U_KLOC2)+1.0
                UOLDDGI2  = UOLDDGI2 + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK2_IPHA )
                VOLDDGI2  = VOLDDGI2 + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK2_IPHA )
                WOLDDGI2  = WOLDDGI2 + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK2_IPHA )
             ENDIF
          END DO

          IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
             DT_I=1.0
             DT_J=1.0
             DTOLD_I=1.0
             DTOLD_J=1.0
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
             DT_I=MAX(1.E-2,T( CV_NODI_IPHA ))
             DT_J=MAX(1.E-2,T( CV_NODJ_IPHA ))
             DTOLD_I=MAX(1.E-2,TOLD( CV_NODI_IPHA ))
             DTOLD_J=MAX(1.E-2,TOLD( CV_NODJ_IPHA ))
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
             DT_I=DEN( CV_NODI_IPHA )*T( CV_NODI_IPHA )
             DT_J=DEN( CV_NODJ_IPHA )*T( CV_NODJ_IPHA )
             DTOLD_I=DENOLD( CV_NODI_IPHA )*TOLD( CV_NODI_IPHA )
             DTOLD_J=DENOLD( CV_NODJ_IPHA )*TOLD( CV_NODJ_IPHA )
          ENDIF
          ! Amend weighting for porosity only across elements...
          IF(ABS(CV_DG_VEL_INT_OPT ) >= 2) THEN 
             IF(ELE /= ELE2) THEN 
                DT_I=VOLFRA_PORE(ELE) *DT_I
                DT_J=VOLFRA_PORE(ELE2)*DT_J
                DTOLD_I=VOLFRA_PORE(ELE) *DTOLD_I
                DTOLD_J=VOLFRA_PORE(ELE2)*DTOLD_J
             ENDIF
          ENDIF


          UDGI_INT = ( DT_I * UDGI + DT_J * UDGI2 ) / (DT_I + DT_J)
          VDGI_INT = ( DT_I * VDGI + DT_J * VDGI2 ) / (DT_I + DT_J)
          WDGI_INT = ( DT_I * WDGI + DT_J * WDGI2 ) / (DT_I + DT_J)

          UOLDDGI_INT = ( DTOLD_I * UOLDDGI + DTOLD_J * UOLDDGI2 ) / (DTOLD_I + DTOLD_J)
          VOLDDGI_INT = ( DTOLD_I * VOLDDGI + DTOLD_J * VOLDDGI2 ) / (DTOLD_I + DTOLD_J)
          WOLDDGI_INT = ( DTOLD_I * WOLDDGI + DTOLD_J * WOLDDGI2 ) / (DTOLD_I + DTOLD_J)


          IF( CV_DG_VEL_INT_OPT < 0 ) THEN

             NDOTQ_INT = CVNORMX( GI ) * UDGI_INT + CVNORMY( GI ) * VDGI_INT + &
                  CVNORMZ( GI ) * WDGI_INT
             NDOTQOLD_INT = CVNORMX( GI ) * UOLDDGI_INT + CVNORMY( GI ) * VOLDDGI_INT + &
                  CVNORMZ( GI ) * WOLDDGI_INT

             IF( NDOTQ_INT <= 0.0 ) THEN  !Incoming
                !   DT_I=1.0
                DT_J=DT_I+DT_J
             ELSE
                DT_I=DT_I+DT_J
                !   DT_J=1.0
             ENDIF
             IF( NDOTQOLD_INT <= 0.0 ) THEN  !Incoming
                !   DTOLD_I=1.0
                DTOLD_J=DTOLD_I+DTOLD_J
             ELSE
                DTOLD_I=DTOLD_I+DTOLD_J
                !   DTOLD_J=1.0
             ENDIF

             UDGI_INT = ( DT_I * UDGI + DT_J * UDGI2 ) / (DT_I + DT_J)
             VDGI_INT = ( DT_I * VDGI + DT_J * VDGI2 ) / (DT_I + DT_J)
             WDGI_INT = ( DT_I * WDGI + DT_J * WDGI2 ) / (DT_I + DT_J)

             UOLDDGI_INT = ( DTOLD_I * UOLDDGI + DTOLD_J * UOLDDGI2 ) / (DTOLD_I + DTOLD_J)
             VOLDDGI_INT = ( DTOLD_I * VOLDDGI + DTOLD_J * VOLDDGI2 ) / (DTOLD_I + DTOLD_J)
             WOLDDGI_INT = ( DTOLD_I * WOLDDGI + DTOLD_J * WOLDDGI2 ) / (DTOLD_I + DTOLD_J)
          ENDIF

          UDGI  = UDGI_INT 
          VDGI  = VDGI_INT 
          WDGI  = WDGI_INT

          UGI_COEF_ELE(:)=DT_I * UGI_COEF_ELE(:) / (DT_I + DT_J)
          VGI_COEF_ELE(:)=DT_I * VGI_COEF_ELE(:) / (DT_I + DT_J)
          WGI_COEF_ELE(:)=DT_I * WGI_COEF_ELE(:) / (DT_I + DT_J)

          UGI_COEF_ELE2(:)=DT_J * UGI_COEF_ELE2(:) / (DT_I + DT_J)
          VGI_COEF_ELE2(:)=DT_J * VGI_COEF_ELE2(:) / (DT_I + DT_J)
          WGI_COEF_ELE2(:)=DT_J * WGI_COEF_ELE2(:) / (DT_I + DT_J)

          UOLDDGI  = UOLDDGI_INT 
          VOLDDGI  = VOLDDGI_INT 
          WOLDDGI  = WOLDDGI_INT

       ENDIF Conditional_ELE2

    ENDIF Conditional_SELE

    NDOTQ =  CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI + CVNORMZ(GI) * WDGI
    NDOTQOLD =  CVNORMX( GI ) * UOLDDGI + CVNORMY( GI ) * VOLDDGI + CVNORMZ(GI) * WOLDDGI

    ! Define whether flux is incoming or outgoing, depending on direction of flow
    IF( NDOTQ >  0.0 ) THEN
       INCOME = 0.0  !Outgoing
    ELSE
       INCOME = 1.0  !Incoming
    ENDIF
    IF( NDOTQOLD >  0.0 ) THEN
       INCOMEOLD = 0.0  !Outgoing
    ELSE
       INCOMEOLD = 1.0  !Incoming
    ENDIF

    RETURN  

  END SUBROUTINE GET_INT_VEL_ORIG




  SUBROUTINE GET_INT_T_DEN( FVT, FVTOLD, FVT2, FVT2OLD, FVD, FVDOLD, LIMD, LIMT, LIMT2, &
       LIMDOLD, LIMTOLD, LIMT2OLD, LIMDT, LIMDTOLD, LIMDTT2, LIMDTT2OLD, &
       FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
       CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
       CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, INCOMEOLD, &
       T, TOLD, T2, T2OLD, DEN, DENOLD, FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, &
       TMIN, TOLDMIN, T2MIN, T2OLDMIN, DENMIN, DENOLDMIN, &
       TMAX, TOLDMAX, T2MAX, T2OLDMAX, DENMAX, DENOLDMAX, &
       SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
       WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
       WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, VOLFRA_PORE, &
       MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
       DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD, &
       T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2 )
    !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
    IMPLICIT NONE
    ! Calculate T and DEN on the CV face at quadrature point GI.
    REAL, intent( inout ) :: FVT,FVTOLD, FVT2, FVT2OLD,FVD,FVDOLD, LIMD,LIMT,LIMT2, &
         LIMDOLD,LIMTOLD,LIMT2OLD,LIMDT,LIMDTOLD,LIMDTT2,LIMDTT2OLD, &
         FEMDGI, FEMTGI, FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI
    REAL, intent( in ) :: INCOME,INCOMEOLD
    INTEGER, intent( in ) :: CV_DISOPT,CV_NONODS,NPHASE,CV_NODI_IPHA,CV_NODJ_IPHA,ELE,ELE2,  &
         CV_NLOC,TOTELE,SCVNGI,GI,IPHASE,SELE,CV_SNLOC,STOTEL, &
         WIC_T_BC_DIRICHLET,WIC_D_BC_DIRICHLET, IGOT_T2
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( CV_NLOC ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( CV_SNLOC ), intent( in ) :: CV_SLOC2LOC
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_T_BC, WIC_D_BC
    INTEGER, DIMENSION( STOTEL * NPHASE * IGOT_T2 ), intent( in ) :: WIC_T2_BC
    REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFEN
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: T, TOLD, DEN, DENOLD, FEMT, FEMTOLD,  &
         FEMDEN, FEMDENOLD
    REAL, DIMENSION( CV_NONODS * NPHASE  * IGOT_T2 ), intent( in ) :: T2, T2OLD, FEMT2, FEMT2OLD
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN, TOLDMIN, DENMIN, DENOLDMIN, TMAX, TOLDMAX, DENMAX, DENOLDMAX
    REAL, DIMENSION( CV_NONODS * NPHASE  * IGOT_T2 ), intent( inout ) :: T2MIN, T2OLDMIN, T2MAX, T2OLDMAX
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: MASS_CV
    INTEGER, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
         TOLDMAX_NOD, DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD
    INTEGER, DIMENSION( CV_NONODS * NPHASE  * IGOT_T2), intent( in ) :: &
         T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD
    REAL, DIMENSION( CV_SNLOC*STOTEL*NPHASE  ), intent( in ) :: SUF_T_BC, SUF_D_BC
    REAL, DIMENSION( CV_SNLOC*STOTEL*NPHASE * IGOT_T2 ), intent( in ) :: SUF_T2_BC
    REAL, DIMENSION( TOTELE  ), intent( in ) :: VOLFRA_PORE
    ! Local variables
    ! If UPWIND then use upwind flux between elements else use central. 
    ! If HI_ORDER_HALF then use high order interpolation when around 
    ! a volume frac of 0.5 and gradually apply limiting near 0 and 1. 
    LOGICAL, PARAMETER :: UPWIND = .TRUE., HI_ORDER_HALF = .FALSE., LIM_VOL_ADJUST = .TRUE.
    LOGICAL :: FIRSTORD, NOLIMI, RESET_STORE
    REAL :: RELAX, RELAXOLD, TMIN_STORE, TMAX_STORE, TOLDMIN_STORE, TOLDMAX_STORE, &
         T2MIN_STORE, T2MAX_STORE, T2OLDMIN_STORE, T2OLDMAX_STORE, &
         DENMIN_STORE, DENMAX_STORE, DENOLDMIN_STORE, DENOLDMAX_STORE
    INTEGER :: CV_KLOC, CV_NODK, CV_NODK_IPHA, CV_KLOC2, CV_NODK2, CV_NODK2_IPHA, CV_STAR_IPHA, &
         CV_SKLOC, CV_SNODK, CV_SNODK_IPHA

    FVT2    = 1.0
    FVT2OLD = 1.0

    IF(SELE == 0) THEN ! Is NOT on boundary of the domain

       FVT    = INCOME * T( CV_NODJ_IPHA ) + ( 1. - INCOME ) * T( CV_NODI_IPHA )
       FVTOLD = INCOMEOLD * TOLD( CV_NODJ_IPHA ) + ( 1. - INCOMEOLD ) * TOLD( CV_NODI_IPHA )

       FVD    = INCOME * DEN( CV_NODJ_IPHA ) + ( 1. - INCOME ) * DEN( CV_NODI_IPHA )
       FVDOLD = INCOMEOLD * DENOLD( CV_NODJ_IPHA ) +  (1. - INCOMEOLD ) * DENOLD( CV_NODI_IPHA )

       IF(IGOT_T2==1) THEN
          FVT2    = INCOME * T2( CV_NODJ_IPHA ) + ( 1. - INCOME ) * T2( CV_NODI_IPHA )
          FVT2OLD = INCOMEOLD * T2OLD( CV_NODJ_IPHA ) + ( 1. - INCOMEOLD ) * T2OLD( CV_NODI_IPHA )
       ENDIF

    ELSE ! Is on boundary of the domain
       IF( WIC_T_BC( SELE + ( IPHASE - 1 ) * STOTEL ) /= WIC_T_BC_DIRICHLET ) THEN 
          ! Dont apply a Dirichlet b.c.
          FVT    = T( CV_NODI_IPHA )
          FVTOLD = TOLD( CV_NODI_IPHA )

          FVD    = DEN( CV_NODI_IPHA )
          FVDOLD = DENOLD( CV_NODI_IPHA )

          IF(IGOT_T2==1) THEN
             FVT2    = T2( CV_NODI_IPHA ) 
             FVT2OLD = T2OLD( CV_NODI_IPHA ) 
          ENDIF
       ELSE
          DO CV_SKLOC = 1, CV_SNLOC
             CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
             CV_NODK_IPHA = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC ) + ( IPHASE - 1 ) * CV_NONODS
             IF(CV_NODK_IPHA == CV_NODI_IPHA) THEN
                CV_SNODK = ( SELE - 1 ) * CV_SNLOC + CV_SKLOC 
                CV_SNODK_IPHA = CV_SNODK + ( IPHASE - 1 ) * STOTEL*CV_SNLOC
             ENDIF
          END DO
          FVT    = INCOME * SUF_T_BC(CV_SNODK_IPHA) + ( 1. - INCOME ) * T( CV_NODI_IPHA )
          FVTOLD = INCOMEOLD * SUF_T_BC(CV_SNODK_IPHA) + ( 1. - INCOMEOLD ) * TOLD( CV_NODI_IPHA )

          FVD    = INCOME * SUF_D_BC(CV_SNODK_IPHA) + ( 1. - INCOME ) * DEN( CV_NODI_IPHA )
          FVDOLD = INCOMEOLD * SUF_D_BC(CV_SNODK_IPHA) +  (1. - INCOMEOLD ) * DENOLD( CV_NODI_IPHA )

          IF(IGOT_T2==1) THEN
             FVT2    = INCOME * SUF_T2_BC(CV_SNODK_IPHA) + ( 1. - INCOME ) * T2( CV_NODI_IPHA )
             FVT2OLD = INCOMEOLD * SUF_T2_BC(CV_SNODK_IPHA) + ( 1. - INCOMEOLD ) * T2OLD( CV_NODI_IPHA )
          ENDIF

       ENDIF
       IF( WIC_D_BC( SELE + ( IPHASE - 1 ) * STOTEL ) /= WIC_D_BC_DIRICHLET ) THEN 
          FVD    = DEN( CV_NODI_IPHA )
          FVDOLD = DENOLD( CV_NODI_IPHA )
       ENDIF
    ENDIF

    ! By default do not use first-order upwinding
    FIRSTORD = .FALSE.

    ! No limiting if CV_DISOPT is 6 or 7  (why not just define limt=femt and skip to assembly?)
    NOLIMI = ( INT( CV_DISOPT / 2 ) == 3 ) 

    ! Make a guess at the CV face value of advected field variable and density 
    ! (Depends on discetisation option, CV_DISOPT)                
    SELECT CASE( CV_DISOPT / 2 )

    CASE( 0 ) ! First-order upwinding
       FEMTGI    = FVT
       FEMTOLDGI = FVTOLD

       FEMT2GI    = FVT2
       FEMT2OLDGI = FVT2OLD

       FEMDGI    = FVD
       FEMDOLDGI = FVDOLD
       FIRSTORD = .TRUE.  

    CASE( 1 ) ! Central differencing [Trapezoidal rule (2 OR 3)]
       FEMTGI    = 0.5 * ( T( CV_NODI_IPHA ) + T( CV_NODJ_IPHA ))   ! need to fix this!
       FEMTOLDGI = 0.5 * ( TOLD( CV_NODI_IPHA ) + TOLD( CV_NODJ_IPHA )) 
       IF(IGOT_T2==1) THEN
          FEMT2GI   = 0.5 * ( T2( CV_NODI_IPHA ) + T2( CV_NODJ_IPHA ))   ! need to fix this!
          FEMT2OLDGI= 0.5 * ( T2OLD( CV_NODI_IPHA ) + T2OLD( CV_NODJ_IPHA )) 
       ELSE
          FEMT2GI    = 1.0
          FEMT2OLDGI = 1.0
       ENDIF
       FEMDGI    = 0.5 * ( DEN( CV_NODI_IPHA ) + DEN( CV_NODJ_IPHA ))   ! need to fix this!
       FEMDOLDGI = 0.5 * ( DENOLD( CV_NODI_IPHA ) + DENOLD( CV_NODJ_IPHA )) 

    CASE DEFAULT ! Finite element approximation (4 OR 5)(6 or 7)(8 or 9)
       FEMTGI    = 0.0
       FEMTOLDGI = 0.0
       FEMT2GI   = 0.0
       FEMT2OLDGI= 0.0
       FEMDGI    = 0.0
       FEMDOLDGI = 0.0

       Conditional_CV_DISOPT_ELE2: IF(SELE /= 0) THEN
          ! Is on boundary of the domain

          IF(WIC_T_BC(SELE+(IPHASE-1)*STOTEL) /= WIC_T_BC_DIRICHLET) THEN ! Dont apply a Dirichlet bc
             DO CV_KLOC = 1, CV_NLOC
                CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                FEMTGI    = FEMTGI     +  SCVFEN( CV_KLOC, GI ) * FEMT( CV_NODK_IPHA )
                FEMTOLDGI = FEMTOLDGI  +  SCVFEN( CV_KLOC, GI ) * FEMTOLD( CV_NODK_IPHA )
                IF(IGOT_T2==1) THEN
                   FEMT2GI    = FEMT2GI     +  SCVFEN( CV_KLOC, GI ) * FEMT2( CV_NODK_IPHA )
                   FEMT2OLDGI = FEMT2OLDGI  +  SCVFEN( CV_KLOC, GI ) * FEMT2OLD( CV_NODK_IPHA )
                ELSE
                   FEMT2GI    = 1.0
                   FEMT2OLDGI = 1.0
                ENDIF
             END DO
          ELSE
             DO CV_SKLOC = 1, CV_SNLOC
                CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
                CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                CV_SNODK = ( SELE - 1 ) * CV_SNLOC + CV_SKLOC 
                CV_SNODK_IPHA = CV_SNODK + ( IPHASE - 1 ) * STOTEL*CV_SNLOC
                FEMTGI = FEMTGI +  SCVFEN( CV_KLOC, GI ) * ( SUF_T_BC( CV_SNODK_IPHA ) & 
                     * INCOME + FEMT( CV_NODK_IPHA ) * ( 1. -INCOME ))
                FEMTOLDGI = FEMTOLDGI + SCVFEN( CV_KLOC, GI ) * ( SUF_T_BC( CV_SNODK_IPHA ) &
                     * INCOMEOLD + FEMTOLD( CV_NODK_IPHA ) * ( 1. - INCOMEOLD ))
                IF(IGOT_T2==1) THEN
                   FEMT2GI = FEMT2GI +  SCVFEN( CV_KLOC, GI ) * ( SUF_T2_BC( CV_SNODK_IPHA ) & 
                        * INCOME + FEMT2( CV_NODK_IPHA ) * ( 1. -INCOME ))
                   FEMT2OLDGI = FEMT2OLDGI + SCVFEN( CV_KLOC, GI ) * ( SUF_T2_BC( CV_SNODK_IPHA ) &
                        * INCOMEOLD + FEMT2OLD( CV_NODK_IPHA ) * ( 1. - INCOMEOLD ))
                ELSE
                   FEMT2GI    = 1.0
                   FEMT2OLDGI = 1.0
                ENDIF
             END DO
          ENDIF

          IF(WIC_D_BC(SELE+(IPHASE-1)*STOTEL) /= WIC_D_BC_DIRICHLET) THEN ! Dont apply a Dirichlet bc
             DO CV_KLOC = 1, CV_NLOC
                CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                FEMDGI    = FEMDGI     +  SCVFEN( CV_KLOC, GI ) * FEMDEN( CV_NODK_IPHA )
                FEMDOLDGI = FEMDOLDGI  +  SCVFEN( CV_KLOC, GI ) * FEMDENOLD( CV_NODK_IPHA )
             END DO
          ELSE
             DO CV_SKLOC = 1, CV_SNLOC
                CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
                CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                CV_SNODK = ( SELE - 1 ) * CV_SNLOC + CV_SKLOC 
                CV_SNODK_IPHA = CV_SNODK + ( IPHASE - 1 ) * STOTEL*CV_SNLOC
                FEMDGI = FEMDGI + SCVFEN( CV_KLOC, GI ) * ( SUF_D_BC( CV_SNODK_IPHA ) &
                     * INCOME + FEMDEN( CV_NODK_IPHA ) * ( 1. - INCOME ))
                FEMDOLDGI = FEMDOLDGI + SCVFEN( CV_KLOC, GI ) * ( SUF_D_BC( CV_SNODK_IPHA ) &
                     * INCOMEOLD + FEMDENOLD( CV_NODK_IPHA ) * ( 1. - INCOMEOLD ))
             END DO
          ENDIF
       ELSE IF( ( ELE2 == 0 ).OR.( ELE2 == ELE ) ) THEN

          DO CV_KLOC = 1, CV_NLOC
             CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
             CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
             FEMTGI    = FEMTGI     +  SCVFEN( CV_KLOC, GI ) * FEMT( CV_NODK_IPHA )
             FEMTOLDGI = FEMTOLDGI  +  SCVFEN( CV_KLOC, GI ) * FEMTOLD( CV_NODK_IPHA )
             FEMDGI    = FEMDGI     +  SCVFEN( CV_KLOC, GI ) * FEMDEN( CV_NODK_IPHA )
             FEMDOLDGI = FEMDOLDGI  +  SCVFEN( CV_KLOC, GI ) * FEMDENOLD( CV_NODK_IPHA )
             IF(IGOT_T2==1) THEN
                FEMT2GI    = FEMT2GI     +  SCVFEN( CV_KLOC, GI ) * FEMT2( CV_NODK_IPHA )
                FEMT2OLDGI = FEMT2OLDGI  +  SCVFEN( CV_KLOC, GI ) * FEMT2OLD( CV_NODK_IPHA )
             ELSE
                FEMT2GI    = 1.0
                FEMT2OLDGI = 1.0
             ENDIF
          END DO

       ELSE  ! DG saturation across elements

          IF(UPWIND) THEN
             DO CV_KLOC = 1, CV_NLOC
                CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                IF(CV_KLOC2 /= 0 ) THEN
                   CV_NODK2 = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 )
                   CV_NODK2_IPHA = CV_NODK2 + ( IPHASE - 1 ) * CV_NONODS
                   CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                   CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                   FEMTGI = FEMTGI +  SCVFEN( CV_KLOC, GI ) * ( FEMT( CV_NODK2_IPHA ) & 
                        * INCOME + FEMT( CV_NODK_IPHA ) * ( 1. -INCOME ))
                   FEMTOLDGI = FEMTOLDGI + SCVFEN( CV_KLOC, GI ) * ( FEMTOLD( CV_NODK2_IPHA ) &
                        * INCOMEOLD + FEMTOLD( CV_NODK_IPHA ) * ( 1. - INCOMEOLD ))
                   FEMDGI = FEMDGI + SCVFEN( CV_KLOC, GI ) * ( FEMDEN( CV_NODK2_IPHA ) &
                        * INCOME + FEMDEN( CV_NODK_IPHA ) * ( 1. - INCOME ))
                   FEMDOLDGI = FEMDOLDGI + SCVFEN( CV_KLOC, GI ) * ( FEMDENOLD( CV_NODK2_IPHA ) &
                        * INCOMEOLD + FEMDENOLD( CV_NODK_IPHA ) * ( 1. - INCOMEOLD ))
                   IF(IGOT_T2==1) THEN
                      FEMT2GI = FEMT2GI +  SCVFEN( CV_KLOC, GI ) * ( FEMT2( CV_NODK2_IPHA ) & 
                           * INCOME + FEMT2( CV_NODK_IPHA ) * ( 1. -INCOME ))
                      FEMT2OLDGI = FEMT2OLDGI + SCVFEN( CV_KLOC, GI ) * ( FEMT2OLD( CV_NODK2_IPHA ) &
                           * INCOMEOLD + FEMT2OLD( CV_NODK_IPHA ) * ( 1. - INCOMEOLD ))
                   ELSE
                      FEMT2GI    = 1.0
                      FEMT2OLDGI = 1.0
                   ENDIF
                ENDIF
             END DO
          ELSE
             DO CV_KLOC = 1, CV_NLOC
                CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                IF(CV_KLOC2 /= 0 ) THEN
                   CV_NODK2 = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 )
                   CV_NODK2_IPHA = CV_NODK2 + ( IPHASE - 1 ) * CV_NONODS
                   CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                   CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                   FEMTGI = FEMTGI +  SCVFEN( CV_KLOC, GI ) * ( FEMT( CV_NODK2_IPHA ) & 
                        * 0.5 + FEMT( CV_NODK_IPHA ) * 0.5 )
                   FEMTOLDGI = FEMTOLDGI + SCVFEN( CV_KLOC, GI ) * ( FEMTOLD( CV_NODK2_IPHA ) &
                        * 0.5 + FEMTOLD( CV_NODK_IPHA ) * 0.5 )
                   FEMDGI = FEMDGI + SCVFEN( CV_KLOC, GI ) * ( FEMDEN( CV_NODK2_IPHA ) &
                        * 0.5 + FEMDEN( CV_NODK_IPHA ) * 0.5 )
                   FEMDOLDGI = FEMDOLDGI + SCVFEN( CV_KLOC, GI ) * ( FEMDENOLD( CV_NODK2_IPHA ) &
                        * 0.5 + FEMDENOLD( CV_NODK_IPHA ) * 0.5 )
                   IF(IGOT_T2==1) THEN
                      FEMT2GI = FEMT2GI +  SCVFEN( CV_KLOC, GI ) * ( FEMT2( CV_NODK2_IPHA ) & 
                           * 0.5 + FEMT2( CV_NODK_IPHA ) * 0.5 )
                      FEMT2OLDGI = FEMT2OLDGI + SCVFEN( CV_KLOC, GI ) * ( FEMT2OLD( CV_NODK2_IPHA ) &
                           * 0.5 + FEMT2OLD( CV_NODK_IPHA ) * 0.5 )
                   ELSE
                      FEMT2GI    = 1.0
                      FEMT2OLDGI = 1.0
                   ENDIF
                ENDIF
             END DO
          ENDIF

       ENDIF Conditional_CV_DISOPT_ELE2

    END SELECT

    ! The original isotropic limiting method for all other cases
    ! Call the original NVD limiting routine to limit the face value of the old 
    ! and new field variable and "density"
    CV_STAR_IPHA = 1 + ( IPHASE - 1 ) * CV_NONODS

    IF(CV_NODI_IPHA == CV_NODJ_IPHA) THEN
       ! On the boundary of the domain and thus use the 1st order flux. 
       LIMT   =FVT
       LIMTOLD=FVTOLD
       LIMD   =FVD
       LIMDOLD=FVDOLD
       LIMT2   =FVT2
       LIMT2OLD=FVT2OLD
    ELSE

       IF(LIM_VOL_ADJUST) THEN
          RESET_STORE=.FALSE. 
          CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
          CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

          CALL CAL_LIM_VOL_ADJUST(TOLDMIN_STORE,TOLDMIN,TOLD,TOLDMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
          CALL CAL_LIM_VOL_ADJUST(TOLDMAX_STORE,TOLDMAX,TOLD,TOLDMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )

          CALL CAL_LIM_VOL_ADJUST(DENMIN_STORE,DENMIN,DEN,DENMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
          CALL CAL_LIM_VOL_ADJUST(DENMAX_STORE,DENMAX,DEN,DENMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

          CALL CAL_LIM_VOL_ADJUST(DENOLDMIN_STORE,DENOLDMIN,DENOLD,DENOLDMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
          CALL CAL_LIM_VOL_ADJUST(DENOLDMAX_STORE,DENOLDMAX,DENOLD,DENOLDMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )

          IF(IGOT_T2==1) THEN
             CALL CAL_LIM_VOL_ADJUST(T2MIN_STORE,T2MIN,T2,T2MIN_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
             CALL CAL_LIM_VOL_ADJUST(T2MAX_STORE,T2MAX,T2,T2MAX_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

             CALL CAL_LIM_VOL_ADJUST(T2OLDMIN_STORE,T2OLDMIN,T2OLD,T2OLDMIN_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
             CALL CAL_LIM_VOL_ADJUST(T2OLDMAX_STORE,T2OLDMAX,T2OLD,T2OLDMAX_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
          ENDIF
       ENDIF

       CALL ONVDLIM( CV_NONODS, &
            LIMT, FEMTGI, INCOME, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
            T( CV_STAR_IPHA ), TMIN( CV_STAR_IPHA ), TMAX( CV_STAR_IPHA ), FIRSTORD, NOLIMI )

       CALL ONVDLIM( CV_NONODS, &
            LIMTOLD, FEMTOLDGI, INCOMEOLD, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
            TOLD( CV_STAR_IPHA ), TOLDMIN( CV_STAR_IPHA ), TOLDMAX( CV_STAR_IPHA ), FIRSTORD, NOLIMI )

       CALL ONVDLIM( CV_NONODS, &
            LIMD, FEMDGI, INCOME, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
            DEN( CV_STAR_IPHA ), DENMIN( CV_STAR_IPHA ), DENMAX( CV_STAR_IPHA ), FIRSTORD, NOLIMI )

       CALL ONVDLIM( CV_NONODS, &
            LIMDOLD,FEMDOLDGI, INCOMEOLD, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
            DENOLD( CV_STAR_IPHA ), DENOLDMIN( CV_STAR_IPHA ), DENOLDMAX( CV_STAR_IPHA ), FIRSTORD, NOLIMI )

       IF(IGOT_T2==1) THEN
          CALL ONVDLIM( CV_NONODS, &
               LIMT2, FEMT2GI, INCOME, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
               T2( CV_STAR_IPHA ), T2MIN( CV_STAR_IPHA ), T2MAX( CV_STAR_IPHA ), FIRSTORD, NOLIMI )

          CALL ONVDLIM( CV_NONODS, &
               LIMT2OLD, FEMT2OLDGI, INCOMEOLD, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
               T2OLD( CV_STAR_IPHA ), T2OLDMIN( CV_STAR_IPHA ), T2OLDMAX( CV_STAR_IPHA ), FIRSTORD, NOLIMI )
       ELSE
          LIMT2   =1.0
          LIMT2OLD=1.0
       ENDIF

       IF(LIM_VOL_ADJUST) THEN
          RESET_STORE=.TRUE. 
          CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
          CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

          CALL CAL_LIM_VOL_ADJUST(TOLDMIN_STORE,TOLDMIN,TOLD,TOLDMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
          CALL CAL_LIM_VOL_ADJUST(TOLDMAX_STORE,TOLDMAX,TOLD,TOLDMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )

          CALL CAL_LIM_VOL_ADJUST(DENMIN_STORE,DENMIN,DEN,DENMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
          CALL CAL_LIM_VOL_ADJUST(DENMAX_STORE,DENMAX,DEN,DENMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

          CALL CAL_LIM_VOL_ADJUST(DENOLDMIN_STORE,DENOLDMIN,DENOLD,DENOLDMIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
          CALL CAL_LIM_VOL_ADJUST(DENOLDMAX_STORE,DENOLDMAX,DENOLD,DENOLDMAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )

          IF(IGOT_T2==1) THEN
             CALL CAL_LIM_VOL_ADJUST(T2MIN_STORE,T2MIN,T2,T2MIN_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
             CALL CAL_LIM_VOL_ADJUST(T2MAX_STORE,T2MAX,T2,T2MAX_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

             CALL CAL_LIM_VOL_ADJUST(T2OLDMIN_STORE,T2OLDMIN,T2OLD,T2OLDMIN_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
             CALL CAL_LIM_VOL_ADJUST(T2OLDMAX_STORE,T2OLDMAX,T2OLD,T2OLDMAX_NOD,RESET_STORE,MASS_CV, &
                  CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOMEOLD )
          ENDIF
       ENDIF

    ENDIF

    ! Amend for porocity...    
    !       if(.false.) then
    IF(ELE2 /=0) THEN 
       FVD    = 0.5*(VOLFRA_PORE(ELE)+VOLFRA_PORE(ELE2))*FVD
       FVDOLD = 0.5*(VOLFRA_PORE(ELE)+VOLFRA_PORE(ELE2))*FVDOLD
       LIMD   = 0.5*(VOLFRA_PORE(ELE)+VOLFRA_PORE(ELE2))*LIMD
       LIMDOLD= 0.5*(VOLFRA_PORE(ELE)+VOLFRA_PORE(ELE2))*LIMDOLD
    ELSE
       FVD    = VOLFRA_PORE(ELE)*FVD
       FVDOLD = VOLFRA_PORE(ELE)*FVDOLD
       LIMD   = VOLFRA_PORE(ELE)*LIMD
       LIMDOLD= VOLFRA_PORE(ELE)*LIMDOLD
    ENDIF
    !      endif

    IF(HI_ORDER_HALF) THEN
       RELAX=MIN(2.*ABS(FEMTGI-0.5),1.0)
       !      relax=0.
       LIMT=RELAX*LIMT+(1.-RELAX)*FEMTGI
       RELAXOLD=MIN(2.*ABS(FEMTOLDGI-0.5),1.0)
       !      RELAXOLD=0.
       LIMTOLD=RELAXOLD*LIMTOLD+(1.-RELAXOLD)*FEMTOLDGI
    ENDIF

    LIMDT = LIMD * LIMT
    LIMDTOLD = LIMDOLD * LIMTOLD

    LIMDTT2 = LIMD * LIMT * LIMT2
    LIMDTT2OLD = LIMDOLD * LIMTOLD * LIMT2OLD

    RETURN

  END SUBROUTINE GET_INT_T_DEN



  SUBROUTINE CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
       CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
    ! Adjust TMIN to take into account different sized CV's. 
    ! if RESET_STORE then reset TMIN to orginal values.
    REAL, intent( in ) :: INCOME
    INTEGER, intent( in ) :: CV_NODI_IPHA, CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE
    LOGICAL, intent( in ) :: RESET_STORE
    REAL, intent( inout ) :: TMIN_STORE
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN
    REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: T
    INTEGER, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: TMIN_NOD
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: MASS_CV
    ! Local variables...
    REAL DX1,   DX2_MIN, COEFF
    INTEGER CV_NODI, CV_NODJ
    LOGICAL, PARAMETER :: DEF_BOUNDED = .false.

    CV_NODI = CV_NODI_IPHA-(IPHASE-1)*CV_NONODS
    CV_NODJ = CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS

    IF(RESET_STORE) THEN

       IF(INCOME<0.5) THEN
          TMIN(CV_NODI_IPHA)=TMIN_STORE
       ELSE
          TMIN(CV_NODJ_IPHA)=TMIN_STORE
       ENDIF

    ELSE

       DX1=0.5*(MASS_CV(CV_NODI) + MASS_CV(CV_NODJ))

       IF(INCOME<0.5) THEN
          DX2_MIN=0.5*(MASS_CV( TMIN_NOD(CV_NODI_IPHA)-(IPHASE-1)*CV_NONODS) &
               +MASS_CV( CV_NODI)) 
          TMIN_STORE=TMIN(CV_NODI_IPHA)
          IF(DEF_BOUNDED) THEN ! This produces strictly bounded always soln
             COEFF=Min(1.0,(DX1/DX2_MIN)) 
          ELSE
             COEFF=(DX1/DX2_MIN)
          ENDIF
          TMIN(CV_NODI_IPHA)=T(CV_NODI_IPHA) + COEFF*(TMIN_STORE-T(CV_NODI_IPHA))
       ELSE
          DX2_MIN=0.5*(MASS_CV( TMIN_NOD(CV_NODJ_IPHA)-(IPHASE-1)*CV_NONODS) &
               +MASS_CV( CV_NODJ)) 
          TMIN_STORE=TMIN(CV_NODJ_IPHA)
          IF(DEF_BOUNDED) THEN ! This produces strictly bounded always soln
             COEFF=Min(1.0,(DX1/DX2_MIN)) 
          ELSE
             COEFF=(DX1/DX2_MIN)
          ENDIF
          TMIN(CV_NODJ_IPHA)=T(CV_NODJ_IPHA) + COEFF*(TMIN_STORE-T(CV_NODJ_IPHA))
       ENDIF

    ENDIF

    RETURN
  END SUBROUTINE CAL_LIM_VOL_ADJUST





  REAL FUNCTION FACE_THETA( DT, CV_THETA, HDC, NDOTQ, LIMDT, DIFF_COEF_DIVDX, &
       T_NODJ_IPHA, T_NODI_IPHA,  &
       NDOTQOLD, LIMDTOLD, DIFF_COEFOLD_DIVDX, TOLD_NODJ_IPHA, TOLD_NODI_IPHA )
    IMPLICIT NONE
    ! Define face value of theta
    REAL :: DT, CV_THETA, HDC, NDOTQ, LIMDT, DIFF_COEF_DIVDX, T_NODJ_IPHA, T_NODI_IPHA,  &
         NDOTQOLD, LIMDTOLD, DIFF_COEFOLD_DIVDX, TOLD_NODJ_IPHA, TOLD_NODI_IPHA
    ! Local variables
    REAL :: FTHETA, HF, HFOLD, GF, PINVTH, QINVTH

    IF( CV_THETA >= 0.0) THEN ! Specified
       FTHETA = CV_THETA
    ELSE ! Non-linear
       HF    = NDOTQ * LIMDT + DIFF_COEF_DIVDX * ( T_NODI_IPHA - T_NODJ_IPHA )
       HFOLD = NDOTQOLD * LIMDTOLD + DIFF_COEFOLD_DIVDX * ( TOLD_NODI_IPHA - TOLD_NODJ_IPHA )
       GF = TOLFUN( DT *  (HF - HFOLD) )
       PINVTH = HDC * ( T_NODI_IPHA - TOLD_NODI_IPHA ) / GF
       QINVTH = HDC * ( T_NODJ_IPHA - TOLD_NODJ_IPHA ) / GF
       ! 0.5 is the original value. 
       !       FTHETA = MAX( 0.5, 1. - 0.5 * MIN( ABS( PINVTH ), ABS( QINVTH ))) 
       !       FTHETA = MAX( 0.5, 1. - 0.25 * MIN( ABS( PINVTH ), ABS( QINVTH ))) 
       FTHETA = MAX( 0.5, 1. - 0.125 * MIN( ABS( PINVTH ), ABS( QINVTH ))) 
    ENDIF

    FACE_THETA = FTHETA

    RETURN

  END FUNCTION FACE_THETA




  SUBROUTINE DGSDETNXLOC2( SNLOC, SNGI, &
       XSL, YSL, ZSL, &
       SN, SNLX, SNLY, SWEIGH, SDETWE, SAREA, &
       D1, D3, DCYL, &
       NORMXN, NORMYN, NORMZN, &
       NORMX, NORMY, NORMZ )
    IMPLICIT NONE

    INTEGER, intent( in ) :: SNLOC, SNGI
    REAL, DIMENSION( SNLOC ), intent( in ) :: XSL, YSL, ZSL
    REAL, DIMENSION( SNLOC, SNGI ), intent( in ) :: SN, SNLX, SNLY
    REAL, DIMENSION( SNGI ), intent( in ) :: SWEIGH
    REAL, DIMENSION( SNGI ), intent( inout ) :: SDETWE 
    REAL, intent( inout ) ::  SAREA
    LOGICAL, intent( in ) ::  D1,D3,DCYL
    REAL, DIMENSION( SNGI ), intent( inout ) :: NORMXN, NORMYN ,NORMZN
    REAL, intent( inout ) :: NORMX, NORMY, NORMZ
    ! Local variables
    real, parameter :: pi = 3.141592654
    INTEGER :: GI, SL, IGLX
    REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY
    REAL :: A, B, C, DETJ, RGI, TWOPI

    SAREA=0.

    IF(D3) THEN
       DO GI=1,SNGI

          DXDLX=0.
          DXDLY=0.
          DYDLX=0.
          DYDLY=0.
          DZDLX=0.
          DZDLY=0.

          DO SL=1,SNLOC
             DXDLX=DXDLX + SNLX(SL,GI)*XSL(SL)
             DXDLY=DXDLY + SNLY(SL,GI)*XSL(SL)
             DYDLX=DYDLX + SNLX(SL,GI)*YSL(SL)
             DYDLY=DYDLY + SNLY(SL,GI)*YSL(SL)
             DZDLX=DZDLX + SNLX(SL,GI)*ZSL(SL)
             DZDLY=DZDLY + SNLY(SL,GI)*ZSL(SL)
          END DO
          A = DYDLX*DZDLY - DYDLY*DZDLX
          B = DXDLX*DZDLY - DXDLY*DZDLX
          C = DXDLX*DYDLY - DXDLY*DYDLX

          DETJ=SQRT( A**2 + B**2 + C**2)
          SDETWE(GI)=DETJ*SWEIGH(GI)
          SAREA=SAREA+SDETWE(GI)

          ! Calculate the normal at the Gauss pts...
          ! Perform x-product. N=T1 x T2
          CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
               DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
               NORMX,NORMY,NORMZ)

       END DO
       ! IF(D3) THEN...
    ELSE IF(.NOT.D1) THEN
       ! ENDOF IF(D3) THEN ELSE...
       TWOPI=1.0
       IF(DCYL) TWOPI=2.*PI

       DO GI=1,SNGI
          RGI=0.
          DXDLX=0.
          DXDLY=0.
          DYDLX=0.
          DYDLY=0.
          DZDLX=0.
          ! DZDLY=1 is to calculate the normal.
          DZDLY=1.
          DO SL=1,SNLOC
             DXDLX=DXDLX + SNLX(SL,GI)*XSL(SL)
             DYDLX=DYDLX + SNLX(SL,GI)*YSL(SL)
             RGI=RGI+SN(SL,GI)*YSL(IGLX)
          END DO
          IF(.NOT.DCYL) RGI=1.0
          DETJ=SQRT( DXDLX**2 + DYDLX**2 )
          SDETWE(GI)=TWOPI*RGI*DETJ*SWEIGH(GI)
          SAREA=SAREA+SDETWE(GI)
          NORMZ=0.
          CALL NORMGI(NORMXN(GI),NORMYN(GI),NORMZN(GI), &
               DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
               NORMX,NORMY,NORMZ)
       END DO
    ELSE ! For 1D...
       DO GI = 1, SNGI
          DXDLX = 0.
          DO SL = 1, SNLOC
             DXDLX = DXDLX + SNLX( SL, GI ) * XSL( SL )
          END DO
          SDETWE( GI ) = SWEIGH( GI )
          SAREA = SAREA + SDETWE( GI )
          NORMXN( GI ) = NORMX
          NORMYN( GI ) = 0.0
          NORMZN( GI ) = 0.0
       END DO
    ENDIF

  END SUBROUTINE DGSDETNXLOC2




  SUBROUTINE DGSIMPLNORM( ELE, SILOC2ILOC, TOTELE, NLOC, SNLOC, XONDGL, &
       X, Y, Z, XNONOD, NORMX, NORMY, NORMZ )
    ! Form approximate surface normal (NORMX,NORMY,NORMZ)
    IMPLICIT NONE
    INTEGER, intent( in ) :: ELE, TOTELE, NLOC, SNLOC, XNONOD
    INTEGER, DIMENSION( SNLOC ), intent( in ) ::  SILOC2ILOC
    INTEGER, DIMENSION( TOTELE * NLOC ), intent( in ) :: XONDGL
    REAL, DIMENSION( XNONOD ), intent( in ) :: X, Y, Z
    REAL, intent( inout ) :: NORMX, NORMY, NORMZ
    ! Local variables
    REAL :: XC, YC, ZC, SXC, SYC, SZC, NORM
    INTEGER :: XNODI, ILOC, SILOC

    XC = 0.0
    YC = 0.0
    ZC = 0.0
    DO ILOC = 1, NLOC
       XNODI = XONDGL(( ELE - 1 ) * NLOC + ILOC )
       XC = XC + X( XNODI ) / REAL( NLOC )
       YC = YC + Y( XNODI ) / REAL( NLOC )
       ZC = ZC + Z( XNODI ) / REAL( NLOC )
    END DO

    SXC = 0.0
    SYC = 0.0
    SZC = 0.0
    DO SILOC = 1, SNLOC
       ILOC = SILOC2ILOC( SILOC )
       XNODI = XONDGL(( ELE - 1 ) * NLOC+ ILOC )
       SXC = SXC + X( XNODI ) / REAL( SNLOC )
       SYC = SYC + Y( XNODI ) / REAL( SNLOC )
       SZC = SZC + Z( XNODI ) / REAL( SNLOC )
    END DO
    NORMX = SXC - XC
    NORMY = SYC - YC
    NORMZ = SZC - ZC

    NORM = SQRT( NORMX**2 + NORMY**2 + NORMZ**2 )

    NORMX = NORMX / NORM
    NORMY = NORMY / NORM
    NORMZ = NORMZ / NORM

    RETURN

  END SUBROUTINE DGSIMPLNORM




  SUBROUTINE CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
       NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
       CV_SLOCLIST, X_NLOC, X_NDGLN)
    ! Calculate FACE_ELE - the list of elements surrounding an 
    ! element and referenced with a face -ve values correspond to surface elements. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: TOTELE, STOTEL, NFACE, NCOLELE, CV_NLOC, CV_SNLOC, CV_NONODS, &
         X_NLOC
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
    INTEGER, DIMENSION( NFACE, CV_SNLOC ), intent( in ) :: CV_SLOCLIST
    INTEGER, DIMENSION( NFACE, TOTELE ), intent( inout ) :: FACE_ELE
    INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
    INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
    ! Local variables
    LOGICAL, DIMENSION( : ), allocatable :: NOD_BELONG_ELE
    INTEGER, DIMENSION( : ), allocatable :: NOD_COUNT_SELE, FIN_ND_SELE, COL_ND_SELE, ELE_ROW
    LOGICAL, DIMENSION( : ), allocatable :: NOD_ON_BOUNDARY
    INTEGER :: ELE, ELE2, SELE, CV_ILOC, CV_INOD, COUNT, IFACE, CV_SILOC, JFACE, &
         CV_ILOC2, CV_INOD2, CV_SILOC2, IFACE2, SELE2 
    LOGICAL :: FOUND, SELE_FOUND, GOT_ALL

    ALLOCATE( NOD_BELONG_ELE( CV_NONODS ))
    ALLOCATE( NOD_COUNT_SELE( CV_NONODS ))
    ALLOCATE( FIN_ND_SELE( CV_NONODS + 1 ))
    ALLOCATE( NOD_ON_BOUNDARY( CV_NONODS )) 
    ALLOCATE( ELE_ROW( NFACE )) 

    FACE_ELE=0
    DO ELE=1,TOTELE
       IFACE=0
       DO COUNT=FINELE(ELE),FINELE(ELE+1)-1
          ELE2=COLELE(COUNT)
          IF(ELE2.NE.ELE) THEN
             IFACE=IFACE+1
             FACE_ELE(IFACE,ELE)=ELE2
          END IF
       END DO
    END DO


    NOD_COUNT_SELE=0
    NOD_ON_BOUNDARY=.FALSE.
    DO SELE=1,STOTEL
       DO CV_SILOC=1,CV_SNLOC
          CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC) 
          !   ewrite(3,*)'sele,CV_INOD:',sele,CV_INOD
          NOD_COUNT_SELE(CV_INOD)=NOD_COUNT_SELE(CV_INOD)+1
          NOD_ON_BOUNDARY(CV_INOD)=.TRUE.
       END DO
    END DO
    !    ewrite(3,*)'NOD_COUNT_SELE:',NOD_COUNT_SELE
    !    ewrite(3,*)'NOD_ON_BOUNDARY:',NOD_ON_BOUNDARY
    !    stop 831

    FIN_ND_SELE(1)=1
    DO CV_INOD=2,CV_NONODS+1
       FIN_ND_SELE(CV_INOD)=FIN_ND_SELE(CV_INOD-1)+NOD_COUNT_SELE(CV_INOD-1)
    END DO
    ALLOCATE(COL_ND_SELE(FIN_ND_SELE(CV_NONODS+1)-1))

    NOD_COUNT_SELE=0
    DO SELE=1,STOTEL
       DO CV_SILOC=1,CV_SNLOC
          CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC) 
          NOD_COUNT_SELE(CV_INOD)=NOD_COUNT_SELE(CV_INOD)+1
          COL_ND_SELE(FIN_ND_SELE(CV_INOD)-1+NOD_COUNT_SELE(CV_INOD))=SELE
       END DO
    END DO

    ! Put surface elements into FACE_ELE    
    NOD_BELONG_ELE=.FALSE.
    DO ELE=1,TOTELE
       DO CV_ILOC=1,CV_NLOC
          CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC) 
          NOD_BELONG_ELE(CV_INOD)=.TRUE.
       END DO
       DO IFACE=1,NFACE
          IF(FACE_ELE(IFACE,ELE)==0) THEN
             DO CV_ILOC=1,CV_NLOC
                CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC) 
                DO COUNT=FIN_ND_SELE(CV_INOD),FIN_ND_SELE(CV_INOD+1)-1
                   SELE=COL_ND_SELE(COUNT) 
                   FOUND=.FALSE.
                   DO JFACE=1,IFACE-1
                      IF(SELE == - FACE_ELE(JFACE,ELE)) FOUND=.TRUE.
                   END DO
                   IF(.NOT.FOUND) THEN ! SELE is a candidate.
                      SELE_FOUND=.TRUE.
                      DO CV_SILOC=1,CV_SNLOC
                         CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC) 
                         IF(.NOT.NOD_BELONG_ELE(CV_INOD)) SELE_FOUND=.FALSE. 
                      END DO
                      IF(SELE_FOUND) FACE_ELE(IFACE,ELE)= - SELE 
                   ENDIF
                END DO
             END DO
          ENDIF
       END DO
       DO CV_ILOC=1,CV_NLOC
          CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC) 
          NOD_BELONG_ELE(CV_INOD)=.FALSE.
       END DO
    END DO

    ! ***********************************************************    
    ! Now re-arrange ordering of elements in FACE_ELE so they are 
    ! consistent with the local nodes in CV_SLOCLIST
    DO ELE=1,TOTELE
       DO IFACE=1,NFACE
          ! Find the suitable element in row ELE of FACE_ELE for face IFACE. 
          DO IFACE2=1,NFACE
             ELE2=FACE_ELE(IFACE2,ELE)
             SELE2 = MAX( 0, - ELE2 )
             ELE2 = MAX( 0, + ELE2 )
             GOT_ALL=.TRUE.
             DO CV_SILOC=1,CV_SNLOC
                CV_ILOC=CV_SLOCLIST(IFACE,CV_SILOC)
                FOUND=.FALSE. 
                IF(SELE2 == 0) THEN ! is a volume element
                   CV_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC)
                   DO CV_ILOC2=1,CV_NLOC
                      CV_INOD2=X_NDGLN((ELE2-1)*X_NLOC+CV_ILOC2) 
                      IF(CV_INOD == CV_INOD2) FOUND=.TRUE.
                   END DO
                ELSE ! is a surface element
                   CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
                   DO CV_SILOC2=1,CV_SNLOC
                      !        ewrite(3,*)'(SELE2-1)*CV_SNLOC+CV_SILOC2,SELE2,CV_SNLOC,CV_SILOC2:', &
                      !                 (SELE2-1)*CV_SNLOC+CV_SILOC2,SELE2,CV_SNLOC,CV_SILOC2
                      CV_INOD2=CV_SNDGLN((SELE2-1)*CV_SNLOC+CV_SILOC2) 
                      IF(CV_INOD == CV_INOD2) FOUND=.TRUE.
                   END DO
                ENDIF
                IF(.NOT.FOUND) GOT_ALL=.FALSE.
             END DO
             IF(GOT_ALL) ELE_ROW(IFACE)=FACE_ELE(IFACE2,ELE)
          END DO
       END DO
       ! Re-order row...
       FACE_ELE(:,ELE)=ELE_ROW(:)
       !     ewrite(3,*)'FACE_ELE(:,ELE):',FACE_ELE(:,ELE)
    END DO
    !   stop 922

    !    do ele=1,totele
    !       ewrite(3,*)'ele',ele,' FACE_ELE(IFACE,ELE):',(FACE_ELE(IFACE,ELE),iface=1,nface) 
    !    end do
    !    stop 2982

    DEALLOCATE( NOD_BELONG_ELE )
    DEALLOCATE( NOD_COUNT_SELE )
    DEALLOCATE( FIN_ND_SELE )
    DEALLOCATE( NOD_ON_BOUNDARY ) 
    DEALLOCATE( COL_ND_SELE )
    DEALLOCATE( ELE_ROW ) 

    RETURN

  END SUBROUTINE CALC_FACE_ELE




  SUBROUTINE SURRO_CV_MINMAX( TMAX, TMIN, TOLDMAX, TOLDMIN, DENMAX, DENMIN, DENOLDMAX, DENOLDMIN, &
       T2MAX, T2MIN, T2OLDMAX, T2OLDMIN, &
       T, TOLD,  T2, T2OLD, DEN, DENOLD, IGOT_T2, NPHASE, CV_NONODS, NCOLACV, FINACV, COLACV, &
       STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC,  SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
       WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
       WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
       T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
       DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD )
    ! For each node, find the largest and smallest value of T and 
    ! DENSITY for both the current and previous timestep, out of 
    ! the node value and all its surrounding nodes including Dirichlet b.c's.
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE,CV_NONODS, NCOLACV,STOTEL,CV_SNLOC, &
         WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
         WIC_D_BC_DIRICHLET, IGOT_T2
    INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC, SUF_D_BC
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2 ), intent( in ) :: SUF_T2_BC
    INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_T_BC, WIC_D_BC
    INTEGER, DIMENSION( STOTEL * NPHASE * IGOT_T2), intent( in ) :: WIC_T2_BC
    INTEGER, DIMENSION( CV_NONODS * NPHASE + 1 ), intent( in ) :: FINACV
    INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: TMAX, TMIN, TOLDMAX, TOLDMIN,  &
         DENMAX, DENMIN, DENOLDMAX, DENOLDMIN
    REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: T2MAX, T2MIN, T2OLDMAX, T2OLDMIN
    REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: T,TOLD,DEN,DENOLD
    REAL, DIMENSION( CV_NONODS*NPHASE * IGOT_T2), intent( in ) :: T2,T2OLD
    REAL, DIMENSION( CV_NONODS ), intent( inout ) :: MASS_CV
    INTEGER, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
         TOLDMAX_NOD, DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD
    INTEGER, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, &
         T2OLDMAX_NOD
    ! Local variables
    INTEGER :: CV_NODI, CV_NODJ, CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, COUNT, CV_SILOC, SELE, &
         CV_INOD, SUF_CV_SI, SUF_CV_SI_IPHA, CV_INOD_IPHA

    Loop_CV_NODI: DO CV_NODI = 1, CV_NONODS

       DO IPHASE = 1, NPHASE
          CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
          TMAX( CV_NODI_IPHA ) = T( CV_NODI_IPHA )
          TMIN( CV_NODI_IPHA ) = T( CV_NODI_IPHA )
          TOLDMAX( CV_NODI_IPHA ) = TOLD( CV_NODI_IPHA )
          TOLDMIN( CV_NODI_IPHA ) = TOLD( CV_NODI_IPHA )
          IF(IGOT_T2==1) THEN
             T2MAX( CV_NODI_IPHA ) = T2( CV_NODI_IPHA )
             T2MIN( CV_NODI_IPHA ) = T2( CV_NODI_IPHA )
             T2OLDMAX( CV_NODI_IPHA ) = T2OLD( CV_NODI_IPHA )
             T2OLDMIN( CV_NODI_IPHA ) = T2OLD( CV_NODI_IPHA )
          ENDIF
          DENMAX( CV_NODI_IPHA ) = DEN( CV_NODI_IPHA )
          DENMIN( CV_NODI_IPHA ) = DEN( CV_NODI_IPHA )
          DENOLDMAX( CV_NODI_IPHA ) = DENOLD( CV_NODI_IPHA )
          DENOLDMIN( CV_NODI_IPHA ) = DENOLD( CV_NODI_IPHA )

          TMAX_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          TMIN_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          TOLDMAX_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          TOLDMIN_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          IF(IGOT_T2==1) THEN
             T2MAX_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
             T2MIN_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
             T2OLDMAX_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
             T2OLDMIN_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          ENDIF
          DENMAX_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          DENMIN_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          DENOLDMAX_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
          DENOLDMIN_NOD( CV_NODI_IPHA ) = CV_NODI_IPHA 
       END DO

       DO COUNT = FINACV( CV_NODI ), FINACV( CV_NODI + 1 ) - 1, 1
          CV_NODJ = COLACV( COUNT )

          IF( CV_NODJ <= CV_NONODS ) THEN

             DO IPHASE = 1,NPHASE
                CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS
                IF( T( CV_NODJ_IPHA ) > TMAX( CV_NODI_IPHA ) ) THEN
                   TMAX( CV_NODI_IPHA ) = T( CV_NODJ_IPHA )
                   TMAX_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF
                IF( T( CV_NODJ_IPHA ) < TMIN( CV_NODI_IPHA ) ) THEN
                   TMIN( CV_NODI_IPHA ) = T( CV_NODJ_IPHA )
                   TMIN_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF

                IF( TOLD( CV_NODJ_IPHA ) > TOLDMAX( CV_NODI_IPHA ) ) THEN
                   TOLDMAX( CV_NODI_IPHA ) = TOLD( CV_NODJ_IPHA )
                   TOLDMAX_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF
                IF( TOLD( CV_NODJ_IPHA ) < TOLDMIN( CV_NODI_IPHA ) ) THEN
                   TOLDMIN( CV_NODI_IPHA ) = TOLD( CV_NODJ_IPHA )
                   TOLDMIN_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF
                ! T2: 
                IF(IGOT_T2==1) THEN
                   IF( T2( CV_NODJ_IPHA ) > T2MAX( CV_NODI_IPHA ) ) THEN
                      T2MAX( CV_NODI_IPHA ) = T2( CV_NODJ_IPHA )
                      T2MAX_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                   ENDIF
                   IF( T2( CV_NODJ_IPHA ) < T2MIN( CV_NODI_IPHA ) ) THEN
                      T2MIN( CV_NODI_IPHA ) = T2( CV_NODJ_IPHA )
                      T2MIN_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                   ENDIF

                   IF( T2OLD( CV_NODJ_IPHA ) > T2OLDMAX( CV_NODI_IPHA ) ) THEN
                      T2OLDMAX( CV_NODI_IPHA ) = T2OLD( CV_NODJ_IPHA )
                      T2OLDMAX_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                   ENDIF
                   IF( T2OLD( CV_NODJ_IPHA ) < T2OLDMIN( CV_NODI_IPHA ) ) THEN
                      T2OLDMIN( CV_NODI_IPHA ) = T2OLD( CV_NODJ_IPHA )
                      T2OLDMIN_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                   ENDIF
                ENDIF
                ! DEN:  
                IF( DEN( CV_NODJ_IPHA ) > DENMAX( CV_NODI_IPHA ) ) THEN
                   DENMAX( CV_NODI_IPHA ) = DEN( CV_NODJ_IPHA )
                   DENMAX_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF
                IF( DEN( CV_NODJ_IPHA ) < DENMIN( CV_NODI_IPHA ) ) THEN
                   DENMIN( CV_NODI_IPHA ) = DEN( CV_NODJ_IPHA )
                   DENMIN_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF

                IF( DENOLD( CV_NODJ_IPHA ) > DENOLDMAX( CV_NODI_IPHA ) ) THEN
                   DENOLDMAX( CV_NODI_IPHA ) = DENOLD( CV_NODJ_IPHA )
                   DENOLDMAX_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF
                IF( DENOLD( CV_NODJ_IPHA ) < DENOLDMIN( CV_NODI_IPHA ) ) THEN
                   DENOLDMIN( CV_NODI_IPHA ) = DENOLD( CV_NODJ_IPHA )
                   DENOLDMIN_NOD( CV_NODI_IPHA ) = CV_NODJ_IPHA 
                ENDIF

             END DO

          ENDIF

       END DO

    END DO Loop_CV_NODI

    ! Take into account the Dirichlet b.c's when working out max and min values.
    Loop_SELE: DO SELE= 1, STOTEL

       Loop_CV_SILOC: DO CV_SILOC = 1, CV_SNLOC

          CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC)
          SUF_CV_SI=(SELE-1)*CV_SNLOC+CV_SILOC

          DO IPHASE=1,NPHASE
             SUF_CV_SI_IPHA = SUF_CV_SI + STOTEL * CV_SNLOC * ( IPHASE - 1 )
             CV_INOD_IPHA=CV_INOD + CV_NONODS*(IPHASE-1)
             IF( (WIC_T_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRICHLET) &
                  .OR.(WIC_T_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRI_ADV_AND_ROBIN)) THEN
                IF(SUF_T_BC( SUF_CV_SI_IPHA ) > TMAX( CV_INOD_IPHA ) ) THEN
                   TMAX( CV_INOD_IPHA ) = SUF_T_BC( SUF_CV_SI_IPHA )
                   TMAX_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF
                IF(SUF_T_BC( SUF_CV_SI_IPHA ) < TMIN( CV_INOD_IPHA ) ) THEN
                   TMIN( CV_INOD_IPHA ) = SUF_T_BC( SUF_CV_SI_IPHA )
                   TMIN_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF

                IF(SUF_T_BC( SUF_CV_SI_IPHA ) > TOLDMAX( CV_INOD_IPHA ) ) THEN
                   TOLDMAX( CV_INOD_IPHA ) = SUF_T_BC( SUF_CV_SI_IPHA )
                   TOLDMAX_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF
                IF(SUF_T_BC( SUF_CV_SI_IPHA ) < TOLDMIN( CV_INOD_IPHA ) ) THEN
                   TOLDMIN( CV_INOD_IPHA ) = SUF_T_BC( SUF_CV_SI_IPHA )
                   TOLDMIN_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF
             ENDIF

             ! T2:
             IF(IGOT_T2==1) THEN
                IF( (WIC_T2_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRICHLET) &
                     .OR.(WIC_T2_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRI_ADV_AND_ROBIN)) THEN
                   IF(SUF_T2_BC( SUF_CV_SI_IPHA ) > T2MAX( CV_INOD_IPHA ) ) THEN
                      T2MAX( CV_INOD_IPHA ) = SUF_T2_BC( SUF_CV_SI_IPHA )
                      T2MAX_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                   ENDIF
                   IF(SUF_T2_BC( SUF_CV_SI_IPHA ) < T2MIN( CV_INOD_IPHA ) ) THEN
                      T2MIN( CV_INOD_IPHA ) = SUF_T2_BC( SUF_CV_SI_IPHA )
                      T2MIN_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                   ENDIF

                   IF(SUF_T2_BC( SUF_CV_SI_IPHA ) > T2OLDMAX( CV_INOD_IPHA ) ) THEN
                      T2OLDMAX( CV_INOD_IPHA ) = SUF_T2_BC( SUF_CV_SI_IPHA )
                      T2OLDMAX_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                   ENDIF
                   IF(SUF_T2_BC( SUF_CV_SI_IPHA ) < T2OLDMIN( CV_INOD_IPHA ) ) THEN
                      T2OLDMIN( CV_INOD_IPHA ) = SUF_T2_BC( SUF_CV_SI_IPHA )
                      T2OLDMIN_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                   ENDIF
                   !                print*, 'T2:',TMAX( CV_INOD_IPHA ), T2MIN( CV_INOD_IPHA ), T2OLDMAX( CV_INOD_IPHA ) , T2OLDMIN( CV_INOD_IPHA )
                ENDIF
             ENDIF
             ! DEN: 
             !          print *,'stotel,nphase=', stotel,nphase
             !          print *,'SELE,SELE+(IPHASE-1)*STOTEL:', SELE,SELE+(IPHASE-1)*STOTEL
             IF( WIC_D_BC(SELE+(IPHASE-1)*STOTEL) == WIC_D_BC_DIRICHLET ) THEN
                IF(SUF_D_BC( SUF_CV_SI_IPHA ) > DENMAX( CV_INOD_IPHA ) ) THEN
                   DENMAX( CV_INOD_IPHA ) = SUF_D_BC( SUF_CV_SI_IPHA )
                   DENMAX_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF
                IF(SUF_D_BC( SUF_CV_SI_IPHA ) < DENMIN( CV_INOD_IPHA ) ) THEN
                   DENMIN( CV_INOD_IPHA ) = SUF_D_BC( SUF_CV_SI_IPHA )
                   DENMIN_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF

                IF(SUF_D_BC( SUF_CV_SI_IPHA ) > DENOLDMAX( CV_INOD_IPHA ) ) THEN
                   DENOLDMAX( CV_INOD_IPHA ) = SUF_D_BC( SUF_CV_SI_IPHA )
                   DENOLDMAX_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF
                IF(SUF_D_BC( SUF_CV_SI_IPHA ) < DENOLDMIN( CV_INOD_IPHA ) ) THEN
                   DENOLDMIN( CV_INOD_IPHA ) = SUF_D_BC( SUF_CV_SI_IPHA )
                   DENOLDMIN_NOD( CV_INOD_IPHA ) =  CV_INOD_IPHA
                ENDIF
                !                print*, 'Dens:',DENMAX( CV_INOD_IPHA ), DENMIN( CV_INOD_IPHA ), DENOLDMAX( CV_INOD_IPHA ),DENOLDMIN( CV_INOD_IPHA )
             ENDIF
          END DO

       END DO Loop_CV_SILOC

    END DO Loop_SELE

    RETURN
  END SUBROUTINE SURRO_CV_MINMAX




  SUBROUTINE CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, CV_NGI, U_SLOC2LOC, CV_SLOC2LOC, &
       FACE_ELE, TOTELE, NFACE, CV_ON_FACE, GI, &
       CV_NONODS, LOG_ON_BOUND, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, STOTEL, &
       CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN ) 
    ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC for a face on the
    ! boundary of the domain
    IMPLICIT NONE
    INTEGER, intent( in ) :: ELE, CV_NGI, CV_NONODS, TOTELE, STOTEL, CV_NLOC, U_NLOC, &
         CV_SNLOC, U_SNLOC, NFACE, GI, CV_ILOC
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
    INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN
    INTEGER, intent( inout ) :: SELE, CV_SILOC
    INTEGER, DIMENSION( NFACE, TOTELE ), intent( in ) :: FACE_ELE
    LOGICAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in )  :: CV_ON_FACE
    INTEGER, DIMENSION( U_SNLOC ), intent( inout ) :: U_SLOC2LOC
    INTEGER, DIMENSION( CV_SNLOC ), intent( inout ) :: CV_SLOC2LOC
    LOGICAL, DIMENSION( CV_NONODS ), intent( inout ) :: LOG_ON_BOUND
    ! local variables
    INTEGER :: IFACE, ELE2, SELE2, CV_JLOC, CV_JNOD, &
         U_JLOC, U_JNOD, CV_KLOC, CV_SKNOD, &
         U_KLOC, U_SKLOC, U_SKNOD, CV_SKLOC
    LOGICAL :: FOUND
    !    ewrite(3,*)'CV_ON_FACE:',CV_ON_FACE
    DO CV_JLOC = 1, CV_NLOC  
       CV_JNOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )
       !       ewrite(3,*)'CV_NONODS,CV_JNOD,CV_NLOC,CV_NGI,CV_JLOC, GI:',  &
       !                CV_NONODS,CV_JNOD,CV_NLOC,CV_NGI,CV_JLOC, GI
       LOG_ON_BOUND( CV_JNOD ) = CV_ON_FACE( CV_JLOC, GI )
    END DO
    !    stop 2821

    ! What face are we on        
    DO IFACE = 1, NFACE
       ELE2=FACE_ELE(IFACE,ELE)
       SELE2 = MAX( 0, - ELE2 )
       ELE2 = MAX( 0, + ELE2 )
       IF(SELE2 /= 0) THEN
          FOUND = .TRUE.
          DO CV_SKLOC = 1, CV_SNLOC  
             CV_SKNOD = CV_SNDGLN(( SELE2 - 1 ) * CV_SNLOC + CV_SKLOC )
             IF( .NOT. LOG_ON_BOUND( CV_SKNOD )) FOUND=.FALSE.
          END DO
          IF( FOUND ) SELE = SELE2
       END IF
    END DO
    !    stop 3983

    ! Calculate CV_SLOC2LOC    
    !    ewrite(3,*)'sele2=',sele2
    DO CV_SKLOC = 1, CV_SNLOC  
       CV_SKNOD = CV_SNDGLN(( SELE - 1 ) * CV_SNLOC + CV_SKLOC )
       DO CV_JLOC = 1, CV_NLOC  
          CV_JNOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )
          IF( CV_SKNOD == CV_JNOD ) CV_KLOC = CV_JLOC
       END DO
       CV_SLOC2LOC( CV_SKLOC ) = CV_KLOC
       IF(CV_KLOC.EQ.CV_ILOC) CV_SILOC=CV_SKLOC
    END DO

    ! Calculate U_SLOC2LOC  
    DO U_SKLOC = 1, U_SNLOC  
       U_SKNOD = U_SNDGLN(( SELE - 1 ) * U_SNLOC + U_SKLOC )
       DO U_JLOC = 1, U_NLOC  
          U_JNOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_JLOC )
          IF( U_SKNOD == U_JNOD ) U_KLOC = U_JLOC
       END DO
       U_SLOC2LOC( U_SKLOC ) = U_KLOC
    END DO

    RETURN   
  END SUBROUTINE CALC_SELE




  !
  SUBROUTINE RE1DN3(NGI,NLOC,WEIGHT,N,NLX )
    IMPLICIT NONE
    ! QUADRATIC VARIATION FOR VELOCITY-2D 9 NODE BRICK ELEMENT.
    ! LINEAR VARIATION FOR PRESSURE-2D 4 NODE BRICK ELEMENT.
    ! NB might have to define surface elements for p and (u,v,w) 
    ! in here as well. 
    ! This is for the 2-D 27node element, which is number as follows
    !   1   2   3
    !      This subroutine defines the shape functions M and N and their
    !      derivatives at the Gauss points
    INTEGER NGI,NLOC
    REAL WEIGHT(NGI)
    REAL N(NLOC,NGI),NLX(NLOC,NGI)
    REAL POSI
    REAL LX(3),LY(3)
    REAL WEI(3)
    REAL XN(3)
    REAL DXN(3)
    INTEGER P,GPOI

    INTEGER NQUAD,ILX,NL
    ! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I

    POSI=0.774596669241483
    LX(1)=-POSI
    LY(1)=-POSI
    LX(2)= 0.
    LY(2)= 0.
    LX(3)= POSI
    LY(3)= POSI
    WEI(1)=0.555555555555556
    WEI(2)=0.888888888888889
    WEI(3)=0.555555555555556
    NQUAD=3
    !
    !
    !  FIND N ETC-----
    ! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I
    !
    do  P=1,NQUAD! Was loop 23
       GPOI=P
       !
       WEIGHT(GPOI)=WEI(P)
       !
       XN(1)=0.5* LX(P) * (LX(P)-1.)
       XN(2)=1. - LX(P)*LX(P)
       XN(3)=0.5* LX(P) * (LX(P)+1.)
       !
       DXN(1)=0.5* (2.*LX(P) -1.)
       DXN(2)= - 2. * LX(P)
       DXN(3)=0.5* (2.*LX(P) + 1.)
       !
       do  ILX=1,3! Was loop 10
          NL=ILX
          N(NL,GPOI)  = XN(ILX) 
          !
          NLX(NL,GPOI)= DXN(ILX) 
       end do ! Was loop 10
       !
    end do ! Was loop 23
  END SUBROUTINE RE1DN3




  SUBROUTINE PUT_IN_CT_RHS(CT, CT_RHS, U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
       CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE, ELE, ELE2, SELE, &
       JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC, U_NDGLN,  NU, NV, NW,  &
       SUFEN, SCVDETWEI, CVNORMX, CVNORMY, CVNORMZ, DENOLD, CV_NODI, CV_NODI_IPHA, &
       UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
       UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
       NDOTQ, NDOTQOLD, LIMDT, LIMDTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD)
    ! This subroutine caculates the discretised cty eqn acting on the velocities i.e. CT, CT_RHS
    IMPLICIT NONE
    INTEGER, intent( in ) :: U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
         CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE,  ELE, ELE2, SELE, &
         CV_NODI, CV_NODI_IPHA
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( U_NLOC ), intent( in ) :: JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC
    REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( inout ) :: CT
    REAL, DIMENSION( CV_NONODS ), intent( inout ) :: CT_RHS
    REAL, DIMENSION( U_NLOC ), intent( in ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( U_NLOC, SCVNGI ), intent( in ) :: SUFEN
    REAL, DIMENSION( SCVNGI ), intent( in ) :: SCVDETWEI, CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DENOLD
    REAL, intent( in ) :: NDOTQ, NDOTQOLD, LIMDT, LIMDTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD

    ! Local variables...
    INTEGER :: U_KLOC, U_KLOC2, JCOUNT_IPHA, IDIM, U_NODK, U_NODK_IPHA, JCOUNT2_IPHA, &
         U_KLOC_LEV, U_NLOC_LEV
    REAL :: RCON,UDGI_IMP,VDGI_IMP,WDGI_IMP,NDOTQ_IMP


    !      U_NLOC_LEV =U_NLOC /CV_NLOC

    !      DO U_KLOC_LEV = 1, U_NLOC_LEV
    !  U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
    DO U_KLOC = 1, U_NLOC

       JCOUNT_IPHA = JCOUNT_KLOC( U_KLOC )  +  (IPHASE - 1 ) * NCOLCT * NDIM
       RCON    = SCVDETWEI( GI ) * FTHETA_T2 * LIMDT  &
            * SUFEN( U_KLOC, GI ) / DENOLD( CV_NODI_IPHA )

       IDIM = 1
       CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) =  CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT) &
            +  RCON *UGI_COEF_ELE(U_KLOC)* CVNORMX( GI )

       IDIM = 2
       IF( NDIM >= 2 ) &
            CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) =  CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) &
            +  RCON * VGI_COEF_ELE(U_KLOC)* CVNORMY( GI )

       IDIM = 3
       IF( NDIM >= 3 ) &
            CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) =  CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) &
            +  RCON * WGI_COEF_ELE(U_KLOC)* CVNORMZ( GI )
    END DO

    IF(SELE /= 0) THEN
       UDGI_IMP=0.0
       VDGI_IMP=0.0
       WDGI_IMP=0.0
       !         DO U_KLOC_LEV = 1, U_NLOC_LEV
       !     U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
       DO U_KLOC = 1, U_NLOC
          U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
          U_NODK_IPHA = U_NODK + (IPHASE-1)*U_NONODS
          UDGI_IMP=UDGI_IMP + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE(U_KLOC) * NU(U_NODK_IPHA) 
          VDGI_IMP=VDGI_IMP + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE(U_KLOC) * NV(U_NODK_IPHA) 
          WDGI_IMP=WDGI_IMP + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE(U_KLOC) * NW(U_NODK_IPHA) 
       END DO
       NDOTQ_IMP=CVNORMX( GI ) * UDGI_IMP + CVNORMY( GI ) * VDGI_IMP + CVNORMZ( GI ) * WDGI_IMP

       CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - SCVDETWEI( GI ) * ( &
            ONE_M_FTHETA_T2OLD * LIMDTOLD * NDOTQOLD &
            + FTHETA_T2  * LIMDT * (NDOTQ-NDOTQ_IMP)        )/ DENOLD( CV_NODI_IPHA ) 
    ELSE
       CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - SCVDETWEI( GI ) * ( &
            ONE_M_FTHETA_T2OLD * LIMDTOLD * NDOTQOLD &
            )/ DENOLD( CV_NODI_IPHA ) 
    ENDIF

    IF((ELE2 /= 0).AND.(ELE2 /= ELE)) THEN
       ! We have a discontinuity between elements so integrate along the face...
       !        DO U_KLOC_LEV = 1, U_NLOC_LEV
       !    U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
       DO U_KLOC = 1, U_NLOC
          U_KLOC2=U_OTHER_LOC(U_KLOC)
          IF(U_KLOC2 /= 0) THEN
             JCOUNT2_IPHA = JCOUNT_KLOC2( U_KLOC2 ) + (IPHASE - 1) * NCOLCT * NDIM

             RCON    = SCVDETWEI( GI ) * FTHETA_T2 * LIMDT  &
                  * SUFEN( U_KLOC, GI ) / DENOLD( CV_NODI_IPHA )

             IDIM = 1
             CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                  =  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT) &
                  +  RCON *UGI_COEF_ELE2(U_KLOC2)* CVNORMX( GI )

             IDIM = 2
             IF( NDIM >= 2 ) &
                  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                  =  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                  +  RCON * VGI_COEF_ELE2(U_KLOC2)* CVNORMY( GI )

             IDIM = 3
             IF( NDIM >= 3 ) &
                  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                  =  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                  +  RCON * WGI_COEF_ELE2(U_KLOC2)* CVNORMZ( GI )
          ENDIF
       END DO
    ENDIF

    RETURN
  END SUBROUTINE PUT_IN_CT_RHS





  SUBROUTINE PRINT_CV_DIST(CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE, &
       SATURA, X_NDGLN, CV_NDGLN, X) 
    IMPLICIT NONE

    INTEGER, intent( in ) :: CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SATURA
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    ! Local Variables
    INTEGER :: ELE, CV_ILOC, IPHASE
    REAL :: x_coord

    ewrite(3,*)'satura :',satura
    do iphase=1,nphase
       ewrite(3,*)'cv REPRESENTATION OF iphase:',iphase
       do ele=1,totele

          IF(CV_NLOC==2) THEN
             cv_iloc=1
             ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
             x_coord=0.5*( x(x_ndgln((ele-1)*x_nloc+cv_iloc)) + x(x_ndgln((ele-1)*x_nloc+cv_iloc+1)) )
             ewrite(3,*)x_coord,SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)

             cv_iloc=cv_nloc
             x_coord=0.5*( x(x_ndgln((ele-1)*x_nloc+cv_iloc)) + x(x_ndgln((ele-1)*x_nloc+cv_iloc-1)) )
             ewrite(3,*)x_coord,SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
             ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
          ENDIF
          IF(CV_NLOC==3) THEN
             cv_iloc=1
             ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
             x_coord=0.5*( x(x_ndgln((ele-1)*x_nloc+cv_iloc)) + x(x_ndgln((ele-1)*x_nloc+cv_iloc+1)) )
             ewrite(3,*)x_coord,SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)

             do cv_iloc=2,cv_nloc-1
                x_coord=0.5*( x(x_ndgln((ele-1)*x_nloc+cv_iloc)) + x(x_ndgln((ele-1)*x_nloc+cv_iloc-1)) )
                ewrite(3,*)x_coord,SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
                x_coord=0.5*( x(x_ndgln((ele-1)*x_nloc+cv_iloc)) + x(x_ndgln((ele-1)*x_nloc+cv_iloc+1)) )
                ewrite(3,*)x_coord,SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
                !        ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),DTX_ELE(cv_iloc,1,ele)
             end do

             cv_iloc=cv_nloc
             x_coord=0.5*( x(x_ndgln((ele-1)*x_nloc+cv_iloc)) + x(x_ndgln((ele-1)*x_nloc+cv_iloc-1)) )
             ewrite(3,*)x_coord,SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
             ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),SATURA(cv_ndgln((ele-1)*cv_nloc+cv_iloc)+(iphase-1)*cv_nonods)
          ENDIF

          ewrite(3,*)'iphase, nphase, ele, totele: ',iphase, nphase, ele, totele
       end do
    end do
    RETURN
  END SUBROUTINE PRINT_CV_DIST





  SUBROUTINE FIND_OPT_INCOME_INTERP( INCOME, W, NDOTQ, NDOTQ2, IPHASE, &
       SAT_CV_NODI_IPHA, SAT_CV_NODJ_IPHA, SAT_FEM_IPHA2, SAT_GEOM_IPHA2, OVER_RELAX, &
       ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
       GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA, CONSERV, MAX_OPER, SAT_BASED )
    ! calculate INCOME & UUNDOTQ for optimal upwinding.
    IMPLICIT NONE
    REAL, intent( inout ) :: INCOME, W
    REAL, intent( in ) :: NDOTQ, NDOTQ2
    INTEGER, intent( in ) :: IPHASE
    REAL, intent( in ) :: SAT_CV_NODI_IPHA, SAT_CV_NODJ_IPHA, SAT_FEM_IPHA2, SAT_GEOM_IPHA2, OVER_RELAX
    REAL, intent( inout ) :: ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA
    REAL, intent( in ) :: GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA
    ! If OVER_RELAX >1 then over relax the upwinding e.g.1 or 2 might be used.
    ! OPT_VEL_UPWIND_COEFS contains the coefficients
    LOGICAL, intent( in ) ::  CONSERV,MAX_OPER,SAT_BASED
    ! local variables...    
    LOGICAL, PARAMETER :: VEL_BASE_INT = .true.
    ! CONSERV then use a highly conservative estimat of upwind parameter. 
    REAL :: ABSC, A, B, NDOTQ_MEAN, SATC, UUNDOTQ, SAT_FEM_IPHA, MAX_SAT, MIN_SAT, &
         SAT_GEOM_IPHA
    REAL :: SAT_CV_NODI, SAT_CV_NODJ, ABS_CV_NODI, ABS_CV_NODJ, ABSC1, ABSC2, WA
    !!  real tolfun

    MIN_SAT=MIN(SAT_CV_NODI_IPHA, SAT_CV_NODJ_IPHA)
    MAX_SAT=MAX(SAT_CV_NODI_IPHA, SAT_CV_NODJ_IPHA)
    SAT_FEM_IPHA=MIN(MAX_SAT,MAX(MIN_SAT,SAT_FEM_IPHA2))
    SAT_GEOM_IPHA=MIN(MAX_SAT,MAX(MIN_SAT,SAT_GEOM_IPHA2))

    IF(VEL_BASE_INT) THEN !****************

       IF((ABS(SAT_CV_NODI_IPHA - SAT_CV_NODJ_IPHA) < 1.E-4) &
            .OR.(ABS(NDOTQ - NDOTQ2) < 1.E-6)) THEN
          INCOME=0.5
          W=0.5
          RETURN
       ENDIF

       ! Find end pt absorptions...
       NDOTQ_MEAN=0.5*(NDOTQ + NDOTQ2)

       Conditional_SATBASED: IF(.NOT.SAT_BASED) THEN
          ! assume a linear variation from downwide node... DA/DS = GRAD_ABS
          Conditional_Conserv: if(.NOT.CONSERV) then
             IF(NDOTQ_MEAN < 0) THEN 
                ABSC=ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA) 
             ELSE
                ABSC=ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)
             ENDIF
             !!        ABSC=0.5*ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA) &
             !!            +0.5*ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA) 
          else
             ABSC1=ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)
             ABSC2=ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA) 
             !        ABSC1=ABS_CV_NODI_IPHA + (SAT_FEM_IPHA - SAT_CV_NODI_IPHA) &
             !     *(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
             !        ABSC2=ABS_CV_NODJ_IPHA + (SAT_FEM_IPHA - SAT_CV_NODJ_IPHA) &
             !     *(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
             IF(NDOTQ_MEAN < 0) THEN 
                if(abs(ABSC1-ABS_CV_NODJ_IPHA) < abs(ABSC2-ABS_CV_NODJ_IPHA)) then
                   ABSC=ABSC1
                else
                   ABSC=ABSC2
                endif
             ELSE
                if(abs(ABSC1-ABS_CV_NODI_IPHA) < abs(ABSC2-ABS_CV_NODI_IPHA)) then
                   ABSC=ABSC1
                else
                   ABSC=ABSC2
                endif
             ENDIF
          endif Conditional_Conserv
          !        ABSC=0.5*(ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)) &
          !            +0.5*(ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA)) 

          if(.false.) then
             if(iphase==2) then
                CALL ABS3P(ABS_CV_NODI_IPHA, 1.0, 1.-SAT_CV_NODI_IPHA, IPHASE)
                CALL ABS3P(ABS_CV_NODJ_IPHA, 1.0, 1.-SAT_CV_NODJ_IPHA, IPHASE)
                CALL ABS3P(ABSC, 1.0, 1.-SAT_FEM_IPHA, IPHASE)
             else
                CALL ABS3P(ABS_CV_NODI_IPHA, 1.0, SAT_CV_NODI_IPHA, IPHASE)
                CALL ABS3P(ABS_CV_NODJ_IPHA, 1.0, SAT_CV_NODJ_IPHA, IPHASE)
                CALL ABS3P(ABSC, 1.0, SAT_FEM_IPHA, IPHASE)
             endif
          endif

          UUNDOTQ = 0.5*(1./tolfun(ABSC))*(NDOTQ*ABS_CV_NODI_IPHA + NDOTQ2*ABS_CV_NODJ_IPHA)
       ENDIF Conditional_SATBASED

       !!  IF(UUNDOTQ < 0.0) THEN
       Conditional_NDOTQMean: IF(NDOTQ_MEAN < 0.0) THEN
          !     W=(UUNDOTQ - NDOTQ)/(NDOTQ2 - NDOTQ)
          !       w=(ABSC-ABS_CV_NODi_IPHA)/(ABS_CV_NODj_IPHA-ABS_CV_NODi_IPHA)
          !       w=(ABS_CV_NODj_IPHA/tolfun(ABSC))*(ABSC-ABS_CV_NODi_IPHA)/(ABS_CV_NODj_IPHA-ABS_CV_NODi_IPHA)
          IF(SAT_BASED) THEN
             W=(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)/(SAT_CV_NODJ_IPHA-SAT_CV_NODI_IPHA)
          ELSE
             IF(MAX_OPER) THEN
                w=max(1.,ABS_CV_NODj_IPHA/tolfun(ABSC))*(ABSC-ABS_CV_NODi_IPHA)/(ABS_CV_NODj_IPHA-ABS_CV_NODi_IPHA)
             ELSE
                w=max(0.,ABS_CV_NODj_IPHA/tolfun(ABSC))*(ABSC-ABS_CV_NODi_IPHA)/(ABS_CV_NODj_IPHA-ABS_CV_NODi_IPHA)
             ENDIF
          ENDIF
          !!     w=0.7
          !!     WA=(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)/(SAT_CV_NODJ_IPHA-SAT_CV_NODI_IPHA)
          !!     W=(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)/(SAT_CV_NODJ_IPHA-SAT_CV_NODI_IPHA)
          !     if(SAT_GEOM_IPHA2<-100.) then
          !       wa=1.0
          !       w=1.0
          !     else
          WA=(SAT_GEOM_IPHA - SAT_CV_NODI_IPHA)/(SAT_CV_NODJ_IPHA-SAT_CV_NODI_IPHA)
          !     endif
          !     if(iphase==2) w=wa
          !     wa=0.5
          !     wa=0.0
          !     if(iphase==2) w=wa
          !     w=wa
          !     wa=min(wa,0.5)
          !!     w=wa
          !!     w=0.8
          W=OVER_RELAX*(w-WA)+WA
          !     W=OVER_RELAX*(w-0.5)+0.5
          !     W=max(w,0.5)
          !     W=max(w,WA,0.5)
          W=max(w,WA,0.)
          !     W=max(w,0.0)
          W=min(w,1.0)
          !     W=min(w,0.9)
          INCOME= W
       ELSE
          !     W=(UUNDOTQ - NDOTQ2)/(NDOTQ - NDOTQ2)
          !       w=(ABSC-ABS_CV_NODJ_IPHA)/(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)
          !       w=(ABS_CV_NODi_IPHA/tolfun(ABSC))*(ABSC-ABS_CV_NODJ_IPHA)/(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)
          IF(SAT_BASED) THEN
             W=(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
          ELSE
             IF(MAX_OPER) THEN
                w=max(1.,ABS_CV_NODi_IPHA/tolfun(ABSC))*(ABSC-ABS_CV_NODJ_IPHA)/(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)
             ELSE
                w=max(0.,ABS_CV_NODi_IPHA/tolfun(ABSC))*(ABSC-ABS_CV_NODJ_IPHA)/(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)
             ENDIF
          ENDIF
          !!     w=0.7
          !!     WA=(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
          !!     W=(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
          !     if(SAT_GEOM_IPHA2<-100.) then
          !       wa=1.0
          !       w=1.0
          !     else
          WA=(SAT_GEOM_IPHA - SAT_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
          !     endif
          !     if(iphase==2) w=wa
          if(.false.) then
             if((iphase==2).and.(abs(SAT_CV_NODI_IPHA-0.6)<0.1)) then
                if(abs(SAT_CV_NODJ_IPHA-0.6)<0.1) then
                   if(abs(wa-0.5)<0.01) then
                      ewrite(3,*)'iphase=',iphase
                      ewrite(3,*)'SAT_GEOM_IPHA,SAT_FEM_IPHA:',SAT_GEOM_IPHA,SAT_FEM_IPHA
                      ewrite(3,*)'SAT_CV_NODI_IPHA,SAT_CV_NODJ_IPHA:',SAT_CV_NODI_IPHA,SAT_CV_NODJ_IPHA
                      ewrite(3,*)'ABSC:',ABSC
                      ewrite(3,*)'GRAD_ABS_CV_NODI_IPHA,GRAD_ABS_CV_NODJ_IPHA:', &
                           GRAD_ABS_CV_NODI_IPHA,GRAD_ABS_CV_NODJ_IPHA
                      ewrite(3,*)'ABS_CV_NODI_IPHA,ABS_CV_NODJ_IPHA:',ABS_CV_NODI_IPHA,ABS_CV_NODJ_IPHA
                      ewrite(3,*)'w,wa:',w,wa
                      ewrite(3,*)'ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA):', &
                           ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)
                      ewrite(3,*)'ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA):', &
                           ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA)
                      ewrite(3,*)'(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA - SAT_CV_NODJ_IPHA):', &
                           (ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA - SAT_CV_NODJ_IPHA)
                      ewrite(3,*)'(ABSC-ABS_CV_NODJ_IPHA)/(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA):', &
                           (ABSC-ABS_CV_NODJ_IPHA)/(ABS_CV_NODI_IPHA-ABS_CV_NODJ_IPHA)
                      stop 2929
                   endif
                endif
             endif
          endif
          !     wa=0.5
          !     wa=0.0
          !     if(iphase==2) w=wa
          !     w=wa
          !     wa=min(wa,0.5)
          !!     w=wa
          !!     w=0.8
          W=OVER_RELAX*(w-WA)+WA
          !     W=OVER_RELAX*(w-0.5)+0.5
          !     W=max(w,0.5)
          !     W=max(w,WA,0.5)
          W=max(w,WA,0.)
          !     W=max(w,0.0)
          W=min(w,1.0)
          !     W=min(w,0.9)
          INCOME= 1.0 - W
       ENDIF Conditional_NDOTQMean
       !  ewrite(3,*)'w,income:',w,income



    ELSE  !****************

       IF(ABS(SAT_CV_NODI_IPHA - SAT_CV_NODJ_IPHA) < 1.E-4) THEN
          INCOME=0.5
          RETURN
       ENDIF

       NDOTQ_MEAN=0.5*(NDOTQ + NDOTQ2)

       ! assume a linear variation from downwide node... DA/DS = GRAD_ABS
       IF(NDOTQ_MEAN < 0) THEN 

          if(.false.) then
             SAT_CV_NODI=SAT_CV_NODI_ipha
             SAT_CV_NODj=SAT_CV_NODj_ipha
             if(iphase==2) then
                SAT_CV_NODI=1.-SAT_CV_NODI_ipha
                SAT_CV_NODj=1.-SAT_CV_NODj_ipha
             endif
             SAT_CV_NODI=min(1.0,max(0.0,SAT_CV_NODI))
             SAT_CV_NODj=min(1.0,max(0.0,SAT_CV_NODj))

             SATC=SAT_CV_NODJ + 0.01*(SAT_CV_NODI-SAT_CV_NODJ) 
             SATC=min(1.0,max(0.0,SATC))
             CALL ABS3P(ABSC, 1.0, SATC, IPHASE)
             CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
             CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
             ewrite(3,*)'ABS_CV_NODJ, SAT_CV_NODJ, IPHASE:',ABS_CV_NODJ, SAT_CV_NODJ, IPHASE
             B=(ABS_CV_NODJ-ABSC)/(SAT_CV_NODJ-SATC)
             A=ABS_CV_NODJ - B*SAT_CV_NODJ

             ABSC=a+b*0.5*(SAT_CV_NODJ+SAT_CV_NODI)
             ewrite(3,*)'guessed valueS ABSC,ABS_CV_NODJ=',ABSC,ABS_CV_NODJ
             ewrite(3,*)'SAT_CV_NODI,SAT_CV_NODj:',SAT_CV_NODI,SAT_CV_NODj
             ewrite(3,*)'SAT_CV_NODI_ipha,SAT_CV_NODj_ipha:',SAT_CV_NODI_ipha,SAT_CV_NODj_ipha
             W=(ABS_CV_NODI-ABSC)/tolfun(ABS_CV_NODI-ABS_CV_NODJ)
             SATC=W*SAT_CV_NODJ + (1.-W)*SAT_CV_NODI
             W=(SATC - SAT_CV_NODI)/(SAT_CV_NODJ - SAT_CV_NODI)
             !       ewrite(3,*)'1w=',w
             !       w=max(w,1.0-w) ! make any non-linear variation subject to upwinding
             !       w=1.-w ! this is correct
             w=0.5 + (w-0.5)*2.
             W=max(w,0.5)
             !       W=max(w,0.0)
             W=min(w,1.0)
             INCOME= W
          endif


          if(.true.) then

             ABSC=ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA) 
             ! ewrite(3,*)'ACTUAL ABSC,ABS_CV_NODi_ipha,ABS_CV_NODj_ipha:', &
             !     ABSC,ABS_CV_NODi_ipha,ABS_CV_NODj_ipha
             !ewrite(3,*)'--------------'

             W=(ABSC - ABS_CV_NODI_IPHA)/tolfun(ABS_CV_NODJ_IPHA-ABS_CV_NODI_IPHA)
             !!     W=(UUNDOTQ - NDOTQ)/(NDOTQ2 - NDOTQ)
             !!   w=1.-w

             !!   SATC=W*SAT_CV_NODJ_IPHA + (1.-W)*SAT_CV_NODI_IPHA
             !!        W=(SATC - SAT_CV_NODI_IPHA)/(SAT_CV_NODJ_IPHA - SAT_CV_NODI_IPHA)

             ! W=the upwind parameter (how much upwind fraction to use)
             w=0.5 + (w-0.5)*OVER_RELAX
             W=max(w,0.5)
             !       W=max(w,0.0)
             W=min(w,1.0)
             INCOME= W
             ewrite(3,*)'w,income:',w,income
          endif
       ELSE

          if(.false.) then
             SAT_CV_NODI=SAT_CV_NODI_ipha
             SAT_CV_NODj=SAT_CV_NODj_ipha
             if(iphase==2) then
                SAT_CV_NODI=1.-SAT_CV_NODI_ipha
                SAT_CV_NODj=1.-SAT_CV_NODj_ipha
             endif
             SAT_CV_NODI=min(1.0,max(0.0,SAT_CV_NODI))
             SAT_CV_NODj=min(1.0,max(0.0,SAT_CV_NODj))
             SATC=SAT_CV_NODI + 0.01*(SAT_CV_NODJ-SAT_CV_NODI) 
             SATC=min(1.0,max(0.0,SATC))
             CALL ABS3P(ABSC, 1.0, SATC, IPHASE)
             CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
             CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
             ewrite(3,*)'ABS_CV_NODI, SAT_CV_NODI, IPHASE:',ABS_CV_NODI, SAT_CV_NODI, IPHASE
             B=(ABS_CV_NODI-ABSC)/(SAT_CV_NODI-SATC)
             A=ABS_CV_NODI - B*SAT_CV_NODI

             ABSC=a+b*0.5*(SAT_CV_NODJ+SAT_CV_NODI)
             ewrite(3,*)'guessed valueS ABSC,ABS_CV_NODi=',ABSC,ABS_CV_NODi
             ewrite(3,*)'SAT_CV_NODI,SAT_CV_NODj:',SAT_CV_NODI,SAT_CV_NODj
             ewrite(3,*)'SAT_CV_NODI_ipha,SAT_CV_NODj_ipha:',SAT_CV_NODI_ipha,SAT_CV_NODj_ipha 
             W=(ABS_CV_NODJ-ABSC)/tolfun(ABS_CV_NODJ-ABS_CV_NODI)
             SATC=W*SAT_CV_NODI + (1.-W)*SAT_CV_NODJ 
             ! W=the upwind parameter (how much upwind fraction to use)
             W=(SATC - SAT_CV_NODJ)/(SAT_CV_NODI - SAT_CV_NODJ)
             !       ewrite(3,*)'2w=',w
             !       w=max(w,1.0-w) ! make any non-linear variation subject to upwinding
             !       w=1.-w ! this is correct
             w=0.5 + (w-0.5)*2.
             W=max(w,0.5)
             !       W=max(w,0.0)
             W=min(w,1.0)
             INCOME=1.0 - W 
          endif


          if(.true.) then

             ABSC=ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA) 
             ewrite(3,*)'ACTUAL ABSC,ABS_CV_NODi_ipha,ABS_CV_NODj_ipha:', &
                  ABSC,ABS_CV_NODi_ipha,ABS_CV_NODj_ipha
             ewrite(3,*)'--------------'

             W=(ABSC - ABS_CV_NODJ_IPHA)/tolfun(ABS_CV_NODI_IPHA - ABS_CV_NODJ_IPHA)
             !!   w=1.-w

             !!   SATC=W*SAT_CV_NODI_IPHA + (1.-W)*SAT_CV_NODJ_IPHA 
             !!        W=(SATC - SAT_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA - SAT_CV_NODJ_IPHA)

             ! W=the upwind parameter (how much upwind fraction to use)
             w=0.5 + (w-0.5)*OVER_RELAX
             W=max(w,0.5)
             !       W=max(w,0.0)
             W=min(w,1.0)
             INCOME=1.0 - W 
             ewrite(3,*)'w,income:',w,income
          endif
       ENDIF

    ENDIF !****************

    RETURN
  END SUBROUTINE FIND_OPT_INCOME_INTERP





  SUBROUTINE FIND_OPT_INCOME2_BET_ELE(INCOME, NDOTQ, NDOTQ2, &
       SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM, IPHASE) 
    ! calculate INCOME & UUNDOTQ for optimal upwinding.
    IMPLICIT NONE
    REAL INCOME, NDOTQ, NDOTQ2, SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM
    INTEGER IPHASE
    ! local variables...
    REAL W,ABS_CV_NODI,ABS_CV_NODJ, &
         ABS_FEM, MAX_SAT, MIN_SAT, SAT_FEM_LIM, UUNDOTQ
    REAL SAT_CV_NODI, SAT_CV_NODJ


    !  real tolfun

    SAT_CV_NODI=MIN(1.0,MAX(0.0,SAT_CV_NODI2))
    SAT_CV_NODJ=MIN(1.0,MAX(0.0,SAT_CV_NODJ2))

    !  ewrite(3,*)'********SAT_CV_NODI,SAT_CV_NODJ,NDOTQ,NDOTQ2,SAT_FEM=', &
    !                   SAT_CV_NODI,SAT_CV_NODJ,NDOTQ,NDOTQ2,SAT_FEM
    IF((ABS(SAT_CV_NODI - SAT_CV_NODJ) < 1.E-4) &
         .OR.(ABS(NDOTQ - NDOTQ2) < 1.E-6)) THEN
       UUNDOTQ=0.5*(NDOTQ+ NDOTQ2)
       INCOME=0.5
       !  ewrite(3,*)'-********SAT_CV_NODI,SAT_CV_NODJ,NDOTQ,NDOTQ2,INCOME=', &
       !                    SAT_CV_NODI,SAT_CV_NODJ,NDOTQ,NDOTQ2,INCOME
       RETURN
    ENDIF

    ! Find end pt absorptions...
    CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
    CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
    MAX_SAT = MAX(SAT_CV_NODI, SAT_CV_NODJ)
    MIN_SAT = MIN(SAT_CV_NODI, SAT_CV_NODJ)
    SAT_FEM_LIM = MAX(MIN(SAT_FEM,MAX_SAT),MIN_SAT)
    CALL ABS3P(ABS_FEM, 1.0, SAT_FEM_LIM, IPHASE)

    !!  UUNDOTQ = ABS_FEM*(NDOTQ/TOLFUN(ABS_CV_NODI) + NDOTQ2/TOLFUN(ABS_CV_NODJ))
    !!  UUNDOTQ = 0.5*ABS_FEM*(NDOTQ/TOLFUN(ABS_CV_NODI) + NDOTQ2/TOLFUN(ABS_CV_NODJ))
    UUNDOTQ = 0.5*(1./TOLFUN(ABS_FEM))*(NDOTQ*ABS_CV_NODI + NDOTQ2*ABS_CV_NODJ)

    !  ewrite(3,*)'UUNDOTQ,ABS_FEM,NDOTQ,ABS_CV_NODI,NDOTQ2,ABS_CV_NODJ:', &
    !           UUNDOTQ,ABS_FEM,NDOTQ,ABS_CV_NODI,NDOTQ2,ABS_CV_NODJ

    IF(UUNDOTQ < 0.0) THEN
       W=(UUNDOTQ - NDOTQ)/TOLFUN(NDOTQ2 - NDOTQ)
       !     W=max(w,0.0)
       W=max(w,0.5)
       W=1.*(max(w,0.5)-0.5)+0.5
       !       W=max(w,0.0)
       W=min(w,1.0)
       INCOME= W
    ELSE
       W=(UUNDOTQ - NDOTQ2)/TOLFUN(NDOTQ - NDOTQ2)
       !     W=max(w,0.0)
       W=max(w,0.5)
       W=1.*(max(w,0.5)-0.5)+0.5
       !       W=max(w,0.0)
       W=min(w,1.0)
       INCOME= 1.0 - W
    ENDIF
    !  ewrite(3,*)'********upwind W=',W

    RETURN
  END SUBROUTINE FIND_OPT_INCOME2_BET_ELE




  SUBROUTINE FIND_OPT_INCOME2(INCOME, NDOTQ, NDOTQ2, &
       SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM, IPHASE) 
    ! calculate INCOME & UUNDOTQ for optimal upwinding.
    IMPLICIT NONE
    REAL INCOME, NDOTQ, NDOTQ2, SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM
    INTEGER IPHASE
    ! local variables...
    REAL W,ABS_CV_NODI,ABS_CV_NODJ, &
         ABS_FEM, MAX_SAT, MIN_SAT, SAT_FEM_LIM, UUNDOTQ, SAT_CV_NODI, SAT_CV_NODJ


    !  real tolfun

    SAT_CV_NODI=MIN(1.0,MAX(0.0,SAT_CV_NODI2))
    SAT_CV_NODJ=MIN(1.0,MAX(0.0,SAT_CV_NODJ2))

    IF((ABS(SAT_CV_NODI - SAT_CV_NODJ) < 1.E-4) &
         .OR.(ABS(NDOTQ - NDOTQ2) < 1.E-6)) THEN
       UUNDOTQ=0.5*(NDOTQ+ NDOTQ2)
       INCOME=0.5
       RETURN
    ENDIF

    ! Find end pt absorptions...
    CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
    CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
    MAX_SAT = MAX(SAT_CV_NODI, SAT_CV_NODJ)
    MIN_SAT = MIN(SAT_CV_NODI, SAT_CV_NODJ)
    SAT_FEM_LIM = MAX(MIN(SAT_FEM,MAX_SAT),MIN_SAT)
    CALL ABS3P(ABS_FEM, 1.0, SAT_FEM_LIM, IPHASE)

    !!  ABS_FEM=2./(1./TOLFUN(ABS_CV_NODI) + 1./TOLFUN(ABS_CV_NODJ))

    !!  UUNDOTQ = ABS_FEM*(NDOTQ/TOLFUN(ABS_CV_NODI) + NDOTQ2/TOLFUN(ABS_CV_NODJ))
    !!  UUNDOTQ = 0.5*ABS_FEM*(NDOTQ/TOLFUN(ABS_CV_NODI) + NDOTQ2/TOLFUN(ABS_CV_NODJ))
    UUNDOTQ = 0.5*(1./tolfun(ABS_FEM))*(NDOTQ*ABS_CV_NODI + NDOTQ2*ABS_CV_NODJ)

    IF(UUNDOTQ < 0.0) THEN
       W=(UUNDOTQ - NDOTQ)/(NDOTQ2 - NDOTQ)
       !!     w=1.-w
       W=max(w,0.5)
       ! W=1.0*(max(w,0.5)-0.5)+0.5
       W=2.5*(max(w,0.5)-0.5)+0.5
       !     W=1.*(max(w,0.5)-0.5)+0.5
       !       W=max(w,0.0)
       W=min(w,1.0)
       !     w=0.5
       INCOME= W
    ELSE
       W=(UUNDOTQ - NDOTQ2)/(NDOTQ - NDOTQ2)
       !!     w=1.-w
       W=max(w,0.5)
       !  W=1.0*(max(w,0.5)-0.5)+0.5
       W=2.5*(max(w,0.5)-0.5)+0.5
       !     W=1.*(max(w,0.5)-0.5)+0.5
       !       W=max(w,0.0)
       W=min(w,1.0)
       !     w=0.5
       INCOME= 1.0 - W
    ENDIF

    RETURN
  END SUBROUTINE FIND_OPT_INCOME2




  SUBROUTINE FIND_OPT_INCOME3(INCOME, NDOTQ, NDOTQ2, &
       SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM, IPHASE) 
    ! calculate INCOME & UUNDOTQ for optimal upwinding.
    IMPLICIT NONE
    REAL INCOME, NDOTQ, NDOTQ2, SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM
    INTEGER IPHASE
    ! local variables...
    REAL W, &
         UUNDOTQ, &
         W1,W2


    !  real tolfun

    CALL FIND_OPT_INCOMEW(INCOME, NDOTQ, NDOTQ2, &
         SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM, IPHASE,W1) 
    CALL FIND_OPT_INCOMEW(INCOME, NDOTQ, NDOTQ2, &
         SAT_CV_NODI2, SAT_CV_NODJ2, 0.5*(SAT_CV_NODI2+SAT_CV_NODJ2), IPHASE,W2) 

    !       W=MAX(W1,W2)
    W=W2

    IF(UUNDOTQ < 0.0) THEN
       INCOME= W
    ELSE
       INCOME= 1.0 - W
    ENDIF

    RETURN
  END SUBROUTINE FIND_OPT_INCOME3




  SUBROUTINE FIND_OPT_INCOMEW(INCOME, NDOTQ, NDOTQ2, &
       SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM, IPHASE,W) 
    ! calculate INCOME & UUNDOTQ for optimal upwinding.
    IMPLICIT NONE
    REAL INCOME, NDOTQ, NDOTQ2, SAT_CV_NODI2, SAT_CV_NODJ2, SAT_FEM,W
    INTEGER IPHASE
    ! local variables...
    REAL ABS_CV_NODI,ABS_CV_NODJ, &
         ABS_FEM, MAX_SAT, MIN_SAT, SAT_FEM_LIM, UUNDOTQ, SAT_CV_NODI, SAT_CV_NODJ


    !  real tolfun

    SAT_CV_NODI=MIN(1.0,MAX(0.0,SAT_CV_NODI2))
    SAT_CV_NODJ=MIN(1.0,MAX(0.0,SAT_CV_NODJ2))

    IF((ABS(SAT_CV_NODI - SAT_CV_NODJ) < 1.E-4) &
         .OR.(ABS(NDOTQ - NDOTQ2) < 1.E-6)) THEN
       UUNDOTQ=0.5*(NDOTQ+ NDOTQ2)
       INCOME=0.5
       RETURN
    ENDIF

    ! Find end pt absorptions...
    CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
    CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
    MAX_SAT = MAX(SAT_CV_NODI, SAT_CV_NODJ)
    MIN_SAT = MIN(SAT_CV_NODI, SAT_CV_NODJ)
    SAT_FEM_LIM = MAX(MIN(SAT_FEM,MAX_SAT),MIN_SAT)
    CALL ABS3P(ABS_FEM, 1.0, SAT_FEM_LIM, IPHASE)

    !!  UUNDOTQ = ABS_FEM*(NDOTQ/TOLFUN(ABS_CV_NODI) + NDOTQ2/TOLFUN(ABS_CV_NODJ))
    UUNDOTQ = 0.5*ABS_FEM*(NDOTQ/TOLFUN(ABS_CV_NODI) + NDOTQ2/TOLFUN(ABS_CV_NODJ))

    IF(UUNDOTQ < 0.0) THEN
       W=(UUNDOTQ - NDOTQ)/(NDOTQ2 - NDOTQ)
       W=max(w,0.5)
       !     W=12.*(max(w,0.5)-0.5)+0.5
       !       W=max(w,0.0)
       W=min(w,1.0)
       INCOME= W
    ELSE
       W=(UUNDOTQ - NDOTQ2)/(NDOTQ - NDOTQ2)
       W=max(w,0.5)
       !     W=12.*(max(w,0.5)-0.5)+0.5
       !       W=max(w,0.0)
       W=min(w,1.0)
       INCOME= 1.0 - W
    ENDIF

    RETURN
  END SUBROUTINE FIND_OPT_INCOMEW




  SUBROUTINE FIND_OPT_INCOME(INCOME, NDOTQ, SAT_CV_NODI, SAT_CV_NODJ, IPHASE) 
    IMPLICIT NONE
    REAL INCOME, NDOTQ, SAT_CV_NODI, SAT_CV_NODJ
    INTEGER IPHASE
    ! local variables...
    REAL ABS_MID,ABS1,ABS2,ABSC,SAT1,SAT2,SATC,W,ABS_CV_NODI,ABS_CV_NODJ,A,B
    INTEGER NMAX,COUNT
    LOGICAL GRADNEG
    !!  real tolfun

    IF(ABS(SAT_CV_NODI - SAT_CV_NODJ) < 1.E-4) THEN
       INCOME=0.5
       RETURN
    ENDIF

    SAT1=MIN(SAT_CV_NODI,SAT_CV_NODJ)
    SAT2=MAX(SAT_CV_NODI,SAT_CV_NODJ)
    SATC=0.5*(SAT1+SAT2)

    ! Find end pt absorptions...
    CALL ABS3P(ABS1, 1.0, SAT1, IPHASE)
    CALL ABS3P(ABS2, 1.0, SAT2, IPHASE)


    ABS_MID=0.5*(ABS1+ABS2)

    IF( (ABS2-ABS1) > 0.0) THEN
       GRADNEG=.TRUE.
    ELSE
       GRADNEG=.FALSE.
    ENDIF
    !  ewrite(3,*)'ABS1,ABS2,ABS_MID,GRADNEG:',ABS1,ABS2,ABS_MID,GRADNEG

    ! perform a binary search for the best point.
    if(.true.) then
       ! assume a linear variation from downwide node...
       !    IF(NDOTQ < 0) THEN 
       if(.false.) then
          IF(NDOTQ > 0) THEN 
             SATC=SAT_CV_NODI + 0.01*(SAT_CV_NODJ-SAT_CV_NODI) 
             CALL ABS3P(ABSC, 1.0, SATC, IPHASE)
             CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
             B=(ABS_CV_NODI-ABSC)/(SAT_CV_NODI-SATC)
             A=ABS_CV_NODI - B*SAT_CV_NODI
          ELSE
             SATC=SAT_CV_NODJ + 0.01*(SAT_CV_NODI-SAT_CV_NODJ) 
             CALL ABS3P(ABSC, 1.0, SATC, IPHASE)
             CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
             B=(ABS_CV_NODJ-ABSC)/(SAT_CV_NODJ-SATC)
             A=ABS_CV_NODJ - B*SAT_CV_NODJ
          ENDIF


          SATC=(ABS_MID-A)/TOLFUN(B)
       else
          IF(NDOTQ < 0) THEN 
             SATC=SAT_CV_NODJ + 0.01*(SAT_CV_NODI-SAT_CV_NODJ) 
             CALL ABS3P(ABSC, 1.0, SATC, IPHASE)
             CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
             B=(ABS_CV_NODJ-ABSC)/(SAT_CV_NODJ-SATC)
             A=ABS_CV_NODJ - B*SAT_CV_NODJ

             ABSC=a+b*0.5*(SAT_CV_NODJ+SAT_CV_NODI)
             !! CALL ABS3P(ABSC, 1.0, 0.5*(SAT_CV_NODJ+SAT_CV_NODI), IPHASE)
             CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
             CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
             !        SATC=(ABS_CV_NODI-ABSC)/tolfun(ABS_CV_NODI-ABS_CV_NODJ) 
             W=(ABS_CV_NODI-ABSC)/tolfun(ABS_CV_NODI-ABS_CV_NODJ)
             SATC=W*SAT_CV_NODJ + (1.-W)*SAT_CV_NODI
             ! ewrite(3,*)'--satc=',satc 
          ELSE
             SATC=SAT_CV_NODI + 0.01*(SAT_CV_NODJ-SAT_CV_NODI) 
             CALL ABS3P(ABSC, 1.0, SATC, IPHASE)
             CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
             B=(ABS_CV_NODI-ABSC)/(SAT_CV_NODI-SATC)
             A=ABS_CV_NODI - B*SAT_CV_NODI

             ABSC=a+b*0.5*(SAT_CV_NODJ+SAT_CV_NODI)
             !! CALL ABS3P(ABSC, 1.0, 0.5*(SAT_CV_NODJ+SAT_CV_NODI), IPHASE)
             CALL ABS3P(ABS_CV_NODI, 1.0, SAT_CV_NODI, IPHASE)
             CALL ABS3P(ABS_CV_NODJ, 1.0, SAT_CV_NODJ, IPHASE)
             !        SATC=(ABS_CV_NODJ-ABSC)/tolfun(ABS_CV_NODJ-ABS_CV_NODI) 
             W=(ABS_CV_NODJ-ABSC)/tolfun(ABS_CV_NODJ-ABS_CV_NODI)
             SATC=W*SAT_CV_NODI + (1.-W)*SAT_CV_NODJ 
          ENDIF
          !    ewrite(3,*)'NDOTQ,ABS_CV_NODI,ABS_CV_NODJ:',NDOTQ,ABS_CV_NODI,ABS_CV_NODJ
          !    ewrite(3,*)'SAT_CV_NODI,SAT_CV_NODJ,satc:',SAT_CV_NODI,SAT_CV_NODJ,satc
          !    ewrite(3,*)'ABSC:',ABSC
       endif
    else
       NMAX=7
       COUNT = 1
       Loop_While: DO WHILE  ( COUNT <= NMAX )
          ! Calculate ABSC from SATC
          CALL ABS3P(ABSC, 1.0, SATC, IPHASE)

          IF(GRADNEG) THEN
             IF(ABSC <= ABS_MID) THEN
                SAT1 = SATC
             ELSE
                SAT2 = SATC
             ENDIF
          ELSE
             IF(ABSC >= ABS_MID) THEN
                SAT1 = SATC
             ELSE
                SAT2 = SATC
             ENDIF
          ENDIF
          SATC=0.5*(SAT1+SAT2)

          COUNT = COUNT + 1

       END DO Loop_While
    endif
    !    EWRITE(3,*)'SATC,ABSC=',SATC,ABSC

    IF(NDOTQ < 0) THEN 
       ! W=the upwind parameter (how much upwind fraction to use)
       W=(SATC - SAT_CV_NODI)/(SAT_CV_NODJ - SAT_CV_NODI)
       !       ewrite(3,*)'1w=',w
       !       w=max(w,1.0-w) ! make any non-linear variation subject to upwinding
       !       w=1.-w ! this is correct
       w=0.5 + (w-0.5)*2.
       W=max(w,0.5)
       !       W=max(w,0.0)
       W=min(w,1.0)
       INCOME= W
    ELSE
       ! W=the upwind parameter (how much upwind fraction to use)
       W=(SATC - SAT_CV_NODJ)/(SAT_CV_NODI - SAT_CV_NODJ)
       !       ewrite(3,*)'2w=',w
       !       w=max(w,1.0-w) ! make any non-linear variation subject to upwinding
       !       w=1.-w ! this is correct
       w=0.5 + (w-0.5)*2.
       W=max(w,0.5)
       !       W=max(w,0.0)
       W=min(w,1.0)
       INCOME=1.0 - W 
    ENDIF

    RETURN
  END SUBROUTINE FIND_OPT_INCOME





  SUBROUTINE ABS3P( ABS, INV_PERM, SAT, IPHASE )
    ! This will need to be changed as it is concerned w a particular case of the BL
    ! Eqn --  should be a function from the ABSORP subrt
    IMPLICIT NONE
    REAL, intent( inout ) :: ABS
    REAL, intent( in ) :: INV_PERM, SAT
    INTEGER, intent( in ) ::  IPHASE
    ! Local variables...
    REAL  :: VISC1, VISC2, S_GC, S_OR, REF_MOBILITY, KR1, KR2, KR, VISC, &
         SATURATION, ABS_SUM, SAT2

    VISC1 = 1.0
    S_GC = 0.1
    S_OR = 0.3
    REF_MOBILITY=10.0

    SATURATION=SAT
    IF( IPHASE == 2 ) SATURATION = 1. - SAT

    IF( SAT < S_GC ) THEN
       KR1 = 0.0
    ELSE IF( SAT > 1. -S_OR ) THEN
       KR1 = 1.0
    ELSE
       KR1 = ( ( SAT - S_GC) / ( 1. - S_GC - S_OR ))**2
    ENDIF

    SAT2=1.0-SAT
    IF( SAT2 < S_OR ) THEN
       KR2 = 0.0
    ELSEIF( SAT2 > 1. - S_GC ) THEN
       KR2 = 1.0
    ELSE
       KR2 = ( ( SAT2 - S_OR ) / ( 1. - S_GC - S_OR ))**2
    ENDIF
    VISC2=REF_MOBILITY

    IF(IPHASE==1) THEN
       KR=KR1
       VISC=VISC1
    ELSE
       KR=KR2
       VISC=VISC2
    ENDIF

    ABS_SUM = KR / MAX(1.e-6,VISC*max(0.01,SATURATION))

    ABS = INV_PERM / MAX( 1.e-6, ABS_SUM )

    if(iphase==1) then
       ABS =  min( 1.e+4, ABS)
       if(saturation < 0.1) then
          ABS = (1. + max(100.*(0.1-saturation),0.0)) * ABS
       endif
    else
       ABS = min( 1.e+5, ABS)
       if(saturation < 0.3) then
          ABS = (1. + max(100.*(0.3-saturation),0.0)) * ABS
       endif
    endif
    RETURN     
  END SUBROUTINE ABS3P



end module cv_advection
      
