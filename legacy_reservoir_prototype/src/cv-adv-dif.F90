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

    use solvers_module
    use spud
    use global_parameters, only: option_path_len
    use futils, only: int2str

    use shape_functions
    use matrix_operations

  contains

    SUBROUTINE CV_ASSEMB( state, &
         CV_RHS, &
         NCOLACV, ACV, FINACV, COLACV, MIDACV, &
         NCOLCT, CT, DIAG_SCALE_PRES, CT_RHS, FINDCT, COLCT, &
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         CV_ELE_TYPE,  &
         NPHASE,  &
         CV_NLOC, U_NLOC, X_NLOC, &
         CV_NDGLN, X_NDGLN, U_NDGLN, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
         X, Y, Z, U, V, W, &
         NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         T, TOLD, DEN, DENOLD, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, &
         CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, SECOND_THETA, CV_BETA, &
         SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_T_BC_ROB1, SUF_T_BC_ROB2,  &
         WIC_T_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, CV_P, &
         SOURCT, ABSORBT, VOLFRA_PORE, & 
         NDIM, GETCV_DISC, GETCT, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         T_FEMT, DEN_FEMT, &
         IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, THETA_GDIFF, &
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         FINDCMC, COLCMC, NCOLCMC, MASS_MN_PRES, THERMAL, &
         MASS_ELE_TRANSP, &
         option_path_spatial_discretisation )

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

      ! Inputs/Outputs
      IMPLICIT NONE
      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_DISOPT, CV_DG_VEL_INT_OPT, NDIM, &
           NCOLM, XU_NLOC, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NCOLCMC
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
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC
      REAL, DIMENSION( NCOLCMC ), intent( inout ) :: MASS_MN_PRES

      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: T, TOLD, DEN, DENOLD
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION( TOTELE * IGOT_THETA_FLUX, CV_NLOC, SCVNGI_THETA, NPHASE ), &
           intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: TDIFFUSION
      REAL, intent( in ) :: DT, CV_THETA, SECOND_THETA, CV_BETA
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE, NDIM ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE * IGOT_T2 ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DERIV
      REAL, DIMENSION( CV_NONODS ), intent( in ) :: CV_P
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SOURCT
      REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE ), intent( in ) :: ABSORBT
      REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE 
      LOGICAL, intent( in ) :: GETCV_DISC, GETCT, GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDM
      INTEGER, DIMENSION( NCOLM ), intent( in ) :: COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( in ) :: MIDM
      INTEGER, DIMENSION( TOTELE + 1 ), intent( in ) :: FINELE
      INTEGER, DIMENSION( NCOLELE ), intent( in ) :: COLELE
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: T_FEMT, DEN_FEMT
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( CV_NONODS ), intent( inout ) :: MEAN_PORE_CV
      REAL, DIMENSION( TOTELE ), intent( inout ) :: MASS_ELE_TRANSP
      character( len = * ), intent( in ), optional :: option_path_spatial_discretisation
!      character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation


      ! Local variables 
      LOGICAL, PARAMETER :: INCLUDE_PORE_VOL_IN_DERIV = .FALSE.
      INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
           WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
           WIC_U_BC_DIRICHLET = 1
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE,LOG_ON_BOUND
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE
      INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, &
           CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
           JCOUNT_KLOC, JCOUNT_KLOC2, COLGPTS, CV_SLOC2LOC, U_SLOC2LOC, &
           TMAX_NOD, TMIN_NOD, TOLDMAX_NOD, &
           TOLDMIN_NOD, DENMAX_NOD, DENMIN_NOD, DENOLDMAX_NOD, DENOLDMIN_NOD, &
           T2MAX_NOD, T2MIN_NOD, T2OLDMAX_NOD, T2OLDMIN_NOD
      INTEGER, DIMENSION( : , : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, &
           FACE_ELE, CV_NEILOC
      REAL, DIMENSION( : ), allocatable :: CVWEIGHT, CVWEIGHT_SHORT, SCVFEWEIGH, SBCVFEWEIGH, &
           TMAX, TMIN, TOLDMAX, &
           TOLDMIN, DENMAX, DENMIN, DENOLDMAX, DENOLDMIN, &
           TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, &
           TOLDMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC, &
           CVNORMX, &
           CVNORMY, CVNORMZ, SCVRA, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
           FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
           SCVDETWEI, SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
           SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
           T2MAX, T2MIN, T2OLDMAX, &
           T2OLDMIN, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, &
           T2OLDMIN_2ND_MC, &
           UP_WIND_NOD, DU, DV, DW
      REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
           UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SCVFENX, SCVFENY, SCVFENZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ, &
           DUMMY_ZERO_NDIM_NDIM
      REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE,  &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
      REAL, DIMENSION( : , :, : ), allocatable :: INV_JAC

      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE,   &
           NCOLGPTS, &
           CV_SILOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, &
           CV_NODI, CV_NODI_IPHA, CV_NODI_JPHA, U_NODK, TIMOPT, &
           JCOUNT_IPHA, IMID_IPHA, &
           NFACE, X_NODI,  &
           CV_INOD, MAT_NODI, FACE_ITS, NFACE_ITS, &
           CVNOD, XNOD, NSMALL_COLM, COUNT2, NOD
      !        ===>  REALS  <===
      REAL :: NDOTQ, NDOTQOLD,  &
           INCOME, INCOMEOLD, HDC, FVT, FVTOLD, FVT2, FVT2OLD, &
           FVD, FVDOLD, LIMT, LIMTOLD, LIMT2, LIMT2OLD,&
           LIMD, LIMDOLD, FTHETA, VTHETA, &
           LIMDT, LIMDTOLD, LIMDTT2, LIMDTT2OLD, &
           FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
           TMID, TOLDMID, &
           DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, BCZERO, ROBIN1, ROBIN2, &
           RSUM, &
           SUM_LIMT, SUM_LIMTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD, THERM_FTHETA, &
           W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW, VOLUME

      REAL, PARAMETER :: W_SUM_ONE = 1.

      integer :: cv_inod_ipha, IGETCT, U_NODK_IPHA, IANISOLIM
      logical :: Have_Temperature_Fields, Have_VolumeFraction_Fields, Have_Components_Fields
      ! Functions...
      !REAL :: R2NORM, FACE_THETA  
      !        ===>  LOGICALS  <===
      LOGICAL :: GETMAT, LIMIT_USE_2ND, &
           D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
           NORMALISE, SUM2ONE, GET_GTHETA, QUAD_OVER_WHOLE_ELE

      character( len = option_path_len ) :: option_path, option_path2, path_temp, path_volf, &
           path_comp, path_spatial_discretisation

      integer, dimension(:), allocatable :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      real, dimension(:), allocatable :: TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, &
              DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT, sol
      INTEGER :: IDUM(1)
      REAL :: RDUM(1)

      IDUM = 0
      RDUM = 0.

      ewrite(3,*) 'In CV_ASSEMB'

      !ewrite(3,*) 'CV_P', CV_P
      !ewrite(3,*) 'DEN', DEN
      !ewrite(3,*) 'DENOLD', DENOLD
      !ewrite(3,*) 'MEAN_PORE_CV', MEAN_PORE_CV

      Have_Temperature_Fields = .false.
      Have_VolumeFraction_Fields = .false.
      Have_Components_Fields = .false.

      path_temp = '/material_phase[0]/scalar_field::Temperature'
      path_volf = '/material_phase[0]/scalar_field::PhaseVolumeFraction'
      path_comp = '/material_phase[' // int2str( nphase ) // ']/scalar_field::ComponentMassFractionPhase1'
      path_spatial_discretisation = '/prognostic/spatial_discretisation/' // &
           'control_volumes/face_value/limit_face_value'

      if( have_option( trim( path_temp ) ) ) Have_Temperature_Fields = .true.
      if( have_option( trim( path_volf ) ) ) Have_VolumeFraction_Fields = .true.
      if( have_option( trim( path_comp ) ) ) Have_Components_Fields = .true.

      limit_use_2nd = .false.

      if ( present( option_path_spatial_discretisation) ) then
         if ( Have_Temperature_Fields .or. Have_VolumeFraction_Fields .or. &
              Have_Components_Fields ) &
              option_path = trim( option_path_spatial_discretisation ) // &
              trim( path_spatial_discretisation )

         if ( have_option( trim( option_path ) // '/limiter::Extrema' ) ) &
              limit_use_2nd = .true.
      end if

      GOT_DIFFUS = ( R2NORM( TDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0 )

      ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS:', &
           CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS

      !ewrite(3,*)'tdiffusion=', tdiffusion
      !ewrite(3,*)'suf_t_bc=', suf_t_bc
      !ewrite(3,*)'wic_t_bc=', wic_t_bc

      ndotq = 0. ; ndotqold = 0.

      QUAD_OVER_WHOLE_ELE=.FALSE. 
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

      ! Allocate memory for the control volume surface shape functions, etc.
      ALLOCATE( JCOUNT_KLOC(  U_NLOC )) ; jcount_kloc = 0
      ALLOCATE( JCOUNT_KLOC2(  U_NLOC )) ; jcount_kloc2 = 0

      ALLOCATE( CVNORMX( SCVNGI ))
      ALLOCATE( CVNORMY( SCVNGI ))
      ALLOCATE( CVNORMZ( SCVNGI ))
      ALLOCATE( SCVRA( SCVNGI ))
      ALLOCATE( COLGPTS( CV_NLOC * SCVNGI )) !The size of this vector is over-estimated
      ALLOCATE( FINDGPTS( CV_NLOC + 1 ))
      ALLOCATE( SNDOTQ( SCVNGI ))
      ALLOCATE( SNDOTQOLD( SCVNGI ))
      ALLOCATE( CV_ON_FACE( CV_NLOC, SCVNGI ))
      ALLOCATE( CVFEM_ON_FACE( CV_NLOC, SCVNGI ))
      ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))
      ALLOCATE( UFEM_ON_FACE( U_NLOC, SCVNGI ))
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
      ALLOCATE( SCVFENX( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENY( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFENZ( CV_NLOC, SCVNGI ))
      ALLOCATE( SCVFEWEIGH( SCVNGI ))

      ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

      ALLOCATE( SCVDETWEI( SCVNGI )) ; SCVDETWEI = 0.
      ALLOCATE( SRA( SCVNGI ))
      ALLOCATE( LOG_ON_BOUND(CV_NONODS))

      ALLOCATE( SBCVN( CV_SNLOC, SBCVNGI ))
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
      ALLOCATE( DUMMY_ZERO_NDIM_NDIM(NDIM,NDIM)) 
      DUMMY_ZERO_NDIM_NDIM=0.0

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

      ALLOCATE( INV_JAC( NDIM, NDIM, SCVNGI ) )

      UP_WIND_NOD = 0.0

      ALLOCATE( ONE_PORE( TOTELE ))
      IF( INCLUDE_PORE_VOL_IN_DERIV ) THEN
         ! solve for actual velocity
         ONE_PORE = VOLFRA_PORE
      ELSE
         ! solve for vel=porosity*actual velocity
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

      ALLOCATE(      TMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(      TMAX_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(   TOLDMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(   TOLDMAX_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(      T2MIN_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(      T2MAX_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(   T2OLDMIN_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(   T2OLDMAX_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(    DENMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(    DENMAX_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE( DENOLDMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE( DENOLDMAX_2ND_MC( CV_NONODS * NPHASE) )

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
      !ewrite(3,*)'in ASSEMB_FORCE_CTY',NCOLGPTS 
      !ewrite(3,*)'in ASSEMB_FORCE_CTY, COLGPTS',size(COLGPTS), COLGPTS
      !ewrite(3,*)'in ASSEMB_FORCE_CTY, FINDGPTS',size(FINDGPTS), FINDGPTS


      !     ======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========

      CALL CV_FEM_SHAPE_FUNS( &
                                ! Volume shape functions...
           NDIM, CV_ELE_TYPE,  & 
           CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
           CVWEIGHT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           UFEN, UFENLX, UFENLY, UFENLZ, &
                                ! Surface of each CV shape functions...
           SCVNGI, CV_NEILOC, CV_ON_FACE, CVFEM_ON_FACE, &  
           SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
           SCVFENLX, SCVFENLY, SCVFENLZ,  &
           SUFEN, SUFENSLX, SUFENSLY,  &
           SUFENLX, SUFENLY, SUFENLZ,  &
                                ! Surface element shape funcs...
           U_ON_FACE, UFEM_ON_FACE, NFACE, & 
           SBCVNGI, SBCVN, SBCVFEN,SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
           CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
           FINDGPTS, COLGPTS, NCOLGPTS, &
           SELE_OVERLAP_SCALE, QUAD_OVER_WHOLE_ELE )  

      !ewrite(3,*)'back in cv-adv-dif'
      !do iphase = 1, nphase
      !   ewrite(3,*) 'Phase', iphase, ',', 'suf_t_bc:', &
      !        suf_t_bc( ( iphase - 1 ) * cv_snloc * stotel + 1 : iphase * cv_snloc * stotel )
      !   ewrite(3,*)''
      !   ewrite(3,*)' Now suf_u_bc:', &
      !        suf_u_bc( ( iphase - 1 ) * u_snloc * stotel + 1 : iphase * u_snloc * stotel ) 
      !   ewrite(3,*)''
      !end do
      !do cv_iloc = 1, cv_nloc
      !ewrite(3,*)'iloc, cv_on_face:', cv_iloc, &
      !( cv_on_face( cv_iloc, gi ), gi = 1, scvngi )
      !ewrite(3,*)'iloc, cvfem_on_face:', cv_iloc, &
      !( cvfem_on_face( cv_iloc, gi ), gi = 1, scvngi )
      !end do


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
      IGETCT=0
      IF(GETCT) IGETCT=1

      CALL PROJ_CV_TO_FEM_4( state, &
           FEMT, FEMTOLD, FEMDEN, FEMDENOLD, T, TOLD, DEN, DENOLD, &
           IGOT_T2, T2,T2OLD, FEMT2,FEMT2OLD, &
           XC_CV, YC_CV, ZC_CV, MASS_CV, MASS_ELE,  &
           NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI_SHORT, CV_NLOC, CVN_SHORT, CVWEIGHT_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
           IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC )

!         femt=t
!         femtold=told
!         FEMT2=t2
!         FEMT2OLD=t2old


      MASS_ELE_TRANSP = MASS_ELE

      NORMALISE = .false.
      IF( NORMALISE ) THEN
         ! make sure the FEM representation sums to unity so we dont get surprising results...
         DO CV_INOD = 1, CV_NONODS
            RSUM = 0.0
            DO IPHASE = 1, NPHASE
               RSUM = RSUM + FEMT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS )
            END DO
            DO IPHASE = 1, NPHASE
               FEMT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMT( CV_INOD + ( IPHASE - 1 ) &
                    * CV_NONODS ) / RSUM
            END DO
         END DO
      ENDIF

      ! Calculate MEAN_PORE_CV   
      SUM_CV = 0.0
      MEAN_PORE_CV = 0.0
      DO ELE = 1, TOTELE
         DO CV_ILOC = 1, CV_NLOC
            CV_INOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
            SUM_CV( CV_INOD ) = SUM_CV( CV_INOD ) + MASS_ELE( ELE )
            MEAN_PORE_CV( CV_INOD ) = MEAN_PORE_CV( CV_INOD ) + & 
                 MASS_ELE( ELE ) * VOLFRA_PORE( ELE )
         END DO
      END DO

      MEAN_PORE_CV = MEAN_PORE_CV / SUM_CV

      !ewrite(3,*) 'MEAN_PORE_CV:', MEAN_PORE_CV

      !sum=0.0
      !do ele = 1, totele
      !   ewrite(3,*) ele, mass_ele(ele)
      !   sum=sum+mass_ele(ele)
      !end do
      !ewrite(3,*)'sum(mass_ele):',sum
      !stop 221

      ! For each node, find the largest and smallest value of T and 
      ! DENSITY for both the current and previous timestep, out of 
      ! the node value and all its surrounding nodes including Dirichlet b.c's.
      CALL SURRO_CV_MINMAX( TMAX, TMIN, TOLDMAX, TOLDMIN, DENMAX, DENMIN, DENOLDMAX, DENOLDMIN, &
           T2MAX, T2MIN, T2OLDMAX, T2OLDMIN, &
           TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, &
           DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, T2OLDMIN_2ND_MC, &
           LIMIT_USE_2ND, &
           T, TOLD, T2, T2OLD, DEN, DENOLD, IGOT_T2, NPHASE, CV_NONODS, NCOLACV, FINACV, COLACV, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC, SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
           WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
           WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
           T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
           DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD )


! **********ANISOTROPIC LIMITING...*******************

      IANISOLIM=0
      IF(CV_DISOPT.GE.8) IANISOLIM=1
      IF (IANISOLIM==0) THEN
         ALLOCATE(TUPWIND_MAT(1), TOLDUPWIND_MAT(1), DENUPWIND_MAT(1), DENOLDUPWIND_MAT(1))
         ALLOCATE(T2UPWIND_MAT(1), T2OLDUPWIND_MAT(1))
      ELSE ! IF (IANISOLIM==1) THEN
! Reduce matrix size...
      COUNT2=0
      DO NOD=1,CV_NONODS
         DO COUNT=FINACV(NOD),FINACV(NOD+1)-1
            IF(COLACV(COUNT).LE.CV_NONODS) COUNT2=COUNT2+1
         END DO
      END DO
      NSMALL_COLM=COUNT2+1

      ALLOCATE(SMALL_FINDRM(CV_NONODS+1),SMALL_COLM(NSMALL_COLM),SMALL_CENTRM(CV_NONODS))
      ALLOCATE(TUPWIND_MAT(NSMALL_COLM*NPHASE), TOLDUPWIND_MAT(NSMALL_COLM*NPHASE), &
           DENUPWIND_MAT(NSMALL_COLM*NPHASE), DENOLDUPWIND_MAT(NSMALL_COLM*NPHASE))
      ALLOCATE(T2UPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2), T2OLDUPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2))

      CALL CALC_ANISOTROP_LIM(&
! Caculate the upwind values stored in matrix form...
           T,TOLD,DEN,DENOLD,T2,T2OLD, &
           TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, &
           T2UPWIND_MAT, T2OLDUPWIND_MAT, &
! Store the upwind element for interpolation and its weights for 
! faster results...
           IGOT_T2,NPHASE,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
           SMALL_FINDRM,SMALL_CENTRM,SMALL_COLM,NSMALL_COLM, &
           X_NDGLN,X_NONODS,NDIM, &
           X,Y,Z, &
           FINACV,COLACV,NCOLACV) 
! endof IF (IANISOLIM==0) THEN ELSE...
       ENDIF

! **********...ANISOTROPIC LIMITING*******************


      ALLOCATE( FACE_ELE( NFACE, TOTELE ) ) ; FACE_ELE = 0
      ! Calculate FACE_ELE
      CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
           NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
           CV_SLOCLIST, X_NLOC, X_NDGLN )

      IF( GOT_DIFFUS ) THEN
         CALL DG_DERIVS( FEMT, FEMTOLD, &
              DTX_ELE, DTY_ELE, DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
              NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
              X_NDGLN, X_NLOC, X_NDGLN, &
              CV_NGI_SHORT, CV_NLOC, CVWEIGHT_SHORT, &
              CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
              CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
              X_NONODS, X, Y, Z,  &
              NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, WIC_T_BC, SUF_T_BC, &
              WIC_T_BC_DIRICHLET, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
              SBCVFEN, SBCVFENSLX, SBCVFENSLY)
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
      !ewrite(1,*)'8 midacv:', midacv( 1 : cv_nonods * nphase )
      !stop 9823

      ! Now we begin the loop over elements to assemble the advection terms
      ! into the matrix (ACV) and the RHS
      Loop_Elements: DO ELE = 1, TOTELE


         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR_INVJAC( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
              CV_NLOC, SCVNGI, &
              SCVFEN, SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH, SCVDETWEI, SCVRA, VOLUME, D1, D3, DCYL, &
              SCVFENX, SCVFENY, SCVFENZ, &
              NDIM, INV_JAC )


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

               Conditional_CheckingNeighbourhood: IF( CV_JLOC == -1 ) THEN

                  ! We are on the boundary or next to another element.  Determine CV_OTHER_LOC
                  ! CVFEM_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
                  ! Look for these nodes on the other elements.
                  CALL FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC,  &
                       MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI,  &
                       TOTELE, X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                       CV_SNLOC, CVFEM_ON_FACE, SCVNGI, GI, X_SHARE, X_NONODS, ELE, ELE2,  &
                       FINELE, COLELE, NCOLELE )

                  !ewrite(3,*)'================================================================================= '
                  !ewrite(3,*)' ele, cv_iloc, cv_nodi, gi, cv_jloc: ', ele, cv_iloc, cv_nodi, gi, cv_jloc
                  !ewrite(3,*)' ele2, integrat_at_gi:', ele2, integrat_at_gi
                  !ewrite(3,*)'================================================================================= '
                  !ewrite(3,*)'cv_other_loc:', cv_other_loc( 1 : cv_nloc )
                  !ewrite(3,*)'u_other_loc:', u_other_loc( 1 : u_nloc )
                  !ewrite(3,*)'mat_other_loc:', mat_other_loc( 1 : mat_nloc )
                  !ewrite(3,*)'INTEGRAT_AT_GI=', INTEGRAT_AT_GI
                  !ewrite(3,*)'================================================================================= '

                  IF(INTEGRAT_AT_GI) THEN
                     CV_JLOC = CV_OTHER_LOC( CV_ILOC )
                     SELE=0

                     IF( CV_JLOC == 0 ) THEN ! We are on the boundary of the domain
                        CV_JLOC = CV_ILOC
                        ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC
                        CALL CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, SCVNGI, U_SLOC2LOC, CV_SLOC2LOC, &
                             FACE_ELE, TOTELE, NFACE, CVFEM_ON_FACE, GI, &
                             CV_NONODS, LOG_ON_BOUND, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, STOTEL, &
                             CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN ) 
                        !EWRITE(3,*)'*****AFTER CALC_SELE SELE,CV_SILOC,CV_SNLOC:',SELE,CV_SILOC,CV_SNLOC
                     ENDIF
                     INTEGRAT_AT_GI=.NOT.((ELE==ELE2).AND.(SELE==0))
                  ENDIF

               ENDIF Conditional_CheckingNeighbourhood


               ! avoid indegrating across the middle of a CV on the boundaries of elements        
               Conditional_integration: IF(INTEGRAT_AT_GI) THEN

                  ! if necessary determine the derivatives between elements ELE and ELE2

                  ! Calculate the control volume normals at the Gauss pts.
                  CALL SCVDETNX( ELE,      GI,          &
                       X_NLOC,  SCVNGI,  TOTELE,  &
                       X_NDGLN,  X_NONODS,         &
                       SCVDETWEI, CVNORMX, CVNORMY, &  
                       CVNORMZ,  SCVFEN,     SCVFENSLX,   &   
                       SCVFENSLY, SCVFEWEIGH, XC_CV(CV_NODI),    & 
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
                  X_NODI = X_NDGLN(( ELE - 1 ) * X_NLOC  + CV_ILOC )

                  Conditional_GETCT1: IF( GETCT ) THEN ! Obtain the CV discretised CT equations plus RHS
                     !ewrite(3,*)'==================================================================='
                     !ewrite(3,*)'CV_NODI, CV_NODJ, ELE, ELE2, GI:', CV_NODI, CV_NODJ, ELE, ELE2, GI
                     !ewrite(3,*)'findct:',FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1
                     !ewrite(3,*)'colct:',colct( FINDCT( CV_NODI ) : FINDCT( CV_NODI + 1 ) - 1 )
                     !ewrite(3,*)'SCVDETWEI:', SCVDETWEI(GI)

                     DO U_KLOC = 1, U_NLOC
                        U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
                        JCOUNT = 0
                        DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1, 1
                           IF(COLCT( COUNT ) == U_NODK) JCOUNT = COUNT
                        END DO
                        JCOUNT_KLOC( U_KLOC ) = JCOUNT
                        !ewrite(3,*)' u_nodk, jcount1:', u_nodk, jcount
                     END DO

                     IF( ( ELE2 /= 0 ) .AND. ( ELE2 /= ELE ) ) THEN
                        DO U_KLOC = 1, U_NLOC
                           U_NODK = U_NDGLN(( ELE2 - 1 ) * U_NLOC + U_KLOC )
                           JCOUNT = 0
                           DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1, 1
                              IF(COLCT( COUNT ) == U_NODK) JCOUNT = COUNT
                           END DO
                           JCOUNT_KLOC2( U_KLOC ) = JCOUNT
                           !ewrite(3,*)' u_nodk, jcount2:', u_nodk, jcount
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
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2, T2OLD, FEMT2, FEMT2OLD, DEN, DENOLD, &
                                U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                                !NU, NV, NW, NUOLD, NVOLD, NWOLD, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT2, LIMT2OLD, FEMDGI, FEMT2GI, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                                T2MIN, T2MAX, T2OLDMIN, T2OLDMAX, T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC, T2OLDMIN_2ND_MC, T2MAX_2ND_MC, T2OLDMAX_2ND_MC, LIMIT_USE_2ND)
                        ELSE
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, &
                                U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                                !NU, NV, NW, NUOLD, NVOLD, NWOLD, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
                                TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC, LIMIT_USE_2ND)
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
                             T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, TOLDMIN_2ND_MC, T2MIN_2ND_MC, T2OLDMIN_2ND_MC, DENMIN_2ND_MC, DENOLDMIN_2ND_MC, &
                             TMAX_2ND_MC, TOLDMAX_2ND_MC, T2MAX_2ND_MC, T2OLDMAX_2ND_MC, DENMAX_2ND_MC, DENOLDMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, NDOTQOLD, DT, &
                             SCVFENX, SCVFENY, SCVFENZ, CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT )

                        SUM_LIMT    = SUM_LIMT    + LIMT
                        SUM_LIMTOLD = SUM_LIMTOLD + LIMTOLD

                     END DO Loop_IPHASE5
                  ENDIF Conditional_SUMLimiting
                  ! get the sum of limiting functions correct...************

                  Loop_IPHASE: DO IPHASE = 1, NPHASE

                     CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                     CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS

                     If_GOT_DIFFUS: IF( GOT_DIFFUS ) THEN
                        ! This sub caculates the effective diffusion 
                        ! coefficient DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
                        CALL DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX,  &
                             CV_NLOC, MAT_NLOC, CV_NONODS, NPHASE, TOTELE, MAT_NONODS, MAT_NDGLN, &
                             SCVFEN, SCVFEN, SCVNGI, GI, IPHASE, NDIM, TDIFFUSION, DUMMY_ZERO_NDIM_NDIM, &
                             HDC, &
                             T(CV_NODJ_IPHA), T(CV_NODI_IPHA), &
                             TOLD(CV_NODJ_IPHA), TOLD(CV_NODI_IPHA), &
                             ELE, ELE2, CVNORMX, CVNORMY, CVNORMZ,  &
                             DTX_ELE, DTY_ELE, DTZ_ELE,DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
                             SELE, STOTEL, WIC_T_BC, WIC_T_BC_DIRICHLET, CV_OTHER_LOC, MAT_OTHER_LOC )
                     ELSE ! IF(GOT_DIFFUS) THEN...
                        DIFF_COEF_DIVDX    = 0.0
                        DIFF_COEFOLD_DIVDX = 0.0
                     END IF If_GOT_DIFFUS

                     NFACE_ITS=1
                     !IF((CV_ELE_TYPE==2).AND.(CV_NONODS==TOTELE*CV_NLOC)) NFACE_ITS=2
                     DO FACE_ITS = 1, NFACE_ITS
                        ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
                        IF(IGOT_T2==1) THEN
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2, T2OLD, FEMT2, FEMT2OLD, DEN, DENOLD, &
                                U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                                !NU, NV, NW, NUOLD, NVOLD, NWOLD, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT2, LIMT2OLD, FEMDGI, FEMT2GI, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                                T2MIN, T2MAX, T2OLDMIN, T2OLDMAX, T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC, T2OLDMIN_2ND_MC, T2MAX_2ND_MC, T2OLDMAX_2ND_MC, LIMIT_USE_2ND)
                        ELSE
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, &
                                U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
                                !NU, NV, NW, NUOLD, NVOLD, NWOLD, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
                                TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC, LIMIT_USE_2ND)
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
                             T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, TOLDMIN_2ND_MC, T2MIN_2ND_MC, T2OLDMIN_2ND_MC, DENMIN_2ND_MC, DENOLDMIN_2ND_MC, &
                             TMAX_2ND_MC, TOLDMAX_2ND_MC, T2MAX_2ND_MC, T2OLDMAX_2ND_MC, DENMAX_2ND_MC, DENOLDMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, NDOTQOLD, DT, &
                             SCVFENX, SCVFENY, SCVFENZ, CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT )

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
                             T(CV_NODJ_IPHA)*DEN(CV_NODJ_IPHA), &
                             T(CV_NODI_IPHA)*DEN(CV_NODI_IPHA), &
                             NDOTQOLD, LIMDTT2OLD, DIFF_COEFOLD_DIVDX, &
                             TOLD(CV_NODJ_IPHA)*DENOLD(CV_NODJ_IPHA), &
                             TOLD(CV_NODI_IPHA)*DENOLD(CV_NODI_IPHA) )
                     ENDIF

                     FTHETA_T2         = FTHETA * LIMT2
                     ONE_M_FTHETA_T2OLD = (1.0-FTHETA) * LIMT2OLD

                     IF(IGOT_THETA_FLUX==1) THEN
                        IF ( GET_THETA_FLUX ) THEN
!                           THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)      =FTHETA*LIMT
!                           ONE_M_THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)=(1.0-FTHETA)*LIMTOLD
                           THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)      =FTHETA*LIMDT/DEN(CV_NODI_IPHA)
                           ONE_M_THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)=(1.0-FTHETA)*LIMDTOLD/DEN(CV_NODI_IPHA)
                        ENDIF
                        IF ( USE_THETA_FLUX ) THEN
                           FTHETA_T2         =THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)
                           ONE_M_FTHETA_T2OLD=ONE_M_THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)
                        ENDIF
                     ENDIF

                     !ewrite(3,*) 'IGOT_THETA_FLUX,GET_THETA_FLUX,USE_THETA_FLUX:', &
                     !     IGOT_THETA_FLUX,GET_THETA_FLUX,USE_THETA_FLUX
                     !ewrite(3,*) 'FTHETA_T2,ONE_M_FTHETA_T2OLD, FTHETA_T2, LIMT2:',FTHETA_T2,ONE_M_FTHETA_T2OLD, FTHETA, LIMT2
                     !ewrite(3,*) 'ndotq, ndotqold', ndotq, ndotqold
                     !stop 2821

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
                             JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC, U_NDGLN, U, V, W,  &
                             !JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC, U_NDGLN, NU, NV, NW,  &
                             SUFEN, SCVDETWEI, CVNORMX, CVNORMY, CVNORMZ, DEN, CV_NODI, CV_NODI_IPHA, &
                             UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
                             UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                             NDOTQNEW, NDOTQOLD, LIMDT, LIMDTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD )


                     ENDIF Conditional_GETCT2

                     Conditional_GETCV_DISC: IF(GETCV_DISC) THEN 
                        ! Obtain the CV discretised advection/diffusion equations

                        IF(GETMAT) THEN
                           ! - Calculate the integration of the limited, high-order flux over a face  
                           ! Conservative discretisation. The matrix (PIVOT ON LOW ORDER SOLN)
                           IF( ( CV_NODI_IPHA /= CV_NODJ_IPHA ) .AND. ( CV_NODJ_IPHA /= 0 ) ) THEN
                              DO COUNT = FINACV( CV_NODI_IPHA ), FINACV( CV_NODI_IPHA + 1 ) - 1
                                 IF( COLACV( COUNT ) == CV_NODJ_IPHA )  JCOUNT_IPHA = COUNT
                              END DO

                              ACV( JCOUNT_IPHA ) =  ACV( JCOUNT_IPHA ) &
                                   + SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW * INCOME * LIMD  & ! advection
                                   - FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX ! Diffusion contribution
                              IF(GET_GTHETA) THEN
                                 THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                      + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX * T(CV_NODJ_IPHA) ! Diffusion contribution
                              ENDIF
                           ELSE IF(SELE/=0) THEN
                              IF(WIC_T_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRICHLET) THEN
                                 CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA )  &
                                      + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &
                                      * SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC+(IPHASE-1)*STOTEL*CV_SNLOC)
                                 IF(GET_GTHETA) THEN
                                    THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                         + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX &
                                         * SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC+(IPHASE-1)*STOTEL*CV_SNLOC)
                                 ENDIF
                              ENDIF
                           ENDIF

                           IMID_IPHA = MIDACV( CV_NODI_IPHA )

                           ACV( IMID_IPHA ) =  ACV( IMID_IPHA ) &
                                +  SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW * ( 1. - INCOME ) * LIMD  & ! advection
                                +  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &  ! Diffusion contribution
                                +  SCVDETWEI( GI ) * ROBIN1  ! Robin bc
                           IF(GET_GTHETA) THEN
                              THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                   -  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX * T( CV_NODI_IPHA ) & !Diffusion contribution
                                   -  SCVDETWEI( GI ) * ROBIN1 * T( CV_NODI_IPHA )  ! Robin bc
                           ENDIF

                           ! CV_BETA=0 for Non-conservative discretisation (CV_BETA=1 for conservative disc)
                           ACV( IMID_IPHA ) =  ACV( IMID_IPHA )  &
                                - SECOND_THETA * FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW * LIMD 
                        ENDIF

                        TMID   =     T( CV_NODI_IPHA )
                        TOLDMID = TOLD( CV_NODI_IPHA )

                        ! Make allowances for no matrix stencil operating from outside the boundary.
                        BCZERO=1.0
                        IF( (SELE /= 0) .AND. (INCOME > 0.5) ) BCZERO=0.0

                        ! Put results into the RHS vector
                        CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA )  &
                                ! subtract 1st order adv. soln.
                             + SECOND_THETA * FTHETA_T2 * NDOTQNEW * SCVDETWEI( GI ) * LIMD * FVT * BCZERO &
                             -  SCVDETWEI( GI ) * ( FTHETA_T2 * NDOTQNEW * LIMDT &
                             + ONE_M_FTHETA_T2OLD * NDOTQOLD * LIMDTOLD ) ! hi order adv

                        !  Subtract out 1st order term non-conservative adv.
                        CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA ) &
                             - FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW * LIMD * TMID

                        ! High-order non-conservative advection contribution 
                        CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA ) &
                             + ( 1. - CV_BETA) * SCVDETWEI( GI ) &
                             * ( FTHETA_T2 * NDOTQNEW * TMID * LIMD  &
                             + ONE_M_FTHETA_T2OLD * NDOTQOLD * LIMDOLD * TOLDMID )

                        ! Diffusion contribution
                        CV_RHS( CV_NODI_IPHA ) =  CV_RHS( CV_NODI_IPHA ) &
                             + (1.-FTHETA) * SCVDETWEI(GI) * DIFF_COEFOLD_DIVDX  &
                             * (TOLD(CV_NODJ_IPHA) - TOLD(CV_NODI_IPHA)) &
                                ! Robin bc
                             + SCVDETWEI(GI) * ROBIN2
                        IF(GET_GTHETA) THEN
                           THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                + (1.-FTHETA) * SCVDETWEI(GI) * DIFF_COEFOLD_DIVDX  &
                                * (TOLD(CV_NODJ_IPHA) - TOLD(CV_NODI_IPHA)) &
                                ! Robin bc
                                + SCVDETWEI(GI) * ROBIN2
                        ENDIF

                        ! this is for the internal energy equation source term..
                        ! - p \div u
                        IF( THERMAL ) THEN
                           THERM_FTHETA = 1.

                           if( igot_t2 /= 0 ) then
                              CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                                   - CV_P( CV_NODI ) * SCVDETWEI( GI ) * ( &
                                   THERM_FTHETA * NDOTQNEW * LIMT2 & 
                                   + ( 1. - THERM_FTHETA ) * NDOTQOLD * LIMT2OLD )
                           else
                              CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                                   - CV_P( CV_NODI ) * SCVDETWEI( GI ) * ( &
                                   THERM_FTHETA * NDOTQNEW & 
                                   + ( 1. - THERM_FTHETA ) * NDOTQOLD )
                           end if

                        END IF

                     ENDIF Conditional_GETCV_DISC

                  END DO Loop_IPHASE

               END IF Conditional_integration

            END DO Loop_GCOUNT

         END DO Loop_CV_ILOC

      END DO Loop_Elements



      if( .false. .and. getct) then
         ewrite(3,*) 'after put_in_ct_rhs'
         ewrite(3,*) 'ct_rhs:', ct_rhs
      end if
      IF( .false. .and. GETCV_DISC ) THEN
         ewrite(3,*) 'after put_in_ct_rhs'
         ewrite(3,*) 'cv_rhs:', cv_rhs
      end if

      IF(GET_GTHETA) THEN
         DO CV_NODI = 1, CV_NONODS
            DO IPHASE = 1, NPHASE    
               CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
               THETA_GDIFF(CV_NODI_IPHA) = THETA_GDIFF(CV_NODI_IPHA) / MASS_CV(CV_NODI)
            END DO
         END DO
      ENDIF

      Conditional_GETCV_DISC2: IF( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion equations

         !ewrite(3,*)'before adding extra bits*****DEN:',DEN
         !ewrite(3,*)'before adding extra bits*****DENOLD:',DENOLD
         !ewrite(3,*)'before adding extra bits*****TOLD:',TOLD
         !ewrite(3,*)'before adding extra bits*****MEAN_PORE_CV:',MEAN_PORE_CV

         Loop_CVNODI2: DO CV_NODI = 1, CV_NONODS ! Put onto the diagonal of the matrix 

            Loop_IPHASE2: DO IPHASE = 1, NPHASE    
               CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
               ! For the gravity term
               !SOURCT2( CV_NODI_IPHA ) = SOURCT( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA ) 
               !SOURCT2( CV_NODI_IPHA ) = SOURCT2( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA ) 
               IMID_IPHA = MIDACV( CV_NODI_IPHA )

               IF(IGOT_T2==1) THEN
                  CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                       + MASS_CV(CV_NODI) * SOURCT(CV_NODI_IPHA) !&
                  !- CV_BETA * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) & ! conservative time term. 
                  !* TOLD(CV_NODI_IPHA)  &
                  !* (DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA) - DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA)) / DT

                  ACV( IMID_IPHA ) =  ACV( IMID_IPHA ) &
                       !+ (CV_BETA * DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA) &
                       + (CV_BETA * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA))  &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) / DT

                  CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                       + (CV_BETA * DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA))  &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) * TOLD(CV_NODI_IPHA) / DT 
               ELSE

                  CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                       + MASS_CV(CV_NODI) * SOURCT(CV_NODI_IPHA) !&
                  !- CV_BETA * MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) & ! conservative time term. 
                  !* TOLD(CV_NODI_IPHA) &
                  !* (DEN(CV_NODI_IPHA) - DENOLD(CV_NODI_IPHA)) / DT

                  ACV( IMID_IPHA ) =  ACV( IMID_IPHA ) &
                       !+ (CV_BETA * DENOLD(CV_NODI_IPHA) + (1.-CV_BETA) * DEN(CV_NODI_IPHA))  &
                       + (CV_BETA * DEN(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA))  &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) / DT

                  CV_RHS( CV_NODI_IPHA ) = CV_RHS( CV_NODI_IPHA ) &
                       + (CV_BETA * DENOLD(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA)) &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) * TOLD(CV_NODI_IPHA) / DT
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

      IF(GETCT) THEN ! Form the rhs of the discretised equations
         DIAG_SCALE_PRES = 0.0 ! Obtain the CV discretised CT eqations plus RHS

         DO IPHASE = 1, NPHASE
            DO CV_NODI = 1, CV_NONODS
               CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS

               !ewrite(3,*) 'MEAN_PORE_CV( CV_NODI ) , TOLD( CV_NODI_IPHA ) , DERIV( CV_NODI_IPHA ), DT , DENOLD( CV_NODI_IPHA ), CV_P( CV_NODI ) ', &
               !     MEAN_PORE_CV( CV_NODI ) , TOLD( CV_NODI_IPHA ) , DERIV( CV_NODI_IPHA ), DT , DENOLD( CV_NODI_IPHA ), CV_P( CV_NODI ) 

               if(.false.) then
                  CT_RHS( CV_NODI ) = CT_RHS( CV_NODI )  -  MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) &
                       * TOLD( CV_NODI_IPHA ) * ( DEN( CV_NODI_IPHA ) - DENOLD( CV_NODI_IPHA )  &
                       - DERIV( CV_NODI_IPHA ) * CV_P( CV_NODI ) ) / ( DT * DENOLD( CV_NODI_IPHA ) )   

                  DIAG_SCALE_PRES( CV_NODI ) = DIAG_SCALE_PRES( CV_NODI ) + &
                       MEAN_PORE_CV( CV_NODI ) * TOLD( CV_NODI_IPHA ) * DERIV( CV_NODI_IPHA )  &
                       / ( DT * DENOLD( CV_NODI_IPHA ) )
               end if

               if(.false.) then

                  CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) * ( &
                       (1.0-W_SUM_ONE) * MEAN_PORE_CV( CV_NODI ) * T( CV_NODI_IPHA ) / DT - &
                       MEAN_PORE_CV( CV_NODI ) *  TOLD( CV_NODI_IPHA ) * DENOLD( CV_NODI_IPHA ) / ( DT * DEN( CV_NODI_IPHA ) ) - &
                       MEAN_PORE_CV( CV_NODI ) * DERIV( CV_NODI_IPHA ) * CV_P( CV_NODI ) * T( CV_NODI_IPHA ) / ( DT * DEN( CV_NODI_IPHA ) ) &
                       )
                  IF(IPHASE==1) THEN ! Add constraint to force sum of volume fracts to be unity... 
                     ! W_SUM_ONE==1 applies the constraint
                     ! W_SUM_ONE==0 does NOT apply the constraint
                     CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) * &
                          W_SUM_ONE * MEAN_PORE_CV( CV_NODI )  / DT
                  ENDIF

               else

                  W_SUM_ONE1 = 1.0 ! =1 Applies constraint to T (if ==1.0)
                  W_SUM_ONE2 = 0.0 ! =1 Applies constraint to Told (if ==1.0)
                  ! the original working code used W_SUM_ONE1 = 1, W_SUM_ONE2 = 1
                  CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) * MEAN_PORE_CV( CV_NODI ) *( &
                       (1.0-W_SUM_ONE1) *  T( CV_NODI_IPHA ) / DT &
                       -(1.0-W_SUM_ONE2) *  TOLD( CV_NODI_IPHA ) / DT  &
                       +TOLD( CV_NODI_IPHA ) * (DEN( CV_NODI_IPHA )-DENOLD( CV_NODI_IPHA )) / ( DT * DEN( CV_NODI_IPHA ) ) - &
                       DERIV( CV_NODI_IPHA ) * CV_P( CV_NODI ) * T( CV_NODI_IPHA ) / ( DT * DEN( CV_NODI_IPHA ) ) &
                       )
                  IF(IPHASE==1) THEN ! Add constraint to force sum of volume fracts to be unity... 
                     ! W_SUM_ONE==1 applies the constraint
                     ! W_SUM_ONE==0 does NOT apply the constraint
                     CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) *   &
                          (W_SUM_ONE1 - W_SUM_ONE2) * MEAN_PORE_CV( CV_NODI )  / DT
                  ENDIF

               endif

               DIAG_SCALE_PRES( CV_NODI ) = DIAG_SCALE_PRES( CV_NODI ) + &
                    MEAN_PORE_CV( CV_NODI ) * T( CV_NODI_IPHA ) * DERIV( CV_NODI_IPHA )  &
                    / ( DT * DEN( CV_NODI_IPHA ) )

               !ewrite(3,*) 'MASS_CV( CV_NODI ), SOURCT( CV_NODI_IPHA ), ABSORBT( CV_NODI, :, : ), T( CV_NODI_IPHA ), DEN( CV_NODI_IPHA ) ', &
               !     MASS_CV( CV_NODI ), SOURCT( CV_NODI_IPHA ), ABSORBT( CV_NODI, :, : ), T( CV_NODI_IPHA ), DEN( CV_NODI_IPHA )

               CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) + MASS_CV( CV_NODI ) * SOURCT( CV_NODI_IPHA ) / DEN( CV_NODI_IPHA )
               !CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) + MASS_CV( CV_NODI ) * SOURCT2( CV_NODI_IPHA ) / DEN( CV_NODI_IPHA )

               DO JPHASE = 1, NPHASE   
                  CV_NODI_JPHA = CV_NODI + ( JPHASE - 1 ) * CV_NONODS
                  CT_RHS( CV_NODI ) = CT_RHS( CV_NODI )  &
                       - MASS_CV( CV_NODI ) * ABSORBT( CV_NODI, IPHASE, JPHASE ) * T( CV_NODI_JPHA ) / DEN( CV_NODI_IPHA )
               END DO
            END DO
         END DO
      END IF

      !ewrite(3,*)'upwind fraction:'
      !ewrite(3,*) 'This is wrong now, but up_wind_nod is not used anyway I think', cv_nonods, x_nonods
      !do iphase=1,nphase
      !   ewrite(3,*)'for phase iphase=',iphase
      !   do ele=1,totele-1
      !      do cv_iloc=1,cv_nloc
      !         cv_nodi = cv_ndgln((ele-1)*cv_nloc+cv_iloc)
      !         mat_nodi = mat_ndgln((ele-1)*cv_nloc+cv_iloc)
      !         cv_nodi_IPHA=cv_nodi +(IPHASE-1)*CV_NONODS
      !
      !         if(cv_nonods==x_nonods) then
      !            !ewrite(3,*)0.5*(x(cv_nodi)+x(cv_nodi+1)),UP_WIND_NOD(cv_nodi_IPHA)
      !            ewrite(3,*)0.5*(x(x_ndgln((ele-1)*x_nloc+cv_iloc))+x(x_ndgln((ele-1)*x_nloc+cv_iloc+1))),  &
      !                 UP_WIND_NOD(cv_nodi_IPHA)
      !         else
      !            if(cv_iloc==cv_nloc) then
      !               ewrite(3,*)x(x_ndgln((ele-1)*x_nloc+cv_iloc)),  &
      !                    UP_WIND_NOD(cv_nodi_IPHA)
      !            else
      !               ewrite(3,*)0.5*(x(x_ndgln((ele-1)*x_nloc+cv_iloc))+x(x_ndgln((ele-1)*x_nloc+cv_iloc+1))),  &
      !                    UP_WIND_NOD(cv_nodi_IPHA)
      !            endif
      !         endif
      !
      !      end do
      !   end do
      !end do

      ! for the output
      T_FEMT = FEMT
      DEN_FEMT = FEMDEN

      ewrite(3,*) '----------sub cv_assemb--------'
      if( .false. .and. getct) then
         ewrite(3,*) 'ct_rhs:', ct_rhs
      end if
      if( .false. .and. GETCV_DISC ) then
         ewrite(3,*) 'cv_rhs:', cv_rhs
      end if

      !if( .false. .and. getct) then
      !   ewrite(3,*)'ct(1:ncolct);',ct(1:ncolct)
      !   ewrite(3,*)'ct(1+ncolct:2*ncolct);',ct(1+ncolct:2*ncolct)
      !end if

      if (.false. .and. getcv_disc) then

         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'
         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'

         ewrite(3,*)'cv_rhs_1:',cv_rhs(1:cv_nonods)
         ewrite(3,*)'cv_rhs_2:',cv_rhs(cv_nonods+1:2*cv_nonods)

         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'
         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'

         ewrite(3,*) '11::',0.5 * mass_cv* t(1:cv_nonods) / dt
         ewrite(3,*) '22::',0.5 * mass_cv * t(cv_nonods+1:2*cv_nonods) / dt

         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'
         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'

         ewrite(3,*) '11::',0.5 * mass_cv* t(1:cv_nonods) / dt - cv_rhs(1:cv_nonods)
         ewrite(3,*) '22::',0.5 * mass_cv * t(cv_nonods+1:2*cv_nonods) / dt - cv_rhs(cv_nonods+1:2*cv_nonods)

         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'
         ewrite(3,*) '-----------------------------------------------------------------------------------------------------------------------------------'

         ewrite(3,*) 'dT, igot_t2, thermal, getmat, cv_beta:', dT, igot_t2, thermal, getmat, cv_beta
         ewrite(3,*) 'MASS:', MASS_CV
         ewrite(3,*) 'T:', T
         ewrite(3,*) 'TOLD:', TOLD
      end if

      ewrite(3,*)'IN cv_adv_dif a CV representation t:'
      CALL PRINT_CV_DIST(CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE, &
           T, X_NDGLN, CV_NDGLN, X) 
      ewrite(3,*) 'just print out - in cv_assemb'

      ! Deallocating temporary working arrays
      DEALLOCATE( JCOUNT_KLOC )
      DEALLOCATE( CVNORMX )
      DEALLOCATE( CVNORMY )
      DEALLOCATE( CVNORMZ )
      DEALLOCATE( COLGPTS ) ! The size of this vector is over-estimated
      DEALLOCATE( FINDGPTS )
      DEALLOCATE( SNDOTQ )
      DEALLOCATE( SNDOTQOLD )
      DEALLOCATE( CV_ON_FACE )
      DEALLOCATE( U_ON_FACE )
      DEALLOCATE( CVFEM_ON_FACE )
      DEALLOCATE( UFEM_ON_FACE )
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
      DEALLOCATE( SCVFENX )
      DEALLOCATE( SCVFENY )
      DEALLOCATE( SCVFENZ )
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

      DEALLOCATE( TMIN_2ND_MC )
      DEALLOCATE( TMAX_2ND_MC )
      DEALLOCATE( TOLDMIN_2ND_MC )
      DEALLOCATE( TOLDMAX_2ND_MC )
      DEALLOCATE( DENMIN_2ND_MC )
      DEALLOCATE( DENMAX_2ND_MC )
      DEALLOCATE( DENOLDMIN_2ND_MC )
      DEALLOCATE( DENOLDMAX_2ND_MC )

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

      RETURN  

    END SUBROUTINE CV_ASSEMB




    SUBROUTINE FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC,  &
         MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI, &
         TOTELE, X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
         CV_SNLOC, CVFEM_ON_FACE, SCVNGI, GI, X_SHARE, X_NONODS, ELE, ELE2,  &
         FINELE, COLELE, NCOLELE) 
      ! We are on the boundary or next to another element. Determine CV_OTHER_LOC,
      ! U_OTHER_LOC. 
      ! CVFEM_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
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
      LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in ) :: CVFEM_ON_FACE
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
         X_SHARE( X_NODK ) = CVFEM_ON_FACE( X_KLOC, GI )
      END DO

      ELE3 = 0
      DO COUNT = FINELE( ELE ), FINELE( ELE + 1 ) - 1, 1
         ELE2 = COLELE( COUNT )
         SUF_COUNT = 0 ! See if we share the same nodes
         IF(ELE2.NE.ELE) THEN
            DO X_KLOC = 1, X_NLOC
               X_NODK = X_NDGLN(( ELE2 - 1 ) * X_NLOC + X_KLOC )
               IF( X_SHARE( X_NODK )) SUF_COUNT = SUF_COUNT + 1
            END DO
         ENDIF
         IF( SUF_COUNT == CV_SNLOC ) ELE3 = ELE2
         !ewrite(3,*)'suf_count:', ele, ele2, suf_count, cv_snloc
      END DO

      DO X_KLOC = 1, X_NLOC
         X_NODK = X_NDGLN(( ELE - 1 ) * X_NLOC + X_KLOC )
         X_SHARE( X_NODK ) = .FALSE.
      END DO

      ELE2 = ELE3
      IF( ELE2 /= 0 ) THEN 
         ! Is CV_NODI in element ELE2 if yes set ELE2=0 as we dont want to integrate in 
         ! the middle of a CV. 
         DO CV_KLOC = 1, CV_NLOC
            CV_NODK = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_KLOC )
            !ewrite(3,*)'cv_nodi, cv_nodk:', cv_nodi, cv_nodk, ele, ele2
            IF( CV_NODK == CV_NODI ) INTEGRAT_AT_GI = .FALSE.
         END DO
      ENDIF

      IF( ( ELE2 /= 0) .AND. INTEGRAT_AT_GI ) THEN ! Determine CV_OTHER_LOC(CV_KLOC)
         CV_OTHER_LOC = 0
         DO CV_KLOC = 1, CV_NLOC
            IF( CVFEM_ON_FACE( CV_KLOC, GI )) THEN ! Find opposite local node
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
                  IF( ( XU_NODK2 == XU_NODK ) .OR. ( XU_NLOC == 1 ) ) THEN
                     DO ILEV = 1, CV_NLOC
                        JLEV = CV_OTHER_LOC( ILEV )
                        IF( JLEV /= 0 ) THEN 
                           U_OTHER_LOC( U_KLOC + ( ILEV - 1 ) * XU_NLOC ) = U_KLOC2 + &
                                ( JLEV - 1 ) * XU_NLOC
                        ENDIF
                     END DO
                  END IF
               END DO
            END DO
         ELSE
            ! Works for non constant and constant (XU_NLOC=1) vel basis functions...
            DO U_KLOC = 1, U_NLOC ! Find opposite local node
               XU_NODK = XU_NDGLN(( ELE - 1 ) * XU_NLOC + U_KLOC )
               DO U_KLOC2 = 1, U_NLOC
                  XU_NODK2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_KLOC2 )
                  IF( ( XU_NODK2 == XU_NODK ) .OR. ( XU_NLOC == 1 ) ) &
                       U_OTHER_LOC( U_KLOC ) = U_KLOC2
               END DO
            END DO
         ENDIF

         MAT_OTHER_LOC = CV_OTHER_LOC
      ELSE
         CV_OTHER_LOC = 0
         U_OTHER_LOC = 0
         MAT_OTHER_LOC = 0
      ENDIF

      RETURN

    END SUBROUTINE FIND_OTHER_SIDE




    SUBROUTINE PROJ_CV_TO_FEM_4( state, &
         FEMT, FEMTOLD, FEMDEN, FEMDENOLD, T, TOLD, DEN, DENOLD, &
         IGOT_T2,T2,T2OLD, FEMT2,FEMT2OLD, &
         XC_CV,YC_CV,ZC_CV, MASS_CV, MASS_ELE, &
         NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI, CV_NLOC, CVN, CVWEIGHT, N, NLX, NLY, NLZ, &
         X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
         IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC )

      ! determine FEMT (finite element wise) etc from T (control volume wise) 
      IMPLICIT NONE
      type( state_type ), dimension( : ), intent( in ) :: state
      INTEGER, intent( in ) :: NDIM, NPHASE, CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
           X_NONODS, NCOLM, IGOT_T2, IGETCT, NCOLCMC
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

      REAL, DIMENSION( IGETCT*NCOLCMC ), intent( inout ) :: MASS_MN_PRES
      INTEGER, DIMENSION( IGETCT*(CV_NONODS + 1) ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( IGETCT*NCOLCMC ), intent( in ) :: COLCMC
      ! Local variables
      REAL, DIMENSION( : ), allocatable :: PSI, FEMPSI, PSI_AVE, PSI_INT
      INTEGER :: NTSOL,NTSOL_AVE,NTSOL_INT,ELE,CV_ILOC,X_INOD,CV_INOD,NL,NFIELD
      CHARACTER(100) :: PATH
      INTEGER :: velocity_max_iterations, nstates, istate
      LOGICAL :: solve_force_balance, have_component

      
     ewrite(3,*) 'In PROJ_CV_TO_FEM_4'

      nstates = option_count( "/material_phase" )
      have_component = .false.
      do istate = 1, nstates
         if( have_option( '/material_phase[' // int2str( istate - 1 ) //']/is_multiphase_component/' ) ) then
            have_component = .true.
         end if
      end do

      solve_force_balance = .false.
      call get_option( "/material_phase[0]/vector_field::Velocity/prognostic/solver/max_iterations", &
           velocity_max_iterations,  default =  500 )
      if( velocity_max_iterations /= 0 ) solve_force_balance = .true.

      if ( solve_force_balance .or. have_component ) then
         path = '/material_phase[0]/scalar_field::Pressure' 
      else
         path = '/material_phase[0]/scalar_field::Temperature' 
      end if

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
           IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC, PATH )

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
         IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC, PATH )

      ! Determine FEMT (finite element wise) etc from T (control volume wise)
      ! Also integrate PSI_INT over each CV and average PSI_AVE over each CV. 
      use shape_functions 
      use shape_functions_Linear_Quadratic
      use solvers_module
      use matrix_operations
      IMPLICIT NONE
      INTEGER, intent( in ) :: NTSOL, NTSOL_AVE, NTSOL_INT, NDIM, CV_NONODS, TOTELE, &
           X_NLOC, CV_NGI, CV_NLOC, X_NONODS, NCOLM, IGETCT, NCOLCMC
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
      CHARACTER(100), intent(in) :: PATH

      REAL, DIMENSION( IGETCT*NCOLCMC ), intent( inout ) :: MASS_MN_PRES
      INTEGER, DIMENSION( IGETCT*(CV_NONODS + 1) ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( IGETCT*NCOLCMC ), intent( in ) :: COLCMC
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

      FEMPSI = 0.0

      PSI_AVE2 = PSI_AVE
      PSI_INT2 = PSI_INT
      PSI_AVE = 0.0
      PSI_INT = 0.0

      FEMPSI_RHS = 0.0
      MAT = 0.0
      MASS_CV = 0.0
      IF(IGETCT.NE.0) MASS_MN_PRES=0.0

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
            !ewrite(3,*)'ele,CV_NODI,CV_ILOC:',ele,CV_NODI,CV_ILOC, x(CV_NODI)

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

               IF(IGETCT.NE.0) THEN
                  CALL POSINMAT( COUNT, CV_NODI, CV_NODJ, CV_NONODS, FINDCMC, COLCMC, NCOLCMC )

                  MASS_MN_PRES( COUNT ) = MASS_MN_PRES( COUNT ) + MN
               ENDIF

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

      ! Form average...
      DO CV_NODI=1,CV_NONODS
         DO IT = 1, NTSOL_AVE
            PSI_AVE( CV_NODI + ( IT - 1 ) * CV_NONODS ) = PSI_AVE( CV_NODI +( IT - 1 ) &
                 * CV_NONODS ) / MASS_CV( CV_NODI )
         END DO
      END DO

      ! Solve...
      DO IT = 1, NTSOL
         CALL SOLVER( MAT,  &
              FEMPSI( 1 + (IT - 1 ) * CV_NONODS : CV_NONODS + (IT - 1 ) * CV_NONODS ),  &
              FEMPSI_RHS( 1 + ( IT - 1 ) * CV_NONODS : CV_NONODS + (IT - 1 ) * CV_NONODS ),  &
              FINDM, &
              COLM, &
              option_path = trim(path) )
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
         DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
         DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
         DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
         NDIM, NDIM_VEL, NPHASE, U_NONODS, TOTELE, U_NDGLN, &
         XU_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI, U_NLOC, CVWEIGHT, &
         N, NLX, NLY, NLZ, &
         CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
         X_NONODS, X, Y, Z, &
         NFACE, FACE_ELE, U_SLOCLIST, CV_SLOCLIST, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC,  &
         SUF_U_BC,SUF_V_BC,SUF_W_BC, &
         WIC_U_BC_DIRICHLET, SBCVNGI, SBUFEN, SBUFENSLX, SBUFENSLY, SBWEIGH, & 
         SBCVFEN, SBCVFENSLX, SBCVFENSLY)   

      ! determine FEMT (finite element wise) etc from T (control volume wise)
      use shape_functions 
      use matrix_operations
      IMPLICIT NONE
      INTEGER, intent( in ) :: NDIM, NDIM_VEL, NPHASE, &
           U_NONODS, TOTELE, X_NLOC, CV_NGI, U_NLOC, &
           X_NONODS, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC_DIRICHLET, SBCVNGI, NFACE
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, UOLD, V, VOLD, W, WOLD
      REAL, DIMENSION( U_NLOC, NPHASE, TOTELE ), intent( inout ) :: DUX_ELE,DUY_ELE,DUZ_ELE, &
           DUOLDX_ELE,DUOLDY_ELE,DUOLDZ_ELE
      REAL, DIMENSION( U_NLOC, NPHASE, TOTELE ), intent( inout ) :: DVX_ELE,DVY_ELE,DVZ_ELE, &
           DVOLDX_ELE,DVOLDY_ELE,DVOLDZ_ELE
      REAL, DIMENSION( U_NLOC, NPHASE, TOTELE ), intent( inout ) :: DWX_ELE,DWY_ELE,DWZ_ELE, &
           DWOLDX_ELE,DWOLDY_ELE,DWOLDZ_ELE
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( STOTEL*NPHASE ), intent( in ) ::  WIC_U_BC
      REAL, DIMENSION( STOTEL*U_SNLOC*NPHASE ), intent( in ) ::  SUF_U_BC,SUF_V_BC,SUF_W_BC
      INTEGER, DIMENSION( NFACE,U_SNLOC ), intent( in ) ::  U_SLOCLIST
      INTEGER, DIMENSION( NFACE,CV_SNLOC ), intent( in ) ::  CV_SLOCLIST
      INTEGER, DIMENSION( NFACE,TOTELE ), intent( in ) ::  FACE_ELE
      REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGHT
      REAL, DIMENSION( U_NLOC, CV_NGI ), intent( in ) :: N, NLX, NLY, NLZ 
      REAL, DIMENSION( X_NLOC, CV_NGI ), intent( in ) :: CVFEN, CVFENLX, CVFENLY, CVFENLZ
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( U_SNLOC, SBCVNGI ), intent( in ) :: SBUFEN, SBUFENSLX, SBUFENSLY
      REAL, DIMENSION( CV_SNLOC, SBCVNGI ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
      REAL, DIMENSION( SBCVNGI ), intent( in ) :: SBWEIGH

      CALL DG_DERIVS( U, UOLD, &
           DUX_ELE, DUY_ELE, DUZ_ELE, DUOLDX_ELE, DUOLDY_ELE, DUOLDZ_ELE, &
           NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, &
           XU_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI, U_NLOC, CVWEIGHT, &
           N, NLX, NLY, NLZ, &
           CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           X_NONODS, X, Y, Z, &
           NFACE, FACE_ELE, U_SLOCLIST, CV_SLOCLIST, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC, SUF_U_BC, &
           WIC_U_BC_DIRICHLET, SBCVNGI, SBUFEN, SBUFENSLX, SBUFENSLY, SBWEIGH, & 
           SBCVFEN, SBCVFENSLX, SBCVFENSLY)   

      IF(NDIM_VEL.GE.2) THEN
         CALL DG_DERIVS( V, VOLD, &
              DVX_ELE, DVY_ELE, DVZ_ELE, DVOLDX_ELE, DVOLDY_ELE, DVOLDZ_ELE, &
              NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, &
              XU_NDGLN, X_NLOC, X_NDGLN, &
              CV_NGI, U_NLOC, CVWEIGHT, &
              N, NLX, NLY, NLZ, &
              CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
              X_NONODS, X, Y, Z, &
              NFACE, FACE_ELE, U_SLOCLIST, CV_SLOCLIST, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC, SUF_V_BC, &
              WIC_U_BC_DIRICHLET, SBCVNGI, SBUFEN, SBUFENSLX, SBUFENSLY, SBWEIGH, & 
              SBCVFEN, SBCVFENSLX, SBCVFENSLY)  
      ELSE
         DVX_ELE=0; DVY_ELE=0; DVZ_ELE=0; DVOLDX_ELE=0; DVOLDY_ELE=0; DVOLDZ_ELE=0
      ENDIF

      IF(NDIM_VEL.GE.3) THEN
         CALL DG_DERIVS( W, WOLD, &
              DWX_ELE, DWY_ELE, DWZ_ELE, DWOLDX_ELE, DWOLDY_ELE, DWOLDZ_ELE, &
              NDIM, NPHASE, U_NONODS, TOTELE, U_NDGLN, &
              XU_NDGLN, X_NLOC, X_NDGLN, &
              CV_NGI, U_NLOC, CVWEIGHT, &
              N, NLX, NLY, NLZ, &
              CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
              X_NONODS, X, Y, Z, &
              NFACE, FACE_ELE, U_SLOCLIST, CV_SLOCLIST, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC, SUF_W_BC, &
              WIC_U_BC_DIRICHLET, SBCVNGI, SBUFEN, SBUFENSLX, SBUFENSLY, SBWEIGH, & 
              SBCVFEN, SBCVFENSLX, SBCVFENSLY) 
      ELSE
         DWX_ELE=0; DWY_ELE=0; DWZ_ELE=0; DWOLDX_ELE=0; DWOLDY_ELE=0; DWOLDZ_ELE=0
      ENDIF

      RETURN
    END SUBROUTINE DG_DERIVS_UVW




    SUBROUTINE DG_DERIVS( FEMT, FEMTOLD, &
         DTX_ELE, DTY_ELE, DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
         NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
         XCV_NDGLN, X_NLOC, X_NDGLN,&
         CV_NGI, CV_NLOC, CVWEIGHT, &
         N, NLX, NLY, NLZ, &
         X_N, X_NLX, X_NLY, X_NLZ, &
         X_NONODS, X, Y, Z, &
         NFACE, FACE_ELE, CV_SLOCLIST, X_SLOCLIST, STOTEL, CV_SNLOC, X_SNLOC, WIC_T_BC, SUF_T_BC, &
         WIC_T_BC_DIRICHLET, SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH, &
         X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY ) 

      ! determine FEMT (finite element wise) etc from T (control volume wise)
      use shape_functions 
      use shape_functions_NDim
      use matrix_operations 
      IMPLICIT NONE
      INTEGER, intent( in ) :: NDIM, NPHASE, CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
           X_NONODS, STOTEL, CV_SNLOC, X_SNLOC, WIC_T_BC_DIRICHLET, SBCVNGI, NFACE
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: FEMT, FEMTOLD
      REAL, DIMENSION( CV_NLOC, NPHASE, TOTELE ), intent( inout ) :: DTX_ELE, DTY_ELE, DTZ_ELE, &
           DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) ::  XCV_NDGLN
      INTEGER, DIMENSION( STOTEL*NPHASE ), intent( in ) ::  WIC_T_BC
      REAL, DIMENSION( STOTEL*CV_SNLOC*NPHASE ), intent( in ) ::  SUF_T_BC
      INTEGER, DIMENSION( NFACE,CV_SNLOC ), intent( in ) ::  CV_SLOCLIST
      INTEGER, DIMENSION( NFACE,X_SNLOC ), intent( in ) ::  X_SLOCLIST
      INTEGER, DIMENSION( NFACE,TOTELE ), intent( in ) ::  FACE_ELE
      REAL, DIMENSION( CV_NGI ), intent( inout ) :: CVWEIGHT
      REAL, DIMENSION( CV_NLOC, CV_NGI ), intent( in ) :: N, NLX, NLY, NLZ 
      REAL, DIMENSION( X_NLOC, CV_NGI ), intent( in ) :: X_N, X_NLX, X_NLY, X_NLZ 
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( CV_SNLOC, SBCVNGI ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
      REAL, DIMENSION( X_SNLOC, SBCVNGI ), intent( in ) :: X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY
      REAL, DIMENSION( SBCVNGI ), intent( in ) :: SBWEIGH
      ! Local variables
      LOGICAL :: D1, D3, DCYL, APPLYBC
      REAL, DIMENSION( : ), allocatable :: DETWEI, RA
      REAL, DIMENSION( :, : ), allocatable :: NX, NY, NZ, X_NX, X_NY, X_NZ
      REAL, DIMENSION( :, :, : ), allocatable :: MASELE
      REAL, DIMENSION( :, :, : ), allocatable :: VTX_ELE, VTY_ELE, VTZ_ELE, VTOLDX_ELE, VTOLDY_ELE, VTOLDZ_ELE
      REAL, DIMENSION( :, : ), allocatable :: MASS, INV_MASS
      REAL, DIMENSION( : ), allocatable :: VTX, VTY, VTZ, VTOLDX, VTOLDY, VTOLDZ, DTX, DTY, DTZ, DTOLDX, DTOLDY, DTOLDZ
      REAL, DIMENSION( : ), allocatable :: XSL, YSL, ZSL, SNORMXN, SNORMYN, SNORMZN, SDETWE
      INTEGER, DIMENSION( : ), allocatable :: SLOC2LOC, X_SLOC2LOC, ILOC_OTHER_SIDE
      REAL :: VOLUME, NN, NNX, NNY, NNZ, NORMX, NORMY, NORMZ, SAREA, NRBC, RNN, RTBC, &
           VLM_NORX, VLM_NORY, VLM_NORZ
      INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_GI, CV_ILOC2, &
           CV_INOD, CV_INOD2, CV_JLOC2, CV_NODJ2, CV_NODJ2_IPHA, CV_NODJ_IPHA, &
           CV_SILOC, CV_SJLOC, ELE2, IFACE, IPHASE, SELE2, SUF_CV_SJ2, SUF_CV_SJ2_IPHA, &
           X_INOD, SGI, X_SILOC, X_ILOC

      ewrite(3,*)'in DG_DERIVS sbrt'
      ALLOCATE( DETWEI( CV_NGI )) 
      ALLOCATE( RA( CV_NGI ))
      ALLOCATE( NX( CV_NLOC, CV_NGI ))
      ALLOCATE( NY( CV_NLOC, CV_NGI ))
      ALLOCATE( NZ( CV_NLOC, CV_NGI ))
      ALLOCATE( X_NX( X_NLOC, CV_NGI ))
      ALLOCATE( X_NY( X_NLOC, CV_NGI ))
      ALLOCATE( X_NZ( X_NLOC, CV_NGI ))
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
      ALLOCATE( X_SLOC2LOC( X_SNLOC ))
      ALLOCATE( ILOC_OTHER_SIDE( CV_SNLOC ))
      ALLOCATE( XSL( X_SNLOC ))
      ALLOCATE( YSL( X_SNLOC ))
      ALLOCATE( ZSL( X_SNLOC ))
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
     ! ewrite(3,*)'****X_NLX:',X_NLX
     ! ewrite(3,*)'****X_NLY:',X_NLY

      Loop_Elements1: DO ELE = 1, TOTELE

         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         !CALL DETNLXR_SUPER( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, X_NLOC, CV_NLOC, CV_NGI, &
         !     X_N, X_NLX, X_NLY, X_NLZ, N, NLX, NLY, NLZ, &
         !     CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
         !     NX, NY, NZ ) 
         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
              X_NLOC, X_NLOC, CV_NGI, &
              X_N, X_NLX, X_NLY, X_NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
              X_NX, X_NY, X_NZ, &
              CV_NLOC, NLX, NLY, NLZ, NX, NY, NZ ) 

         !ewrite(3,*)'N',N
         !ewrite(3,*)'nlx:',nlx
         !ewrite(3,*)'nx:',nx
         !ewrite(3,*)'CVWEIGHT:',CVWEIGHT
         !ewrite(3,*)'DETWEI:',DETWEI
         !ewrite(3,*)'volume=',volume
         !stop 12


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

     ! ewrite(3,*)'totele=',totele

      Loop_Elements2: DO ELE=1,TOTELE

         Between_Elements_And_Boundary: DO IFACE = 1, NFACE
            ELE2 = FACE_ELE( IFACE, ELE )
            SELE2 = MAX( 0, - ELE2 )
            ELE2 = MAX( 0, + ELE2 )
            !ewrite(3,*)'FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE ):',FACE_ELE( 1, ELE ),FACE_ELE( 2, ELE )

            ! The surface nodes on element face IFACE.  
            SLOC2LOC( : ) = CV_SLOCLIST( IFACE, : )
            X_SLOC2LOC( : ) = X_SLOCLIST( IFACE, : )

            ! Form approximate surface normal (NORMX,NORMY,NORMZ)
            CALL DGSIMPLNORM( ELE, X_SLOC2LOC, TOTELE, X_NLOC, X_SNLOC, X_NDGLN, &
                 X, Y, Z, X_NONODS, NORMX, NORMY, NORMZ )

            ! Recalculate the normal...
            DO X_SILOC = 1, X_SNLOC
               X_ILOC = X_SLOC2LOC( X_SILOC )
               X_INOD = X_NDGLN(( ELE - 1 ) * X_NLOC + X_ILOC )
               XSL( X_SILOC ) = X( X_INOD )
               YSL( X_SILOC ) = Y( X_INOD )
               ZSL( X_SILOC ) = Z( X_INOD )
            END DO

            CALL DGSDETNXLOC2(X_SNLOC, SBCVNGI, &
                 XSL, YSL, ZSL, &
                 X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY, SBWEIGH, SDETWE, SAREA, &
                 (NDIM==1), (NDIM==3), (NDIM==-2), &
                 SNORMXN, SNORMYN, SNORMZN, &
                 NORMX, NORMY, NORMZ )

            !ewrite(3,*)'*********************'
            !ewrite(3,*)'ele,ele2,sele2:',ele,ele2,sele2
            !ewrite(3,*)'iface=',iface
            IF(SELE2 == 0) THEN
               ! Calculate the nodes on the other side of the face:
               !ewrite(3,*)'X_NLOC,CV_SNLOC,CV_NLOC,ele,ele2:',X_NLOC,CV_SNLOC,CV_NLOC,ele,ele2
               !ewrite(3,*)'SLOC2LOC:',SLOC2LOC
               DO CV_SILOC = 1, CV_SNLOC
                  CV_ILOC = SLOC2LOC( CV_SILOC )
                  CV_INOD = XCV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
                  ! ewrite(3,*)'CV_SILOC,CV_ILOC,CV_INOD:',CV_SILOC,CV_ILOC,CV_INOD
                  DO CV_ILOC2 = 1, CV_NLOC
                     CV_INOD2 = XCV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )
                 ! ewrite(3,*)'CV_INOD2,CV_INOD=',CV_INOD2,CV_INOD
                     IF( CV_INOD2 == CV_INOD ) ILOC_OTHER_SIDE( CV_SILOC ) = CV_ILOC2
                  END DO
               END DO
               ! ewrite(3,*)'ILOC_OTHER_SIDE:',ILOC_OTHER_SIDE
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
                       ! ewrite(3,*)'(ELE2-1)*CV_NLOC+CV_JLOC2,ELE2,CV_NLOC,CV_JLOC2:', &
                       !      (ELE2-1)*CV_NLOC+CV_JLOC2,ELE2,CV_NLOC,CV_JLOC2
                       ! ewrite(3,*)'ILOC_OTHER_SIDE:',ILOC_OTHER_SIDE
                       ! ewrite(3,*)'ELE,ELE2,SELE2,CV_JLOC2=',ELE,ELE2,SELE2,CV_JLOC2
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
                     END DO
                  END DO
               END DO
            ENDIF


         END DO Between_Elements_And_Boundary


      END DO Loop_Elements2



      ! Solve local system for the gradients DTX_ELE etc:
!        ewrite(3,*)'masele:', masele   
!        ewrite(3,*)'ndim=',ndim 

      Loop_Elements3: DO ELE=1,TOTELE
        ! ewrite(3,*)'ele=',ele

         MASS(:,:)=MASELE(:,:,ELE)
        ! ewrite(3,*)'mass=',mass
        ! ewrite(3,*)'MASELE(:,:,ELE):',MASELE(:,:,ELE)
         CALL MATDMATINV( MASS, INV_MASS, CV_NLOC)
        ! ewrite(3,*)'here 1'
        ! ewrite(3,*)'inv_mass=',inv_mass

         !INV_MASS = 0.0
         !INV_MASS(1,1) = 1./sum( MASS(1,:))
         !INV_MASS(2,2) =  1./sum( MASS(2,:))
         !INV_MASS(3,3) =  1./sum( MASS(3,:))


         Loop_IPHASE: DO IPHASE=1,NPHASE

            VTX(:)=VTX_ELE(:, IPHASE,ELE )
            VTY(:)=VTY_ELE(:, IPHASE,ELE )
            VTZ(:)=VTZ_ELE(:, IPHASE,ELE )
            VTOLDX(:)=VTOLDX_ELE(:, IPHASE,ELE )
            VTOLDY(:)=VTOLDY_ELE(:, IPHASE,ELE )
            VTOLDZ(:)=VTOLDZ_ELE(:, IPHASE,ELE )
            ! ewrite(3,*)'heree 2'
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
            ! ewrite(3,*)'heree 3'
            DTX_ELE(:, IPHASE,ELE )=DTX(:)
            ! ewrite(3,*)'heree 3.01'
            DTY_ELE(:, IPHASE,ELE )=DTY(:)
            ! ewrite(3,*)'heree 3.02'
            DTZ_ELE(:, IPHASE,ELE )=DTZ(:)
            ! ewrite(3,*)'heree 3.1'
            DTOLDX_ELE(:, IPHASE,ELE )=DTOLDX(:)
            DTOLDY_ELE(:, IPHASE,ELE )=DTOLDY(:)
            DTOLDZ_ELE(:, IPHASE,ELE )=DTOLDZ(:)
            ! ewrite(3,*)'heree 3.2'

         END DO Loop_IPHASE

      END DO Loop_Elements3


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
      ! ewrite(3,*)'about to leave DG_DERIVS'

      RETURN

    END SUBROUTINE DG_DERIVS




    SUBROUTINE ONVDLIM_ALL( TOTELE, &
         TDLIM, TDCEN, INCOME, PELE, PELEOT, &
         ETDNEW, TDMIN, TDMAX, &
         TDMIN_2nd_mc, TDMAX_2nd_mc, FIRORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE, &
         IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, TUPWIND_MAT )
      implicit none 
      ! This sub calculates the limited face values TDADJ(1...SNGI) from the central 
      ! difference face values TDCEN(1...SNGI) using a NVD shceme.  
      ! The method is based on a new approach to limiting based on applying 
      ! limiting only when we have a local min or max.
      ! INCOME(1...SNGI)=1 for incomming to element ELE  else =0.
      ! LIBETA is the flux limiting parameter. 
      ! TDMAX(PELE)=maximum of the surrounding element values of element PELE.
      ! TDMIN(PELE)=minimum of the surrounding element values of element PELE.
      ! TDMIN_2nd_mc(PELE), TDMAX_2nd_mc(PELE) are the second minima and max 
      ! excluding (minus c) these values surrounding them. 
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
      REAL, intent( in ) :: TDCEN, INCOME, COURANT_OR_MINUS_ONE
      INTEGER, intent( in ) :: PELE, PELEOT
      REAL, DIMENSION( TOTELE ), intent( in ) :: ETDNEW, TDMIN, TDMAX, TDMIN_2nd_mc, TDMAX_2nd_mc
      LOGICAL, intent( in ) :: FIRORD, NOLIMI, LIMIT_USE_2ND

      INTEGER, intent( in ) :: IANISOTROPIC
      INTEGER, intent( in ) :: NSMALL_COLM
      INTEGER, DIMENSION( (TOTELE+1)*IANISOTROPIC ), intent( in ) :: SMALL_FINDRM
      INTEGER, DIMENSION( NSMALL_COLM*IANISOTROPIC ), intent( in ) :: SMALL_COLM
      REAL, DIMENSION( NSMALL_COLM*IANISOTROPIC), intent( in ) :: TUPWIND_MAT 

      IF(IANISOTROPIC==1) THEN ! limit based on 2 largest and 2 minima values of T:
         CALL ONVDLIM_ANO( TOTELE, &
         TDLIM, TDCEN, INCOME, PELE, PELEOT, &
         ETDNEW, TDMIN, TDMAX, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE, &
         SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, TUPWIND_MAT )
      ELSE IF(LIMIT_USE_2ND) THEN ! limit based on 2 largest and 2 minima values of T:
         CALL ONVDLIM_2nd( TOTELE, &
         TDLIM, TDCEN, INCOME, PELE, PELEOT, &
         ETDNEW, TDMIN, TDMAX, TDMIN_2nd_mc, TDMAX_2nd_mc, FIRORD, NOLIMI )
      ELSE
         CALL ONVDLIM( TOTELE, &
         TDLIM, TDCEN, INCOME, PELE, PELEOT, &
         ETDNEW, TDMIN, TDMAX, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE )
      ENDIF

      RETURN

      END SUBROUTINE ONVDLIM_ALL




    SUBROUTINE ONVDLIM_2nd( TOTELE, &
         TDLIM, TDCEN, INCOME, PELE, PELEOT, &
         ETDNEW, TDMIN, TDMAX, TDMIN_2nd_mc, TDMAX_2nd_mc, FIRORD, NOLIMI )
      implicit none 
      ! This sub calculates the limited face values TDADJ(1...SNGI) from the central 
      ! difference face values TDCEN(1...SNGI) using a NVD shceme.  
      ! The method is based on a new approach to limiting based on applying 
      ! limiting only when we have a local min or max.
      ! INCOME(1...SNGI)=1 for incomming to element ELE  else =0.
      ! LIBETA is the flux limiting parameter. 
      ! TDMAX(PELE)=maximum of the surrounding element values of element PELE.
      ! TDMIN(PELE)=minimum of the surrounding element values of element PELE.
      ! TDMIN_2nd_mc(PELE), TDMAX_2nd_mc(PELE) are the second minima and max 
      ! excluding (minus c) these values surrounding them. 
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
      REAL, DIMENSION( TOTELE ), intent( in ) :: ETDNEW, TDMIN, TDMAX, TDMIN_2nd_mc, TDMAX_2nd_mc
      LOGICAL, intent( in ) :: FIRORD, NOLIMI
      ! Local variables   
      REAL, PARAMETER :: TOLER=1.0E-10
      REAL :: TC_in, TD_in, TC_out, TD_out
      REAL :: TUPWIN_in, TUPWIN_out, TD_MAX_in, TD_MAX_out, TD_MIN_in, TD_MIN_out
      REAL :: TDLIM_in, TDLIM_out

      IF( NOLIMI ) THEN
         TDLIM = TDCEN
         RETURN
      ENDIF
      Conditional_FIRORD: IF( .NOT.FIRORD ) THEN ! Velocity is pointing into element

      Conditional_PELEOT: IF( PELEOT /= PELE ) THEN

         TC_in =ETDNEW( PELEOT )! upwind value for incomming (to cell PELE) vel
         TD_in =ETDNEW( PELE )  ! downwind value for incomming vel
         TC_out=ETDNEW( PELE )  ! upwind value for outgoing vel
         TD_out=ETDNEW( PELEOT )! downwind value for outgoing vel

         IF( ETDNEW( PELEOT ) > ETDNEW( PELE )) THEN
            TUPWIN_out = TDMIN( PELE )
            TD_MIN_out=TDMIN_2nd_mc( PELE )
            TD_MAX_out=TDMAX( PELE )   !1.E+10

            TUPWIN_in  = TDMAX( PELEOT )
            TD_MAX_in =TDMAX_2nd_mc( PELEOT ) 
            TD_MIN_in =TDMIN( PELEOT )!-1.E+10 used to bound everything
         ELSE
            TUPWIN_out = TDMAX( PELE )
            TD_MAX_out=TDMAX_2nd_mc( PELE )
            TD_MIN_out=TDMIN( PELE )   !-1.e+10

            TUPWIN_in  = TDMIN( PELEOT )
            TD_MIN_in =TDMIN_2nd_mc( PELEOT )
            TD_MAX_in =TDMAX( PELEOT ) !+1.e+10 used to bound everything
         ENDIF
! Calculate normalisation parameters for out going velocities *******
            CALL LIM_TC_2MAX2MIN(TDLIM_out, TC_out,TD_out, TUPWIN_out, TDCEN, TD_MIN_out,TD_MAX_out)

! Calculate normalisation parameters for incomming velocities *******
            CALL LIM_TC_2MAX2MIN(TDLIM_in, TC_in,TD_in, TUPWIN_in, TDCEN, TD_MIN_in,TD_MAX_in)

! overall TDLIM*******:
            TDLIM=INCOME*TDLIM_in + (1.0-INCOME)*TDLIM_out

      ELSE
! Next to boundary...
         TDLIM = ETDNEW( PELE )

      ENDIF Conditional_PELEOT

! Conditional_FIRORD: IF( .NOT.FIRORD ) THEN...
      ELSE

         ! Velocity is going out of element
!         TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 
         TDLIM = INCOME * ETDNEW( PELEOT ) + ( 1.0 - INCOME ) * ETDNEW( PELE )

      ENDIF Conditional_FIRORD

      RETURN

    END SUBROUTINE ONVDLIM_2nd



    SUBROUTINE LIM_TC_2MAX2MIN(TDLIM, TC,TD, TUPWIN, TDCEN, TD_MIN,TD_MAX)  
! Calculate normalisation parameters to obtain a 
! limited value of the face variable TDLIM based on out going velocities *******
! TC=value of the upwind cell.
! TD=value of the downwind cell.
! TUPWIN=value of the far field upwind cell.
! tdcen=central diff value (suggested) of the face value. 
! This sub calculates the downwind cell value TDELE_int and limits 
! it between [TD_MIN,TD_MAX].
      REAL, intent( inout ) :: TDLIM  
      REAL, intent( in ) :: TC,TD, TUPWIN, TDCEN, TD_MIN,TD_MAX
! local variables...
       REAL :: TDELE_int,tdcen_w,DENOOU,FTILOU,CTILOU,TDLIM_w,w

       TDELE_int = 2.*TC - TUPWIN
       TDELE_int = max(min(TDELE_int, TD_MAX),TD_MIN)
       w=min(max((TD-tdcen)/tolfun(TD - TC),0.0),1.0)
! value to be limited along new edge
       tdcen_w=w*TC +(1.-w)*tdele_int
       DENOOU = tolfun(TDELE_int - TUPWIN)
       FTILOU = ( TDCEN_w        - TUPWIN ) / DENOOU
       CTILOU = ( TC             - TUPWIN ) / DENOOU

       TDLIM_w=  TUPWIN + NVDFUNNEW( FTILOU, CTILOU, -1.0 ) * DENOOU
! project TDLIM_w back onto original edge.
       w=min(max((TDELE_int-TDLIM_w)/tolfun(TDELE_int - TC),0.0),1.0)
       TDLIM=w*TC + (1.-w)*TD
    RETURN

    END SUBROUTINE LIM_TC_2MAX2MIN




    SUBROUTINE ONVDLIM( TOTELE, &
         TDLIM, TDCEN, INCOME, PELE, PELEOT, &
         ETDNEW, TDMIN, TDMAX, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE )
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
      REAL, intent( in ) :: TDCEN, INCOME, COURANT_OR_MINUS_ONE
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
         !tdele = min( max( tdcen + 1.*( TDELE - tdcen ) , 0.) , 1. )
         DENOIN = TDELE - TUPWIN 

         IF( ABS( DENOIN ) < TOLER ) DENOIN = SIGN( TOLER, DENOIN )

         UCIN = ETDNEW( PELEOT )
         CTILIN = ( UCIN - TUPWIN ) / DENOIN

         ! Calculate normalisation parameters for out going velocities 
         TDELE = ETDNEW( PELEOT )
         !tdele =  min( max( tdcen + 1.*( TDELE - tdcen ) , 0.) , 1. )
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
         TDLIM= INCOME*( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, COURANT_OR_MINUS_ONE ) * DENOIN ) &
              + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, COURANT_OR_MINUS_ONE ) &
              * DENOOU )

      ENDIF Conditional_FIRORD

      RETURN

    END SUBROUTINE ONVDLIM




    SUBROUTINE ONVDLIM_ANO( TOTELE, &
         TDLIM, TDCEN, INCOME, PELE, PELEOT, &
         ETDNEW, TDMIN, TDMAX, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE, &
         SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, TUPWIND_MAT )
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
      INTEGER, intent( in ) :: TOTELE, NSMALL_COLM
      REAL, intent( inout ) :: TDLIM  
      REAL, intent( in ) :: TDCEN, INCOME, COURANT_OR_MINUS_ONE
      INTEGER, intent( in ) :: PELE, PELEOT
      REAL, DIMENSION( TOTELE ), intent( in ) :: ETDNEW, TDMIN, TDMAX
      LOGICAL, intent( in ) :: FIRORD, NOLIMI
      INTEGER, DIMENSION( TOTELE+1 ), intent( in ) :: SMALL_FINDRM
      INTEGER, DIMENSION( NSMALL_COLM ), intent( in ) :: SMALL_COLM
      REAL, DIMENSION( NSMALL_COLM ), intent( in ) :: TUPWIND_MAT 
      ! Local variables   
      REAL, PARAMETER :: TOLER=1.0E-10
      REAL :: UCIN, UCOU, TUPWIN, TUPWI2, TDELE, DENOIN, CTILIN, DENOOU, &
           CTILOU, FTILIN, FTILOU
      INTEGER :: COUNT

      IF( NOLIMI ) THEN
         TDLIM = TDCEN
         RETURN
      ENDIF

      Conditional_PELEOT: IF( PELEOT /= PELE ) THEN

!         IF( ETDNEW( PELEOT ) > ETDNEW( PELE )) THEN
!            TUPWIN = TDMAX( PELEOT )
!            TUPWI2 = TDMIN( PELE )
!         ELSE
!            TUPWIN = TDMIN( PELEOT )
!            TUPWI2 = TDMAX( PELE ) 
!         ENDIF

          DO COUNT=SMALL_FINDRM(PELEOT),SMALL_FINDRM(PELEOT+1)-1
             IF(SMALL_COLM(COUNT)==PELE) TUPWIN = TUPWIND_MAT(COUNT)
 !            IF(SMALL_COLM(COUNT)==PELE) TUPWI2 = TUPWIND_MAT(COUNT)
          END DO
          DO COUNT=SMALL_FINDRM(PELE),SMALL_FINDRM(PELE+1)-1
             IF(SMALL_COLM(COUNT)==PELEOT) TUPWI2 = TUPWIND_MAT(COUNT)
 !            IF(SMALL_COLM(COUNT)==PELEOT) TUPWIN = TUPWIND_MAT(COUNT)
          END DO


         ! Calculate normalisation parameters for incomming velocities 
         TDELE = ETDNEW( PELE )
         !tdele = min( max( tdcen + 1.*( TDELE - tdcen ) , 0.) , 1. )
         DENOIN = TDELE - TUPWIN 

         IF( ABS( DENOIN ) < TOLER ) DENOIN = SIGN( TOLER, DENOIN )

         UCIN = ETDNEW( PELEOT )
         CTILIN = ( UCIN - TUPWIN ) / DENOIN

         ! Calculate normalisation parameters for out going velocities 
         TDELE = ETDNEW( PELEOT )
         !tdele =  min( max( tdcen + 1.*( TDELE - tdcen ) , 0.) , 1. )
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
         TDLIM= INCOME*( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, COURANT_OR_MINUS_ONE ) * DENOIN ) &
              + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, COURANT_OR_MINUS_ONE ) &
              * DENOOU )

      ENDIF Conditional_FIRORD

      RETURN

    END SUBROUTINE ONVDLIM_ANO





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
      LOGICAL, PARAMETER :: DOWNWIND_EXTRAP =.TRUE.
      REAL :: TILDEUF, MAXUF

      ! For the region 0 < UC < 1 on the NVD, define the limiter
      IF( ( UC > 0.0 ) .AND. ( UC < 1.0 ) ) THEN
         ! For the interface tracking scheme, use Hyper-C (NOTE: at high Courant numbers
         !  limit is defined by XI, not the Hyper-C scheme.)

         IF( COURAT > 0.0 ) THEN
            IF(DOWNWIND_EXTRAP) THEN
! new method based on downwind extrapolation...
!               MAXUF = MAX( 0.0, UF )
               MAXUF = MAX( 0.0, UF, UC )
!               TILDEUF = MIN( 1.0, UC/ (10.0 * COURAT), MAXUF )
               TILDEUF = MIN( 1.0, UC * 10.0, MAXUF )
            ELSE
               !TILDEUF = MIN( 1.0, max( UC / COURAT, XI * UC ))
               ! halve the slope for now...
               TILDEUF = MIN( 1.0, max( UC / (2.0 * COURAT), XI * UC ))
            ENDIF
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
      use shape_functions_NDim
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
      IMPLICIT NONE
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

      !ewrite(3,*)' In SCVDETNX'

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

         !ewrite(3,*) 'CVNORMX(GI), CVNORMY(GI), CVNORMZ(GI):',CVNORMX(GI), CVNORMY(GI), CVNORMZ(GI)
         !ewrite(3,*) 'DXDLX,       DYDLX,       DZDLX:',DXDLX,       DYDLX,       DZDLX
         !ewrite(3,*) 'DXDLY,       DYDLY,       DZDLY:',DXDLY,       DYDLY,       DZDLY
         !ewrite(3,*) 'POSVGIX,     POSVGIY,     POSVGIZ:',POSVGIX,     POSVGIY,     POSVGIZ
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
      use shape_functions_NDim
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
      use shape_functions_NDim
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




    SUBROUTINE BETWEEN_ELE_SOLVE_DIF(UDIFF_SUF_STAB, &
                       DIFF_FOR_BETWEEN_U_ELE, DIFF_FOR_BETWEEN_U_ELE2, &
                       MAT_ELE, MAT_ELE2, U_SLOC2LOC,U_ILOC_OTHER_SIDE,  &
                       SBUFEN,SBCVNGI,U_NLOC,U_SNLOC,NDIM,NPHASE,GOT_OTHER_ELE) 
! Calculate the between element diffusion coefficients for stabaization scheme UDIFF_SUF_STAB.
      use matrix_operations
      implicit none      
      LOGICAL, PARAMETER :: FAST_AND_SIMP=.TRUE.
! If FAST_AND_SIMP use a simple mean to calculate the between element diffusion.
      INTEGER, intent( in ) :: SBCVNGI,U_NLOC,U_SNLOC,NDIM,NPHASE
      LOGICAL, intent( in ) :: GOT_OTHER_ELE
      REAL, DIMENSION( NPHASE,SBCVNGI,NDIM,NDIM  ), intent( out ) :: UDIFF_SUF_STAB
      REAL, DIMENSION( NPHASE,U_NLOC ), intent( in ) :: DIFF_FOR_BETWEEN_U_ELE,DIFF_FOR_BETWEEN_U_ELE2
      REAL, DIMENSION( U_NLOC,U_NLOC ), intent( in ) :: MAT_ELE, MAT_ELE2
      INTEGER, DIMENSION( U_SNLOC ), intent( in ) :: U_SLOC2LOC,U_ILOC_OTHER_SIDE
      REAL, DIMENSION( U_SNLOC,SBCVNGI ), intent( in ) :: SBUFEN
! Local variables...
      REAL, DIMENSION( : , : ), allocatable :: MAT_LOC_2ELES,MAT
      REAL, DIMENSION( : ), allocatable :: VECRHS_2ELES, DIFF
      REAL, DIMENSION( SBCVNGI ) :: DIFF_ADD_STAB
      INTEGER, DIMENSION( U_NLOC ) :: OTHER_SI2
      INTEGER, DIMENSION( 2*U_NLOC ) :: GLOB_NO
      REAL, DIMENSION( U_SNLOC ) :: DIFF_SUF
      INTEGER :: U_SILOC,U_ILOC,U_ILOC2,U_JLOC2,IGL,JGL,SGI,IPHASE,NLEN,IDIM
      LOGICAL :: GOTDEC


! initialize to zero to set the off diagonal terms to 0.0
      UDIFF_SUF_STAB=0.0

      IF(FAST_AND_SIMP) THEN
! use a simple average either side of the interface...
        IF(GOT_OTHER_ELE) THEN
           DO IPHASE=1,NPHASE
              DIFF_ADD_STAB( : )=0.0
              DO U_SILOC = 1, U_SNLOC
                 U_ILOC = U_SLOC2LOC( U_SILOC )
                 U_ILOC2= U_ILOC_OTHER_SIDE( U_SILOC ) 
! Calculate added stabilization diffusion DIFF_ADD_STAB( SGI,IDIM,IPHASE )
            DIFF_ADD_STAB( : )=DIFF_ADD_STAB( : )+SBUFEN(U_SILOC,:) &
       *0.5*(DIFF_FOR_BETWEEN_U_ELE(IPHASE,U_ILOC)+DIFF_FOR_BETWEEN_U_ELE2(IPHASE,U_ILOC2))
              END DO
              DIFF_ADD_STAB( : )=MAX(0.0,DIFF_ADD_STAB( : ))
              DO IDIM=1,NDIM
                 UDIFF_SUF_STAB(IPHASE,:,IDIM,IDIM)=DIFF_ADD_STAB( : )
              END DO
           END DO
        ELSE
           DO IPHASE=1,NPHASE
              DIFF_ADD_STAB( : )=0.0
              DO U_SILOC = 1, U_SNLOC
                 U_ILOC = U_SLOC2LOC( U_SILOC )
! Calculate added stabilization diffusion DIFF_ADD_STAB( SGI,IDIM,IPHASE )
            DIFF_ADD_STAB( : )=DIFF_ADD_STAB( : )+SBUFEN(U_SILOC,:)*DIFF_FOR_BETWEEN_U_ELE(IPHASE,U_ILOC)
              END DO
              DIFF_ADD_STAB( : )=MAX(0.0,DIFF_ADD_STAB( : ))
              DO IDIM=1,NDIM
                 UDIFF_SUF_STAB(IPHASE,:,IDIM,IDIM)=DIFF_ADD_STAB( : )
              END DO
           END DO
        ENDIF

! ENDOF IF(FAST_AND_SIMP) THEN...
      ELSE

      DO U_ILOC=1,U_NLOC
        GLOB_NO(U_ILOC)=U_ILOC
      END DO

      IF(GOT_OTHER_ELE) THEN
            OTHER_SI2(1:U_NLOC)=0
            DO U_SILOC = 1, U_SNLOC
                U_ILOC = U_SLOC2LOC( U_SILOC )
                U_ILOC2= U_ILOC_OTHER_SIDE( U_SILOC ) 
                OTHER_SI2(U_ILOC2)=U_ILOC
            END DO

            IGL=U_NLOC
            DO U_ILOC2=1,U_NLOC
              IF(OTHER_SI2(U_ILOC2)==0) THEN
                 IGL=IGL+1
                 GLOB_NO(U_ILOC2+U_NLOC)=IGL
              ELSE
                 GLOB_NO(U_ILOC2+U_NLOC)=OTHER_SI2(U_ILOC2)
              ENDIF
            END DO

            NLEN=2*U_NLOC-U_SNLOC
            ALLOCATE(VECRHS_2ELES(NLEN))
   
            ALLOCATE(MAT_LOC_2ELES(NLEN,NLEN))
            MAT_LOC_2ELES(:,:)=0.0
            MAT_LOC_2ELES(1:U_NLOC,1:U_NLOC)=MAT_ELE(1:U_NLOC,1:U_NLOC)
            DO U_ILOC2=1,U_NLOC 
               IGL=GLOB_NO(U_ILOC2+U_NLOC)
               DO U_JLOC2=1,U_NLOC 
                  JGL=GLOB_NO(U_JLOC2+U_NLOC)
                  MAT_LOC_2ELES(IGL,JGL)=MAT_LOC_2ELES(IGL,JGL)+MAT_ELE2(U_ILOC2,U_JLOC2)
                END DO
            END DO
      ELSE
            NLEN=U_NLOC
            ALLOCATE(VECRHS_2ELES(NLEN))

            ALLOCATE(MAT_LOC_2ELES(NLEN,NLEN))
            MAT_LOC_2ELES(:,:)=0.0
            MAT_LOC_2ELES(1:U_NLOC,1:U_NLOC)=MAT_ELE(1:U_NLOC,1:U_NLOC)
      ENDIF

      ALLOCATE(DIFF(NLEN))
      ALLOCATE(MAT(NLEN,NLEN))
      MAT=MAT_LOC_2ELES ! MAT is overwritten

      GOTDEC =.FALSE.
      DO IPHASE=1,NPHASE

         VECRHS_2ELES(:)   =0.0
         VECRHS_2ELES(1:U_NLOC)=DIFF_FOR_BETWEEN_U_ELE(IPHASE,1:U_NLOC)
         IF(GOT_OTHER_ELE) THEN

            DO U_ILOC2=1,U_NLOC 
               IGL=GLOB_NO(U_ILOC2+U_NLOC)
               VECRHS_2ELES(IGL)=VECRHS_2ELES(IGL)+DIFF_FOR_BETWEEN_U_ELE2(IPHASE,U_ILOC2)
            END DO
         ENDIF

            ! Solve MAT_LOC_2ELES *DIFF = VECRHS_2ELES  
         ! MAT is overwritten by decomposition
         CALL SMLINNGOT( MAT, DIFF, VECRHS_2ELES, NLEN, NLEN, GOTDEC)
         GOTDEC =.TRUE.
         DO U_SILOC=1,U_SNLOC
            U_ILOC = U_SLOC2LOC( U_SILOC )
            DIFF_SUF(U_SILOC)=DIFF(U_ILOC)
         END DO
! Calculate added stabilization diffusion DIFF_ADD_STAB( SGI,IDIM,IPHASE )
         DIFF_ADD_STAB( : )=0.0
         DO U_SILOC=1,U_SNLOC
            DIFF_ADD_STAB( : )=DIFF_ADD_STAB( : )+SBUFEN(U_SILOC,:)*DIFF_SUF(U_SILOC)
         END DO
! Make certain the diffusion is positive between the elements...
         DIFF_ADD_STAB( : )=MAX(0.0,DIFF_ADD_STAB( : ))
         DO IDIM=1,NDIM
            UDIFF_SUF_STAB(IPHASE,1:SBCVNGI,IDIM,IDIM  )=UDIFF_SUF_STAB(IPHASE,1:SBCVNGI,IDIM,IDIM  ) &
               + DIFF_ADD_STAB( 1:SBCVNGI )
         END DO
! ENDOF DO IPHASE=1,NPHASE...
      END DO

! ENDOF IF(FAST_AND_SIMP) THEN ELSE...
      ENDIF
      
    END SUBROUTINE BETWEEN_ELE_SOLVE_DIF




    SUBROUTINE DIFFUS_CAL_COEFF(DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX,  &
         CV_NLOC, MAT_NLOC, CV_NONODS, NPHASE, TOTELE, MAT_NONODS, MAT_NDGLN, &
         SMATFEN, SCVFEN, SCVNGI, GI, IPHASE, NDIM, TDIFFUSION, DIFF_GI_ADDED, &
         HDC, &
         T_CV_NODJ_IPHA, T_CV_NODI_IPHA, &
         TOLD_CV_NODJ_IPHA, TOLD_CV_NODI_IPHA, &
         ELE, ELE2, CVNORMX, CVNORMY, CVNORMZ, &
         DTX_ELE, DTY_ELE, DTZ_ELE, DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE, &
         SELE, STOTEL, WIC_T_BC, WIC_T_BC_DIRICHLET, CV_OTHER_LOC, MAT_OTHER_LOC )
      ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
      ! based on a non-linear method and a non-oscillating scheme.
      IMPLICIT NONE
      INTEGER, intent( in ) :: CV_NLOC, MAT_NLOC, CV_NONODS,NPHASE, TOTELE, MAT_NONODS, &
           &                   SCVNGI, GI, IPHASE, NDIM, ELE, ELE2, &
           &                   SELE, STOTEL, WIC_T_BC_DIRICHLET
      REAL, intent( in ) :: HDC, T_CV_NODJ_IPHA, T_CV_NODI_IPHA, &
           &                TOLD_CV_NODJ_IPHA, TOLD_CV_NODI_IPHA
      REAL, intent( inout ) :: DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
      INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) ::MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::WIC_T_BC
      INTEGER, DIMENSION( CV_NLOC ), intent( in ) ::CV_OTHER_LOC
      INTEGER, DIMENSION( MAT_NLOC ), intent( in ) ::MAT_OTHER_LOC
      REAL, DIMENSION( MAT_NLOC,SCVNGI ), intent( in ) :: SMATFEN
      REAL, DIMENSION( CV_NLOC,SCVNGI ), intent( in ) :: SCVFEN
      REAL, DIMENSION( MAT_NONODS,NDIM,NDIM,NPHASE ), intent( in ) :: TDIFFUSION
      REAL, DIMENSION( NDIM, NDIM ), intent( in ) :: DIFF_GI_ADDED
      REAL, DIMENSION( CV_NLOC, NPHASE, TOTELE ), intent( in ) :: DTX_ELE,DTY_ELE,DTZ_ELE, &
           &                                                      DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
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
             !ewrite(3,*)'MAT_KLOC,cv_nloc,mat_nloc,MAT_NODK,mat_nonods,cv_nonods:', &
             !         MAT_KLOC,cv_nloc,mat_nloc,MAT_NODK,mat_nonods,cv_nonods
            DIFF_GI( 1:NDIM , 1:NDIM ) = DIFF_GI( 1:NDIM , 1:NDIM ) &
                 + SMATFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK, 1:NDIM , 1:NDIM , IPHASE )
!                 + SCVFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK, 1:NDIM , 1:NDIM , IPHASE )
         END DO
         DIFF_GI( 1:NDIM , 1:NDIM ) = DIFF_GI( 1:NDIM , 1:NDIM )+DIFF_GI_ADDED( 1:NDIM , 1:NDIM )

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
!                  DIFF_GI2( 1:NDIM, 1:NDIM )=DIFF_GI2( 1:NDIM, 1:NDIM ) +SCVFEN( MAT_KLOC, GI ) &
                  DIFF_GI2( 1:NDIM, 1:NDIM )=DIFF_GI2( 1:NDIM, 1:NDIM ) +SMATFEN( MAT_KLOC, GI ) &
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


!          EWRITE(3,*)'CV_NODI_IPHA,CV_NODJ_IPHA:',CV_NODI_IPHA,CV_NODJ_IPHA
         DIFF_COEF_DIVDX    = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX, N_DOT_DKDT / &
              TOLFUN( T_CV_NODJ_IPHA  - T_CV_NODI_IPHA )  )
!              TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA )) )
         DIFF_COEFOLD_DIVDX = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX, N_DOT_DKDTOLD /  &
              TOLFUN( TOLD_CV_NODJ_IPHA  - TOLD_CV_NODI_IPHA )  )
!              TOLFUN( TOLD( CV_NODJ_IPHA ) - TOLD( CV_NODI_IPHA )))

         ! Make sure the diffusion has an upper bound...       
         DIFF_COEF_DIVDX    = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX, DIFF_COEF_DIVDX )
         DIFF_COEFOLD_DIVDX = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX, DIFF_COEFOLD_DIVDX )

      END IF Cond_ZerDiff

      !ewrite(3,*)'HDC, TOLFUN, DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX, DIFF_GI, DIFF_GI2:', &
      !     HDC, TOLFUN( T_CV_NODJ_IPHA - T_CV_NODI_IPHA ), DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX, &
      !     '::', DIFF_GI, '::', DIFF_GI2

      DEALLOCATE( DIFF_GI, DIFF_GI2 )

      RETURN            

    END SUBROUTINE DIFFUS_CAL_COEFF




    SUBROUTINE DIFFUS_CAL_COEFF_SURFACE(DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX,  &
         CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE, TOTELE, MAT_NONODS,MAT_NDGLN, &
         SBCVFEN,SBCVNGI,SGI,IPHASE,NDIM,TDIFFUSION,DIFF_GI_ADDED, &
         HDC, &
         T_CV_NODJ_IPHA, T_CV_NODI_IPHA, &
         TOLD_CV_NODJ_IPHA, TOLD_CV_NODI_IPHA, &
         ELE,ELE2,SNORMXN,SNORMYN,SNORMZN, &
         DTX_ELE,DTY_ELE,DTZ_ELE,DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE, &
         SELE,STOTEL,WIC_T_BC,WIC_T_BC_DIRICHLET, MAT_OTHER_LOC,CV_SLOC2LOC )
      ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
      ! based on a non-linear method and a non-oscillating scheme.
      IMPLICIT NONE
      INTEGER, intent( in ) :: CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE, TOTELE, MAT_NONODS, &
           SBCVNGI, SGI, IPHASE, NDIM,ELE, ELE2, &
           SELE, STOTEL, WIC_T_BC_DIRICHLET
      REAL, intent( in ) :: HDC,T_CV_NODJ_IPHA, T_CV_NODI_IPHA, &
                                TOLD_CV_NODJ_IPHA, TOLD_CV_NODI_IPHA
      REAL, intent( inout ) :: DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
      INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) ::MAT_NDGLN
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) ::WIC_T_BC
      INTEGER, DIMENSION( MAT_NLOC ), intent( in ) ::MAT_OTHER_LOC
      INTEGER, DIMENSION( CV_SNLOC ), intent( in ) ::CV_SLOC2LOC
      REAL, DIMENSION( CV_SNLOC, SBCVNGI  ), intent( in ) :: SBCVFEN
      REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE  ), intent( in ) :: TDIFFUSION
      REAL, DIMENSION( NDIM, NDIM ), intent( in ) :: DIFF_GI_ADDED
      REAL, DIMENSION( CV_NLOC, NPHASE, TOTELE ), intent( in ) :: DTX_ELE,DTY_ELE,DTZ_ELE, &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
      REAL, DIMENSION( SBCVNGI ), intent( in ) :: SNORMXN,SNORMYN,SNORMZN

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
      INTEGER :: CV_SKLOC,CV_KLOC,CV_KLOC2,MAT_KLOC,MAT_KLOC2,MAT_NODK,MAT_NODK2
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

         DO CV_SKLOC = 1, CV_SNLOC
            CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
            
            DTDX_GI = DTDX_GI + SBCVFEN(CV_SKLOC,SGI) * DTX_ELE(CV_KLOC,IPHASE,ELE)
            DTDY_GI = DTDY_GI + SBCVFEN(CV_SKLOC,SGI) * DTY_ELE(CV_KLOC,IPHASE,ELE)
            DTDZ_GI = DTDZ_GI + SBCVFEN(CV_SKLOC,SGI) * DTZ_ELE(CV_KLOC,IPHASE,ELE)

            DTOLDDX_GI = DTOLDDX_GI + SBCVFEN(CV_SKLOC,SGI) * DTOLDX_ELE(CV_KLOC,IPHASE,ELE)
            DTOLDDY_GI = DTOLDDY_GI + SBCVFEN(CV_SKLOC,SGI) * DTOLDY_ELE(CV_KLOC,IPHASE,ELE)
            DTOLDDZ_GI = DTOLDDZ_GI + SBCVFEN(CV_SKLOC,SGI) * DTOLDZ_ELE(CV_KLOC,IPHASE,ELE)
         END DO

         DIFF_GI(:,:) = 0.0
         DO CV_SKLOC = 1, CV_SNLOC
            MAT_KLOC = CV_SLOC2LOC( CV_SKLOC )
            MAT_NODK = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_KLOC )
             !ewrite(3,*)'MAT_KLOC,cv_nloc,mat_nloc,MAT_NODK,mat_nonods,cv_nonods:', &
             !         MAT_KLOC,cv_nloc,mat_nloc,MAT_NODK,mat_nonods,cv_nonods
            DIFF_GI( 1:NDIM , 1:NDIM ) = DIFF_GI( 1:NDIM , 1:NDIM ) &
                 + SBCVFEN(CV_SKLOC,SGI) * TDIFFUSION( MAT_NODK, 1:NDIM , 1:NDIM , IPHASE )
         END DO
         DIFF_GI(:,:) = max(0.0, DIFF_GI(:,:) ) 
         DIFF_GI( 1:NDIM , 1:NDIM ) = DIFF_GI( 1:NDIM , 1:NDIM )+DIFF_GI_ADDED( 1:NDIM , 1:NDIM )


         N_DOT_DKDT=SNORMXN(SGI)*(DIFF_GI(1,1)*DTDX_GI+DIFF_GI(1,2)*DTDY_GI+DIFF_GI(1,3)*DTDZ_GI) &
                   +SNORMYN(SGI)*(DIFF_GI(2,1)*DTDX_GI+DIFF_GI(2,2)*DTDY_GI+DIFF_GI(2,3)*DTDZ_GI) &
                   +SNORMZN(SGI)*(DIFF_GI(3,1)*DTDX_GI+DIFF_GI(3,2)*DTDY_GI+DIFF_GI(3,3)*DTDZ_GI)

         N_DOT_DKDTOLD=SNORMXN(SGI)*(DIFF_GI(1,1)*DTOLDDX_GI  &
              +DIFF_GI(1,2)*DTOLDDY_GI+DIFF_GI(1,3)*DTOLDDZ_GI) &
              +SNORMYN(SGI)*(DIFF_GI(2,1)*DTOLDDX_GI  &
              +DIFF_GI(2,2)*DTOLDDY_GI+DIFF_GI(2,3)*DTOLDDZ_GI) &
              +SNORMZN(SGI)*(DIFF_GI(3,1)*DTOLDDX_GI  &
              +DIFF_GI(3,2)*DTOLDDY_GI+DIFF_GI(3,3)*DTOLDDZ_GI)

         ! This is the minimum diffusion...
         DIFF_STAND_DIVDX=( ABS(SNORMXN(SGI))*DIFF_GI(1,1) &
              +ABS(SNORMYN(SGI))*DIFF_GI(2,2) + ABS(SNORMZN(SGI))*DIFF_GI(3,3) ) /HDC

         Conditional_MAT_DISOPT_ELE2: IF( ( ELE2 /= 0 ).AND.( ELE2 /= ELE) ) THEN
            DIFF_GI2(:,:) = 0.0
            DO CV_SKLOC = 1, CV_SNLOC
               MAT_KLOC = CV_SLOC2LOC( CV_SKLOC )
               MAT_KLOC2 = MAT_OTHER_LOC( MAT_KLOC )
               IF( MAT_KLOC2 /= 0 ) THEN
                  MAT_NODK2 = MAT_NDGLN(( ELE2 - 1 ) * MAT_NLOC + MAT_KLOC2 )
                  DIFF_GI2( 1:NDIM, 1:NDIM )=DIFF_GI2( 1:NDIM, 1:NDIM ) +SBCVFEN(CV_SKLOC,SGI) &
                       *TDIFFUSION(MAT_NODK2, 1:NDIM, 1:NDIM ,IPHASE)
               ENDIF
            END DO
            DIFF_GI2(:,:) = max(0.0, DIFF_GI2(:,:) ) 
            DIFF_GI2( 1:NDIM, 1:NDIM ) = DIFF_GI2( 1:NDIM , 1:NDIM ) + DIFF_GI_ADDED( 1:NDIM , 1:NDIM )

            DTDX_GI2 = 0.0
            DTDY_GI2 = 0.0
            DTDZ_GI2 = 0.0
            DTOLDDX_GI2 = 0.0
            DTOLDDY_GI2 = 0.0
            DTOLDDZ_GI2 = 0.0
            DO CV_SKLOC = 1, CV_SNLOC
               CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
               CV_KLOC2 = MAT_OTHER_LOC( CV_KLOC )
               IF(CV_KLOC2 /= 0 )THEN
                  DTDX_GI2 = DTDX_GI2 + SBCVFEN(CV_SKLOC,SGI) * DTX_ELE(CV_KLOC2,IPHASE,ELE2)
                  DTDY_GI2 = DTDY_GI2 + SBCVFEN(CV_SKLOC,SGI) * DTY_ELE(CV_KLOC2,IPHASE,ELE2)
                  DTDZ_GI2 = DTDZ_GI2 + SBCVFEN(CV_SKLOC,SGI) * DTZ_ELE(CV_KLOC2,IPHASE,ELE2)

                  DTOLDDX_GI2 = DTOLDDX_GI2 + SBCVFEN(CV_SKLOC,SGI) * DTOLDX_ELE(CV_KLOC2,IPHASE,ELE2)
                  DTOLDDY_GI2 = DTOLDDY_GI2 + SBCVFEN(CV_SKLOC,SGI) * DTOLDY_ELE(CV_KLOC2,IPHASE,ELE2)
                  DTOLDDZ_GI2 = DTOLDDZ_GI2 + SBCVFEN(CV_SKLOC,SGI) * DTOLDZ_ELE(CV_KLOC2,IPHASE,ELE2)
               ENDIF
            END DO
            N_DOT_DKDT2=SNORMXN(SGI)*(DIFF_GI2(1,1)*DTDX_GI2  &
                 +DIFF_GI2(1,2)*DTDY_GI2+DIFF_GI2(1,3)*DTDZ_GI2) &
                 +SNORMYN(SGI)*(DIFF_GI2(2,1)*DTDX_GI2  &
                 +DIFF_GI2(2,2)*DTDY_GI2+DIFF_GI2(2,3)*DTDZ_GI2) &
                 +SNORMZN(SGI)*(DIFF_GI2(3,1)*DTDX_GI2  &
                 +DIFF_GI2(3,2)*DTDY_GI2+DIFF_GI2(3,3)*DTDZ_GI2) 

            N_DOT_DKDTOLD2=SNORMXN(SGI)*(DIFF_GI2(1,1)*DTOLDDX_GI2  &
                 +DIFF_GI2(1,2)*DTOLDDY_GI2+DIFF_GI2(1,3)*DTOLDDZ_GI2) &
                 +SNORMYN(SGI)*(DIFF_GI2(2,1)*DTOLDDX_GI2  &
                 +DIFF_GI2(2,2)*DTOLDDY_GI2+DIFF_GI2(2,3)*DTOLDDZ_GI2) &
                 +SNORMZN(SGI)*(DIFF_GI2(3,1)*DTOLDDX_GI2  &
                 +DIFF_GI2(3,2)*DTOLDDY_GI2+DIFF_GI2(3,3)*DTOLDDZ_GI2) 

            ! This is the minimum diffusion...
            DIFF_STAND_DIVDX2 = ( ABS(SNORMXN(SGI))*DIFF_GI2(1,1) &
                 +ABS(SNORMYN(SGI))*DIFF_GI2(2,2) + ABS(SNORMZN(SGI))*DIFF_GI2(3,3) ) /HDC

            N_DOT_DKDT = 0.5*( N_DOT_DKDT + N_DOT_DKDT2 )
            N_DOT_DKDTOLD= 0.5*( N_DOT_DKDTOLD + N_DOT_DKDTOLD2 )

            !   ewrite(3,*)'DToldDX_GI,DToldDX_GI2,N_DOT_DKDTOLD,N_DOT_DKDTOLD2:',  &
            !            DToldDX_GI,DToldDX_GI2,N_DOT_DKDTOLD,N_DOT_DKDTOLD2

            ! This is the minimum diffusion...
            DIFF_STAND_DIVDX = MIN( DIFF_STAND_DIVDX, DIFF_STAND_DIVDX2 ) 

         ENDIF Conditional_MAT_DISOPT_ELE2


!          EWRITE(3,*)'CV_NODI_IPHA,CV_NODJ_IPHA:',CV_NODI_IPHA,CV_NODJ_IPHA
         DIFF_COEF_DIVDX    = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX, N_DOT_DKDT / &
              TOLFUN( T_CV_NODJ_IPHA  - T_CV_NODI_IPHA )  )
!              TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA )) )
         DIFF_COEFOLD_DIVDX = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX, N_DOT_DKDTOLD /  &
              TOLFUN( TOLD_CV_NODJ_IPHA  - TOLD_CV_NODI_IPHA )  )
!              TOLFUN( TOLD( CV_NODJ_IPHA ) - TOLD( CV_NODI_IPHA )))

         ! Make sure the diffusion has an upper bound...       
         DIFF_COEF_DIVDX    = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX, DIFF_COEF_DIVDX )
         DIFF_COEFOLD_DIVDX = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX, DIFF_COEFOLD_DIVDX )

      END IF Cond_ZerDiff

      !    ewrite(3,*)'HDC,DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX:', &
      !             HDC,DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
      !    ewrite(3,*)'TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA )):',  &
      !             TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA ))

      DEALLOCATE( DIFF_GI, DIFF_GI2 )

      RETURN

    END SUBROUTINE DIFFUS_CAL_COEFF_SURFACE





    SUBROUTINE GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ,INCOME, NDOTQOLD, INCOMEOLD, &
         HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
         T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, &
         U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
         CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
         SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
         SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
         UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
         VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
         MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
         FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
         TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC,  LIMIT_USE_2ND)
      ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
      IMPLICIT NONE
      INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
           CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
           SELE, U_SNLOC, STOTEL, WIC_U_BC_DIRICHLET, CV_ELE_TYPE, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, &
           NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
           IN_ELE_UPWIND, DG_ELE_UPWIND
      REAL, intent( in ) :: HDC, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI
      REAL, intent( inout ) :: NDOTQNEW,NDOTQ, INCOME, NDOTQOLD, INCOMEOLD
      INTEGER, DIMENSION( TOTELE*U_NLOC ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( U_NLOC ), intent( in ) :: U_OTHER_LOC
      INTEGER, DIMENSION( U_SNLOC ), intent( in ) :: U_SLOC2LOC
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_U_BC
      INTEGER, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
           TOLDMAX_NOD
      REAL, DIMENSION( U_NLOC, SCVNGI  ), intent( in ) :: SUFEN
      REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( in ) :: T, TOLD, FEMT, FEMTOLD, DEN, DENOLD
      REAL, DIMENSION( U_NONODS * NPHASE  ), intent( in ) :: U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( SCVNGI ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
      REAL, DIMENSION( STOTEL * U_SNLOC * NPHASE ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE, NDIM ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( U_NLOC ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
      REAL, DIMENSION( TOTELE  ), intent( in ) :: VOLFRA_PORE
      REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFEN
      INTEGER, DIMENSION( TOTELE*CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( CV_NLOC ), intent( in ) :: CV_OTHER_LOC
      INTEGER, DIMENSION( CV_SNLOC ), intent( in ) :: CV_SLOC2LOC
      REAL, DIMENSION( CV_NONODS  ), intent( in ) :: MASS_CV
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) :: MAT_NDGLN
      REAL, DIMENSION( CV_NONODS*NPHASE  ), intent( inout ) :: UP_WIND_NOD
      REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN, TOLDMIN, TMAX, TOLDMAX
      logical, intent( in ) :: LIMIT_USE_2ND
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC
      ! local variables
      character( len = option_path_len ) :: overlapping_path 
      logical :: is_overlapping   
      INTEGER :: U_NLOC_LEV,U_KLOC_LEV,U_KLOC,U_NODK_IPHA, U_KLOC2, U_NODK2_IPHA

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.

      IF( is_overlapping ) THEN
         ! For overlapping basis function approach.
         CALL GET_INT_VEL_OVERLAP( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
              HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
              T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
              CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
              CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
              SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
              SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
              UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
              VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
              MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
              FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
              TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
              IN_ELE_UPWIND, DG_ELE_UPWIND, &
              TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC,  LIMIT_USE_2ND)
      ELSE
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
              IN_ELE_UPWIND, DG_ELE_UPWIND, &
              TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC,  LIMIT_USE_2ND)
      ENDIF

      ! Calculate NDOTQNEW from NDOTQ
      NDOTQNEW=NDOTQ ! initialize it like this so that it contains the b.c's
      U_NLOC_LEV =U_NLOC /CV_NLOC
      DO U_KLOC = 1, U_NLOC
               U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) + (IPHASE-1)*U_NONODS
               NDOTQNEW=NDOTQNEW &
               + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE(U_KLOC) * ( U( U_NODK_IPHA ) - NU( U_NODK_IPHA ) ) * CVNORMX(GI) &
               + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE(U_KLOC) * ( V( U_NODK_IPHA ) - NV( U_NODK_IPHA ) ) * CVNORMY(GI) &
               + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE(U_KLOC) * ( W( U_NODK_IPHA ) - NW( U_NODK_IPHA ) ) * CVNORMZ(GI)
      END DO

      IF((ELE2 /= 0).AND.(ELE2 /= ELE)) THEN
         ! We have a discontinuity between elements so integrate along the face...
         DO U_KLOC = 1, U_NLOC
            U_KLOC2=U_OTHER_LOC(U_KLOC)
            IF(U_KLOC2 /= 0) THEN
               U_NODK2_IPHA = U_NDGLN(( ELE2 - 1 ) * U_NLOC + U_KLOC2 ) +(IPHASE-1)*U_NONODS
               NDOTQNEW=NDOTQNEW &
                    + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE2(U_KLOC2) * ( U( U_NODK2_IPHA ) - NU( U_NODK2_IPHA ) ) * CVNORMX(GI) &
                    + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE2(U_KLOC2) * ( V( U_NODK2_IPHA ) - NV( U_NODK2_IPHA ) ) * CVNORMY(GI) &
                    + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE2(U_KLOC2) * ( W( U_NODK2_IPHA ) - NW( U_NODK2_IPHA ) ) * CVNORMZ(GI)
            END IF
         END DO
      END IF

      RETURN
    END SUBROUTINE GET_INT_VEL





    SUBROUTINE GET_INT_VEL_OVERLAP( NPHASE, NDOTQ, INCOME, NDOTQOLD, INCOMEOLD, &
         HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
         T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, NU, NV, NW, NUOLD, NVOLD, NWOLD, &
         CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
         CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
         SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
         SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
         UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
         VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
         MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
         FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
         TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC,  LIMIT_USE_2ND)
      !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===
      ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
      ! it assumes an overlapping decomposition approach for velocity. 
      IMPLICIT NONE
      INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, CV_SNLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
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
      REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE, NDIM ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( U_NLOC ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
      REAL, DIMENSION( TOTELE  ), intent( in ) :: VOLFRA_PORE
      REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFEN
      INTEGER, DIMENSION( TOTELE*CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( CV_NLOC ), intent( in ) :: CV_OTHER_LOC
      INTEGER, DIMENSION( CV_SNLOC ), intent( in ) :: CV_SLOC2LOC
      REAL, DIMENSION( CV_NONODS  ), intent( in ) :: MASS_CV
      REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, DIMENSION( TOTELE*MAT_NLOC ), intent( in ) :: MAT_NDGLN
      REAL, DIMENSION( CV_NONODS*NPHASE  ), intent( inout ) :: UP_WIND_NOD
      REAL, DIMENSION( CV_NONODS * NPHASE  ), intent( inout ) :: TMIN, TOLDMIN, TMAX, TOLDMAX
      logical, intent( in ) :: LIMIT_USE_2ND
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC

      ! Local variables
      REAL :: UDGI,VDGI,WDGI,UOLDDGI,VOLDDGI,WOLDDGI,  &
           UDGI2,VDGI2,WDGI2,UOLDDGI2,VOLDDGI2,WOLDDGI2, DT_I,DT_J,DTOLD_I,DTOLD_J, &
           UDGI_INT,VDGI_INT,WDGI_INT,UOLDDGI_INT,VOLDDGI_INT,WOLDDGI_INT, &
           NDOTQ_INT,NDOTQOLD_INT,FEMTOLDGI2,NDOTQ2,NDOTQOLD2, &
           FEMTOLDGI_IPHA, OVER_RELAX, ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
           GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA, V_NODI, G_NODI, V_NODJ, G_NODJ, &
           GEOMTOLDGI_IPHA, W_UPWIND, INCOME3, INCOME4, INCOMEOLD3, INCOMEOLD4, &
           LIMT3, LIMTOLD3, FEMTGI_IPHA, GEOMTGI_IPHA, UPWIND_FRAC
      REAL :: NVEC(3),SUF_SIG_DIAGTEN_BC_GI(3), UGI_TMP(3)
      INTEGER :: U_KLOC,U_NODK,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_SKLOC, &
           U_SNODK,U_SNODK_IPHA, II,  &
           U_KLOC_LEV, U_NLOC_LEV, U_SKLOC_LEV, U_SNLOC_LEV, CV_KLOC, CV_KNOD, &
           CV_KLOC2, CV_KNOD2, CV_NODI, CV_NODJ, IDIM, JDIM, IJ, MAT_NODI, MAT_NODJ, &
           CV_SKLOC, CV_SNODK, CV_SNODK_IPHA
      LOGICAL :: CONSERV, MAX_OPER, SAT_BASED
      ! IN_ELE_UPWIND=1 switches on upwinding within and element (=3 recommended).
      !    INTEGER, PARAMETER :: IN_ELE_UPWIND = 3
      !    INTEGER, PARAMETER :: IN_ELE_UPWIND = 2
      ! DG_ELE_UPWIND=1 switches on upwinding between elements (=3 recommended).
      !    INTEGER, PARAMETER :: DG_ELE_UPWIND = 3
      !    INTEGER, PARAMETER :: DG_ELE_UPWIND = 2
      LOGICAL, PARAMETER :: LIM_VOL_ADJUST2 = .true.
      LOGICAL :: RESET_STORE, LIM_VOL_ADJUST
      REAL :: TMIN_STORE, TMAX_STORE, TOLDMIN_STORE, TOLDMAX_STORE
      ! coefficients for this element ELE

      ! The adjustment method for variable CV volumes is not ready for new limiter method...
      LIM_VOL_ADJUST=LIM_VOL_ADJUST2.AND.(.NOT.LIMIT_USE_2ND)

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
               U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) + (IPHASE-1)*U_NONODS
               UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
               VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
               WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )
               UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK_IPHA )
               VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK_IPHA )
               WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK_IPHA )
            END DO
! Here we assume that sigma_out/sigma_in is a diagonal matrix 
! which effectively assumes that the anisotropy just inside the domain 
! is the same as just outside the domain. 
! Multiply by a normalized sigma tensor so that we use the 
! sigma from just outside the boundary:

            SUF_SIG_DIAGTEN_BC_GI=0.0
            DO CV_SKLOC = 1, CV_SNLOC
               CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
               IF(CV_KLOC==CV_ILOC) THEN
                  CV_SNODK = ( SELE - 1 ) * CV_SNLOC + CV_SKLOC
                  CV_SNODK_IPHA = CV_SNODK + ( IPHASE - 1 ) * STOTEL*CV_SNLOC
                  SUF_SIG_DIAGTEN_BC_GI(1:NDIM)=SUF_SIG_DIAGTEN_BC( CV_SNODK_IPHA,1:NDIM)
               ENDIF
            END DO
! Only modify boundary velocity for incomming velocity...
            IF(UDGI*CVNORMX(GI)+VDGI*CVNORMY(GI)+WDGI*CVNORMZ(GI).LT.0.0) THEN ! Incomming...
               UGI_TMP = SUF_SIG_DIAGTEN_BC_GI(1:3) * (/UDGI, VDGI, WDGI/)
               UDGI=UGI_TMP(1) ; VDGI=UGI_TMP(2) ; WDGI=UGI_TMP(3)
            ENDIF
            
            IF(UOLDDGI*CVNORMX(GI)+VOLDDGI*CVNORMY(GI)+WOLDDGI*CVNORMZ(GI).LT.0.0) THEN ! Incomming...
               UGI_TMP = SUF_SIG_DIAGTEN_BC_GI(1:3) * (/UOLDDGI, VOLDDGI, WOLDDGI/)
               UOLDDGI=UGI_TMP(1) ; VOLDDGI=UGI_TMP(2) ; WOLDDGI=UGI_TMP(3)
            ENDIF

            UGI_COEF_ELE=0.0
            VGI_COEF_ELE=0.0
            WGI_COEF_ELE=0.0
            DO U_KLOC_LEV = 1, U_NLOC_LEV
               U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
               IF(UDGI*CVNORMX(GI)+VDGI*CVNORMY(GI)+WDGI*CVNORMZ(GI).LT.0.0) THEN ! Incomming...
                  UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+1.0*SUF_SIG_DIAGTEN_BC_GI(1)
                  VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+1.0*SUF_SIG_DIAGTEN_BC_GI(2)
                  WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+1.0*SUF_SIG_DIAGTEN_BC_GI(3)
               ELSE
                  UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+1.0
                  VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+1.0
                  WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+1.0
               ENDIF
            END DO

         ELSE ! Specified vel bc.
            UDGI = 0.0
            VDGI = 0.0
            WDGI = 0.0
            UOLDDGI = 0.0
            VOLDDGI = 0.0
            WOLDDGI = 0.0
            UGI_COEF_ELE=0.0
            VGI_COEF_ELE=0.0
            WGI_COEF_ELE=0.0 
            DO U_SKLOC_LEV = 1, U_SNLOC_LEV
               U_SKLOC = (CV_ILOC-1)*U_SNLOC_LEV + U_SKLOC_LEV
               U_KLOC = U_SLOC2LOC( U_SKLOC )
               U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
               U_NODK_IPHA = U_NODK + ( IPHASE - 1 ) * U_NONODS
               U_SNODK = ( SELE - 1 ) * U_SNLOC + U_SKLOC 
               U_SNODK_IPHA = U_SNODK + ( IPHASE - 1 ) * STOTEL*U_SNLOC
               IF(WIC_U_BC(SELE+(IPHASE-1)*STOTEL) == 10) THEN
                  UDGI = UDGI + SUFEN( U_KLOC, GI ) * 0.5 *(NU( U_NODK_IPHA )+SUF_U_BC( U_SNODK_IPHA ))
                  VDGI = VDGI + SUFEN( U_KLOC, GI ) * 0.5 *(NV( U_NODK_IPHA )+SUF_V_BC( U_SNODK_IPHA ))
                  WDGI = WDGI + SUFEN( U_KLOC, GI ) * 0.5 *(NW( U_NODK_IPHA )+SUF_W_BC( U_SNODK_IPHA ))
                  UGI_COEF_ELE(U_KLOC)=UGI_COEF_ELE(U_KLOC)+0.5
                  VGI_COEF_ELE(U_KLOC)=VGI_COEF_ELE(U_KLOC)+0.5
                  WGI_COEF_ELE(U_KLOC)=WGI_COEF_ELE(U_KLOC)+0.5
                  UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * 0.5 * (NUOLD( U_NODK_IPHA )+SUF_U_BC( U_SNODK_IPHA ))
                  VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * 0.5 * (NVOLD( U_NODK_IPHA )+SUF_V_BC( U_SNODK_IPHA ))
                  WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * 0.5 * (NWOLD( U_NODK_IPHA )+SUF_W_BC( U_SNODK_IPHA ))
               ELSE
                  UDGI = UDGI + SUFEN( U_KLOC, GI ) * SUF_U_BC( U_SNODK_IPHA )
                  VDGI = VDGI + SUFEN( U_KLOC, GI ) * SUF_V_BC( U_SNODK_IPHA )
                  WDGI = WDGI + SUFEN( U_KLOC, GI ) * SUF_W_BC( U_SNODK_IPHA )
                  UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * SUF_U_BC( U_SNODK_IPHA ) 
                  VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * SUF_V_BC( U_SNODK_IPHA ) 
                  WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * SUF_W_BC( U_SNODK_IPHA )
               END IF
            END DO
         END IF

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

            NDOTQ = 0.5 * ( CVNORMX(GI) * (UDGI+UDGI2) + CVNORMY(GI) * (VDGI+VDGI2)  &
                 + CVNORMZ(GI) * (WDGI+WDGI2) )
            NDOTQOLD = 0.5 * ( CVNORMX(GI) * (UOLDDGI+UOLDDGI2) + CVNORMY(GI) * (VOLDDGI+VOLDDGI2) &
                 + CVNORMZ(GI) * (WOLDDGI+WOLDDGI2) )

            CV_NODI = CV_NODI_IPHA - (IPHASE-1)*CV_NONODS
            CV_NODJ = CV_NODJ_IPHA - (IPHASE-1)*CV_NONODS

            IF(IN_ELE_UPWIND==1) THEN

               IF(NDOTQ < 0.0) THEN
                  INCOME=1.0
               ELSE
                  INCOME=0.0
               END IF

               IF(NDOTQOLD < 0.0) THEN
                  INCOMEOLD=1.0
               ELSE
                  INCOMEOLD=0.0
               END IF

            ELSE IF(IN_ELE_UPWIND==2) THEN ! the best

               UPWIND_FRAC=0.8
               IF(NDOTQ < 0.0) THEN
                  INCOME=0.8
                  !INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               ELSE
                  INCOME=0.8
                  !INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               END IF
               IF(NDOTQOLD < 0.0) THEN
                  INCOMEOLD=0.8
                  !INCOMEOLD=1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               ELSE
                  INCOMEOLD=0.2
                  !INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               END IF

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
                       + SCVFEN(CV_KLOC,GI) * FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
                  FEMTOLDGI_IPHA = FEMTOLDGI_IPHA &
                       + SCVFEN(CV_KLOC,GI) * FEMTOLD(CV_KNOD+(IPHASE-1)*CV_NONODS)
               END DO
               !FEMTGI_IPHA = ( MASS_CV(CV_NODJ) * T(CV_NODI_IPHA) + &
               !     MASS_CV(CV_NODI) * T(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               !FEMTOLDGI_IPHA = ( MASS_CV(CV_NODJ) * TOLD(CV_NODI_IPHA) + &
               !     MASS_CV(CV_NODI) * TOLD(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))

               !GEOMTGI_IPHA = ( MASS_CV(CV_NODJ) * T(CV_NODI_IPHA) + &
               !     MASS_CV(CV_NODI) * T(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               !GEOMTOLDGI_IPHA = ( MASS_CV(CV_NODJ) * TOLD(CV_NODI_IPHA) + &
               !     MASS_CV(CV_NODI) * TOLD(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))

               ! Central is fine as its within an element with equally spaced nodes. 
               FEMTGI_IPHA = 0.5 * ( T(CV_NODI_IPHA) + T(CV_NODJ_IPHA) )
               FEMTOLDGI_IPHA = 0.5 * ( TOLD(CV_NODI_IPHA) + TOLD(CV_NODJ_IPHA) )
               GEOMTGI_IPHA = 0.5 * ( T(CV_NODI_IPHA) + T(CV_NODJ_IPHA) )
               GEOMTOLDGI_IPHA = 0.5*(TOLD(CV_NODI_IPHA) + TOLD(CV_NODJ_IPHA) )
               !GEOMTGI_IPHA = FEMTGI_IPHA
               !GEOMTOLDGI_IPHA = FEMTOLDGI_IPHA

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
               OVER_RELAX=1.0
               MAX_OPER=.TRUE.
               CONSERV=.TRUE.
               SAT_BASED=.FALSE.
               IF(.false.) THEN
                  !FEMTOLDGI_IPHA=(LIMT-FEMTGI) + GEOMTOLDGI_IPHA
                  FEMTOLDGI_IPHA=LIMT
                  !FEMTOLDGI_IPHA=sqrt(LIMT-FEMTOLDGI_IPHA) + FEMTOLDGI_IPHA
                  !GEOMTOLDGI_IPHA=FEMTGI

                  OVER_RELAX=1.0
                  MAX_OPER=.FALSE.
                  SAT_BASED=.TRUE.
               END IF

               CALL FIND_OPT_INCOME_INTERP(INCOME, W_UPWIND, NDOTQ, NDOTQ2, IPHASE, &
                    T(CV_NODI_IPHA),T(CV_NODJ_IPHA), FEMTGI_IPHA, GEOMTGI_IPHA,OVER_RELAX, &
                    ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                    GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

               CALL FIND_OPT_INCOME_INTERP(INCOMEOLD, W_UPWIND, NDOTQOLD, NDOTQOLD2, IPHASE, &
                    TOLD(CV_NODI_IPHA),TOLD(CV_NODJ_IPHA), FEMTOLDGI_IPHA, GEOMTOLDGI_IPHA,OVER_RELAX, &
                    ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                    GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

               if (.true.) then ! Chris look at this...

                  IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                     INCOME3=1.0
                  ELSE
                     INCOME3=0.0
                  END IF
                  IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                     INCOMEOLD3=1.0
                  ELSE
                     INCOMEOLD3=0.0
                  END IF

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
                  END IF


                  CALL ONVDLIMsqrt( CV_NONODS, &
                       LIMT3, FEMTGI_IPHA, INCOME3, CV_NODI, CV_NODJ, &
                       T( (IPHASE-1)*CV_NONODS+1 ), &
                       TMIN( (IPHASE-1)*CV_NONODS+1 ), TMAX( (IPHASE-1)*CV_NONODS+1 ), .FALSE. , .FALSE.)


                  CALL ONVDLIMsqrt( CV_NONODS, &
                       LIMTOLD3, FEMTOLDGI_IPHA, INCOMEOLD3, CV_NODI, CV_NODJ, &
                       TOLD( (IPHASE-1)*CV_NONODS+1 ), &
                       TOLDMIN( (IPHASE-1)*CV_NONODS+1 ), TOLDMAX( (IPHASE-1)*CV_NONODS+1 ), .FALSE., .FALSE. )

                  INCOME4 = ( LIMT3 - T( CV_NODI_IPHA ) ) &
                       / TOLFUN( T( CV_NODJ_IPHA ) - T( CV_NODI_IPHA ) )

                  INCOMEOLD4 = ( LIMTOLD3 - TOLD( CV_NODI_IPHA ) ) &
                       / TOLFUN( TOLD( CV_NODJ_IPHA ) - TOLD( CV_NODI_IPHA ) )

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
                  END IF

               end if

            ELSE

               UPWIND_FRAC=0.5
               IF(NDOTQ < 0.0) THEN
                  !INCOME=0.8
                  INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               ELSE
                  !INCOME=0.2
                  INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               END IF
               IF(NDOTQOLD < 0.0) THEN
                  !INCOMEOLD=0.8
                  INCOMEOLD=1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               ELSE
                  !INCOMEOLD=0.2
                  INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
               END IF
               !INCOME=0.5
               !INCOMEOLD=0.5
            END IF

            UDGI = 0.0
            VDGI = 0.0
            WDGI = 0.0
            UOLDDGI = 0.0
            VOLDDGI = 0.0
            WOLDDGI = 0.0
            UGI_COEF_ELE=0.0
            VGI_COEF_ELE=0.0
            WGI_COEF_ELE=0.0 
            DO U_KLOC_LEV = 1, U_NLOC_LEV
               U_KLOC =(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
               U_KLOC2=(CV_JLOC-1)*U_NLOC_LEV + U_KLOC_LEV
               U_NODK_IPHA  = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
               U_NODK2_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC2 )+(IPHASE-1)*U_NONODS

               UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA ) * (1.0-INCOME) &
                    + SUFEN( U_KLOC2, GI ) * NU( U_NODK2_IPHA ) * INCOME 
               VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA ) * (1.0-INCOME) &
                    + SUFEN( U_KLOC2, GI ) * NV( U_NODK2_IPHA ) * INCOME
               WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA ) * (1.0-INCOME) &
                    + SUFEN( U_KLOC2, GI ) * NW( U_NODK2_IPHA ) * INCOME

               UOLDDGI = UOLDDGI + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK_IPHA ) * (1.0-INCOMEOLD) &
                    + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK2_IPHA ) * INCOMEOLD
               VOLDDGI = VOLDDGI + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK_IPHA ) * (1.0-INCOMEOLD) &
                    + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK2_IPHA ) * INCOMEOLD
               WOLDDGI = WOLDDGI + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK_IPHA ) * (1.0-INCOMEOLD) &
                    + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK2_IPHA ) * INCOMEOLD

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
            UGI_COEF_ELE = 0.0
            VGI_COEF_ELE = 0.0
            WGI_COEF_ELE = 0.0
            DO U_KLOC_LEV = 1, U_NLOC_LEV
               U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
               U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) + (IPHASE-1)*U_NONODS

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
            UGI_COEF_ELE2 = 0.0
            VGI_COEF_ELE2 = 0.0
            WGI_COEF_ELE2 = 0.0
            DO U_KLOC_LEV = 1, U_NLOC_LEV
               U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
               U_KLOC2 = U_OTHER_LOC( U_KLOC )
               IF( U_KLOC2 /= 0 ) THEN
                  U_NODK2_IPHA = U_NDGLN((ELE2-1)*U_NLOC+U_KLOC2)+(IPHASE-1)*U_NONODS

                  UDGI2  = UDGI2 + SUFEN( U_KLOC, GI ) * NU( U_NODK2_IPHA )
                  VDGI2  = VDGI2 + SUFEN( U_KLOC, GI ) * NV( U_NODK2_IPHA )
                  WDGI2  = WDGI2 + SUFEN( U_KLOC, GI ) * NW( U_NODK2_IPHA )

                  UGI_COEF_ELE2(U_KLOC2)=UGI_COEF_ELE2(U_KLOC2)+1.0
                  VGI_COEF_ELE2(U_KLOC2)=VGI_COEF_ELE2(U_KLOC2)+1.0
                  WGI_COEF_ELE2(U_KLOC2)=WGI_COEF_ELE2(U_KLOC2)+1.0

                  UOLDDGI2  = UOLDDGI2 + SUFEN( U_KLOC, GI ) * NUOLD( U_NODK2_IPHA )
                  VOLDDGI2  = VOLDDGI2 + SUFEN( U_KLOC, GI ) * NVOLD( U_NODK2_IPHA )
                  WOLDDGI2  = WOLDDGI2 + SUFEN( U_KLOC, GI ) * NWOLD( U_NODK2_IPHA )
               END IF
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
                  END IF
                  IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                     INCOMEOLD=1.0
                  ELSE
                     INCOMEOLD=0.0
                  END IF

               ELSE IF(DG_ELE_UPWIND==2) THEN ! the best

                  UPWIND_FRAC=0.8
                  IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                     !INCOME=0.8
                     INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  ELSE
                     !INCOME=0.2
                     INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  END IF
                  IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                     !INCOMEOLD=0.8
                     INCOMEOLD=1.0 -2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  ELSE
                     !INCOMEOLD=0.2
                     INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  END IF

               ELSE IF(DG_ELE_UPWIND==3) THEN ! the best optimal upwind frac.

                  MAT_NODI=MAT_NDGLN((ELE-1)*MAT_NLOC+CV_ILOC)
                  MAT_NODJ=MAT_NDGLN((ELE2-1)*MAT_NLOC+CV_JLOC)

                  FEMTGI_IPHA = 0.0
                  FEMTOLDGI_IPHA = 0.0
                  DO CV_KLOC = 1, CV_NLOC
                     CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                     IF(CV_KLOC2 /= 0 ) THEN
                        if (.false.) then
                           CV_KNOD = CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                           CV_KNOD2 = CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
                           FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) & 
                                * 0.5 * ( FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS) &
                                + FEMT(CV_KNOD2+(IPHASE-1)*CV_NONODS) )
                           FEMTOLDGI_IPHA = FEMTOLDGI_IPHA + SCVFEN(CV_KLOC,GI) & 
                                * 0.5 * ( FEMTOLD(CV_KNOD+(IPHASE-1)*CV_NONODS) &
                                + FEMTOLD(CV_KNOD2+(IPHASE-1)*CV_NONODS) )
                        end if

                        if (.true.) then
                           if (0.5*(NDOTQ+NDOTQ2)>0.0) then
                              CV_KNOD = CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                              FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                   * FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
                           else
                              CV_KNOD2 = CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
                              FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                   * FEMT(CV_KNOD2+(IPHASE-1)*CV_NONODS)
                           end if
                           if (0.5*(NDOTQOLD+NDOTQOLD2)>0.0) then
                              CV_KNOD = CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                              FEMTOLDGI_IPHA = FEMTOLDGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                   * FEMTOLD(CV_KNOD+(IPHASE-1)*CV_NONODS)
                           else
                              CV_KNOD2 = CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
                              FEMTOLDGI_IPHA = FEMTOLDGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                   * FEMTOLD(CV_KNOD2+(IPHASE-1)*CV_NONODS)
                           end if
                        end if

                        if (.false.) then
                           CV_KNOD = CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                           FEMTGI_IPHA = FEMTGI_IPHA &
                                + SCVFEN(CV_KLOC,GI) * FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
                           FEMTOLDGI_IPHA = FEMTOLDGI_IPHA &
                                + SCVFEN(CV_KLOC,GI) * FEMTOLD(CV_KNOD+(IPHASE-1)*CV_NONODS)
                        end if
                     END IF
                  END DO
                  !FEMTOLDGI_IPHA = (MASS_CV(CV_NODJ)* T(CV_NODI_IPHA) + &
                  !   MASS_CV(CV_NODI)*T(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  !FEMTOLDGI_IPHA = (MASS_CV(CV_NODi)* T(CV_NODI_IPHA) + &
                  !   MASS_CV(CV_NODj)*T(CV_NODJ_IPHA) )/ (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  !FEMTOLDGI_IPHA = 0.5*(T(CV_NODI_IPHA) + T(CV_NODJ_IPHA) )
                  GEOMTGI_IPHA = ( MASS_CV(CV_NODJ) * T(CV_NODI_IPHA) + &
                       MASS_CV(CV_NODI) * T(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  GEOMTOLDGI_IPHA = ( MASS_CV(CV_NODJ) * TOLD(CV_NODI_IPHA) + &
                       MASS_CV(CV_NODI) * TOLD(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  !GEOMTGI_IPHA =FEMTGI_IPHA
                  !GEOMTOLDGI_IPHA =FEMTOLDGI_IPHA

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
                  CONSERV=.TRUE.
                  SAT_BASED=.FALSE.

                  CALL FIND_OPT_INCOME_INTERP(INCOME, W_UPWIND, NDOTQ, NDOTQ2, IPHASE, &
                       T(CV_NODI_IPHA),T(CV_NODJ_IPHA), FEMTGI_IPHA, GEOMTGI_IPHA,OVER_RELAX, &
                       ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                       GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

                  CALL FIND_OPT_INCOME_INTERP(INCOMEOLD, W_UPWIND, NDOTQOLD, NDOTQOLD2, IPHASE,  &
                       TOLD(CV_NODI_IPHA),TOLD(CV_NODJ_IPHA), FEMTOLDGI_IPHA, GEOMTOLDGI_IPHA,OVER_RELAX, &
                       ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                       GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA,CONSERV,MAX_OPER,SAT_BASED)

                  if(.true.) then

                     IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN
                        INCOME3=1.0
                     ELSE
                        INCOME3=0.0
                     END IF
                     IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                        INCOMEOLD3=1.0
                     ELSE
                        INCOMEOLD3=0.0
                     END IF

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
                     END IF


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
                     END IF
 
                  end if

                  OVER_RELAX=1.
                  MAX_OPER=.TRUE.
                  CONSERV=.TRUE.
                  SAT_BASED=.FALSE.
                  IF(.true.) THEN
                     FEMTGI_IPHA=LIMT3
                     FEMTOLDGI_IPHA=LIMTOLD3
                     OVER_RELAX=1.0
                     MAX_OPER=.FALSE.
                     SAT_BASED=.TRUE.
                  END IF

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
                     !INCOME=0.8
                     INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  ELSE
                     !INCOME=0.2
                     INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  END IF
                  IF(0.5*(NDOTQOLD+NDOTQOLD2) < 0.0) THEN
                     !INCOMEOLD=0.8
                     INCOMEOLD=1.0 -2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  ELSE
                     !INCOMEOLD=0.2
                     INCOMEOLD=2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                  END IF
               END IF

               DT_I = 1.0 - INCOME
               DT_J = 1.0 - DT_I

               DTOLD_I = 1.0 - INCOMEOLD
               DTOLD_J = 1.0 - DTOLD_I
            END IF

            ! Amend weighting for porosity only across elements...
            IF((ABS(CV_DG_VEL_INT_OPT ) == 2).OR.(ABS(CV_DG_VEL_INT_OPT ) == 3)) THEN 
               IF(ELE /= ELE2) THEN 
                  DT_I=VOLFRA_PORE(ELE) *DT_I
                  DT_J=VOLFRA_PORE(ELE2)*DT_J
                  DTOLD_I=VOLFRA_PORE(ELE) *DTOLD_I
                  DTOLD_J=VOLFRA_PORE(ELE2)*DTOLD_J
               END IF
            END IF

            !dt_i=0.5
            !dt_j=0.5
            !dtold_i=0.5
            !dtold_j=0.5

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

               IF( NDOTQ_INT <= 0.0 ) THEN ! Incoming
                  !DT_I=1.0
                  DT_J=DT_I+DT_J
               ELSE
                  DT_I=DT_I+DT_J
                  !DT_J=1.0
               END IF
               IF( NDOTQOLD_INT <= 0.0 ) THEN ! Incoming
                  !DTOLD_I=1.0
                  DTOLD_J=DTOLD_I+DTOLD_J
               ELSE
                  DTOLD_I=DTOLD_I+DTOLD_J
                  !DTOLD_J=1.0
               END IF

               UDGI_INT = ( DT_I * UDGI + DT_J * UDGI2 ) / (DT_I + DT_J)
               VDGI_INT = ( DT_I * VDGI + DT_J * VDGI2 ) / (DT_I + DT_J)
               WDGI_INT = ( DT_I * WDGI + DT_J * WDGI2 ) / (DT_I + DT_J)

               UOLDDGI_INT = ( DTOLD_I * UOLDDGI + DTOLD_J * UOLDDGI2 ) / (DTOLD_I + DTOLD_J)
               VOLDDGI_INT = ( DTOLD_I * VOLDDGI + DTOLD_J * VOLDDGI2 ) / (DTOLD_I + DTOLD_J)
               WOLDDGI_INT = ( DTOLD_I * WOLDDGI + DTOLD_J * WOLDDGI2 ) / (DTOLD_I + DTOLD_J)
            END IF

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

         END IF Conditional_ELE2

      END IF Conditional_SELE

      NDOTQ =  CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI + CVNORMZ(GI) * WDGI
      NDOTQOLD =  CVNORMX( GI ) * UOLDDGI + CVNORMY( GI ) * VOLDDGI + CVNORMZ(GI) * WOLDDGI

      ! Define whether flux is incoming or outgoing, depending on direction of flow
      INCOME = 1.
      IF( NDOTQ >= 0. ) INCOME = 0.

      INCOMEOLD = 1.
      IF( NDOTQOLD >= 0. ) INCOMEOLD = 0.

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
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC,  LIMIT_USE_2ND)

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
      logical, intent( in ) :: LIMIT_USE_2ND
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC

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
         T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2, &
         TMIN_2ND_MC, TOLDMIN_2ND_MC, T2MIN_2ND_MC, T2OLDMIN_2ND_MC, DENMIN_2ND_MC, DENOLDMIN_2ND_MC, &
         TMAX_2ND_MC, TOLDMAX_2ND_MC, T2MAX_2ND_MC, T2OLDMAX_2ND_MC, DENMAX_2ND_MC, DENOLDMAX_2ND_MC, LIMIT_USE_2ND, &
         HDC, NDOTQ, NDOTQOLD, DT, &
         SCVFENX, SCVFENY, SCVFENZ, CVNORMX, CVNORMY, CVNORMZ, &
         U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
         IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
         TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT )
      !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
      IMPLICIT NONE
      ! Calculate T and DEN on the CV face at quadrature point GI.
      REAL, intent( inout ) :: FVT,FVTOLD, FVT2, FVT2OLD,FVD,FVDOLD, LIMD,LIMT,LIMT2, &
           LIMDOLD,LIMTOLD,LIMT2OLD,LIMDT,LIMDTOLD,LIMDTT2,LIMDTT2OLD, &
           FEMDGI, FEMTGI, FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI
      REAL, intent( in ) :: INCOME,INCOMEOLD,HDC,NDOTQ,NDOTQOLD,DT
      logical, intent( in ) :: LIMIT_USE_2ND
      INTEGER, intent( in ) :: CV_DISOPT,CV_NONODS,NPHASE,CV_NODI_IPHA,CV_NODJ_IPHA,ELE,ELE2,  &
           CV_NLOC,TOTELE,SCVNGI,GI,IPHASE,SELE,CV_SNLOC,STOTEL, &
           WIC_T_BC_DIRICHLET,WIC_D_BC_DIRICHLET, IGOT_T2, U_NLOC,U_NONODS,NDIM
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( CV_NLOC ), intent( in ) :: CV_OTHER_LOC
      INTEGER, DIMENSION( CV_SNLOC ), intent( in ) :: CV_SLOC2LOC
      INTEGER, DIMENSION( STOTEL * NPHASE ), intent( in ) :: WIC_T_BC, WIC_D_BC
      INTEGER, DIMENSION( STOTEL * NPHASE * IGOT_T2 ), intent( in ) :: WIC_T2_BC
      REAL, DIMENSION( U_NLOC, SCVNGI  ), intent( in ) :: SUFEN
      REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFEN
      REAL, DIMENSION( CV_NLOC, SCVNGI  ), intent( in ) :: SCVFENX, SCVFENY, SCVFENZ
      REAL, DIMENSION( SCVNGI  ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
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
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC,  &
           DENMAX_2ND_MC, DENMIN_2ND_MC, DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, T2OLDMIN_2ND_MC
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
      REAL, DIMENSION( NDIM, NDIM, SCVNGI ), intent( in ) :: INV_JAC

      INTEGER, intent( in ) :: IANISOTROPIC
      INTEGER, intent( in ) :: NSMALL_COLM
      INTEGER, DIMENSION( (TOTELE+1)*IANISOTROPIC ), intent( in ) :: SMALL_FINDRM
      INTEGER, DIMENSION( NSMALL_COLM*IANISOTROPIC ), intent( in ) :: SMALL_COLM
      REAL, DIMENSION( NSMALL_COLM*NPHASE*IANISOTROPIC), intent( in ) :: TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT
      REAL, DIMENSION( NSMALL_COLM*NPHASE*IANISOTROPIC*IGOT_T2), intent( in ) :: T2UPWIND_MAT, T2OLDUPWIND_MAT
      ! Local variables
      ! If UPWIND then use upwind flux between elements else use central. 
      ! If HI_ORDER_HALF then use high order interpolation when around 
      ! a volume frac of 0.5 and gradually apply limiting near 0 and 1. 
      LOGICAL, PARAMETER :: UPWIND = .TRUE., HI_ORDER_HALF = .FALSE., LIM_VOL_ADJUST2 = .TRUE.
      LOGICAL, PARAMETER :: DOWNWIND_EXTRAP = .TRUE. ! Extrapolate a downwind value for interface tracking.
      
      ! Scaling to reduce the downwind bias(=1downwind, =0central)
      LOGICAL, PARAMETER :: SCALE_DOWN_WIND = .true.

      LOGICAL :: FIRSTORD, NOLIMI, RESET_STORE, LIM_VOL_ADJUST
      REAL :: RELAX, RELAXOLD, TMIN_STORE, TMAX_STORE, TOLDMIN_STORE, TOLDMAX_STORE, &
           T2MIN_STORE, T2MAX_STORE, T2OLDMIN_STORE, T2OLDMAX_STORE, &
           DENMIN_STORE, DENMAX_STORE, DENOLDMIN_STORE, DENOLDMAX_STORE, &
           COURANT_OR_MINUS_ONE_NEW, COURANT_OR_MINUS_ONE_OLD
      INTEGER :: CV_KLOC, CV_NODK, CV_NODK_IPHA, CV_KLOC2, CV_NODK2, CV_NODK2_IPHA, CV_STAR_IPHA, &
           CV_SKLOC, CV_SNODK, CV_SNODK_IPHA, U_KLOC,U_NODK,U_NODK_IPHA, IDIM

      LOGICAL :: VOF_METHOD, VOF_INTER, VOF_INTER_OLD
      REAL :: T_BETWEEN_MIN, T_BETWEEN_MAX, TOLD_BETWEEN_MIN, TOLD_BETWEEN_MAX
      REAL :: T_AVE_EDGE, T_AVE_ELE, TOLD_AVE_EDGE, TOLD_AVE_ELE
      REAL :: T_MIDVAL, TOLD_MIDVAL
      REAL :: T_UPWIND,TOLD_UPWIND,TMIN_UPWIND,TMAX_UPWIND,TOLDMIN_UPWIND,TOLDMAX_UPWIND
      REAL :: RSHAPE,RSHAPE_OLD,RGRAY, UDGI,VDGI,WDGI, RSCALE, TXGI,TYGI,TZGI
      REAL :: VEC_VEL(3), VEC_VEL2(3), ELE_LENGTH_SCALE

      ! The adjustment method is not ready for the LIMIT_USE_2ND - the new limiting method.
      LIM_VOL_ADJUST =LIM_VOL_ADJUST2.AND.(.NOT.LIMIT_USE_2ND)

      FVT2    = 1.0
      FVT2OLD = 1.0

      if ( cv_disopt>=8 ) then
         courant_or_minus_one_new = abs ( dt * ndotq / hdc )
         courant_or_minus_one_old = abs ( dt * ndotqold / hdc )
      else
         courant_or_minus_one_new = -1.0
         courant_or_minus_one_old = -1.0
      end if


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

            RSCALE=1.0 ! Scaling to reduce the downwind bias(=1downwind, =0central)
            IF(SCALE_DOWN_WIND) THEN
               IF(DOWNWIND_EXTRAP.AND.(courant_or_minus_one_new.GE.0.0)) THEN
                  TXGI=0.0
                  TYGI=0.0
                  TZGI=0.0
                  DO CV_KLOC = 1, CV_NLOC
                    CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                    CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                    TXGI=TXGI+SCVFENX( CV_KLOC, GI )*FEMT(CV_NODK_IPHA)
      IF(NDIM.GE.2) TYGI=TYGI+SCVFENY( CV_KLOC, GI )*FEMT(CV_NODK_IPHA)
      IF(NDIM.GE.3) TZGI=TZGI+SCVFENZ( CV_KLOC, GI )*FEMT(CV_NODK_IPHA)
                  END DO

                  UDGI=0.0
                  VDGI=0.0
                  WDGI=0.0
                  DO U_KLOC = 1, U_NLOC
                      U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
                      U_NODK_IPHA = U_NODK + ( IPHASE - 1 ) * U_NONODS
                      UDGI=UDGI+SUFEN( U_KLOC, GI )*U(U_NODK_IPHA)
        IF(NDIM.GE.2) VDGI=VDGI+SUFEN( U_KLOC, GI )*V(U_NODK_IPHA)
        IF(NDIM.GE.3) WDGI=WDGI+SUFEN( U_KLOC, GI )*W(U_NODK_IPHA)
                  END DO
 
!                  RSCALE=ABS(CVNORMX(GI)*UDGI+CVNORMY(GI)*VDGI+CVNORMZ(GI)*WDGI) &
!                        /TOLFUN(UDGI**2+VDGI**2+WDGI**2)
! cosine rule:
                  RSCALE=ABS(TXGI*UDGI+TYGI*VDGI+TZGI*WDGI) &
                        /TOLFUN((UDGI**2+VDGI**2+WDGI**2)*SQRT(TXGI**2+TYGI**2+TZGI**2))
! no cosine rule:
!                  RSCALE=1.0 &
!                        /TOLFUN(sqrt(UDGI**2+VDGI**2+WDGI**2))

                  VEC_VEL(1)=UDGI
                  VEC_VEL(2)=VDGI
                  VEC_VEL(3)=WDGI
                  VEC_VEL2=0.0
                  DO IDIM=1,NDIM
                     VEC_VEL2(IDIM)=SUM( INV_JAC(IDIM, 1:NDIM, GI)*VEC_VEL(1:NDIM) )
                  END DO
! normalize the velocity in here: 
!                  VEC_VEL2=VEC_VEL2/TOLFUN(SQRT( UDGI**2+VDGI**2+WDGI**2))

                  ELE_LENGTH_SCALE=0.5*SQRT( (UDGI**2+VDGI**2+WDGI**2)/TOLFUN( SUM( VEC_VEL2(1:NDIM)**2 ))  )
!                  ELE_LENGTH_SCALE=1.0/TOLFUN( SQRT(SUM( VEC_VEL2(1:NDIM)**2 )) )  
!                  ELE_LENGTH_SCALE=0.5*HDC
! For discontinuous elements half the length scale...
                  IF(U_NONODS==CV_NONODS) ELE_LENGTH_SCALE=0.5*ELE_LENGTH_SCALE

               ENDIF
            ENDIF
!          print *,'ele_length_scale=',ele_length_scale
            DO CV_KLOC = 1, CV_NLOC
               CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
               CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
               IF(DOWNWIND_EXTRAP.AND.(courant_or_minus_one_new.GE.0.0)) THEN ! Extrapolate to the downwind value...
!                  RGRAY=RSCALE*ELE_LENGTH_SCALE*( CVNORMX(GI)*SCVFENX( CV_KLOC, GI ) &
!                       + CVNORMY(GI)*SCVFENY( CV_KLOC, GI )+CVNORMZ(GI)*SCVFENZ( CV_KLOC, GI ) )
                  RGRAY=RSCALE*ELE_LENGTH_SCALE*( UDGI*SCVFENX( CV_KLOC, GI ) &
                       + VDGI*SCVFENY( CV_KLOC, GI )+WDGI*SCVFENZ( CV_KLOC, GI ) )
                  RSHAPE    =SCVFEN( CV_KLOC, GI ) + RGRAY
                  RSHAPE_OLD=SCVFEN( CV_KLOC, GI ) + RGRAY
                  FEMTGI    = FEMTGI     +  RSHAPE     * FEMT( CV_NODK_IPHA )
                  FEMTOLDGI = FEMTOLDGI  +  RSHAPE_OLD * FEMTOLD( CV_NODK_IPHA )
               ELSE
                  FEMTGI    = FEMTGI     +  SCVFEN( CV_KLOC, GI ) * FEMT( CV_NODK_IPHA )
                  FEMTOLDGI = FEMTOLDGI  +  SCVFEN( CV_KLOC, GI ) * FEMTOLD( CV_NODK_IPHA )
               ENDIF
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

      VOF_INTER = .TRUE.
      VOF_INTER_OLD = .TRUE.

      ! use vof type approach
      VOF_METHOD= .false.

      IF( VOF_METHOD ) THEN 
         ! Find interface value pts and see if 0.5 lies within them - then  VOF_INTER=.TRUE.
         ! if VOF_INTER then use normal limiting options 
         ! else assume 1.0 or 0.0 depending on the closest value. 
         ! VOF_INTER_OLD is for the old variable

         T_AVE_EDGE = 0.5 * ( T( CV_NODI_IPHA ) + T( CV_NODJ_IPHA ) )
         TOLD_AVE_EDGE = 0.5 * ( TOLD( CV_NODI_IPHA ) + TOLD( CV_NODJ_IPHA ) )

         T_AVE_ELE = 0.
         TOLD_AVE_ELE = 0. 
         DO CV_KLOC = 1, CV_NLOC
            CV_NODK_IPHA = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC ) + ( IPHASE - 1 ) * CV_NONODS
            T_AVE_ELE = T_AVE_ELE + T( CV_NODK_IPHA ) / CV_NLOC
            TOLD_AVE_ELE = TOLD_AVE_ELE + TOLD( CV_NODK_IPHA ) / CV_NLOC
         END DO

         T_UPWIND=T( CV_NODI_IPHA ) * (1. - INCOME) + T( CV_NODJ_IPHA ) * INCOME
         TOLD_UPWIND=TOLD( CV_NODI_IPHA ) * (1. - INCOMEOLD) + TOLD( CV_NODJ_IPHA ) * INCOMEOLD
         TMIN_UPWIND=TMIN( CV_NODI_IPHA ) * (1. - INCOME) + TMIN( CV_NODJ_IPHA ) * INCOME
         TMAX_UPWIND=TMAX( CV_NODI_IPHA ) * (1. - INCOME) + TMAX( CV_NODJ_IPHA ) * INCOME 
         TOLDMIN_UPWIND=TOLDMIN( CV_NODI_IPHA ) * (1. - INCOMEOLD) + TOLDMIN( CV_NODJ_IPHA ) * INCOMEOLD
         TOLDMAX_UPWIND=TOLDMAX( CV_NODI_IPHA ) * (1. - INCOMEOLD) + TOLDMAX( CV_NODJ_IPHA ) * INCOMEOLD

         IF( ( T_UPWIND >= TMAX_UPWIND-1.E-6 ) .OR. ( T_UPWIND <= TMIN_UPWIND+1.E-6 ) ) THEN
           ! Extrema...
           T_MIDVAL = 0.5 * ( TMIN_UPWIND + TMAX_UPWIND )
           !IF( T_UPWIND >= 0.95 ) T_MIDVAL = 0.5  
         ELSE
           T_MIDVAL = 0.5 
         ENDIF

         IF( ( TOLD_UPWIND >= TOLDMAX_UPWIND-1.E-6 ) .OR. ( TOLD_UPWIND <= TOLDMIN_UPWIND+1.E-6 ) ) THEN
           ! Extrema...
           TOLD_MIDVAL = 0.5 * ( TOLDMIN_UPWIND + TOLDMAX_UPWIND )
           !IF( TOLD_UPWIND >= 0.95 ) TOLD_MIDVAL = 0.5  
         ELSE
           TOLD_MIDVAL = 0.5
         ENDIF


         ! hard wired for interfaces...
         T_MIDVAL = 0.5 
         TOLD_MIDVAL = 0.5


         T_BETWEEN_MIN = MIN( T_AVE_EDGE, T_AVE_ELE )
         T_BETWEEN_MAX = MAX( T_AVE_EDGE, T_AVE_ELE )

         TOLD_BETWEEN_MIN = MIN( TOLD_AVE_EDGE, TOLD_AVE_ELE )
         TOLD_BETWEEN_MAX = MAX( TOLD_AVE_EDGE, TOLD_AVE_ELE )

         femtgi = ( t_between_min + t_between_max ) / 2 
         femtoldgi = ( told_between_min + told_between_max ) / 2 

         VOF_INTER = .FALSE.
         VOF_INTER_OLD = .FALSE.
         IF( T_BETWEEN_MIN < T_MIDVAL .AND. T_BETWEEN_MAX > T_MIDVAL ) THEN
            VOF_INTER = .TRUE.
            femtgi = (t_between_min+t_between_max)/2 !max(min(femtgi, t_between_max), t_between_min)
         END IF

         IF( TOLD_BETWEEN_MIN < TOLD_MIDVAL .AND. TOLD_BETWEEN_MAX > TOLD_MIDVAL ) then
            VOF_INTER_OLD = .TRUE. 
            femtoldgi = (told_between_min+told_between_max)/2 !max(min(femtoldgi, told_between_max), told_between_min)
         END IF


!if( ( T(CV_NODI_IPHA) >= TMAX(CV_NODI_IPHA)-1.E-4 .or. T(CV_NODJ_IPHA) >= TMAX(CV_NODJ_IPHA)-1.E-4 ) .and. (T(CV_NODI_IPHA) <= 0.5+0.1 .OR. T(CV_NODJ_IPHA)<= 0.5+0.1) ) VOF_INTER=.true.
!if( ( T(CV_NODI_IPHA) <= TMIN(CV_NODI_IPHA)+1.E-4 .or. T(CV_NODJ_IPHA) <= TMIN(CV_NODJ_IPHA)+1.E-4 ) .and. (T(CV_NODI_IPHA) >= 0.5-0.1 .OR. T(CV_NODJ_IPHA)>= 0.5-0.1) ) VOF_INTER=.true.

!if( ( Told(CV_NODI_IPHA) >= ToldMAX(CV_NODI_IPHA)-1.E-4 .or. Told(CV_NODJ_IPHA) >= ToldMAX(CV_NODJ_IPHA)-1.E-4 ) .and. (Told(CV_NODI_IPHA) <= 0.5+0.1 .OR. Told(CV_NODJ_IPHA)<= 0.5+0.1) ) VOF_INTER_OLD=.true.
!if( ( Told(CV_NODI_IPHA) <= ToldMIN(CV_NODI_IPHA)+1.E-4 .or. Told(CV_NODJ_IPHA) <= ToldMIN(CV_NODJ_IPHA)+1.E-4 ) .and. (Told(CV_NODI_IPHA) >= 0.5-0.1 .OR. Told(CV_NODJ_IPHA)>= 0.5-0.1) ) VOF_INTER_OLD=.true.

!if( ( T(CV_NODI_IPHA) >= TMAX(CV_NODI_IPHA)-1.E-4 .or. T(CV_NODJ_IPHA) >= TMAX(CV_NODJ_IPHA)-1.E-4 ) .and. T_UPWIND <= 0.5+0.1 ) VOF_INTER=.true.
!if( ( T(CV_NODI_IPHA) <= TMIN(CV_NODI_IPHA)+1.E-4 .or. T(CV_NODJ_IPHA) <= TMIN(CV_NODJ_IPHA)+1.E-4 ) .and. T_UPWIND >= 0.5-0.1 ) VOF_INTER=.true.

!if( ( Told(CV_NODI_IPHA) >= ToldMAX(CV_NODI_IPHA)-1.E-4 .or. Told(CV_NODJ_IPHA) >= ToldMAX(CV_NODJ_IPHA)-1.E-4 ) .and. Told_UPWIND <= 0.5+0.1 ) VOF_INTER_OLD=.true.
!if( ( Told(CV_NODI_IPHA) <= ToldMIN(CV_NODI_IPHA)+1.E-4 .or. Told(CV_NODJ_IPHA) <= ToldMIN(CV_NODJ_IPHA)+1.E-4 ) .and. Told_UPWIND >= 0.5-0.1 ) VOF_INTER_OLD=.true.

!if( T_UPWIND >= TMAX_UPWIND-1.E-4 .and. T_UPWIND <= 0.5+0.1 ) VOF_INTER=.true.
!if( T_UPWIND <= TMIN_UPWIND+1.E-4 .and. T_UPWIND >= 0.5-0.1 ) VOF_INTER=.true.
!if( T_UPWIND >= 0.5+0.3 .or. T_UPWIND <= 0.5-0.3 ) VOF_INTER=.true.

!if( TOLD_UPWIND >= TOLDMAX_UPWIND-1.E-4 .and. TOLD_UPWIND <= 0.5+0.1 ) VOF_INTER_OLD=.true.
!if( TOLD_UPWIND <= TOLDMIN_UPWIND+1.E-4 .and. TOLD_UPWIND >= 0.5-0.1 ) VOF_INTER_OLD=.true.
!if( TOLD_UPWIND >= 0.5+0.3 .or. TOLD_UPWIND <= 0.5-0.3 ) VOF_INTER_OLD=.true.

!if( FEMTGI < 0.5 ) VOF_INTER=.true.
!if( FEMToldGI < 0.5 ) VOF_INTER_OLD=.true.




         IF( .NOT.VOF_INTER ) THEN
            IF( FEMTGI > T_MIDVAL ) THEN
               FEMTGI = 1.!MIN(1.0, TMAX_UPWIND )!1.0
            ELSE 
               FEMTGI = 0.!MAX(0.0, TMIN_UPWIND )!0.0
            END IF
         END IF

         IF( .NOT.VOF_INTER_OLD ) THEN
            IF( FEMTOLDGI > TOLD_MIDVAL ) THEN
              FEMTOLDGI = 1.!MIN(1.0, TOLDMAX_UPWIND )!1.0
            ELSE 
              FEMTOLDGI = 0.!MAX(0.0, TOLDMIN_UPWIND )!0.0
            END IF
         END IF

      END IF

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

         CALL ONVDLIM_ALL( CV_NONODS, &
              LIMT, FEMTGI, INCOME, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
              T( CV_STAR_IPHA ), TMIN( CV_STAR_IPHA ), TMAX( CV_STAR_IPHA ), &
              TMIN_2ND_MC( CV_STAR_IPHA ), TMAX_2ND_MC( CV_STAR_IPHA ), FIRSTORD, .not.VOF_INTER, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_NEW, &
              IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, TUPWIND_MAT(1+(IPHASE-1)*NSMALL_COLM*IANISOTROPIC) )

         CALL ONVDLIM_ALL( CV_NONODS, &
              LIMTOLD, FEMTOLDGI, INCOMEOLD, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
              TOLD( CV_STAR_IPHA ), TOLDMIN( CV_STAR_IPHA ), TOLDMAX( CV_STAR_IPHA ), &
              TOLDMIN_2ND_MC( CV_STAR_IPHA ), TOLDMAX_2ND_MC( CV_STAR_IPHA ), FIRSTORD, .not.VOF_INTER_OLD, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_OLD, &
              IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, TOLDUPWIND_MAT(1+(IPHASE-1)*NSMALL_COLM*IANISOTROPIC) )

         CALL ONVDLIM_ALL( CV_NONODS, &
              LIMD, FEMDGI, INCOME, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
              DEN( CV_STAR_IPHA ), DENMIN( CV_STAR_IPHA ), DENMAX( CV_STAR_IPHA ), &
              DENMIN_2ND_MC( CV_STAR_IPHA ), DENMAX_2ND_MC( CV_STAR_IPHA ), FIRSTORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_NEW, &
              IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, DENUPWIND_MAT(1+(IPHASE-1)*NSMALL_COLM*IANISOTROPIC) )

         CALL ONVDLIM_ALL( CV_NONODS, &
              LIMDOLD,FEMDOLDGI, INCOMEOLD, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
              DENOLD( CV_STAR_IPHA ), DENOLDMIN( CV_STAR_IPHA ), DENOLDMAX( CV_STAR_IPHA ), &
              DENOLDMIN_2ND_MC( CV_STAR_IPHA ), DENOLDMAX_2ND_MC( CV_STAR_IPHA ), FIRSTORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_OLD, &
              IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, DENOLDUPWIND_MAT(1+(IPHASE-1)*NSMALL_COLM*IANISOTROPIC) )

         IF(IGOT_T2==1) THEN
            CALL ONVDLIM_ALL( CV_NONODS, &
                 LIMT2, FEMT2GI, INCOME, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
                 T2( CV_STAR_IPHA ), T2MIN( CV_STAR_IPHA ), T2MAX( CV_STAR_IPHA ), &
                 T2MIN_2ND_MC( CV_STAR_IPHA ), T2MAX_2ND_MC( CV_STAR_IPHA ),FIRSTORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_NEW, &
                 IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, T2UPWIND_MAT(1+(IPHASE-1)*NSMALL_COLM*IANISOTROPIC) )

            CALL ONVDLIM_ALL( CV_NONODS, &
                 LIMT2OLD, FEMT2OLDGI, INCOMEOLD, CV_NODI_IPHA-(IPHASE-1)*CV_NONODS, CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS, &
                 T2OLD( CV_STAR_IPHA ), T2OLDMIN( CV_STAR_IPHA ), T2OLDMAX( CV_STAR_IPHA ), &
                 T2OLDMIN_2ND_MC( CV_STAR_IPHA ), T2OLDMAX_2ND_MC( CV_STAR_IPHA ), FIRSTORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_OLD, &
                 IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, T2OLDUPWIND_MAT(1+(IPHASE-1)*NSMALL_COLM*IANISOTROPIC) )
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
      !if(.false.) then
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
      !endif

      IF(HI_ORDER_HALF) THEN
         RELAX=MIN(2.*ABS(FEMTGI-0.5),1.0)
         !relax=0.
         LIMT=RELAX*LIMT+(1.-RELAX)*FEMTGI
         RELAXOLD=MIN(2.*ABS(FEMTOLDGI-0.5),1.0)
         !relaxold=0.
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
      implicit none
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
            !ewrite(3,*)'sele,CV_INOD:',sele,CV_INOD
            NOD_COUNT_SELE(CV_INOD)=NOD_COUNT_SELE(CV_INOD)+1
            NOD_ON_BOUNDARY(CV_INOD)=.TRUE.
         END DO
      END DO
      !ewrite(3,*)'NOD_COUNT_SELE:',NOD_COUNT_SELE
      !ewrite(3,*)'NOD_ON_BOUNDARY:',NOD_ON_BOUNDARY
      !stop 831

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

      !do ele=1,totele
      !   ewrite(3,*)'ele',ele,' FACE_ELE(IFACE,ELE):',(FACE_ELE(IFACE,ELE),iface=1,nface) 
      !end do

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
                        ! ewrite(3,*)'(SELE2-1)*CV_SNLOC+CV_SILOC2,SELE2,CV_SNLOC,CV_SILOC2:', &
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

      !do ele=1,totele
      !   ewrite(3,*)'ele',ele,' FACE_ELE(IFACE,ELE):',(FACE_ELE(IFACE,ELE),iface=1,nface) 
      !end do
      !stop 2982

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
         TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC, &
         T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, T2OLDMIN_2ND_MC, &
         LIMIT_USE_2ND, &
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
      LOGICAL, intent( in ) :: LIMIT_USE_2ND
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

      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC,  &
           DENMAX_2ND_MC, DENMIN_2ND_MC, DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC
      REAL, DIMENSION( CV_NONODS * NPHASE * IGOT_T2 ), intent( inout ) :: T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, T2OLDMIN_2ND_MC

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
                     !ewrite(3,*) 'T2:',TMAX( CV_INOD_IPHA ), T2MIN( CV_INOD_IPHA ), T2OLDMAX( CV_INOD_IPHA ) , T2OLDMIN( CV_INOD_IPHA )
                  ENDIF
               ENDIF
               ! DEN: 
               !ewrite(3,*) 'stotel,nphase=', stotel,nphase
               !ewrite(3,*) 'SELE,SELE+(IPHASE-1)*STOTEL:', SELE,SELE+(IPHASE-1)*STOTEL
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
                  !ewrite(3,*) 'Dens:',DENMAX( CV_INOD_IPHA ), DENMIN( CV_INOD_IPHA ), DENOLDMAX( CV_INOD_IPHA ),DENOLDMIN( CV_INOD_IPHA )
               ENDIF
            END DO

         END DO Loop_CV_SILOC

      END DO Loop_SELE


      IF(LIMIT_USE_2ND) THEN
      Loop2_CV_NODI: DO CV_NODI = 1, CV_NONODS

         DO IPHASE = 1,NPHASE
             CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
             !ewrite(3,*)'CV_NODI, IPHASE, CV_NODI_IPHA=',CV_NODI, IPHASE, CV_NODI_IPHA
             TMAX_2ND_MC(CV_NODI_IPHA)  =-1.E+10
             TMIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
             TOLDMAX_2ND_MC( CV_NODI_IPHA )=-1.E+10
             TOLDMIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
            if(IGOT_T2==1) then
             T2MAX_2ND_MC( CV_NODI_IPHA )=-1.E+10
             T2MIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
             T2OLDMAX_2ND_MC( CV_NODI_IPHA )=-1.E+10
             T2OLDMIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
            endif
             DENMAX_2ND_MC( CV_NODI_IPHA )=-1.E+10
             DENMIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
             DENOLDMAX_2ND_MC( CV_NODI_IPHA )=-1.E+10
             DENOLDMIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
         END DO

         DO COUNT = FINACV( CV_NODI ), FINACV( CV_NODI + 1 ) - 1, 1
            CV_NODJ = COLACV( COUNT )

            IF( CV_NODJ <= CV_NONODS ) THEN

               DO IPHASE = 1,NPHASE
                  CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                  CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS
                IF(CV_NODI_IPHA.NE.CV_NODJ_IPHA) THEN

                  IF(TMAX_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA) THEN
                  IF( T( CV_NODJ_IPHA ) > TMAX_2ND_MC( CV_NODI_IPHA ) ) THEN
                     TMAX_2ND_MC( CV_NODI_IPHA ) = T( CV_NODJ_IPHA )
                  ENDIF
                  ENDIF
                  IF(TMIN_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA) THEN
                  IF( T( CV_NODJ_IPHA ) < TMIN_2ND_MC( CV_NODI_IPHA ) ) THEN
                     TMIN_2ND_MC( CV_NODI_IPHA ) = T( CV_NODJ_IPHA )
                  ENDIF
                  ENDIF

                  IF( TOLDMAX_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                  IF( TOLD( CV_NODJ_IPHA ) > TOLDMAX_2ND_MC( CV_NODI_IPHA ) ) THEN
                     TOLDMAX_2ND_MC( CV_NODI_IPHA ) = TOLD( CV_NODJ_IPHA )
                  ENDIF
                  ENDIF
                  IF( TOLDMIN_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                  IF( TOLD( CV_NODJ_IPHA ) < TOLDMIN_2ND_MC( CV_NODI_IPHA ) ) THEN
                     TOLDMIN_2ND_MC( CV_NODI_IPHA ) = TOLD( CV_NODJ_IPHA )
                  ENDIF
                  ENDIF
                  ! T2: 
                  IF(IGOT_T2==1) THEN
                     IF( T2MAX_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                     IF( T2( CV_NODJ_IPHA ) > T2MAX_2ND_MC( CV_NODI_IPHA ) ) THEN
                        T2MAX_2ND_MC( CV_NODI_IPHA ) = T2( CV_NODJ_IPHA )
                     ENDIF
                     ENDIF
                     IF( T2MIN_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                     IF( T2( CV_NODJ_IPHA ) < T2MIN_2ND_MC( CV_NODI_IPHA ) ) THEN
                        T2MIN_2ND_MC( CV_NODI_IPHA ) = T2( CV_NODJ_IPHA )
                     ENDIF
                     ENDIF

                     IF( T2OLDMAX_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                     IF( T2OLD( CV_NODJ_IPHA ) > T2OLDMAX_2ND_MC( CV_NODI_IPHA ) ) THEN
                        T2OLDMAX_2ND_MC( CV_NODI_IPHA ) = T2OLD( CV_NODJ_IPHA )
                     ENDIF
                     ENDIF
                     IF( T2OLDMIN_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                     IF( T2OLD( CV_NODJ_IPHA ) < T2OLDMIN_2ND_MC( CV_NODI_IPHA ) ) THEN
                        T2OLDMIN_2ND_MC( CV_NODI_IPHA ) = T2OLD( CV_NODJ_IPHA )
                     ENDIF
                     ENDIF
                  ENDIF
                  ! DEN:  
                  IF( DENMAX_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                  IF( DEN( CV_NODJ_IPHA ) > DENMAX_2ND_MC( CV_NODI_IPHA ) ) THEN
                     DENMAX_2ND_MC( CV_NODI_IPHA ) = DEN( CV_NODJ_IPHA ) 
                  ENDIF
                  ENDIF
                  IF( DENMIN_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                  IF( DEN( CV_NODJ_IPHA ) < DENMIN_2ND_MC( CV_NODI_IPHA ) ) THEN
                     DENMIN_2ND_MC( CV_NODI_IPHA ) = DEN( CV_NODJ_IPHA )
                  ENDIF
                  ENDIF

                  IF( DENOLDMAX_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                  IF( DENOLD( CV_NODJ_IPHA ) > DENOLDMAX_2ND_MC( CV_NODI_IPHA ) ) THEN
                     DENOLDMAX_2ND_MC( CV_NODI_IPHA ) = DENOLD( CV_NODJ_IPHA )
                  ENDIF
                  ENDIF
                  IF( DENOLDMIN_NOD( CV_NODI_IPHA ) .NE. CV_NODJ_IPHA ) THEN
                  IF( DENOLD( CV_NODJ_IPHA ) < DENOLDMIN_2ND_MC( CV_NODI_IPHA ) ) THEN
                     DENOLDMIN_2ND_MC( CV_NODI_IPHA ) = DENOLD( CV_NODJ_IPHA )
                  ENDIF
                  ENDIF

                ENDIF

               END DO
            ENDIF
         END DO

      END DO Loop2_CV_NODI


      ENDIF

      RETURN
    END SUBROUTINE SURRO_CV_MINMAX




    SUBROUTINE CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, SCVNGI, U_SLOC2LOC, CV_SLOC2LOC, &
         FACE_ELE, TOTELE, NFACE, CVFEM_ON_FACE, GI, &
         CV_NONODS, LOG_ON_BOUND, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, STOTEL, &
         CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN ) 
      ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC for a face on the
      ! boundary of the domain
      IMPLICIT NONE
      INTEGER, intent( in ) :: ELE, SCVNGI, CV_NONODS, TOTELE, STOTEL, CV_NLOC, U_NLOC, &
           CV_SNLOC, U_SNLOC, NFACE, GI, CV_ILOC
      INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( STOTEL * CV_SNLOC ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION( STOTEL * U_SNLOC ), intent( in ) :: U_SNDGLN
      INTEGER, intent( inout ) :: SELE, CV_SILOC
      INTEGER, DIMENSION( NFACE, TOTELE ), intent( in ) :: FACE_ELE
      LOGICAL, DIMENSION( CV_NLOC, SCVNGI ), intent( in )  :: CVFEM_ON_FACE
      INTEGER, DIMENSION( U_SNLOC ), intent( inout ) :: U_SLOC2LOC
      INTEGER, DIMENSION( CV_SNLOC ), intent( inout ) :: CV_SLOC2LOC
      LOGICAL, DIMENSION( CV_NONODS ), intent( inout ) :: LOG_ON_BOUND
      ! local variables
      INTEGER :: IFACE, ELE2, SELE2, CV_JLOC, CV_JNOD, &
           U_JLOC, U_JNOD, CV_KLOC, CV_SKNOD, &
           U_KLOC, U_SKLOC, U_SKNOD, CV_SKLOC, ngi, igi
      LOGICAL :: FOUND

      !ewrite(3,*)'In Calc_Sele'

      DO CV_JLOC = 1, CV_NLOC  
         CV_JNOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )
         !ewrite(3,*)'cv_jloc, gi, cvfem_on_face:', cv_jloc, gi, cvfem_on_face( cv_jloc, gi )
         LOG_ON_BOUND( CV_JNOD ) = CVFEM_ON_FACE( CV_JLOC, GI )
      END DO

      SELE = 0
      ! What face are we on        
      DO IFACE = 1, NFACE
         ELE2 = FACE_ELE( IFACE, ELE )
         SELE2 = MAX( 0, - ELE2 )
         !!SELE = SELE2
         ELE2 = MAX( 0, + ELE2 )
         IF( SELE2 /= 0 ) THEN
            FOUND = .TRUE.
            DO CV_SKLOC = 1, CV_SNLOC 
               CV_SKNOD = CV_SNDGLN(( SELE2 - 1 ) * CV_SNLOC + CV_SKLOC )
               !ewrite(3,*)'CV_SKNOD, LOG_ON_BOUND:',CV_SKNOD, LOG_ON_BOUND(CV_SKNOD)
               IF( .NOT. LOG_ON_BOUND( CV_SKNOD )) FOUND=.FALSE.
            END DO
            IF( FOUND ) SELE = SELE2
         END IF
      END DO

      ! Calculate CV_SLOC2LOC  
      Conditional_Sele: IF ( SELE /= 0 ) THEN   
         DO CV_SKLOC = 1, CV_SNLOC  
            CV_SKNOD = CV_SNDGLN(( SELE - 1 ) * CV_SNLOC + CV_SKLOC )
            DO CV_JLOC = 1, CV_NLOC  
               CV_JNOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_JLOC )
               IF( CV_SKNOD == CV_JNOD ) CV_KLOC = CV_JLOC
            END DO
            CV_SLOC2LOC( CV_SKLOC ) = CV_KLOC
            IF( CV_KLOC == CV_ILOC ) CV_SILOC = CV_SKLOC
         END DO
         !ewrite(3,*) 'cv_sloc2loc: ', cv_sloc2loc

         ! Calculate U_SLOC2LOC 
         DO U_SKLOC = 1, U_SNLOC  
            U_SKNOD = U_SNDGLN(( SELE - 1 ) * U_SNLOC + U_SKLOC )
            DO U_JLOC = 1, U_NLOC  
               U_JNOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_JLOC )
               IF( U_SKNOD == U_JNOD ) U_KLOC = U_JLOC
            END DO
            U_SLOC2LOC( U_SKLOC ) = U_KLOC
         END DO
      END IF Conditional_Sele

      RETURN   
    END SUBROUTINE CALC_SELE



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
         SUFEN, SCVDETWEI, CVNORMX, CVNORMY, CVNORMZ, DEN, CV_NODI, CV_NODI_IPHA, &
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
      REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DEN
      REAL, intent( in ) :: NDOTQ, NDOTQOLD, LIMDT, LIMDTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD

      ! Local variables...
      INTEGER :: U_KLOC, U_KLOC2, JCOUNT_IPHA, IDIM, U_NODK, U_NODK_IPHA, JCOUNT2_IPHA, &
           U_KLOC_LEV, U_NLOC_LEV
      REAL :: RCON,UDGI_IMP,VDGI_IMP,WDGI_IMP,NDOTQ_IMP

      !ewrite(3,*)' In PUT_IN_CT_RHS CVNORMX/Y/Z:',CVNORMX( GI ),CVNORMY( GI ),CVNORMZ( GI ), &
      !     ':', CVNORMX( GI )**2+CVNORMY( GI )**2+CVNORMZ( GI )**2
      !ewrite(3,*)' SCVDETWEI( GI ):',SCVDETWEI( GI )
      !ewrite(3,*)' SUFEN( :, GI ):',SUFEN( :, GI )
      !ewrite(3,*)' jcount_kloc, limdt, theta:', jcount_kloc, limdt, ftheta_t2

      !ewrite(3,*)' ugi_coef_ele:', ugi_coef_ele
      !ewrite(3,*)' vgi_coef_ele:', vgi_coef_ele
      !ewrite(3,*)' wgi_coef_ele:', wgi_coef_ele
      !ewrite(3,*)' ugi_coef_ele2:', ugi_coef_ele2
      !ewrite(3,*)' vgi_coef_ele2:', vgi_coef_ele2
      !ewrite(3,*)' wgi_coef_ele2:', wgi_coef_ele2

      !ewrite(3,*)' DEN, DETWEI', DEN( CV_NODI_IPHA ), SCVDETWEI( GI ) 
      !ewrite(3,*)' FTHETA_T2,  FTHETA_T2OLD, LIMDT, LIMDTOLD', FTHETA_T2, ONE_M_FTHETA_T2OLD, LIMDT, LIMDTOLD
      !ewrite(3,*)' ndotq, ndotqold', ndotq, ndotqold


      DO U_KLOC = 1, U_NLOC

         JCOUNT_IPHA = JCOUNT_KLOC( U_KLOC )  +  (IPHASE - 1 ) * NCOLCT * NDIM

         RCON    = SCVDETWEI( GI ) * FTHETA_T2 * LIMDT  &
              * SUFEN( U_KLOC, GI ) / DEN( CV_NODI_IPHA )

         !ewrite(3,*) 'U_KLOC, CV_NODI_IPHA, JCOUNT_IPHA, RCON:', &
         !     U_KLOC, CV_NODI_IPHA, JCOUNT_IPHA, RCON

         IDIM = 1
         CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) &
              =  CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT) &
              +  RCON * UGI_COEF_ELE(U_KLOC) * CVNORMX( GI )

         IDIM = 2
         IF( NDIM >= 2 ) &
              CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) &
              =  CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) &
              +  RCON * VGI_COEF_ELE(U_KLOC) * CVNORMY( GI )

         IDIM = 3
         IF( NDIM >= 3 ) &
              CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) &
              =  CT( JCOUNT_IPHA + ( IDIM - 1 ) * NCOLCT ) &
              +  RCON * WGI_COEF_ELE(U_KLOC) * CVNORMZ( GI )
      END DO

      IF(SELE /= 0) THEN
         UDGI_IMP=0.0
         VDGI_IMP=0.0
         WDGI_IMP=0.0
         !DO U_KLOC_LEV = 1, U_NLOC_LEV
         !U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
         DO U_KLOC = 1, U_NLOC
            U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
            U_NODK_IPHA = U_NODK + (IPHASE-1)*U_NONODS
            UDGI_IMP=UDGI_IMP + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE(U_KLOC) * NU(U_NODK_IPHA) 
            VDGI_IMP=VDGI_IMP + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE(U_KLOC) * NV(U_NODK_IPHA) 
            WDGI_IMP=WDGI_IMP + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE(U_KLOC) * NW(U_NODK_IPHA) 
         END DO

         NDOTQ_IMP=CVNORMX( GI ) * UDGI_IMP + CVNORMY( GI ) * VDGI_IMP + CVNORMZ( GI ) * WDGI_IMP

         !ewrite(3,*) 'CVNORMX( GI ),UDGI_IMP,CVNORMY( GI ),VDGI_IMP,CVNORMZ( GI ),WDGI_IMP:', &
         !     CVNORMX( GI ),UDGI_IMP,CVNORMY( GI ),VDGI_IMP,CVNORMZ( GI ),WDGI_IMP
         !ewrite(3,*) 'NDOTQOLD,NDOTQ,NDOTQ_IMP:',NDOTQOLD,NDOTQ,NDOTQ_IMP
         !ewrite(3,*)' DEN, DETWEI', DEN( CV_NODI_IPHA ), SCVDETWEI( GI ) 
         !ewrite(3,*)' FTHETA_T2,  ONE_M_FTHETA_T2OLD, LIMDT, LIMDTOLD', FTHETA_T2, ONE_M_FTHETA_T2OLD, LIMDT, LIMDTOLD
         !ewrite(3,*)' NDOTQOLD, NDOTQ, NDOTQ_IMP', NDOTQOLD, NDOTQ, NDOTQ_IMP 

         CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - SCVDETWEI( GI ) * ( &
              ONE_M_FTHETA_T2OLD * LIMDTOLD * NDOTQOLD &
              + FTHETA_T2  * LIMDT * (NDOTQ-NDOTQ_IMP) &
              ) / DEN( CV_NODI_IPHA )

      ELSE

         !ewrite(3,*)' DEN, DETWEI', DEN( CV_NODI_IPHA ), SCVDETWEI( GI ) 
         !ewrite(3,*)' ONE_M_FTHETA_T2OLD, LIMDTOLD, NDOTQOLD', FTHETA_T2, ONE_M_FTHETA_T2OLD, LIMDTOLD, NDOTQOLD

         CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - SCVDETWEI( GI ) * ( &
              ONE_M_FTHETA_T2OLD * LIMDTOLD * NDOTQOLD &
              ) / DEN( CV_NODI_IPHA ) 
      ENDIF

      IF((ELE2 /= 0).AND.(ELE2 /= ELE)) THEN
         ! We have a discontinuity between elements so integrate along the face...
         !DO U_KLOC_LEV = 1, U_NLOC_LEV
         !U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
         DO U_KLOC = 1, U_NLOC
            U_KLOC2=U_OTHER_LOC(U_KLOC)
            IF(U_KLOC2 /= 0) THEN
               JCOUNT2_IPHA = JCOUNT_KLOC2( U_KLOC2 ) + (IPHASE - 1) * NCOLCT * NDIM

               RCON = SCVDETWEI( GI ) * FTHETA_T2 * LIMDT  &
                    * SUFEN( U_KLOC, GI ) / DEN( CV_NODI_IPHA )

               !ewrite(3,*) 'U_KLOC2, CV_NODI_IPHA, JCOUNT2_IPHA, RCON UGI:', &
               !     U_KLOC2, CV_NODI_IPHA, JCOUNT2_IPHA, RCON, UGI_COEF_ELE2( U_KLOC2 )

               IDIM = 1
               CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                    =  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                    +  RCON * UGI_COEF_ELE2( U_KLOC2 ) * CVNORMX( GI )

               IDIM = 2
               IF( NDIM >= 2 ) &
                    CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                    =  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                    +  RCON * VGI_COEF_ELE2( U_KLOC2 ) * CVNORMY( GI )

               IDIM = 3
               IF( NDIM >= 3 ) &
                    CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                    =  CT( JCOUNT2_IPHA + ( IDIM - 1 ) * NCOLCT ) &
                    +  RCON * WGI_COEF_ELE2( U_KLOC2 ) * CVNORMZ( GI )
            ENDIF
         END DO
      ENDIF

      !ewrite(3,*)'check ct:,', NCOLCT * NDIM * NPHASE, size(ct), ct(1:NCOLCT * NDIM * NPHASE)

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

      !ewrite(3,*)'satura :',satura

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
               !        ABSC=0.5*ABS_CV_NODI_IPHA + GRAD_ABS_CV_NODI_IPHA*(SAT_FEM_IPHA - SAT_CV_NODI_IPHA) &
               !            +0.5*ABS_CV_NODJ_IPHA + GRAD_ABS_CV_NODJ_IPHA*(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA) 
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
            !     w=0.7
            !     WA=(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)/(SAT_CV_NODJ_IPHA-SAT_CV_NODI_IPHA)
            !     W=(SAT_FEM_IPHA - SAT_CV_NODI_IPHA)/(SAT_CV_NODJ_IPHA-SAT_CV_NODI_IPHA)
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
            !     w=wa
            !     w=0.8
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
            !     w=0.7
            !     WA=(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
            !     W=(SAT_FEM_IPHA - SAT_CV_NODJ_IPHA)/(SAT_CV_NODI_IPHA-SAT_CV_NODJ_IPHA)
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
            !     w=wa
            !     w=0.8
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
               !ewrite(3,*)'ABS_CV_NODJ, SAT_CV_NODJ, IPHASE:',ABS_CV_NODJ, SAT_CV_NODJ, IPHASE
               B=(ABS_CV_NODJ-ABSC)/(SAT_CV_NODJ-SATC)
               A=ABS_CV_NODJ - B*SAT_CV_NODJ

               ABSC=a+b*0.5*(SAT_CV_NODJ+SAT_CV_NODI)
               !ewrite(3,*)'guessed valueS ABSC,ABS_CV_NODJ=',ABSC,ABS_CV_NODJ
               !ewrite(3,*)'SAT_CV_NODI,SAT_CV_NODj:',SAT_CV_NODI,SAT_CV_NODj
               !ewrite(3,*)'SAT_CV_NODI_ipha,SAT_CV_NODj_ipha:',SAT_CV_NODI_ipha,SAT_CV_NODj_ipha
               W=(ABS_CV_NODI-ABSC)/tolfun(ABS_CV_NODI-ABS_CV_NODJ)
               SATC=W*SAT_CV_NODJ + (1.-W)*SAT_CV_NODI
               W=(SATC - SAT_CV_NODI)/(SAT_CV_NODJ - SAT_CV_NODI)
               !ewrite(3,*)'1w=',w
               !w=max(w,1.0-w) ! make any non-linear variation subject to upwinding
               !w=1.-w ! this is correct
               w=0.5 + (w-0.5)*2.
               W=max(w,0.5)
               !W=max(w,0.0)
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
               !ewrite(3,*)'w,income:',w,income
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
      !EWRITE(3,*)'SATC,ABSC=',SATC,ABSC

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




    SUBROUTINE CALC_ANISOTROP_LIM(&
         ! Caculate the upwind values stored in matrix form...
         T,TOLD,DEN,DENOLD,T2,T2OLD, &
         TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, &
         T2UPWIND_MAT, T2OLDUPWIND_MAT, &
         ! Store the upwind element for interpolation and its weights for 
         ! faster results...
         IGOT_T2,NPHASE,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
         SMALL_FINDRM,SMALL_CENTRM,SMALL_COLM,NSMALL_COLM, &
         X_NDGLN,X_NONODS,NDIM, &
         X,Y,Z, &
         FINACV,COLACV,NCOLACV) 
      ! For the anisotropic limiting scheme we find the upwind values
      ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
      ! value for each node pair is stored in the matrices TUPWIND AND
      IMPLICIT NONE
      INTEGER, intent( in ) :: CV_NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, &
           NSMALL_COLM, NDIM,NCOLACV,IGOT_T2,NPHASE
      REAL, DIMENSION( CV_NONODS*NPHASE ), intent( in ) :: T,TOLD,DEN,DENOLD
      REAL, DIMENSION( CV_NONODS*NPHASE*IGOT_T2), intent( in ) :: T2,T2OLD
      REAL, DIMENSION( NSMALL_COLM*NPHASE ), intent( inout ) :: TUPWIND_MAT, TOLDUPWIND_MAT, &
           DENUPWIND_MAT, DENOLDUPWIND_MAT
      REAL, DIMENSION( NSMALL_COLM*NPHASE*IGOT_T2 ), intent( inout ) :: T2UPWIND_MAT, T2OLDUPWIND_MAT
      INTEGER, DIMENSION( TOTELE*X_NLOC ), intent( in ) :: X_NDGLN
      INTEGER, DIMENSION( TOTELE*CV_NLOC ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( CV_NONODS+1 ), intent( inout ) :: SMALL_FINDRM
      INTEGER, DIMENSION( NSMALL_COLM ), intent( inout ) :: SMALL_COLM
      INTEGER, DIMENSION( CV_NONODS ), intent( inout ) :: SMALL_CENTRM
      REAL, DIMENSION( X_NONODS ), intent( in ) :: X,Y,Z
      INTEGER, DIMENSION( CV_NONODS*NPHASE+1 ), intent( in ) :: FINACV
      INTEGER, DIMENSION( NCOLACV ), intent( in ) :: COLACV
      REAL, DIMENSION(:), ALLOCATABLE :: SOL

      ! Allocate memory 
      INTEGER :: COUNT, COUNT2, CV_NOD
      INTEGER :: IDUM(1)
      REAL :: RDUM(1)

      COUNT2=0
      DO CV_NOD=1,CV_NONODS
         SMALL_FINDRM(CV_NOD)=COUNT2+1
         DO COUNT=FINACV(CV_NOD),FINACV(CV_NOD+1)-1
            IF(COLACV(COUNT).LE.CV_NONODS) THEN
               COUNT2=COUNT2+1
               SMALL_COLM(COUNT2)=COLACV(COUNT)
               IF(SMALL_COLM(COUNT2)==CV_NOD) SMALL_CENTRM(CV_NOD)=COUNT2
            ENDIF
         END DO
      END DO
      SMALL_FINDRM(CV_NONODS+1)=COUNT2+1

      ! Allocate memory and find upwind field values for limiting...
      IF(IGOT_T2.NE.0) THEN

         allocate( sol( 6*nsmall_colm*nphase) )
         sol = (/TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT/)

         ! Obtain the weights
         CALL CALC_ANISOTROP_LIM_VALS( &
              ! Caculate the upwind values stored in matrix form...
              (/T,TOLD,DEN,DENOLD,T2,T2OLD/), &
              SOL,  &
              ! Store the upwind element for interpolation and its weights for 
              ! faster results...
              IDUM,RDUM, 0,  &
              NPHASE*6,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
              SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
              X_NDGLN,X_NONODS,NDIM, &
              X,Y,Z, &
              .FALSE., .FALSE.)

         TUPWIND_MAT = sol( 1 : nsmall_colm*nphase )
         TOLDUPWIND_MAT= sol( 1+nsmall_colm*nphase : 2*nsmall_colm*nphase )
         DENUPWIND_MAT = sol( 1+2*nsmall_colm*nphase : 3*nsmall_colm*nphase )
         DENOLDUPWIND_MAT = sol(1+3*nsmall_colm*nphase : 4*nsmall_colm*nphase )
         T2UPWIND_MAT = sol( 1+4*nsmall_colm*nphase : 5*nsmall_colm*nphase )
         T2OLDUPWIND_MAT = sol( 1+5*nsmall_colm*nphase : 6*nsmall_colm*nphase )

         deallocate( sol )

      ELSE

         allocate( sol( 4*nsmall_colm*nphase ) )
         sol = (/TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT/)

         CALL CALC_ANISOTROP_LIM_VALS( &
              ! Caculate the upwind values stored in matrix form...
              (/T,TOLD,DEN,DENOLD/),&
              SOL,  &
              ! Store the upwind element for interpolation and its weights for 
              ! faster results...
              IDUM,RDUM, 0,  &
              NPHASE*4,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
              SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
              X_NDGLN,X_NONODS,NDIM, &
              X,Y,Z, &
              .FALSE., .FALSE.)

         TUPWIND_MAT = sol( 1 : nsmall_colm*nphase )
         TOLDUPWIND_MAT = sol( 1 + nsmall_colm*nphase : 2*nsmall_colm*nphase )
         DENUPWIND_MAT = sol( 1+2*nsmall_colm*nphase : 3*nsmall_colm*nphase )
         DENOLDUPWIND_MAT = sol( 1+3*nsmall_colm*nphase : 4*nsmall_colm*nphase )

         deallocate( sol )

      ENDIF

    END SUBROUTINE CALC_ANISOTROP_LIM




    SUBROUTINE CALC_ANISOTROP_LIM_VALS( &
         ! Caculate the upwind values stored in matrix form...
         T, &
         TUPWIND,  &
         ! Store the upwind element for interpolation and its weights for 
         ! faster results...
         ELEMATPSI,ELEMATWEI, IELEMATPSI,  &
         NFIELD,NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
         SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
         X_NDGLN,X_NONODS,NDIM, &
         X,Y,Z, &
         STORE_ELE, RET_STORE_ELE)
      ! For the anisotropic limiting scheme we find the upwind values
      ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
      ! value for each node pair is stored in the matrices TUPWIND AND
      IMPLICIT NONE
      INTEGER :: NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, NSMALL_COLM, NFIELD,NDIM
      REAL, DIMENSION( NONODS*NFIELD ), intent( in ) :: T
      REAL, DIMENSION( NSMALL_COLM*NFIELD ), intent( inout ) :: TUPWIND
      INTEGER :: X_NDGLN(TOTELE*X_NLOC), CV_NDGLN(TOTELE*CV_NLOC) 
      INTEGER :: SMALL_FINDRM(NONODS+1),SMALL_COLM(NSMALL_COLM),SMALL_CENTRM(NONODS)
      REAL    :: X(X_NONODS),Y(X_NONODS),Z(X_NONODS)

      INTEGER :: ELEMATPSI( NSMALL_COLM * IELEMATPSI ), IELEMATPSI      
      REAL :: ELEMATWEI( NSMALL_COLM * CV_NLOC * IELEMATPSI )
      LOGICAL, INTENT(IN) :: STORE_ELE, RET_STORE_ELE

      LOGICAL :: D3,DCYL
      ! Allocate memory for the interpolated upwind values
      real, dimension( :, : ), allocatable :: N,NLX,NLY,NLZ
      real, dimension( :, : ), allocatable :: UN, UNLX, UNLY, UNLZ, CVN
      real, dimension( : ), allocatable :: WEIGHT
      integer, dimension( : ), allocatable :: SUB_NDGLNO
      INTEGER :: COUNT, COUNT2, NOD, SUB_TOTELE, NGI,NLOC
      REAL L1(50),L2(50),L3(50),L4(50)

! **********************Calculate linear shape functions...
      IF(NDIM==1) THEN
         NLOC=1
         NGI=2
      ELSE IF(NDIM==2) THEN
         NLOC=3
         NGI=3
      ELSE IF(NDIM==3) THEN 
         NLOC=4
         NGI=4
      ENDIF
      ALLOCATE(N(NLOC,NGI), NLX(NLOC,NGI), NLY(NLOC,NGI), NLZ(NLOC,NGI))
      ALLOCATE(WEIGHT(NGI))
!
! Shape functions for triangles and tets...
      CALL TRIQUAold(L1, L2, L3, L4, WEIGHT, ndim==3,NGI)
      ! Work out the shape functions and there derivatives...
      CALL SHATRIold(L1, L2, L3, L4, WEIGHT, ndim==3,&
           &              NLOC,NGI,&
           &              N,NLX,NLY,NLZ)

! ******************************************************************
! Calculate the sub elements for quadratic element SUB_NDGLNO ... 
      IF(CV_NLOC==NLOC) THEN 
         SUB_TOTELE=TOTELE
      ELSE
         IF(NDIM==1) THEN
            SUB_TOTELE=2*TOTELE
         ELSE IF(NDIM==2) THEN
            SUB_TOTELE=4*TOTELE
         ELSE IF(NDIM==3) THEN 
            SUB_TOTELE=10*TOTELE
         ENDIF
      ENDIF

      ALLOCATE(SUB_NDGLNO(SUB_TOTELE*NLOC))

      !print *,'CV_nLOC,NLOC:',CV_NLOC,NLOC
      IF(CV_NLOC==NLOC) THEN 
         SUB_NDGLNO=CV_NDGLN
      ELSE
         stop 2921
         ! do something...
      ENDIF

! Calculate the sub elements for quadratic element SUB_NDGLNO ... 
! ******************************************************************


      CALL CALC_ANISOTROP_LIM_VALS2( &
! Caculate the upwind values stored in matrix form...
           T, &
           TUPWIND,  &
! Store the upwind element for interpolation and its weights for 
! faster results...
           ELEMATPSI,ELEMATWEI, IELEMATPSI,  &
           NFIELD,NONODS,NLOC,NGI,SUB_TOTELE,SUB_NDGLNO, &
           SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
           X_NDGLN,X_NONODS,NDIM, &
           X,Y,Z, &
           N,NLX,NLY,NLZ, WEIGHT, &
           STORE_ELE, RET_STORE_ELE)

      RETURN
    END SUBROUTINE CALC_ANISOTROP_LIM_VALS



       SUBROUTINE CALC_ANISOTROP_LIM_VALS2( &
! Caculate the upwind values stored in matrix form...
           T, &
           TUPWIND,  &
! Store the upwind element for interpolation and its weights for 
! faster results...
           ELEMATPSI,ELEMATWEI, IELEMATPSI,  &
           NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
           FINDRM,COLM,NCOLM, &
           X_NDGLN,X_NONODS,NDIM, &
           X,Y,Z, &
           N,NLX,NLY,NLZ, WEIGHT, &
           STORE_ELE, RET_STORE_ELE)
        ! For the anisotropic limiting scheme we find the upwind values
        ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
        ! value for each node pair is stored in the matrices TUPWIND AND
      IMPLICIT NONE
      INTEGER :: NONODS,X_NONODS,TOTELE,NLOC,NGI,NCOLM,NFIELD,NDIM
       REAL, DIMENSION( NONODS*NFIELD ), intent( in ) :: T
       REAL, DIMENSION( NCOLM*NFIELD ), intent( inout ) :: TUPWIND
      INTEGER :: NDGLNO(TOTELE*NLOC),X_NDGLN(TOTELE*NLOC)
      INTEGER :: FINDRM(NONODS+1),COLM(NCOLM),CENTRM(NONODS)
      INTEGER :: ELEMATPSI( NCOLM * IELEMATPSI ), IELEMATPSI      
      REAL :: ELEMATWEI( NCOLM * NLOC * IELEMATPSI )

      REAL    :: X(X_NONODS),Y(X_NONODS),Z(X_NONODS)
      REAL, INTENT(IN) :: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL, INTENT(IN) :: WEIGHT(NGI) 
      LOGICAL, INTENT(IN) :: STORE_ELE, RET_STORE_ELE


      LOGICAL :: D3,DCYL
        ! Allocate memory for the interpolated upwind values
      LOGICAL, PARAMETER :: BOUND   = .TRUE.,  REFLECT = .TRUE. ! limiting options
      INTEGER, DIMENSION( : ), allocatable :: NOD_FINDELE,NOD_COLELE, NLIST, INLIST, DUMMYINT
      REAL, DIMENSION( : ), allocatable :: DUMMYREAL
      INTEGER MXNCOLEL,NCOLEL

        ! Over-estimate the size of the COLELE array
        MXNCOLEL=20*TOTELE+500

        ALLOCATE( NOD_FINDELE(NONODS+1) )
        ALLOCATE( NOD_COLELE(MXNCOLEL) )
        
        ALLOCATE(   NLIST(NONODS) )
        ALLOCATE(  INLIST(NONODS) )
        
        ! Calculate node element list - moved from (I)FINPTS
        CALL PHILNODELE(NONODS,NOD_FINDELE,NOD_COLELE, &
                        NCOLEL,MXNCOLEL, &
                        TOTELE,NLOC,NDGLNO, &
                        NLIST,INLIST)
        
        IF(STORE_ELE) THEN
    
          CALL FINPTSSTORE(T,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                  TUPWIND,FINDRM,COLM,NCOLM,NDIM, &
                  X_NDGLN,X_NONODS, &
                  X,Y,Z, &
                  N,NLX,NLY,NLZ, WEIGHT, &
                  NOD_FINDELE,NOD_COLELE,NCOLEL, &
                  ELEMATPSI,ELEMATWEI,1, &
                  BOUND, REFLECT)
    
        ELSE IF(RET_STORE_ELE) THEN 
! Find the weights for the interpolation
! This does depend on the solns T when BOUND...
          
          CALL GETSTOREELEWEI(T,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
                    TUPWIND,FINDRM,COLM,NCOLM,BOUND, &
                    ELEMATPSI,ELEMATWEI)
    
        ELSE ! assume we have not stored anything (elements or weights)...
           ALLOCATE(DUMMYINT(NCOLM))
           ALLOCATE(DUMMYREAL(NCOLM*NLOC))

          CALL FINPTSSTORE(T,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                    TUPWIND,FINDRM,COLM,NCOLM,NDIM, &
                    X_NDGLN,X_NONODS, &
                    X,Y,Z, &
                    N,NLX,NLY,NLZ, WEIGHT, &
                    NOD_FINDELE,NOD_COLELE,NCOLEL, &
                    DUMMYINT,DUMMYREAL,0, &
                    BOUND, REFLECT)
 
          DEALLOCATE(DUMMYINT,DUMMYREAL) 

        ENDIF

        DEALLOCATE( NOD_FINDELE, NOD_COLELE, NLIST, INLIST )

      END SUBROUTINE CALC_ANISOTROP_LIM_VALS2
!
!     
!
!
        SUBROUTINE GETSTOREELEWEI(PSI,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
     &     MATPSI,FINDRM,COLM,NCOLM,BOUND,&
     &     ELEMATPSI,ELEMATWEI)
! use the stored interpolation coeffs to caclulate MATPSI.
!     This sub finds the matrix values MATPSI for a given point on the 
!     stencil 
      IMPLICIT NONE
        REAL FRALINE
        LOGICAL BOUND
        PARAMETER(FRALINE=0.001)
        INTEGER NFIELD,NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
        REAL PSI(NONODS*NFIELD)
        INTEGER NCOLM
        INTEGER FINDRM(NONODS+1),COLM(NCOLM)
        REAL MATPSI(NCOLM*NFIELD)
        INTEGER ELEMATPSI(NCOLM)
        REAL ELEMATWEI(NCOLM*NLOC)
!  LOCAL VARIABLES...
        INTEGER NOD,COUNT,ELEWIC,ILOC,INOD,IFIELD
        INTEGER KNOD,COUNT2,JNOD
        REAL RMATPSI,RMATPSIOLD
        REAL, ALLOCATABLE, DIMENSION(:)::MINPSI
        REAL, ALLOCATABLE, DIMENSION(:)::MAXPSI
        
        if(bound) then
          ALLOCATE(MINPSI(TOTELE*NFIELD))
          ALLOCATE(MAXPSI(TOTELE*NFIELD))
        
! find the max and min local to each element...
        !print *,'going into minmaxelewic'
          CALL MINMAXELEWIC(PSI,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
     &     FINDRM,COLM,NCOLM,&
     &     MINPSI,MAXPSI)
        !print *,'out of minmaxelewic'
        endif
        
      do NOD=1,NONODS! Was loop 
!         print *,'nod=',nod
      do COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 
!            print *,'count=',count
            IF(NOD.NE.COLM(COUNT)) THEN
              ELEWIC=ELEMATPSI(COUNT)
!              print *,'elewic=',elewic
              DO IFIELD=1,NFIELD
                 RMATPSI=0.0
                 RMATPSIOLD=0.0
                 DO ILOC=1,NLOC! Was loop 
                    INOD=NDGLNO((ELEWIC-1)*NLOC+ILOC)
                    RMATPSI=RMATPSI+ELEMATWEI((COUNT-1)*NLOC+ILOC)*PSI(INOD+(IFIELD-1)*NONODS)
                 END DO
              
                 RMATPSI   =PSI(NOD+(IFIELD-1)*NONODS)   &
                      +(1./FRALINE)*(RMATPSI   -PSI(NOD+(IFIELD-1)*NONODS))

! make locally bounded...
                 if(bound) then
                    MATPSI(COUNT+(IFIELD-1)*NCOLM)   &
      =MAX(MIN(RMATPSI,   MAXPSI(ELEWIC+(IFIELD-1)*TOTELE)),   &
     &                            MINPSI(ELEWIC+(IFIELD-1)*TOTELE))
                 else
                    MATPSI(COUNT+(IFIELD-1)*NCOLM)   =RMATPSI
                 endif
              END DO
            ENDIF
          END DO
        END DO
        RETURN
        
  end subroutine getstoreelewei
        
!
!   
!
!
        SUBROUTINE MINMAXELEWIC(PSI,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
     &     FINDRM,COLM,NCOLM,&
     &     MINPSI,MAXPSI)
! This sub calculates the max and min values of PSI in local vacinity of 
! an element. 
      IMPLICIT NONE
        REAL FRALINE
        PARAMETER(FRALINE=0.001)
        INTEGER NFIELD,NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
        REAL PSI(NONODS*NFIELD)
        INTEGER NCOLM
        INTEGER FINDRM(NONODS+1),COLM(NCOLM)
        REAL MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
!  LOCAL VARIABLES...
        INTEGER NOD,COUNT,ELEWIC,ILOC,INOD,IFIELD
        INTEGER KNOD,COUNT2,JNOD
        REAL RMATPSI
        
        MINPSI   =1.E+20
        MAXPSI   =-1.E+20
! find the max and min local to each element...
      do ELEWIC=1,TOTELE! Was loop 
      do ILOC=1,NLOC! Was loop 
              KNOD=NDGLNO((ELEWIC-1)*NLOC+ILOC)
!     Search around node KNOD for max and min PSI...
      do COUNT2=FINDRM(KNOD),FINDRM(KNOD+1)-1! Was loop 
                  JNOD=COLM(COUNT2)
                  DO IFIELD=1,NFIELD
                     MINPSI(ELEWIC+(IFIELD-1)*TOTELE)  &
        =MIN(PSI(JNOD+(IFIELD-1)*NONODS),   MINPSI(ELEWIC+(IFIELD-1)*TOTELE))
                     MAXPSI(ELEWIC+(IFIELD-1)*TOTELE)  &
        =MAX(PSI(JNOD+(IFIELD-1)*NONODS),   MAXPSI(ELEWIC+(IFIELD-1)*TOTELE))
                  END DO
              END DO
           END DO
        END DO
        RETURN
        
  end subroutine minmaxelewic
        
!     
!     
!     
!     
  SUBROUTINE FINPTSSTORE(PSI,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
       &     MATPSI,FINDRM,COLM,NCOLM,NDIM, &
       &     X_NDGLN,X_NONODS, &
       &     X,Y,Z,&
       &     N,NLX,NLY,NLZ, WEIGHT,&
       !     work space...
       &     FINDELE,COLELE,NCOLEL,&
       &     ELEMATPSI,ELEMATWEI,IGETSTOR,&
       &     BOUND, REFLECT)
    !     This sub finds the matrix values MATPSI for a given point on the 
    !     stencil 
    ! IF IGETSTOR=1 then get ELEMATPSI,ELEMATWEI.
    IMPLICIT NONE
    LOGICAL BOUND,REFLECT
    ! IF REFLECT then use a reflection condition at boundary to 
    ! do limiting. 
    INTEGER NFIELD,NONODS,NLOC,NGI,TOTELE,NDIM,NDGLNO(TOTELE*NLOC)
    REAL PSI(NONODS*NFIELD)
    INTEGER NCOLM,NCOLEL
    INTEGER FINDRM(NONODS+1),COLM(NCOLM)
    REAL MATPSI(NCOLM*NFIELD)
    INTEGER X_NDGLN(TOTELE*NLOC),X_NONODS
    REAL X(X_NONODS),Y(X_NONODS),Z(X_NONODS)
    REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL WEIGHT(NGI)
    !     work space...
    INTEGER FINDELE(NONODS+1),COLELE(NCOLEL)
    INTEGER IGETSTOR
    INTEGER ELEMATPSI(NCOLM*IGETSTOR)
    REAL ELEMATWEI(NCOLM*NLOC*IGETSTOR)
    ! ELEWIC is the element to do interpolation from
    ! LOCCORDSK contains the weights. 
    !     Local variables...
    INTEGER NOD,COUNT,NODI,NODJ,ILOC,GI,ELE
    INTEGER ELEWIC,XNOD,XNODJ,IFIELD
    REAL LOCCORDSK(NLOC)
    REAL NORMX1,NORMY1,NORMZ1
    REAL VOLUME,INVH,LENG
    !     work space...
    REAL DETWEI(NGI),RA(NGI)
    REAL NX(NLOC,NGI),NY(NLOC,NGI),NZ(NLOC,NGI)
    REAL, ALLOCATABLE, DIMENSION(:)::NORMX
    REAL, ALLOCATABLE, DIMENSION(:)::NORMY
    REAL, ALLOCATABLE, DIMENSION(:)::NORMZ
    REAL, ALLOCATABLE, DIMENSION(:)::MLUM
    REAL, ALLOCATABLE, DIMENSION(:)::MINPSI
    REAL, ALLOCATABLE, DIMENSION(:)::MAXPSI
    INTEGER, ALLOCATABLE, DIMENSION(:)::NOD2XNOD

    NORMX1=0.0
    NORMY1=0.0
    NORMZ1=0.0
    IF(REFLECT) THEN
       !     calculate normals...********************
       ALLOCATE(NORMX(NONODS))
       ALLOCATE(NORMY(NONODS))
       ALLOCATE(NORMZ(NONODS))
       ALLOCATE(MLUM(NONODS))
       NORMX(1:NONODS) = 0.0
       NORMY(1:NONODS) = 0.0
       NORMZ(1:NONODS) = 0.0
       MLUM(1:NONODS) = 0.0
       DO ELE=1,TOTELE! Was loop 
          !     Calculate DETWEI,RA,NX,NY,NZ for element ELE
          CALL DETNLXR(ELE, X,Y,Z, X_NDGLN, TOTELE,X_NONODS,NLOC,NGI, &
               &        N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME,NDIM==1,NDIM==3,.FALSE., &
               &        NX,NY,NZ) 
          !     
          DO ILOC=1,NLOC! Was loop 
             NODI=NDGLNO((ELE-1)*NLOC+ILOC)
             DO GI=1,NGI! Was loop 
                NORMX(NODI)=NORMX(NODI)+NX(ILOC,GI)*DETWEI(GI)
                NORMY(NODI)=NORMY(NODI)+NY(ILOC,GI)*DETWEI(GI)
                IF(NDIM==3)NORMZ(NODI)=NORMZ(NODI)+NZ(ILOC,GI)*DETWEI(GI)
                MLUM(NODI) =MLUM(NODI) +N(ILOC,GI) *DETWEI(GI)
             END DO
          END DO
       END DO
       !     Renormalise
       DO NODI=1,NONODS! Was loop 
          INVH=(ABS(NORMX(NODI))+ABS(NORMY(NODI))+ABS(NORMZ(NODI)))&
               &          /MLUM(NODI)
          IF(INVH.GT.1.E-5) THEN
             LENG=SQRT(NORMX(NODI)**2+NORMY(NODI)**2+NORMZ(NODI)**2)
             NORMX(NODI)=NORMX(NODI)/LENG
             NORMY(NODI)=NORMY(NODI)/LENG
             NORMZ(NODI)=NORMZ(NODI)/LENG
          ELSE
             NORMX(NODI)=0.0
             NORMY(NODI)=0.0
             NORMZ(NODI)=0.0
          ENDIF
       END DO
       ! In parallel need to distribute NORMX,NORMY,NORMZ
       !     calculate normals...************************
       ! ENDOF IF(REFLECT) THEN...    
    ENDIF

    ALLOCATE(MINPSI(TOTELE*NFIELD))
    ALLOCATE(MAXPSI(TOTELE*NFIELD))

    IF(BOUND) THEN
       ! find the max and min local to each element...
       CALL MINMAXELEWIC(PSI,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
            &     FINDRM,COLM,NCOLM,&
            &     MINPSI,MAXPSI)
    ENDIF

    !     
    !     Calculate node element list. 

    ALLOCATE(NOD2XNOD(NONODS))
    DO ELE=1,TOTELE! Was loop 
       DO ILOC=1,NLOC! Was loop 
          NOD =NDGLNO((ELE-1)*NLOC+ILOC)
          XNOD=X_NDGLN((ELE-1)*NLOC+ILOC)
          NOD2XNOD(NOD)=XNOD
       END DO
    END DO
    !     
    DO NOD=1,NONODS! Was loop 10

       XNOD=NOD2XNOD(NOD)
       !     
       DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1! Was loop 20

          NODJ=COLM(COUNT)
          XNODJ=NOD2XNOD(NODJ)
          DO IFIELD=1,NFIELD
             MATPSI(COUNT+(IFIELD-1)*NCOLM)=0.
          END DO
          !     
          IF(NOD.NE.NODJ) THEN
             IF(REFLECT) THEN
                NORMX1=NORMX(NOD)
                NORMY1=NORMY(NOD)
                NORMZ1=NORMZ(NOD)
             ENDIF
             CALL MATPTSSTORE(MATPSI,COUNT,NFIELD,NOD,&
                  &              PSI,NONODS,X_NONODS,&
                  &              NLOC,TOTELE,X_NDGLN,NDGLNO,&
                  &              FINDRM,COLM,NCOLM,&
                  &              X(XNOD),Y(XNOD),Z(XNOD),&
                  &              X(XNODJ),Y(XNODJ),Z(XNODJ),&
                  &              NORMX1,NORMY1,NORMZ1,&
                  &              X,Y,Z,&
                  !     work space...
                  &              FINDELE,COLELE,NCOLEL, &
                  &              MINPSI,MAXPSI, &
                  &              ELEWIC,LOCCORDSK,BOUND,REFLECT,NDIM)
             IF(IGETSTOR.EQ.1) THEN
                ELEMATPSI(COUNT)=ELEWIC
                DO ILOC=1,NLOC! Was loop 
                   ELEMATWEI((COUNT-1)*NLOC+ILOC)=LOCCORDSK(ILOC)
                END DO
             ENDIF
          ENDIF

       END DO ! Was loop 20
    END DO ! Was loop 10

    RETURN

  end subroutine finptsstore
!     
!     
!     
!     
      SUBROUTINE MATPTSSTORE(MATPSI,COUNT,NFIELD,NOD,&
     &     PSI,NONODS,X_NONODS,&
     &     NLOC,TOTELE,X_NDGLN,NDGLNO,&
     &     FINDRM,COLM,NCOLM,&
     &     X1,Y1,Z1,&
     &     X2,Y2,Z2,&
     &     NORMX1,NORMY1,NORMZ1,&
     &     X,Y,Z,&
!     work space...
     &     FINDELE,COLELE,NCOLEL,&
     &     MINPSI,MAXPSI,  &
     &     ELEWIC,LOCCORDSK,BOUND,REFLECT,NDIM)
!     This sub calculates the value of PSI that would be at the 
!     other side of the stencil if we had a linear variation and within 
!     a single element.     
! IF BOUND then make locally bounded.
      IMPLICIT NONE
      REAL INFINY,FRALINE
      LOGICAL REFLECT
! IF REFLECT then use a reflection condition at boundary to 
! do limiting. 
      PARAMETER(INFINY=1.E+20,FRALINE=0.001)
      LOGICAL BOUND
      INTEGER COUNT,NFIELD,NOD,NONODS,X_NONODS,NLOC,TOTELE,NDIM
      REAL MATPSI(NCOLM*NFIELD),PSI(NONODS*NFIELD)
      INTEGER X_NDGLN(NLOC*TOTELE),NDGLNO(NLOC*TOTELE)
      INTEGER NCOLM,FINDRM(NONODS+1),COLM(NCOLM)
      REAL X1,Y1,Z1,X2,Y2,Z2,NORMX1,NORMY1,NORMZ1
      REAL X(X_NONODS),Y(X_NONODS),Z(X_NONODS)
      INTEGER NCOLEL
      INTEGER FINDELE(NONODS+1),COLELE(NCOLEL)
      REAL MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
      INTEGER ELEWIC
      REAL LOCCORDSK(NLOC)
!     
!     Local variables...
      REAL XC,YC,ZC
      REAL LOCCORDS(4)
      INTEGER LOCNODS(4),LOCNODSK(4)
      INTEGER NLOCNODS(4),NLOCNODSK(4)
      INTEGER ELE,ILOC,KNOD,JNOD,IFIELD, COUNT2
      REAL MINCOR,MINCORK,SUM
      REAL VX,VY,VZ,T2X,T2Y,T2Z,T1X,T1Y,T1Z,DIST12,RN,RMATPSI
      REAL REFX,REFY,REFZ,REFX2,REFY2,REFZ2
!     
      XC=X1 - FRALINE*(X2-X1)
      YC=Y1 - FRALINE*(Y2-Y1)
      ZC=Z1 - FRALINE*(Z2-Z1)
      
      IF(REFLECT) THEN
         IF(ABS(NORMX1)+ABS(NORMY1)+ABS(NORMZ1).NE.0.0) THEN
!  if (XC,YC,ZC) is outside the domain
!     The rotation matrix in 3-D is R=  
!     NX    NY    NZ
!     T1X   T1Y   T1Z
!     T2X   T2Y   T2Z
!     
            VX=X1-X2
            VY=Y1-Y2
            VZ=Z1-Z2
!     
            CALL XPROD(T2X,T2Y,T2Z, NORMX1,NORMY1,NORMZ1, VX,VY,VZ)
!     
            DIST12=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
            RN=SQRT(T2X**2+T2Y**2+T2Z**2)
            IF(RN.LT.(1.E-5)*DIST12) THEN
!     Simply have VX,VY,VZ going in the opposite direction...
               XC=X1 - VX*FRALINE
               YC=Y1 - VY*FRALINE
               ZC=Z1 - VZ*FRALINE
            ELSE
               T2X= T2X/RN
               T2Y= T2Y/RN
               T2Z= T2Z/RN
!     T1=Nx (-T2)
               CALL XPROD(T1X,T1Y,T1Z, NORMX1,NORMY1,NORMZ1, -T2X,-T2Y,-T2Z)
!     
               REFX2= NORMX1*VX + NORMY1*VY + NORMZ1*VZ
               REFY2= T1X   *VX + T1Y   *VY + T1Z   *VZ
               REFZ2= T2X   *VX + T2Y   *VY + T2Z   *VZ
!     Reflect...
               REFX2=-REFX2
!     MAP BACK USING R^T
               REFX = NORMX1*REFX2 + T1X*REFY2 + T2X*REFZ2
               REFY = NORMY1*REFX2 + T1Y*REFY2 + T2Y*REFZ2
               REFZ = NORMZ1*REFX2 + T1Z*REFY2 + T2Z*REFZ2
!     
!     (REFX,REFY,REFZ) is the reflected direction...
               XC=X1 + REFX*FRALINE
               YC=Y1 + REFY*FRALINE
               ZC=Z1 + REFZ*FRALINE
            ENDIF

         ENDIF
      ENDIF
!     
      MINCORK=-INFINY
!
      DO COUNT2=FINDELE(NOD),FINDELE(NOD+1)-1! Was loop 10
         ELE=COLELE(COUNT2)
!     
         NLOCNODS(1)=NDGLNO((ELE-1)*NLOC+1)
         NLOCNODS(2)=NDGLNO((ELE-1)*NLOC+2)
         NLOCNODS(3)=NDGLNO((ELE-1)*NLOC+3)
         if(ndim==3) NLOCNODS(4)=NDGLNO((ELE-1)*NLOC+4)
!     
         LOCNODS(1)=X_NDGLN((ELE-1)*NLOC+1)
         LOCNODS(2)=X_NDGLN((ELE-1)*NLOC+2)
         LOCNODS(3)=X_NDGLN((ELE-1)*NLOC+3)
         if(ndim==3) LOCNODS(4)=X_NDGLN((ELE-1)*NLOC+4)
!     
! Calculate the local coord but with 4th point replaced by INOD...
! Find local coords LOCCORDS of point INOD corresponding to these nodes LOCNODS...
 
         IF (NDIM==3) THEN
         CALL TRILOCCORDS(XC,YC,ZC, &
     &        LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),LOCCORDS(4),&
!     The 4 corners of the tet...
     &        X(LOCNODS(1)),Y(LOCNODS(1)),Z(LOCNODS(1)),&
     &        X(LOCNODS(2)),Y(LOCNODS(2)),Z(LOCNODS(2)),&
     &        X(LOCNODS(3)),Y(LOCNODS(3)),Z(LOCNODS(3)),&
     &        X(LOCNODS(4)),Y(LOCNODS(4)),Z(LOCNODS(4))  )
         ELSE
        CALL TRILOCCORDS2D(XC,YC, &
     &        LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),&
!     The 3 corners of the tri...
     &        X(LOCNODS(1)),Y(LOCNODS(1)),&
     &        X(LOCNODS(2)),Y(LOCNODS(2)),&
     &        X(LOCNODS(3)),Y(LOCNODS(3))  )

         END IF

!     
!         MINCOR=MIN(LOCCORDS(1),LOCCORDS(2), LOCCORDS(3),LOCCORDS(4))
         MINCOR=MINVAL( LOCCORDS(1:NLOC) )

         IF(MINCOR.GT.MINCORK) THEN 
            MINCORK=MINCOR
            DO ILOC=1,NLOC! Was loop 
               LOCCORDSK(ILOC)=LOCCORDS(ILOC)
               LOCNODSK(ILOC)=LOCNODS(ILOC)
               NLOCNODSK(ILOC)=NLOCNODS(ILOC)
            END DO
            ELEWIC=ELE
         ENDIF
      END DO ! Was loop 10

!     Set all the negative basis to zero and re-normalise 
!     to put on the face of an element...
      SUM=0.0
      DO ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=MAX(0.0,LOCCORDSK(ILOC))
         SUM=SUM+LOCCORDSK(ILOC)
      END DO
      DO ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=LOCCORDSK(ILOC)/SUM
      END DO
      DO IFIELD=1,NFIELD
         RMATPSI=0.0
         DO ILOC=1,NLOC! Was loop 
            RMATPSI   =RMATPSI  +LOCCORDSK(ILOC)*PSI(NLOCNODSK(ILOC)+(IFIELD-1)*NONODS)
!         XC=XC+LOCCORDSK(ILOC)*X(LOCNODSK(ILOC))
!         YC=YC+LOCCORDSK(ILOC)*Y(LOCNODSK(ILOC))
!         ZC=ZC+LOCCORDSK(ILOC)*Z(LOCNODSK(ILOC))
         END DO
!     Exaduate difference by a factor of 100.
         RMATPSI   =PSI(NOD+(IFIELD-1)*NONODS)  &
         +(1./FRALINE)*(RMATPSI   -PSI(NOD+(IFIELD-1)*NONODS))

!     Now correct to make sure that we get a bounded soln...
         IF(BOUND) THEN
           RMATPSI   =MAX(MIN(RMATPSI,   MAXPSI(ELEWIC+(IFIELD-1)*TOTELE)),   MINPSI(ELEWIC+(IFIELD-1)*TOTELE))
         ENDIF
         MATPSI(COUNT+(IFIELD-1)*NCOLM)   =RMATPSI
      END DO
!
      RETURN

      end subroutine matptsstore
!     
!     
!     
!     
      SUBROUTINE PHILNODELE(NONODS,FINDELE,COLELE, &
     &     NCOLEL,MXNCOLEL, &
     &     TOTELE,NLOC,NDGLNO, &
     &     NLIST,INLIST)
      !=================================================================
      ! This sub calculates the node to element list FINDELE,COLELE
      ! 
      ! Note NLIST and INLIST are only used locally but are passed 
      ! down from parent routine where they are dynamically allocated.
      !
      ! INPUTS:
      ! ------
      ! NDGLNO  - List of global node numbers
      !
      ! OUTPUTS: 
      ! -------
      ! COLELE  - This is a list of the element numbers that each node
      !           belongs to.  So it lists all elements for node 1, then
      !           all elements for node 2, and so on...
      ! FINDELE - is the pointer to the place in COLELE that gives the 
      !           first element associated with a given global node
      ! 
      ! Called from subroutines IFINPTS and FINPTS, which are
      ! subroutines of CONSTRUCT_ADVECTION_DIFFUSION_CV
      ! 
      ! Description                                   Programmer      Date
      ! ==================================================================
      ! Original version..................................CCP   2013-28-01
      !
      !================================================================ 
      IMPLICIT NONE
      integer, intent( in ) :: NONODS,MXNCOLEL,TOTELE,NLOC
      integer, intent( inout ) :: NCOLEL
      integer, dimension( NONODS+1 ), intent( inout ) :: FINDELE
      integer, dimension( MXNCOLEL ), intent( inout ) :: COLELE
      integer, dimension( TOTELE*NLOC ), intent( in ) :: NDGLNO
      integer, dimension( NONODS ), intent( inout ) :: NLIST,INLIST
!     Local variables...
      INTEGER NOD,ELE,ILOC,COUNT, INOD
!     
      do NOD=1,NONODS! Was loop 
         NLIST(NOD)=0
         INLIST(NOD)=0
      END DO

      ! NLIST is the number of elements each node belongs to...
      do ELE=1,TOTELE! Was loop 
      do ILOC=1,NLOC! Was loop 
!          print *,'iloc,nloc,totele,ele:', iloc,nloc,totele,ele
!          print *,'NDGLNO((ELE-1)*NLOC+ILOC):',NDGLNO((ELE-1)*NLOC+ILOC)
            INOD=NDGLNO((ELE-1)*NLOC+ILOC)
            NLIST(INOD)=NLIST(INOD)+1
         END DO
      END DO

      ! FINDELE is a pointer to the first element 
      ! associated with a given global node (NOD)
      COUNT=0
      do NOD=1,NONODS! Was loop 
         FINDELE(NOD)=COUNT+1
         COUNT=COUNT+NLIST(NOD)
      END DO
      FINDELE(NONODS+1)=COUNT+1
      NCOLEL=COUNT

      ! COLELE is a list of the element numbers each node belongs
      ! to stored in the order of the global nodes...
      ! INLIST is the element number the node belongs to.
      do ELE=1,TOTELE! Was loop 
      do ILOC=1,NLOC! Was loop 
            INOD=NDGLNO((ELE-1)*NLOC+ILOC)
            INLIST(INOD)=INLIST(INOD)+1
            IF (FINDELE(INOD)-1+INLIST(INOD).GT.MXNCOLEL) THEN
               STOP 'COLELE ARRAY OUT OF BOUNDS--SUB:PHILNODELE'
            ENDIF
            COLELE(FINDELE(INOD)-1+INLIST(INOD))=ELE ! 
         END DO
      END DO
      RETURN

  end subroutine philnodele
!     
!     
!     
!     
      Subroutine TRILOCCORDS(Xp,Yp,Zp, &
     &     N1, N2, N3, N4, &
     &     X1,Y1,Z1, &
     &     X2,Y2,Z2, &
     &     X3,Y3,Z3, &
     &     X4,Y4,Z4  )

      IMPLICIT NONE
      Real Xp, Yp, Zp

      Real N1, N2, N3, N4
      
      Real X1,Y1,Z1
      Real X2,Y2,Z2
      Real X3,Y3,Z3
      Real X4,Y4,Z4

      Real Volume

!     calculate element volume...

      Volume = TetVolume(X1, Y1, Z1, &
     &     X2, Y2, Z2, &
     &     X3, Y3, Z3, &
     &     X4, Y4, Z4)
      
      Volume = Volume /6.0


!     vol coords...

      N1 = TetVolume(Xp, Yp, Zp, &
     &     X2, Y2, Z2, &
     &     X3, Y3, Z3, &
     &     X4, Y4, Z4) 

      N1 = N1/(6.0*Volume)

      

      N2 = TetVolume(X1, Y1, Z1, &
     &     Xp, Yp, Zp, &
     &     X3, Y3, Z3, &
     &     X4, Y4, Z4) 
      
      N2 = N2/(6.0*Volume)
      


      N3 = TetVolume(X1, Y1, Z1, &
     &     X2, Y2, Z2, &
     &     Xp, Yp, Zp, &
     &     X4, Y4, Z4) 

      N3 = N3/(6.0*Volume)

      
      N4 = TetVolume(X1, Y1, Z1, &
     &     X2, Y2, Z2, &
     &     X3, Y3, Z3, &
     &     Xp, Yp, Zp) 

      N4 = N4/(6.0*Volume)


      Return

  end subroutine triloccords


pure function tetvolume(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3)
      IMPLICIT NONE
   
     real, intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3
     
     real :: tetvolume
      
     ! tetvolume = 1.0 / 6.0 * det |three tet edge vectors|
     ! Chris' tets have a clockwise base, hence the sign change in the det
     tetvolume = &
       (  &
         & - (x1 - x0) * ((y2 - y0) * (z3 - z0) - (y3 - y0) * (z2 - z0)) &
         & + (y1 - y0) * ((x2 - x0) * (z3 - z0) - (x3 - x0) * (z2 - z0)) &
         & - (z1 - z0) * ((x2 - x0) * (y3 - y0) - (x3 - x0) * (y2 - y0)) &
       & ) / 6.0
       
   end function tetvolume

 
!
! 
! 
!     
      Subroutine TRILOCCORDS2D(Xp,Yp, &
     &     N1, N2, N3,  &
     &     X1,Y1, &
     &     X2,Y2, &
     &     X3,Y3 )

      IMPLICIT NONE
      Real Xp,Yp, &
     &     N1, N2, N3,  &
     &     X1,Y1, &
     &     X2,Y2, &
     &     X3,Y3 

      Real AREA

      AREA = TRIAREAF( X1, Y1, X2, Y2, X3, Y3)
!     area coords...

      N1 = TRIAREAF(Xp, Yp,  &
     &     X2, Y2,  &
     &     X3, Y3 ) 

      N1 = N1/AREA

      

      N2 = TRIAREAF(X1, Y1, &
     &     Xp, Yp,  &
     &     X3, Y3 ) 
      
      N2 = N2/AREA
      


      N3 = TRIAREAF(X1, Y1,  &
     &     X2, Y2,  &
     &     Xp, Yp ) 

      N3 = N3/AREA
       

      Return

      end subroutine triloccords2d
! 
! 
! 
! 


!
!
! 
! 

! -----------------------------------------------------------------------------

  end module cv_advection
      
