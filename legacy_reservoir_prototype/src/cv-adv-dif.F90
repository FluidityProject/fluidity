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
  use global_parameters, only: option_path_len, timestep
  use futils, only: int2str
  use adapt_state_prescribed_module

  use shape_functions
  use matrix_operations

contains

    SUBROUTINE CV_ASSEMB( state, &
         CV_RHS, &
         NCOLACV, CSR_ACV, dense_acv, FINACV, COLACV, MIDACV, &
         SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM,&
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
         SUF_T_BC_ROB1, SUF_T_BC_ROB2, &
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
         option_path_spatial_discretisation,&
         StorageIndexes )

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
      !     CSR_ACV   - Matrix for assembling the advection terms (empty on input)
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
      !     CSR_ACV   - Matrix updated to include the advection terms
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
      type( state_type ), dimension( : ), intent( inout ) :: state
      INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_DISOPT, CV_DG_VEL_INT_OPT, NDIM, &
           NCOLM, XU_NLOC, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NCOLCMC
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(: ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T2_BC
      REAL, DIMENSION( : ), intent( inout ) :: CV_RHS
      !REAL, DIMENSION( NCOLACV ), intent( inout ) :: ACV
      real, dimension(:), intent(inout) :: CSR_ACV
      real, dimension(:,:,:), intent(inout) :: dense_acv
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
      REAL, DIMENSION( : ), intent( inout ) :: CT
      ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
      REAL, DIMENSION( : ), intent( inout ) :: DIAG_SCALE_PRES
      REAL, DIMENSION( :  ), intent( inout ) :: CT_RHS
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES

      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ), intent( in ) :: U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( : ), intent( in ) :: T, TOLD, DEN, DENOLD
      REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( : ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION( :, :, :, : ), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
      REAL, intent( in ) :: DT, CV_THETA, SECOND_THETA, CV_BETA
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION(: ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: CV_P
      REAL, DIMENSION( : ), intent( in ) :: SOURCT
      REAL, DIMENSION( :, :, : ), intent( in ) :: ABSORBT
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      LOGICAL, intent( in ) :: GETCV_DISC, GETCT, GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION(: ), intent( inout ) :: T_FEMT, DEN_FEMT
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
      REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE_TRANSP
      character( len = * ), intent( in ), optional :: option_path_spatial_discretisation
      integer, dimension(:), intent(in) :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      integer, dimension(:), intent(inout) :: StorageIndexes
      !character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation


      ! Local variables
      INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
           WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
           WIC_U_BC_DIRICHLET = 1
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE
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
           CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
           FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
           SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
           SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
           T2MAX, T2MIN, T2OLDMAX, &
           T2OLDMIN, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, &
           T2OLDMIN_2ND_MC, &
           UP_WIND_NOD, DU, DV, DW, PERM_ELE
      REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
           UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ, &
           DUMMY_ZERO_NDIM_NDIM
      REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE,  &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
      REAL, pointer, DIMENSION( : , :, : ) :: INV_JAC
      real, pointer, dimension(:) :: SCVDETWEI, SCVRA
      real, pointer, dimension(:,:,:) :: SCVFENX_ALL
      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE,   &
           NCOLGPTS, &
           CV_SILOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, rhs_nodj_ipha,rhs_nodi_ipha,&
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
           W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW

      REAL, PARAMETER :: W_SUM_ONE = 1.
      real, pointer :: VOLUME
      integer :: cv_inod_ipha, IGETCT, U_NODK_IPHA, IANISOLIM
      logical :: Have_Temperature_Fields, Have_VolumeFraction_Fields, Have_Components_Fields
      logical :: overlapping
      ! Functions...
      !REAL :: R2NORM, FACE_THETA
      !        ===>  LOGICALS  <===
      character( len = option_path_len ) :: overlapping_path
      LOGICAL :: GETMAT, LIMIT_USE_2ND, &
           D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
           NORMALISE, SUM2ONE, GET_GTHETA, QUAD_OVER_WHOLE_ELE,is_overlapping

      character( len = option_path_len ) :: option_path, option_path2, path_temp, path_volf, &
           path_comp, path_spatial_discretisation


      real, dimension(:), allocatable :: TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, &
              DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT
      INTEGER :: IDUM(1)
      REAL :: RDUM(1),n1,n2,n3
      type( scalar_field ), pointer :: perm
     !Reals to store the irresidual water and irreducible oil values, used in GET_INT_T_DEN
    real ::  s_gc, s_or

      overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) overlapping = .true.

!      ALLOCATE( PERM_ELE( TOTELE ) ) ; PERM_ELE = 0.0
!      if ( overlapping ) then
!         perm => extract_scalar_field( state(1), "Permeability" )
!         perm_ele = perm % val
!      end if

      IDUM = 0
      RDUM = 0.

!       call TRILOCCORDS2D(Xp,Yp, &
!     &     N1, N2, N3,  &
!     &     X1,Y1, &
!     &     X2,Y2, &
!     &     X3,Y3 )

!       call TRILOCCORDS2D(1.,1., &
!     &     N1, N2, N3,  &
!     &     0.,0., &
!     &     1.0,0.0, &
!     &     0.0,1.0 )

!        print *,'n1,n2,n3:',n1,n2,n3
!        stop 722

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

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.


      GOT_DIFFUS = ( R2NORM( TDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0 )

      ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS:', &
           CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS
      ewrite(3,*)'GETCV_DISC, GETCT', GETCV_DISC, GETCT

      !ewrite(3,*)'tdiffusion=', tdiffusion
      !ewrite(3,*)'suf_t_bc=', suf_t_bc
      !ewrite(3,*)'wic_t_bc=', wic_t_bc
      !ewrite(3,*)'nu=', nu
      !ewrite(3,*)'nv=', nv
      !ewrite(3,*)'nw=', w
      !ewrite(3,*)'den=', den
      !ewrite(3,*)'denold=', denold

      ndotq = 0. ; ndotqold = 0.

      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

      ! Allocate memory for the control volume surface shape functions, etc.
      ALLOCATE( JCOUNT_KLOC( U_NLOC )) ; jcount_kloc = 0
      ALLOCATE( JCOUNT_KLOC2( U_NLOC )) ; jcount_kloc2 = 0

      ALLOCATE( CVNORMX( SCVNGI ))
      ALLOCATE( CVNORMY( SCVNGI ))
      ALLOCATE( CVNORMZ( SCVNGI ))
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
      ALLOCATE( SCVFEWEIGH( SCVNGI ))

      ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

      ALLOCATE( SRA( SCVNGI ))

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


      UP_WIND_NOD = 0.0

      ALLOCATE( ONE_PORE( TOTELE ))
      IF( have_option( '/porous_media/actual_velocity' ) ) THEN
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
           T, TOLD, T2, T2OLD, DEN, DENOLD, IGOT_T2, NPHASE, CV_NONODS, size(small_colm), SMALL_FINDRM, SMALL_COLM, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC, SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
           WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
           WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
           T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
           DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD )

! **********ANISOTROPIC LIMITING...*******************
      IANISOLIM=0
     ! print *,'CV_DISOPT=',CV_DISOPT
     ! stop 72
      IF(CV_DISOPT.GE.5) IANISOLIM=1
! limiting not ready yet for P2 DG:
!      IF(
!      IF( (NDIM==2).AND.(CV_NLOC==6) ) IANISOLIM=0
!      IF( (NDIM==3).AND.(CV_NLOC==10) ) IANISOLIM=0
!IANISOLIM=0
      IF (IANISOLIM==0) THEN
         ALLOCATE(TUPWIND_MAT(1), TOLDUPWIND_MAT(1), DENUPWIND_MAT(1), DENOLDUPWIND_MAT(1))
         ALLOCATE(T2UPWIND_MAT(1), T2OLDUPWIND_MAT(1))
         NSMALL_COLM=1
      ELSE ! IF (IANISOLIM==1) THEN
         ! Reduce matrix size...
         NSMALL_COLM=size(small_colm)

!        !We need to allocate CSR_ADV
!        allocate(CSR_ACV(NPHASE*NSMALL_COLM))

         !ALLOCATE(SMALL_FINDRM(CV_NONODS+1),SMALL_COLM(NSMALL_COLM),SMALL_CENTRM(CV_NONODS))
         ALLOCATE(TUPWIND_MAT(NSMALL_COLM*NPHASE), TOLDUPWIND_MAT(NSMALL_COLM*NPHASE), &
              DENUPWIND_MAT(NSMALL_COLM*NPHASE), DENOLDUPWIND_MAT(NSMALL_COLM*NPHASE))
         ALLOCATE(T2UPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2), T2OLDUPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2))

         CALL CALC_ANISOTROP_LIM(&
              ! Caculate the upwind values stored in matrix form...
              T,TOLD,DEN,DENOLD,T2,T2OLD, &
              FEMT,FEMTOLD,FEMDEN,FEMDENOLD,FEMT2,FEMT2OLD, (CV_NONODS.NE.X_NONODS), &
              TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, &
              T2UPWIND_MAT, T2OLDUPWIND_MAT, &
              ! Store the upwind element for interpolation and its weights for
              ! faster results...
              IGOT_T2,NPHASE,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
              SMALL_FINDRM,SMALL_CENTRM,SMALL_COLM,NSMALL_COLM, &
              X_NDGLN,X_NONODS,NDIM, &
              X,Y,Z, XC_CV, YC_CV, ZC_CV)

         if ( .false. ) then
            ewrite(3,*) 'TUPWIND_MAT', TUPWIND_MAT(1:nsmall_colm)
            ewrite(3,*) 'TOLDUPWIND_MAT', TOLDUPWIND_MAT(1:nsmall_colm)
            if (igot_t2==1)then
               ewrite(3,*) 'T2UPWIND_MAT', T2UPWIND_MAT(1:nsmall_colm)
               ewrite(3,*) 'T2OLDUPWIND_MAT', T2OLDUPWIND_MAT(1:nsmall_colm)
            end if
         end if

      END IF
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
              SBCVFEN, SBCVFENSLX, SBCVFENSLY&
              ,  state, "CV", StorageIndexes(16))
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
            CSR_ACV = 0.0
         ENDIF
      ENDIF


      GET_GTHETA=.FALSE.
      IF( IGOT_THETA_FLUX == 1 ) THEN
         IF( GET_THETA_FLUX ) THEN
            THETA_FLUX = 0.0
            ONE_M_THETA_FLUX = 0.0
            IF( IGOT_T2 == 1 ) THEN
               GET_GTHETA = .TRUE.
               THETA_GDIFF = 0.0
            END IF
         ENDIF
      ENDIF

     !Get the irreducible water and residual oil before entering the loop
     !this value are used in GET_INT_T_DEN to keep the phases between realistic values
     !Default value has to be the same as in subroutine get_corey_options in multi_eos.F90
      !Default value of  S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
          s_gc, default=0.1)
      !Default value of    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
           s_or, default=0.3)


      ! Now we begin the loop over elements to assemble the advection terms
      ! into the matrix (ACV) and the RHS
      Loop_Elements: DO ELE = 1, TOTELE


         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR_INVJAC( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
              CV_NLOC, SCVNGI, &
              SCVFEN, SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH, SCVDETWEI, SCVRA, VOLUME, D1, D3, DCYL, &
              SCVFENX_ALL, &
              NDIM, INV_JAC, state, "CVI", StorageIndexes(29))

         Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC ! Loop over the nodes of the element

            ! Global node number of the local node
            CV_NODI = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )

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
                       X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                       CV_SNLOC, CVFEM_ON_FACE(:,GI), X_SHARE, X_NONODS, ELE, ELE2,  &
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

                  IF ( INTEGRAT_AT_GI ) THEN
                     CV_JLOC = CV_OTHER_LOC( CV_ILOC )
                     SELE = 0

                     IF ( CV_JLOC == 0 ) THEN ! We are on the boundary of the domain
                        CV_JLOC = CV_ILOC
                        ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC
                        CALL CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                             FACE_ELE, NFACE, CVFEM_ON_FACE( :, GI ), &
                             CV_NONODS, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
                             CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )
                        !EWRITE(3,*)'*****AFTER CALC_SELE SELE,CV_SILOC,CV_SNLOC:',SELE,CV_SILOC,CV_SNLOC
                     END IF
                     INTEGRAT_AT_GI = .NOT.( (ELE==ELE2) .AND. (SELE==0) )
                  END IF

               END IF Conditional_CheckingNeighbourhood

               ! avoid indegrating across the middle of a CV on the boundaries of elements
               Conditional_integration: IF ( INTEGRAT_AT_GI ) THEN

                  ! if necessary determine the derivatives between elements ELE and ELE2

                  ! Calculate the control volume normals at the Gauss pts.
                  CALL SCVDETNX( ELE, GI, &
                       X_NLOC, SCVNGI, TOTELE, &
                       X_NDGLN, X_NONODS, &
                       SCVDETWEI, CVNORMX, CVNORMY, &
                       CVNORMZ, SCVFEN, SCVFENSLX, &
                       SCVFENSLY, SCVFEWEIGH, XC_CV(CV_NODI), &
                       YC_CV(CV_NODI), ZC_CV(CV_NODI), &
                       X, Y, Z, &
                       D1, D3, DCYL )

                  ! ================ COMPUTE THE FLUX ACROSS SUB-CV FACE ===============

                  ! Find its global node number
                  IF ( ELE2 == 0 ) THEN
                     CV_NODJ = CV_NDGLN( ( ELE - 1 )  * CV_NLOC + CV_JLOC )
                  ELSE
                     CV_NODJ = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_JLOC )
                  END IF
                  X_NODI = X_NDGLN( ( ELE - 1 ) * X_NLOC  + CV_ILOC )

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
                           IF(COLCT( COUNT ) == U_NODK) then
                              JCOUNT = COUNT
                              exit
                           end if
                        END DO
                        JCOUNT_KLOC( U_KLOC ) = JCOUNT
                        !ewrite(3,*)' u_nodk, jcount1:', u_nodk, jcount
                     END DO

                     IF( ( ELE2 /= 0 ) .AND. ( ELE2 /= ELE ) ) THEN
                        DO U_KLOC = 1, U_NLOC
                           U_NODK = U_NDGLN(( ELE2 - 1 ) * U_NLOC + U_KLOC )
                           JCOUNT = 0
                           DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1, 1
                              IF(COLCT( COUNT ) == U_NODK) then
                                 JCOUNT = COUNT
                                 exit
                              end if
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
                        rhs_NODI_IPHA = IPHASE + (cv_nodi-1 ) * nphase
                        rhs_NODJ_IPHA = IPHASE + (cv_nodj-1 ) * nphase

                        ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
                        IF(IGOT_T2==1) THEN
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW,  NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2OLD, FEMT2OLD, DENOLD, &
                                U, V, W, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT2OLD,  FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                                T2OLDMIN, T2OLDMAX, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2OLDMIN_2ND_MC, T2OLDMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2OLDUPWIND_MAT )
                            CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2, FEMT2, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT2, FEMDGI, FEMT2GI,  UP_WIND_NOD, &
                                T2MIN, T2MAX,  T2MIN_NOD, T2MAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC, T2MAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2UPWIND_MAT )
                        ELSE
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                TOLD, FEMTOLD, DENOLD, &
                                U, V, W, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1,  LIMTOLD, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
                                TOLDMIN, TOLDMAX, TOLDMIN_NOD, TOLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TOLDMIN_2ND_MC, TOLDMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TOLDUPWIND_MAT )
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T, FEMT, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
                                TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TMIN_2ND_MC, TMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TUPWIND_MAT )
                        ENDIF


                        !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                        ! Calculate T and DEN on the CV face at quadrature point GI.
                        CALL GET_INT_T_DEN( FVTOLD, FVT2OLD, FVDOLD, &
                             LIMDOLD, LIMTOLD, LIMT2OLD, LIMDTOLD, LIMDTT2OLD,&
                             FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOMEOLD, &
                             TOLD, T2OLD, DENOLD, FEMTOLD, FEMT2OLD, FEMDENOLD, &
                             TOLDMIN, T2OLDMIN, DENOLDMIN, &
                             TOLDMAX, T2OLDMAX, DENOLDMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TOLDMIN_NOD, TOLDMAX_NOD, &
                             DENOLDMIN_NOD, DENOLDMAX_NOD, &
                             T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2, &
                             TOLDMIN_2ND_MC, T2OLDMIN_2ND_MC, DENOLDMIN_2ND_MC, &
                             TOLDMAX_2ND_MC, T2OLDMAX_2ND_MC, DENOLDMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQOLD, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TOLDUPWIND_MAT, DENOLDUPWIND_MAT, T2OLDUPWIND_MAT , &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )
                        
                        CALL GET_INT_T_DEN( FVT, FVT2, FVD, LIMD, LIMT, LIMT2, &
                             LIMDT, LIMDTT2,&
                             FEMDGI, FEMTGI,FEMT2GI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, &
                             T, T2, DEN, FEMT, FEMT2, FEMDEN, &
                             TMIN, T2MIN, DENMIN, &
                             TMAX, T2MAX, DENMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TMIN_NOD, TMAX_NOD, &
                             DENMIN_NOD, DENMAX_NOD, &
                             T2MIN_NOD, T2MAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, T2MIN_2ND_MC, DENMIN_2ND_MC, &
                             TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT, &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )

                        SUM_LIMT    = SUM_LIMT    + LIMT
                        SUM_LIMTOLD = SUM_LIMTOLD + LIMTOLD

                     END DO Loop_IPHASE5
                  ENDIF Conditional_SUMLimiting
                  ! get the sum of limiting functions correct...************

                  DO COUNT = SMALL_FINDRM( CV_NODI ), SMALL_FINDRM( CV_NODI + 1 ) - 1
                     IF( SMALL_COLM( COUNT ) == CV_NODJ ) THEN 
                        JCOUNT_IPHA = COUNT
                        EXIT
                     END IF!An exit may improve the performance!!!
                  END DO

                  Loop_IPHASE: DO IPHASE = 1, NPHASE

                     CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                     CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS
                     rhs_NODI_IPHA = IPHASE + (cv_nodi-1 ) * nphase
                     rhs_NODJ_IPHA = IPHASE + (cv_nodj-1 ) * nphase

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
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW,  NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2OLD, FEMT2OLD, DENOLD, &
                                U, V, W, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT2OLD, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                                T2OLDMIN, T2OLDMAX, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2OLDMIN_2ND_MC, T2OLDMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2OLDUPWIND_MAT )
                            CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2, FEMT2, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT2, FEMDGI, FEMT2GI,  UP_WIND_NOD, &
                                T2MIN, T2MAX,  T2MIN_NOD, T2MAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC, T2MAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2UPWIND_MAT )
                        ELSE
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                TOLD, FEMTOLD, DENOLD, &
                                U, V, W, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1,  LIMTOLD, FEMDOLdGI, FEMTOLDGI, UP_WIND_NOD, &
                                TOLDMIN, TOLDMAX, TOLDMIN_NOD, TOLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TOLDMIN_2ND_MC, TOLDMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TOLDUPWIND_MAT )
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T, FEMT, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
                                TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TMIN_2ND_MC, TMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TUPWIND_MAT )
                        ENDIF

                        !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                        ! Calculate T and DEN on the CV face at quadrature point GI.
                        CALL GET_INT_T_DEN( FVTOLD, FVT2OLD, FVDOLD, &
                             LIMDOLD, LIMTOLD, LIMT2OLD, LIMDTOLD, LIMDTT2OLD,&
                             FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOMEOLD, &
                             TOLD, T2OLD, DENOLD, FEMTOLD, FEMT2OLD, FEMDENOLD, &
                             TOLDMIN, T2OLDMIN, DENOLDMIN, &
                             TOLDMAX, T2OLDMAX, DENOLDMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TOLDMIN_NOD, TOLDMAX_NOD, &
                             DENOLDMIN_NOD, DENOLDMAX_NOD, &
                             T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2, &
                             TOLDMIN_2ND_MC, T2OLDMIN_2ND_MC, DENOLDMIN_2ND_MC, &
                             TOLDMAX_2ND_MC, T2OLDMAX_2ND_MC, DENOLDMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQOLD, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TOLDUPWIND_MAT, DENOLDUPWIND_MAT, T2OLDUPWIND_MAT , &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )
                        
                        CALL GET_INT_T_DEN( FVT, FVT2, FVD, LIMD, LIMT, LIMT2, &
                             LIMDT, LIMDTT2,&
                             FEMDGI, FEMTGI,FEMT2GI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, &
                             T, T2, DEN, FEMT, FEMT2, FEMDEN, &
                             TMIN, T2MIN, DENMIN, &
                             TMAX, T2MAX, DENMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TMIN_NOD, TMAX_NOD, &
                             DENMIN_NOD, DENMAX_NOD, &
                             T2MIN_NOD, T2MAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, T2MIN_2ND_MC, DENMIN_2ND_MC, &
                             TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT, &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )

                     END DO

                     IF( SUM2ONE ) THEN
                        LIMT = LIMT / SUM_LIMT
                        LIMTOLD = LIMTOLD / SUM_LIMTOLD
                        LIMDT = LIMT * LIMD
                        LIMDTOLD = LIMTOLD * LIMDOLD
                     END IF

                     ! Define face value of theta
                     IF(IGOT_T2==1) THEN
                        FTHETA=FACE_THETA( DT, CV_THETA, ( cv_disopt>=8 ),HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                             T(CV_NODJ_IPHA)*DEN(CV_NODJ_IPHA)*T2(CV_NODJ_IPHA), &
                             T(CV_NODI_IPHA)*DEN(CV_NODI_IPHA)*T2(CV_NODI_IPHA), &
                             NDOTQOLD, LIMDTT2OLD, DIFF_COEFOLD_DIVDX, &
                             TOLD(CV_NODJ_IPHA)*DENOLD(CV_NODJ_IPHA)*T2OLD(CV_NODJ_IPHA), &
                             TOLD(CV_NODI_IPHA)*DENOLD(CV_NODI_IPHA)*T2OLD(CV_NODI_IPHA) )
                     ELSE
                        FTHETA=FACE_THETA( DT, CV_THETA, ( cv_disopt>=8 ),HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
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

! ********FOR RELATIVE PERM B.C WHEN SPECIFYING PRESSURE*********
      IF( .false. ) THEN ! still bugs in the original better code but this is ok...
!      IF( is_overlapping ) THEN ! still bugs in the original better code but this is ok...

              !       RBC_OUT(1:NDIM)=1.0
              !       IF(SELE.NE.0) RBC_OUT(1:NDIM) &
              ! =SUF_SIG_DIAGTEN_BC(( SELE - 1 ) * CV_SNLOC + 1 +( IPHASE - 1 ) * STOTEL*CV_SNLOC,1:NDIM)

!                     IF((SELE.NE.0).AND.(NDOTQ.GE.0.0))  &
!                     IF((SELE.NE.0).AND.(NDOTQ.LE.0.0))  &
                     IF(SELE.NE.0)  &
                   SCVDETWEI( GI )=SCVDETWEI( GI ) &
               *SUF_SIG_DIAGTEN_BC(( SELE - 1 ) * CV_SNLOC + 1 +( IPHASE - 1 ) * STOTEL*CV_SNLOC,1)
      ENDIF
! ********FOR RELATIVE PERM B.C WHEN SPECIFYING PRESSURE*********

                     !====================== ACV AND RHS ASSEMBLY ===================
                     Conditional_GETCT2 : IF( GETCT ) THEN ! Obtain the CV discretised CT eqations plus RHS

                        CALL PUT_IN_CT_RHS(CT, CT_RHS, U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
                             CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE, ELE, ELE2, SELE, &
                             JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC, U_NDGLN, U, V, W, &
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
                           IF( ( CV_NODI_IPHA /= CV_NODJ_IPHA ) .AND. ( CV_NODJ /= 0 ) ) THEN
                              CSR_ACV( IPHASE+(JCOUNT_IPHA-1)*NPHASE ) =  CSR_ACV( IPHASE+(JCOUNT_IPHA-1)*NPHASE ) &
                                   + SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW * INCOME * LIMD  & ! advection
                                   - FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX ! Diffusion contribution

                              IF(GET_GTHETA) THEN
                                 THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                      + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX * T(CV_NODJ_IPHA) ! Diffusion contribution
                              ENDIF
                           ELSE IF(SELE/=0) THEN
                              IF(WIC_T_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRICHLET) THEN
                                 CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA )  &
                                      + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &
                                      * SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC+(IPHASE-1)*STOTEL*CV_SNLOC)
                                 IF(GET_GTHETA) THEN
                                    THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                         + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX &
                                         * SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC+(IPHASE-1)*STOTEL*CV_SNLOC)
                                 ENDIF
                              ENDIF
                           ENDIF

                           IMID_IPHA =  IPHASE+(SMALL_CENTRM(CV_NODI)-1)*NPHASE

                           CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA ) &
                                +  SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW * ( 1. - INCOME ) * LIMD & ! advection
                                +  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &  ! Diffusion contribution
                                +  SCVDETWEI( GI ) * ROBIN1  ! Robin bc

                           IF(GET_GTHETA) THEN
                              THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                                   -  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX * T( CV_NODI_IPHA ) & ! Diffusion contribution
                                   -  SCVDETWEI( GI ) * ROBIN1 * T( CV_NODI_IPHA )  ! Robin bc
                           ENDIF

                           ! CV_BETA=0 for Non-conservative discretisation (CV_BETA=1 for conservative disc)
                              CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA )  &
                                - SECOND_THETA * FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW * LIMD
                        ENDIF

                        TMID   =     T( CV_NODI_IPHA )
                        TOLDMID = TOLD( CV_NODI_IPHA )

                        ! Make allowances for no matrix stencil operating from outside the boundary.
                        BCZERO=1.0
                        IF( (SELE /= 0) .AND. (INCOME > 0.5) ) BCZERO=0.0

                        ! Put results into the RHS vector
                        CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA )  &
                                ! subtract 1st order adv. soln.
                             + SECOND_THETA * FTHETA_T2 * NDOTQNEW * SCVDETWEI( GI ) * LIMD * FVT * BCZERO &
                             -  SCVDETWEI( GI ) * ( FTHETA_T2 * NDOTQNEW * LIMDT &
                             + ONE_M_FTHETA_T2OLD * NDOTQOLD * LIMDTOLD ) ! hi order adv

                        ! Subtract out 1st order term non-conservative adv.
                        CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA ) &
                             - FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW * LIMD * TMID

                        ! High-order non-conservative advection contribution
                        CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA ) &
                             + ( 1. - CV_BETA) * SCVDETWEI( GI ) &
                             * ( FTHETA_T2 * NDOTQNEW * TMID * LIMD  &
                             + ONE_M_FTHETA_T2OLD * NDOTQOLD * LIMDOLD * TOLDMID )

                        ! Diffusion contribution
                        CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA ) &
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
                              CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                                   - CV_P( CV_NODI ) * SCVDETWEI( GI ) * ( &
                                   THERM_FTHETA * NDOTQNEW * LIMT2 &
                                   + ( 1. - THERM_FTHETA ) * NDOTQOLD * LIMT2OLD )
                           else
                              CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
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
         ewrite(3,*) 'cv_rhs1:', cv_rhs(1:cv_nonods)
         !ewrite(3,*) 'cv_rhs2:', cv_rhs(1+cv_nonods:2*cv_nonods)
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
         !ewrite(3,*)'before adding extra bits*****SOURCT:',SOURCT

         Loop_CVNODI2: DO CV_NODI = 1, CV_NONODS ! Put onto the diagonal of the matrix

            DO COUNT = SMALL_FINDRM( CV_NODI ), SMALL_FINDRM( CV_NODI + 1 ) - 1, 1
               IF( SMALL_COLM( COUNT ) == CV_NODI )  then
                  JCOUNT_IPHA = COUNT !An exit may improve the performance!!!
                  exit
               end IF
            END DO


            Loop_IPHASE2: DO IPHASE = 1, NPHASE
               CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
               rhs_NODI_IPHA = iphase + ( cv_nodi -1 ) * nphase
               ! For the gravity term
               !SOURCT2( CV_NODI_IPHA ) = SOURCT( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA )
               !SOURCT2( CV_NODI_IPHA ) = SOURCT2( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA )
!               IMID_IPHA = MIDACV( CV_NODI_IPHA )
                IMID_IPHA = IPHASE+ (SMALL_CENTRM(CV_NODI)-1)*NPHASE


               IF(IGOT_T2==1) THEN
                  CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                       + MASS_CV(CV_NODI) * SOURCT(CV_NODI_IPHA) !&
                  !- CV_BETA * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) & ! conservative time term.
                  !* TOLD(CV_NODI_IPHA)  &
                  !* (DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA) - DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA)) / DT


                  CSR_ACV( IMID_IPHA ) = CSR_ACV( IMID_IPHA ) &
                       !+ (CV_BETA * DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA) &
                       + (CV_BETA * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA))  &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) / DT


                  CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                       + (CV_BETA * DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA))  &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) * TOLD(CV_NODI_IPHA) / DT
               ELSE

                  CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                       + MASS_CV(CV_NODI) * SOURCT(CV_NODI_IPHA) !&
                  !- CV_BETA * MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) & ! conservative time term.
                  !* TOLD(CV_NODI_IPHA) &
                  !* (DEN(CV_NODI_IPHA) - DENOLD(CV_NODI_IPHA)) / DT

                  CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA ) &
                       !+ (CV_BETA * DENOLD(CV_NODI_IPHA) + (1.-CV_BETA) * DEN(CV_NODI_IPHA))  &
                       + (CV_BETA * DEN(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA))  &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) / DT

                  CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                       + (CV_BETA * DENOLD(CV_NODI_IPHA) &
                       + (1.-CV_BETA) * DEN(CV_NODI_IPHA)) &
                       * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) * TOLD(CV_NODI_IPHA) / DT
               ENDIF


               Conditional_GETMAT2: IF( GETMAT ) THEN

                  DO JPHASE = 1, NPHASE
                     !CV_NODI_JPHA = CV_NODI + ( JPHASE - 1 ) * CV_NONODS
!here we should put the full_ACV
                     dense_acv(JPHASE,IPHASE,CV_NODI)  = dense_acv(JPHASE,IPHASE,CV_NODI)  &
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
               rhs_NODI_IPHA = IPHASE + (cv_nodi-1 ) * nphase

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
!                  print *,'den=',den
!                  stop 578

                  W_SUM_ONE1 = 1.0 ! =1 Applies constraint to T (if ==1.0)
!                  W_SUM_ONE1 = 0.0 ! =1 Applies constraint to T (if ==1.0)
                  W_SUM_ONE2 = 0.0 ! =1 Applies constraint to Told (if ==1.0)
                  ! the original working code used W_SUM_ONE1 = 1, W_SUM_ONE2 = 1
                  CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) * MEAN_PORE_CV( CV_NODI ) *( &
                       (1.0-W_SUM_ONE1) *  T( CV_NODI_IPHA ) / DT &
                       -(1.0-W_SUM_ONE2) *  TOLD( CV_NODI_IPHA ) / DT  &
                       +TOLD( CV_NODI_IPHA ) * (DEN( CV_NODI_IPHA )-DENOLD( CV_NODI_IPHA )) / ( DT * DEN( CV_NODI_IPHA ) ) &
                       -DERIV( CV_NODI_IPHA ) * CV_P( CV_NODI ) * T( CV_NODI_IPHA ) / ( DT * DEN( CV_NODI_IPHA ) ) &
                       )
                  IF(IPHASE==1) THEN ! Add constraint to force sum of volume fracts to be unity...
                     ! W_SUM_ONE==1 applies the constraint
                     ! W_SUM_ONE==0 does NOT apply the constraint
                     CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) * &
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

      if( .false. .and. getct) then
         print *,'******'
         do CV_NODI=1,cv_nonods
            print *,'cv_nodi,ct:',cv_nodi,(ct(count),count=findct(cv_nodi),findct(cv_nodi+1)-1)
         end do
!         print *,'ct(1:ncolct);',ct(1:ncolct)
!          stop 1761
!         ewrite(3,*)'ct(1:ncolct);',ct(1:ncolct)
!         ewrite(3,*)'ct(1+ncolct:2*ncolct);',ct(1+ncolct:2*ncolct)
      end if

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
!      DEALLOCATE( PERM_ELE )

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
      DEALLOCATE( SCVFEWEIGH )

      DEALLOCATE( SUFEN )
      DEALLOCATE( SUFENSLX )
      DEALLOCATE( SUFENSLY )
      DEALLOCATE( SUFENLX )
      DEALLOCATE( SUFENLY )
      DEALLOCATE( SUFENLZ )

      DEALLOCATE( SRA )

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

    SUBROUTINE CV_ASSEMB_ADV_DIF( state, &
         LIMTOLD,LIMT2OLD,LIMDOLD,LIMDTOLD,LIMDTT2OLD,NDOTQOLD,&
         CV_RHS, &
         NCOLACV, CSR_ACV, dense_acv, FINACV, COLACV, MIDACV, &
         SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM,&
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
         SUF_T_BC_ROB1, SUF_T_BC_ROB2, &
         WIC_T_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, CV_P, &
         SOURCT, ABSORBT, VOLFRA_PORE, &
         NDIM, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         T_FEMT, DEN_FEMT, &
         IGOT_T2, T2, T2OLD,SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
         THETA_GDIFF, &
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         FINDCMC, COLCMC, NCOLCMC, MASS_MN_PRES, THERMAL, &
         MASS_ELE_TRANSP, &
         option_path_spatial_discretisation, &
         THETA_FLUX, ONE_M_THETA_FLUX,&
         StorageIndexes)

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
      !     CSR_ACV   - Matrix for assembling the advection terms (empty on input)
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
      !     CSR_ACV   - Matrix updated to include the advection terms
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
      type( state_type ), dimension( : ), intent( inout ) :: state
      real, DIMENSION( :, : ), intent( in ) :: LIMTOLD,LIMDOLD,LIMT2OLD,LIMDTOLD,LIMDTT2OLD, ndotqold
      INTEGER, intent( in ) :: NCOLACV, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_DISOPT, CV_DG_VEL_INT_OPT, NDIM, &
           NCOLM, XU_NLOC, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NCOLCMC
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(: ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T2_BC
      REAL, DIMENSION( : ), intent( inout ) :: CV_RHS
      !REAL, DIMENSION( NCOLACV ), intent( inout ) :: ACV
      real, dimension(:), intent(inout) :: CSR_ACV
      real, dimension(:,:,:), intent(inout) :: dense_acv
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES

      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ), intent( in ) :: U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( : ), intent( in ) :: T, TOLD, DEN, DENOLD
      REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( : ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION( :, : ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
      REAL, intent( in ) :: DT, CV_THETA, SECOND_THETA, CV_BETA
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION(: ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: CV_P
      REAL, DIMENSION( : ), intent( in ) :: SOURCT
      REAL, DIMENSION( :, :, : ), intent( in ) :: ABSORBT
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      LOGICAL, intent( in ) ::  GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION(: ), intent( inout ) :: T_FEMT, DEN_FEMT
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
      REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE_TRANSP
      character( len = * ), intent( in ), optional :: option_path_spatial_discretisation
      integer, dimension(:), intent(in) :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      !character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation
      integer, dimension(:), intent(inout) :: StorageIndexes

      ! Local variables
      INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
           WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
           WIC_U_BC_DIRICHLET = 1
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE
      INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, &
           CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
           COLGPTS, CV_SLOC2LOC, U_SLOC2LOC, &
           TMAX_NOD, TMIN_NOD, DENMAX_NOD, DENMIN_NOD,&
           T2MAX_NOD, T2MIN_NOD
      INTEGER, DIMENSION( : , : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, &
           FACE_ELE, CV_NEILOC
      REAL, DIMENSION( : ), allocatable :: CVWEIGHT, CVWEIGHT_SHORT, SCVFEWEIGH, SBCVFEWEIGH, &
           TMAX, TMIN, DENMAX, DENMIN, &
           TMAX_2ND_MC, TMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, &
           CVNORMX, &
           CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
           FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
           SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
           SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
           T2MAX, T2MIN,&
           T2MAX_2ND_MC, T2MIN_2ND_MC,&
           UP_WIND_NOD
      REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
           UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ, &
           DUMMY_ZERO_NDIM_NDIM
      REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE,  &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
      REAL, pointer,  DIMENSION( : , :, : ) :: INV_JAC
      real, pointer, dimension(:) :: SCVDETWEI, SCVRA
      real, pointer, dimension(:,:,:) :: SCVFENX_ALL
      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE,   &
           NCOLGPTS, &
           CV_SILOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, rhs_nodj_ipha,rhs_nodi_ipha,&
           CV_NODI, CV_NODI_IPHA, U_NODK, TIMOPT, &
           JCOUNT_IPHA, IMID_IPHA, &
           NFACE, X_NODI,  &
           CV_INOD, MAT_NODI, FACE_ITS, NFACE_ITS, &
           CVNOD, XNOD, NSMALL_COLM, NOD
      !        ===>  REALS  <===
      REAL :: NDOTQ,  &
           INCOME, HDC, FVT, FVT2,&
           FVD, LIMT, LIMT2,&
           LIMD, FTHETA, VTHETA, &
           LIMDT, LIMDTT2, &
           FEMDGI, FEMTGI,FEMT2GI, &
           TMID, TOLDMID, &
           DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, BCZERO, ROBIN1, ROBIN2, &
           RSUM, &
           FTHETA_T2, ONE_M_FTHETA_T2OLD, THERM_FTHETA, &
           W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW

      REAL, PARAMETER :: W_SUM_ONE = 1.
      real, pointer :: VOLUME
      integer :: U_NODK_IPHA, IANISOLIM
      logical :: Have_Temperature_Fields, Have_VolumeFraction_Fields, Have_Components_Fields
      logical :: overlapping
      ! Functions...
      !REAL :: R2NORM, FACE_THETA
      !        ===>  LOGICALS  <===
      character( len = option_path_len ) :: overlapping_path
      LOGICAL :: GETMAT, LIMIT_USE_2ND, &
           D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
           NORMALISE, GET_GTHETA, QUAD_OVER_WHOLE_ELE,is_overlapping

      character( len = option_path_len ) :: option_path, option_path2, path_temp, path_volf, &
           path_comp, path_spatial_discretisation


      real, dimension(:), allocatable :: TUPWIND_MAT,  DENUPWIND_MAT, &
              T2UPWIND_MAT
      INTEGER :: IDUM(1)
      REAL :: RDUM(1)
     !Reals to store the irresidual water and irreducible oil values, used in GET_INT_T_DEN
    real ::  s_gc, s_or

    logical got_theta_flux
    integer :: global_face

    got_theta_flux=present(theta_flux)
    global_face=0

      overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) overlapping = .true.

      IDUM = 0
      RDUM = 0.

      ewrite(3,*) 'In CV_ASSEMB_ADV_DIF'

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

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.


      GOT_DIFFUS = ( R2NORM( TDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0 )

      ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS:', &
           CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS

      ndotq = 0. 

      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

      ! Allocate memory for the control volume surface shape functions, etc.

      ALLOCATE( CVNORMX( SCVNGI ))
      ALLOCATE( CVNORMY( SCVNGI ))
      ALLOCATE( CVNORMZ( SCVNGI ))
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
      ALLOCATE( SCVFEWEIGH( SCVNGI ))

      ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

      ALLOCATE( SRA( SCVNGI ))

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
      ! The porocity mapped to the CV nodes
      ALLOCATE( SUM_CV( CV_NONODS ))
      ALLOCATE( UP_WIND_NOD( CV_NONODS * NPHASE ))


      UP_WIND_NOD = 0.0

      ALLOCATE( ONE_PORE( TOTELE ))
      IF( have_option( '/porous_media/actual_velocity' ) ) THEN
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
      ALLOCATE(      T2MIN( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(      T2MAX( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(    DENMIN( CV_NONODS * NPHASE) )
      ALLOCATE(    DENMAX( CV_NONODS * NPHASE) )

      ALLOCATE(      TMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(      TMAX_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(      T2MIN_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(      T2MAX_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(    DENMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE(    DENMAX_2ND_MC( CV_NONODS * NPHASE) )

      ALLOCATE(      TMIN_NOD( CV_NONODS * NPHASE) )
      ALLOCATE(      TMAX_NOD( CV_NONODS * NPHASE) )
      ALLOCATE(      T2MIN_NOD( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(      T2MAX_NOD( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE(    DENMIN_NOD( CV_NONODS * NPHASE) )
      ALLOCATE(    DENMAX_NOD( CV_NONODS * NPHASE) )

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

      CALL PROJ_CV_TO_FEM_4( state, &
           FEMT, FEMTOLD, FEMDEN, FEMDENOLD, T, TOLD, DEN, DENOLD, &
           IGOT_T2, T2,T2OLD, FEMT2,FEMT2OLD, &
           XC_CV, YC_CV, ZC_CV, MASS_CV, MASS_ELE,  &
           NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI_SHORT, CV_NLOC, CVN_SHORT, CVWEIGHT_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
           0, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC )

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
      CALL SURRO_CV_MINMAX_1time( TMAX, TMIN, DENMAX, DENMIN, &
           T2MAX, T2MIN, &
           TMAX_2ND_MC, TMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, &
           LIMIT_USE_2ND, &
           T, T2, DEN, IGOT_T2, NPHASE, CV_NONODS, size(small_colm), SMALL_FINDRM, SMALL_COLM, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC, SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
           WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
           WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, &
           T2MIN_NOD, T2MAX_NOD, &
           DENMIN_NOD, DENMAX_NOD )

! **********ANISOTROPIC LIMITING...*******************
      IANISOLIM=0
     ! print *,'CV_DISOPT=',CV_DISOPT
     ! stop 72
      IF(CV_DISOPT.GE.5) IANISOLIM=1
! limiting not ready yet for P2 DG:
!      IF(
!      IF( (NDIM==2).AND.(CV_NLOC==6) ) IANISOLIM=0
!      IF( (NDIM==3).AND.(CV_NLOC==10) ) IANISOLIM=0
!IANISOLIM=0
      IF (IANISOLIM==0) THEN
         ALLOCATE(TUPWIND_MAT(1), DENUPWIND_MAT(1))
         ALLOCATE(T2UPWIND_MAT(1))
         NSMALL_COLM=1
      ELSE ! IF (IANISOLIM==1) THEN
         ! Reduce matrix size...
         NSMALL_COLM=size(small_colm)

!        !We need to allocate CSR_ADV
!        allocate(CSR_ACV(NPHASE*NSMALL_COLM))

         !ALLOCATE(SMALL_FINDRM(CV_NONODS+1),SMALL_COLM(NSMALL_COLM),SMALL_CENTRM(CV_NONODS))
         ALLOCATE(TUPWIND_MAT(NSMALL_COLM*NPHASE), &
              DENUPWIND_MAT(NSMALL_COLM*NPHASE))
         ALLOCATE(T2UPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2))

         CALL CALC_ANISOTROP_LIM_1time(&
              ! Caculate the upwind values stored in matrix form...
              T,DEN,T2, &
              FEMT,FEMDEN,FEMT2,(CV_NONODS.NE.X_NONODS), &
              TUPWIND_MAT, DENUPWIND_MAT, &
              T2UPWIND_MAT, &
              ! Store the upwind element for interpolation and its weights for
              ! faster results...
              IGOT_T2,NPHASE,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
              SMALL_FINDRM,SMALL_CENTRM,SMALL_COLM,NSMALL_COLM, &
              X_NDGLN,X_NONODS,NDIM, &
              X,Y,Z, XC_CV, YC_CV, ZC_CV )

      END IF
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
              SBCVFEN, SBCVFENSLX, SBCVFENSLY&
              ,  state,"ADVDIF", StorageIndexes(17) )
      ENDIF

      !     =============== DEFINE THETA FOR TIME-STEPPING ===================
      ! Define the type of time integration:
      ! Timopt is 0 if CV_DISOPT is even (theta specified);
      ! Timopt is 1 if CV_DISOPT is odd (non-linear theta scheme)
      TIMOPT = MOD( CV_DISOPT, 2 )

      VTHETA = 1.0
      FTHETA = CV_THETA

      CV_RHS = 0.0
      IF( GETMAT ) THEN
         CSR_ACV = 0.0
      ENDIF


      GET_GTHETA=.FALSE.
      IF( got_theta_flux ) THEN
         IF( GET_THETA_FLUX ) THEN
!            THETA_FLUX = 0.0
!            ONE_M_THETA_FLUX = 0.0
            IF( IGOT_T2 == 1 ) THEN
               GET_GTHETA = .TRUE.
               THETA_GDIFF = 0.0
            END IF
         ENDIF
      ENDIF

     !Get the irreducible water and residual oil before entering the loop
     !this value are used in GET_INT_T_DEN to keep the phases between realistic values
     !Default value has to be the same as in subroutine get_corey_options in multi_eos.F90
      !Default value of  S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
          s_gc, default=0.1)
      !Default value of    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
           s_or, default=0.3)


      ! Now we begin the loop over elements to assemble the advection terms
      ! into the matrix (ACV) and the RHS
      Loop_Elements: DO ELE = 1, TOTELE


         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR_INVJAC( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
              CV_NLOC, SCVNGI, &
              SCVFEN, SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH, SCVDETWEI, SCVRA, VOLUME, D1, D3, DCYL, &
              SCVFENX_ALL, &
              NDIM, INV_JAC, state, "advI", StorageIndexes(27))


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
                       X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                       CV_SNLOC, CVFEM_ON_FACE(:,GI), X_SHARE, X_NONODS, ELE, ELE2,  &
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
                        CALL CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                             FACE_ELE, NFACE, CVFEM_ON_FACE( :, GI ), &
                             CV_NONODS, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
                             CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )
                        !EWRITE(3,*)'*****AFTER CALC_SELE SELE,CV_SILOC,CV_SNLOC:',SELE,CV_SILOC,CV_SNLOC
                     END IF
                     INTEGRAT_AT_GI=.NOT.((ELE==ELE2).AND.(SELE==0))
                  END IF

               END IF Conditional_CheckingNeighbourhood

               ! avoid indegrating across the middle of a CV on the boundaries of elements
               Conditional_integration: IF(INTEGRAT_AT_GI) THEN

                  global_face=global_face+1

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

                  JCOUNT_IPHA = -99999
                  DO COUNT = SMALL_FINDRM( CV_NODI ), SMALL_FINDRM( CV_NODI + 1 ) - 1
                     IF( SMALL_COLM( COUNT ) == CV_NODJ ) THEN 
                        JCOUNT_IPHA = COUNT
                        EXIT
                     END IF!An exit may improve the performance!!!
                  END DO

                  Loop_IPHASE: DO IPHASE = 1, NPHASE

                     CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                     CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS
                     rhs_NODI_IPHA = IPHASE + (cv_nodi-1 ) * nphase
                     rhs_NODJ_IPHA = IPHASE + (cv_nodj-1 ) * nphase

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
                            CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2, FEMT2, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT2, FEMDGI, FEMT2GI,  UP_WIND_NOD, &
                                T2MIN, T2MAX,  T2MIN_NOD, T2MAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC,  T2MAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2UPWIND_MAT )
                        ELSE
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T, FEMT, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
                                TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TMIN_2ND_MC, TMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TUPWIND_MAT )
                        ENDIF

                        !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                        ! Calculate T and DEN on the CV face at quadrature point GI.
                        
                        CALL GET_INT_T_DEN( FVT, FVT2, FVD, LIMD, LIMT, LIMT2, &
                             LIMDT, LIMDTT2,&
                             FEMDGI, FEMTGI,FEMT2GI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, &
                             T, T2, DEN, FEMT, FEMT2, FEMDEN, &
                             TMIN, T2MIN, DENMIN, &
                             TMAX, T2MAX, DENMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TMIN_NOD, TMAX_NOD, &
                             DENMIN_NOD, DENMAX_NOD, &
                             T2MIN_NOD, T2MAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, T2MIN_2ND_MC, DENMIN_2ND_MC, &
                             TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT, &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )

                     END DO

                     ! Define face value of theta
                     IF(IGOT_T2==1) THEN
                        FTHETA=FACE_THETA( DT, CV_THETA, ( cv_disopt>=8 ),HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                             T(CV_NODJ_IPHA)*DEN(CV_NODJ_IPHA)*T2(CV_NODJ_IPHA), &
                             T(CV_NODI_IPHA)*DEN(CV_NODI_IPHA)*T2(CV_NODI_IPHA), &
                             NDOTQOLD(iphase,global_face), LIMDTT2OLD(iphase,global_face), DIFF_COEFOLD_DIVDX, &
                             TOLD(CV_NODJ_IPHA)*DENOLD(CV_NODJ_IPHA)*T2OLD(CV_NODJ_IPHA), &
                             TOLD(CV_NODI_IPHA)*DENOLD(CV_NODI_IPHA)*T2OLD(CV_NODI_IPHA) )
                     ELSE
                        FTHETA=FACE_THETA( DT, CV_THETA, ( cv_disopt>=8 ),HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                             T(CV_NODJ_IPHA)*DEN(CV_NODJ_IPHA), &
                             T(CV_NODI_IPHA)*DEN(CV_NODI_IPHA), &
                             NDOTQOLD(iphase,global_face), LIMDTT2OLD(iphase,global_face), DIFF_COEFOLD_DIVDX, &
                             TOLD(CV_NODJ_IPHA)*DENOLD(CV_NODJ_IPHA), &
                             TOLD(CV_NODI_IPHA)*DENOLD(CV_NODI_IPHA) )
                     ENDIF

                     FTHETA_T2         = FTHETA * LIMT2
                     ONE_M_FTHETA_T2OLD = (1.0-FTHETA) * LIMT2OLD(iphase,global_face)

                     IF( got_theta_flux ) THEN ! inner loop
                        IF ( GET_THETA_FLUX ) THEN
!                           THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)      =FTHETA*LIMT
!                           ONE_M_THETA_FLUX(ELE, CV_ILOC, GI, IPHASE)=(1.0-FTHETA)*LIMTOLD
                           THETA_FLUX(IPHASE,global_face)      =FTHETA*LIMDT/DEN(CV_NODI_IPHA)
                           ONE_M_THETA_FLUX(IPHASE,global_face)=(1.0-FTHETA)*LIMDTOLD(iphase,global_face)/DEN(CV_NODI_IPHA)
                        ENDIF
                        IF ( USE_THETA_FLUX ) THEN
                           FTHETA_T2         =THETA_FLUX(IPHASE,global_face)
                           ONE_M_FTHETA_T2OLD=ONE_M_THETA_FLUX(IPHASE,global_face)
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


      
      IF(GETMAT) THEN
         
         ! - Calculate the integration of the limited, high-order flux over a face
         ! Conservative discretisation. The matrix (PIVOT ON LOW ORDER SOLN)
         IF( ( CV_NODI_IPHA /= CV_NODJ_IPHA ) .AND. ( CV_NODJ /= 0 ) ) THEN
            CSR_ACV( IPHASE+(JCOUNT_IPHA-1)*NPHASE ) =  CSR_ACV( IPHASE+(JCOUNT_IPHA-1)*NPHASE ) &
                 + SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW * INCOME * LIMD  & ! advection
                 - FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX ! Diffusion contribution
            
            IF(GET_GTHETA) THEN
               THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                    + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX * T(CV_NODJ_IPHA) ! Diffusion contribution
            ENDIF
         ELSE IF(SELE/=0) THEN
            IF(WIC_T_BC(SELE+(IPHASE-1)*STOTEL) == WIC_T_BC_DIRICHLET) THEN
               CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA )  &
                    + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &
                    * SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC+(IPHASE-1)*STOTEL*CV_SNLOC)
               IF(GET_GTHETA) THEN
                  THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                       + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX &
                       * SUF_T_BC(CV_SILOC+(SELE-1)*CV_SNLOC+(IPHASE-1)*STOTEL*CV_SNLOC)
               ENDIF
            ENDIF
         ENDIF
                           
         IMID_IPHA =  IPHASE+(SMALL_CENTRM(CV_NODI)-1)*NPHASE
         
         CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA ) &
              +  SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW * ( 1. - INCOME ) * LIMD & ! advection
              +  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX  &  ! Diffusion contribution
              +  SCVDETWEI( GI ) * ROBIN1  ! Robin bc
         
         IF(GET_GTHETA) THEN
            THETA_GDIFF( CV_NODI_IPHA ) =  THETA_GDIFF( CV_NODI_IPHA ) &
                 -  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX * T( CV_NODI_IPHA ) & ! Diffusion contribution
                 -  SCVDETWEI( GI ) * ROBIN1 * T( CV_NODI_IPHA )  ! Robin bc
         ENDIF
         
         ! CV_BETA=0 for Non-conservative discretisation (CV_BETA=1 for conservative disc)
         CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA )  &
              - SECOND_THETA * FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW * LIMD
      ENDIF
      
      TMID   =     T( CV_NODI_IPHA )
      TOLDMID = TOLD( CV_NODI_IPHA )
      
      ! Make allowances for no matrix stencil operating from outside the boundary.
      BCZERO=1.0
      IF( (SELE /= 0) .AND. (INCOME > 0.5) ) BCZERO=0.0
      
      ! Put results into the RHS vector
      CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA )  &
           ! subtract 1st order adv. soln.
           + SECOND_THETA * FTHETA_T2 * NDOTQNEW * SCVDETWEI( GI ) * LIMD * FVT * BCZERO &
           -  SCVDETWEI( GI ) * ( FTHETA_T2 * NDOTQNEW * LIMDT &
           + ONE_M_FTHETA_T2OLD * NDOTQOLD(iphase,global_face) * LIMDTOLD(iphase,global_face) ) ! hi order adv
      
      ! Subtract out 1st order term non-conservative adv.
      CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA ) &
           - FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW * LIMD * TMID
      
      ! High-order non-conservative advection contribution
      CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA ) &
           + ( 1. - CV_BETA) * SCVDETWEI( GI ) &
           * ( FTHETA_T2 * NDOTQNEW * TMID * LIMD  &
           + ONE_M_FTHETA_T2OLD * NDOTQOLD (iphase,global_face) * LIMDOLD(iphase,global_face) * TOLDMID )

                        ! Diffusion contribution
      CV_RHS( rhs_NODI_IPHA ) =  CV_RHS( rhs_NODI_IPHA ) &
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
            CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                 - CV_P( CV_NODI ) * SCVDETWEI( GI ) * ( &
                 THERM_FTHETA * NDOTQNEW * LIMT2 &
                 + ( 1. - THERM_FTHETA ) * NDOTQOLD(iphase,global_face) * LIMT2OLD(iphase,global_face) )
         else
            CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                 - CV_P( CV_NODI ) * SCVDETWEI( GI ) * ( &
                 THERM_FTHETA * NDOTQNEW &
                 + ( 1. - THERM_FTHETA ) * NDOTQOLD(iphase,global_face) )
         end if

      END IF


                  END DO Loop_IPHASE

               END IF Conditional_integration

            END DO Loop_GCOUNT

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
      
      Loop_CVNODI2: DO CV_NODI = 1, CV_NONODS ! Put onto the diagonal of the matrix
         
         DO COUNT = SMALL_FINDRM( CV_NODI ), SMALL_FINDRM( CV_NODI + 1 ) - 1, 1
            IF( SMALL_COLM( COUNT ) == CV_NODI )  then
               JCOUNT_IPHA = COUNT !An exit may improve the performance!!!
               exit
            end IF
         END DO


         Loop_IPHASE2: DO IPHASE = 1, NPHASE
            CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
            rhs_NODI_IPHA = iphase + ( cv_nodi -1 ) * nphase
            ! For the gravity term
            !SOURCT2( CV_NODI_IPHA ) = SOURCT( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA )
            !SOURCT2( CV_NODI_IPHA ) = SOURCT2( CV_NODI_IPHA ) * DEN( CV_NODI_IPHA )
            !               IMID_IPHA = MIDACV( CV_NODI_IPHA )
            IMID_IPHA = IPHASE+ (SMALL_CENTRM(CV_NODI)-1)*NPHASE
            

            IF(IGOT_T2==1) THEN
               CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                    + MASS_CV(CV_NODI) * SOURCT(CV_NODI_IPHA) !&
               !- CV_BETA * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) & ! conservative time term.
               !* TOLD(CV_NODI_IPHA)  &
               !* (DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA) - DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA)) / DT


               CSR_ACV( IMID_IPHA ) = CSR_ACV( IMID_IPHA ) &
                    !+ (CV_BETA * DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA) &
                    + (CV_BETA * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA) &
                    + (1.-CV_BETA) * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA))  &
                    * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) / DT


               CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                    + (CV_BETA * DENOLD(CV_NODI_IPHA) * T2OLD(CV_NODI_IPHA) &
                    + (1.-CV_BETA) * DEN(CV_NODI_IPHA) * T2(CV_NODI_IPHA))  &
                    * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) * TOLD(CV_NODI_IPHA) / DT
            ELSE

               CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                    + MASS_CV(CV_NODI) * SOURCT(CV_NODI_IPHA) !&
               !- CV_BETA * MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) & ! conservative time term.
               !* TOLD(CV_NODI_IPHA) &
               !* (DEN(CV_NODI_IPHA) - DENOLD(CV_NODI_IPHA)) / DT

               CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA ) &
                    !+ (CV_BETA * DENOLD(CV_NODI_IPHA) + (1.-CV_BETA) * DEN(CV_NODI_IPHA))  &
                    + (CV_BETA * DEN(CV_NODI_IPHA) &
                    + (1.-CV_BETA) * DEN(CV_NODI_IPHA))  &
                    * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) / DT

               CV_RHS( rhs_NODI_IPHA ) = CV_RHS( rhs_NODI_IPHA ) &
                    + (CV_BETA * DENOLD(CV_NODI_IPHA) &
                    + (1.-CV_BETA) * DEN(CV_NODI_IPHA)) &
                    * MEAN_PORE_CV(CV_NODI) * MASS_CV(CV_NODI) * TOLD(CV_NODI_IPHA) / DT
            ENDIF


            Conditional_GETMAT2: IF( GETMAT ) THEN
                  
               DO JPHASE = 1, NPHASE
                  !CV_NODI_JPHA = CV_NODI + ( JPHASE - 1 ) * CV_NONODS
                  !here we should put the full_ACV
                  dense_acv(JPHASE,IPHASE,CV_NODI)  = dense_acv(JPHASE,IPHASE,CV_NODI)  &
                       +  MASS_CV( CV_NODI ) * ABSORBT( CV_NODI, IPHASE, JPHASE )
               END DO

            end IF Conditional_GETMAT2
            
         END DO Loop_IPHASE2
         
      END DO Loop_CVNODI2

      ! for the output
      T_FEMT = FEMT
      DEN_FEMT = FEMDEN

      ewrite(3,*) '----------sub cv_assemb--------'
      ewrite(3,*)'IN cv_adv_dif a CV representation t:'
      CALL PRINT_CV_DIST(CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE, &
           T, X_NDGLN, CV_NDGLN, X)
      ewrite(3,*) 'just print out - in cv_assemb_adv_adif'

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
      DEALLOCATE( SCVFEWEIGH )

      DEALLOCATE( SUFEN )
      DEALLOCATE( SUFENSLX )
      DEALLOCATE( SUFENSLY )
      DEALLOCATE( SUFENLX )
      DEALLOCATE( SUFENLY )
      DEALLOCATE( SUFENLZ )

      DEALLOCATE( SRA )

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
      DEALLOCATE( DENMIN )
      DEALLOCATE( DENMAX )

      DEALLOCATE( TMIN_2ND_MC )
      DEALLOCATE( TMAX_2ND_MC )
      DEALLOCATE( DENMIN_2ND_MC )
      DEALLOCATE( DENMAX_2ND_MC )

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

    END SUBROUTINE CV_ASSEMB_ADV_DIF

    SUBROUTINE CV_ASSEMB_CT( state, &
         SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM,&
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
         SUF_T_BC_ROB1, SUF_T_BC_ROB2, &
         WIC_T_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, CV_P, &
         SOURCT, ABSORBT, VOLFRA_PORE, &
         NDIM, &
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
         option_path_spatial_discretisation,&
         StorageIndexes )

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
      !     CSR_ACV   - Matrix for assembling the advection terms (empty on input)
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
      !     CSR_ACV   - Matrix updated to include the advection terms
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
      type( state_type ), dimension( : ), intent( inout ) :: state
      INTEGER, intent( in ) :: NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_DISOPT, CV_DG_VEL_INT_OPT, NDIM, &
           NCOLM, XU_NLOC, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NCOLCMC
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(: ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T2_BC
      REAL, DIMENSION( : ), intent( inout ) :: CT
      ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
      REAL, DIMENSION( : ), intent( inout ) :: DIAG_SCALE_PRES
      REAL, DIMENSION( :  ), intent( inout ) :: CT_RHS
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES

      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ), intent( in ) :: U, V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
      REAL, DIMENSION( : ), intent( in ) :: T, TOLD, DEN, DENOLD
      REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( : ), intent( inout ) :: THETA_GDIFF
      REAL, DIMENSION( :, :), intent( inout ) :: THETA_FLUX, ONE_M_THETA_FLUX
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
      REAL, intent( in ) :: DT, CV_THETA, SECOND_THETA, CV_BETA
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION(: ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: CV_P
      REAL, DIMENSION( : ), intent( in ) :: SOURCT
      REAL, DIMENSION( :, :, : ), intent( in ) :: ABSORBT
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      LOGICAL, intent( in ) :: GET_THETA_FLUX, USE_THETA_FLUX, THERMAL
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION(: ), intent( inout ) :: T_FEMT, DEN_FEMT
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
      REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE_TRANSP
      character( len = * ), intent( in ), optional :: option_path_spatial_discretisation
      integer, dimension(:), intent(in) :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      !character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation
      integer, dimension(:), intent(inout) :: StorageIndexes

      ! Local variables
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
           CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
           FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
           SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
           SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
           T2MAX, T2MIN, T2OLDMAX, &
           T2OLDMIN, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, &
           T2OLDMIN_2ND_MC, &
           UP_WIND_NOD, DU, DV, DW, PERM_ELE
      REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
           UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ, &
           DUMMY_ZERO_NDIM_NDIM
      REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE,  &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
      REAL, pointer, DIMENSION( : , :, : ) :: INV_JAC
      real, pointer, dimension(:) :: SCVDETWEI, SCVRA
      real, pointer, dimension(:,:,:) :: SCVFENX_ALL
      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE,   &
           NCOLGPTS, &
           CV_SILOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, rhs_nodj_ipha,rhs_nodi_ipha,&
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
           W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW

      REAL, PARAMETER :: W_SUM_ONE = 1.
      real, pointer :: VOLUME
      integer :: cv_inod_ipha, IGETCT, U_NODK_IPHA, IANISOLIM
      logical :: Have_Temperature_Fields, Have_VolumeFraction_Fields, Have_Components_Fields
      logical :: overlapping
      ! Functions...
      !REAL :: R2NORM, FACE_THETA
      !        ===>  LOGICALS  <===
      character( len = option_path_len ) :: overlapping_path
      LOGICAL :: GETMAT, LIMIT_USE_2ND, &
           D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
           NORMALISE, SUM2ONE, GET_GTHETA, QUAD_OVER_WHOLE_ELE,is_overlapping

      character( len = option_path_len ) :: option_path, option_path2, path_temp, path_volf, &
           path_comp, path_spatial_discretisation


      real, dimension(:), allocatable :: TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, &
              DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT
      INTEGER :: IDUM(1)
      REAL :: RDUM(1),n1,n2,n3
      type( scalar_field ), pointer :: perm
     !Reals to store the irresidual water and irreducible oil values, used in GET_INT_T_DEN
    real ::  s_gc, s_or

    integer :: global_face
    global_face=0

      overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) overlapping = .true.

      IDUM = 0
      RDUM = 0.


      ewrite(3,*) 'In CV_ASSEMB_CT'

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

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.


      GOT_DIFFUS = ( R2NORM( TDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0 )

      ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS:', &
           CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA, GOT_DIFFUS

      ndotq = 0. ; ndotqold = 0.

      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

      ! Allocate memory for the control volume surface shape functions, etc.
      ALLOCATE( JCOUNT_KLOC( U_NLOC )) ; jcount_kloc = 0
      ALLOCATE( JCOUNT_KLOC2( U_NLOC )) ; jcount_kloc2 = 0

      ALLOCATE( CVNORMX( SCVNGI ))
      ALLOCATE( CVNORMY( SCVNGI ))
      ALLOCATE( CVNORMZ( SCVNGI ))
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
      ALLOCATE( SCVFEWEIGH( SCVNGI ))

      ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

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

      UP_WIND_NOD = 0.0

      ALLOCATE( ONE_PORE( TOTELE ))
      IF( have_option( '/porous_media/actual_velocity' ) ) THEN
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
      IGETCT=1

      CALL PROJ_CV_TO_FEM_4( state, &
           FEMT, FEMTOLD, FEMDEN, FEMDENOLD, T, TOLD, DEN, DENOLD, &
           IGOT_T2, T2,T2OLD, FEMT2,FEMT2OLD, &
           XC_CV, YC_CV, ZC_CV, MASS_CV, MASS_ELE,  &
           NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI_SHORT, CV_NLOC, CVN_SHORT, CVWEIGHT_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
           IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC )

      MASS_ELE_TRANSP = MASS_ELE

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

      ! For each node, find the largest and smallest value of T and
      ! DENSITY for both the current and previous timestep, out of
      ! the node value and all its surrounding nodes including Dirichlet b.c's.
      CALL SURRO_CV_MINMAX( TMAX, TMIN, TOLDMAX, TOLDMIN, DENMAX, DENMIN, DENOLDMAX, DENOLDMIN, &
           T2MAX, T2MIN, T2OLDMAX, T2OLDMIN, &
           TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, &
           DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, T2OLDMIN_2ND_MC, &
           LIMIT_USE_2ND, &
           T, TOLD, T2, T2OLD, DEN, DENOLD, IGOT_T2, NPHASE, CV_NONODS, size(small_colm), SMALL_FINDRM, SMALL_COLM, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC, SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
           WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
           WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
           T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, T2OLDMAX_NOD, &
           DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD )

! **********ANISOTROPIC LIMITING...*******************
      IANISOLIM=0
     ! print *,'CV_DISOPT=',CV_DISOPT
     ! stop 72
      IF(CV_DISOPT.GE.5) IANISOLIM=1
! limiting not ready yet for P2 DG:
!      IF(
!      IF( (NDIM==2).AND.(CV_NLOC==6) ) IANISOLIM=0
!      IF( (NDIM==3).AND.(CV_NLOC==10) ) IANISOLIM=0
!IANISOLIM=0
      IF (IANISOLIM==0) THEN
         ALLOCATE(TUPWIND_MAT(1), TOLDUPWIND_MAT(1), DENUPWIND_MAT(1), DENOLDUPWIND_MAT(1))
         ALLOCATE(T2UPWIND_MAT(1), T2OLDUPWIND_MAT(1))
         NSMALL_COLM=1
      ELSE ! IF (IANISOLIM==1) THEN
         ! Reduce matrix size...
         NSMALL_COLM=size(small_colm)

!        !We need to allocate CSR_ADV
!        allocate(CSR_ACV(NPHASE*NSMALL_COLM))

         !ALLOCATE(SMALL_FINDRM(CV_NONODS+1),SMALL_COLM(NSMALL_COLM),SMALL_CENTRM(CV_NONODS))
         ALLOCATE(TUPWIND_MAT(NSMALL_COLM*NPHASE), TOLDUPWIND_MAT(NSMALL_COLM*NPHASE), &
              DENUPWIND_MAT(NSMALL_COLM*NPHASE), DENOLDUPWIND_MAT(NSMALL_COLM*NPHASE))
         ALLOCATE(T2UPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2), T2OLDUPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2))

         CALL CALC_ANISOTROP_LIM(&
              ! Caculate the upwind values stored in matrix form...
              T,TOLD,DEN,DENOLD,T2,T2OLD, &
              FEMT,FEMTOLD,FEMDEN,FEMDENOLD,FEMT2,FEMT2OLD, (CV_NONODS.NE.X_NONODS), &
              TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, &
              T2UPWIND_MAT, T2OLDUPWIND_MAT, &
              ! Store the upwind element for interpolation and its weights for
              ! faster results...
              IGOT_T2,NPHASE,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
              SMALL_FINDRM,SMALL_CENTRM,SMALL_COLM,NSMALL_COLM, &
              X_NDGLN,X_NONODS,NDIM, &
              X,Y,Z, XC_CV, YC_CV, ZC_CV)

         if ( .false. ) then
            ewrite(3,*) 'TUPWIND_MAT', TUPWIND_MAT(1:nsmall_colm)
            ewrite(3,*) 'TOLDUPWIND_MAT', TOLDUPWIND_MAT(1:nsmall_colm)
            if (igot_t2==1)then
               ewrite(3,*) 'T2UPWIND_MAT', T2UPWIND_MAT(1:nsmall_colm)
               ewrite(3,*) 'T2OLDUPWIND_MAT', T2OLDUPWIND_MAT(1:nsmall_colm)
            end if
         end if

      END IF
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
              SBCVFEN, SBCVFENSLX, SBCVFENSLY&
              ,  state, "CT", StorageIndexes(18))
      ENDIF

      !     =============== DEFINE THETA FOR TIME-STEPPING ===================
      ! Define the type of time integration:
      ! Timopt is 0 if CV_DISOPT is even (theta specified);
      ! Timopt is 1 if CV_DISOPT is odd (non-linear theta scheme)
      TIMOPT = MOD( CV_DISOPT, 2 )

      VTHETA = 1.0
      FTHETA = CV_THETA

      
      CT_RHS = 0.0
      CT = 0.0

      GET_GTHETA=.FALSE.
      IF( IGOT_THETA_FLUX == 1 ) THEN
         IF( GET_THETA_FLUX ) THEN
            THETA_FLUX = 0.0
            ONE_M_THETA_FLUX = 0.0
            IF( IGOT_T2 == 1 ) THEN
               GET_GTHETA = .TRUE.
               THETA_GDIFF = 0.0
            END IF
         ENDIF
      end IF
     !Get the irreducible water and residual oil before entering the loop
     !this value are used in GET_INT_T_DEN to keep the phases between realistic values
     !Default value has to be the same as in subroutine get_corey_options in multi_eos.F90
      !Default value of  S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
          s_gc, default=0.1)
      !Default value of    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
           s_or, default=0.3)


      ! Now we begin the loop over elements to assemble the advection terms
      ! into the matrix (ACV) and the RHS
      Loop_Elements: DO ELE = 1, TOTELE


         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR_INVJAC( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
              CV_NLOC, SCVNGI, &
              SCVFEN, SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH, SCVDETWEI, SCVRA, VOLUME, D1, D3, DCYL, &
              SCVFENX_ALL,  &
              NDIM, INV_JAC, state, "ctI", StorageIndexes(28) )


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
                       X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                       CV_SNLOC, CVFEM_ON_FACE(:,GI), X_SHARE, X_NONODS, ELE, ELE2,  &
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
                        CALL CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                             FACE_ELE, NFACE, CVFEM_ON_FACE( :, GI ), &
                             CV_NONODS, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
                             CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )    
                        !EWRITE(3,*)'*****AFTER CALC_SELE SELE,CV_SILOC,CV_SNLOC:',SELE,CV_SILOC,CV_SNLOC
                     END IF
                     INTEGRAT_AT_GI=.NOT.((ELE==ELE2).AND.(SELE==0))
                  END IF

               END IF Conditional_CheckingNeighbourhood

               ! avoid indegrating across the middle of a CV on the boundaries of elements
               Conditional_integration: IF(INTEGRAT_AT_GI) THEN

                  ! if necessary determine the derivatives between elements ELE and ELE2

                  global_face=glObal_face+1

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

                  ! Obtain the CV discretised CT equations plus RHS
                     !ewrite(3,*)'==================================================================='
                     !ewrite(3,*)'CV_NODI, CV_NODJ, ELE, ELE2, GI:', CV_NODI, CV_NODJ, ELE, ELE2, GI
                     !ewrite(3,*)'findct:',FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1
                     !ewrite(3,*)'colct:',colct( FINDCT( CV_NODI ) : FINDCT( CV_NODI + 1 ) - 1 )
                     !ewrite(3,*)'SCVDETWEI:', SCVDETWEI(GI)

                     DO U_KLOC = 1, U_NLOC
                        U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
                        JCOUNT = 0
                        DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1, 1
                           IF(COLCT( COUNT ) == U_NODK) then
                              JCOUNT = COUNT
                              exit
                           end if
                        END DO
                        JCOUNT_KLOC( U_KLOC ) = JCOUNT
                        !ewrite(3,*)' u_nodk, jcount1:', u_nodk, jcount
                     END DO

                     IF( ( ELE2 /= 0 ) .AND. ( ELE2 /= ELE ) ) THEN
                        DO U_KLOC = 1, U_NLOC
                           U_NODK = U_NDGLN(( ELE2 - 1 ) * U_NLOC + U_KLOC )
                           JCOUNT = 0
                           DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1, 1
                              IF(COLCT( COUNT ) == U_NODK) then
                                 JCOUNT = COUNT
                                 exit
                              end if
                           END DO
                           JCOUNT_KLOC2( U_KLOC ) = JCOUNT
                           !ewrite(3,*)' u_nodk, jcount2:', u_nodk, jcount
                        END DO

                     ENDIF

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
                        rhs_NODI_IPHA = IPHASE + (cv_nodi-1 ) * nphase
                        rhs_NODJ_IPHA = IPHASE + (cv_nodj-1 ) * nphase

                        ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
                        IF(IGOT_T2==1) THEN
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW,  NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2OLD, FEMT2OLD, DENOLD, &
                                U, V, W, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT2OLD, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                                T2OLDMIN, T2OLDMAX, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2OLDMIN_2ND_MC, T2OLDMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2OLDUPWIND_MAT )
                            CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2, FEMT2, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1, LIMT2, FEMDGI, FEMT2GI,  UP_WIND_NOD, &
                                T2MIN, T2MAX,  T2MIN_NOD, T2MAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC, T2MAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2UPWIND_MAT )
                        ENDIF


                        !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                        ! Calculate T and DEN on the CV face at quadrature point GI.
                         CALL GET_INT_T_DEN( FVTOLD, FVT2OLD, FVDOLD, &
                             LIMDOLD, LIMTOLD, LIMT2OLD, LIMDTOLD, LIMDTT2OLD,&
                             FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOMEOLD, &
                             TOLD, T2OLD, DENOLD, FEMTOLD, FEMT2OLD, FEMDENOLD, &
                             TOLDMIN, T2OLDMIN, DENOLDMIN, &
                             TOLDMAX, T2OLDMAX, DENOLDMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TOLDMIN_NOD, TOLDMAX_NOD, &
                             DENOLDMIN_NOD, DENOLDMAX_NOD, &
                             T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2, &
                             TOLDMIN_2ND_MC, T2OLDMIN_2ND_MC, DENOLDMIN_2ND_MC, &
                             TOLDMAX_2ND_MC, T2OLDMAX_2ND_MC, DENOLDMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQOLD, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TOLDUPWIND_MAT, DENOLDUPWIND_MAT, T2OLDUPWIND_MAT , &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )
                        
                        CALL GET_INT_T_DEN( FVT, FVT2, FVD, LIMD, LIMT, LIMT2, &
                             LIMDT, LIMDTT2,&
                             FEMDGI, FEMTGI,FEMT2GI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, &
                             T, T2, DEN, FEMT, FEMT2, FEMDEN, &
                             TMIN, T2MIN, DENMIN, &
                             TMAX, T2MAX, DENMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TMIN_NOD, TMAX_NOD, &
                             DENMIN_NOD, DENMAX_NOD, &
                             T2MIN_NOD, T2MAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, T2MIN_2ND_MC, DENMIN_2ND_MC, &
                             TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:),CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT, &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )

                        SUM_LIMT    = SUM_LIMT    + LIMT
                        SUM_LIMTOLD = SUM_LIMTOLD + LIMTOLD

                     END DO Loop_IPHASE5
                  ENDIF Conditional_SUMLimiting
                  ! get the sum of limiting functions correct...************

                  DO COUNT = SMALL_FINDRM( CV_NODI ), SMALL_FINDRM( CV_NODI + 1 ) - 1
                     IF( SMALL_COLM( COUNT ) == CV_NODJ ) THEN 
                        JCOUNT_IPHA = COUNT
                        EXIT
                     END IF!An exit may improve the performance!!!
                  END DO

                  Loop_IPHASE: DO IPHASE = 1, NPHASE

                     CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                     CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS
                     rhs_NODI_IPHA = IPHASE + (cv_nodi-1 ) * nphase
                     rhs_NODJ_IPHA = IPHASE + (cv_nodj-1 ) * nphase

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
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW,  NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2OLD, FEMT2OLD, DENOLD, &
                                U, V, W, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT2OLD, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                                T2OLDMIN, T2OLDMAX, T2OLDMIN_NOD, T2OLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2OLDMIN_2ND_MC, T2OLDMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2OLDUPWIND_MAT )
                            CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2, FEMT2, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT2, FEMDGI, FEMT2GI,  UP_WIND_NOD, &
                                T2MIN, T2MAX,  T2MIN_NOD, T2MAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC, T2MAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2UPWIND_MAT )
                        ELSE
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                TOLD, FEMTOLD, DENOLD, &
                                U, V, W, NUOLD, NVOLD, NWOLD, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                1,  LIMTOLD, FEMDOLDGI, FEMT2OLDGI, UP_WIND_NOD, &
                                TOLDMIN, TOLDMAX, TOLDMIN_NOD, TOLDMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TOLDMIN_2ND_MC, TOLDMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TOLDUPWIND_MAT )
                           CALL GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T, FEMT, DEN, &
                                U, V, W, NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
                                TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TMIN_2ND_MC, TMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TUPWIND_MAT )
                        ENDIF

                        !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                        ! Calculate T and DEN on the CV face at quadrature point GI.
                        CALL GET_INT_T_DEN( FVTOLD, FVT2OLD, FVDOLD, &
                             LIMDOLD, LIMTOLD, LIMT2OLD, LIMDTOLD, LIMDTT2OLD,&
                             FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOMEOLD, &
                             TOLD, T2OLD, DENOLD, FEMTOLD, FEMT2OLD, FEMDENOLD, &
                             TOLDMIN, T2OLDMIN, DENOLDMIN, &
                             TOLDMAX, T2OLDMAX, DENOLDMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TOLDMIN_NOD, TOLDMAX_NOD, &
                             DENOLDMIN_NOD, DENOLDMAX_NOD, &
                             T2OLDMIN_NOD, T2OLDMAX_NOD, IGOT_T2, &
                             TOLDMIN_2ND_MC, T2OLDMIN_2ND_MC, DENOLDMIN_2ND_MC, &
                             TOLDMAX_2ND_MC, T2OLDMAX_2ND_MC, DENOLDMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQOLD, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TOLDUPWIND_MAT, DENOLDUPWIND_MAT, T2OLDUPWIND_MAT , &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )
                        
                        CALL GET_INT_T_DEN( FVT, FVT2, FVD, LIMD, LIMT, LIMT2, &
                             LIMDT, LIMDTT2,&
                             FEMDGI, FEMTGI,FEMT2GI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, &
                             T, T2, DEN, FEMT, FEMT2, FEMDEN, &
                             TMIN, T2MIN, DENMIN, &
                             TMAX, T2MAX, DENMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TMIN_NOD, TMAX_NOD, &
                             DENMIN_NOD, DENMAX_NOD, &
                             T2MIN_NOD, T2MAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, T2MIN_2ND_MC, DENMIN_2ND_MC, &
                             TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:), CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT, &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )

                     END DO

                     IF( SUM2ONE ) THEN
                        LIMT   =LIMT/SUM_LIMT
                        LIMTOLD=LIMTOLD/SUM_LIMTOLD
                        LIMDT=LIMT*LIMD
                        LIMDTOLD=LIMTOLD*LIMDOLD
                     ENDIF

                     ! Define face value of theta
                     IF(IGOT_T2==1) THEN
                        FTHETA=FACE_THETA( DT, CV_THETA, ( cv_disopt>=8 ),HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
                             T(CV_NODJ_IPHA)*DEN(CV_NODJ_IPHA)*T2(CV_NODJ_IPHA), &
                             T(CV_NODI_IPHA)*DEN(CV_NODI_IPHA)*T2(CV_NODI_IPHA), &
                             NDOTQOLD, LIMDTT2OLD, DIFF_COEFOLD_DIVDX, &
                             TOLD(CV_NODJ_IPHA)*DENOLD(CV_NODJ_IPHA)*T2OLD(CV_NODJ_IPHA), &
                             TOLD(CV_NODI_IPHA)*DENOLD(CV_NODI_IPHA)*T2OLD(CV_NODI_IPHA) )
                     ELSE
                        FTHETA=FACE_THETA( DT, CV_THETA, ( cv_disopt>=8 ),HDC, NDOTQ, LIMDTT2, DIFF_COEF_DIVDX, &
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
                           THETA_FLUX(IPHASE,global_face) =FTHETA*LIMDT/DEN(CV_NODI_IPHA)
                           ONE_M_THETA_FLUX(IPHASE,global_face)=(1.0-FTHETA)*LIMDTOLD/DEN(CV_NODI_IPHA)
                        ENDIF
                        IF ( USE_THETA_FLUX ) THEN
                           FTHETA_T2         =THETA_FLUX(IPHASE,global_face)
                           ONE_M_FTHETA_T2OLD=ONE_M_THETA_FLUX(IPHASE,global_face)
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

! ********FOR RELATIVE PERM B.C WHEN SPECIFYING PRESSURE*********
      IF( .false. ) THEN ! still bugs in the original better code but this is ok...
!      IF( is_overlapping ) THEN ! still bugs in the original better code but this is ok...

              !       RBC_OUT(1:NDIM)=1.0
              !       IF(SELE.NE.0) RBC_OUT(1:NDIM) &
              ! =SUF_SIG_DIAGTEN_BC(( SELE - 1 ) * CV_SNLOC + 1 +( IPHASE - 1 ) * STOTEL*CV_SNLOC,1:NDIM)

!                     IF((SELE.NE.0).AND.(NDOTQ.GE.0.0))  &
!                     IF((SELE.NE.0).AND.(NDOTQ.LE.0.0))  &
                     IF(SELE.NE.0)  &
                   SCVDETWEI( GI )=SCVDETWEI( GI ) &
               *SUF_SIG_DIAGTEN_BC(( SELE - 1 ) * CV_SNLOC + 1 +( IPHASE - 1 ) * STOTEL*CV_SNLOC,1)
      ENDIF
! ********FOR RELATIVE PERM B.C WHEN SPECIFYING PRESSURE*********

      !====================== CT AND RHS ASSEMBLY ===================
      
      CALL PUT_IN_CT_RHS(CT, CT_RHS, U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
           CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE, ELE, ELE2, SELE, &
           JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC, U_NDGLN, U, V, W, &
           SUFEN, SCVDETWEI, CVNORMX, CVNORMY, CVNORMZ, DEN, CV_NODI, CV_NODI_IPHA, &
           UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
           NDOTQNEW, NDOTQOLD, LIMDT, LIMDTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD )


        

   END DO Loop_IPHASE
   
               END IF Conditional_integration

            END DO Loop_GCOUNT

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

   
      DIAG_SCALE_PRES = 0.0 ! Obtain the CV discretised CT eqations plus RHS

      DO IPHASE = 1, NPHASE
         DO CV_NODI = 1, CV_NONODS
            CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
            rhs_NODI_IPHA = IPHASE + (cv_nodi-1 ) * nphase
            
            W_SUM_ONE1 = 1.0 ! =1 Applies constraint to T (if ==1.0)
            !                  W_SUM_ONE1 = 0.0 ! =1 Applies constraint to T (if ==1.0)
            W_SUM_ONE2 = 0.0 ! =1 Applies constraint to Told (if ==1.0)
                  ! the original working code used W_SUM_ONE1 = 1, W_SUM_ONE2 = 1
            CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) * MEAN_PORE_CV( CV_NODI ) *( &
                 (1.0-W_SUM_ONE1) *  T( CV_NODI_IPHA ) / DT &
                 -(1.0-W_SUM_ONE2) *  TOLD( CV_NODI_IPHA ) / DT  &
                 +TOLD( CV_NODI_IPHA ) * (DEN( CV_NODI_IPHA )-DENOLD( CV_NODI_IPHA )) / ( DT * DEN( CV_NODI_IPHA ) ) &
                 -DERIV( CV_NODI_IPHA ) * CV_P( CV_NODI ) * T( CV_NODI_IPHA ) / ( DT * DEN( CV_NODI_IPHA ) ) &
                 )
            IF(IPHASE==1) THEN ! Add constraint to force sum of volume fracts to be unity...
                     ! W_SUM_ONE==1 applies the constraint
               ! W_SUM_ONE==0 does NOT apply the constraint
               CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - MASS_CV( CV_NODI ) * &
                    (W_SUM_ONE1 - W_SUM_ONE2) * MEAN_PORE_CV( CV_NODI )  / DT
            ENDIF

        

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
 
      ! for the output
      T_FEMT = FEMT
      DEN_FEMT = FEMDEN

      ewrite(3,*) '----------sub cv_assemb--------'

      ewrite(3,*)'IN cv_adv_dif a CV representation t:'
      CALL PRINT_CV_DIST(CV_NONODS,X_NONODS,TOTELE,CV_NLOC,X_NLOC,NPHASE, &
           T, X_NDGLN, CV_NDGLN, X)
      ewrite(3,*) 'just print out - in cv_assemb_ct'

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
      DEALLOCATE( SCVFEWEIGH )

      DEALLOCATE( SUFEN )
      DEALLOCATE( SUFENSLX )
      DEALLOCATE( SUFENSLY )
      DEALLOCATE( SUFENLX )
      DEALLOCATE( SUFENLY )
      DEALLOCATE( SUFENLZ )

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

    END SUBROUTINE CV_ASSEMB_CT

    SUBROUTINE CV_GET_ALL_LIMITED_VALS( state, LIMT,LIMT2,LIMD,LIMDT,LIMDTT2,NDOTQOLD,&
         SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM,&
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         CV_ELE_TYPE,  &
         NPHASE,  &
         CV_NLOC, U_NLOC, X_NLOC, &
         CV_NDGLN, X_NDGLN, U_NDGLN, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
         X, Y, Z, U, V, W, &
         NU, NV, NW, &
         T, DEN, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
         CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, SECOND_THETA, CV_BETA, &
         SUF_T_BC, SUF_D_BC, SUF_U_BC, SUF_V_BC, SUF_W_BC, SUF_SIG_DIAGTEN_BC, &
         SUF_T_BC_ROB1, SUF_T_BC_ROB2, &
         WIC_T_BC, WIC_D_BC, WIC_U_BC, &
         DERIV, CV_P, &
         SOURCT, ABSORBT, VOLFRA_PORE, &
         NDIM, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         IGOT_T2, T2, SCVNGI_THETA,&
         SUF_T2_BC, SUF_T2_BC_ROB1, SUF_T2_BC_ROB2, WIC_T2_BC, IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         FINDCMC, COLCMC, NCOLCMC, &
         MASS_ELE_TRANSP, &
         option_path_spatial_discretisation, StorageIndexes)

      !  =====================================================================
      !     This subroutine gets the limited face values for a variable
      !
      !     This routine uses a Control Volume (CV) formulation to compute
      !     the advection terms. The general procedure is as follows:
      !
      !     This procedure is implemented by considering the CV to be made up
      !     of a number of sub-control-volumes, which represent the part of
      !     the control volume within a given element.  The assembly of terms
      !     considers each of these sub-CVs in turn, calculating (and limiting)
      !     the flux across sub-CV faces that are external to the CV...
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
      !     CSR_ACV   - Matrix for assembling the advection terms (empty on input)
      !     CV_RHS      - Right-hand side vector for advection-diffusion terms
      !     X,Y,Z    - Node co-ordinates
      !     NU       - Nodal velocity component
      !     T   - advected field values at nodes
      !     CV_DISOPT   - The discretisation/limiting option (see above)
      !     DT       - The time step
      !     ELE_TYP   - Integer flag definining element type
      !
      !     IMPORTANT OUTPUTS:
      !     -----------------
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
      type( state_type ), dimension( : ), intent( inout ) :: state
      real, dimension(:,:), intent(inout) :: LIMT,LIMT2,LIMD,LIMDT,LIMDTT2,NDOTQOLD
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_DISOPT, CV_DG_VEL_INT_OPT, NDIM, &
           NCOLM, XU_NLOC, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NCOLCMC, igot_t2
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(: ), intent( in ) :: U_SNDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC, WIC_D_BC, WIC_U_BC
      INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T2_BC
      ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC

      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      REAL, DIMENSION( : ), intent( in ) :: U, V, W
      REAL, DIMENSION( : ), intent( in ) :: NU, NV, NW
      REAL, DIMENSION( : ), intent( in ) :: T, DEN
      REAL, DIMENSION( : ), intent( in ) :: T2 
      REAL, intent( in ) :: DT, CV_THETA, SECOND_THETA, CV_BETA
      REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
      REAL, DIMENSION( :  ), intent( in ) :: SUF_T2_BC
      REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION(: ), intent( in ) :: SUF_T_BC_ROB1, SUF_T_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC_ROB1, SUF_T2_BC_ROB2
      REAL, DIMENSION( : ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: CV_P
      REAL, DIMENSION( : ), intent( in ) :: SOURCT
      REAL, DIMENSION( :, :, : ), intent( in ) :: ABSORBT
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
      REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE_TRANSP
      character( len = * ), intent( in ), optional :: option_path_spatial_discretisation
      integer, dimension(:), intent(in) :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      integer, dimension(:), intent(inout) :: StorageIndexes
      !character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation

      ! Local variables
      INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
           WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
           WIC_U_BC_DIRICHLET = 1
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE,LOG_ON_BOUND
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE
      INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, &
           CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
           COLGPTS, CV_SLOC2LOC, U_SLOC2LOC, &
           TMAX_NOD, TMIN_NOD, &
           DENMAX_NOD, DENMIN_NOD, &
           T2MAX_NOD, T2MIN_NOD
      INTEGER, DIMENSION( : , : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, &
           FACE_ELE, CV_NEILOC
      REAL, DIMENSION( : ), allocatable :: CVWEIGHT, CVWEIGHT_SHORT, SCVFEWEIGH, SBCVFEWEIGH, &
           TMAX, TMIN, TOLDMAX, &
           TOLDMIN, DENMAX, DENMIN, DENOLDMAX, DENOLDMIN, &
           TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, &
           TOLDMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC, &
           CVNORMX, &
           CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
           FEMT, FEMTOLD, FEMT2, FEMT2OLD, FEMDEN, FEMDENOLD, XC_CV, YC_CV, ZC_CV, &
           SRA, UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2,  &
           SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
           T2MAX, T2MIN, T2OLDMAX, &
           T2OLDMIN, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, &
           T2OLDMIN_2ND_MC, &
           UP_WIND_NOD, DU, DV, DW, PERM_ELE
      REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
           UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ, &
           DUMMY_ZERO_NDIM_NDIM
      REAL, DIMENSION( : , :, : ), allocatable :: DTX_ELE,DTY_ELE,DTZ_ELE,  &
           DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
      REAL, pointer, DIMENSION( : , :, : ) :: INV_JAC
      real, pointer, dimension(:) :: SCVRA, SCVDETWEI
      real, pointer, dimension(:,:,:) :: SCVFENX_ALL
      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE,   &
           NCOLGPTS, &
           CV_SILOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, rhs_nodj_ipha,rhs_nodi_ipha,&
           CV_NODI, CV_NODI_IPHA, CV_NODI_JPHA, U_NODK, TIMOPT, &
           JCOUNT_IPHA, IMID_IPHA, &
           NFACE, X_NODI,  &
           CV_INOD, MAT_NODI, FACE_ITS, NFACE_ITS, &
           CVNOD, XNOD, NSMALL_COLM, COUNT2, NOD
      !        ===>  REALS  <===
      REAL :: NDOTQ,& 
           INCOME, INCOMEOLD, HDC, FVT, FVTOLD, FVT2, FVT2OLD, &
           FVD, FVDOLD, FTHETA, VTHETA, &
           FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
           TMID, TOLDMID, &
           DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, BCZERO, ROBIN1, ROBIN2, &
           RSUM, &
           SUM_LIMT, SUM_LIMTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD, THERM_FTHETA, &
           W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW

      REAL, DIMENSION( NCOLCMC)  :: MASS_MN_PRES

      REAL, PARAMETER :: W_SUM_ONE = 1.
      real, pointer :: VOLUME
      integer :: cv_inod_ipha, U_NODK_IPHA, IANISOLIM
      logical :: Have_Temperature_Fields, Have_VolumeFraction_Fields, Have_Components_Fields
      logical :: overlapping
      ! Functions...
      !REAL :: R2NORM, FACE_THETA
      !        ===>  LOGICALS  <===
      character( len = option_path_len ) :: overlapping_path
      LOGICAL :: GETMAT, LIMIT_USE_2ND, &
           D1, D3, DCYL, INTEGRAT_AT_GI, &
           NORMALISE, SUM2ONE, GET_GTHETA, QUAD_OVER_WHOLE_ELE,is_overlapping

      character( len = option_path_len ) :: option_path, option_path2, path_temp, path_volf, &
           path_comp, path_spatial_discretisation


      real, dimension(:), allocatable :: TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, &
              DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT
      INTEGER :: IDUM(1)
      REAL :: RDUM(1),n1,n2,n3
      type( scalar_field ), pointer :: perm
     !Reals to store the irresidual water and irreducible oil values, used in GET_INT_T_DEN
    real ::  s_gc, s_or

    INTEGER :: GLOBAL_FACE

    GLOBAL_FACE=0

      overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) overlapping = .true.

      IDUM = 0
      RDUM = 0.


      ewrite(3,*) 'In CV_ASSEMB_CT'

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

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.

      ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA:', &
           CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, LIMIT_USE_2ND, SECOND_THETA

      ndotq = 0. 

      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

  

      ! Allocate memory for the control volume surface shape functions, etc.

      ALLOCATE( CVNORMX( SCVNGI ))
      ALLOCATE( CVNORMY( SCVNGI ))
      ALLOCATE( CVNORMZ( SCVNGI ))
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
      ALLOCATE( SCVFEWEIGH( SCVNGI ))

      ALLOCATE( SUFEN( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENSLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLX( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLY( U_NLOC, SCVNGI ))
      ALLOCATE( SUFENLZ( U_NLOC, SCVNGI ))

      ALLOCATE( SRA( SCVNGI ))

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
      IF( have_option( '/porous_media/actual_velocity' ) ) THEN
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
      ALLOCATE( TMIN( CV_NONODS * NPHASE) )
      ALLOCATE( TMAX( CV_NONODS * NPHASE) )
      ALLOCATE( T2MIN( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE( T2MAX( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE( DENMIN( CV_NONODS * NPHASE) )
      ALLOCATE( DENMAX( CV_NONODS * NPHASE) )

      ALLOCATE( TMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE( TMAX_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE( T2MIN_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE( T2MAX_2ND_MC( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE( DENMIN_2ND_MC( CV_NONODS * NPHASE) )
      ALLOCATE( DENMAX_2ND_MC( CV_NONODS * NPHASE) )

      ALLOCATE( TMIN_NOD( CV_NONODS * NPHASE) )
      ALLOCATE( TMAX_NOD( CV_NONODS * NPHASE) )
      ALLOCATE( T2MIN_NOD( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE( T2MAX_NOD( CV_NONODS * NPHASE* IGOT_T2) )
      ALLOCATE( DENMIN_NOD( CV_NONODS * NPHASE) )
      ALLOCATE( DENMAX_NOD( CV_NONODS * NPHASE) )

      NCOLGPTS = 0
      COLGPTS = 0
      FINDGPTS = 0

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

      ! Determine FEMT (finite element wise) etc from T (control volume wise)
      ! Also determine the CV mass matrix MASS_CV and centre of the CV's XC_CV,YC_CV,ZC_CV.
      ! This is for projecting to finite element basis functions...
      ALLOCATE( FEMT( CV_NONODS * NPHASE ))
      ALLOCATE( FEMDEN( CV_NONODS * NPHASE ))
      ALLOCATE( FEMT2( CV_NONODS * NPHASE * IGOT_T2 ))
      ALLOCATE( MASS_CV( CV_NONODS ))
      ALLOCATE( MASS_ELE( TOTELE ))
      ALLOCATE( XC_CV( CV_NONODS ))
      ALLOCATE( YC_CV( CV_NONODS ))
      ALLOCATE( ZC_CV( CV_NONODS ))
      ALLOCATE( DTX_ELE( CV_NLOC, NPHASE,TOTELE ))
      ALLOCATE( DTY_ELE( CV_NLOC, NPHASE,TOTELE ))
      ALLOCATE( DTZ_ELE( CV_NLOC, NPHASE,TOTELE ))

      ewrite(3,*)'here2'

      CALL PROJ_CV_TO_FEM_2( state, &
           FEMT, FEMDEN, T, DEN, &
           IGOT_T2, T2, FEMT2, &
           XC_CV, YC_CV, ZC_CV, MASS_CV, MASS_ELE,  &
           NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
           CV_NGI_SHORT, CV_NLOC, CVN_SHORT, CVWEIGHT_SHORT, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, &
           X_NONODS, X, Y, Z, NCOLM, FINDM, COLM, MIDM, &
           0 , MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC )

      MASS_ELE_TRANSP = MASS_ELE

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

      ! For each node, find the largest and smallest value of T and
      ! DENSITY for both the current and previous timestep, out of
      ! the node value and all its surrounding nodes including Dirichlet b.c's.

    CALL SURRO_CV_MINMAX_1Time( TMAX, TMIN, DENMAX, DENMIN, &
           T2MAX, T2MIN, &
           TMAX_2ND_MC, TMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, &
           T2MAX_2ND_MC, T2MIN_2ND_MC, &
           LIMIT_USE_2ND, &
           T, T2, DEN, IGOT_T2, NPHASE, CV_NONODS, size(small_colm), SMALL_FINDRM, SMALL_COLM, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC, SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
           WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
           WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, &
           T2MIN_NOD, T2MAX_NOD, &
           DENMIN_NOD, DENMAX_NOD )

! **********ANISOTROPIC LIMITING...*******************
      IANISOLIM=0
     ! print *,'CV_DISOPT=',CV_DISOPT
     ! stop 72
      IF(CV_DISOPT.GE.5) IANISOLIM=1
! limiting not ready yet for P2 DG:
!      IF(
!      IF( (NDIM==2).AND.(CV_NLOC==6) ) IANISOLIM=0
!      IF( (NDIM==3).AND.(CV_NLOC==10) ) IANISOLIM=0
!IANISOLIM=0
      IF (IANISOLIM==0) THEN
         ALLOCATE(TUPWIND_MAT(1),  DENUPWIND_MAT(1))
         ALLOCATE(T2UPWIND_MAT(1))
         NSMALL_COLM=1
      ELSE ! IF (IANISOLIM==1) THEN
         ! Reduce matrix size...
         NSMALL_COLM=size(small_colm)

         ALLOCATE(TUPWIND_MAT(NSMALL_COLM*NPHASE), &
              DENUPWIND_MAT(NSMALL_COLM*NPHASE))
         ALLOCATE(T2UPWIND_MAT(NSMALL_COLM*NPHASE*IGOT_T2))

         
         CALL CALC_ANISOTROP_LIM_1time(&
              ! Caculate the upwind values stored in matrix form...
              T,DEN,T2, &
              FEMT,FEMDEN,FEMT2, (CV_NONODS.NE.X_NONODS), &
              TUPWIND_MAT, DENUPWIND_MAT, &
              T2UPWIND_MAT, &
              ! Store the upwind element for interpolation and its weights for
              ! faster results...
              IGOT_T2,NPHASE,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
              SMALL_FINDRM,SMALL_CENTRM,SMALL_COLM,NSMALL_COLM, &
              X_NDGLN,X_NONODS,NDIM, &
              X,Y,Z, XC_CV, YC_CV, ZC_CV )


      END IF
! **********...ANISOTROPIC LIMITING*******************

      ALLOCATE( FACE_ELE( NFACE, TOTELE ) ) ; FACE_ELE = 0
      ! Calculate FACE_ELE
      CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
           NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
           CV_SLOCLIST, X_NLOC, X_NDGLN )

      !     =============== DEFINE THETA FOR TIME-STEPPING ===================
      ! Define the type of time integration:
      ! Timopt is 0 if CV_DISOPT is even (theta specified);
      ! Timopt is 1 if CV_DISOPT is odd (non-linear theta scheme)
      TIMOPT = MOD( CV_DISOPT, 2 )

      VTHETA = 1.0
      FTHETA = CV_THETA

     !Get the irreducible water and residual oil before entering the loop
     !this value are used in GET_INT_T_DEN to keep the phases between realistic values
     !Default value has to be the same as in subroutine get_corey_options in multi_eos.F90
      !Default value of  S_GC = 0.1
      call get_option("/material_phase[0]/multiphase_properties/immobile_fraction", &
          s_gc, default=0.1)
      !Default value of    S_OR = 0.3
      call get_option("/material_phase[1]/multiphase_properties/immobile_fraction", &
           s_or, default=0.3)


      ! Now we begin the loop over elements to assemble the advection terms
      ! into the matrix (ACV) and the RHS
      Loop_Elements: DO ELE = 1, TOTELE


         ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
         CALL DETNLXR_INVJAC( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
              CV_NLOC, SCVNGI, &
              SCVFEN, SCVFENLX, SCVFENLY, SCVFENLZ, SCVFEWEIGH, SCVDETWEI, SCVRA, VOLUME, D1, D3, DCYL, &
              SCVFENX_ALL, &
              NDIM, INV_JAC,state, "ali", StorageIndexes(30) )


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
                       X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                       CV_SNLOC, CVFEM_ON_FACE(:,GI), X_SHARE, X_NONODS, ELE, ELE2,  &
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
                        CALL CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                             FACE_ELE, NFACE, CVFEM_ON_FACE( :, GI ), &
                             CV_NONODS, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
                             CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )
                        !EWRITE(3,*)'*****AFTER CALC_SELE SELE,CV_SILOC,CV_SNLOC:',SELE,CV_SILOC,CV_SNLOC
                     END IF
                     INTEGRAT_AT_GI=.NOT.((ELE==ELE2).AND.(SELE==0))
                  END IF

               END IF Conditional_CheckingNeighbourhood

               ! avoid indegrating across the middle of a CV on the boundaries of elements
               Conditional_integration: IF(INTEGRAT_AT_GI) THEN

                  global_face=global_face+1

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

                  Loop_IPHASE: DO IPHASE = 1, NPHASE

                     CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                     CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS

                     NFACE_ITS=1
                     !IF((CV_ELE_TYPE==2).AND.(CV_NONODS==TOTELE*CV_NLOC)) NFACE_ITS=2
                     DO FACE_ITS = 1, NFACE_ITS
                        ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
                        if (igot_t2==1) then
                           CALL GET_INT_VEL_2time( NPHASE, NDOTQNEW,NDOTQ,INCOME,  NDOTQOLD(iphase,global_face), INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T2,T2, FEMT2, FEMT2, DEN,DEN, &
                                U, V, W, NU, NV, NW,NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT2(iphase,global_face),LIMT2(iphase,global_face), FEMDGI,FEMDGI, FEMT2GI, FEMT2GI, UP_WIND_NOD, &
                                T2MIN, T2MAX, T2MIN, T2MAX, T2MIN_NOD, T2MAX_NOD, T2MIN_NOD, T2MAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                T2MIN_2ND_MC,T2MIN_2ND_MC, T2MAX_2ND_MC,T2MAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                T2UPWIND_MAT,T2UPWIND_MAT )
                        else
                            CALL GET_INT_VEL_2time( NPHASE, NDOTQNEW,NDOTQ,INCOME,  NDOTQOLD(iphase,global_face), INCOMEOLD, &
                                HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
                                T,T, FEMT, FEMT, DEN,DEN, &
                                U, V, W, NU, NV, NW,NU, NV, NW, &
                                CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC,WIC_U_BC, &
                                SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
                                UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
                                ONE_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
                                MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
                                FACE_ITS, LIMT(iphase,global_face),LIMT(iphase,global_face), FEMDGI,FEMDGI, FEMTGI, FEMTGI, UP_WIND_NOD, &
                                TMIN, TMAX, TMIN, TMAX, TMIN_NOD, TMAX_NOD, TMIN_NOD, TMAX_NOD, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                TMIN_2ND_MC,TMIN_2ND_MC, TMAX_2ND_MC,TMAX_2ND_MC, LIMIT_USE_2ND,&
                                overlapping, &
                                IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                                TUPWIND_MAT,TUPWIND_MAT )
                         end if
                        !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                        ! Calculate T and DEN on the CV face at quadrature point GI.
                        CALL GET_INT_T_DEN( FVT, FVT2, FVD, LIMD(iphase,global_face), LIMT(iphase,global_face), LIMT2(iphase,global_face), &
                             LIMDT(iphase,global_face), LIMDTT2(iphase,global_face),&
                             FEMDGI, FEMTGI,FEMT2GI, &
                             CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
                             CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, &
                             T, T2, DEN, FEMT, FEMT2, FEMDEN, &
                             TMIN, T2MIN, DENMIN, &
                             TMAX, T2MAX, DENMAX, &
                             SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
                             WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, ONE_PORE, &
                             MASS_CV, TMIN_NOD, TMAX_NOD, &
                             DENMIN_NOD, DENMAX_NOD, &
                             T2MIN_NOD, T2MAX_NOD, IGOT_T2, &
                             TMIN_2ND_MC, T2MIN_2ND_MC, DENMIN_2ND_MC, &
                             TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, &
                             LIMIT_USE_2ND, HDC, NDOTQ, DT, &
                             SCVFENX_ALL(1,:,:), SCVFENX_ALL(2,:,:), SCVFENX_ALL(3,:,:),CVNORMX, CVNORMY, CVNORMZ, &
                             U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             IANISOLIM, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
                             TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT, &
                             !Values to limit the flow when reaching the irreducible  saturation for a phase
                             s_gc, s_or )
                     end DO
                  end DO Loop_IPHASE
               end IF Conditional_integration
            end DO Loop_GCOUNT
         end DO Loop_CV_ILOC
      end DO Loop_Elements

      DEALLOCATE( FACE_ELE )
      deallocate(tupwind_mat,denupwind_mat)
      deallocate(t2upwind_mat) 

      DEALLOCATE( FEMT )
      DEALLOCATE( FEMDEN)
      DEALLOCATE( FEMT2)
      DEALLOCATE( MASS_CV)
      DEALLOCATE( MASS_ELE)
      DEALLOCATE( XC_CV)
      DEALLOCATE( YC_CV)
      DEALLOCATE( ZC_CV)
      DEALLOCATE( DTX_ELE)
      DEALLOCATE( DTY_ELE)
      DEALLOCATE( DTZ_ELE)   

      DEALLOCATE( TMIN )
      DEALLOCATE( TMAX )
      DEALLOCATE( T2MIN )
      DEALLOCATE( T2MAX )
      DEALLOCATE( DENMIN )
      DEALLOCATE( DENMAX )

      DEALLOCATE( TMIN_2ND_MC )
      DEALLOCATE( TMAX_2ND_MC )
      DEALLOCATE( T2MIN_2ND_MC )
      DEALLOCATE( T2MAX_2ND_MC )
      DEALLOCATE( DENMIN_2ND_MC )
      DEALLOCATE( DENMAX_2ND_MC )

      DEALLOCATE( TMIN_NOD )
      DEALLOCATE( TMAX_NOD )
      DEALLOCATE( T2MIN_NOD )
      DEALLOCATE( T2MAX_NOD )
      DEALLOCATE( DENMIN_NOD )
      DEALLOCATE( DENMAX_NOD )

      DEALLOCATE( ONE_PORE)

      DEALLOCATE( CVNORMX)
      DEALLOCATE( CVNORMY )
      DEALLOCATE( CVNORMZ )
      DEALLOCATE( COLGPTS ) !The size of this vector is over-estimated
      DEALLOCATE( FINDGPTS )
      DEALLOCATE( SNDOTQ )
      DEALLOCATE( SNDOTQOLD )
      DEALLOCATE( CV_ON_FACE )
      DEALLOCATE( CVFEM_ON_FACE )
      DEALLOCATE( U_ON_FACE )
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

      DEALLOCATE( SCVFEN)
      DEALLOCATE( SCVFENSLX)
      DEALLOCATE( SCVFENSLY)
      DEALLOCATE( SCVFENLX)
      DEALLOCATE( SCVFENLY)
      DEALLOCATE( SCVFENLZ)
      DEALLOCATE( SCVFEWEIGH)

      DEALLOCATE( SUFEN )
      DEALLOCATE( SUFENSLX )
      DEALLOCATE( SUFENSLY )
      DEALLOCATE( SUFENLX )
      DEALLOCATE( SUFENLY )
      DEALLOCATE( SUFENLZ )

      DEALLOCATE( SRA)

      DEALLOCATE( SBCVN)
      DEALLOCATE( SBCVFEN)
      DEALLOCATE( SBCVFENSLX)
      DEALLOCATE( SBCVFENSLY)
      DEALLOCATE( SBCVFEWEIGH)
      DEALLOCATE( SBCVFENLX)
      DEALLOCATE( SBCVFENLY)
      DEALLOCATE( SBCVFENLZ)
      DEALLOCATE( SBUFEN)
      DEALLOCATE( SBUFENSLX)
      DEALLOCATE( SBUFENSLY)
      DEALLOCATE( SBUFENLX)
      DEALLOCATE( SBUFENLY)
      DEALLOCATE( SBUFENLZ)
      DEALLOCATE( DUMMY_ZERO_NDIM_NDIM)

      DEALLOCATE( CV_SLOC2LOC)
      DEALLOCATE( U_SLOC2LOC)
      DEALLOCATE( CV_SLOCLIST)
      DEALLOCATE( U_SLOCLIST)
      DEALLOCATE( CV_NEILOC)

      DEALLOCATE( SELE_OVERLAP_SCALE )

      DEALLOCATE( UGI_COEF_ELE,  VGI_COEF_ELE,  WGI_COEF_ELE )
      DEALLOCATE( UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2 )
  
      DEALLOCATE( SUM_CV)
      DEALLOCATE( UP_WIND_NOD)



     
    end SUBROUTINE CV_GET_ALL_LIMITED_VALS


    function CV_count_faces( SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM,&
         CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
         CV_ELE_TYPE,  &
         NPHASE,  &
         CV_NLOC, U_NLOC, X_NLOC, &
         CV_NDGLN, X_NDGLN, U_NDGLN, &
         CV_SNLOC, U_SNLOC, STOTEL, CV_SNDGLN, U_SNDGLN, &
         X, Y, Z,&
         MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
         NDIM, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         FINDCMC, COLCMC, NCOLCMC ) result(global_face)

      !  =====================================================================
      !     This subroutine gets the limited face values for a variable
      !
      !     This routine uses a Control Volume (CV) formulation to compute
      !     the advection terms. The general procedure is as follows:
      !
      !     This procedure is implemented by considering the CV to be made up
      !     of a number of sub-control-volumes, which represent the part of
      !     the control volume within a given element.  The assembly of terms
      !     considers each of these sub-CVs in turn, calculating (and limiting)
      !     the flux across sub-CV faces that are external to the CV...
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
      !     CSR_ACV   - Matrix for assembling the advection terms (empty on input)
      !     CV_RHS      - Right-hand side vector for advection-diffusion terms
      !     X,Y,Z    - Node co-ordinates
      !     NU       - Nodal velocity component
      !     T   - advected field values at nodes
      !     CV_DISOPT   - The discretisation/limiting option (see above)
      !     DT       - The time step
      !     ELE_TYP   - Integer flag definining element type
      !
      !     IMPORTANT OUTPUTS:
      !     -----------------
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
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, NDIM, &
           NCOLM, XU_NLOC, NCOLELE, &
           NCOLCMC
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(: ), intent( in ) :: U_SNDGLN
      ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC

      REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      integer, dimension(:), intent(in) :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      !character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation


      ! Local variables
      INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
           WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
           WIC_U_BC_DIRICHLET = 1
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE,LOG_ON_BOUND
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE
      INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, &
           CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
           COLGPTS, CV_SLOC2LOC, U_SLOC2LOC, &
           TMAX_NOD, TMIN_NOD, &
           DENMAX_NOD, DENMIN_NOD, &
           T2MAX_NOD, T2MIN_NOD
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
           UP_WIND_NOD, DU, DV, DW, PERM_ELE
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

      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE,   &
           NCOLGPTS, &
           CV_SILOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, rhs_nodj_ipha,rhs_nodi_ipha,&
           CV_NODI, CV_NODI_IPHA, CV_NODI_JPHA, U_NODK, TIMOPT, &
           JCOUNT_IPHA, IMID_IPHA, &
           NFACE, X_NODI,  &
           CV_INOD, MAT_NODI, FACE_ITS, NFACE_ITS, &
           CVNOD, XNOD, NSMALL_COLM, COUNT2, NOD
      !        ===>  REALS  <===
      REAL :: NDOTQ,& 
           INCOME, INCOMEOLD, HDC, FVT, FVTOLD, FVT2, FVT2OLD, &
           FVD, FVDOLD, FTHETA, VTHETA, &
           FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
           TMID, TOLDMID, &
           DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, BCZERO, ROBIN1, ROBIN2, &
           RSUM, &
           SUM_LIMT, SUM_LIMTOLD, FTHETA_T2, ONE_M_FTHETA_T2OLD, THERM_FTHETA, &
           W_SUM_ONE1, W_SUM_ONE2, NDOTQNEW, VOLUME

      REAL, PARAMETER :: W_SUM_ONE = 1.

      integer :: cv_inod_ipha, U_NODK_IPHA, IANISOLIM
      logical :: Have_Temperature_Fields, Have_VolumeFraction_Fields, Have_Components_Fields
      logical :: overlapping
      ! Functions...
      !REAL :: R2NORM, FACE_THETA
      !        ===>  LOGICALS  <===
      character( len = option_path_len ) :: overlapping_path
      LOGICAL :: QUAD_OVER_WHOLE_ELE,is_overlapping,integrat_at_gi

    INTEGER :: GLOBAL_FACE

    GLOBAL_FACE=0

      overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) overlapping = .true.

      ewrite(3,*) 'In CV_FACE_COUNT'

      is_overlapping = .false.
      call get_option( '/geometry/mesh::VelocityMesh/from_mesh/mesh_shape/element_type', &
           overlapping_path )
      if( trim( overlapping_path ) == 'overlapping' ) is_overlapping = .true.

      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

  

      ! Allocate memory for the control volume surface shape functions, etc.

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

      ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
      ALLOCATE( U_SLOC2LOC( U_SNLOC ))
      ALLOCATE( CV_SLOCLIST( NFACE, CV_SNLOC ))
      ALLOCATE( U_SLOCLIST( NFACE, U_SNLOC ))
      ALLOCATE( CV_NEILOC( CV_NLOC, SCVNGI ))

      ALLOCATE( SELE_OVERLAP_SCALE(CV_NLOC) )


      ewrite(3,*)'here1'


      X_SHARE = .FALSE.

      ! If using the original limiting scheme, the first step is to estimate
      ! the upwind field value from the surrounding nodes

      ! Allocate memory for terms needed by GETGXYZ OR ONVDLIM

      NCOLGPTS = 0
      COLGPTS = 0
      FINDGPTS = 0

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

      ALLOCATE( XC_CV( CV_NONODS ))
      ALLOCATE( YC_CV( CV_NONODS ))
      ALLOCATE( ZC_CV( CV_NONODS ))

! **********...ANISOTROPIC LIMITING*******************

      ALLOCATE( FACE_ELE( NFACE, TOTELE ) ) ; FACE_ELE = 0
      ! Calculate FACE_ELE
      CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
           NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
           CV_SLOCLIST, X_NLOC, X_NDGLN )

      !     =============== DEFINE THETA FOR TIME-STEPPING ===================
      ! Define the type of time integration:
      ! Timopt is 0 if CV_DISOPT is even (theta specified);
      ! Timopt is 1 if CV_DISOPT is odd (non-linear theta scheme)

      ! Now we begin the loop over elements to assemble the advection terms
      ! into the matrix (ACV) and the RHS
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

               Conditional_CheckingNeighbourhood: IF( CV_JLOC == -1 ) THEN

                  ! We are on the boundary or next to another element.  Determine CV_OTHER_LOC
                  ! CVFEM_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
                  ! Look for these nodes on the other elements.
                  CALL FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC,  &
                       MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI,  &
                       X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                       CV_SNLOC, CVFEM_ON_FACE(:,GI), X_SHARE, X_NONODS, ELE, ELE2,  &
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
                        CALL CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                             FACE_ELE, NFACE, CVFEM_ON_FACE( :, GI ), &
                             CV_NONODS, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
                             CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )
                        !EWRITE(3,*)'*****AFTER CALC_SELE SELE,CV_SILOC,CV_SNLOC:',SELE,CV_SILOC,CV_SNLOC
                     END IF
                     INTEGRAT_AT_GI=.NOT.((ELE==ELE2).AND.(SELE==0))
                  END IF

               END IF Conditional_CheckingNeighbourhood

               ! avoid indegrating across the middle of a CV on the boundaries of elements
               Conditional_integration: IF(INTEGRAT_AT_GI) THEN

                  global_face=global_face+1

               end IF Conditional_integration
            end DO Loop_GCOUNT
         end DO Loop_CV_ILOC
      end DO Loop_Elements

      DEALLOCATE( FACE_ELE )
 
      DEALLOCATE( XC_CV)
      DEALLOCATE( YC_CV)
      DEALLOCATE( ZC_CV)

      DEALLOCATE( SCVRA )
      DEALLOCATE( COLGPTS ) !The size of this vector is over-estimated
      DEALLOCATE( FINDGPTS )
      DEALLOCATE( SNDOTQ )
      DEALLOCATE( SNDOTQOLD )
      DEALLOCATE( CV_ON_FACE )
      DEALLOCATE( CVFEM_ON_FACE )
      DEALLOCATE( U_ON_FACE )
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

      DEALLOCATE( SCVFEN)
      DEALLOCATE( SCVFENSLX)
      DEALLOCATE( SCVFENSLY)
      DEALLOCATE( SCVFENLX)
      DEALLOCATE( SCVFENLY)
      DEALLOCATE( SCVFENLZ)
      DEALLOCATE( SCVFENX)
      DEALLOCATE( SCVFENY)
      DEALLOCATE( SCVFENZ)
      DEALLOCATE( SCVFEWEIGH)

      DEALLOCATE( SUFEN )
      DEALLOCATE( SUFENSLX )
      DEALLOCATE( SUFENSLY )
      DEALLOCATE( SUFENLX )
      DEALLOCATE( SUFENLY )
      DEALLOCATE( SUFENLZ )

      DEALLOCATE( SCVDETWEI) 
      DEALLOCATE( SRA)

      DEALLOCATE( SBCVN)
      DEALLOCATE( SBCVFEN)
      DEALLOCATE( SBCVFENSLX)
      DEALLOCATE( SBCVFENSLY)
      DEALLOCATE( SBCVFEWEIGH)
      DEALLOCATE( SBCVFENLX)
      DEALLOCATE( SBCVFENLY)
      DEALLOCATE( SBCVFENLZ)
      DEALLOCATE( SBUFEN)
      DEALLOCATE( SBUFENSLX)
      DEALLOCATE( SBUFENSLY)
      DEALLOCATE( SBUFENLX)
      DEALLOCATE( SBUFENLY)
      DEALLOCATE( SBUFENLZ)

      DEALLOCATE( CV_SLOC2LOC)
      DEALLOCATE( U_SLOC2LOC)
      DEALLOCATE( CV_SLOCLIST)
      DEALLOCATE( U_SLOCLIST)
      DEALLOCATE( CV_NEILOC)

      DEALLOCATE( SELE_OVERLAP_SCALE )

     
    end function CV_COUNT_FACES


    SUBROUTINE FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC,  &
         MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI, &
         X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
         CV_SNLOC, CVFEM_ON_FACE, X_SHARE, X_NONODS, ELE, ELE2,  &
         FINELE, COLELE, NCOLELE) 
      ! We are on the boundary or next to another element. Determine CV_OTHER_LOC,
      ! U_OTHER_LOC. 
      ! CVFEM_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
      ! Look for these nodes on the other elements. 
      ! ELE2=0 also when we are between elements but are trying to integrate across 
      ! the middle of a CV. 
      IMPLICIT NONE
      INTEGER, intent( in ) :: CV_NLOC, CV_NODI, U_NLOC, X_NONODS, X_NLOC, XU_NLOC, &
           &                   CV_SNLOC, MAT_NLOC, ELE, NCOLELE 
      INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN, CV_NDGLN, XU_NDGLN
      LOGICAL, DIMENSION( : ), intent( in ) :: CVFEM_ON_FACE
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE

      INTEGER, DIMENSION( : ), intent( inout ) :: CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC
      LOGICAL, DIMENSION( : ), intent( inout ) :: X_SHARE
      INTEGER, intent( inout ) :: ELE2
      LOGICAL, intent( inout ) :: INTEGRAT_AT_GI
      ! Local variables
      INTEGER :: X_KLOC, X_NODK, X_NODK2, COUNT, ELE3, SUF_COUNT, CV_KLOC, CV_KLOC2, &
           &     U_KLOC, U_KLOC2, CV_NODK, XU_NODK, XU_NODK2, ILEV, JLEV

      !ewrite(3,*) 'In FIND_OTHER_SIDE'

      DO X_KLOC = 1, X_NLOC
         X_NODK = X_NDGLN( ( ELE - 1) * X_NLOC + X_KLOC )
         X_SHARE( X_NODK ) = CVFEM_ON_FACE( X_KLOC )
      END DO

      ELE3 = 0
      DO COUNT = FINELE( ELE ), FINELE( ELE + 1 ) - 1, 1
         ELE2 = COLELE( COUNT )
         SUF_COUNT = 0 ! See if we share the same nodes
         IF ( ELE2 /= ELE ) THEN
            DO X_KLOC = 1, X_NLOC
               X_NODK = X_NDGLN( ( ELE2 - 1 ) * X_NLOC + X_KLOC )
               IF ( X_SHARE( X_NODK ) ) SUF_COUNT = SUF_COUNT + 1
            END DO
         END IF
         IF ( SUF_COUNT == CV_SNLOC ) ELE3 = ELE2
         !ewrite(3,*)'suf_count:', ele, ele2, suf_count, cv_snloc
      END DO

      DO X_KLOC = 1, X_NLOC
         X_NODK = X_NDGLN( ( ELE - 1 ) * X_NLOC + X_KLOC )
         X_SHARE( X_NODK ) = .FALSE.
      END DO

      ELE2 = ELE3
      IF ( ELE2 /= 0 ) THEN 
         ! Is CV_NODI in element ELE2 if yes set ELE2=0 as we don't want to integrate in 
         ! the middle of a CV. 
         DO CV_KLOC = 1, CV_NLOC
            CV_NODK = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_KLOC )
            !ewrite(3,*)'cv_nodi, cv_nodk:', cv_nodi, cv_nodk, ele, ele2
            IF ( CV_NODK == CV_NODI ) INTEGRAT_AT_GI = .FALSE.
         END DO
      END IF

      IF ( ( ELE2 /= 0 ) .AND. INTEGRAT_AT_GI ) THEN ! Determine CV_OTHER_LOC(CV_KLOC)
         CV_OTHER_LOC = 0
         DO CV_KLOC = 1, CV_NLOC
            IF ( CVFEM_ON_FACE( CV_KLOC ) ) THEN ! Find opposite local node
               X_NODK = X_NDGLN( ( ELE - 1 ) * X_NLOC + CV_KLOC )
               DO CV_KLOC2 = 1, CV_NLOC
                  X_NODK2 = X_NDGLN( ( ELE2 - 1 ) * X_NLOC + CV_KLOC2 )
                  IF ( X_NODK2 == X_NODK ) CV_OTHER_LOC( CV_KLOC ) = CV_KLOC2
               END DO
            END IF
         END DO

         U_OTHER_LOC = 0 ! Determine U_OTHER_LOC(U_KLOC)
         IF ( XU_NLOC /= U_NLOC ) THEN ! This is for the overlapping approach
            DO U_KLOC = 1, XU_NLOC ! Find opposite local node
               XU_NODK = XU_NDGLN( ( ELE - 1 ) * XU_NLOC + U_KLOC )
               DO U_KLOC2 = 1, XU_NLOC
                  XU_NODK2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_KLOC2 )
                  ! XU_NLOC==1 is a special case...
                  IF ( ( XU_NODK2 == XU_NODK ) .OR. ( XU_NLOC == 1 ) ) THEN
                     DO ILEV = 1, CV_NLOC
                        JLEV = CV_OTHER_LOC( ILEV )
                        IF ( JLEV /= 0 ) THEN 
                           U_OTHER_LOC( U_KLOC + ( ILEV - 1 ) * XU_NLOC ) = U_KLOC2 + &
                                ( JLEV - 1 ) * XU_NLOC
                        END IF
                     END DO
                  END IF
               END DO
            END DO
         ELSE
            ! Works for non constant and constant (XU_NLOC=1) vel basis functions...
            DO U_KLOC = 1, U_NLOC ! Find opposite local node
               XU_NODK = XU_NDGLN( ( ELE - 1 ) * XU_NLOC + U_KLOC )
               DO U_KLOC2 = 1, U_NLOC
                  XU_NODK2 = XU_NDGLN(( ELE2 - 1 ) * XU_NLOC + U_KLOC2 )
                  IF ( ( XU_NODK2 == XU_NODK ) .OR. ( XU_NLOC == 1 ) ) &
                       U_OTHER_LOC( U_KLOC ) = U_KLOC2
               END DO
            END DO
         END IF

         MAT_OTHER_LOC = CV_OTHER_LOC
      ELSE
         CV_OTHER_LOC = 0
         U_OTHER_LOC = 0
         MAT_OTHER_LOC = 0
      END IF

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
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    REAL, DIMENSION( : ), intent( inout ) :: FEMT, FEMTOLD, FEMDEN, FEMDENOLD
    REAL, DIMENSION( : ), intent( inout ) :: FEMT2, FEMT2OLD
    REAL, DIMENSION( : ), intent( in ) :: T, TOLD, DEN, DENOLD
    REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
    REAL, DIMENSION( : ), intent( inout ) :: MASS_CV, XC_CV,YC_CV,ZC_CV
    REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE
    REAL, DIMENSION( :, : ), intent( in ) :: CVN
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
    INTEGER, DIMENSION( : ), intent( in ) :: FINDM
    INTEGER, DIMENSION( : ), intent( in ) :: COLM
    INTEGER, DIMENSION( : ), intent( in ) :: MIDM

    REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
    INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
    ! Local variables
    REAL, DIMENSION( : ), allocatable :: PSI, FEMPSI, PSI_AVE, PSI_INT
    INTEGER :: NTSOL,NTSOL_AVE,NTSOL_INT,ELE,CV_ILOC,X_INOD,CV_INOD,NL,NFIELD
    CHARACTER(len=100) :: PATH
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

  SUBROUTINE PROJ_CV_TO_FEM_2( state, &
       FEMT, FEMDEN, T, DEN, &
       IGOT_T2,T2, FEMT2, &
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
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    REAL, DIMENSION( : ), intent( inout ) :: FEMT, FEMDEN
    REAL, DIMENSION( : ), intent( inout ) :: FEMT2
    REAL, DIMENSION( : ), intent( in ) :: T,  DEN
    REAL, DIMENSION( : ), intent( in ) :: T2
    REAL, DIMENSION( : ), intent( inout ) :: MASS_CV, XC_CV,YC_CV,ZC_CV
    REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE
    REAL, DIMENSION( :, : ), intent( in ) :: CVN
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
    INTEGER, DIMENSION( : ), intent( in ) :: FINDM
    INTEGER, DIMENSION( : ), intent( in ) :: COLM
    INTEGER, DIMENSION( : ), intent( in ) :: MIDM

    REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
    INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
    ! Local variables
    REAL, DIMENSION( : ), allocatable :: PSI, FEMPSI, PSI_AVE, PSI_INT
    INTEGER :: NTSOL,NTSOL_AVE,NTSOL_INT,ELE,CV_ILOC,X_INOD,CV_INOD,NL,NFIELD
    CHARACTER(len=100) :: PATH
    INTEGER :: velocity_max_iterations, nstates, istate
    LOGICAL :: solve_force_balance, have_component


    ewrite(3,*) 'In PROJ_CV_TO_FEM_2'

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

    NFIELD=2 + IGOT_T2
    ALLOCATE( PSI( NFIELD * CV_NONODS*NPHASE ))
    ALLOCATE( FEMPSI( NFIELD * CV_NONODS*NPHASE ))
    ALLOCATE( PSI_AVE( 3 * CV_NONODS ))
    ALLOCATE( PSI_INT( 1 * CV_NONODS ))

    NTSOL = NFIELD*NPHASE
    NTSOL_AVE = 3
    NTSOL_INT = 1

    NL=CV_NONODS*NPHASE
    PSI( 1 + 0 * NL : NL + 0 * NL )  =      T( 1 : NL )  
    PSI( 1 + 1 * NL : NL + 1 * NL )  =    DEN( 1 : NL ) 
    IF(IGOT_T2==1) THEN
       PSI( 1 + 2 * NL : NL + 2 * NL )  =      T2( 1 : NL )  
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
    FEMDEN( 1 : NL ) = FEMPSI( 1 + 1 * NL : NL + 1 * NL ) 
    IF(IGOT_T2==1) THEN 
       FEMT2( 1 : NL ) = FEMPSI( 1 + 2 * NL : NL + 2 * NL ) 
    ENDIF

    XC_CV( 1 : CV_NONODS ) = PSI_AVE( 1 : CV_NONODS )
    YC_CV( 1 : CV_NONODS ) = PSI_AVE( 1 +CV_NONODS:   2*CV_NONODS )
    ZC_CV( 1 : CV_NONODS ) = PSI_AVE( 1 +2*CV_NONODS: 3*CV_NONODS )
    MASS_CV( 1 : CV_NONODS ) = PSI_INT( 1 : CV_NONODS )

    DEALLOCATE( PSI )
    DEALLOCATE( FEMPSI )
    DEALLOCATE( PSI_AVE )
    DEALLOCATE( PSI_INT )

    ewrite(3,*) 'Leaving PROJ_CV_TO_FEM_2'

    RETURN

  END SUBROUTINE PROJ_CV_TO_FEM_2


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
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    REAL, DIMENSION( : ), intent( inout ) :: FEMPSI
    REAL, DIMENSION( : ), intent( in ) :: PSI
    REAL, DIMENSION( :, : ), intent( in ) :: CVN
    REAL, DIMENSION( : ), intent( inout ) :: PSI_AVE
    REAL, DIMENSION( : ), intent( inout ) :: PSI_INT
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE
    INTEGER, DIMENSION( : ), intent( in ) :: FINDM
    INTEGER, DIMENSION( : ), intent( in ) :: COLM
    INTEGER, DIMENSION( : ), intent( in ) :: MIDM
    CHARACTER(100), intent(in) :: PATH

    REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
    INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
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
       SBCVFEN, SBCVFENSLX, SBCVFENSLY, state, StorName, Indexes)

    ! determine FEMT (finite element wise) etc from T (control volume wise)
    use shape_functions 
    use matrix_operations
    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, NDIM_VEL, NPHASE, &
         U_NONODS, TOTELE, X_NLOC, CV_NGI, U_NLOC, &
         X_NONODS, STOTEL, U_SNLOC, CV_SNLOC, WIC_U_BC_DIRICHLET, SBCVNGI, NFACE
    REAL, DIMENSION( :), intent( in ) :: U, UOLD, V, VOLD, W, WOLD
    REAL, DIMENSION( : , : , : ), intent( inout ) :: DUX_ELE,DUY_ELE,DUZ_ELE, &
         DUOLDX_ELE,DUOLDY_ELE,DUOLDZ_ELE
    REAL, DIMENSION( : , : , : ), intent( inout ) :: DVX_ELE,DVY_ELE,DVZ_ELE, &
         DVOLDX_ELE,DVOLDY_ELE,DVOLDZ_ELE
    REAL, DIMENSION( : , : , : ), intent( inout ) :: DWX_ELE,DWY_ELE,DWZ_ELE, &
         DWOLDX_ELE,DWOLDY_ELE,DWOLDZ_ELE
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  WIC_U_BC
    REAL, DIMENSION( : ), intent( in ) ::  SUF_U_BC,SUF_V_BC,SUF_W_BC
    INTEGER, DIMENSION( :,: ), intent( in ) ::  U_SLOCLIST
    INTEGER, DIMENSION( :,: ), intent( in ) ::  CV_SLOCLIST
    INTEGER, DIMENSION( :,: ), intent( in ) ::  FACE_ELE
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( :, : ), intent( in ) :: CVFEN, CVFENLX, CVFENLY, CVFENLZ
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( :, : ), intent( in ) :: SBUFEN, SBUFENSLX, SBUFENSLY
    REAL, DIMENSION( :, : ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
    REAL, DIMENSION( : ), intent( in ) :: SBWEIGH
    type( state_type ), dimension( : ), intent( inout ) :: state
    character(len=*), intent(in) :: StorName
    integer, dimension(:), intent(inout) :: Indexes

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
         SBCVFEN, SBCVFENSLX, SBCVFENSLY,&
         state, StorName//"U", Indexes(1))

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
            SBCVFEN, SBCVFENSLX, SBCVFENSLY,&
            state, StorName//"V", Indexes(2))
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
            SBCVFEN, SBCVFENSLX, SBCVFENSLY,&
            state,StorName//"W", Indexes(3))
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
       X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY,&
        state, StorName,indx )

    ! determine FEMT (finite element wise) etc from T (control volume wise)
    use shape_functions 
    use shape_functions_NDim
    use matrix_operations 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NDIM, NPHASE, CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
         X_NONODS, STOTEL, CV_SNLOC, X_SNLOC, WIC_T_BC_DIRICHLET, SBCVNGI, NFACE
    REAL, DIMENSION( : ), intent( in ) :: FEMT, FEMTOLD
    REAL, DIMENSION( : , : , : ), intent( inout ) :: DTX_ELE, DTY_ELE, DTZ_ELE, &
         DTOLDX_ELE, DTOLDY_ELE, DTOLDZ_ELE
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  XCV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  WIC_T_BC
    REAL, DIMENSION( : ), intent( in ) ::  SUF_T_BC
    INTEGER, DIMENSION( :,: ), intent( in ) ::  CV_SLOCLIST
    INTEGER, DIMENSION( :,: ), intent( in ) ::  X_SLOCLIST
    INTEGER, DIMENSION( :,: ), intent( in ) ::  FACE_ELE
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( :, : ), intent( in ) :: X_N, X_NLX, X_NLY, X_NLZ
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( :, : ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
    REAL, DIMENSION( :, : ), intent( in ) :: X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY
    REAL, DIMENSION( : ), intent( in ) :: SBWEIGH
    type( state_type ), dimension( : ), intent( inout ) :: state
    character(len=*), intent(in) :: StorName
    integer, intent(inout) :: indx
    ! Local variables
    LOGICAL :: D1, D3, DCYL, APPLYBC
    REAL, pointer, DIMENSION( : ) :: DETWEI, RA
    REAL, pointer, DIMENSION( :, :,: ) :: NX_ALL, X_NX_ALL
    real, pointer :: VOLUME
    REAL, DIMENSION( :, :, : ), allocatable :: MASELE
    REAL, DIMENSION( :, :, : ), allocatable :: VTX_ELE, VTY_ELE, VTZ_ELE, VTOLDX_ELE, VTOLDY_ELE, VTOLDZ_ELE
    REAL, DIMENSION( :, : ), allocatable :: MASS, INV_MASS
    REAL, DIMENSION( : ), allocatable :: VTX, VTY, VTZ, VTOLDX, VTOLDY, VTOLDZ, DTX, DTY, DTZ, DTOLDX, DTOLDY, DTOLDZ
    REAL, DIMENSION( : ), allocatable :: XSL, YSL, ZSL, SNORMXN, SNORMYN, SNORMZN, SDETWE
    INTEGER, DIMENSION( : ), allocatable :: SLOC2LOC, X_SLOC2LOC, ILOC_OTHER_SIDE
    REAL :: NN, NNX, NNY, NNZ, NORMX, NORMY, NORMZ, SAREA, NRBC, RNN, RTBC, &
         VLM_NORX, VLM_NORY, VLM_NORZ
    INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_GI, CV_ILOC2, &
         CV_INOD, CV_INOD2, CV_JLOC2, CV_NODJ2, CV_NODJ2_IPHA, CV_NODJ_IPHA, &
         CV_SILOC, CV_SJLOC, ELE2, IFACE, IPHASE, SELE2, SUF_CV_SJ2, SUF_CV_SJ2_IPHA, &
         X_INOD, SGI, X_SILOC, X_ILOC

    ewrite(3,*)'in DG_DERIVS sbrt'
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
            X_NX_ALL, &
            CV_NLOC, NLX, NLY, NLZ, NX_ALL,&
            state,StorName, indx )

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
                NNX = NNX + N( CV_ILOC, CV_GI )  * NX_ALL(1, CV_JLOC, CV_GI ) * DETWEI( CV_GI )
                NNY = NNY + N( CV_ILOC, CV_GI )  * NX_ALL(2, CV_JLOC, CV_GI ) * DETWEI( CV_GI )
                NNZ = NNZ + N( CV_ILOC, CV_GI )  * NX_ALL(3, CV_JLOC, CV_GI ) * DETWEI( CV_GI )
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



  SUBROUTINE DG_DERIVS_ALL( FEMT, FEMTOLD, &
       DTX_ELE, DTOLDX_ELE, &
       NDIM, NPHASE, NCOMP, CV_NONODS, TOTELE, CV_NDGLN, &
       XCV_NDGLN, X_NLOC, X_NDGLN,&
       CV_NGI, CV_NLOC, CVWEIGHT, &
       N, NLX, NLY, NLZ, &
       X_N, X_NLX, X_NLY, X_NLZ, &
       X_NONODS, X, Y, Z, &
       NFACE, FACE_ELE, CV_SLOCLIST, X_SLOCLIST, STOTEL, CV_SNLOC, X_SNLOC, WIC_T_BC, SUF_T_BC, &
       SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBWEIGH, &
       X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY,&
      state, StorName, StorageIndexes  )

    ! determine FEMT (finite element wise) etc from T (control volume wise)
    use shape_functions
    use shape_functions_NDim
    use matrix_operations
    IMPLICIT NONE

    INTEGER, intent( in ) :: NDIM, NPHASE, NCOMP,CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
         X_NONODS, STOTEL, CV_SNLOC, X_SNLOC, SBCVNGI, NFACE
    REAL, DIMENSION( :, :, : ), intent( in ) :: FEMT, FEMTOLD
    REAL, DIMENSION( :, :, :, :, : ), intent( inout ) :: DTX_ELE, DTOLDX_ELE
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  XCV_NDGLN
    INTEGER, DIMENSION( :, :, : ), intent( in ) ::  WIC_T_BC
    REAL, DIMENSION( :, :, :, : ), intent( in ) ::  SUF_T_BC
    INTEGER, DIMENSION( :, : ), intent( in ) ::  CV_SLOCLIST
    INTEGER, DIMENSION( :, : ), intent( in ) ::  X_SLOCLIST
    INTEGER, DIMENSION( :, : ), intent( in ) ::  FACE_ELE
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( :, : ), intent( in ) :: X_N, X_NLX, X_NLY, X_NLZ
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( :, : ), intent( in ) :: SBCVFEN, SBCVFENSLX, SBCVFENSLY
    REAL, DIMENSION( :, : ), intent( in ) :: X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY
    REAL, DIMENSION( : ), intent( in ) :: SBWEIGH
    type( state_type ), dimension( : ), intent( inout ) :: state
    character(len=*), intent(in) :: StorName
    integer, dimension(:), intent(inout) :: StorageIndexes
    ! Local variables
    REAL, DIMENSION( :, :, : ), ALLOCATABLE :: MASELE
    REAL, DIMENSION( :, :, :, :, : ), ALLOCATABLE :: VTX_ELE, VTOLDX_ELE
    LOGICAL :: D1, D3, DCYL, APPLYBC( NCOMP, NPHASE )
    REAL, pointer, dimension( : ) :: DETWEI, RA
    REAL, pointer, DIMENSION( :,:,:):: NX_ALL
    REAL, pointer, DIMENSION( :, :, : ) :: X_NX_ALL
    REAL, pointer :: VOLUME
    REAL, DIMENSION( CV_NLOC, CV_NLOC )  :: MASS, INV_MASS
    REAL, DIMENSION( NDIM, X_SNLOC ) :: XSL( 3, X_SNLOC ), SNORMXN( NDIM, SBCVNGI ), SDETWE( SBCVNGI )
    INTEGER  :: SLOC2LOC( CV_SNLOC ), X_SLOC2LOC( X_SNLOC ), ILOC_OTHER_SIDE( CV_SNLOC )
    REAL :: NN, NNX( NDIM ), NORMX( 3 ), SAREA, NRBC, RNN, RTBC, VLM_NORX( NDIM )
    INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_GI, CV_ILOC2, &
         CV_INOD, CV_INOD2, CV_JLOC2, CV_NODJ2, CV_NODJ2_IPHA, CV_NODJ_IPHA, &
         CV_SILOC, CV_SJLOC, CV_SJLOC2, ELE2, IFACE, IPHASE, SELE2, SUF_CV_SJ2, SUF_CV_SJ2_IPHA, &
         X_INOD, SGI, X_SILOC, X_ILOC, ICOMP, IDIM
    INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1

    ewrite(3,*)'in DG_DERIVS'

    DTX_ELE = 0.0 ; DTOLDX_ELE = 0.0

    ALLOCATE( MASELE( CV_NLOC, CV_NLOC, TOTELE ) )
    ALLOCATE( VTX_ELE( NDIM, NCOMP, NPHASE, CV_NLOC, TOTELE ) )
    ALLOCATE( VTOLDX_ELE( NDIM, NCOMP, NPHASE, CV_NLOC, TOTELE ) )

    MASELE = 0.0
    VTX_ELE = 0.0

    VTOLDX_ELE = 0.0

    D1 = ( NDIM == 1 )
    D3 = ( NDIM == 3 )
    DCYL = .FALSE.

    Loop_Elements1: DO ELE = 1, TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, &
            X_NLOC, X_NLOC, CV_NGI, &
            X_N, X_NLX, X_NLY, X_NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
            X_NX_ALL, &
            CV_NLOC, NLX, NLY, NLZ, NX_ALL&
            , state,StorName , StorageIndexes(23) )

       Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC

          CV_NODI = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )

          Loop_CV_JLOC: DO CV_JLOC = 1, CV_NLOC

             CV_NODJ = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_JLOC )

             NN  = SUM( N( CV_ILOC, : ) * N(  CV_JLOC, : ) * DETWEI )
             NNX = MATMUL( NX_ALL( :, CV_JLOC, : ), N( CV_ILOC, : )  * DETWEI )

             MASELE( CV_ILOC, CV_JLOC, ELE) = MASELE( CV_ILOC, CV_JLOC, ELE ) + NN

             DO IPHASE = 1, NPHASE
                DO ICOMP = 1, NCOMP
                   VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                        VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                        + NNX (:) * FEMT( ICOMP, IPHASE, CV_NODJ )

                   VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                        VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                        + NNX (:) * FEMTOLD( ICOMP, IPHASE, CV_NODJ )
                END DO
             END DO

          END DO Loop_CV_JLOC

       END DO Loop_CV_ILOC

    END DO Loop_Elements1


    Loop_Elements2: DO ELE = 1, TOTELE

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
               X, Y, Z, X_NONODS, NORMX( 1 ), NORMX( 2 ), NORMX( 3 ) )

          ! Recalculate the normal...
          DO X_SILOC = 1, X_SNLOC
             X_ILOC = X_SLOC2LOC( X_SILOC )
             X_INOD = X_NDGLN(( ELE - 1 ) * X_NLOC + X_ILOC )
             XSL( 1, X_SILOC ) = X( X_INOD )
             XSL( 2, X_SILOC ) = Y( X_INOD )
             XSL( 3, X_SILOC ) = Z( X_INOD )
          END DO

          CALL DGSDETNXLOC2(X_SNLOC, SBCVNGI, &
               XSL( 1, : ), XSL( 2, : ), XSL( 3, : ), &
               X_SBCVFEN, X_SBCVFENSLX, X_SBCVFENSLY, SBWEIGH, SDETWE, SAREA, &
               (NDIM==1), (NDIM==3), (NDIM==-2), &
               SNORMXN( 1, : ), SNORMXN( 2, : ), SNORMXN( 3, : ), &
               NORMX( 1 ), NORMX( 2 ), NORMX( 3 ) )

          IF ( SELE2 == 0 ) THEN
             ! Calculate the nodes on the other side of the face:

             DO CV_SILOC = 1, CV_SNLOC
                CV_ILOC = SLOC2LOC( CV_SILOC )
                CV_INOD = XCV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )

                DO CV_ILOC2 = 1, CV_NLOC
                   CV_INOD2 = XCV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_ILOC2 )

                   IF( CV_INOD2 == CV_INOD ) ILOC_OTHER_SIDE( CV_SILOC ) = CV_ILOC2
                END DO
             END DO

             APPLYBC = (ELE /= ELE2) .AND. (ELE2 /= 0)

          ELSE
             APPLYBC = ( WIC_T_BC( :, :, SELE2 ) == WIC_T_BC_DIRICHLET )
          END IF

          DO CV_SILOC = 1, CV_SNLOC
             CV_ILOC = SLOC2LOC( CV_SILOC )
             DO CV_SJLOC = 1, CV_SNLOC
                CV_JLOC = SLOC2LOC( CV_SJLOC )
                CV_NODJ = CV_NDGLN( (ELE-1)*CV_NLOC + CV_JLOC )
                IF ( SELE2 /= 0 ) THEN
                   CV_JLOC2 = CV_JLOC
                   CV_SJLOC2 = CV_SJLOC
                   CV_NODJ2 = CV_NODJ
                   SUF_CV_SJ2 = CV_SJLOC + CV_SNLOC * ( SELE2 - 1 )
                   NRBC = 0.0
                ELSE
                   CV_JLOC2 = ILOC_OTHER_SIDE( CV_SJLOC )
                   CV_NODJ2 = CV_NDGLN( (ELE2-1)*CV_NLOC + CV_JLOC2 )
                   NRBC = 1.0
                END IF

                ! Have a surface integral on element boundary... 
                VLM_NORX(:) = MATMUL( SNORMXN( :, : ), &
                     SDETWE(:) * SBCVFEN( CV_SILOC, : ) * SBCVFEN( CV_SJLOC, : ) )

                ! add diffusion term...
                DO IPHASE = 1, NPHASE
                   DO ICOMP = 1, NCOMP
                      IF ( APPLYBC( ICOMP, IPHASE ) ) THEN
                  
                         VTX_ELE(:, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                              VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                              - VLM_NORX(:) * 0.5 * ( FEMT( ICOMP, IPHASE, CV_NODJ ) - FEMT( ICOMP, IPHASE, CV_NODJ2 ) * NRBC )
                         VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                              VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                              - VLM_NORX(:) * 0.5 * ( FEMTOLD( ICOMP, IPHASE, CV_NODJ ) - FEMTOLD( ICOMP, IPHASE, CV_NODJ2 ) * NRBC )

                         IF ( SELE2 /= 0 ) THEN
                            IF ( WIC_T_BC( ICOMP, IPHASE, SELE2 ) == WIC_T_BC_DIRICHLET ) THEN

                               RTBC = SUF_T_BC( ICOMP, IPHASE, CV_SJLOC2, SELE2 )

                               VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                    + VLM_NORX(:) * 0.5 * RTBC
                               VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                                    + VLM_NORX(:) * 0.5 * RTBC

                            END IF
                         END IF
                      END IF
                   END DO
                END DO
             END DO

          END DO

       END DO Between_Elements_And_Boundary


    END DO Loop_Elements2


    Loop_Elements3: DO ELE = 1, TOTELE

       MASS( :, : ) = MASELE( :, :, ELE )
       CALL MATDMATINV( MASS, INV_MASS, CV_NLOC )

       FORALL ( IDIM = 1:NDIM, ICOMP = 1:NCOMP, IPHASE = 1:NPHASE )

          DTX_ELE( ICOMP, IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ), VTX_ELE( IDIM, ICOMP, IPHASE, :, ELE ) )
          DTOLDX_ELE( ICOMP, IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ) , VTOLDX_ELE( IDIM, ICOMP, IPHASE, :, ELE ) )

          !DTX_ELE( ICOMP, IDIM, :, IPHASE, ELE ) = MATMUL( INV_MASS( :, : ), VTX_ELE( IDIM, ICOMP, IPHASE, :, ELE ) )
          !DTOLDX_ELE( ICOMP, IDIM, :, IPHASE, ELE ) = MATMUL( INV_MASS( :, : ) , VTOLDX_ELE( IDIM, ICOMP, IPHASE, :, ELE ) )
       END FORALL

    END DO Loop_Elements3

    DEALLOCATE( MASELE, VTX_ELE, VTOLDX_ELE )

    ewrite(3,*)'about to leave DG_DERIVS'

    RETURN

  END SUBROUTINE DG_DERIVS_ALL

 PURE SUBROUTINE ONVDLIM_ALL( TOTELE, &
       TDLIM, TDCEN, INCOME, PELE, PELEOT, &
       ETDNEW, TDMIN, TDMAX, &
       TDMIN_2nd_mc, TDMAX_2nd_mc, FIRORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE, &
       IANISOTROPIC, TUPWIN, TUPWI2 )
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
    INTEGER, intent( in ) :: TOTELE, IANISOTROPIC
    REAL, intent( inout ) :: TDLIM  
    REAL, intent( in ) :: TDCEN, INCOME, COURANT_OR_MINUS_ONE, TUPWIN, TUPWI2 
    INTEGER, intent( in ) :: PELE, PELEOT
    REAL, DIMENSION( : ), intent( in ) :: ETDNEW, TDMIN, TDMAX, TDMIN_2nd_mc, TDMAX_2nd_mc
    LOGICAL, intent( in ) :: FIRORD, NOLIMI, LIMIT_USE_2ND
    REAL :: ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE, TDMAX_PELE, TDMIN_PELEOT, TDMAX_PELEOT, &
         &  TDMIN_2ND_MC_PELE, TDMAX_2ND_MC_PELE, TDMIN_2ND_MC_PELEOT, TDMAX_2ND_MC_PELEOT  

    ETDNEW_PELE = ETDNEW( PELE )
    ETDNEW_PELEOT = ETDNEW( PELEOT )
    IF( PELEOT /= PELE ) THEN
       IF ( IANISOTROPIC /= 1 ) THEN
          TDMIN_PELE   = TDMIN( PELE )
          TDMAX_PELE   = TDMAX( PELE )
          TDMIN_PELEOT = TDMIN( PELEOT )
          TDMAX_PELEOT = TDMAX( PELEOT )
          IF ( LIMIT_USE_2ND ) THEN
             TDMIN_2ND_MC_PELE   = TDMIN_2ND_MC( PELE )
             TDMAX_2ND_MC_PELE   = TDMAX_2ND_MC( PELE )
             TDMIN_2ND_MC_PELEOT = TDMIN_2ND_MC( PELEOT )
             TDMAX_2ND_MC_PELEOT = TDMAX_2ND_MC( PELEOT )
          END IF
       END IF
    END IF

    IF ( IANISOTROPIC == 1 ) THEN ! limit based on 2 largest and 2 minima values of T:
       CALL ONVDLIM_ANO( TOTELE, &
            TDLIM, TDCEN, INCOME, PELE, PELEOT, &
            ETDNEW_PELE, ETDNEW_PELEOT, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE, &
            TUPWIN, TUPWI2 )
    ELSE IF ( LIMIT_USE_2ND ) THEN ! limit based on 2 largest and 2 minima values of T:
       CALL ONVDLIM_2nd( TOTELE, &
            TDLIM, TDCEN, INCOME, PELE, PELEOT, &
            ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE, TDMAX_PELE, TDMIN_PELEOT, TDMAX_PELEOT, &
            TDMIN_2ND_MC_PELE, TDMAX_2ND_MC_PELE, TDMIN_2ND_MC_PELEOT, TDMAX_2ND_MC_PELEOT, &
            FIRORD, NOLIMI )
    ELSE
       CALL ONVDLIM( TOTELE, &
            TDLIM, TDCEN, INCOME, PELE, PELEOT, &
            ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE,TDMAX_PELE, &
            TDMIN_PELEOT, TDMAX_PELEOT, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE )
    END IF

    RETURN

    contains 

        PURE SUBROUTINE ONVDLIM_ANO( TOTELE, &
       TDLIM, TDCEN, INCOME, PELE, PELEOT, &
       ETDNEW_PELE, ETDNEW_PELEOT, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE, &
       TUPWIN2, TUPWI22 )
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
    REAL, intent( in ) :: TDCEN, INCOME, COURANT_OR_MINUS_ONE, TUPWIN2, TUPWI22 
    INTEGER, intent( in ) :: PELE, PELEOT
    REAL, intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
    LOGICAL, intent( in ) :: FIRORD, NOLIMI
    ! Local variables   
    REAL, PARAMETER :: TOLER=1.0E-10
    INTEGER, PARAMETER :: POWER=1
    REAL :: UCIN, UCOU, TDELE, DENOIN, CTILIN, DENOOU, &
         CTILOU, FTILIN, FTILOU, TUPWIN, TUPWI2
    INTEGER :: COUNT

    IF( NOLIMI ) THEN
       TDLIM = TDCEN
       RETURN
    ENDIF

    Conditional_PELEOT: IF( PELEOT /= PELE ) THEN

       TUPWIN = TUPWIN2 ** POWER
       TUPWI2 = TUPWI22 ** POWER

       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE ** POWER
       DENOIN = TDELE - TUPWIN

       IF( ABS( DENOIN ) < TOLER ) DENOIN = SIGN( TOLER, DENOIN )

       UCIN = ETDNEW_PELEOT ** POWER
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT ** POWER
       DENOOU = TDELE - TUPWI2

       IF( ABS( DENOOU ) < TOLER ) DENOOU = SIGN( TOLER, DENOOU )
       UCOU = ETDNEW_PELE ** POWER
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU

    ELSE

       ! Calculate normalisation parameters for incomming velocities 
       TUPWIN = ETDNEW_PELE ** POWER
       UCIN = ETDNEW_PELE ** POWER
       DENOIN = 1.
       CTILIN = 0.

       ! Calculate normalisation parameters for out going velocities 
       TUPWI2 = ETDNEW_PELE ** POWER
       UCOU = ETDNEW_PELE ** POWER
       DENOOU = 1.
       CTILOU = 0.

    END IF Conditional_PELEOT

    Conditional_FIRORD: IF( FIRORD ) THEN ! Velocity is pointing into element

       ! Velocity is going out of element
       TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 

    ELSE

       FTILIN = ( TDCEN - TUPWIN ) / DENOIN
       FTILOU = ( TDCEN - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, COURANT_OR_MINUS_ONE ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, COURANT_OR_MINUS_ONE ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD

    RETURN

  END SUBROUTINE ONVDLIM_ANO


       PURE SUBROUTINE ONVDLIM_2nd( TOTELE, &
       TDLIM, TDCEN, INCOME, PELE, PELEOT, &
       ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE, TDMAX_PELE, TDMIN_PELEOT, TDMAX_PELEOT, &
       TDMIN_2ND_MC_PELE, TDMAX_2ND_MC_PELE, TDMIN_2ND_MC_PELEOT, TDMAX_2ND_MC_PELEOT, &
       FIRORD, NOLIMI )
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
    REAL, intent( in ) :: TDCEN, INCOME, ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE, TDMAX_PELE, &
         &                TDMIN_PELEOT, TDMAX_PELEOT, TDMIN_2ND_MC_PELE, TDMAX_2ND_MC_PELE, &
         &                TDMIN_2ND_MC_PELEOT, TDMAX_2ND_MC_PELEOT
    INTEGER, intent( in ) :: PELE, PELEOT
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

          TC_in  = ETDNEW_PELEOT ! upwind value for incomming (to cell PELE) vel
          TD_in  = ETDNEW_PELE   ! downwind value for incomming vel
          TC_out = ETDNEW_PELE   ! upwind value for outgoing vel
          TD_out = ETDNEW_PELEOT ! downwind value for outgoing vel

          IF( ETDNEW_PELEOT > ETDNEW_PELE ) THEN
             TUPWIN_out = TDMIN_PELE
             TD_MIN_out = TDMIN_2nd_mc_PELE
             TD_MAX_out = TDMAX_PELE !1.E+10

             TUPWIN_in = TDMAX_PELEOT
             TD_MAX_in = TDMAX_2nd_mc_PELEOT
             TD_MIN_in = TDMIN_PELEOT !-1.E+10 used to bound everything
          ELSE
             TUPWIN_out = TDMAX_PELE
             TD_MAX_out = TDMAX_2nd_mc_PELE
             TD_MIN_out = TDMIN_PELE !-1.e+10

             TUPWIN_in  = TDMIN_PELEOT
             TD_MIN_in  = TDMIN_2nd_mc_PELEOT
             TD_MAX_in  = TDMAX_PELEOT !+1.e+10 used to bound everything
          END IF
          ! Calculate normalisation parameters for out going velocities *******
          CALL LIM_TC_2MAX2MIN(TDLIM_out, TC_out,TD_out, TUPWIN_out, TDCEN, TD_MIN_out, TD_MAX_out)

          ! Calculate normalisation parameters for incomming velocities *******
          CALL LIM_TC_2MAX2MIN(TDLIM_in, TC_in,TD_in, TUPWIN_in, TDCEN, TD_MIN_in, TD_MAX_in)

          ! overall TDLIM*******:
          TDLIM = INCOME * TDLIM_in + (1.0-INCOME) * TDLIM_out

       ELSE
          ! Next to boundary...
          TDLIM = ETDNEW_PELE

       ENDIF Conditional_PELEOT

       ! Conditional_FIRORD: IF( .NOT.FIRORD ) THEN...
    ELSE

       ! Velocity is going out of element
       TDLIM = INCOME * ETDNEW_PELEOT + ( 1.0 - INCOME ) * ETDNEW_PELE

    ENDIF Conditional_FIRORD

    RETURN

  END SUBROUTINE ONVDLIM_2nd



 PURE SUBROUTINE ONVDLIM( TOTELE, &
       TDLIM, TDCEN, INCOME, PELE, PELEOT, &
       ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE, TDMAX_PELE, &
       TDMIN_PELEOT, TDMAX_PELEOT, FIRORD, NOLIMI, COURANT_OR_MINUS_ONE )
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
    REAL, intent( in ) :: TDCEN, INCOME, COURANT_OR_MINUS_ONE, ETDNEW_PELE, ETDNEW_PELEOT, &
         &                TDMIN_PELE, TDMAX_PELE, TDMIN_PELEOT, TDMAX_PELEOT
    INTEGER, intent( in ) :: PELE, PELEOT
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

       IF( ETDNEW_PELEOT > ETDNEW_PELE ) THEN
          TUPWIN = TDMAX_PELEOT
          TUPWI2 = TDMIN_PELE
       ELSE
          TUPWIN = TDMIN_PELEOT
          TUPWI2 = TDMAX_PELE 
       END IF

       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE
       DENOIN = TDELE - TUPWIN 

       IF( ABS( DENOIN ) < TOLER ) DENOIN = SIGN( TOLER, DENOIN )

       UCIN = ETDNEW_PELEOT
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT
       DENOOU = TDELE - TUPWI2

       IF( ABS( DENOOU ) < TOLER) DENOOU = SIGN( TOLER, DENOOU )
       UCOU = ETDNEW_PELE
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU

    ELSE

       ! Calculate normalisation parameters for incomming velocities 
       TUPWIN = ETDNEW_PELE
       UCIN = ETDNEW_PELE
       DENOIN = 1.
       CTILIN = 0.

       ! Calculate normalisation parameters for out going velocities 
       TUPWI2 = ETDNEW_PELE
       UCOU = ETDNEW_PELE
       DENOOU = 1.
       CTILOU = 0.

    END IF Conditional_PELEOT

    Conditional_FIRORD: IF( FIRORD ) THEN ! Velocity is pointing into element

       ! Velocity is going out of element
       TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 

    ELSE

       FTILIN = ( TDCEN - TUPWIN ) / DENOIN
       FTILOU = ( TDCEN - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, COURANT_OR_MINUS_ONE ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, COURANT_OR_MINUS_ONE ) * DENOOU )

    END IF Conditional_FIRORD

    RETURN

  END SUBROUTINE ONVDLIM


  END SUBROUTINE ONVDLIM_ALL





 PURE  SUBROUTINE LIM_TC_2MAX2MIN(TDLIM, TC,TD, TUPWIN, TDCEN, TD_MIN,TD_MAX)  
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






 PURE  SUBROUTINE ONVDLIMsqrt( TOTELE, &
       TDLIM, TDCEN, INCOME, PELE, PELEOT, &
       ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE, TDMAX_PELE, TDMIN_PELEOT, TDMAX_PELEOT, FIRORD, NOLIMI )
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
    REAL, intent( in ) :: TDCEN, INCOME, ETDNEW_PELE, ETDNEW_PELEOT, TDMIN_PELE, &
         &                TDMAX_PELE, TDMIN_PELEOT, TDMAX_PELEOT
    INTEGER, intent( in ) :: PELE, PELEOT
    LOGICAL, intent( in ) :: FIRORD, NOLIMI
    ! Local variables   
    REAL, PARAMETER :: TOLER=1.0E-10
    INTEGER, PARAMETER :: POWER=2
    REAL :: UCIN, UCOU, TUPWIN, TUPWI2, TDELE, DENOIN, CTILIN, DENOOU, &
         CTILOU, FTILIN, FTILOU

    IF( NOLIMI ) THEN
       TDLIM = TDCEN
       RETURN
    ENDIF

    Conditional_PELEOT: IF( PELEOT /= PELE ) THEN

       IF( ETDNEW_PELEOT > ETDNEW_PELE ) THEN
          TUPWIN = TDMAX_PELEOT ** POWER
          TUPWI2 = TDMIN_PELE ** POWER
       ELSE
          TUPWIN = TDMIN_PELEOT ** POWER
          TUPWI2 = TDMAX_PELE ** POWER
       ENDIF

       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE ** POWER
       DENOIN = TDELE - TUPWIN

       IF( ABS( DENOIN ) < TOLER ) DENOIN = SIGN( TOLER, DENOIN )

       UCIN = ETDNEW_PELEOT ** POWER
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT ** POWER
       DENOOU = TDELE - TUPWI2

       IF( ABS( DENOOU ) < TOLER) DENOOU = SIGN( TOLER, DENOOU )
       UCOU = ETDNEW_PELE ** POWER
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU

    ELSE

       ! Calculate normalisation parameters for incomming velocities 
       TUPWIN = ETDNEW_PELE ** POWER
       UCIN = ETDNEW_PELE ** POWER
       DENOIN = 1.
       CTILIN = 0.

       ! Calculate normalisation parameters for out going velocities 
       TUPWI2 = ETDNEW_PELE ** POWER
       UCOU = ETDNEW_PELE ** POWER
       DENOOU = 1.
       CTILOU = 0.

    ENDIF Conditional_PELEOT

    Conditional_FIRORD: IF( FIRORD ) THEN ! Velocity is pointing into element
       !     Conditional_FIRORD: IF( .TRUE. ) THEN ! Velocity is pointing into element

       ! Velocity is going out of element
       TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 
       TDLIM = MAX(TDLIM,0.0) ** (1.0/POWER)

    ELSE

       FTILIN = ( MAX(0.0,TDCEN) ** POWER - TUPWIN ) / DENOIN
       FTILOU = ( MAX(0.0,TDCEN) ** POWER - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEWsqrt( FTILIN, CTILIN, -1.0 ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEWsqrt( FTILOU, CTILOU, -1.0 ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD

    RETURN

  END SUBROUTINE ONVDLIMsqrt




  PURE REAL FUNCTION NVDFUNNEWsqrt( UF, UC, COURAT )
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
    REAL, intent(in) :: UC, UF, COURAT
    ! Local variables
    REAL, PARAMETER :: XI = 2. ! XI = 1.
    !REAL, PARAMETER :: XI = 10. ! XI = 1.
    !REAL, PARAMETER :: XI = 5. ! XI = 1.
    REAL :: TILDEUF, MAXUF

    ! For the region 0 < UC < 1 on the NVD, define the limiter
    IF( ( UC > 0.0 ) .AND. ( UC < 1.0 ) ) THEN
       ! For the interface tracking scheme, use Hyper-C (NOTE: at high Courant numbers
       !  limit is defined by XI, not the Hyper-C scheme.)

       IF( COURAT > 0.0 ) THEN
          TILDEUF = MIN( 1.0, max( UC / COURAT, XI * UC ))
       ELSE !For the normal limiting
          MAXUF = MAX( 0.0, UF )
          !MAXUF = MAX( UC, UF )
          TILDEUF = MIN( 1.0, XI * UC, MAXUF )
          !TILDEUF = MIN( 1.0, 1.5* UC, MAXUF )
          !TILDEUF = MIN( 1.0, SQRT(UC), MAXUF )
       ENDIF

    ELSE ! Outside the region 0<UC<1 on the NVD, use first-order upwinding
       TILDEUF = UC
    ENDIF

    NVDFUNNEWsqrt = TILDEUF

  end function nvdfunnewsqrt




  PURE REAL FUNCTION NVDFUNNEW( UF, UC, COURAT )
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
    REAL, intent(in)  :: UC, UF, COURAT
    ! Local variables
    REAL, PARAMETER :: XI = 2., TOLER=1.0E-10
    LOGICAL, PARAMETER :: DOWNWIND_EXTRAP = .TRUE.
    REAL :: TILDEUF, MAXUF

    ! For the region 0 < UC < 1 on the NVD, define the limiter
    IF( ( UC > 0.0 ) .AND. ( UC < 1.0 ) ) THEN
       ! For the interface tracking scheme, use Hyper-C (NOTE: at high Courant numbers
       !  limit is defined by XI, not the Hyper-C scheme.)

       IF( COURAT > TOLER ) THEN
          IF(DOWNWIND_EXTRAP) THEN
             ! new method based on downwind extrapolation...
             !               MAXUF = MAX( 0.0, UF )
             MAXUF = MAX( 0.0, UF, UC )
             !               TILDEUF = MIN( 1.0, UC/ (2.0 * COURAT), MAXUF )
             !               TILDEUF = MIN( 1.0, max(UC/ (2.0 * COURAT), XI * UC ), MAXUF )
             !               TILDEUF = MIN( 1.0, max(UC/ (3.0 * COURAT), XI * UC ), MAXUF )
             !               TILDEUF = MIN( 1.0, UC * 10.0, MAXUF )
             IF ( UF < UC ) THEN
                TILDEUF = UC
             ELSE
                TILDEUF = MIN( 1.0, max(UC / (3.0 * COURAT), XI * UC ) )
             END IF
          ELSE
             !TILDEUF = MIN( 1.0, max( UC / COURAT, XI * UC ))
             ! halve the slope for now...
             TILDEUF = MIN( 1.0, max( UC / (2.0 * COURAT), XI * UC ))
          ENDIF
       ELSE !For the normal limiting
          MAXUF = MAX( 0.0, UF )
          !MAXUF = MAX( UC, UF )
          TILDEUF = MIN( 1.0, XI * UC, MAXUF )
          !TILDEUF = MIN( 1.0, 3.0 * UC, MAXUF )
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
    REAL, DIMENSION( : ), intent( inout ) ::  MASS_CV
    REAL, DIMENSION( : ), intent( in ) :: X
    INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN
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
    INTEGER, DIMENSION( : ), intent( in ) :: XNDGLN
    REAL, DIMENSION( : ), intent( inout ) :: CVDETWEI, CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( :, : ), intent( in ) :: SVN, SVNLX, SVNLY
    REAL, DIMENSION( : ), intent( in ) :: SVWEIGH
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
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
    INTEGER, DIMENSION( : ), intent( in ) :: SNDGLN
    REAL, DIMENSION( :), intent( in ) :: X, Y, Z
    INTEGER, DIMENSION( : ), intent( in ) :: FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: COLM
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
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
    REAL, DIMENSION( :, : ), intent( in ) :: SBUFEN, SBUFENSLX, SBUFENSLY
    REAL, DIMENSION( : ), intent( in ) :: SBWEIGH
    REAL, DIMENSION( : ), intent( inout ) :: SBDETWEI
    REAL, DIMENSION( : ), intent( inout ) :: SNORMXN, SNORMYN, SNORMZN
    REAL, intent( in ) :: SNORMX, SNORMY, SNORMZ
    REAL, DIMENSION( : ), intent( in ) :: X,Y,Z
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
    INTEGER, DIMENSION( : ), intent( in ) :: SNDGLN
    REAL, DIMENSION( :, : ), intent( in ) :: SN, SNLX, SNLY
    REAL, DIMENSION( : ), intent( in ) :: SWEIGH
    REAL, DIMENSION( : ), intent( inout ) :: DETWEI
    REAL, intent( inout ) :: SAREA
    LOGICAL, intent( in ) :: D1, D3, DCYL
    REAL, DIMENSION( : ), intent( inout ) :: NORMXN, NORMYN, NORMZN
    REAL, intent( inout ) :: NORMX, NORMY, NORMZ
    REAL, DIMENSION( : ), intent( in ) :: X,Y,Z
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

    !RSUM = 0.0
    !DO I = 1, NVEC
    !   RSUM = RSUM + VEC( I )**2
    !END DO

    R2NORM = SQRT( SUM( VEC**2 ) )

    RETURN
  END FUNCTION R2NORM



pure real function ptolfun(value)
    ! This function is a tolerance function for strictly positive values used as a denominator.
    ! If the value of VALUE less than 1E-10, then it returns TOLERANCE otherwise VALUE.

    implicit none
    real, intent(in) :: value
    ! Local
    real, parameter :: tolerance = 1.e-10

    ptolfun = max( tolerance, value )

!    if( value > tolerance ) then
!       ptolfun = value
!    else
!       ptolfun = tolerance
!    endif

    return

  end function ptolfun




  pure real function tolfun(value)
    ! This function is a tolerance function for a value which is used as a denominator.
    ! If the absolute value of VALUE less than 1E-10, then it returns SIGN(A,B) i.e. 
    ! the absolute value of A times the sign of B where A is TOLERANCE and B is VALUE.

    implicit none
    real, intent(in) :: value
    ! Local
    real, parameter :: tolerance = 1.e-10

    tolfun = sign( 1.0, value ) * max( tolerance, abs(value) )

!    if( abs( value ) < tolerance ) then
!       tolfun = sign( tolerance, value )
!    else
!       tolfun = value
!    endif

    return

  end function tolfun





  SUBROUTINE BETWEEN_ELE_SOLVE_DIF(UDIFF_SUF_STAB, &
       DIFF_FOR_BETWEEN_U_ELE, DIFF_FOR_BETWEEN_U_ELE2, &
       MAT_ELE, MAT_ELE2, U_SLOC2LOC,U_ILOC_OTHER_SIDE,  &
       SBUFEN,SBCVNGI,U_NLOC,U_SNLOC,NDIM,NPHASE,GOT_OTHER_ELE) 
    ! Calculate the between element diffusion coefficients for stabaization scheme UDIFF_SUF_STAB.
    use matrix_operations
    implicit none     
    !  FAST_AND_SIMP is a good option and is accurate...
    LOGICAL, PARAMETER :: FAST_AND_SIMP=.TRUE.
    !      LOGICAL, PARAMETER :: FAST_AND_SIMP=.FALSE.
    ! If FAST_AND_SIMP use a simple mean to calculate the between element diffusion.
    INTEGER, intent( in ) :: SBCVNGI,U_NLOC,U_SNLOC,NDIM,NPHASE
    LOGICAL, intent( in ) :: GOT_OTHER_ELE
    REAL, DIMENSION( :,:,:,:  ), intent( out ) :: UDIFF_SUF_STAB
    REAL, DIMENSION( :,: ), intent( in ) :: DIFF_FOR_BETWEEN_U_ELE,DIFF_FOR_BETWEEN_U_ELE2
    REAL, DIMENSION( :,: ), intent( in ) :: MAT_ELE, MAT_ELE2
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC,U_ILOC_OTHER_SIDE
    REAL, DIMENSION( :,: ), intent( in ) :: SBUFEN
    ! Local variables...
    REAL, DIMENSION( : , : ), allocatable :: MAT_LOC_2ELES,MAT
    REAL, DIMENSION( : ), allocatable :: VECRHS_2ELES, DIFF
    REAL, DIMENSION( SBCVNGI ) :: DIFF_ADD_STAB
    REAL, DIMENSION( NPHASE, SBCVNGI ) :: DIFF_ADD_STAB2
    INTEGER, DIMENSION( U_NLOC ) :: OTHER_SI2
    INTEGER, DIMENSION( 2*U_NLOC ) :: GLOB_NO
    REAL, DIMENSION( U_SNLOC ) :: DIFF_SUF
    INTEGER :: U_SILOC,U_ILOC,U_ILOC2,U_JLOC2,IGL,JGL,SGI,IPHASE,NLEN,IDIM
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: IPIV
    LOGICAL :: GOTDEC


    ! initialize to zero to set the off diagonal terms to 0.0
    UDIFF_SUF_STAB=0.0

    IF(FAST_AND_SIMP) THEN
       ! use a simple average either side of the interface...
       IF(GOT_OTHER_ELE) THEN
          DIFF_ADD_STAB2( :, : )=0.0
          DO U_SILOC = 1, U_SNLOC
             U_ILOC = U_SLOC2LOC( U_SILOC )
             U_ILOC2= U_ILOC_OTHER_SIDE( U_SILOC ) 
             DO IPHASE=1,NPHASE
                ! Calculate added stabilization diffusion DIFF_ADD_STAB( SGI,IDIM,IPHASE )
                DIFF_ADD_STAB2( IPHASE, : )=DIFF_ADD_STAB2( IPHASE, : )+SBUFEN(U_SILOC,:) &
                     *0.5*(DIFF_FOR_BETWEEN_U_ELE(IPHASE,U_ILOC)+DIFF_FOR_BETWEEN_U_ELE2(IPHASE,U_ILOC2))
             END DO
          END DO
          DIFF_ADD_STAB2=MAX(0.0,DIFF_ADD_STAB2)

          DO IPHASE=1,NPHASE
             DO IDIM=1,NDIM
                UDIFF_SUF_STAB(IDIM,IDIM, IPHASE,:)=DIFF_ADD_STAB2( IPHASE, : )
             END DO
          END DO
       ELSE
          DIFF_ADD_STAB2( :, : )=0.0
          DO U_SILOC = 1, U_SNLOC
             U_ILOC = U_SLOC2LOC( U_SILOC )
             DO IPHASE=1,NPHASE
                ! Calculate added stabilization diffusion DIFF_ADD_STAB2( SGI,IDIM,IPHASE )
                 DIFF_ADD_STAB2( IPHASE, : )=DIFF_ADD_STAB2( IPHASE, : )+SBUFEN(U_SILOC,:)*DIFF_FOR_BETWEEN_U_ELE(IPHASE,U_ILOC)
             END DO
          END DO
          DIFF_ADD_STAB2=MAX(0.0,DIFF_ADD_STAB2)

          DO IPHASE=1,NPHASE
             DO IDIM=1,NDIM
                UDIFF_SUF_STAB(IDIM,IDIM, IPHASE,:)=DIFF_ADD_STAB2( IPHASE, : )
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
       ALLOCATE(IPIV(NLEN,NLEN))
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
          CALL SMLINNGOT( MAT, DIFF, VECRHS_2ELES, NLEN, NLEN, IPIV, GOTDEC)
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
             UDIFF_SUF_STAB(IDIM,IDIM,  IPHASE,1:SBCVNGI )=UDIFF_SUF_STAB(IDIM,IDIM,  IPHASE,1:SBCVNGI )  &
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
    INTEGER, DIMENSION( : ), intent( in ) ::MAT_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::WIC_T_BC
    INTEGER, DIMENSION( : ), intent( in ) ::CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) ::MAT_OTHER_LOC
    REAL, DIMENSION( :,: ), intent( in ) :: SMATFEN
    REAL, DIMENSION( :,: ), intent( in ) :: SCVFEN
    REAL, DIMENSION( :,:,:,: ), intent( in ) :: TDIFFUSION
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
         DIFF_STAND_DIVDX,DIFF_STAND_DIVDX2,COEF
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
       COEF=&
            CVNORMX(GI)*(DIFF_GI(1,1)*CVNORMX(GI)+ DIFF_GI(1,2)*CVNORMY(GI)+DIFF_GI(1,3)*CVNORMZ(GI))&
            +CVNORMY(GI)*(DIFF_GI(2,1)*CVNORMX(GI)+ DIFF_GI(2,2)*CVNORMY(GI)+DIFF_GI(2,3)*CVNORMZ(GI))&
            +CVNORMZ(GI)*(DIFF_GI(3,1)*CVNORMX(GI)+ DIFF_GI(3,2)*CVNORMY(GI)+DIFF_GI(3,3)*CVNORMZ(GI))
       DIFF_STAND_DIVDX= COEF  /HDC
       !         DIFF_STAND_DIVDX=( ABS(CVNORMX(GI))*DIFF_GI(1,1) &
       !              +ABS(CVNORMY(GI))*DIFF_GI(2,2) + ABS(CVNORMZ(GI))*DIFF_GI(3,3) ) /HDC

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
          COEF=&
               CVNORMX(GI)*(DIFF_GI2(1,1)*CVNORMX(GI)+ DIFF_GI2(1,2)*CVNORMY(GI)+DIFF_GI2(1,3)*CVNORMZ(GI))&
               +CVNORMY(GI)*(DIFF_GI2(2,1)*CVNORMX(GI)+ DIFF_GI2(2,2)*CVNORMY(GI)+DIFF_GI2(2,3)*CVNORMZ(GI))&
               +CVNORMZ(GI)*(DIFF_GI2(3,1)*CVNORMX(GI)+ DIFF_GI2(3,2)*CVNORMY(GI)+DIFF_GI2(3,3)*CVNORMZ(GI))
          DIFF_STAND_DIVDX2 = COEF  /HDC
          !            DIFF_STAND_DIVDX2 = ( ABS(CVNORMX(GI))*DIFF_GI2(1,1) &
          !                 +ABS(CVNORMY(GI))*DIFF_GI2(2,2) + ABS(CVNORMZ(GI))*DIFF_GI2(3,3) ) /HDC

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





  SUBROUTINE DIFFUS_CAL_COEFF_STRESS_OR_TENSOR( DIFF_COEF_DIVDX, &
       DIFF_COEFOLD_DIVDX, STRESS_FORM, ZERO_OR_TWO_THIRDS, &
       U_SNLOC, U_NLOC, CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE,  &
       SBCVFEN,SBCVNGI, NDIM_VEL, NDIM, SLOC_UDIFFUSION, SLOC2_UDIFFUSION, DIFF_GI_ADDED, &
       HDC, &
       U_CV_NODJ_IPHA_ALL, U_CV_NODI_IPHA_ALL, &
       UOLD_CV_NODJ_IPHA_ALL, UOLD_CV_NODI_IPHA_ALL, &
       ELE, ELE2, SNORMXN_ALL, &
       SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL,   SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL,  &
       SELE, STOTEL, WIC_U_BC, WIC_U_BC_DIRICHLET )
    ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
    ! based on a non-linear method and a non-oscillating scheme.
! This implements the stress and tensor form of diffusion and calculates a jump conidition. 
! which is in DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
! The coefficient are in N_DOT_DKDU, N_DOT_DKDUOLD. 
! look at the manual DG treatment of viscocity. 
    IMPLICIT NONE
    LOGICAL, intent( in ) :: STRESS_FORM
    INTEGER, intent( in ) :: U_SNLOC, U_NLOC, CV_SNLOC,CV_NLOC, MAT_NLOC, NPHASE,  &
         &                   SBCVNGI, NDIM_VEL, NDIM, ELE, ELE2, &
         &                   SELE, STOTEL, WIC_U_BC_DIRICHLET
    REAL, intent( in ) :: HDC
    REAL, DIMENSION(NDIM_VEL,NPHASE,SBCVNGI), intent( in ) :: U_CV_NODJ_IPHA_ALL, U_CV_NODI_IPHA_ALL, &
                                                          UOLD_CV_NODJ_IPHA_ALL, UOLD_CV_NODI_IPHA_ALL
    REAL, intent( in ) :: ZERO_OR_TWO_THIRDS
    REAL, DIMENSION( NDIM,NPHASE,SBCVNGI ), intent( inout ) :: DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
    INTEGER, DIMENSION( NPHASE,STOTEL ), intent( in ) ::WIC_U_BC
    REAL, DIMENSION( U_SNLOC, SBCVNGI  ), intent( in ) :: SBCVFEN
    REAL, DIMENSION( NDIM,NDIM,NPHASE,U_SNLOC ), intent( in ) :: SLOC_UDIFFUSION, SLOC2_UDIFFUSION
    ! DIFF_GI_ADDED( IDIM, :,:) is for dimension IDIM e.g IDIM=1 corresponds to U 
    ! the rest is for the diffusion tensor. 
    REAL, DIMENSION( NDIM_VEL, NDIM,NDIM, NPHASE, SBCVNGI), intent( in ) :: DIFF_GI_ADDED
    REAL, DIMENSION( NDIM_VEL, NDIM , NPHASE, U_SNLOC ), intent( in ) :: SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL,   SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL
    REAL, DIMENSION( NDIM, SBCVNGI ), intent( in ) :: SNORMXN_ALL

    ! local variables
    !        ===>  REALS  <===
    ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use 
    ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction. 
    ! If SIMPLE_DIFF_CALC then use a simple and fast diffusion calculation.
    LOGICAL, PARAMETER :: SIMPLE_DIFF_CALC2 = .FALSE.
    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.005, DIFF_MAX_FRAC = 200.0

    REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI, DIFF_GI2

    REAL, DIMENSION( : , :, :, :,  : ), allocatable :: DIFF_GI_BOTH
    REAL, DIMENSION( :, :, : ), allocatable :: N_DOT_DKDU, N_DOT_DKDUOLD, N_DOT_DKDU2, N_DOT_DKDUOLD2
    REAL, DIMENSION( :, :, : ), allocatable :: DIFF_STAND_DIVDX_U, DIFF_STAND_DIVDX2_U, &
             DIFF_COEF_DIVDX_U, DIFF_COEFOLD_DIVDX_U
    REAL, DIMENSION( :, : ), allocatable :: IDENT
    REAL, DIMENSION( : ), allocatable :: RZER_DIFF_ALL
    REAL :: COEF
    INTEGER :: CV_KLOC,CV_KLOC2,MAT_KLOC,MAT_KLOC2,MAT_NODK,MAT_NODK2,IDIM,JDIM,CV_SKLOC
    INTEGER :: SGI,IPHASE
    LOGICAL :: ZER_DIFF,SIMPLE_DIFF_CALC

    SIMPLE_DIFF_CALC=SIMPLE_DIFF_CALC2

    ALLOCATE( RZER_DIFF_ALL(NPHASE) )


    ZER_DIFF=.FALSE.
    RZER_DIFF_ALL=1.0
    IF(SELE /= 0) THEN
       ZER_DIFF=.TRUE.
       RZER_DIFF_ALL=0.0
       DO IPHASE=1,NPHASE
          IF(WIC_U_BC( IPHASE, SELE) == WIC_U_BC_DIRICHLET) THEN
             ZER_DIFF=.FALSE.
             RZER_DIFF_ALL(IPHASE)=1.0
          ENDIF
       END DO
    ENDIF

!    ZER_DIFF=.FALSE.
!    IF(SELE /= 0) THEN
!       IF(WIC_U_BC(SELE+(IPHASE-1)*STOTEL) /= WIC_U_BC_DIRICHLET) THEN
!          ZER_DIFF=.TRUE.
!       ELSE
!          SIMPLE_DIFF_CALC=.FALSE.
!       ENDIF
!    ENDIF



    Cond_ZerDiff: IF(ZER_DIFF) THEN

       DIFF_COEF_DIVDX    = 0.0 
       DIFF_COEFOLD_DIVDX = 0.0

    ELSE


       ALLOCATE( DIFF_GI_BOTH(NDIM_VEL, NDIM,NDIM, NPHASE,SBCVNGI) )

       ALLOCATE( N_DOT_DKDU( NDIM_VEL,NPHASE,SBCVNGI )  )
       ALLOCATE( N_DOT_DKDUOLD( NDIM_VEL,NPHASE,SBCVNGI )  )
       ALLOCATE( N_DOT_DKDU2( NDIM_VEL,NPHASE,SBCVNGI )  )
       ALLOCATE( N_DOT_DKDUOLD2( NDIM_VEL,NPHASE,SBCVNGI )  )

       ALLOCATE( DIFF_STAND_DIVDX_U( NDIM_VEL,NPHASE,SBCVNGI )  )
       ALLOCATE( DIFF_STAND_DIVDX2_U( NDIM_VEL,NPHASE,SBCVNGI )  )

       ALLOCATE( DIFF_COEF_DIVDX_U( NDIM_VEL,NPHASE,SBCVNGI )  )
       ALLOCATE( DIFF_COEFOLD_DIVDX_U( NDIM_VEL,NPHASE,SBCVNGI )  )



       IF(SIMPLE_DIFF_CALC) THEN ! The simplest method we can think of... 

          ALLOCATE( DIFF_GI(NDIM,NDIM,NPHASE,SBCVNGI) )
          ALLOCATE( DIFF_GI2(NDIM,NDIM,NPHASE,SBCVNGI) )

          ALLOCATE( IDENT(NDIM,NDIM) )
          IDENT=0.0
          DO IDIM=1,NDIM
            IDENT(IDIM,IDIM)=1.0
          END DO

          DIFF_GI = 0.0
          DO CV_SKLOC = 1, CV_SNLOC
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   DIFF_GI( 1:NDIM , 1:NDIM, IPHASE, SGI ) = DIFF_GI( 1:NDIM , 1:NDIM, IPHASE, SGI  ) &
                     + SBCVFEN(CV_SKLOC,SGI) * SLOC_UDIFFUSION( 1:NDIM , 1:NDIM , IPHASE, CV_SKLOC )
                END DO
             END DO
          END DO
          DIFF_GI=MAX(0.0, DIFF_GI) 

          Conditional_MAT_DISOPT_ELE2_2: IF( ( ELE2 /= 0 ).AND.( ELE2 /= ELE) ) THEN
             DIFF_GI2 = 0.0
             DO CV_SKLOC = 1, CV_SNLOC
                DO SGI=1,SBCVNGI
                   DO IPHASE=1, NPHASE
                      DIFF_GI2( 1:NDIM, 1:NDIM, IPHASE, SGI )=DIFF_GI2( 1:NDIM, 1:NDIM, IPHASE, SGI ) +SBCVFEN(CV_SKLOC,SGI) &
                        *SLOC2_UDIFFUSION(1:NDIM, 1:NDIM ,IPHASE, CV_SKLOC)
                   END DO
                END DO
             END DO
             DIFF_GI2=MAX(0.0, DIFF_GI2) 
             DIFF_GI=0.5*(DIFF_GI+DIFF_GI2)
          ENDIF Conditional_MAT_DISOPT_ELE2_2

          IF(STRESS_FORM) THEN
             
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   DO IDIM=1, NDIM_VEL
                      DIFF_COEF_DIVDX(IDIM,IPHASE,SGI)=8.*( SUM( (1.+IDENT(IDIM,:))*SNORMXN_ALL(:,SGI)**2*DIFF_GI(IDIM,:,IPHASE,SGI) ) &
                        +DIFF_GI_ADDED(IDIM, 1,1, IPHASE,SGI) ) /HDC
!                        +DIFF_GI_ADDED(IDIM, IPHASE,SGI,1,1) ) /HDC
                   END DO
                END DO
             END DO
          ELSE
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   COEF=0.0
                   DO IDIM=1,NDIM
                      COEF=COEF + SNORMXN_ALL(IDIM,SGI)*( SUM( DIFF_GI(IDIM,:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )  )
                   END DO
                   DIFF_COEF_DIVDX(:,IPHASE,SGI)=8.*( COEF + DIFF_GI_ADDED(:, 1,1, IPHASE,SGI) ) /HDC
                END DO
             END DO
          ENDIF

          DIFF_COEFOLD_DIVDX=DIFF_COEF_DIVDX

          ! END OF IF(SIMPLE_DIFF_CALC) THEN...
       ELSE



! Calculate DIFF_COEF_DIVDX, N_DOT_DKDU, N_DOT_DKDUOLD
          CALL FOR_TENS_DERIVS_NDOTS(DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD,  &
                 DIFF_GI_ADDED, SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL, SLOC_UDIFFUSION, &
                 NDIM_VEL, NDIM, NPHASE, U_SNLOC, SBCVNGI, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, STRESS_FORM )



          Conditional_MAT_DISOPT_ELE2: IF( ( ELE2 /= 0 ).AND.( ELE2 /= ELE) ) THEN



! Calculate DIFF_COEF_DIVDX, N_DOT_DKDU, N_DOT_DKDUOLD
             CALL FOR_TENS_DERIVS_NDOTS(DIFF_STAND_DIVDX2_U, N_DOT_DKDU2, N_DOT_DKDUOLD2,  &  
                    DIFF_GI_ADDED, SLOC2_DUX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL, SLOC2_UDIFFUSION, &
                    NDIM_VEL, NDIM, NPHASE, U_SNLOC, SBCVNGI, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, STRESS_FORM )




               N_DOT_DKDU = 0.5*( N_DOT_DKDU + N_DOT_DKDU2 )
               N_DOT_DKDUOLD= 0.5*( N_DOT_DKDUOLD + N_DOT_DKDUOLD2 )

            ! This is the minimum diffusion...
               DIFF_STAND_DIVDX_U    = 0.5*( DIFF_STAND_DIVDX_U + DIFF_STAND_DIVDX2_U ) 

          ENDIF Conditional_MAT_DISOPT_ELE2



          DO SGI=1,SBCVNGI
             DO IPHASE=1, NPHASE
                DO IDIM=1,NDIM_VEL

                   DIFF_COEF_DIVDX_U(IDIM,IPHASE,SGI)    = N_DOT_DKDU(IDIM,IPHASE,SGI) / &
                      TOLFUN( U_CV_NODJ_IPHA_ALL(IDIM,IPHASE,SGI)  - U_CV_NODI_IPHA_ALL(IDIM,IPHASE,SGI) )  
                   DIFF_COEFOLD_DIVDX_U(IDIM,IPHASE,SGI) = N_DOT_DKDUOLD(IDIM,IPHASE,SGI) /  &
                      TOLFUN( UOLD_CV_NODJ_IPHA_ALL(IDIM,IPHASE,SGI)  - UOLD_CV_NODI_IPHA_ALL(IDIM,IPHASE,SGI) )  

                END DO
             END DO
          END DO

         ! Make sure the diffusion has an lower bound...       
            DIFF_COEF_DIVDX_U    = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEF_DIVDX_U )
            DIFF_COEFOLD_DIVDX_U = MAX( DIFF_MIN_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEFOLD_DIVDX_U )
         ! Make sure the diffusion has an upper bound...       
            DIFF_COEF_DIVDX_U    = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEF_DIVDX_U )
            DIFF_COEFOLD_DIVDX_U = MIN( DIFF_MAX_FRAC*DIFF_STAND_DIVDX_U, DIFF_COEFOLD_DIVDX_U )

! Redfine for output...
            DIFF_COEF_DIVDX   = DIFF_COEF_DIVDX_U
            DIFF_COEFOLD_DIVDX = DIFF_COEFOLD_DIVDX_U

          ! END OF IF(SIMPLE_DIFF_CALC) THEN ELSE...
       ENDIF

    END IF Cond_ZerDiff

! 
! Zero if we are on boundary and applying Dirichlet b.c's
      DO IPHASE=1, NPHASE
          DIFF_COEF_DIVDX(:,IPHASE,:)    =  RZER_DIFF_ALL(IPHASE)*DIFF_COEF_DIVDX(:,IPHASE,:)
          DIFF_COEFOLD_DIVDX(:,IPHASE,:) =  RZER_DIFF_ALL(IPHASE)*DIFF_COEFOLD_DIVDX(:,IPHASE,:)
      END DO


    RETURN            

  END SUBROUTINE DIFFUS_CAL_COEFF_STRESS_OR_TENSOR





        SUBROUTINE FOR_TENS_DERIVS_NDOTS( DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD,  &
                 DIFF_GI_ADDED, SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL, SLOC_UDIFFUSION, &
                 NDIM_VEL, NDIM, NPHASE, U_SNLOC, SBCVNGI, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, STRESS_FORM )

! Calculate DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD
! This implements the stress and tensor form of diffusion and calculates a jump conidition. 
! DIFF_STAND_DIVDX_U is the minimal amount of diffusion. 
! The coefficient are in N_DOT_DKDU, N_DOT_DKDUOLD. 
! look at the manual DG treatment of viscocity. 
    IMPLICIT NONE
      INTEGER, intent( in )  :: NDIM_VEL, NDIM, NPHASE, U_SNLOC, SBCVNGI
      REAL, intent( in )  :: HDC, ZERO_OR_TWO_THIRDS
      LOGICAL, intent( in )  :: STRESS_FORM
    REAL, DIMENSION( NDIM,NPHASE,SBCVNGI ), intent( inout ) :: DIFF_STAND_DIVDX_U
    REAL, DIMENSION( NDIM_VEL,NPHASE,SBCVNGI ), intent( inout ) :: N_DOT_DKDU, N_DOT_DKDUOLD
    ! DIFF_GI_ADDED( IDIM, :,:) is for dimension IDIM e.g IDIM=1 corresponds to U 
    ! the rest is for the diffusion tensor. 
    REAL, DIMENSION( NDIM_VEL, NDIM,NDIM, NPHASE, SBCVNGI), intent( in ) :: DIFF_GI_ADDED
    REAL, DIMENSION( NDIM_VEL, NDIM , NPHASE, U_SNLOC ), intent( in ) :: SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL 
    REAL, DIMENSION( U_SNLOC, SBCVNGI  ), intent( in ) :: SBCVFEN
    REAL, DIMENSION( NDIM,NDIM,NPHASE,U_SNLOC ), intent( in ) :: SLOC_UDIFFUSION
    REAL, DIMENSION( NDIM, SBCVNGI ), intent( in ) :: SNORMXN_ALL

    ! local variables
    REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI, STRESS_INDEX, STRESS_INDEXOLD
!    REAL, DIMENSION( : , :, :, :,  : ), allocatable :: DIFF_GI_BOTH
    REAL, DIMENSION( :, :, :, : ), allocatable :: DUDX_ALL_GI, DUOLDDX_ALL_GI
    REAL, DIMENSION( :, : ), allocatable :: IDENT
    REAL :: COEF, DIVU, DIVUOLD
    INTEGER :: U_KLOC,U_KLOC2,MAT_KLOC,MAT_KLOC2,IDIM,JDIM,IDIM_VEL,U_SKLOC
    INTEGER :: SGI,IPHASE
    LOGICAL :: ZER_DIFF,SIMPLE_DIFF_CALC


       ALLOCATE( DIFF_GI(NDIM,NDIM,NPHASE,SBCVNGI) )
       ALLOCATE( STRESS_INDEX(NDIM,NDIM,NPHASE,SBCVNGI) )
       ALLOCATE( STRESS_INDEXOLD(NDIM,NDIM,NPHASE,SBCVNGI) )
!       ALLOCATE( DIFF_GI_BOTH(NDIM_VEL, NDIM,NDIM, NPHASE,SBCVNGI) )

       ALLOCATE( DUDX_ALL_GI( NDIM_VEL,NDIM,NPHASE,SBCVNGI )  )
       ALLOCATE( DUOLDDX_ALL_GI( NDIM_VEL,NDIM,NPHASE,SBCVNGI )  )

       ALLOCATE( IDENT(NDIM,NDIM) )


          IDENT=0.0
          DO IDIM=1,NDIM
            IDENT(IDIM,IDIM)=1.0
          END DO


          DUDX_ALL_GI = 0.0
          DUOLDDX_ALL_GI = 0.0

          DO U_SKLOC = 1, U_SNLOC
             DO SGI=1,SBCVNGI
             ! U, V & W: 
                   DUDX_ALL_GI(:,:,:,SGI)    = DUDX_ALL_GI(:,:,:,SGI)    + SBCVFEN(U_SKLOC,SGI) * SLOC_DUX_ELE_ALL(:,:,:,U_SKLOC)
                   DUOLDDX_ALL_GI(:,:,:,SGI) = DUOLDDX_ALL_GI(:,:,:,SGI) + SBCVFEN(U_SKLOC,SGI) * SLOC_DUOLDX_ELE_ALL(:,:,:,U_SKLOC)
             END DO
          END DO

          DIFF_GI = 0.0
          DO U_SKLOC = 1, U_SNLOC
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   DIFF_GI( 1:NDIM , 1:NDIM, IPHASE,SGI ) = DIFF_GI( 1:NDIM , 1:NDIM, IPHASE,SGI ) &
                  + SBCVFEN(U_SKLOC,SGI) * SLOC_UDIFFUSION( 1:NDIM , 1:NDIM , IPHASE, U_SKLOC )
                END DO
             END DO
          END DO
          DIFF_GI=MAX(0.0, DIFF_GI) 


          ! U:
!          DIFF_GI_BOTH=DIFF_GI_ADDED
!          IF(.NOT.STRESS_FORM) THEN
!             DO IDIM=1,NDIM_VEL
!                DIFF_GI_BOTH(IDIM,:,:,:,:) = DIFF_GI_BOTH(IDIM,:,:,:,:) + DIFF_GI(:,:,:,:)
!             END DO
!          ENDIF



          IF(STRESS_FORM) THEN 
! FOR STRESS FORM...
! BUT 1st tensor form for added diffusion from stabilization say...
             N_DOT_DKDU=0.0
             N_DOT_DKDUOLD=0.0
             DIFF_STAND_DIVDX_U=0.0
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   DO IDIM_VEL=1,NDIM_VEL
                      DO IDIM=1,NDIM
! tensor form...
                         N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   =  N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   &
                          +  SNORMXN_ALL(IDIM,SGI)*SUM( DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI) * DUDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) ) 
! tensor form...
                         N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)= N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)  &
                          +  SNORMXN_ALL(IDIM,SGI)*SUM( DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI) * DUOLDDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) ) 
! for minimal amount of diffusion calc...
                         DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   =  DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   &
                          + SNORMXN_ALL(IDIM,SGI)*SUM( DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI) * SNORMXN_ALL(:,SGI) )  /HDC

                     END DO
                  END DO
                END DO
             END DO

! stress form needs to add this...
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE

                   DIVU=0.0
                   DIVUOLD=0.0
                   DO IDIM=1,NDIM
                      DIVU=DIVU+DUDX_ALL_GI(IDIM,IDIM,IPHASE,SGI)
                      DIVUOLD=DIVUOLD+DUOLDDX_ALL_GI(IDIM,IDIM,IPHASE,SGI)
                   END DO

                   DO IDIM_VEL=1,NDIM_VEL
! Stress form...
                         N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   =  N_DOT_DKDU(IDIM_VEL,IPHASE,SGI) &
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI(IDIM_VEL,:,IPHASE,SGI)*DUDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )  & 
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI(IDIM_VEL,:,IPHASE,SGI)*DUDX_ALL_GI(:,IDIM_VEL,IPHASE,SGI) ) & 
! stress form addition...  
                                                           - ZERO_OR_TWO_THIRDS*SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*DIVU

! Stress form...
                         N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)= N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI) &
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI(IDIM_VEL,:,IPHASE,SGI)*DUOLDDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )  & 
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI(IDIM_VEL,:,IPHASE,SGI)*DUOLDDX_ALL_GI(:,IDIM_VEL,IPHASE,SGI) ) & 
! stress form addition...  
                                                           - ZERO_OR_TWO_THIRDS*SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*DIVUOLD
                                   
! This is for the minimum & max. diffusion...
                         DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   =  DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI) &
                                                           + (   SUM( SNORMXN_ALL(:,SGI)*DIFF_GI(IDIM_VEL,:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )  & 
                                                               + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI(IDIM_VEL,:,IPHASE,SGI)*SNORMXN_ALL(IDIM_VEL,SGI) )    )/HDC

                   END DO
                END DO
             END DO


          ELSE  ! IF(STRESS_FORM) THEN ELSE
! tensor form...
! tensor form for added diffusion from stabilization as well...
             N_DOT_DKDU=0.0
             N_DOT_DKDUOLD=0.0
             DIFF_STAND_DIVDX_U=0.0
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   DO IDIM_VEL=1,NDIM_VEL
                      DO IDIM=1,NDIM
! tensor form...
                         N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   =  N_DOT_DKDU(IDIM_VEL,IPHASE,SGI)   &
                          +  SNORMXN_ALL(IDIM,SGI)*SUM( (DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI)+DIFF_GI(IDIM,:,IPHASE,SGI)) * DUDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) ) 
! tensor form...
                         N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)= N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)  &
                          +  SNORMXN_ALL(IDIM,SGI)*SUM( (DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI)+DIFF_GI(IDIM,:,IPHASE,SGI)) * DUOLDDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) ) 
! This is for the minimum & max. diffusion...
                         DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   =  DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   &
                          +  SNORMXN_ALL(IDIM,SGI)*SUM( (DIFF_GI_ADDED(IDIM_VEL,IDIM,:,IPHASE,SGI)+DIFF_GI(IDIM,:,IPHASE,SGI)) * SNORMXN_ALL(:,SGI) )   /HDC

                     END DO
                  END DO
                END DO
             END DO


             ! ENDOF IF(STRESS_FORM) THEN ELSE...
          ENDIF
! just in case...
! the factor of 8 is there to take into account that HD is measured between centres of elements...
          DIFF_STAND_DIVDX_U=abs( 8.*DIFF_STAND_DIVDX_U )

        RETURN

        END SUBROUTINE FOR_TENS_DERIVS_NDOTS





    SUBROUTINE CALC_STRESS_TEN(STRESS_IJ, ZERO_OR_TWO_THIRDS, NDIM,    &
                 UFENX_ILOC, UFENX_JLOC,  TEN_XX )
! determine stress form of viscocity...
      IMPLICIT NONE
      INTEGER, intent( in )  :: NDIM
      REAL, DIMENSION( NDIM, NDIM  ), intent( inOUT ) :: STRESS_IJ
      REAL, DIMENSION( NDIM ), intent( in ) :: UFENX_ILOC
      REAL, DIMENSION( NDIM,NDIM ), intent( in ) :: TEN_XX
      REAL, DIMENSION( NDIM ), intent( in ) :: UFENX_JLOC
      REAL, intent( in ) :: ZERO_OR_TWO_THIRDS
! Local variables...
      REAL :: FEN_TEN_XX(NDIM,NDIM)
      INTEGER :: IDIM,JDIM


         DO IDIM=1,NDIM
            FEN_TEN_XX(IDIM,:)=UFENX_ILOC(:) * TEN_XX(IDIM,:)
         END DO

         DO IDIM=1,NDIM
               STRESS_IJ( IDIM,IDIM ) = STRESS_IJ( IDIM,IDIM ) + SUM( FEN_TEN_XX(IDIM,:) * UFENX_JLOC(:) ) 
         END DO

         DO IDIM=1,NDIM
            DO JDIM=1,NDIM
               STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) + FEN_TEN_XX(IDIM,JDIM) * UFENX_JLOC(JDIM) 
            END DO
         END DO

         DO IDIM=1,NDIM
            DO JDIM=1,NDIM
               STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) &
                      - ZERO_OR_TWO_THIRDS * FEN_TEN_XX(IDIM,IDIM) * UFENX_JLOC(JDIM) 
            END DO
         END DO

      RETURN            

    END SUBROUTINE CALC_STRESS_TEN






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
    INTEGER, DIMENSION( : ), intent( in ) ::MAT_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::WIC_T_BC
    INTEGER, DIMENSION( : ), intent( in ) ::MAT_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) ::CV_SLOC2LOC
    REAL, DIMENSION( :, :  ), intent( in ) :: SBCVFEN
    REAL, DIMENSION( :, :, :, :  ), intent( in ) :: TDIFFUSION
    REAL, DIMENSION( :, : ), intent( in ) :: DIFF_GI_ADDED
    REAL, DIMENSION(CV_NLOC, NPHASE, TOTELE), intent( in ) :: DTX_ELE,DTY_ELE,DTZ_ELE, &
         DTOLDX_ELE,DTOLDY_ELE,DTOLDZ_ELE
    REAL, DIMENSION( : ), intent( in ) :: SNORMXN,SNORMYN,SNORMZN

    ! local variables
    !        ===>  REALS  <===
    ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use 
    ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction. 
    !REAL, PARAMETER :: DIFF_MIN_FRAC = 0.05, DIFF_MAX_FRAC = 20.0
    REAL, PARAMETER :: DIFF_MIN_FRAC = 1., DIFF_MAX_FRAC = 1.
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





  SUBROUTINE GET_INT_VEL( NPHASE, NDOTQNEW, NDOTQ,INCOME, &
       HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
       T, FEMT, DEN, &
       U, V, W, NU, NV, NW, &
       CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
       SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
       UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
       VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
       MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
       FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
       TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       TMIN_2ND_MC, TMAX_2ND_MC,  LIMIT_USE_2ND,&
       is_overlapping,  &
       IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
       TUPWIND_MAT )
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, WIC_U_BC_DIRICHLET, CV_ELE_TYPE, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, &
         NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, LIMT, FEMDGI, FEMTGI
    REAL, intent( inout ) :: NDOTQNEW,NDOTQ, INCOME
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_U_BC
    INTEGER, DIMENSION( : ), intent( in ) :: TMIN_NOD, TMAX_NOD 
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :  ), intent( in ) :: T, FEMT,DEN
    REAL, DIMENSION( :  ), intent( in ) :: U, V, W, NU, NV, NW 
    REAL, DIMENSION( : ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
    REAL, DIMENSION( : ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( :  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    REAL, DIMENSION( :  ), intent( in ) :: MASS_CV
    REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( :  ), intent( inout ) :: UP_WIND_NOD
    REAL, DIMENSION( :  ), intent( inout ) :: TMIN, TMAX
    logical, intent( in ) :: LIMIT_USE_2ND
    REAL, DIMENSION( : ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC

    INTEGER, intent( in ) :: IANISOTROPIC
    INTEGER, intent( in ) :: NSMALL_COLM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT
    ! local variables
    logical, INTENT (IN) :: is_overlapping   
    INTEGER :: U_NLOC_LEV,U_KLOC_LEV,U_KLOC,U_NODK_IPHA, U_KLOC2, U_NODK2_IPHA

    IF( is_overlapping ) THEN
       ! For overlapping basis function approach.
       CALL GET_INT_VEL_OVERLAP( NPHASE, NDOTQ, INCOME, &
            HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            T, FEMT, DEN, NU, NV, NW, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
            SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
            TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            TMIN_2ND_MC, TMAX_2ND_MC,  LIMIT_USE_2ND, &
            IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
            TUPWIND_MAT)

    ELSE
       CALL GET_INT_VEL_ORIG( NPHASE, NDOTQ, INCOME, &
            HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            T, FEMT, DEN, NU, NV, NW, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
            WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
            TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            TMIN_2ND_MC, TMAX_2ND_MC,  LIMIT_USE_2ND )
    END IF 

    ! Calculate NDOTQNEW from NDOTQ
    NDOTQNEW = NDOTQ ! initialize it like this so that it contains the b.c's
    U_NLOC_LEV = U_NLOC / CV_NLOC
    DO U_KLOC = 1, U_NLOC
       U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) + (IPHASE-1)*U_NONODS
       NDOTQNEW=NDOTQNEW &
            + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE(U_KLOC) * ( U( U_NODK_IPHA ) - NU( U_NODK_IPHA ) ) * CVNORMX(GI) &
            + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE(U_KLOC) * ( V( U_NODK_IPHA ) - NV( U_NODK_IPHA ) ) * CVNORMY(GI) &
            + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE(U_KLOC) * ( W( U_NODK_IPHA ) - NW( U_NODK_IPHA ) ) * CVNORMZ(GI)
    END DO
    !     endif

    IF( (ELE2 /= 0) .AND. (ELE2 /= ELE) ) THEN
       ! We have a discontinuity between elements so integrate along the face...
       DO U_KLOC = 1, U_NLOC
          U_KLOC2 = U_OTHER_LOC( U_KLOC )
          IF( U_KLOC2 /= 0 ) THEN
             U_NODK2_IPHA = U_NDGLN(( ELE2 - 1 ) * U_NLOC + U_KLOC2 ) + (IPHASE-1)*U_NONODS
             NDOTQNEW = NDOTQNEW &
                  + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE2(U_KLOC2) * ( U( U_NODK2_IPHA ) - NU( U_NODK2_IPHA ) ) * CVNORMX(GI) &
                  + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE2(U_KLOC2) * ( V( U_NODK2_IPHA ) - NV( U_NODK2_IPHA ) ) * CVNORMY(GI) &
                  + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE2(U_KLOC2) * ( W( U_NODK2_IPHA ) - NW( U_NODK2_IPHA ) ) * CVNORMZ(GI)
          END IF
       END DO
    END IF

    RETURN


  end SUBROUTINE GET_INT_VEL


 SUBROUTINE GET_INT_VEL_2time( NPHASE, NDOTQNEW, NDOTQ,INCOME, NDOTQOLD, INCOMEOLD, &
       HDC, GI, IPHASE,SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
       T, TOLD, FEMT, FEMTOLD, DEN, DENOLD, &
       U, V, W, NU, NV, NW,  NUOLD, NVOLD, NWOLD,&
       CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
       SUF_SIG_DIAGTEN_BC,&
       WIC_U_BC_DIRICHLET, &
       UGI_COEF_ELE, vGI_COEF_ELE, wGI_COEF_ELE, &
       UGI_COEF_ELE2, vGI_COEF_ELE2, WGI_COEF_ELE2, &
       VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
       MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
       FACE_ITS, LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
       TMIN, TMAX, TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, TOLDMAX_NOD, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       TMIN_2ND_MC, TOLDMIN_2ND_MC, TMAX_2ND_MC, TOLDMAX_2ND_MC,  LIMIT_USE_2ND,&
       is_overlapping,  &
       IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
       TUPWIND_MAT, TOLDUPWIND_MAT )
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, CV_ELE_TYPE, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, &
         NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, WIC_U_BC_DIRICHLET
    REAL, intent( in ) :: HDC
    REAL, intent(in) :: LIMT, LIMTOLD, FEMDGI, FEMTGI, FEMDOLDGI, FEMTOLDGI
    REAL, intent( inout ) :: NDOTQNEW,NDOTQ, INCOME, NDOTQOLD, INCOMEOLD
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_U_BC
    INTEGER, DIMENSION( : ), intent( in ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
         TOLDMAX_NOD
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :  ), intent( in ) :: T, TOLD, FEMT, FEMTOLD, DEN, DENOLD
    REAL, DIMENSION( :  ), intent( in ) :: U,V, W, NU, NV, NW, NUOLD, NVOLD, NWOLD
    REAL, DIMENSION( : ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC

    REAL, DIMENSION( : ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( :  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    REAL, DIMENSION( :  ), intent( in ) :: MASS_CV
    REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( :  ), intent( inout ) :: UP_WIND_NOD
    REAL, DIMENSION( : ), intent( inout ) :: TMIN, TOLDMIN, TMAX, TOLDMAX
    logical, intent( in ) :: LIMIT_USE_2ND
    REAL, DIMENSION( : ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC

    INTEGER, intent( in ) :: IANISOTROPIC
    INTEGER, intent( in ) :: NSMALL_COLM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT, TOLDUPWIND_MAT
    ! local variables
    logical, INTENT (IN) :: is_overlapping   
    INTEGER :: U_NLOC_LEV,U_KLOC_LEV,U_KLOC,U_NODK_IPHA, U_KLOC2, U_NODK2_IPHA, u_nodk2
  
    IF( is_overlapping ) THEN
       ! For overlapping basis function approach.
       CALL GET_INT_VEL_OVERLAP( NPHASE, NDOTQOLD, INCOMEOLD, &
            HDC, GI, iphase, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            TOLD, FEMTOLD, DENOLD, NUOLD,NVOLD,NWOLD, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC,SUF_V_BC, SUF_W_BC,   WIC_U_BC, &
            SUF_SIG_DIAGTEN_BC, &
            WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE,&
            UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMTOLD, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
            TOLDMIN, TOLDMAX, TMIN_NOD, TMAX_NOD,&
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            TOLDMIN_2ND_MC, TOLDMAX_2ND_MC, LIMIT_USE_2ND, &
            IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
            TOLDUPWIND_MAT )
       CALL GET_INT_VEL_OVERLAP( NPHASE, NDOTQ, INCOME, &
            HDC, GI, iphase,SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            T, FEMT, DEN, NU,NV,NW, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC,SUF_V_BC,SUF_W_BC, WIC_U_BC, &
            SUF_SIG_DIAGTEN_BC, &
            WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE,&
            UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV, OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
            TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            TMIN_2ND_MC, TMAX_2ND_MC,  LIMIT_USE_2ND, &
            IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
            TUPWIND_MAT)

    ELSE
       CALL GET_INT_VEL_ORIG( NPHASE, NDOTQOLD, INCOMEOLD, &
            HDC, GI, iphase,SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            TOLD, FEMTOLD, DENOLD, NUOLD,NVOLD,NWOLD, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC,SUF_V_BC,SUF_W_BC, WIC_U_BC, &
            WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE,&
            UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMTOLD, FEMDOLDGI, FEMTOLDGI, UP_WIND_NOD, &
            TMIN, TMAX, TMIN_NOD, TMAX_NOD,  &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            TOLDMIN_2ND_MC, TOLDMAX_2ND_MC,  LIMIT_USE_2ND )
       CALL GET_INT_VEL_ORIG( NPHASE, NDOTQ, INCOME, &
            HDC, GI, IPHASE,SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
            T, FEMT, DEN, NU,NV,NW, &
            CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC,SUF_V_BC,SUF_W_BC, WIC_U_BC, &
            WIC_U_BC_DIRICHLET, &
            UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE,&
            UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
            VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
            MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
            FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
            TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            TMIN_2ND_MC, TMAX_2ND_MC,  LIMIT_USE_2ND )
    END IF

    ! Calculate NDOTQNEW from NDOTQ
    NDOTQNEW = NDOTQ ! initialize it like this so that it contains the b.c's
    U_NLOC_LEV = U_NLOC / CV_NLOC
    DO U_KLOC = 1, U_NLOC
       U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) + (IPHASE-1)*U_NONODS
       NDOTQNEW=NDOTQNEW &
            + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE(U_KLOC) * ( U( U_NODK_IPHA ) - NU( U_NODK_IPHA ) ) * CVNORMX(GI) &
            + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE(U_KLOC) * ( V( U_NODK_IPHA ) - NV( U_NODK_IPHA ) ) * CVNORMY(GI) &
            + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE(U_KLOC) * ( W( U_NODK_IPHA ) - NW( U_NODK_IPHA ) ) * CVNORMZ(GI)
    END DO
    !     endif

    IF( (ELE2 /= 0) .AND. (ELE2 /= ELE) ) THEN
       ! We have a discontinuity between elements so integrate along the face...
       DO U_KLOC = 1, U_NLOC
          U_KLOC2 = U_OTHER_LOC( U_KLOC )
          IF( U_KLOC2 /= 0 ) THEN
             U_NODK2_IPHA = U_NDGLN(( ELE2 - 1 ) * U_NLOC + U_KLOC2 ) + (IPHASE-1)*U_NONODS
             NDOTQNEW = NDOTQNEW &
                  + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE2(U_KLOC2) * ( U( U_NODK2_IPHA ) - NU( U_NODK2_IPHA ) ) * CVNORMX(GI) &
                  + SUFEN( U_KLOC, GI ) * VGI_COEF_ELE2(U_KLOC2) * ( V( U_NODK2_IPHA ) - NV( U_NODK2_IPHA ) ) * CVNORMY(GI) &
                  + SUFEN( U_KLOC, GI ) * WGI_COEF_ELE2(U_KLOC2) * ( W( U_NODK2_IPHA ) - NW( U_NODK2_IPHA ) ) * CVNORMZ(GI)
          END IF
       END DO
    END IF

    RETURN


end SUBROUTINE GET_INT_VEL_2TIME


      PURE SUBROUTINE GET_INT_VEL_OVERLAP( NPHASE, NDOTQ,INCOME, &
       HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
       T, FEMT, DEN, NU, NV, NW,&
       CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
       SUF_SIG_DIAGTEN_BC, WIC_U_BC_DIRICHLET, &
       UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
       VOLFRA_PORE, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
       MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
       FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
       TMIN, TMAX, TMIN_NOD, TMAX_NOD, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       TMIN_2ND_MC, TMAX_2ND_MC,  LIMIT_USE_2ND, &
       IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
       TUPWIND_MAT)
    !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    ! it assumes an overlapping decomposition approach for velocity. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, CV_SNLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, WIC_U_BC_DIRICHLET, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, LIMT, FEMDGI, FEMTGI
    REAL, intent( inout ) :: NDOTQ, INCOME
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_U_BC
    INTEGER, DIMENSION( :  ), intent( in ) :: TMIN_NOD, TMAX_NOD
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :  ), intent( in ) :: T, FEMT, DEN
    REAL, DIMENSION( :  ), intent( in ) :: NU, NV, NW
    REAL, DIMENSION( : ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
    REAL, DIMENSION( : ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( :  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( : , :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    REAL, DIMENSION( :  ), intent( in ) :: MASS_CV
    REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( :  ), intent( inout ) :: UP_WIND_NOD
    REAL, DIMENSION( :  ), intent( inout ) :: TMIN, TMAX
    logical, intent( in ) :: LIMIT_USE_2ND
    REAL, DIMENSION( : ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC

      INTEGER, intent( in ) :: IANISOTROPIC
      INTEGER, intent( in ) :: NSMALL_COLM
      INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
      INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
      REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT

    ! Local variables
    REAL :: UDGI,VDGI,WDGI,  &
         UDGI2,VDGI2,WDGI2, DT_I,DT_J, &
         UDGI_INT,VDGI_INT,WDGI_INT, &
         NDOTQ_INT,NDOTQ2, &
         FEMTOLDGI_IPHA, OVER_RELAX, ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
         GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA, V_NODI, G_NODI, V_NODJ, G_NODJ, &
         GEOMTOLDGI_IPHA, W_UPWIND, W_UPWINDOLD, INCOME3, INCOME4, &
         LIMT3, FEMTGI_IPHA, GEOMTGI_IPHA, UPWIND_FRAC, TUPWIN, TUPWI2
    REAL :: NVEC(3),SUF_SIG_DIAGTEN_BC_GI(3), UGI_TMP(3)
    INTEGER :: U_KLOC,U_NODK,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_SKLOC, &
         U_SNODK,U_SNODK_IPHA, II,  COUNT, &
         U_KLOC_LEV, U_NLOC_LEV, U_SKLOC_LEV, U_SNLOC_LEV, CV_KLOC, CV_KNOD, &
         CV_KLOC2, CV_KNOD2, CV_NODI, CV_NODJ, IDIM, JDIM, IJ, MAT_NODI, MAT_NODJ, &
         CV_SKLOC, CV_SNODK, CV_SNODK_IPHA, CV_STAR_IPHA
    LOGICAL :: CONSERV, MAX_OPER, SAT_BASED, got_dt_ij
    ! IN_ELE_UPWIND=1 switches on upwinding within and element (=3 recommended).
    !    INTEGER, PARAMETER :: IN_ELE_UPWIND = 3
    !    INTEGER, PARAMETER :: IN_ELE_UPWIND = 2
    ! DG_ELE_UPWIND=1 switches on upwinding between elements (=3 recommended).
    !    INTEGER, PARAMETER :: DG_ELE_UPWIND = 3
    !    INTEGER, PARAMETER :: DG_ELE_UPWIND = 2
    ! FORCE_UPWIND_VEL forces the use of upwinding for velocity - good for testing as high order can be complex. 
    ! FORCE_UPWIND_VEL_DG_ELE forces upwind vel between elements
    LOGICAL, PARAMETER :: FORCE_UPWIND_VEL = .false.
    !      LOGICAL, PARAMETER :: FORCE_UPWIND_VEL = .true.
    LOGICAL, PARAMETER :: FORCE_UPWIND_VEL_DG_ELE = .true.
    LOGICAL, PARAMETER :: LIM_VOL_ADJUST2 = .true.
    ! If ROW_AVE then use a Roe averaged flux between the elements or CV's when ROE_AVE =.true. 
    LOGICAL, PARAMETER :: ROE_AVE = .false. !.false. !.true.
    ! limit the volume fraction based on net letting flux go into certain elements for DG between the elements...
    LOGICAL, PARAMETER :: LIMIT_SAT_BASED_INTERP = .false.
    ! option for between element velocity calculation...
    integer, PARAMETER :: between_ele_dg_opt = 5
    ! limit the volume fraction based on switching to the 1st order scheme for DG between the elements...
    ! which should be ued with LIMIT_WITHIN_REAL=.true.
    ! LIMIT_WITHIN_REAL makes sure the saturations are bounded can be used with DG 
    LOGICAL, PARAMETER :: LIMIT_SAT_BASED_UPWIND = .true.
    LOGICAL, PARAMETER :: LIMIT_WITHIN_REAL=.false.
    ! The fracture modelling with the upwind based method for DG... 
    LOGICAL, PARAMETER :: high_order_upwind_vel_for_dg = .true.
    LOGICAL :: RESET_STORE, LIM_VOL_ADJUST, enforce_abs
    REAL :: TMIN_STORE, TMAX_STORE
    REAL :: PERM_TILDE,NDOTQ_TILDE, NDOTQ2_TILDE, NDOTQOLD_TILDE, NDOTQOLD2_TILDE, rden_ave, Q_UNDERLY
    REAL :: NDOTQ_KEEP_IN,  NDOTQ_KEEP, NDOTQ2_KEEP
    ! coefficients for this element ELE
    real :: gamma, abs_tilde, abs_tilde_i, abs_tilde_j, abs_max, abs_min, w_relax, grad2nd
    real :: abs_tilde_2nd, abs_tildeold_2nd, max_nodtq_keep, min_nodtq_keep
    real :: w_weight, relax

    real :: DT_I_upwind, DT_J_upwind
    real :: dt_max, dt_min
    real :: abs_tilde1, abs_tilde2
    real :: wrelax, wrelax1, wrelax2
    real :: T_PELE, T_PELEOT, TMIN_PELE, TMAX_PELE, TMIN_PELEOT, TMAX_PELEOT


    ! print *,'IN_ELE_UPWIND,CV_DG_VEL_INT_OPT,  DG_ELE_UPWIND, IANISOTROPIC:',IN_ELE_UPWIND,CV_DG_VEL_INT_OPT,  DG_ELE_UPWIND, IANISOTROPIC
    !  stop 8721

    ! The adjustment method for variable CV volumes is not ready for new limiter method...
    LIM_VOL_ADJUST = ( LIM_VOL_ADJUST2 .AND. (.NOT.LIMIT_USE_2ND) ) .AND. IANISOTROPIC==0

    UGI_COEF_ELE=0.0
    VGI_COEF_ELE=0.0
    WGI_COEF_ELE=0.0 

    ! coefficients for this element ELE2
    UGI_COEF_ELE2=0.0
    VGI_COEF_ELE2=0.0
    WGI_COEF_ELE2=0.0 

    U_NLOC_LEV = U_NLOC / CV_NLOC
    U_SNLOC_LEV = U_SNLOC / CV_NLOC


    got_dt_ij=.false.


    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain. 
       IF( WIC_U_BC( SELE + ( IPHASE - 1 ) * STOTEL) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_NODK_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) + (IPHASE-1)*U_NONODS
             UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
             VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
             WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )
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
                SUF_SIG_DIAGTEN_BC_GI( 1:NDIM ) = SUF_SIG_DIAGTEN_BC( CV_SNODK_IPHA, 1:NDIM )
             ENDIF
          END DO

          ! Only modify boundary velocity for incoming velocity...
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

          IF(UDGI*CVNORMX(GI)+VDGI*CVNORMY(GI)+WDGI*CVNORMZ(GI).LT.0.0) THEN ! Incomming...
             UGI_TMP = SUF_SIG_DIAGTEN_BC_GI(1:3) * (/UDGI, VDGI, WDGI/)
             UDGI=UGI_TMP(1) ; VDGI=UGI_TMP(2) ; WDGI=UGI_TMP(3)
          ENDIF

       ELSE ! Specified vel bc.
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
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
             ELSE
                UDGI = UDGI + SUFEN( U_KLOC, GI ) * SUF_U_BC( U_SNODK_IPHA )
                VDGI = VDGI + SUFEN( U_KLOC, GI ) * SUF_V_BC( U_SNODK_IPHA )
                WDGI = WDGI + SUFEN( U_KLOC, GI ) * SUF_W_BC( U_SNODK_IPHA )
             END IF
          END DO
       END IF

    ELSE ! Conditional_SELE. Not on the boundary of the domain.
       Conditional_ELE2: IF(( ELE2 == 0 ).OR.( ELE2 == ELE)) THEN
          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
          UDGI2 = 0.0
          VDGI2 = 0.0
          WDGI2 = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC =(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_KLOC2=(CV_JLOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_NODK_IPHA  = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC ) +(IPHASE-1)*U_NONODS
             U_NODK2_IPHA = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC2 )+(IPHASE-1)*U_NONODS

             UDGI = UDGI + SUFEN( U_KLOC, GI ) * NU( U_NODK_IPHA )
             VDGI = VDGI + SUFEN( U_KLOC, GI ) * NV( U_NODK_IPHA )
             WDGI = WDGI + SUFEN( U_KLOC, GI ) * NW( U_NODK_IPHA )

             UDGI2 = UDGI2 + SUFEN( U_KLOC2, GI ) * NU( U_NODK2_IPHA )
             VDGI2 = VDGI2 + SUFEN( U_KLOC2, GI ) * NV( U_NODK2_IPHA )
             WDGI2 = WDGI2 + SUFEN( U_KLOC2, GI ) * NW( U_NODK2_IPHA )

          END DO

          NDOTQ = CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI  &
               + CVNORMZ(GI) * WDGI 
          NDOTQ2 =  CVNORMX( GI ) * UDGI2 + CVNORMY( GI ) * VDGI2  &
               + CVNORMZ(GI) * WDGI2 

          NDOTQ_TILDE = 0.5 * ( NDOTQ +  NDOTQ2 )
          NDOTQ2_TILDE =   NDOTQ_TILDE

          CV_NODI = CV_NODI_IPHA - (IPHASE-1)*CV_NONODS
          CV_NODJ = CV_NODJ_IPHA - (IPHASE-1)*CV_NONODS

          IF((IN_ELE_UPWIND==1).OR.FORCE_UPWIND_VEL) THEN

             IF(NDOTQ_TILDE < 0.0) THEN
                INCOME=1.0
             ELSE
                INCOME=0.0
             END IF

             IF(abs(NDOTQ2 - NDOTQ).lt. 1.e-7) INCOME=0.5



          ELSE IF(IN_ELE_UPWIND==2) THEN ! the best

             UPWIND_FRAC=0.8
             IF(NDOTQ_TILDE < 0.0) THEN
                INCOME=0.8
                !INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ELSE
                INCOME=0.8
                !INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             END IF

          ELSE IF(IN_ELE_UPWIND==3) THEN ! the best optimal upwind frac.


             MAT_NODI=MAT_NDGLN((ELE-1)*MAT_NLOC+CV_ILOC)
             MAT_NODJ=MAT_NDGLN((ELE-1)*MAT_NLOC+CV_JLOC)

             FEMTGI_IPHA = 0.0
             DO CV_KLOC=1,CV_NLOC
                CV_KNOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                FEMTGI_IPHA = FEMTGI_IPHA &
                     + SCVFEN(CV_KLOC,GI) * FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
             END DO
             ! Central is fine as its within an element with equally spaced nodes. 
             !FEMTGI_IPHA = 0.5*( t(cv_nodi_ipha)+t(cv_nodj_ipha) )
             !FEMTOLDGI_IPHA = 0.5*( told(cv_nodi_ipha)+told(cv_nodj_ipha) )
             if(IANISOTROPIC==0) then ! this is the only needed for isotropic limiting for velocity...
                FEMTGI_IPHA = ( MASS_CV(CV_NODJ) * T(CV_NODI_IPHA) + &
                     MASS_CV(CV_NODI) * T(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             endif
             GEOMTGI_IPHA = ( MASS_CV(CV_NODJ) * T(CV_NODI_IPHA) + &
                  MASS_CV(CV_NODI) * T(CV_NODJ_IPHA) ) / (MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))

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
                   G_NODI = G_NODI + OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM) * NVEC(JDIM)
                   IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODJ-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                   V_NODJ = V_NODJ + OPT_VEL_UPWIND_COEFS(IJ) * NVEC(JDIM)
                   G_NODJ = G_NODJ + OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM) * NVEC(JDIM)
                END DO
                ABS_CV_NODI_IPHA      = ABS_CV_NODI_IPHA + NVEC(IDIM)*V_NODI
                GRAD_ABS_CV_NODI_IPHA = GRAD_ABS_CV_NODI_IPHA + NVEC(IDIM)*G_NODI
                ABS_CV_NODJ_IPHA      = ABS_CV_NODJ_IPHA + NVEC(IDIM)*V_NODJ
                GRAD_ABS_CV_NODJ_IPHA = GRAD_ABS_CV_NODJ_IPHA + NVEC(IDIM)*G_NODJ
             END DO

             !        print *,'-iphase,cv_nodi,cv_nodj,ABS_CV_NODI_IPHA,ABS_CV_NODJ_IPHA:',iphase,cv_nodi,cv_nodj,ABS_CV_NODI_IPHA,ABS_CV_NODJ_IPHA





             ! Make sure we have some sort of velocity (only needed between elements)...
             ! take the mean of the underlying velocity...
             !                        Q_UNDERLY=0.5*( NDOTQ*ABS_CV_NODI_IPHA + NDOTQ2*ABS_CV_NODJ_IPHA )
             !                        NDOTQ_TILDE =Q_UNDERLY/ABS_CV_NODI_IPHA
             !                        NDOTQ2_TILDE=Q_UNDERLY/ABS_CV_NODJ_IPHA
             !
             !                        QOLD_UNDERLY=0.5*( NDOTQOLD*ABS_CV_NODI_IPHA + NDOTQOLD2*ABS_CV_NODJ_IPHA )
             !                        NDOTQOLD_TILDE =QOLD_UNDERLY/ABS_CV_NODI_IPHA
             !                        NDOTQOLD2_TILDE=QOLD_UNDERLY/ABS_CV_NODJ_IPHA

             ! These are the new limits of the velocities...
             NDOTQ_KEEP  = NDOTQ   ! this is associated with saturation at NODI
             NDOTQ2_KEEP = NDOTQ2  ! this is associated with saturation at NODJ

             ! between these limits work out which are associated with the flux limited saturation at this interface. 

             if(ROE_AVE) then

                ! do the Roe average of the rest of the velocity...
                NDOTQ_TILDE  = ( DEN(CV_NODI_IPHA) * T(CV_NODI_IPHA) * NDOTQ -  &
                     &           DEN(CV_NODJ_IPHA) * T(CV_NODJ_IPHA) * NDOTQ2 ) & 
                     / tolfun( VOLFRA_PORE(ELE)*DEN(CV_NODI_IPHA) * T(CV_NODI_IPHA) - VOLFRA_PORE(ELE)*DEN(CV_NODJ_IPHA) * T(CV_NODJ_IPHA) )
                NDOTQ2_TILDE = NDOTQ_TILDE

                ! Make sure we have some sort of velocity (only needed between elements)...

             endif



             ! high order order (low order is an option above)



             IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) THEN
                INCOME3=1.0
                ! Make sure the fem description is not biased to the downwind...         
             ELSE
                INCOME3=0.0
             END IF

             IF ( LIM_VOL_ADJUST ) THEN
                RESET_STORE=.FALSE. 
                CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
                     CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
                CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
                     CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
             END IF



             if(IANISOTROPIC==0) then ! this is the only good isotropic limiting for velocity...

                T_PELE = T( CV_NODI_IPHA )
                IF( CV_NODJ /= CV_NODI ) THEN
                   T_PELEOT    = T( CV_NODJ_IPHA )
                   TMIN_PELE   = TMIN( CV_NODI_IPHA )
                   TMAX_PELE   = TMAX( CV_NODI_IPHA )
                   TMIN_PELEOT = TMIN( CV_NODJ_IPHA )
                   TMAX_PELEOT = TMAX( CV_NODJ_IPHA )
                END IF
                CALL ONVDLIMsqrt( CV_NONODS, &
                     LIMT3, FEMTGI_IPHA, INCOME3, CV_NODI, CV_NODJ, &
                     T_PELE, T_PELEOT, TMIN_PELE, TMAX_PELE, TMIN_PELEOT, TMAX_PELEOT, .FALSE. , .FALSE. )


             else

                CV_STAR_IPHA = 1 + ( IPHASE - 1 ) * CV_NONODS
                CV_NODI = CV_NODI_IPHA-(IPHASE-1)*CV_NONODS
                CV_NODJ = CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS

                IF ( IANISOTROPIC == 1 ) THEN
                   IF( CV_NODI /= CV_NODJ ) THEN
                      DO COUNT = SMALL_FINDRM(CV_NODJ), SMALL_FINDRM(CV_NODJ+1)-1
                         IF ( SMALL_COLM(COUNT) == CV_NODI ) THEN
                            TUPWIN = TUPWIND_MAT( COUNT )
                            EXIT
                         END IF
                      END DO
                      DO COUNT = SMALL_FINDRM(CV_NODI), SMALL_FINDRM(CV_NODI+1)-1
                         IF ( SMALL_COLM(COUNT) == CV_NODJ ) THEN
                            TUPWI2 = TUPWIND_MAT( COUNT )
                            EXIT
                         END IF
                      END DO
                   END IF
                END IF

                CALL ONVDLIM_ALL( CV_NONODS, &
                     LIMT3, FEMTGI_IPHA, INCOME3, CV_NODI, CV_NODJ, &
                     T( CV_STAR_IPHA:CV_STAR_IPHA +CV_NONODS-1 ), TMIN( CV_STAR_IPHA:CV_STAR_IPHA +CV_NONODS-1 ), TMAX( CV_STAR_IPHA:CV_STAR_IPHA +CV_NONODS-1 ), &
                     TMIN_2ND_MC( CV_STAR_IPHA:CV_STAR_IPHA +CV_NONODS-1 ), TMAX_2ND_MC( CV_STAR_IPHA:CV_STAR_IPHA +CV_NONODS-1 ), .FALSE., .FALSE., LIMIT_USE_2ND, -1.0, &
                     IANISOTROPIC, TUPWIN, TUPWI2 )               

             endif




             IF ( LIM_VOL_ADJUST ) THEN
                RESET_STORE=.TRUE. 
                CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
                     CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
                CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
                     CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME3 )
             END IF




             abs_tilde_i = ABS_CV_NODI_IPHA
             abs_tilde_j = ABS_CV_NODJ_IPHA




             ! amend the absorption between abs_tilde and the upwind value...
             ! ***WITHIN AN ELEMENT****



             abs_tilde =  0.5*(  ABS_CV_NODI_IPHA  + ( limt3   -  T(CV_NODI_IPHA)  ) * GRAD_ABS_CV_NODI_IPHA   +   &
                  ABS_CV_NODJ_IPHA  + ( limt3   -  T(CV_NODJ_IPHA)  ) * GRAD_ABS_CV_NODJ_IPHA )


             abs_max=max(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)
             abs_min=min(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)

             abs_tilde    = min(abs_max, max(abs_min,  abs_tilde ))

             NDOTQ_KEEP_IN    =  0.5*( NDOTQ*ABS_CV_NODI_IPHA    + NDOTQ2*ABS_CV_NODJ_IPHA )    /abs_tilde

             INCOME    = MIN(1.0, MAX(0.0,  (NDOTQ_KEEP_IN - NDOTQ)/TOLFUN( NDOTQ2 - NDOTQ ) ))
             IF(abs(NDOTQ2 - NDOTQ).lt. 1.e-7) INCOME=0.5

             !                     INCOME =0.5*ABS_CV_NODI_IPHA* MASS_CV(CV_NODI) /(0.5*(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI) +ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ) ))
             !                     INCOMEOLD =0.5*ABS_CV_NODI_IPHA* MASS_CV(CV_NODI) /(0.5*(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI) +ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ) ))

             !               stop 821




             ! based on switching to the 1st order scheme...
             if(LIMIT_SAT_BASED_UPWIND) then
                if( min(ABS_CV_NODI_IPHA,ABS_CV_NODJ_IPHA)/ (ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA).lt.1.e-3) then ! resort to using the velocity from the highest absorption cell. 
                   !                 if( ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA.gt.1.0e+3 ) then ! resort to using the velocity from the highest absorption cell. 
                   !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) THEN 
                   IF(0.5*(NDOTQ+NDOTQ2) < 0.0) THEN 
                      INCOME=1.0
                   ELSE
                      INCOME=0.0
                   END IF
                endif
             endif





          ELSE

             UPWIND_FRAC=0.5
             IF(NDOTQ_TILDE < 0.0) THEN
                !INCOME=0.8
                INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODJ)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             ELSE
                !INCOME=0.2
                INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV(CV_NODI)/(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
             END IF
             !INCOME=0.5
          END IF

          UDGI = 0.0
          VDGI = 0.0
          WDGI = 0.0
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
          END DO

          UDGI2 = 0.0
          VDGI2 = 0.0
          WDGI2 = 0.0
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

             END IF
          END DO

          NDOTQ = CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI  &
               + CVNORMZ(GI) * WDGI

          NDOTQ2 =  CVNORMX( GI ) * UDGI2 + CVNORMY( GI ) * VDGI2  &
               + CVNORMZ(GI) * WDGI2 

          CV_NODI=CV_NODI_IPHA -(IPHASE-1)*CV_NONODS
          CV_NODJ=CV_NODJ_IPHA -(IPHASE-1)*CV_NONODS

          IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
             DT_I=1.0
             DT_J=1.0
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
             DT_I=MAX(1.E-2,T( CV_NODI_IPHA ))
             DT_J=MAX(1.E-2,T( CV_NODJ_IPHA ))
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
             DT_I=DEN( CV_NODI_IPHA )*T( CV_NODI_IPHA )
             DT_J=DEN( CV_NODJ_IPHA )*T( CV_NODJ_IPHA )
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 4) THEN


             ! DG_ELE_UPWIND==1: the upwind method. 
             ! DG_ELE_UPWIND==3: the best optimal upwind frac.


             MAT_NODI=MAT_NDGLN((ELE-1)*MAT_NLOC+CV_ILOC)
             MAT_NODJ=MAT_NDGLN((ELE2-1)*MAT_NLOC+CV_JLOC)


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


             ! Make sure we have some sort of velocity (only needed between elements)...
             ! take the mean of the underlying velocity...
             if(.false.) then
                Q_UNDERLY=( NDOTQ*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ NDOTQ2*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ) ) &
                     /(MASS_CV(CV_NODI)+MASS_CV(CV_NODJ))
                NDOTQ_TILDE =Q_UNDERLY/ABS_CV_NODI_IPHA
                NDOTQ2_TILDE=Q_UNDERLY/ABS_CV_NODJ_IPHA

             else ! better tested option...

                Q_UNDERLY=0.5*( NDOTQ*ABS_CV_NODI_IPHA + NDOTQ2*ABS_CV_NODJ_IPHA )
                NDOTQ_TILDE =Q_UNDERLY/ABS_CV_NODI_IPHA
                NDOTQ2_TILDE=Q_UNDERLY/ABS_CV_NODJ_IPHA

             endif

             !   Q_UNDERLY=0.5*( NDOTQ/ABS_CV_NODI_IPHA + NDOTQ2/ABS_CV_NODJ_IPHA )
             !   NDOTQ_TILDE =Q_UNDERLY*ABS_CV_NODI_IPHA
             !   NDOTQ2_TILDE=Q_UNDERLY*ABS_CV_NODJ_IPHA

             !   QOLD_UNDERLY=0.5*( NDOTQOLD/ABS_CV_NODI_IPHA + NDOTQOLD2/ABS_CV_NODJ_IPHA )
             !   NDOTQOLD_TILDE =QOLD_UNDERLY*ABS_CV_NODI_IPHA
             !   NDOTQOLD2_TILDE=QOLD_UNDERLY*ABS_CV_NODJ_IPHA

             ! These are the new limits of the velocities...
             NDOTQ_KEEP  = NDOTQ_TILDE   ! this is associated with saturation at NODI
             NDOTQ2_KEEP = NDOTQ2_TILDE  ! this is associated with saturation at NODJ

             ! between these limits work out which are associated with the flux limited saturation at this interface. 

             if(ROE_AVE) then

                ! do the Roe average of the rest of the velocity...
                NDOTQ_TILDE  = ( DEN(CV_NODI_IPHA) * T(CV_NODI_IPHA) * NDOTQ -  &
                     &           DEN(CV_NODJ_IPHA) * T(CV_NODJ_IPHA) * NDOTQ2 ) & 
                     / tolfun( VOLFRA_PORE(ELE)*DEN(CV_NODI_IPHA) * T(CV_NODI_IPHA) - VOLFRA_PORE(ELE2)*DEN(CV_NODJ_IPHA) * T(CV_NODJ_IPHA) )
                NDOTQ2_TILDE = NDOTQ_TILDE

                ! Make sure we have some sort of velocity (only needed between elements)...

             endif


             ! high order dg...
             IF ( high_order_upwind_vel_for_dg ) THEN ! high order


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
                      end if

                      if (.true.) then
                         if (0.5*(NDOTQ_TILDE+NDOTQ2_TILDE)>0.0) then
                            CV_KNOD = CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 * FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
                         else
                            CV_KNOD2 = CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 * FEMT(CV_KNOD2+(IPHASE-1)*CV_NONODS)
                         end if
                      end if

                      if (.false.) then
                         CV_KNOD = CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
                         FEMTGI_IPHA = FEMTGI_IPHA &
                              + SCVFEN(CV_KLOC,GI) * FEMT(CV_KNOD+(IPHASE-1)*CV_NONODS)
                      end if
                   END IF
                END DO




                IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) THEN
                   INCOME3=1.0
                ELSE
                   INCOME3=0.0
                END IF

                ! Between elements...




                ! redefine so that it detects oscillations...
                abs_tilde1 =  ABS_CV_NODI_IPHA  + 0.5*( T(CV_NODJ_IPHA)   -  T(CV_NODI_IPHA)  ) * GRAD_ABS_CV_NODI_IPHA   
                abs_tilde2 =  ABS_CV_NODJ_IPHA  + 0.5*( T(CV_NODi_IPHA)   -  T(CV_NODJ_IPHA)  ) * GRAD_ABS_CV_NODJ_IPHA 


                abs_max=max(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)
                abs_min=min(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)

                abs_tilde1    = min(1000.*abs_max, max(0.001*abs_min,  abs_tilde1 )) 

                abs_tilde2    = min(1000.*abs_max, max(0.001*abs_min,  abs_tilde2 )) 



                got_dt_ij=.true.
                ! new high order method - works well...

                wrelax= min(  &
                     ABS_CV_NODI_IPHA/abs_tilde1, abs_tilde1/ABS_CV_NODI_IPHA,   &
                     ABS_CV_NODJ_IPHA/abs_tilde2, abs_tilde2/ABS_CV_NODJ_IPHA    )

                !     wrelax= min(  &
                !       ABS_CV_NODI_IPHA/ABS_CV_NODJ_IPHA, ABS_CV_NODJ_IPHA/ABS_CV_NODI_IPHA   )

                if(.true.) then
                   wrelax1= min( 1.0, abs_tilde1/ABS_CV_NODI_IPHA )
                   wrelax2= min( 1.0, abs_tilde2/ABS_CV_NODJ_IPHA  )
                else
                   wrelax1= min( ABS_CV_NODI_IPHA/abs_tilde1, abs_tilde1/ABS_CV_NODI_IPHA )
                   wrelax2= min( ABS_CV_NODJ_IPHA/abs_tilde2, abs_tilde2/ABS_CV_NODJ_IPHA  )

                endif


                if(income3.lt.0.5) then ! flux limit
                   wrelax2=wrelax1
                else
                   wrelax1=wrelax2
                endif


                ! no 4 (more stable than 2 and 3):
                if(between_ele_dg_opt==1) then
                   if ( income3 < 0.5 ) then ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*0.5
                   else
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*0.5
                   endif

                else if(between_ele_dg_opt==2) then
                   ! no 4 (more stable than 2 and 3) (also make sure we dont violate bounds):
                   if ( income3 < 0.5 ) then ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*0.5
                      !              IF(T(CV_NODI_IPHA).LT.0.2) THEN
                      IF(ABS_CV_NODI_IPHA.GT.1.0E+10) THEN
                         DT_I = 0.0
                         DT_J = 0.0
                      ENDIF
                   else
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*0.5
                      !               IF(T(CV_NODJ_IPHA).LT.0.2) THEN
                      IF(ABS_CV_NODJ_IPHA.GT.1.0E+10) THEN
                         DT_I = 0.0
                         DT_J = 0.0
                      ENDIF
                   endif

                   ! no 2(again):
                else if(between_ele_dg_opt==3) then
                   if ( income3 < 0.5 ) then ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   else
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   endif


                   ! no 2(again):
                else if(between_ele_dg_opt==4) then
                   if ( income3 < 0.5 ) then ! upwind
                      DT_I = ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      IF(ABS_CV_NODI_IPHA.GT.1.0E+10) THEN
                         DT_I = 0.0
                         DT_J = 0.0
                      ENDIF
                   else
                      DT_I = ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      IF(ABS_CV_NODJ_IPHA.GT.1.0E+10) THEN
                         DT_I = 0.0
                         DT_J = 0.0
                      ENDIF
                   endif

                   ! no 2(again):
                else if(between_ele_dg_opt==5) then
                   if ( income3 < 0.5 ) then ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      IF(ABS_CV_NODI_IPHA.GT.1.0E+10) THEN
                         DT_I = 0.0
                         DT_J = 0.0
                      ENDIF
                   else
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      IF(ABS_CV_NODJ_IPHA.GT.1.0E+10) THEN
                         DT_I = 0.0
                         DT_J = 0.0
                      ENDIF
                   endif

                   ! no 3(again): XXX failed gravity problem...
                else if(between_ele_dg_opt==6) then
                   if ( income3 < 0.5 ) then ! upwind
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde1)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   else
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde2)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde2)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   endif

                   ! no 2(again,again) XXXXfor the gravity problem:
                else if(between_ele_dg_opt==7) then
                   if ( income3 < 0.5 ) then ! upwind
                      DT_I = (1.-wrelax1)*(1.-min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1))  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   else
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde2)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                      DT_J = (1.-wrelax2)*(1.-min(1.0,ABS_CV_NODI_IPHA/abs_tilde2))  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                           /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   endif

                   ! no 3 soln...
                else if(between_ele_dg_opt==8) then
                   DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA  ) &
                        /(ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA)   + wrelax1*0.5
                   DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA  ) &
                        /(ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA)   + wrelax2*0.5

                   ! no 2 soln...
                   !        DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  ) &
                   !             /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))   + wrelax1*0.5
                   !        DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  ) &
                   !             /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))   + wrelax2*0.5

                   !        DTOLD_I = ( (1.-wrelaxold1)*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  ) &
                   !                /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))  + wrelaxold1*0.5
                   !        DTOLD_J = ( (1.-wrelaxold2)*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  ) &
                   !                /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))  + wrelaxold2*0.5

                else if(between_ele_dg_opt==9) then
                   ! no 0 soln (more stable than 2 & 3)...
                   DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI) + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ) ) &
                        /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ) + wrelax2*ABS_CV_NODI_IPHA*MASS_CV(CV_NODI) ) &
                        /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   ! no 0 soln (more stable than 2 & 3)...
                else if(between_ele_dg_opt==10) then
                   DT_I = ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ)  &
                        /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                   DT_J = ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)  &
                        /(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI)+ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ))
                endif





                ! END OF IF ( DG_ELE_UPWIND==3 ) THEN ...
             ELSE  ! Low order...
                !if ( abs(T(CV_NODI_IPHA)-T(CV_NODJ_IPHA)).lt.1.e-4 ) then ! low order


                INCOME =0.5*ABS_CV_NODI_IPHA* MASS_CV(CV_NODI) /(0.5*(ABS_CV_NODI_IPHA*MASS_CV(CV_NODI) +ABS_CV_NODJ_IPHA*MASS_CV(CV_NODJ) ))


                ! limit the volume fraction based on net letting flux go into certain elements...
                !     if(.false.) then
                if(LIMIT_SAT_BASED_INTERP) then
                   enforce_abs = .false.

                   if(iphase==1) then
                      !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) >= 0.0) THEN 
                      IF(0.5*(NDOTQ+NDOTQ2) >= 0.0) THEN 
                         IF(T(CV_NODi_IPHA) < 0.5) then
                            enforce_abs = .true.
                            relax=min( 10.*(0.5-T(CV_NODi_IPHA)), 1.0)
                         ENDIF
                      ELSE
                         IF(T(CV_NODj_IPHA) < 0.5) THEN
                            enforce_abs = .true. 
                            relax=min( 10.*(0.5-T(CV_NODJ_IPHA)), 1.0)
                         ENDIF
                      ENDIF
                   endif

                   if(iphase==2) then
                      !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) >= 0.0) THEN 
                      IF(0.5*(NDOTQ+NDOTQ2) >= 0.0) THEN 
                         IF(T(CV_NODi_IPHA) < 0.5) then
                            enforce_abs = .true.
                            relax=min( 10.*(0.5-T(CV_NODi_IPHA)), 1.0)
                         ENDIF
                      ELSE
                         IF(T(CV_NODj_IPHA) < 0.5) THEN
                            enforce_abs = .true. 
                            relax=min( 10.*(0.5-T(CV_NODJ_IPHA)), 1.0)
                         ENDIF
                      ENDIF
                   endif



                   if(enforce_abs) then
                      INCOME =  (1.-RELAX)*INCOME  + RELAX*( (1.-INCOME)**2 *(1. - 3.*INCOME) +  INCOME*(1.-INCOME)*(1.+3.*INCOME)  ) ! non-lumped
                   endif

                   income=min(1.0, max(0.0, income))
                endif





                ! endof IF ( DG_ELE_UPWIND==3 ) THEN ELSE...
             ENDIF








             if(.not.got_dt_ij) then

                DT_I = 1.0 - INCOME
                DT_J = 1.0 - DT_I

             endif


          END IF



          ! Amend weighting for porosity only across elements...
          IF((ABS(CV_DG_VEL_INT_OPT ) == 2).OR.(ABS(CV_DG_VEL_INT_OPT ) == 3)) THEN 
             IF(ELE /= ELE2) THEN 
                DT_I=VOLFRA_PORE(ELE) *DT_I
                DT_J=VOLFRA_PORE(ELE2)*DT_J
             END IF
          END IF

          !dt_i=0.5
          !dt_j=0.5
          !dtold_i=0.5
          !dtold_j=0.5

          if(.not.got_dt_ij) then ! normalize...

             DT_I_upwind = DT_I
             DT_J_upwind = DT_J

             DT_I =  DT_I_upwind/(DT_I_upwind+DT_J_upwind)
             DT_J =  DT_J_upwind/(DT_I_upwind+DT_J_upwind)

          else
             income=income3
          endif





          UDGI_INT = ( DT_I * UDGI + DT_J * UDGI2 ) 
          VDGI_INT = ( DT_I * VDGI + DT_J * VDGI2 ) 
          WDGI_INT = ( DT_I * WDGI + DT_J * WDGI2 ) 



          IF( CV_DG_VEL_INT_OPT < 0 ) THEN
   !          FLAbort('2821') ! should not be going in here??
   !          if(got_dt_ij) FLAbort('2611') ! we can not use these CV_DG_VEL_INT_OPT < 0, got_dt_ij=.true. options together.

             NDOTQ_INT = CVNORMX( GI ) * UDGI_INT + CVNORMY( GI ) * VDGI_INT + &
                  CVNORMZ( GI ) * WDGI_INT

             IF( NDOTQ_INT <= 0.0 ) THEN ! Incoming
                !DT_I=1.0
                DT_J=DT_I+DT_J
             ELSE
                DT_I=DT_I+DT_J
                !DT_J=1.0
             END IF

             UDGI_INT = ( DT_I * UDGI + DT_J * UDGI2 ) / (DT_I + DT_J)
             VDGI_INT = ( DT_I * VDGI + DT_J * VDGI2 ) / (DT_I + DT_J)
             WDGI_INT = ( DT_I * WDGI + DT_J * WDGI2 ) / (DT_I + DT_J)

          END IF

          UDGI  = UDGI_INT 
          VDGI  = VDGI_INT 
          WDGI  = WDGI_INT

          UGI_COEF_ELE(:)=DT_I * UGI_COEF_ELE(:) 
          VGI_COEF_ELE(:)=DT_I * VGI_COEF_ELE(:) 
          WGI_COEF_ELE(:)=DT_I * WGI_COEF_ELE(:) 

          UGI_COEF_ELE2(:)=DT_J * UGI_COEF_ELE2(:) 
          VGI_COEF_ELE2(:)=DT_J * VGI_COEF_ELE2(:) 
          WGI_COEF_ELE2(:)=DT_J * WGI_COEF_ELE2(:) 


       END IF Conditional_ELE2

    END IF Conditional_SELE

    NDOTQ =  CVNORMX( GI ) * UDGI + CVNORMY( GI ) * VDGI + CVNORMZ(GI) * WDGI

    ! Define whether flux is incoming or outgoing, depending on direction of flow
    INCOME = 1.
    IF( NDOTQ >= 0. ) INCOME = 0.
    !  IF( NDOTQ_tilde >= 0. ) INCOME = 0.
    !IF( NDOTQ > 0. ) INCOME = 0.

    RETURN  

  END SUBROUTINE GET_INT_VEL_OVERLAP




    PURE SUBROUTINE GET_INT_VEL_ORIG( NPHASE, NDOTQ, INCOME, &
       HDC, GI, IPHASE, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, U_NDGLN, &
       T, FEMT,DEN, NU, NV, NW, &
       CV_NODI_IPHA, CV_NODJ_IPHA, CVNORMX, CVNORMY, CVNORMZ,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC, SUF_V_BC, SUF_W_BC, WIC_U_BC, &
       WIC_U_BC_DIRICHLET, &
       UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2, &
       VOLFRA_PORE, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_NDGLN, CV_OTHER_LOC, &
       MASS_CV,OPT_VEL_UPWIND_COEFS,NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NDGLN, MAT_NONODS, &
       FACE_ITS, LIMT, FEMDGI, FEMTGI, UP_WIND_NOD, &
       TMIN, TMAX, TMIN_NOD, TMAX_NOD,&
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       TMIN_2ND_MC, TMAX_2ND_MC, LIMIT_USE_2ND )

    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, IPHASE, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ_IPHA, CV_NODI_IPHA, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, WIC_U_BC_DIRICHLET, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NOPT_VEL_UPWIND_COEFS, NDIM, MAT_NLOC, MAT_NONODS, FACE_ITS, &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, LIMT, FEMDGI, FEMTGI
    REAL, intent( inout ) :: NDOTQ, INCOME
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_U_BC
    INTEGER, DIMENSION( :  ), intent( in ) :: TMIN_NOD, TMAX_NOD
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :  ), intent( in ) :: T, FEMT, DEN
    REAL, DIMENSION( :  ), intent( in ) :: NU, NV, NW
    REAL, DIMENSION( : ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( : ), intent( in ) :: SUF_U_BC, SUF_V_BC, SUF_W_BC
    REAL, DIMENSION( : ), intent( inout ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( :  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    REAL, DIMENSION( :  ), intent( in ) :: MASS_CV
    REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
    INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( :  ), intent( inout ) :: UP_WIND_NOD
    REAL, DIMENSION( :  ), intent( inout ) :: TMIN, TMAX
    logical, intent( in ) :: LIMIT_USE_2ND
    REAL, DIMENSION( : ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC
    integer :: s,e

    ! Local variables
    REAL :: UDGI,VDGI,WDGI,  &
         UDGI2,VDGI2,WDGI2, DT_I,DT_J, &
         UDGI_INT,VDGI_INT,WDGI_INT, &
         NDOTQ_INT
    INTEGER :: U_KLOC,U_NODK,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_NODK2,U_SKLOC, &
         U_SNODK,U_SNODK_IPHA, II, CV_NODI, CV_NODJ, IDIM
    integer, dimension(U_SNLOC) ::  U_NODK_IPHA_V, U_SNODK_IPHA_V
    REAL, DIMENSION(NDIM) :: CVNORMX_ALL

    ! Local variable for indirect addressing
    REAL, DIMENSION ( :, :, : ), allocatable :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
    REAL, DIMENSION ( :, : ), allocatable :: UDGI_ALL, UDGI2_ALL, UDGI_INT_ALL
    REAL, DIMENSION ( :, :, : ), allocatable :: LOC_NU, SLOC_NU, SLOC2_NU, SUF_U_BC_ALL
    REAL, DIMENSION ( :, : ), allocatable :: GLOB_T, GLOB_DEN


    ALLOCATE( UGI_COEF_ELE_ALL( NDIM, NPHASE, U_NLOC) )
    ALLOCATE( UGI_COEF_ELE2_ALL( NDIM, NPHASE, U_NLOC) )
    ALLOCATE( UDGI_ALL( NDIM, NPHASE ), UDGI2_ALL( NDIM, NPHASE ), UDGI_INT_ALL( NDIM, NPHASE ) )
    ALLOCATE( LOC_NU(NDIM, NPHASE, U_NLOC) )
    ALLOCATE( SLOC_NU(NDIM, NPHASE, U_SNLOC), SLOC2_NU(NDIM, NPHASE, U_SNLOC) )
    ALLOCATE( SUF_U_BC_ALL(NDIM, NPHASE, U_SNLOC) )
    ALLOCATE( GLOB_T( NPHASE, CV_NONODS), GLOB_DEN( NPHASE, CV_NONODS) )

    CVNORMX_ALL(1)= CVNORMX(GI)
    IF (NDIM>=2 ) CVNORMX_ALL(2)=CVNORMY(GI)
    IF (NDIM==3 ) CVNORMX_ALL(3)=CVNORMZ(GI)

    ! coefficients for this element ELE
    UGI_COEF_ELE_ALL=0.0

    ! coefficients for this element ELE2
    UGI_COEF_ELE2_ALL=0.0


    ! Copy local variables----Start
    DO U_KLOC = 1, U_NLOC
       U_NODK = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_KLOC )
       U_NODK_IPHA = U_NODK+(IPHASE-1)*U_NONODS
       !DO IPHASE=1,NPHASE
          DO IDIM = 1, NDIM
             IF ( IDIM==1 ) THEN
                LOC_NU(IDIM,IPHASE,U_KLOC)=NU(U_NODK_IPHA)
             END IF
             IF ( IDIM==2 ) THEN
                LOC_NU(IDIM,IPHASE,U_KLOC)=NV(U_NODK_IPHA)
             END IF
             IF ( IDIM==3 ) THEN
                LOC_NU(IDIM,IPHASE,U_KLOC)=NW(U_NODK_IPHA)
             END IF
          END DO
       !END DO
    END DO

    DO U_SKLOC = 1, U_SNLOC
       U_KLOC = U_SLOC2LOC( U_SKLOC )
       U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
       U_NODK_IPHA = U_NODK + ( IPHASE - 1 ) * U_NONODS
       IF ( ELE2==0 ) THEN
          U_KLOC2 = U_KLOC
          U_NODK2 = U_NODK
       ELSE
          U_KLOC2 = U_OTHER_LOC( U_SKLOC ) 
          U_NODK2 = U_NDGLN( (ELE2-1)*U_NLOC+U_KLOC2 )
       END IF
       U_NODK2_IPHA = U_NODK2 + ( IPHASE - 1 ) * U_NONODS
       U_SNODK = ( SELE - 1 ) * U_SNLOC + U_SKLOC 
       U_SNODK_IPHA= U_SNODK + ( IPHASE - 1 ) * STOTEL*U_SNLOC
       !DO IPHASE=1,NPHASE
          DO IDIM = 1, NDIM
             IF ( IDIM==1 ) THEN
                SLOC_NU(IDIM,IPHASE,U_SKLOC)=NU(U_NODK_IPHA)
                SLOC2_NU(IDIM,IPHASE,U_SKLOC)=NU(U_NODK2_IPHA)
                IF( SELE /= 0 ) SUF_U_BC_ALL(IDIM,IPHASE,U_SKLOC)=SUF_U_BC(U_SNODK_IPHA)
             END IF
             IF ( IDIM==2 ) THEN
                SLOC_NU(IDIM,IPHASE,U_SKLOC)=NV(U_NODK_IPHA)
                SLOC2_NU(IDIM,IPHASE,U_SKLOC)=NV(U_NODK2_IPHA)
                IF( SELE /= 0 ) SUF_U_BC_ALL(IDIM,IPHASE,U_SKLOC)=SUF_V_BC(U_SNODK_IPHA)
             END IF
             IF ( IDIM==3 ) THEN
                SLOC_NU(IDIM,IPHASE,U_SKLOC)=NW(U_NODK_IPHA)
                SLOC2_NU(IDIM,IPHASE,U_SKLOC)=NW(U_NODK2_IPHA)
                IF( SELE /= 0 ) SUF_U_BC_ALL(IDIM,IPHASE,U_SKLOC)=SUF_W_BC(U_SNODK_IPHA)
             END IF
          END DO
       !END DO
    END DO

    CV_NODI = CV_NODI_IPHA - CV_NONODS*(IPHASE-1)
    CV_NODJ = CV_NODJ_IPHA - CV_NONODS*(IPHASE-1)
    GLOB_T( IPHASE, CV_NODI ) = T(CV_NODI_IPHA)
    GLOB_T( IPHASE, CV_NODJ ) = T(CV_NODJ_IPHA)
    GLOB_DEN( IPHASE, CV_NODI ) = DEN(CV_NODI_IPHA)
    GLOB_DEN( IPHASE, CV_NODJ ) = DEN(CV_NODJ_IPHA)
    ! Copy local variables----End

    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain. 
       IF( WIC_U_BC( SELE + ( IPHASE - 1 ) * STOTEL) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
          UDGI_ALL = 0.0
          DO U_KLOC = 1, U_NLOC
                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * LOC_NU( :, IPHASE, U_KLOC )
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) = UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) + 1.0
          END DO
       ELSE ! Specified vel bc.
          UDGI_ALL = 0.0
          DO U_SKLOC = 1, U_SNLOC
             U_KLOC = U_SLOC2LOC( U_SKLOC )

             IF (WIC_U_BC(SELE+(IPHASE-1)*STOTEL) == 10) THEN
                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * 0.5 * &
                                     ( SLOC_NU( :, IPHASE, U_KLOC ) + SUF_U_BC_ALL( :, IPHASE, U_SKLOC ) )
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) = UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) + 0.5   
             ELSE
                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * SUF_U_BC_ALL(:, IPHASE, U_SKLOC)
             END IF

          END DO
       END IF

    ELSE ! Conditional_SELE. Not on the boundary of the domain.
       
       UDGI_ALL = 0.0
       DO U_KLOC = 1, U_NLOC
          UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * LOC_NU( :, IPHASE, U_KLOC )
          UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) = UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) + 1.0
       END DO

       Conditional_ELE2: IF( ELE2 /= 0 ) THEN
          UDGI2_ALL = 0.0
          DO U_SKLOC = 1, U_SNLOC
             U_KLOC = U_SLOC2LOC( U_SKLOC )
             U_KLOC2 = U_OTHER_LOC( U_KLOC )
!             IF( U_KLOC2 /= 0 ) THEN ! MAKE SURE WE DONT NEED THIS...
             UDGI2_ALL(:, IPHASE) = UDGI2_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * SLOC2_NU(:, IPHASE, U_SKLOC)
             UGI_COEF_ELE2_ALL( :, IPHASE, U_KLOC2) = UGI_COEF_ELE2_ALL( :, IPHASE, U_KLOC2) + 1.0
!             ENDIF
          END DO

          IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
             DT_I=1.0
             DT_J=1.0
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
             DT_I=MAX(1.E-2,GLOB_T( IPHASE, CV_NODI ))
             DT_J=MAX(1.E-2,GLOB_T( IPHASE, CV_NODJ ))
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
             DT_I=GLOB_DEN( IPHASE, CV_NODI )*GLOB_T( IPHASE, CV_NODI )
             DT_J=GLOB_DEN( IPHASE, CV_NODJ )*GLOB_T( IPHASE, CV_NODJ )
          ENDIF
          ! Amend weighting for porosity only across elements...
          IF(ABS(CV_DG_VEL_INT_OPT ) >= 2) THEN 
             IF(ELE /= ELE2) THEN 
                DT_I=VOLFRA_PORE(ELE) *DT_I
                DT_J=VOLFRA_PORE(ELE2)*DT_J
             ENDIF
          ENDIF

          UDGI_INT_ALL(:,IPHASE) = (DT_I * UDGI_ALL(:,IPHASE) + DT_J * UDGI2_ALL(:, IPHASE)) / (DT_I + DT_J)

          IF( CV_DG_VEL_INT_OPT < 0 ) THEN

             NDOTQ_INT = DOT_PRODUCT( CVNORMX_ALL, UDGI_INT_ALL( :, IPHASE ) )

             IF( NDOTQ_INT <= 0.0 ) THEN  !Incoming
                !   DT_I=1.0
                DT_J=DT_I+DT_J
             ELSE
                DT_I=DT_I+DT_J
                !   DT_J=1.0
             ENDIF

             UDGI_INT_ALL(:,IPHASE) = (DT_I * UDGI_ALL(:,IPHASE) + DT_J * UDGI2_ALL(:, IPHASE)) / (DT_I + DT_J)

          ENDIF

          UDGI_ALL( :, IPHASE ) = UDGI_INT_ALL( :, IPHASE )

          UGI_COEF_ELE_ALL(:,IPHASE,:)=DT_I * UGI_COEF_ELE_ALL(:,IPHASE,:) / (DT_I + DT_J)

          UGI_COEF_ELE2_ALL(:,IPHASE,:)=DT_J * UGI_COEF_ELE2_ALL(:,IPHASE,:) / (DT_I + DT_J)

       ENDIF Conditional_ELE2

    ENDIF Conditional_SELE

    NDOTQ =  DOT_PRODUCT( CVNORMX_ALL, UDGI_ALL( :, IPHASE ) )

    ! Define whether flux is incoming or outgoing, depending on direction of flow
    !IF( NDOTQ >  0.0 ) THEN
    !   INCOME = 0.0  !Outgoing
    !ELSE
    !   INCOME = 1.0  !Incoming
    !ENDIF
    INCOME = 0.5*( 1. + SIGN(1.0, -NDOTQ) )

    ! Passing back at the moment
    UGI_COEF_ELE = UGI_COEF_ELE_ALL(1,IPHASE,:)
    IF (NDIM>=2 ) VGI_COEF_ELE = UGI_COEF_ELE_ALL(2,IPHASE,:)
    IF (NDIM==3 ) WGI_COEF_ELE = UGI_COEF_ELE_ALL(3,IPHASE,:)
    UGI_COEF_ELE2 = UGI_COEF_ELE2_ALL(1,IPHASE,:)
    IF (NDIM>=2 ) VGI_COEF_ELE2 = UGI_COEF_ELE2_ALL(2,IPHASE,:)
    IF (NDIM==3 ) WGI_COEF_ELE2 = UGI_COEF_ELE2_ALL(3,IPHASE,:)

    DEALLOCATE( UGI_COEF_ELE_ALL )
    DEALLOCATE( UGI_COEF_ELE2_ALL )
    DEALLOCATE( UDGI_ALL, UDGI2_ALL, UDGI_INT_ALL )
    DEALLOCATE( LOC_NU )
    DEALLOCATE( SLOC_NU, SLOC2_NU )
    DEALLOCATE( SUF_U_BC_ALL )
    DEALLOCATE( GLOB_T, GLOB_DEN )

    RETURN  

  END SUBROUTINE GET_INT_VEL_ORIG


  Pure SUBROUTINE GET_INT_T_DEN(FVT, FVT2, FVD, LIMD, LIMT, LIMT2, &
       LIMDT, LIMDTT2,  &
       FEMDGI, FEMTGI,FEMT2GI, &
       CV_DISOPT, CV_NONODS, NPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, ELE, ELE2, GI, IPHASE,  &
       CV_NLOC, TOTELE, CV_NDGLN, CV_OTHER_LOC, SCVNGI, SCVFEN, INCOME, &
       T, T2, DEN, FEMT, FEMT2, FEMDEN, &
       TMIN, T2MIN, DENMIN, &
       TMAX, T2MAX, DENMAX, &
       SELE, CV_SNLOC, STOTEL, CV_SLOC2LOC, SUF_T_BC, SUF_T2_BC, SUF_D_BC, &
       WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
       WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, VOLFRA_PORE, &
       MASS_CV, TMIN_NOD, TMAX_NOD, &
       DENMIN_NOD, DENMAX_NOD, &
       T2MIN_NOD, T2MAX_NOD, IGOT_T2, &
       TMIN_2ND_MC, T2MIN_2ND_MC, DENMIN_2ND_MC, &
       TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, LIMIT_USE_2ND, &
       HDC, NDOTQ, DT, &
       SCVFENX, SCVFENY, SCVFENZ, CVNORMX, CVNORMY, CVNORMZ, &
       U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
       IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
       TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT, &
       !Values to limit the flow when reaching the irreducible  saturation for a phase
       s_gc, s_or )
    !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
    IMPLICIT NONE
    ! Calculate T and DEN on the CV face at quadrature point GI.
    REAL, intent( inout ) :: FVT, FVT2, FVD, LIMD,LIMT,LIMT2, &
         LIMDT,LIMDTT2, &
         FEMDGI, FEMTGI, FEMT2GI
    REAL, intent( in ) :: INCOME,HDC,NDOTQ,DT, s_gc, s_or
    logical, intent( in ) :: LIMIT_USE_2ND
    INTEGER, intent( in ) :: CV_DISOPT,CV_NONODS,NPHASE,CV_NODI_IPHA,CV_NODJ_IPHA,ELE,ELE2,  &
         CV_NLOC,TOTELE,SCVNGI,GI,IPHASE,SELE,CV_SNLOC,STOTEL, &
         WIC_T_BC_DIRICHLET,WIC_D_BC_DIRICHLET, IGOT_T2, U_NLOC,U_NONODS,NDIM
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_T_BC, WIC_D_BC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_T2_BC
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFEN
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFENX, SCVFENY, SCVFENZ
    REAL, DIMENSION( :  ), intent( in ) :: CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( :  ), intent( in ) :: T, DEN, FEMT, FEMDEN
    REAL, DIMENSION( : ), intent( in ) :: T2, FEMT2
    REAL, DIMENSION( :  ), intent( inout ) :: TMIN, DENMIN, TMAX, DENMAX
    REAL, DIMENSION( : ), intent( inout ) :: T2MIN, T2MAX
    REAL, DIMENSION( : ), intent( in ) :: MASS_CV
    INTEGER, DIMENSION( :  ), intent( in ) :: TMIN_NOD, TMAX_NOD, DENMIN_NOD, DENMAX_NOD
    INTEGER, DIMENSION( :), intent( in ) :: T2MIN_NOD, T2MAX_NOD
    REAL, DIMENSION( :  ), intent( in ) :: SUF_T_BC, SUF_D_BC
    REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC
    REAL, DIMENSION( :  ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( : ), intent( in ) :: TMAX_2ND_MC, TMIN_2ND_MC,&
         DENMAX_2ND_MC, DENMIN_2ND_MC
    REAL, DIMENSION( : ), intent( in ) :: T2MAX_2ND_MC, T2MIN_2ND_MC
    REAL, DIMENSION( : ), intent( in ) :: U, V, W
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    REAL, DIMENSION( : , : , : ), intent( in ) :: INV_JAC

    INTEGER, intent( in ) :: IANISOTROPIC
    INTEGER, intent( in ) :: NSMALL_COLM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT, DENUPWIND_MAT
    REAL, DIMENSION( : ), intent( in ) :: T2UPWIND_MAT
    ! Local variables
    ! If UPWIND then use upwind flux between elements else use central. 
    ! If HI_ORDER_HALF then use high order interpolation when around 
    ! a volume frac of 0.5 and gradually apply limiting near 0 and 1. 
    LOGICAL, PARAMETER :: UPWIND = .TRUE., HI_ORDER_HALF = .FALSE., LIM_VOL_ADJUST2 = .TRUE.
    LOGICAL :: DOWNWIND_EXTRAP ! Extrapolate a downwind value for interface tracking.

    ! Scaling to reduce the downwind bias(=1downwind, =0central)
    LOGICAL, PARAMETER :: SCALE_DOWN_WIND = .true.
    ! Non-linear Petrov-Galerkin option for interface tracking...
    ! =4 is anisotropic downwind diffusion based on a projected 1D system (1st recommend)
    ! =0 is anisotropic downwind diffusion based on a velocity projection like SUPG 
    ! (2nd recommend, most compressive)
    ! =2 is isotropic downwind diffusion  (3rd recommend,least compressive)
    ! =5 is isotropic downwind diffusion with magnitude of =0 option. 
    ! In tests they all produce similar results.
    !      INTEGER, PARAMETER :: NON_LIN_PETROV_INTERFACE = 5
    INTEGER, PARAMETER :: NON_LIN_PETROV_INTERFACE = 3

    LOGICAL :: FIRSTORD, NOLIMI, RESET_STORE, LIM_VOL_ADJUST
    REAL :: RELAX, TMIN_STORE, TMAX_STORE, TOLDMIN_STORE, &
         T2MIN_STORE, T2MAX_STORE, &
         DENMIN_STORE, DENMAX_STORE, &
         COURANT_OR_MINUS_ONE_NEW
    INTEGER :: CV_KLOC, CV_NODK, CV_NODK_IPHA, CV_KLOC2, CV_NODK2, CV_NODK2_IPHA, CV_STAR_IPHA, &
         CV_SKLOC, CV_SNODK, CV_SNODK_IPHA, U_KLOC,U_NODK,U_NODK_IPHA, IDIM, ELE_DOWN
    INTEGER, DIMENSION( CV_NLOC ) :: E_CV_NODK_IPHA
    INTEGER, DIMENSION( U_NLOC )  ::E_U_NODK_IPHA
    integer, DIMENSION( : ), ALLOCATABLE :: SE_CV_NODK_IPHA, SE_CV_SNODK_IPHA
    REAL :: T_BETWEEN_MIN, T_BETWEEN_MAX
    REAL :: T_AVE_EDGE, T_AVE_ELE
    REAL :: T_MIDVAL
    REAL :: T_UPWIND, TMIN_UPWIND,TMAX_UPWIND
    REAL :: RSHAPE,RGRAY, UDGI,VDGI,WDGI, RSCALE, TXGI,TYGI,TZGI
    REAL :: VEC_VEL(3), VEC_VEL2(3), ELE_LENGTH_SCALE
    REAL :: TGI,TDTGI, U_DOT_GRADT_GI, COEF, A_STAR_T, A_STAR_X, A_STAR_Y, A_STAR_Z
    REAL :: RESIDGI, P_STAR, DIFF_COEF, COEF2, FEMTGI_DDG, income2

    ! The adjustment method is not ready for the LIMIT_USE_2ND - the new limiting method.
    LIM_VOL_ADJUST = ( LIM_VOL_ADJUST2 .AND. (.NOT.LIMIT_USE_2ND) ) .AND. IANISOTROPIC==0

    FVT2    = 1.0

    ! Extrapolate a downwind value for interface tracking.
    DOWNWIND_EXTRAP = ( cv_disopt>=8 )

    if ( DOWNWIND_EXTRAP ) then
       courant_or_minus_one_new = abs ( dt * ndotq / hdc )
    else
       courant_or_minus_one_new = -1.0
    end if

    ! Figure out global node numbers
    E_CV_NODK_IPHA = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + 1: ELE*CV_NLOC ) + ( IPHASE - 1 ) * CV_NONODS
    E_U_NODK_IPHA = U_NDGLN( ( ELE - 1 ) * U_NLOC + 1: ELE*U_NLOC ) + ( IPHASE - 1 ) * U_NONODS

    IF ( SELE /=0 ) THEN
       ALLOCATE( SE_CV_NODK_IPHA( CV_SNLOC ), SE_CV_SNODK_IPHA( CV_SNLOC ) )
       SE_CV_NODK_IPHA  = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_SLOC2LOC ) + ( IPHASE - 1 ) * CV_NONODS
       DO CV_SKLOC = 1, CV_SNLOC
          SE_CV_SNODK_IPHA( CV_SKLOC ) = ( SELE - 1 ) * CV_SNLOC + CV_SKLOC + ( IPHASE - 1 ) * STOTEL * CV_SNLOC
       END DO
    END IF

    IF ( SELE == 0 ) THEN ! Is NOT on boundary of the domain

       FVT    = INCOME * T( CV_NODJ_IPHA ) + ( 1. - INCOME ) * T( CV_NODI_IPHA )

       FVD    = INCOME * DEN( CV_NODJ_IPHA ) + ( 1. - INCOME ) * DEN( CV_NODI_IPHA )

       IF ( IGOT_T2 == 1 ) THEN
          FVT2    = INCOME * T2( CV_NODJ_IPHA ) + ( 1. - INCOME ) * T2( CV_NODI_IPHA )
       END IF

    ELSE ! Is on boundary of the domain

       IF( WIC_T_BC( SELE + ( IPHASE - 1 ) * STOTEL ) /= WIC_T_BC_DIRICHLET ) THEN 
          ! Dont apply a Dirichlet b.c.
          FVT    = T( CV_NODI_IPHA )

          FVD    = DEN( CV_NODI_IPHA )

          IF ( IGOT_T2 == 1 ) THEN
             FVT2    = T2( CV_NODI_IPHA ) 
          END IF
       ELSE
          !DO CV_SKLOC = 1, CV_SNLOC
          !   CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
          !   CV_NODK_IPHA = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_KLOC ) + ( IPHASE - 1 ) * CV_NONODS
          !   IF ( CV_NODK_IPHA == CV_NODI_IPHA ) THEN
          !      CV_SNODK = ( SELE - 1 ) * CV_SNLOC + CV_SKLOC 
          !      CV_SNODK_IPHA = CV_SNODK + ( IPHASE - 1 ) * STOTEL*CV_SNLOC
          !      EXIT
          !   END IF
          !END DO
          DO CV_SKLOC = 1, CV_SNLOC
             IF ( SE_CV_NODK_IPHA( CV_SKLOC ) == CV_NODI_IPHA ) EXIT
          END DO
          CV_SNODK_IPHA = SE_CV_SNODK_IPHA( CV_SKLOC )

          FVT    = INCOME * SUF_T_BC(CV_SNODK_IPHA) + ( 1. - INCOME ) * T( CV_NODI_IPHA )
          FVD    = INCOME * SUF_D_BC(CV_SNODK_IPHA) + ( 1. - INCOME ) * DEN( CV_NODI_IPHA )
          IF ( IGOT_T2 == 1 ) THEN
             FVT2    = INCOME * SUF_T2_BC(CV_SNODK_IPHA) + ( 1. - INCOME ) * T2( CV_NODI_IPHA )
          END IF

       END IF
       IF( WIC_D_BC( SELE + ( IPHASE - 1 ) * STOTEL ) /= WIC_D_BC_DIRICHLET ) THEN 
          FVD    = DEN( CV_NODI_IPHA )
       END IF
    END IF

    ! By default do not use first-order upwinding
    FIRSTORD = .FALSE.

    ! No limiting if CV_DISOPT is 6 or 7  (why not just define limt=femt and skip to assembly?)
    NOLIMI = ( INT( CV_DISOPT / 2 ) == 3 ) 

    ! Make a guess at the CV face value of advected field variable and density 
    ! (Depends on discetisation option, CV_DISOPT)                
    SELECT CASE( CV_DISOPT / 2 )

    CASE( 0 ) ! First-order upwinding
       FEMTGI    = FVT
       FEMT2GI    = FVT2
       FEMDGI    = FVD

       FIRSTORD = .TRUE.  

    CASE( 1 ) ! Central differencing [Trapezoidal rule (2 OR 3)]
       FEMTGI    = 0.5 * ( T( CV_NODI_IPHA ) + T( CV_NODJ_IPHA ) )        ! need to fix this!
       IF ( IGOT_T2 == 1 ) THEN
          FEMT2GI   = 0.5 * ( T2( CV_NODI_IPHA ) + T2( CV_NODJ_IPHA ) )   ! need to fix this!
       ELSE
          FEMT2GI    = 1.0
       END IF
       FEMDGI    = 0.5 * ( DEN( CV_NODI_IPHA ) + DEN( CV_NODJ_IPHA ) )    ! need to fix this!

    CASE DEFAULT ! Finite element approximation (4 OR 5)(6 or 7)(8 or 9)
       FEMTGI    = 0.0
       FEMT2GI   = 0.0
       FEMDGI    = 0.0

       Conditional_CV_DISOPT_ELE2: IF ( SELE /= 0 ) THEN
          ! Is on boundary of the domain

          IF ( WIC_T_BC(SELE+(IPHASE-1)*STOTEL) /= WIC_T_BC_DIRICHLET ) THEN ! Don't apply a Dirichlet bc
             !CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
             !CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
             FEMTGI    = dot_product( SCVFEN(:, GI ) , FEMT( E_CV_NODK_IPHA ))
             IF ( IGOT_T2 == 1 ) THEN
                FEMT2GI    = dot_product( SCVFEN(:, GI ) , FEMT2( E_CV_NODK_IPHA ))
             ELSE
                FEMT2GI    = 1.0
             END IF
          ELSE
             DO CV_SKLOC = 1, CV_SNLOC
                CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
                CV_NODK_IPHA = SE_CV_NODK_IPHA( CV_SKLOC )
                CV_SNODK_IPHA = SE_CV_SNODK_IPHA( CV_SKLOC )
                FEMTGI = FEMTGI +  SCVFEN( CV_KLOC, GI ) * ( SUF_T_BC( CV_SNODK_IPHA ) & 
                     * INCOME + FEMT( CV_NODK_IPHA ) * ( 1. - INCOME ) )
                IF ( IGOT_T2 == 1 ) THEN
                   FEMT2GI = FEMT2GI +  SCVFEN( CV_KLOC, GI ) * ( SUF_T2_BC( CV_SNODK_IPHA ) & 
                        * INCOME + FEMT2( CV_NODK_IPHA ) * ( 1. -INCOME ) )
                ELSE
                   FEMT2GI    = 1.0
                END IF
             END DO
          END IF

          IF ( WIC_D_BC(SELE+(IPHASE-1)*STOTEL) /= WIC_D_BC_DIRICHLET) THEN ! Dont apply a Dirichlet bc
             FEMDGI    = dot_product( SCVFEN( :, GI ) , FEMDEN( E_CV_NODK_IPHA ) )
          ELSE
             DO CV_SKLOC = 1, CV_SNLOC
                CV_KLOC = CV_SLOC2LOC( CV_SKLOC )
                CV_NODK_IPHA = SE_CV_NODK_IPHA( CV_SKLOC )
                CV_SNODK_IPHA = SE_CV_SNODK_IPHA( CV_SKLOC )
                FEMDGI = FEMDGI + SCVFEN( CV_KLOC, GI ) * ( SUF_D_BC( CV_SNODK_IPHA ) &
                     * INCOME + FEMDEN( CV_NODK_IPHA ) * ( 1. - INCOME ))
             END DO
          END IF

       ELSE IF( ( ELE2 == 0 ) .OR. ( ELE2 == ELE ) ) THEN

          RSCALE = 1.0 ! Scaling to reduce the downwind bias(=1downwind, =0central)
          IF ( SCALE_DOWN_WIND ) THEN

             IF ( DOWNWIND_EXTRAP .AND. (courant_or_minus_one_new.GE.0.0) ) THEN

                IF ( NON_LIN_PETROV_INTERFACE == 0 ) THEN ! NOT non-linear Petrov-Galerkin Interface

                   TXGI=0.0
                   TYGI=0.0
                   TZGI=0.0
                   DO CV_KLOC = 1, CV_NLOC
                      CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )
                      TXGI = TXGI + SCVFENX( CV_KLOC, GI ) * FEMT(CV_NODK_IPHA)
                      IF(NDIM.GE.2) TYGI = TYGI + SCVFENY( CV_KLOC, GI ) * FEMT(CV_NODK_IPHA)
                      IF(NDIM.GE.3) TZGI = TZGI + SCVFENZ( CV_KLOC, GI ) * FEMT(CV_NODK_IPHA)
                   END DO

                   UDGI=0.0
                   VDGI=0.0
                   WDGI=0.0
                   DO U_KLOC = 1, U_NLOC
                      U_NODK_IPHA = E_U_NODK_IPHA( U_KLOC )
                      UDGI = UDGI + SUFEN( U_KLOC, GI ) * U(U_NODK_IPHA)
                      IF(NDIM.GE.2) VDGI = VDGI+SUFEN( U_KLOC, GI ) * V(U_NODK_IPHA)
                      IF(NDIM.GE.3) WDGI = WDGI+SUFEN( U_KLOC, GI ) * W(U_NODK_IPHA)
                   END DO

                   ! no cosine rule:
                   RSCALE = 1.0 / PTOLFUN( SQRT( UDGI**2 + VDGI**2 + WDGI**2 ) )

                   VEC_VEL(1) = UDGI
                   VEC_VEL(2) = VDGI
                   VEC_VEL(3) = WDGI
                   VEC_VEL2 = 0.0
                   DO IDIM = 1, NDIM
                      VEC_VEL2(IDIM) = SUM( INV_JAC(IDIM, 1:NDIM, GI) * VEC_VEL(1:NDIM) )
                   END DO
                   ! normalize the velocity in here: 
                   !VEC_VEL2=VEC_VEL2/TOLFUN(SQRT( UDGI**2+VDGI**2+WDGI**2))

                   ELE_LENGTH_SCALE = 0.5 * SQRT( (UDGI**2+VDGI**2+WDGI**2) / PTOLFUN( SUM( VEC_VEL2(1:NDIM)**2 ) ) )

                   ! For discontinuous elements half the length scale...
                   IF(U_NONODS==CV_NONODS) ELE_LENGTH_SCALE=0.5*ELE_LENGTH_SCALE
                   ! For quadratic elements...
                   IF( ((NDIM==2).AND.(CV_NLOC==6)).or.((NDIM==3).AND.(CV_NLOC==10)) ) &
                        ELE_LENGTH_SCALE=0.5*ELE_LENGTH_SCALE

                ELSE ! Petrov-Galerkin end of IF(NON_LIN_PETROV_INTERFACE==0) THEN 


                      !CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                      !CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                      !CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )
                      !if ( cv_nonods == u_nonods ) then ! DG
                      if ( .true. ) then 
                         TXGI = dot_product(SCVFENX( : , GI ) , FEMT(E_CV_NODK_IPHA))
                         IF(NDIM.GE.2) then
                            TYGI = dot_product( SCVFENY( : , GI ) , FEMT(E_CV_NODK_IPHA))
                         else
                            TYGI=0.0
                         end IF
                         IF(NDIM.GE.3) then
                            TZGI = dot_product( SCVFENZ( : , GI ) , FEMT(E_CV_NODK_IPHA))
                         else
                            TZGI=0.0
                         end IF
                         TGI = dot_product(SCVFEN( : , GI ) , FEMT(E_CV_NODK_IPHA))
!!$                      else
!!$                         TXGI = TXGI + SCVFENX( CV_KLOC, GI ) * T(CV_NODK_IPHA)
!!$                         IF(NDIM.GE.2) TYGI=TYGI + SCVFENY( CV_KLOC, GI ) * T(CV_NODK_IPHA)
!!$                         IF(NDIM.GE.3) TZGI=TZGI + SCVFENZ( CV_KLOC, GI ) * T(CV_NODK_IPHA)
!!$                         TGI = TGI + SCVFEN( CV_KLOC, GI ) * T(CV_NODK_IPHA)
                      endif


                   TDTGI = 0.0

                   

                   UDGI = dot_product(SUFEN( : , GI ) , U(E_U_NODK_IPHA))
                   IF(NDIM.GE.2) then
                      VDGI = dot_product(SUFEN( : , GI ) , V(E_U_NODK_IPHA))
                   else
                      VDGI=0.0
                   end IF
                   IF(NDIM.GE.3) then
                      WDGI = dot_product(SUFEN( : , GI ) , W(E_U_NODK_IPHA))
                   else
                      WDGI=0.0
                   end IF

                   U_DOT_GRADT_GI = TDTGI + UDGI*TXGI + VDGI*TYGI + WDGI*TZGI

                   IF ( NON_LIN_PETROV_INTERFACE == 5 ) THEN 
                      COEF = 1.0 / PTOLFUN( SQRT( TDTGI**2 + TXGI**2 + TYGI**2 + TZGI**2 ) )
                   ELSE
                      COEF = U_DOT_GRADT_GI / PTOLFUN( TDTGI*TDTGI + TXGI*TXGI + TYGI*TYGI + TZGI*TZGI )
                   END IF

                   A_STAR_T = COEF * TDTGI
                   A_STAR_X = COEF * TXGI
                   A_STAR_Y = COEF * TYGI
                   A_STAR_Z = COEF * TZGI

                   RESIDGI = SQRT ( UDGI*UDGI + VDGI*VDGI + WDGI*WDGI ) / HDC

                   VEC_VEL(1) = A_STAR_X
                   VEC_VEL(2) = A_STAR_Y
                   VEC_VEL(3) = A_STAR_Z
                   VEC_VEL2=0.0

                   VEC_VEL2(1:NDIM) = matmul( INV_JAC(:,:,GI),VEC_VEL(1:NDIM) )

                   P_STAR = 0.5 * HDC / PTOLFUN( SQRT( A_STAR_X*A_STAR_X + A_STAR_Y*A_STAR_Y + A_STAR_Z*A_STAR_Z ) )

                   ! For discontinuous elements half the length scale...
                   !IF(U_NONODS==CV_NONODS) P_STAR=0.5*P_STAR 
                   ! For quadratic elements...
                   !IF( ((NDIM==2).AND.(CV_NLOC==6)).or.((NDIM==3).AND.(CV_NLOC==10)) ) &
                   !     P_STAR=0.5*P_STAR

                   select case (NON_LIN_PETROV_INTERFACE)
                   case ( 1 )     ! standard approach
                      DIFF_COEF = COEF * P_STAR * RESIDGI
                   case ( 2 )     ! standard approach making it +ve
                      DIFF_COEF = MAX( 0.0, COEF * P_STAR * RESIDGI )
                   case ( 3 )     ! residual squared approach
                      DIFF_COEF = P_STAR * RESIDGI**2 / PTOLFUN( TDTGI*TDTGI + TXGI*TXGI + TYGI*TYGI + TZGI*TZGI ) 
                   case ( 4 )     ! anisotropic diffusion in the A* direction.
                      COEF2 = CVNORMX(GI)*A_STAR_X + CVNORMY(GI)*A_STAR_Y + CVNORMZ(GI)*A_STAR_Z
                   case default   ! isotropic diffusion with u magnitide
                      DIFF_COEF = SQRT( UDGI**2 + VDGI**2 + WDGI**2 ) * P_STAR
                   END select
                   ! Make the diffusion coefficient negative (compressive)
                   DIFF_COEF = -DIFF_COEF
                   RSCALE = 1. / TOLFUN( CVNORMX(GI)*UDGI + CVNORMY(GI)*VDGI + CVNORMZ(GI)*WDGI )
 
                END IF ! Petrov-Galerkin end of IF(NON_LIN_PETROV_INTERFACE==0) THEN 
            
             END IF ! DOWNWIND_EXTRAP .AND. CFL>=0
         
          END IF ! SCALE_DOWN_WIND

          DO CV_KLOC = 1, CV_NLOC
             CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )

             IF( DOWNWIND_EXTRAP .AND. (courant_or_minus_one_new.GE.0.0) ) THEN ! Extrapolate to the downwind value...
             
                IF ( NON_LIN_PETROV_INTERFACE.NE.0 ) THEN
                   if ( .false. ) then
                      RGRAY = 0.5 * (udgi**2+vdgi**2+wdgi**2) * &
                           ( CVNORMX(GI)*SCVFENX( CV_KLOC, GI ) + CVNORMY(GI)*SCVFENY( CV_KLOC, GI ) + CVNORMZ(GI)*SCVFENZ( CV_KLOC, GI ) ) &
                           / tolfun( hdc * abs(u_dot_gradt_gi) * sqrt( TXGI**2 + TYGI**2 + TZGI**2 ) &
                           * ( CVNORMX(GI)*UDGI + CVNORMY(GI)*VDGI + CVNORMZ(GI)*WDGI )  )
                   else
                      IF ( NON_LIN_PETROV_INTERFACE == 4 ) THEN ! anisotropic diffusion...
                         RGRAY = RSCALE * COEF2 * P_STAR * ( UDGI*SCVFENX( CV_KLOC, GI ) &
                              + VDGI*SCVFENY( CV_KLOC, GI ) + WDGI*SCVFENZ( CV_KLOC, GI ) )
                      ELSE
                         RGRAY = - DIFF_COEF * RSCALE * ( CVNORMX(GI)*SCVFENX( CV_KLOC, GI ) &
                              + CVNORMY(GI)*SCVFENY( CV_KLOC, GI ) + CVNORMZ(GI)*SCVFENZ( CV_KLOC, GI ) )
                      END IF
                   end if
                ELSE
                   RGRAY = RSCALE * ELE_LENGTH_SCALE * ( UDGI*SCVFENX( CV_KLOC, GI ) &
                        + VDGI*SCVFENY( CV_KLOC, GI ) + WDGI*SCVFENZ( CV_KLOC, GI ) )
                END IF
                
                RSHAPE    = SCVFEN( CV_KLOC, GI ) + RGRAY
                FEMTGI    = FEMTGI     +  RSHAPE     * FEMT( CV_NODK_IPHA )
             ELSE

                FEMTGI    = FEMTGI     +  SCVFEN( CV_KLOC, GI ) * FEMT( CV_NODK_IPHA )             
             END IF

             FEMDGI    = FEMDGI     +  SCVFEN( CV_KLOC, GI ) * FEMDEN( CV_NODK_IPHA )
             IF(IGOT_T2==1) THEN
                FEMT2GI    = FEMT2GI     +  SCVFEN( CV_KLOC, GI ) * FEMT2( CV_NODK_IPHA )
             ELSE
                FEMT2GI    = 1.0
             END IF
          END DO

       ELSE  ! DG saturation across elements

          IF ( UPWIND ) THEN

             ! Interface tracking...
             RSCALE=1.0 ! Scaling to reduce the downwind bias(=1downwind, =0central)
             IF( DOWNWIND_EXTRAP .AND. (courant_or_minus_one_new.GE.0.0) ) THEN

                ELE_DOWN = ELE
                IF ( INCOME < 0.5 ) ELE_DOWN = ELE2

                IF ( NON_LIN_PETROV_INTERFACE == 0 ) THEN ! NOT non-linear Petrov-Galerkin Interface

                   TXGI=0.0
                   TYGI=0.0
                   TZGI=0.0
                   DO CV_KLOC = 1, CV_NLOC
                      !CV_NODK = CV_NDGLN( ( ELE_DOWN - 1 ) * CV_NLOC + CV_KLOC )
                      !CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                      CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )
                      TXGI = TXGI + SCVFENX( CV_KLOC, GI ) * FEMT(CV_NODK_IPHA)
                      IF(NDIM.GE.2) TYGI = TYGI + SCVFENY( CV_KLOC, GI ) * FEMT(CV_NODK_IPHA)
                      IF(NDIM.GE.3) TZGI = TZGI + SCVFENZ( CV_KLOC, GI ) * FEMT(CV_NODK_IPHA)
                   END DO

                   UDGI=0.0
                   VDGI=0.0
                   WDGI=0.0
                   DO U_KLOC = 1, U_NLOC
                      !U_NODK = U_NDGLN( ( ELE_DOWN - 1 ) * U_NLOC + U_KLOC )
                      !U_NODK_IPHA = U_NODK + ( IPHASE - 1 ) * U_NONODS
                      U_NODK_IPHA = E_U_NODK_IPHA( U_KLOC )
                      UDGI = UDGI + SUFEN( U_KLOC, GI )*U(U_NODK_IPHA)
                      IF(NDIM.GE.2) VDGI = VDGI + SUFEN( U_KLOC, GI ) * V(U_NODK_IPHA)
                      IF(NDIM.GE.3) WDGI = WDGI + SUFEN( U_KLOC, GI ) * W(U_NODK_IPHA)
                   END DO

                   ! cosine rule with velocity and normal:
                   !RSCALE=ABS(CVNORMX(GI)*UDGI+CVNORMY(GI)*VDGI+CVNORMZ(GI)*WDGI) &
                   !      /TOLFUN(UDGI**2+VDGI**2+WDGI**2)
                   ! cosine rule with velocity and concentration gradient:
                   !RSCALE=ABS(TXGI*UDGI+TYGI*VDGI+TZGI*WDGI) &
                   !      /TOLFUN((UDGI**2+VDGI**2+WDGI**2)*SQRT(TXGI**2+TYGI**2+TZGI**2))
                   ! no cosine rule:
                   RSCALE = 1.0 &
                        / PTOLFUN( SQRT( UDGI**2 + VDGI**2 + WDGI**2 ) )

                   VEC_VEL(1) = UDGI
                   VEC_VEL(2) = VDGI
                   VEC_VEL(3) = WDGI
                   VEC_VEL2 = 0.0
                   DO IDIM=1,NDIM
                      VEC_VEL2(IDIM) = SUM( INV_JAC(IDIM, 1:NDIM, GI) * VEC_VEL(1:NDIM) )
                   END DO
                   ! normalize the velocity in here: 
                   !VEC_VEL2=VEC_VEL2/TOLFUN(SQRT( UDGI**2+VDGI**2+WDGI**2))

                   ELE_LENGTH_SCALE = 0.5 * SQRT( (UDGI**2 + VDGI**2 + WDGI**2 ) / PTOLFUN( SUM( VEC_VEL2(1:NDIM)**2 ) ) )

                   ! For discontinuous elements half the length scale...
                   IF(U_NONODS==CV_NONODS) ELE_LENGTH_SCALE=0.5*ELE_LENGTH_SCALE
                   ! For quadratic elements...
                   IF( ((NDIM==2).AND.(CV_NLOC==6)).or.((NDIM==3).AND.(CV_NLOC==10)) ) &
                        ELE_LENGTH_SCALE = 0.5 * ELE_LENGTH_SCALE
                
                ELSE ! Petrov-Galerkin end of IF(NON_LIN_PETROV_INTERFACE==0) THEN 

                   TXGI = dot_product( SCVFENX( : , GI ) , FEMT(E_CV_NODK_IPHA))
                   IF(NDIM.GE.2) then
                      TYGI = dot_product( SCVFENY( : , GI ) , FEMT(E_CV_NODK_IPHA))
                   else
                      TYGI=0.0
                   end IF
                   IF(NDIM.GE.3) then
                      TZGI = dot_product( SCVFENZ( : , GI ) , FEMT(E_CV_NODK_IPHA))
                   else
                      TZGI=0.0
                   end IF
                   TGI = dot_product(SCVFEN( : , GI ) , FEMT(E_CV_NODK_IPHA) )

                   TDTGI=0.0

                   UDGI = dot_product( SUFEN( : , GI ) , U(E_U_NODK_IPHA))
                   IF(NDIM.GE.2) then
                      VDGI = dot_product(  SUFEN( : , GI ) , V(E_U_NODK_IPHA))
                   else
                      VDGI=0.0
                   end IF
                   IF(NDIM.GE.3) then
                      WDGI = dot_product(  SUFEN( : , GI)  , W(E_U_NODK_IPHA) )
                   else
                      WDGI =0.0
                   end IF

                   U_DOT_GRADT_GI = TDTGI + UDGI*TXGI + VDGI*TYGI + WDGI*TZGI

                   COEF = U_DOT_GRADT_GI / PTOLFUN( TDTGI*TDTGI + TXGI*TXGI + TYGI*TYGI + TZGI*TZGI )
                   IF ( NON_LIN_PETROV_INTERFACE == 5 ) THEN 
                      COEF = 1.0 / TOLFUN( SQRT( TDTGI**2 + TXGI**2 + TYGI**2 + TZGI**2 ) )
                   ELSE
                      COEF = U_DOT_GRADT_GI / TOLFUN( TDTGI**2 + TXGI**2 + TYGI**2 + TZGI**2 )
                   ENDIF

                   A_STAR_T = COEF * TDTGI
                   A_STAR_X = COEF * TXGI
                   A_STAR_Y = COEF * TYGI
                   A_STAR_Z = COEF * TZGI

                   RESIDGI = SQRT( A_STAR_X**2 + A_STAR_Y**2 + A_STAR_Z**2) / HDC

                   VEC_VEL(1) = A_STAR_X
                   VEC_VEL(2) = A_STAR_Y
                   VEC_VEL(3) = A_STAR_Z
                   VEC_VEL2 = 0.0
                   
                   VEC_VEL2(:) = matmul(INV_JAC( : , : , GI) , VEC_VEL(:) )

                   P_STAR = 0.5 * HDC / PTOLFUN( SQRT( A_STAR_X*A_STAR_X + A_STAR_Y*A_STAR_Y + A_STAR_Z*A_STAR_Z ) )

                   ! For discontinuous elements half the length scale...
                   !IF(U_NONODS==CV_NONODS) P_STAR=0.5*P_STAR 
                   ! For quadratic elements...
                   !IF( ((NDIM==2).AND.(CV_NLOC==6)).or.((NDIM==3).AND.(CV_NLOC==10)) ) P_STAR=0.5*P_STAR

                   IF ( NON_LIN_PETROV_INTERFACE == 1 ) THEN              ! standard approach
                      DIFF_COEF = COEF * P_STAR * RESIDGI
                   ELSE IF( NON_LIN_PETROV_INTERFACE == 2 ) THEN          ! standard approach making it +ve
                      DIFF_COEF = MAX( 0.0, COEF * P_STAR * RESIDGI )
                   ELSE IF( NON_LIN_PETROV_INTERFACE == 3 ) THEN          ! residual squared approach
                      DIFF_COEF = P_STAR * RESIDGI**2 / PTOLFUN( TDTGI**2 + TXGI**2 + TYGI**2 + TZGI**2 ) 
                   ELSE IF( NON_LIN_PETROV_INTERFACE == 4 ) THEN          ! anisotropic diffusion in the A* direction.
                      COEF2 = CVNORMX(GI)*A_STAR_X + CVNORMY(GI)*A_STAR_Y + CVNORMZ(GI)*A_STAR_Z
                   ELSE                                                   ! isotropic diffusion with u magnitide
                      DIFF_COEF = SQRT( UDGI**2 + VDGI**2 + WDGI**2 ) * P_STAR
                   END IF
                  
                   ! Make the diffusion coefficient negative (compressive)
                   DIFF_COEF = -DIFF_COEF
                   RSCALE = 1. / TOLFUN( CVNORMX(GI)*UDGI + CVNORMY(GI)*VDGI + CVNORMZ(GI)*WDGI )

                END IF ! Petrov-Galerkin end of IF(NON_LIN_PETROV_INTERFACE==0) THEN 

                FEMTGI = 0.0
                FEMDGI = 0.0
                FEMT2GI = 0.0
                DO CV_KLOC = 1, CV_NLOC
                   !CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                   !CV_NODK = CV_NDGLN( ( ELE_DOWN - 1 ) * CV_NLOC + CV_KLOC )
                   !CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                   CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )
                   !CV_NODK2 = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 )
                   !CV_NODK2_IPHA = CV_NODK2 + ( IPHASE - 1 ) * CV_NONODS
                   ! Extrapolate to the downwind value...
                   IF ( NON_LIN_PETROV_INTERFACE .NE. 0 ) THEN 
                      IF ( NON_LIN_PETROV_INTERFACE == 4 ) THEN ! anisotropic diffusion...
                         RGRAY = RSCALE * COEF2 *P_STAR * ( UDGI*SCVFENX( CV_KLOC, GI ) &
                              + VDGI*SCVFENY( CV_KLOC, GI ) + WDGI*SCVFENZ( CV_KLOC, GI ) )
                      ELSE
                         RGRAY = -DIFF_COEF * RSCALE * ( CVNORMX(GI)*SCVFENX( CV_KLOC, GI ) &
                              + CVNORMY(GI)*SCVFENY( CV_KLOC, GI ) + CVNORMZ(GI)*SCVFENZ( CV_KLOC, GI ) )
                      END IF
                   ELSE
                      RGRAY = RSCALE * ELE_LENGTH_SCALE * ( UDGI*SCVFENX( CV_KLOC, GI ) &
                           + VDGI*SCVFENY( CV_KLOC, GI ) + WDGI*SCVFENZ( CV_KLOC, GI ) )
                   END IF
                   RSHAPE     = SCVFEN( CV_KLOC, GI ) + RGRAY

                   !IF ( NON_LIN_PETROV_INTERFACE .NE. 0 ) THEN 
                   !   FEMTGI    = FEMTGI     +  RSHAPE     * 0.5 * ( FEMT( CV_NODK_IPHA ) + FEMT( CV_NODK2_IPHA ) )
                   !   FEMTOLDGI = FEMTOLDGI  +  RSHAPE_OLD * 0.5 * ( FEMTOLD( CV_NODK_IPHA ) + FEMTOLD( CV_NODK2_IPHA ) )
                   !ELSE
                   FEMTGI    = FEMTGI     +  RSHAPE     * FEMT( CV_NODK_IPHA )
                   !END IF

                   FEMDGI    = FEMDGI     +  SCVFEN( CV_KLOC, GI ) * FEMDEN( CV_NODK_IPHA )
                   IF ( IGOT_T2 == 1 ) THEN
                      FEMT2GI    = FEMT2GI     +  SCVFEN( CV_KLOC, GI ) * FEMT2( CV_NODK_IPHA )

                   ELSE
                      FEMT2GI    = 1.0
                   END IF
                END DO

                IF( .TRUE. ) THEN ! Downwinding for DG...
                   FEMTGI_DDG = 0.0

                   DO CV_KLOC = 1, CV_NLOC
                      CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                      IF ( CV_KLOC2 /= 0 ) THEN
                         CV_NODK2 = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 )
                         CV_NODK2_IPHA = CV_NODK2 + ( IPHASE - 1 ) * CV_NONODS
                         !CV_NODK = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_KLOC )
                         !CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                         CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )
                         ! Extrapolate to the downwind value...
                         ! USE DOWNWINDING...
                         !FEMTGI_DDG = FEMTGI_DDG +  SCVFEN( CV_KLOC, GI ) * ( FEMT( CV_NODK2_IPHA ) & 
                         !     * (1.-INCOME) + FEMT( CV_NODK_IPHA ) * INCOME )
                         !FEMTOLDGI_DDG = FEMTOLDGI_DDG + SCVFEN( CV_KLOC, GI ) * ( FEMTOLD( CV_NODK2_IPHA ) &
                         !     * (1.-INCOMEOLD) + FEMTOLD( CV_NODK_IPHA ) * INCOMEOLD )
                         ! Central...
                         FEMTGI_DDG = FEMTGI_DDG +  SCVFEN( CV_KLOC, GI ) * ( FEMT( CV_NODK2_IPHA ) & 
                              * 0.5 + FEMT( CV_NODK_IPHA ) * 0.5 )
                        
                         if ( .false. ) then
                            FEMTGI = FEMTGI +  SCVFEN( CV_KLOC, GI ) * ( FEMT( CV_NODK2_IPHA ) & 
                                 * (1.-INCOME) + FEMT( CV_NODK_IPHA ) * INCOME )
                            FEMDGI = FEMDGI + SCVFEN( CV_KLOC, GI ) * ( FEMDEN( CV_NODK2_IPHA ) &
                                 * (1.-INCOME) + FEMDEN( CV_NODK_IPHA ) * INCOME )
                            IF(IGOT_T2==1) THEN
                               FEMT2GI = FEMT2GI +  SCVFEN( CV_KLOC, GI ) * ( FEMT2( CV_NODK2_IPHA ) & 
                                    * (1.-INCOME) + FEMT2( CV_NODK_IPHA ) * INCOME )
                            ELSE
                               FEMT2GI    = 1.0
                            ENDIF
                         end if

                      END IF
                   END DO
                   
                   ! use the method with the max difference like ENO but opposite...
                   IF ( ABS(FEMTGI_DDG-FVT) .GT. ABS(FEMTGI-FVT) ) FEMTGI = FEMTGI_DDG
                   !FEMTGI   = 0.5 * ( FEMTGI_DDG + FEMTGI )  
                   FEMTGI = FEMTGI_DDG
                   !IF ( ABS(FEMTOLDGI_DDG-FVTOLD) .GT. ABS(FEMTOLDGI-FVTOLD) ) FEMTOLDGI=FEMTOLDGI_DDG
                   !FEMTOLDGI = 0.5 * ( FEMTOLDGI_DDG + FEMTOLDGI )

                END IF ! ENDOF DOWNWINDING FOR DG
                ! END OF IF ( DOWNWIND_EXTRAP .AND. (courant_or_minus_one_new.GE.0.0) ) THEN ...

             ELSE ! Standard DG upwinding...

                DO CV_KLOC = 1, CV_NLOC
                   CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                   IF ( CV_KLOC2 /= 0 ) THEN
                      CV_NODK2 = CV_NDGLN(( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 )
                      CV_NODK2_IPHA = CV_NODK2 + ( IPHASE - 1 ) * CV_NONODS
                      !CV_NODK = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_KLOC )
                      !CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                      CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )
                      FEMTGI = FEMTGI +  SCVFEN( CV_KLOC, GI ) * ( FEMT( CV_NODK2_IPHA ) & 
                           * INCOME + FEMT( CV_NODK_IPHA ) * ( 1. -INCOME ) )
                      FEMDGI = FEMDGI + SCVFEN( CV_KLOC, GI ) * ( FEMDEN( CV_NODK2_IPHA ) &
                           * INCOME + FEMDEN( CV_NODK_IPHA ) * ( 1. - INCOME ) )

                      IF ( IGOT_T2 == 1 ) THEN
                         FEMT2GI = FEMT2GI +  SCVFEN( CV_KLOC, GI ) * ( FEMT2( CV_NODK2_IPHA ) & 
                              * INCOME + FEMT2( CV_NODK_IPHA ) * ( 1. -INCOME ) )
                      ELSE
                         FEMT2GI    = 1.0
                      ENDIF

                   END IF
                END DO
             ENDIF ! END OF IF(DOWNWIND_EXTRAP.AND.(courant_or_minus_one_new.GE.0.0)) THEN ELSE ...

          ELSE ! Central DG...

             DO CV_KLOC = 1, CV_NLOC
                CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                IF ( CV_KLOC2 /= 0 ) THEN
                   CV_NODK2 = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 )
                   CV_NODK2_IPHA = CV_NODK2 + ( IPHASE - 1 ) * CV_NONODS
                   !CV_NODK = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_KLOC )
                   !CV_NODK_IPHA = CV_NODK + ( IPHASE - 1 ) * CV_NONODS
                   CV_NODK_IPHA = E_CV_NODK_IPHA( CV_KLOC )
                   FEMTGI = FEMTGI +  SCVFEN( CV_KLOC, GI ) * ( FEMT( CV_NODK2_IPHA ) & 
                        * 0.5 + FEMT( CV_NODK_IPHA ) * 0.5 )
                   FEMDGI = FEMDGI + SCVFEN( CV_KLOC, GI ) * ( FEMDEN( CV_NODK2_IPHA ) &
                        * 0.5 + FEMDEN( CV_NODK_IPHA ) * 0.5 )
                   IF ( IGOT_T2 == 1 ) THEN
                      FEMT2GI = FEMT2GI +  SCVFEN( CV_KLOC, GI ) * ( FEMT2( CV_NODK2_IPHA ) & 
                           * 0.5 + FEMT2( CV_NODK_IPHA ) * 0.5 )
                   ELSE
                      FEMT2GI    = 1.0
                   ENDIF
                ENDIF
             END DO

          END IF ! IF UPWIND / Central DG

       ENDIF Conditional_CV_DISOPT_ELE2

    END SELECT

    IF ( CV_NODI_IPHA == CV_NODJ_IPHA ) THEN

       ! On the boundary of the domain and thus use the 1st order flux. 
       LIMT    = FVT
       LIMD    = FVD
       LIMT2   = FVT2

    ELSE

       CALL LIMITERS( &
            IGOT_T2, CV_NONODS, NPHASE, IPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, IANISOTROPIC, &
            LIM_VOL_ADJUST, FIRSTORD, NOLIMI, LIMIT_USE_2ND, &
            FEMTGI , FEMDGI, FEMT2GI,INCOME, &
            COURANT_OR_MINUS_ONE_NEW, &
            T, DEN, T2, &
            TMIN, DENMIN, TMAX, DENMAX, &
            T2MIN, T2MAX, MASS_CV, TMIN_NOD, TMAX_NOD, &
            DENMIN_NOD, DENMAX_NOD, &
            T2MIN_NOD, T2MAX_NOD, &
            TUPWIND_MAT, DENUPWIND_MAT, &
            T2UPWIND_MAT, &
            NSMALL_COLM, SMALL_FINDRM, SMALL_COLM, &
            TMAX_2ND_MC, TMIN_2ND_MC, &
            DENMAX_2ND_MC, DENMIN_2ND_MC, &
            T2MAX_2ND_MC, T2MIN_2ND_MC, &
            LIMT, LIMD, LIMT2 )

    END IF

    ! Amend for porosity...
    IF ( ELE2 /= 0 ) THEN 
       FVD    = 0.5 * ( VOLFRA_PORE(ELE) + VOLFRA_PORE(ELE2) ) * FVD
       LIMD   = 0.5 * ( VOLFRA_PORE(ELE) + VOLFRA_PORE(ELE2) ) * LIMD
    ELSE
       FVD    = VOLFRA_PORE(ELE) * FVD
       LIMD   = VOLFRA_PORE(ELE) * LIMD
    END IF

    IF ( HI_ORDER_HALF ) THEN
       RELAX = MIN ( 2.*ABS(FEMTGI-0.5), 1.0 )
       LIMT = RELAX * LIMT + (1.-RELAX) * FEMTGI
    END IF

    LIMDT = LIMD * LIMT

    LIMDTT2 = LIMD * LIMT * LIMT2

    IF ( SELE /=0 ) DEALLOCATE( SE_CV_NODK_IPHA, SE_CV_SNODK_IPHA )

    RETURN


    contains
      
      PURE SUBROUTINE LIMITERS( &
       IGOT_T2, CV_NONODS, NPHASE, IPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, IANISOTROPIC, &
       LIM_VOL_ADJUST, FIRSTORD, NOLIMI, LIMIT_USE_2ND, &
       FEMTGI, FEMDGI, FEMT2GI, INCOME, &
       COURANT_OR_MINUS_ONE_NEW, &
       T, DEN, T2, &
       TMIN, DENMIN, TMAX, DENMAX, &
       T2MIN, T2MAX, MASS_CV, TMIN_NOD, TMAX_NOD, &
       DENMIN_NOD, DENMAX_NOD, &
       T2MIN_NOD, T2MAX_NOD, &
       TUPWIND_MAT, DENUPWIND_MAT, &
       T2UPWIND_MAT, &
       NSMALL_COLM, SMALL_FINDRM, SMALL_COLM, &
       TMAX_2ND_MC, TMIN_2ND_MC, &
       DENMAX_2ND_MC, DENMIN_2ND_MC, &
       T2MAX_2ND_MC, T2MIN_2ND_MC, &
       LIMT, LIMD, LIMT2 )

    implicit none


    INTEGER, intent( in ) :: IGOT_T2, CV_NONODS, NPHASE, IPHASE, CV_NODI_IPHA, CV_NODJ_IPHA, IANISOTROPIC
    LOGICAL, intent( in ) :: LIM_VOL_ADJUST, FIRSTORD, NOLIMI, LIMIT_USE_2ND 
    REAL, intent( in ) :: FEMTGI, FEMDGI, FEMT2GI,INCOME  
    REAL, intent( in ) :: COURANT_OR_MINUS_ONE_NEW
   
    REAL, DIMENSION( : ), intent( in ) :: T, DEN
    REAL, DIMENSION( : ), intent( in ) :: T2 
    REAL, DIMENSION( : ), intent( inout ) :: TMIN, DENMIN, TMAX, DENMAX
    REAL, DIMENSION( : ), intent( inout ) :: T2MIN, T2MAX

    REAL, DIMENSION( : ), intent( in ) :: MASS_CV
    INTEGER, DIMENSION( :  ), intent( in ) :: TMIN_NOD, TMAX_NOD, DENMIN_NOD, DENMAX_NOD
    INTEGER, DIMENSION( : ), intent( in ) :: T2MIN_NOD, T2MAX_NOD
    REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT, DENUPWIND_MAT
    REAL, DIMENSION( : ), intent( in ) :: T2UPWIND_MAT
    INTEGER, intent( in ) :: NSMALL_COLM
    INTEGER, DIMENSION( :), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( : ), intent( in ) :: TMAX_2ND_MC, TMIN_2ND_MC, &
         DENMAX_2ND_MC, DENMIN_2ND_MC
    REAL, DIMENSION( : ), intent( in ) :: T2MAX_2ND_MC, T2MIN_2ND_MC
    REAL, intent( inout ) :: LIMT, LIMD, LIMT2

    LOGICAL :: RESET_STORE
    REAL :: TMIN_STORE, TMAX_STORE, &
         T2MIN_STORE, T2MAX_STORE, &
         DENMIN_STORE, DENMAX_STORE, &
         TUPWIN, TUPWI2, DENUPWIN, DENUPWI2, T2UPWIN, T2UPWI2
    INTEGER :: CV_NODI, CV_NODJ, COUNT, S, E

    CV_NODI = CV_NODI_IPHA - (IPHASE-1)*CV_NONODS
    CV_NODJ = CV_NODJ_IPHA - (IPHASE-1)*CV_NONODS

    IF ( LIM_VOL_ADJUST ) THEN
       RESET_STORE = .FALSE.
       CALL CAL_LIM_VOL_ADJUST( TMIN_STORE, TMIN, T, TMIN_NOD, RESET_STORE, MASS_CV, &
            CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, NPHASE, INCOME )
       CALL CAL_LIM_VOL_ADJUST( TMAX_STORE, TMAX, T, TMAX_NOD, RESET_STORE, MASS_CV, &
            CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, NPHASE, INCOME )
       CALL CAL_LIM_VOL_ADJUST( DENMIN_STORE, DENMIN, DEN, DENMIN_NOD, RESET_STORE, MASS_CV, &
            CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, NPHASE, INCOME )
       CALL CAL_LIM_VOL_ADJUST( DENMAX_STORE, DENMAX, DEN, DENMAX_NOD, RESET_STORE, MASS_CV, &
            CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, NPHASE, INCOME )

       IF(IGOT_T2==1) THEN
          CALL CAL_LIM_VOL_ADJUST( T2MIN_STORE, T2MIN, T2, T2MIN_NOD, RESET_STORE, MASS_CV, &
               CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, NPHASE, INCOME )
          CALL CAL_LIM_VOL_ADJUST( T2MAX_STORE, T2MAX, T2, T2MAX_NOD, RESET_STORE, MASS_CV, &
               CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, NPHASE, INCOME )
       END IF

    END IF

    IF ( IANISOTROPIC == 1 ) THEN
       IF( CV_NODJ /= CV_NODI ) THEN
          DO COUNT = SMALL_FINDRM(CV_NODJ), SMALL_FINDRM(CV_NODJ+1)-1
             IF ( SMALL_COLM(COUNT) == CV_NODI ) THEN
                TUPWIN = TUPWIND_MAT( COUNT )
                DENUPWIN = DENUPWIND_MAT( COUNT )
                IF ( IGOT_T2 == 1 ) T2UPWIN = T2UPWIND_MAT( COUNT )
                EXIT
             END IF
          END DO
          DO COUNT = SMALL_FINDRM(CV_NODI), SMALL_FINDRM(CV_NODI+1)-1
             IF ( SMALL_COLM(COUNT) == CV_NODJ ) THEN
                TUPWI2 = TUPWIND_MAT( COUNT )
                DENUPWI2 = DENUPWIND_MAT( COUNT )
                IF ( IGOT_T2 == 1 ) T2UPWI2 = T2UPWIND_MAT( COUNT )
                EXIT
             END IF
          END DO
       END IF
    END IF

    S = ( IPHASE - 1 ) * CV_NONODS + 1
    E = IPHASE * CV_NONODS

    CALL ONVDLIM_ALL( CV_NONODS, &
         LIMT, FEMTGI, INCOME, CV_NODI, CV_NODJ, &
         T( S:E ), TMIN( S:E ), TMAX( S:E ), &
         TMIN_2ND_MC( S:E ), TMAX_2ND_MC( S:E ), FIRSTORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_NEW, &
         IANISOTROPIC, TUPWIN, TUPWI2 )

    CALL ONVDLIM_ALL( CV_NONODS, &
         LIMD, FEMDGI, INCOME, CV_NODI, CV_NODJ, &
         DEN( S:E ), DENMIN( S:E ), DENMAX( S:E ), &
         DENMIN_2ND_MC( S:E ), DENMAX_2ND_MC( S:E ), FIRSTORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_NEW, &
         IANISOTROPIC, DENUPWIN, DENUPWI2 )

    IF ( IGOT_T2 == 1 ) THEN
       CALL ONVDLIM_ALL( CV_NONODS, &
            LIMT2, FEMT2GI, INCOME, CV_NODI, CV_NODJ, &
            T2( S:E ), T2MIN( S:E ), T2MAX( S:E ), &
            T2MIN_2ND_MC( S:E ), T2MAX_2ND_MC( S:E ), FIRSTORD, NOLIMI, LIMIT_USE_2ND, COURANT_OR_MINUS_ONE_NEW, &
            IANISOTROPIC, T2UPWIN, T2UPWI2 )

    ELSE
       LIMT2   =1.0
    END IF

    IF ( LIM_VOL_ADJUST ) THEN
       RESET_STORE = .TRUE. 
       CALL CAL_LIM_VOL_ADJUST(TMIN_STORE,TMIN,T,TMIN_NOD,RESET_STORE,MASS_CV, &
            CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
       CALL CAL_LIM_VOL_ADJUST(TMAX_STORE,TMAX,T,TMAX_NOD,RESET_STORE,MASS_CV, &
            CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

       CALL CAL_LIM_VOL_ADJUST(DENMIN_STORE,DENMIN,DEN,DENMIN_NOD,RESET_STORE,MASS_CV, &
            CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
       CALL CAL_LIM_VOL_ADJUST(DENMAX_STORE,DENMAX,DEN,DENMAX_NOD,RESET_STORE,MASS_CV, &
            CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )

       IF ( IGOT_T2 == 1 ) THEN
          CALL CAL_LIM_VOL_ADJUST(T2MIN_STORE,T2MIN,T2,T2MIN_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
          CALL CAL_LIM_VOL_ADJUST(T2MAX_STORE,T2MAX,T2,T2MAX_NOD,RESET_STORE,MASS_CV, &
               CV_NODI_IPHA,CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE, INCOME )
       END IF

    END IF

    RETURN

  END SUBROUTINE LIMITERS

  END SUBROUTINE GET_INT_T_DEN


  PURE SUBROUTINE CAL_LIM_VOL_ADJUST( TMIN_STORE, TMIN, T, TMIN_NOD, RESET_STORE, MASS_CV, &
       &                         CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, CV_NONODS, NPHASE, INCOME )
    ! Adjust TMIN to take into account different sized CV's. 
    ! if RESET_STORE then reset TMIN to orginal values.
    implicit none
    REAL, intent( in ) :: INCOME
    INTEGER, intent( in ) :: CV_NODI_IPHA, CV_NODJ_IPHA,IPHASE, CV_NONODS, NPHASE
    LOGICAL, intent( in ) :: RESET_STORE
    REAL, intent( inout ) :: TMIN_STORE
    REAL, DIMENSION( : ), intent( inout ) :: TMIN
    REAL, DIMENSION( : ), intent( in ) :: T
    INTEGER, DIMENSION( : ), intent( in ) :: TMIN_NOD
    REAL, DIMENSION( : ), intent( in ) :: MASS_CV
    ! Local variables...
    REAL DX1, DX2_MIN, COEFF
    INTEGER CV_NODI, CV_NODJ
    LOGICAL, PARAMETER :: DEF_BOUNDED = .FALSE.

    CV_NODI = CV_NODI_IPHA-(IPHASE-1)*CV_NONODS
    CV_NODJ = CV_NODJ_IPHA-(IPHASE-1)*CV_NONODS

    IF ( RESET_STORE ) THEN

       IF ( INCOME < 0.5 ) THEN
          TMIN( CV_NODI_IPHA ) = TMIN_STORE
       ELSE
          TMIN( CV_NODJ_IPHA ) = TMIN_STORE
       ENDIF

    ELSE

       DX1 = 0.5 * ( MASS_CV(CV_NODI) + MASS_CV(CV_NODJ) )

       IF ( INCOME < 0.5 ) THEN
          DX2_MIN = 0.5 * ( MASS_CV( TMIN_NOD(CV_NODI_IPHA)-(IPHASE-1)*CV_NONODS) &
               + MASS_CV( CV_NODI) )
          TMIN_STORE = TMIN(CV_NODI_IPHA)
          IF ( DEF_BOUNDED ) THEN ! This produces strictly bounded always soln
             COEFF = MIN( 1.0, (DX1/DX2_MIN) ) 
          ELSE
             COEFF = (DX1/DX2_MIN)
          ENDIF
          TMIN( CV_NODI_IPHA ) = T(CV_NODI_IPHA) + COEFF * ( TMIN_STORE - T(CV_NODI_IPHA) )
       ELSE
          DX2_MIN = 0.5 * ( MASS_CV( TMIN_NOD(CV_NODJ_IPHA)-(IPHASE-1)*CV_NONODS) &
               + MASS_CV( CV_NODJ)) 
          TMIN_STORE = TMIN(CV_NODJ_IPHA)
          IF ( DEF_BOUNDED ) THEN ! This produces strictly bounded always soln
             COEFF = MIN( 1.0, (DX1/DX2_MIN) ) 
          ELSE
             COEFF = (DX1/DX2_MIN)
          END IF
          TMIN(CV_NODJ_IPHA) = T(CV_NODJ_IPHA) + COEFF * ( TMIN_STORE - T(CV_NODJ_IPHA) )
       END IF

    END IF

    RETURN
  END SUBROUTINE CAL_LIM_VOL_ADJUST





  PURE REAL FUNCTION FACE_THETA( DT, CV_THETA, INTERFACE_TRACK, HDC, NDOTQ, LIMDT, DIFF_COEF_DIVDX, &
       T_NODJ_IPHA, T_NODI_IPHA,  &
       NDOTQOLD, LIMDTOLD, DIFF_COEFOLD_DIVDX, TOLD_NODJ_IPHA, TOLD_NODI_IPHA )
    IMPLICIT NONE
    ! Define face value of theta
    REAL, intent(in) :: DT, CV_THETA, HDC, NDOTQ, LIMDT, DIFF_COEF_DIVDX, T_NODJ_IPHA, T_NODI_IPHA,  &
         NDOTQOLD, LIMDTOLD, DIFF_COEFOLD_DIVDX, TOLD_NODJ_IPHA, TOLD_NODI_IPHA
    LOGICAL, intent(in) :: INTERFACE_TRACK
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
       IF(INTERFACE_TRACK) THEN ! For interface tracking use forward Euler as much as possible...
          !            FTHETA = MAX( 0.0, 1. - 0.125 * MIN( ABS( PINVTH ), ABS( QINVTH ))) 
          FTHETA = MAX( 0.0, 1. - 0.5 * MIN( ABS( PINVTH ), ABS( QINVTH ))) 
       ELSE ! for Crank Nickolson time stepping base scheme... 
          FTHETA = MAX( 0.5, 1. - 0.125 * MIN( ABS( PINVTH ), ABS( QINVTH ))) 
       ENDIF
    ENDIF

    FACE_THETA = FTHETA

    RETURN

  END FUNCTION FACE_THETA




    SUBROUTINE DGSIMPLNORM_ALL( NLOC, SNLOC, NDIM,  &
         XL_ALL, XSL_ALL, NORMX_ALL )
      ! Form approximate surface normal (NORMX_ALL(1),NORMX_ALL(2),NORMX_ALL(3))
      IMPLICIT NONE
      INTEGER, intent( in ) :: NLOC, SNLOC, NDIM
      REAL, DIMENSION( NDIM, NLOC ), intent( in ) :: XL_ALL
      REAL, DIMENSION( NDIM, SNLOC ), intent( in ) :: XSL_ALL
      REAL, DIMENSION( NDIM ), intent( inout ) :: NORMX_ALL
      ! Local variables
      REAL :: XC(NDIM), SXC(NDIM), NORM
      INTEGER :: ILOC, IDIM

      DO IDIM = 1, NDIM
         XC(IDIM) = SUM( XL_ALL( IDIM, : ) )

         SXC(IDIM) = SUM( XSL_ALL( IDIM, : ) )
      END DO

      NORMX_ALL = SXC / REAL( SNLOC )  - XC /  REAL( NLOC ) 

      NORM = SQRT( SUM( NORMX_ALL**2 ) )

      NORMX_ALL = NORMX_ALL / NORM

      RETURN

    END SUBROUTINE DGSIMPLNORM_ALL




  SUBROUTINE DGSIMPLNORM( ELE, SILOC2ILOC, TOTELE, NLOC, SNLOC, XONDGL, &
       X, Y, Z, XNONOD, NORMX, NORMY, NORMZ )
    ! Form approximate surface normal (NORMX,NORMY,NORMZ)
    IMPLICIT NONE
    INTEGER, intent( in ) :: ELE, TOTELE, NLOC, SNLOC, XNONOD
    INTEGER, DIMENSION( : ), intent( in ) ::  SILOC2ILOC
    INTEGER, DIMENSION( : ), intent( in ) :: XONDGL
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
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
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
    INTEGER, DIMENSION( :, : ), intent( in ) :: CV_SLOCLIST
    INTEGER, DIMENSION( :, : ), intent( inout ) :: FACE_ELE
    INTEGER, DIMENSION( : ), intent( in ) :: FINELE
    INTEGER, DIMENSION( : ), intent( in ) :: COLELE
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
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
    REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
    REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_T_BC, WIC_D_BC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_T2_BC
    INTEGER, DIMENSION( : ), intent( in ) :: FINACV
    INTEGER, DIMENSION( : ), intent( in ), target :: COLACV
    REAL, DIMENSION( : ), intent( inout ) :: TMAX, TMIN, TOLDMAX, TOLDMIN,  &
         DENMAX, DENMIN, DENOLDMAX, DENOLDMIN
    REAL, DIMENSION( : ), intent( inout ) :: T2MAX, T2MIN, T2OLDMAX, T2OLDMIN

    REAL, DIMENSION( : ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC, TOLDMAX_2ND_MC, TOLDMIN_2ND_MC,  &
         DENMAX_2ND_MC, DENMIN_2ND_MC, DENOLDMAX_2ND_MC, DENOLDMIN_2ND_MC
    REAL, DIMENSION( : ), intent( inout ) :: T2MAX_2ND_MC, T2MIN_2ND_MC, T2OLDMAX_2ND_MC, T2OLDMIN_2ND_MC

    REAL, DIMENSION( : ), intent( in ) :: T,TOLD,DEN,DENOLD
    REAL, DIMENSION( :), intent( in ) :: T2,T2OLD
    REAL, DIMENSION( : ), intent( inout ) :: MASS_CV
    INTEGER, DIMENSION( : ), intent( inout ) :: TMIN_NOD, TMAX_NOD, TOLDMIN_NOD, &
         TOLDMAX_NOD, DENMIN_NOD, DENMAX_NOD, DENOLDMIN_NOD, DENOLDMAX_NOD
    INTEGER, DIMENSION( : ), intent( inout ) :: T2MIN_NOD, T2MAX_NOD, T2OLDMIN_NOD, &
         T2OLDMAX_NOD
    ! Local variables
    INTEGER :: CV_NODI, CV_NODJ, CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, COUNT, CV_SILOC, SELE, &
         CV_INOD, SUF_CV_SI, SUF_CV_SI_IPHA, CV_INOD_IPHA
    integer, dimension(:), pointer :: cv_neigh_ptr

    Loop_CV_NODI: DO CV_NODI = 1, CV_NONODS

       cv_neigh_ptr=>colacv(finacv(cv_nodi):finacv(cv_nodi+1)-1)

       DO IPHASE = 1, NPHASE
          CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
          TMAX( CV_NODI_IPHA ) = maxval(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          TMAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(maxloc(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
          TMIN( CV_NODI_IPHA ) = minval(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          TMIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(minloc(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
          TOLDMAX( CV_NODI_IPHA ) =&
               maxval(TOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          TOLDMAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(maxloc(TOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
          TOLDMIN( CV_NODI_IPHA ) =&
               minval(TOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          TOLDMIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(minloc(TOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
          IF(IGOT_T2==1) THEN
             T2MAX( CV_NODI_IPHA ) =&
                  maxval(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
             T2MAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                  cv_neigh_ptr(maxloc(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
             T2MIN( CV_NODI_IPHA ) =&
                  minval(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
             T2MIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                  cv_neigh_ptr(minloc(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
             T2OLDMAX( CV_NODI_IPHA ) =&
                  maxval(T2OLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
             T2OLDMAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                  cv_neigh_ptr(maxloc(T2OLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
             T2OLDMIN( CV_NODI_IPHA ) =&
                  minval(T2OLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
             T2OLDMIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                  cv_neigh_ptr(minloc(T2OLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))
          ENDIF
          DENMAX( CV_NODI_IPHA ) =&
               maxval(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          DENMAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(maxloc(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
          DENMIN( CV_NODI_IPHA ) =&
               minval(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          DENMIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(minloc(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
          DENOLDMAX( CV_NODI_IPHA ) =&
               maxval(DENOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          DENOLDMAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(maxloc(DENOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
          DENOLDMIN( CV_NODI_IPHA ) =&
               minval(DENOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          DENOLDMIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
               cv_neigh_ptr(minloc(DENOLD( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS  
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
          END DO

       END DO Loop2_CV_NODI


    ENDIF

    RETURN
  END SUBROUTINE SURRO_CV_MINMAX
 
  SUBROUTINE SURRO_CV_MINMAX_1time( TMAX, TMIN, DENMAX, DENMIN, &
       T2MAX, T2MIN, &
       TMAX_2ND_MC, TMIN_2ND_MC, DENMAX_2ND_MC, DENMIN_2ND_MC, &
       T2MAX_2ND_MC, T2MIN_2ND_MC, &
       LIMIT_USE_2ND, &
       T,  T2, DEN, IGOT_T2, NPHASE, CV_NONODS, NCOLACV, FINACV, COLACV, &
       STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC,  SUF_T2_BC, SUF_D_BC, WIC_T_BC, WIC_T2_BC, WIC_D_BC, &
       WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
       WIC_D_BC_DIRICHLET, MASS_CV, TMIN_NOD, TMAX_NOD, &
       T2MIN_NOD, T2MAX_NOD, &
       DENMIN_NOD, DENMAX_NOD)
    ! For each node, find the largest and smallest value of T and 
    ! DENSITY for both the current and previous timestep, out of 
    ! the node value and all its surrounding nodes including Dirichlet b.c's.
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE,CV_NONODS, NCOLACV,STOTEL,CV_SNLOC, &
         WIC_T_BC_DIRICHLET, WIC_T_BC_DIRI_ADV_AND_ROBIN, &
         WIC_D_BC_DIRICHLET, IGOT_T2
    LOGICAL, intent( in ) :: LIMIT_USE_2ND
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
    REAL, DIMENSION( : ), intent( in ) :: SUF_T_BC, SUF_D_BC
    REAL, DIMENSION( : ), intent( in ) :: SUF_T2_BC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_T_BC, WIC_D_BC
    INTEGER, DIMENSION( : ), intent( in ) :: WIC_T2_BC
    INTEGER, DIMENSION( : ), intent( in ) :: FINACV
    INTEGER, DIMENSION( : ), intent( in ), target :: COLACV
    REAL, DIMENSION( : ), intent( inout ) :: TMAX, TMIN,  &
         DENMAX, DENMIN
    REAL, DIMENSION( : ), intent( inout ) :: T2MAX, T2MIN

    REAL, DIMENSION( : ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC,  &
         DENMAX_2ND_MC, DENMIN_2ND_MC
    REAL, DIMENSION( : ), intent( inout ) :: T2MAX_2ND_MC, T2MIN_2ND_MC

    REAL, DIMENSION( : ), intent( in ) :: T,DEN
    REAL, DIMENSION( :), intent( in ) :: T2
    REAL, DIMENSION( : ), intent( inout ) :: MASS_CV
    INTEGER, DIMENSION( : ), intent( inout ) :: TMIN_NOD, TMAX_NOD, DENMIN_NOD, DENMAX_NOD
    INTEGER, DIMENSION( : ), intent( inout ) :: T2MIN_NOD, T2MAX_NOD
    ! Local variables
    INTEGER :: CV_NODI, CV_NODJ, CV_NODI_IPHA, CV_NODJ_IPHA, IPHASE, COUNT, CV_SILOC, SELE, &
         CV_INOD, SUF_CV_SI, SUF_CV_SI_IPHA, CV_INOD_IPHA
    integer, dimension(:), pointer :: cv_neigh_ptr

    Loop_CV_NODI: DO CV_NODI = 1, CV_NONODS

       cv_neigh_ptr=>colacv(finacv(cv_nodi):finacv(cv_nodi+1)-1)

       DO IPHASE = 1, NPHASE
          CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
          TMAX( CV_NODI_IPHA ) = maxval(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          if (TMAX( CV_NODI_IPHA )==T(CV_NODI_IPHA )) THEN
             TMAX_NOD=CV_NODI_IPHA
          else
             TMAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                  cv_neigh_ptr(maxloc(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS
          end if
          TMIN( CV_NODI_IPHA ) = minval(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          TMIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                cv_neigh_ptr(minloc(T( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS

          IF(IGOT_T2==1) THEN
             T2MAX( CV_NODI_IPHA ) =&
                  maxval(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
             T2MAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                   cv_neigh_ptr(maxloc(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS
             T2MIN( CV_NODI_IPHA ) =&
                  minval(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
             T2MIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                   cv_neigh_ptr(minloc(T2( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS
          ENDIF
          DENMAX( CV_NODI_IPHA ) =&
               maxval(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          DENMAX_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                cv_neigh_ptr(maxloc(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS
          DENMIN( CV_NODI_IPHA ) =&
               minval(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  ))
          DENMIN_NOD( CV_NODI_IPHA:CV_NODI_IPHA ) =&
                cv_neigh_ptr(minloc(DEN( cv_neigh_ptr+ ( IPHASE - 1 ) * CV_NONODS  )))+ ( IPHASE - 1 ) * CV_NONODS
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

             if(IGOT_T2==1) then
                T2MAX_2ND_MC( CV_NODI_IPHA )=-1.E+10
                T2MIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
             endif
             DENMAX_2ND_MC( CV_NODI_IPHA )=-1.E+10
             DENMIN_2ND_MC( CV_NODI_IPHA )=+1.E+10
          END DO

          DO COUNT = FINACV( CV_NODI ), FINACV( CV_NODI + 1 ) - 1, 1
             CV_NODJ = COLACV( COUNT )

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
                   
                ENDIF

             END DO
          END DO

       END DO Loop2_CV_NODI


    ENDIF

    RETURN
  END SUBROUTINE SURRO_CV_MINMAX_1time



  SUBROUTINE CALC_SELE( ELE, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
       FACE_ELE, NFACE, CVFEM_ON_FACE, &
       CV_NONODS, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
       CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN ) 
    ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC for a face on the
    ! boundary of the domain
    IMPLICIT NONE
    INTEGER, intent( in ) :: ELE, CV_NONODS, CV_NLOC, U_NLOC, &
         CV_SNLOC, U_SNLOC, NFACE, CV_ILOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: U_SNDGLN
    INTEGER, DIMENSION( :, : ), intent( in ) :: FACE_ELE
    LOGICAL, DIMENSION( : ), intent( in )  :: CVFEM_ON_FACE
    INTEGER, intent( inout ) :: SELE, CV_SILOC
    INTEGER, DIMENSION( : ), intent( inout ) :: U_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( inout ) :: CV_SLOC2LOC
    ! local variables
    INTEGER :: IFACE, ELE2, SELE2, CV_JLOC, CV_JNOD, &
         U_JLOC, U_JNOD, CV_KLOC, CV_SKNOD, &
         U_KLOC, U_SKLOC, U_SKNOD, CV_SKLOC, CV_SKLOC2, I
    LOGICAL :: FOUND
    INTEGER, DIMENSION( : ), ALLOCATABLE :: LOG_ON_BOUND

    !ewrite(3,*)'In Calc_Sele'

    ALLOCATE( LOG_ON_BOUND( CV_SNLOC ) ) ; I = 1
    DO CV_JLOC = 1, CV_NLOC  
       CV_JNOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_JLOC )
       IF ( .NOT.CVFEM_ON_FACE( CV_JLOC ) ) THEN
          LOG_ON_BOUND( I ) = CV_JNOD
          I = I + 1
       END IF
    END DO

    SELE = 0
    ! What face are we on        
    DO IFACE = 1, NFACE
       ELE2 = FACE_ELE( IFACE, ELE )
       SELE2 = MAX( 0, - ELE2 )
       ELE2 = MAX( 0, + ELE2 )
       IF ( SELE2 /= 0 ) THEN
          FOUND = .TRUE.
          DO CV_SKLOC = 1, CV_SNLOC
             CV_SKNOD = CV_SNDGLN( ( SELE2 - 1 ) * CV_SNLOC + CV_SKLOC )
             DO CV_SKLOC2 = 1, CV_SNLOC
                IF ( CV_SKNOD == LOG_ON_BOUND( CV_SKLOC2 ) ) FOUND = .FALSE.
             END DO
          END DO
          IF( FOUND ) SELE = SELE2
       END IF
    END DO

    ! Calculate CV_SLOC2LOC  
    Conditional_Sele: IF ( SELE /= 0 ) THEN   
       DO CV_SKLOC = 1, CV_SNLOC  
          CV_SKNOD = CV_SNDGLN( ( SELE - 1 ) * CV_SNLOC + CV_SKLOC )
          DO CV_JLOC = 1, CV_NLOC  
             CV_JNOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_JLOC )
             IF( CV_SKNOD == CV_JNOD ) EXIT
          END DO
          CV_SLOC2LOC( CV_SKLOC ) = CV_JLOC
          IF( CV_JLOC == CV_ILOC ) CV_SILOC = CV_SKLOC
       END DO

       ! Calculate U_SLOC2LOC 
       DO U_SKLOC = 1, U_SNLOC  
          U_SKNOD = U_SNDGLN( ( SELE - 1 ) * U_SNLOC + U_SKLOC )
          DO U_JLOC = 1, U_NLOC  
             U_JNOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_JLOC )
             IF( U_SKNOD == U_JNOD ) EXIT
          END DO
          U_SLOC2LOC( U_SKLOC ) = U_JLOC
       END DO
    END IF Conditional_Sele

    DEALLOCATE( LOG_ON_BOUND )

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
    REAL, dimension(:) :: WEIGHT
    REAL , dimension (:,:) :: N,NLX
    !Local variables
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
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: JCOUNT_KLOC, JCOUNT_KLOC2, U_OTHER_LOC
    REAL, DIMENSION( : ), intent( inout ) :: CT
    REAL, DIMENSION( : ), intent( inout ) :: CT_RHS
    REAL, DIMENSION( : ), intent( in ) :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
         UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
    REAL, DIMENSION( :, : ), intent( in ) :: SUFEN
    REAL, DIMENSION( : ), intent( in ) :: SCVDETWEI, CVNORMX, CVNORMY, CVNORMZ
    REAL, DIMENSION( : ), intent( in ) :: NU, NV, NW
    REAL, DIMENSION( : ), intent( in ) :: DEN
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
    REAL, DIMENSION( : ), intent( in ) :: SATURA
    REAL, DIMENSION( : ), intent( in ) :: X
    INTEGER, DIMENSION( : ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
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

       !         IF((ABS(SAT_CV_NODI_IPHA - SAT_CV_NODJ_IPHA) < 1.E-4) &
       !              .OR.(ABS(NDOTQ - NDOTQ2) < 1.E-6)) THEN
       IF(ABS(SAT_CV_NODI_IPHA - SAT_CV_NODJ_IPHA) < 1.E-4) THEN
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
       T, TOLD, DEN, DENOLD, T2, T2OLD, &
       FEMT, FEMTOLD, FEMDEN, FEMDENOLD, FEMT2, FEMT2OLD, USE_FEMT, &
       TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, &
       T2UPWIND_MAT, T2OLDUPWIND_MAT, &
       IGOT_T2, NPHASE, CV_NONODS,CV_NLOC, X_NLOC,TOTELE, CV_NDGLN, &
       SMALL_FINDRM, SMALL_CENTRM, SMALL_COLM,NSMALL_COLM, &
       X_NDGLN, X_NONODS, NDIM, &
       X, Y, Z, XC_CV, YC_CV, ZC_CV) 
    ! For the anisotropic limiting scheme we find the upwind values
    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
    ! value for each node pair is stored in the matrices TUPWIND AND
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, &
         NSMALL_COLM, NDIM,IGOT_T2,NPHASE
    REAL, DIMENSION( : ), intent( in ) :: T,TOLD,DEN,DENOLD
    REAL, DIMENSION( :), intent( in ) :: T2,T2OLD
    REAL, DIMENSION( :), intent( in ) :: FEMT,FEMTOLD,FEMDEN,FEMDENOLD
    REAL, DIMENSION( :), intent( in ) :: FEMT2,FEMT2OLD
    LOGICAL, intent( in ) :: USE_FEMT ! Use the FEM solns rather than CV's when interpolating soln
    REAL, DIMENSION( : ), intent( inout ) :: TUPWIND_MAT, TOLDUPWIND_MAT, &
         DENUPWIND_MAT, DENOLDUPWIND_MAT
    REAL, DIMENSION( : ), intent( inout ) :: T2UPWIND_MAT, T2OLDUPWIND_MAT
    INTEGER, DIMENSION(: ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    INTEGER, DIMENSION( : ), intent( in) :: SMALL_CENTRM
    REAL, DIMENSION( : ), intent( in ) :: X,Y,Z
    REAL, DIMENSION( : ), intent( in ) :: XC_CV, YC_CV, ZC_CV
    REAL, DIMENSION(:), ALLOCATABLE :: SOL

    ! Allocate memory 
    INTEGER :: COUNT, COUNT2, CV_NOD


    ! Allocate memory and find upwind field values for limiting...
    IF(IGOT_T2.NE.0) THEN

       allocate( sol( 6*nsmall_colm*nphase) )
       sol = (/TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT/)

       ! Obtain the weights
       CALL CALC_ANISOTROP_LIM_VALS( &
            ! Caculate the upwind values stored in matrix form...
            (/T,TOLD,DEN,DENOLD,T2,T2OLD/), &
            (/FEMT,FEMTOLD,FEMDEN,FEMDENOLD,FEMT2,FEMT2OLD/), USE_FEMT, &
            SOL,  &
            NPHASE*6,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
            SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
            X_NDGLN,X_NONODS,NDIM, &
            X,Y,Z, XC_CV, YC_CV, ZC_CV )

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
            (/FEMT,FEMTOLD,FEMDEN,FEMDENOLD/), USE_FEMT, &
            SOL,  &
            NPHASE*4,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
            SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
            X_NDGLN,X_NONODS,NDIM, &
            X,Y,Z, XC_CV, YC_CV, ZC_CV )

       TUPWIND_MAT = sol( 1 : nsmall_colm*nphase )
       TOLDUPWIND_MAT = sol( 1 + nsmall_colm*nphase : 2*nsmall_colm*nphase )
       DENUPWIND_MAT = sol( 1+2*nsmall_colm*nphase : 3*nsmall_colm*nphase )
       DENOLDUPWIND_MAT = sol( 1+3*nsmall_colm*nphase : 4*nsmall_colm*nphase )

       deallocate( sol )

    ENDIF

  END SUBROUTINE CALC_ANISOTROP_LIM


  SUBROUTINE CALC_ANISOTROP_LIM_1time(&
       ! Caculate the upwind values stored in matrix form...
       T, DEN, T2, &
       FEMT, FEMDEN, FEMT2, USE_FEMT, &
       TUPWIND_MAT, DENUPWIND_MAT, &
       T2UPWIND_MAT, &
       IGOT_T2, NPHASE, CV_NONODS,CV_NLOC, X_NLOC,TOTELE, CV_NDGLN, &
       SMALL_FINDRM, SMALL_CENTRM, SMALL_COLM,NSMALL_COLM, &
       X_NDGLN, X_NONODS, NDIM, &
       X, Y, Z, XC_CV, YC_CV, ZC_CV) 
    ! For the anisotropic limiting scheme we find the upwind values
    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
    ! value for each node pair is stored in the matrices TUPWIND AND
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, &
         NSMALL_COLM, NDIM,IGOT_T2,NPHASE
    REAL, DIMENSION( : ), intent( in ) :: T,DEN
    REAL, DIMENSION( :), intent( in ) :: T2
    REAL, DIMENSION( :), intent( in ) :: FEMT,FEMDEN
    REAL, DIMENSION( :), intent( in ) :: FEMT2
    LOGICAL, intent( in ) :: USE_FEMT ! Use the FEM solns rather than CV's when interpolating soln
    REAL, DIMENSION( : ), intent( inout ) :: TUPWIND_MAT, &
         DENUPWIND_MAT
    REAL, DIMENSION( : ), intent( inout ) :: T2UPWIND_MAT
    INTEGER, DIMENSION(: ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    INTEGER, DIMENSION( : ), intent( in) :: SMALL_CENTRM
    REAL, DIMENSION( : ), intent( in ) :: X,Y,Z
    REAL, DIMENSION( : ), intent( in ) :: XC_CV, YC_CV, ZC_CV
    REAL, DIMENSION(:), ALLOCATABLE :: SOL

    ! Allocate memory 
    INTEGER :: COUNT, COUNT2, CV_NOD


    ! Allocate memory and find upwind field values for limiting...
    IF(IGOT_T2.NE.0) THEN

       allocate( sol( 3*nsmall_colm*nphase) )
       sol = (/TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT/)

       ! Obtain the weights
       CALL CALC_ANISOTROP_LIM_VALS( &
            ! Caculate the upwind values stored in matrix form...
            (/T,DEN,T2/), &
            (/FEMT,FEMDEN,FEMT2/), USE_FEMT, &
            SOL,  &
            NPHASE*3,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
            SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
            X_NDGLN,X_NONODS,NDIM, &
            X,Y,Z, XC_CV, YC_CV, ZC_CV )

       TUPWIND_MAT = sol( 1 : nsmall_colm*nphase )
       DENUPWIND_MAT = sol( 1+nsmall_colm*nphase : 2*nsmall_colm*nphase )
       T2UPWIND_MAT = sol( 1+2*nsmall_colm*nphase : 3*nsmall_colm*nphase )

       deallocate( sol )

    ELSE

       allocate( sol( 2*nsmall_colm*nphase ) )
       sol = (/TUPWIND_MAT, DENUPWIND_MAT/)

       CALL CALC_ANISOTROP_LIM_VALS( &
            ! Caculate the upwind values stored in matrix form...
            (/T,DEN/),&
            (/FEMT,FEMDEN/), USE_FEMT, &
            SOL,  &
            NPHASE*2,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
            SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
            X_NDGLN,X_NONODS,NDIM, &
            X,Y,Z, XC_CV, YC_CV, ZC_CV )

       TUPWIND_MAT = sol( 1 : nsmall_colm*nphase )
       DENUPWIND_MAT = sol( 1+nsmall_colm*nphase : 2*nsmall_colm*nphase )

       deallocate( sol )

    ENDIF

  END SUBROUTINE CALC_ANISOTROP_LIM_1time

  SUBROUTINE CALC_ANISOTROP_LIM_VALS( &
       ! Caculate the upwind values stored in matrix form...
       T, &
       FEMT, USE_FEMT, &
       TUPWIND, &
       NFIELD,NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
       SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
       X_NDGLN,X_NONODS,NDIM, &
       X,Y,Z, XC_CV, YC_CV, ZC_CV )
    ! For the anisotropic limiting scheme we find the upwind values
    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
    ! value for each node pair is stored in the matrices TUPWIND AND
    IMPLICIT NONE
    INTEGER, intent(in) :: NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, NSMALL_COLM, NFIELD,NDIM
    REAL, DIMENSION( : ), intent( in ) :: T
    REAL, DIMENSION( : ), intent( in ) :: FEMT
    LOGICAL, intent( in ) :: USE_FEMT
    REAL, DIMENSION( : ), intent( inout ) :: TUPWIND
    INTEGER, DIMENSION( :  ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( :  ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( : ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( : ), intent( in ) :: XC_CV, YC_CV, ZC_CV
    ! the centre of each CV is: XC_CV, YC_CV, ZC_CV

    ! Allocate memory for the interpolated upwind values
    real, dimension( :, : ), allocatable :: N, NLX, NLY, NLZ
    real, dimension( : ), allocatable :: WEIGHT, L1, L2, L3, L4
    integer, dimension( : ), allocatable :: SUB_NDGLNO, SUB_XNDGLNO, ndgln_p2top1
    INTEGER :: COUNT, COUNT2, NOD, SUB_TOTELE, NGI,NLOC, ELE, IL_LOC, IQ_LOC, &
         LOC_ELE, SUB_ELE, SUB_LIN_TOTELE

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
    ALLOCATE( N(NLOC,NGI), NLX(NLOC,NGI), NLY(NLOC,NGI), NLZ(NLOC,NGI) )
    ALLOCATE( WEIGHT(NGI) )
    ALLOCATE( L1(NGI), L2(NGI), L3(NGI), L4(NGI) )
    !
    ! Shape functions for triangles and tets...
    CALL TRIQUAold( L1, L2, L3, L4, WEIGHT, ndim==3, NGI )
    ! Work out the shape functions and there derivatives...
    CALL SHATRIold( L1, L2, L3, L4, WEIGHT, ndim==3, &
         &          NLOC,NGI,&
         &          N,NLX,NLY,NLZ)

    ! ******************************************************************
    ! Calculate the sub elements for quadratic element SUB_NDGLNO ... 
    IF(CV_NLOC==NLOC) THEN 
       SUB_TOTELE=TOTELE
    ELSE
       IF(NDIM==1) THEN
          sub_lin_totele=2
       ELSE IF(NDIM==2) THEN
          sub_lin_totele=4
       ELSE IF(NDIM==3) THEN 
          sub_lin_totele=8
       ENDIF
       SUB_TOTELE= sub_lin_totele * totele

       allocate( ndgln_p2top1( sub_lin_totele*nloc ) ) ; ndgln_p2top1 = 0
       call conv_quad_to_lin_tri_tet( ndgln_p2top1, nloc, cv_nloc, sub_lin_totele )

    ENDIF

    ALLOCATE( SUB_NDGLNO( SUB_TOTELE*NLOC ) )
    ALLOCATE( SUB_XNDGLNO( SUB_TOTELE*NLOC ) )

    IF ( CV_NLOC==NLOC ) THEN
       SUB_NDGLNO = CV_NDGLN
       SUB_XNDGLNO = X_NDGLN
    ELSE

       SUB_ELE=0
       DO ELE = 1, TOTELE
          DO LOC_ELE = 1, SUB_LIN_TOTELE

             SUB_ELE = SUB_ELE + 1

             DO IL_LOC = 1, NLOC
                IQ_LOC = ndgln_p2top1( (loc_ELE-1)*NLOC + IL_LOC )
                SUB_NDGLNO( (sub_ele-1)*nloc + il_loc ) = cv_ndgln( (ele-1)*cv_nloc + iq_loc )
                SUB_XNDGLNO( (sub_ele-1)*nloc + il_loc ) = x_ndgln( (ele-1)*cv_nloc + iq_loc )
             END DO

          END DO
       END DO
       deallocate( ndgln_p2top1 )
    END IF

    ! Calculate the sub elements for quadratic element SUB_NDGLNO ... 
    ! ******************************************************************
    CALL CALC_ANISOTROP_LIM_VALS2( &
         ! Caculate the upwind values stored in matrix form...
         T, &
         FEMT, USE_FEMT, &
         TUPWIND,  &
         NFIELD, NONODS, NLOC, NGI, SUB_TOTELE, SUB_NDGLNO, &
         SMALL_FINDRM,SMALL_COLM, NSMALL_COLM, &
         SUB_XNDGLNO, X_NONODS, NDIM, &
         X, Y, Z, XC_CV, YC_CV, ZC_CV, &
         N, NLX, NLY, NLZ, WEIGHT )

    DEALLOCATE( N, NLX, NLY, NLZ, L1, L2, L3, L4, &
         WEIGHT, SUB_NDGLNO, SUB_XNDGLNO )

    RETURN
  END SUBROUTINE CALC_ANISOTROP_LIM_VALS


  SUBROUTINE CALC_ANISOTROP_LIM_VALS2( &
       ! Caculate the upwind values stored in matrix form...
       T, &
       FEMT, USE_FEMT, &
       TUPWIND,  &
       NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
       FINDRM,COLM,NCOLM, &
       X_NDGLN,X_NONODS,NDIM, &
       X,Y,Z, XC_CV, YC_CV, ZC_CV, &
       N,NLX,NLY,NLZ, WEIGHT )
    ! For the anisotropic limiting scheme we find the upwind values
    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
    ! value for each node pair is stored in the matrices TUPWIND AND
    IMPLICIT NONE
    INTEGER, intent(in) :: NONODS,X_NONODS,TOTELE,NLOC,NGI,NCOLM,NFIELD,NDIM
    REAL, DIMENSION( : ), intent( in ) :: T
    REAL, DIMENSION( : ), intent( in ) :: FEMT
    LOGICAL, intent( in ) :: USE_FEMT
    REAL, DIMENSION( : ), intent( inout ) :: TUPWIND
    INTEGER, INTENT(IN) :: NDGLNO(TOTELE*NLOC),X_NDGLN(TOTELE*NLOC)
    INTEGER, INTENT(IN) :: FINDRM(NONODS+1),COLM(NCOLM)

    REAL, DIMENSION(:), intent( in ) :: X,Y,Z
    REAL, DIMENSION( : ), intent( in ) :: XC_CV, YC_CV, ZC_CV
    !Local variables
    REAL, INTENT(IN) :: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL, INTENT(IN) :: WEIGHT(NGI) 

    INTEGER, DIMENSION( : ), ALLOCATABLE, SAVE :: ELEMATPSI
    REAL, DIMENSION( :  ), ALLOCATABLE, SAVE :: ELEMATWEI
    LOGICAL, SAVE :: STORE_ELE=.TRUE., RET_STORE_ELE=.FALSE.

    LOGICAL :: D3,DCYL
    ! Allocate memory for the interpolated upwind values
    LOGICAL, PARAMETER :: BOUND  = .TRUE., REFLECT = .FALSE. ! limiting options
    INTEGER, DIMENSION( : ), allocatable :: NOD_FINDELE,NOD_COLELE, NLIST, INLIST, DUMMYINT
    REAL, DIMENSION( : ), allocatable :: DUMMYREAL
    INTEGER MXNCOLEL,NCOLEL,adapt_time_steps
    REAL current_time

    ! Over-estimate the size of the COLELE array
    MXNCOLEL=20*TOTELE+500

    ALLOCATE( NOD_FINDELE(X_NONODS+1) )
    ALLOCATE( NOD_COLELE(MXNCOLEL) )

    ALLOCATE( NLIST(X_NONODS) )
    ALLOCATE( INLIST(X_NONODS) )

    ! Calculate node element list - moved from (I)FINPTS
    CALL PHILNODELE(X_NONODS,NOD_FINDELE,NOD_COLELE, &
         NCOLEL,MXNCOLEL, &
         TOTELE,NLOC,X_NDGLN, &
         NLIST,INLIST)

    IF( STORE_ELE ) THEN

       ALLOCATE( ELEMATPSI( NCOLM ) )
       ALLOCATE( ELEMATWEI( NCOLM * NLOC ) )

       CALL FINPTSSTORE(T,FEMT,USE_FEMT,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
            TUPWIND,FINDRM,COLM,NCOLM,NDIM, &
            X_NDGLN,X_NONODS, &
            X,Y,Z, XC_CV, YC_CV, ZC_CV,&
            N,NLX,NLY,NLZ, WEIGHT, &
            NOD_FINDELE,NOD_COLELE,NCOLEL, &
            ELEMATPSI,ELEMATWEI,1, &
            BOUND, REFLECT)

    ELSE IF( RET_STORE_ELE ) THEN 
       
       ! Find the weights for the interpolation
       ! This does depend on the solns T when BOUND...
       CALL GETSTOREELEWEI(T,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
            TUPWIND,FINDRM,COLM,NCOLM,BOUND, &
            ELEMATPSI,ELEMATWEI)

    ELSE 

       ! Assume we have not stored anything (elements or weights)...
       ALLOCATE(DUMMYINT(NCOLM))
       ALLOCATE(DUMMYREAL(NCOLM*NLOC))

       CALL FINPTSSTORE(T,FEMT,USE_FEMT,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
            TUPWIND,FINDRM,COLM,NCOLM,NDIM, &
            X_NDGLN,X_NONODS, &
            X,Y,Z, XC_CV, YC_CV, ZC_CV,&
            N,NLX,NLY,NLZ, WEIGHT, &
            NOD_FINDELE,NOD_COLELE,NCOLEL, &
            DUMMYINT,DUMMYREAL,0, &
            BOUND, REFLECT)

       DEALLOCATE(DUMMYINT,DUMMYREAL) 

    ENDIF

    store_ele = .false. ; ret_store_ele = .true.
    if( have_option( '/mesh_adaptivity/hr_adaptivity') ) then
       if( have_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps') ) then
          call get_option( '/mesh_adaptivity/hr_adaptivity/period_in_timesteps', &
               adapt_time_steps )
       end if
       if( mod( timestep, adapt_time_steps ) == 0 ) store_ele = .true.
    elseif( have_option( '/mesh_adaptivity/prescribed_adaptivity' ) ) then
       call get_option( '/timestepping/current_time', current_time )
       if( do_adapt_state_prescribed( current_time ) ) store_ele = .true.
    end if
    if ( store_ele ) then
       ret_store_ele = .false.
       deallocate( elematpsi, elematwei )
    end if

    DEALLOCATE( NOD_FINDELE, NOD_COLELE, NLIST, INLIST )

  END SUBROUTINE CALC_ANISOTROP_LIM_VALS2


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
    INTEGER, intent(in) :: NFIELD,NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
    REAL, INTENT(IN) :: PSI(NONODS*NFIELD)
    INTEGER, INTENT(IN) :: NCOLM
    INTEGER, INTENT(IN) :: FINDRM(NONODS+1),COLM(NCOLM)
    REAL, INTENT(INOUT) :: MATPSI(NCOLM*NFIELD)
    INTEGER, INTENT(IN) :: ELEMATPSI(NCOLM)
    REAL, INTENT(IN) ::  ELEMATWEI(NCOLM*NLOC)
    !  LOCAL VARIABLES...
    INTEGER NOD,COUNT,ELEWIC,ILOC,INOD,IFIELD
    INTEGER KNOD,COUNT2,JNOD
    REAL RMATPSI
    REAL, ALLOCATABLE, DIMENSION(:)::MINPSI
    REAL, ALLOCATABLE, DIMENSION(:)::MAXPSI

    if ( bound ) then
       ALLOCATE( MINPSI( TOTELE*NFIELD ) )
       ALLOCATE( MAXPSI( TOTELE*NFIELD ) )

       ! find the max and min local to each element...
       CALL MINMAXELEWIC( PSI,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
            &     FINDRM,COLM,NCOLM,&
            &     MINPSI,MAXPSI )
    end if

    do NOD = 1, NONODS
       do COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
          IF(NOD.NE.COLM(COUNT)) THEN
             ELEWIC = ELEMATPSI( COUNT )
             DO IFIELD = 1, NFIELD
                RMATPSI=0.0
                DO ILOC = 1, NLOC
                   INOD = NDGLNO( (ELEWIC-1)*NLOC + ILOC )
                   RMATPSI = RMATPSI + ELEMATWEI( (COUNT-1)*NLOC+ILOC) * PSI(INOD+(IFIELD-1)*NONODS)
                END DO

                RMATPSI   =PSI(NOD+(IFIELD-1)*NONODS)   &
                     +(1./FRALINE)*(RMATPSI   -PSI(NOD+(IFIELD-1)*NONODS))

                ! make locally bounded...
                if ( bound ) then
                   MATPSI(COUNT+(IFIELD-1)*NCOLM)   &
                        =MAX(MIN(RMATPSI,   MAXPSI(ELEWIC+(IFIELD-1)*TOTELE)),   &
                        &                            MINPSI(ELEWIC+(IFIELD-1)*TOTELE))
                else
                   MATPSI(COUNT+(IFIELD-1)*NCOLM)   =RMATPSI
                end if
             END DO
          END IF
       END DO
    END DO

    RETURN

  end subroutine getstoreelewei

  SUBROUTINE MINMAXELEWIC(PSI,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
       &     FINDRM,COLM,NCOLM,&
       &     MINPSI,MAXPSI)
    ! This sub calculates the max and min values of PSI in local vacinity of 
    ! an element. 
    IMPLICIT NONE
    INTEGER, intent(in) :: NFIELD,NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
    REAL, INTENT(IN) :: PSI(NONODS*NFIELD)
    INTEGER, INTENT(IN) :: NCOLM
    INTEGER, INTENT(IN) :: FINDRM(NONODS+1),COLM(NCOLM)
    REAL, INTENT(INOUT) :: MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
    !  LOCAL VARIABLES...
    INTEGER NOD,COUNT,ELEWIC,ILOC,INOD,IFIELD
    INTEGER KNOD,COUNT2,JNOD
    REAL RMATPSI

    MINPSI   =1.E+20
    MAXPSI   =-1.E+20
    ! find the max and min local to each element...
    DO ELEWIC=1,TOTELE! Was loop 
       DO ILOC=1,NLOC! Was loop 
          KNOD=NDGLNO((ELEWIC-1)*NLOC+ILOC)
          ! Search around node KNOD for max and min PSI...
          DO COUNT2 = FINDRM(KNOD), FINDRM(KNOD+1)-1
             JNOD = COLM( COUNT2 )
             DO IFIELD = 1, NFIELD
                MINPSI( ELEWIC+(IFIELD-1)*TOTELE )  &
                     = MIN( PSI(JNOD+(IFIELD-1)*NONODS), MINPSI(ELEWIC+(IFIELD-1)*TOTELE) )
                MAXPSI( ELEWIC+(IFIELD-1)*TOTELE )  &
                     = MAX( PSI(JNOD+(IFIELD-1)*NONODS), MAXPSI(ELEWIC+(IFIELD-1)*TOTELE) )
             END DO
          END DO
       END DO
    END DO

    !ewrite(3,*) '***M-m', MAXPSI-MINPSI

    RETURN

  end subroutine minmaxelewic
        
!     
!     
!     
!     
  SUBROUTINE FINPTSSTORE(PSI,FEMPSI,USE_FEMPSI,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
       &     MATPSI,FINDRM,COLM,NCOLM,NDIM, &
       &     X_NDGLN,X_NONODS, &
       &     X,Y,Z, XC_CV, YC_CV, ZC_CV,&
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
    INTEGER, intent(in) :: NFIELD,NONODS,NLOC,NGI,TOTELE,NDIM,X_NONODS
    INTEGER, intent(in) :: NDGLNO(TOTELE*NLOC)
    REAL, intent(in) :: PSI(NONODS*NFIELD)
    REAL, intent(in) :: FEMPSI(NONODS*NFIELD)
    LOGICAL, intent(in) :: USE_FEMPSI
    INTEGER, intent(in) :: NCOLM,NCOLEL
    INTEGER, intent(in) :: FINDRM(NONODS+1),COLM(NCOLM)
    REAL, intent(inout) :: MATPSI(NCOLM*NFIELD)
    INTEGER, intent(in) :: X_NDGLN(TOTELE*NLOC)
    REAL, intent(in) :: X(X_NONODS),Y(X_NONODS),Z(X_NONODS)
    REAL, DIMENSION( : ), intent( in ) :: XC_CV, YC_CV, ZC_CV
    REAL, intent(in) :: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
    REAL, intent(in) :: WEIGHT(NGI)
    !     work space...
    INTEGER, intent(in) :: FINDELE(X_NONODS+1),COLELE(NCOLEL)
    INTEGER, intent(in) :: IGETSTOR
    INTEGER, intent(inout) :: ELEMATPSI(NCOLM*IGETSTOR)
    REAL, intent(inout) :: ELEMATWEI(NCOLM*NLOC*IGETSTOR)
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
    REAL :: X1,Y1,Z1, X2,Y2,Z2

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
          END IF
       END DO
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

    MATPSI=0.
    DO NOD=1,NONODS! Was loop 10

       XNOD=NOD2XNOD(NOD)
       !     
       DO COUNT=FINDRM(NOD ),FINDRM(NOD+1)-1! Was loop 20

          NODJ=COLM(COUNT)
          XNODJ=NOD2XNOD(NODJ)
          !     
          IF(NOD.NE.NODJ) THEN

             IF(REFLECT) THEN
                NORMX1=NORMX(NOD)
                NORMY1=NORMY(NOD)
                NORMZ1=NORMZ(NOD)
             ENDIF
             IF(NONODS.NE.X_NONODS) THEN ! Its a DG soln field...
               X1=XC_CV(NOD)
               Y1=YC_CV(NOD)
               Z1=ZC_CV(NOD)

               X2=XC_CV(NODJ)
               Y2=YC_CV(NODJ)
               Z2=ZC_CV(NODJ)
             ELSE
               X1=X(XNOD)
               Y1=Y(XNOD)
               Z1=Z(XNOD)

               X2=X(XNODJ)
               Y2=Y(XNODJ)
               Z2=Z(XNODJ)
             ENDIF
             CALL MATPTSSTORE(MATPSI,COUNT,NFIELD,NOD,XNOD,&
                  &              PSI,FEMPSI,USE_FEMPSI,NONODS,X_NONODS,&
                  &              NLOC,TOTELE,X_NDGLN,NDGLNO,&
                  &              NCOLM,&
                  &              X1,Y1,Z1,&
                  &              X2,Y2,Z2,&
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
!    stop 67

    RETURN

    end subroutine finptsstore
!     
!     
!     
!     
      SUBROUTINE MATPTSSTORE(MATPSI,COUNT,NFIELD,NOD,XNOD,&
     &     PSI,FEMPSI,USE_FEMPSI,NONODS,X_NONODS,&
     &     NLOC,TOTELE,X_NDGLN,NDGLNO,&
     &     NCOLM,&
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
      REAL INFINY,FRALINE2
      LOGICAL, intent(in) :: REFLECT
! IF REFLECT then use a reflection condition at boundary to 
! do limiting. 
      PARAMETER(INFINY=1.E+20,FRALINE2=0.001)
      LOGICAL, intent(in) :: BOUND
      INTEGER, intent(in) :: COUNT,NFIELD,NOD,XNOD,NONODS,X_NONODS,NLOC,TOTELE,NDIM
      REAL, intent(in) :: PSI(NONODS*NFIELD)
      REAL, intent(in) :: FEMPSI(NONODS*NFIELD)
      LOGICAL, intent(in) :: USE_FEMPSI
      REAL, intent(inout) :: MATPSI(NCOLM*NFIELD)
      INTEGER, intent(in) :: X_NDGLN(NLOC*TOTELE),NDGLNO(NLOC*TOTELE)
      INTEGER, intent(in) :: NCOLM
      REAL, intent(in) :: X1,Y1,Z1,X2,Y2,Z2,NORMX1,NORMY1,NORMZ1
      REAL, intent(in) :: X(X_NONODS),Y(X_NONODS),Z(X_NONODS)
      INTEGER, intent(in) :: NCOLEL
      INTEGER, intent(in) :: FINDELE(X_NONODS+1),COLELE(NCOLEL)
      REAL, intent(in) :: MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
      INTEGER, intent(inout) :: ELEWIC
      REAL, intent(inout) :: LOCCORDSK(NLOC)
!     
!     Local variables...
      REAL XC,YC,ZC
      REAL LOCCORDS(4)
      INTEGER LOCNODS(4),LOCNODSK(4)
      INTEGER NLOCNODS(4),NLOCNODSK(4)
      INTEGER ELE,ILOC,KNOD,JNOD,IFIELD, COUNT2
      REAL MINCOR,MINCORK,RSUM
      REAL VX,VY,VZ,T2X,T2Y,T2Z,T1X,T1Y,T1Z,DIST12,RN,RMATPSI
      REAL REFX,REFY,REFZ,REFX2,REFY2,REFZ2, FRALINE
      LOGICAL IS_DG

      IS_DG=NONODS.NE.X_NONODS
!     
      FRALINE=FRALINE2
      IF(IS_DG) FRALINE=1.0

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
      DO COUNT2=FINDELE(XNOD),FINDELE(XNOD+1)-1! Was loop 10
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
                 &        X(LOCNODS(4)),Y(LOCNODS(4)),Z(LOCNODS(4)) )
         ELSE
            CALL TRILOCCORDS2D(XC,YC, &
                 &        LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),&
                 !     The 3 corners of the tri...
                 &        X(LOCNODS(1)),Y(LOCNODS(1)),&
                 &        X(LOCNODS(2)),Y(LOCNODS(2)),&
                 &        X(LOCNODS(3)),Y(LOCNODS(3)) )
         END IF

         MINCOR=MINVAL( LOCCORDS(1:NLOC) )
!          print *,'ele,LOCCORDS(1:NLOC):',ele,LOCCORDS(1:NLOC)

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
!        stop 677


!     Set all the negative basis to zero and re-normalise 
!     to put on the face of an element...
      RSUM=0.0
      DO ILOC=1,NLOC! Was loop 
         LOCCORDSK(ILOC)=MAX(0.0,LOCCORDSK(ILOC))
         RSUM=RSUM+LOCCORDSK(ILOC)
      END DO
      IF(RSUM.LT.1.E-5) THEN ! Just in case RSUM=0.0
            LOCCORDSK(1:NLOC)=1.0/REAL(NLOC)
      ELSE
         DO ILOC=1,NLOC! Was loop 
            LOCCORDSK(ILOC)=LOCCORDSK(ILOC)/RSUM
         END DO
      ENDIF
!         print *,'nod,ELEWIC,LOCCORDSk(1:NLOC)=',nod,ELEWIC,LOCCORDSk(1:NLOC)
      DO IFIELD=1,NFIELD
         RMATPSI=0.0
         DO ILOC=1,NLOC! Was loop 
            IF(USE_FEMPSI) THEN
               RMATPSI   =RMATPSI  +LOCCORDSK(ILOC)*FEMPSI(NLOCNODSK(ILOC)+(IFIELD-1)*NONODS)
            ELSE
               RMATPSI   =RMATPSI  +LOCCORDSK(ILOC)*PSI(NLOCNODSK(ILOC)+(IFIELD-1)*NONODS)
            ENDIF
!         XC=XC+LOCCORDSK(ILOC)*X(LOCNODSK(ILOC))
!         YC=YC+LOCCORDSK(ILOC)*Y(LOCNODSK(ILOC))
!         ZC=ZC+LOCCORDSK(ILOC)*Z(LOCNODSK(ILOC))
         END DO
!     Exaduate difference by a factor of 100.
         IF(USE_FEMPSI) THEN
            RMATPSI   = FEMPSI( NOD + (IFIELD-1)*NONODS )  &
         + (1./FRALINE) * ( RMATPSI - FEMPSI( NOD + (IFIELD-1)*NONODS) )
         ELSE
            RMATPSI   = PSI( NOD + (IFIELD-1)*NONODS )  &
         + (1./FRALINE) * ( RMATPSI - PSI( NOD + (IFIELD-1)*NONODS) )
         ENDIF

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
      integer, dimension( : ), intent( inout ) :: FINDELE
      integer, dimension( : ), intent( inout ) :: COLELE
      integer, dimension( : ), intent( in ) :: NDGLNO
      integer, dimension( : ), intent( inout ) :: NLIST,INLIST
!     Local variables...
      INTEGER NOD,ELE,ILOC,COUNT, INOD
!     
      NLIST=0
      INLIST=0

      ! NLIST is the number of elements each node belongs to...
    !  print *,'NONODS,totele,MXNCOLEL:',NONODS,totele,MXNCOLEL
      do ELE=1,TOTELE! Was loop 
      do ILOC=1,NLOC! Was loop 
    !      print *,'iloc,nloc,totele,ele:', iloc,nloc,totele,ele
    !      print *,'NDGLNO((ELE-1)*NLOC+ILOC):',NDGLNO((ELE-1)*NLOC+ILOC)
            INOD=NDGLNO((ELE-1)*NLOC+ILOC)
            NLIST(INOD)=NLIST(INOD)+1
         END DO
      END DO
    !  stop 771

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
      DO ELE=1,TOTELE! Was loop 
         DO ILOC=1,NLOC! Was loop 
            INOD=NDGLNO((ELE-1)*NLOC+ILOC)
            INLIST(INOD)=INLIST(INOD)+1
            IF (FINDELE(INOD)-1+INLIST(INOD).GT.MXNCOLEL) THEN
               STOP 'COLELE ARRAY OUT OF BOUNDS--SUB:PHILNODELE'
            ENDIF
            COLELE(FINDELE(INOD)-1+INLIST(INOD))=ELE
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

      AREA = TRIAREAF_SIGN( X1, Y1, X2, Y2, X3, Y3)
!     area coords...

      N1 = TRIAREAF_SIGN(Xp, Yp,  &
     &     X2, Y2,  &
     &     X3, Y3 ) 

      N1 = N1/AREA

      

      N2 = TRIAREAF_SIGN(X1, Y1, &
     &     Xp, Yp,  &
     &     X3, Y3 ) 
      
      N2 = N2/AREA
      


      N3 = TRIAREAF_SIGN(X1, Y1,  &
     &     X2, Y2,  &
     &     Xp, Yp ) 

      N3 = N3/AREA
       

      Return

      end subroutine triloccords2d
! 
! 
! 
! 
     
      subroutine conv_quad_to_lin_tri_tet( ndgln_p2top1, nloc_lin, cv_nloc, sub_lin_totele )
        ! convert quadratic element into a series of linear elements...
        integer, intent( in ) :: nloc_lin, cv_nloc, sub_lin_totele
        integer, intent( inout ) :: ndgln_p2top1(sub_lin_totele*nloc_lin)
        ! local variables...
        integer :: sub_ele

        if(cv_nloc==6) then ! quadratic triangle...
           sub_ele = 1
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 1
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4

           sub_ele = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 5

           sub_ele = 3
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 3
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 5

           sub_ele = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 5
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 6

        else if(cv_nloc==10) then ! quadratic triangle...

           sub_ele = 1
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 7
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 8
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 9
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 10

           sub_ele = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 1
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 7

           sub_ele = 3
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 7
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 8
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 4

           sub_ele = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 2
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 3
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 8

           sub_ele = 5
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 3
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 5
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 8

           sub_ele = 6
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 5
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 9
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 8

           sub_ele = 7
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 5
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 6
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 4
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 9

           sub_ele = 8
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 1 ) = 7
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 2 ) = 9
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 3 ) = 8
           ndgln_p2top1( ( sub_ele - 1 ) * nloc_lin + 4 ) = 4

        else 
           ewrite(3,*) 'not a viable option for calc_sub_lin_tri_tet'
        end if

        return

      end subroutine conv_quad_to_lin_tri_tet




!
! 
! 

! -----------------------------------------------------------------------------

  end module cv_advection
      
