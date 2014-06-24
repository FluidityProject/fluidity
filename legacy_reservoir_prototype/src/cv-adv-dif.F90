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

  use fields
  use solvers

  use solvers_module
  use spud
  use global_parameters, only: option_path_len, field_name_len, timestep, is_overlapping, is_compact_overlapping
  use futils, only: int2str
  use adapt_state_prescribed_module
  use sparse_tools
  use sparsity_patterns

  use shape_functions
  use matrix_operations
  use Copy_Outof_State
  use boundary_conditions

  INTEGER, PARAMETER :: WIC_T_BC_DIRICHLET = 1, WIC_T_BC_ROBIN = 2, &
       WIC_T_BC_DIRI_ADV_AND_ROBIN = 3, WIC_D_BC_DIRICHLET = 1, &
       WIC_U_BC_DIRICHLET = 1, &
       WIC_U_BC_ROBIN = 2, &
       WIC_U_BC_DIRI_ADV_AND_ROBIN = 3, &
       WIC_U_BC_DIRICHLET_INOUT = 2, &
       WIC_P_BC_DIRICHLET = 1

contains

    SUBROUTINE CV_ASSEMB( state, packed_state, &
         tracer, velocity, density, &
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
         DEN_ALL, DENOLD_ALL, &
         MAT_NLOC, MAT_NDGLN, MAT_NONODS, TDIFFUSION, IGOT_THERM_VIS, THERM_U_DIFFUSION, &
         CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, SECOND_THETA, CV_BETA, &
         SUF_SIG_DIAGTEN_BC, &
         DERIV, CV_P, &
         SOURCT, ABSORBT, VOLFRA_PORE, &
         NDIM, GETCV_DISC, GETCT, &
         NCOLM, FINDM, COLM, MIDM, &
         XU_NLOC, XU_NDGLN, FINELE, COLELE, NCOLELE, &
         OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS, &
         DEN_FEMT, &
         IGOT_T2, T2, T2OLD, IGOT_THETA_FLUX, SCVNGI_THETA, GET_THETA_FLUX, USE_THETA_FLUX, &
         THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J, THETA_GDIFF, &
         IN_ELE_UPWIND, DG_ELE_UPWIND, &
         NOIT_DIM, &
         MEAN_PORE_CV, &
         FINDCMC, COLCMC, NCOLCMC, MASS_MN_PRES, THERMAL, RETRIEVE_SOLID_CTY, &
         MASS_ELE_TRANSP, &
         StorageIndexes, Field_selector, icomp,&
         option_path_spatial_discretisation,T_input,TOLD_input, FEMT_input,&
         saturation)

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
      type( state_type ), intent( inout ) :: packed_state
      type(tensor_field), intent(inout), target :: tracer
      type(tensor_field), intent(in), target :: density
      type(tensor_field), intent(in) :: velocity

      INTEGER, intent( in ) :: NCOLACV, NCOLCT, CV_NONODS, U_NONODS, X_NONODS, MAT_NONODS, &
           TOTELE, &
           CV_ELE_TYPE, &
           NPHASE, CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           CV_SNLOC, U_SNLOC, STOTEL, CV_DISOPT, CV_DG_VEL_INT_OPT, NDIM, &
           NCOLM, XU_NLOC, NCOLELE, NOPT_VEL_UPWIND_COEFS, &
           IGOT_T2, IGOT_THETA_FLUX, SCVNGI_THETA, IN_ELE_UPWIND, DG_ELE_UPWIND, &
           NCOLCMC, Field_selector
      INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: XU_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      INTEGER, DIMENSION(: ), intent( in ) :: U_SNDGLN
      REAL, DIMENSION( : ), intent( inout ) :: CV_RHS
      real, dimension(:), intent(inout) :: CSR_ACV
      real, dimension(:,:,:), intent(inout) :: dense_acv
      INTEGER, DIMENSION( : ), intent( in ) :: FINACV
      INTEGER, DIMENSION( : ), intent( in ) :: COLACV
      INTEGER, DIMENSION( : ), intent( in ) :: MIDACV
      REAL, DIMENSION( :, :, : ), intent( inout ), allocatable :: CT
      ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
      REAL, DIMENSION( : ), intent( inout ), allocatable :: DIAG_SCALE_PRES
      REAL, DIMENSION( :  ), intent( inout ) :: CT_RHS
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
      REAL, DIMENSION( : ), optional, intent( in ) :: T_input, TOLD_input, FEMT_input !<========TEMPORARY UNTIL ALL THE VARIABLES ARE INSIDE PACKED_STATE!!!
      REAL, DIMENSION( :, : ), intent( in ) :: DEN_ALL, DENOLD_ALL
      REAL, DIMENSION( : ), intent( in ) :: T2, T2OLD
      REAL, DIMENSION( :, : ), intent( inout ) :: THETA_GDIFF ! (NPHASE,CV_NONODS)
      REAL, DIMENSION( :, : ), intent( inout ), optional :: THETA_FLUX, ONE_M_THETA_FLUX, THETA_FLUX_J, ONE_M_THETA_FLUX_J
      REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
      INTEGER, intent( in ) :: IGOT_THERM_VIS
      !REAL, DIMENSION(NDIM,NDIM,NPHASE,MAT_NONODS*IGOT_THERM_VIS), intent( in ) :: THERM_U_DIFFUSION
      REAL, DIMENSION(:,:,:,:), intent( in ) :: THERM_U_DIFFUSION

      REAL, intent( in ) :: DT, CV_THETA, SECOND_THETA, CV_BETA
      REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
      REAL, DIMENSION( NPHASE, CV_NONODS ), intent( in ) :: DERIV
      REAL, DIMENSION( : ), intent( in ) :: CV_P
      REAL, DIMENSION( : ), intent( in ) :: SOURCT
      REAL, DIMENSION( :, :, : ), intent( in ) :: ABSORBT
      REAL, DIMENSION( : ), intent( in ) :: VOLFRA_PORE
      LOGICAL, intent( in ) :: GETCV_DISC, GETCT, GET_THETA_FLUX, USE_THETA_FLUX, THERMAL, RETRIEVE_SOLID_CTY
      INTEGER, DIMENSION( : ), intent( in ) :: FINDM
      INTEGER, DIMENSION( : ), intent( in ) :: COLM
      INTEGER, DIMENSION( : ), intent( in ) :: MIDM
      INTEGER, DIMENSION( : ), intent( in ) :: FINELE
      INTEGER, DIMENSION( : ), intent( in ) :: COLELE
      REAL, DIMENSION( : ), intent( inout ) :: DEN_FEMT!, T_FEMT
      REAL, DIMENSION( : ), intent( in ) :: OPT_VEL_UPWIND_COEFS
      INTEGER, INTENT( IN ) :: NOIT_DIM
      REAL, DIMENSION( : ), intent( inout ) :: MEAN_PORE_CV
      REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE_TRANSP
      character( len = * ), intent( in ), optional :: option_path_spatial_discretisation
      integer, dimension(:), intent(in) :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      integer, dimension(:), intent(inout) :: StorageIndexes
     integer, optional, intent(in) :: icomp
     type(tensor_field), intent(in), optional, target :: saturation
      !character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation


      ! Local variables
      REAL :: ZERO_OR_TWO_THIRDS
! if integrate_other_side then just integrate over a face when cv_nodj>cv_nodi
      logical, PARAMETER :: integrate_other_side= .true.
! if .not.correct_method_petrov_method then we can compare our results directly with previous code...
      logical, PARAMETER :: correct_method_petrov_method= .true.
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE
      LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE, U_ON_FACE, &
           CVFEM_ON_FACE, UFEM_ON_FACE
      INTEGER, DIMENSION( : ), allocatable :: &
           CV_OTHER_LOC, U_OTHER_LOC, MAT_OTHER_LOC, &
           JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, CV_SLOC2LOC, U_SLOC2LOC
      INTEGER, DIMENSION( : , : ), allocatable :: FACE_ELE
      REAL, DIMENSION( : ), allocatable ::  &
           CVNORMX, &
           CVNORMY, CVNORMZ, MASS_CV, MASS_ELE, SNDOTQ, SNDOTQOLD,  &
           SRA,   &
           SUM_CV, ONE_PORE, &
           DU, DV, DW, PERM_ELE
      REAL, DIMENSION( :, : ), allocatable :: CVNORMX_ALL, XC_CV_ALL
!      REAL, DIMENSION( :, : ), allocatable :: UGI_COEF_ELE, VGI_COEF_ELE, WGI_COEF_ELE, &
!           UGI_COEF_ELE2, VGI_COEF_ELE2, WGI_COEF_ELE2
      REAL, DIMENSION( :, :, : ), allocatable :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL


    !***Pointers for Shape function calculation***
      integer, pointer :: NCOLGPTS
      integer, dimension(:), pointer ::  FINDGPTS, COLGPTS
      integer, dimension(:,:), pointer ::  CV_NEILOC, CV_SLOCLIST, U_SLOCLIST

      real, dimension(:), pointer :: CVWEIGHT, CVWEIGHT_SHORT, SCVFEWEIGH, SBCVFEWEIGH,&
        SELE_OVERLAP_SCALE
      REAL, DIMENSION( : , :, : ), pointer :: CVFENLX_ALL, CVFENLX_SHORT_ALL, UFENLX_ALL,&
      SCVFENLX_ALL, SUFENLX_ALL, SBCVFENLX_ALL, SBUFENLX_ALL
      REAL, DIMENSION( : , : ), pointer :: CVN, CVN_SHORT, CVFEN, CVFEN_SHORT, &
           UFEN, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SUFEN, SUFENSLX, SUFENSLY,  &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBUFEN, SBUFENSLX, SBUFENSLY, &
           DUMMY_ZERO_NDIM_NDIM
    !###Pointers for Shape function calculation###

      REAL, DIMENSION( : ), allocatable :: SHAPE_CV_SNL
      REAL, DIMENSION( :, :, : ), allocatable :: DUMMY_ZERO_NDIM_NDIM_NPHASE
      REAL, DIMENSION( :, :, :, : ), allocatable :: DTX_ELE_ALL, DTOLDX_ELE_ALL
      REAL, pointer, DIMENSION( :, :, : ) :: INV_JAC
      real, pointer, dimension( : ) :: SCVDETWEI, SCVRA
      real, pointer, dimension( :, :, : ) :: SCVFENX_ALL

! Variables used to calculate CV face values: 
      REAL, DIMENSION( :, : ), allocatable :: LOC_F, LOC_FEMF
      REAL, DIMENSION( :, : ), allocatable :: SLOC_F, SLOC_FEMF, SLOC2_F, SLOC2_FEMF
      REAL, DIMENSION( :, :, : ), allocatable :: LOC_UF, SLOC_UF, SLOC2_UF
      INTEGER, DIMENSION( : ), allocatable :: SELE_LOC_WIC_F_BC
      REAL, DIMENSION( :, : ), allocatable :: SLOC_SUF_F_BC
      REAL, DIMENSION( : ), allocatable :: FUPWIND_IN, FUPWIND_OUT, FVF, LIMF, F_INCOME, F_NDOTQ

! Variables used in GET_INT_VEL_NEW: 
      REAL, DIMENSION ( :, :, : ), allocatable :: LOC_U, LOC2_U
      REAL, DIMENSION ( :, :, : ), allocatable :: LOC_NU, LOC2_NU, SLOC_NU, LOC_NUOLD, LOC2_NUOLD, SLOC_NUOLD
      REAL, DIMENSION ( :, : ), allocatable :: LOC_U_HAT, LOC2_U_HAT
      INTEGER :: CV_KNOD, CV_KNOD2, U_SNODK
      REAL, DIMENSION ( :, : ), allocatable :: LOC_FEMT, LOC2_FEMT, LOC_FEMTOLD, LOC2_FEMTOLD
      REAL, DIMENSION ( :, : ), allocatable :: LOC_FEMT2, LOC2_FEMT2, LOC_FEMT2OLD, LOC2_FEMT2OLD
      REAL, DIMENSION ( : ), allocatable :: LOC_T_I, LOC_T_J, LOC_DEN_I, LOC_DEN_J
      REAL, DIMENSION ( : ), allocatable :: LOC_TOLD_I, LOC_TOLD_J, LOC_DENOLD_I, LOC_DENOLD_J
      REAL, DIMENSION ( : ), allocatable :: LOC_T2_I, LOC_T2_J, LOC_T2OLD_I, LOC_T2OLD_J

! NPHASE Variables: 
      REAL, DIMENSION( : ), allocatable :: DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX, &
               NDOTQNEW,  NDOTQOLD, INCOMEOLD,  INCOME_J, INCOMEOLD_J, LIMT2OLD, LIMDTOLD, &
               NDOTQ, INCOME, LIMT2, LIMTOLD, LIMT, LIMT_HAT, &
               FVTOLD, FVT2OLD, FVDOLD, &
               LIMDOLD, LIMDTT2OLD,&
               FVT, FVT2, FVD, LIMD,  &
               LIMDT, LIMDTT2, SUM_LIMT, SUM_LIMTOLD
      LOGICAL :: DISTCONTINUOUS_METHOD, QUAD_ELEMENTS


      !        ===> INTEGERS <===
      INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, COUNT, ICOUNT, JCOUNT, &
           ELE, ELE2, GI, GCOUNT, SELE,   &
           U_SILOC, CV_SILOC, U_KLOC, &
           CV_ILOC, CV_JLOC, IPHASE, JPHASE, &
           CV_NODJ, CV_NODJ_IPHA, rhs_nodj_ipha,rhs_nodi_ipha,&
           CV_NODI, CV_NODI_IPHA, CV_NODI_JPHA, U_NODK, TIMOPT, &
           ICOUNT_IPHA, JCOUNT_IPHA, IMID_IPHA, JMID_IPHA, &
           NFACE, X_NODI,  &
           CV_INOD, MAT_NODI,  MAT_NODJ, FACE_ITS, NFACE_ITS, &
           CVNOD, XNOD, NSMALL_COLM, COUNT2, NOD
      !        ===>  REALS  <===
      REAL :: HDC, &
           FTHETA, VTHETA, &
           FEMDGI, FEMTGI,FEMT2GI, FEMDOLDGI, FEMTOLDGI, FEMT2OLDGI, &
           TMID, TOLDMID, TMID_J, TOLDMID_J,&
           ROBIN1, ROBIN2, &
           RSUM, &
           FTHETA_T2, ONE_M_FTHETA_T2OLD, FTHETA_T2_J, ONE_M_FTHETA_T2OLD_J, THERM_FTHETA, &
           W_SUM_ONE1, W_SUM_ONE2

      real, pointer :: VOLUME
      integer :: cv_inod_ipha, IGETCT, U_NODK_IPHA, IANISOLIM, global_face
      logical :: Have_Temperature_Fields, Have_VolumeFraction_Fields, Have_Components_Fields
      ! Functions...
      !REAL :: R2NORM, FACE_THETA
      !        ===>  LOGICALS  <===
      LOGICAL :: GETMAT, &
           D1, D3, DCYL, GOT_DIFFUS, INTEGRAT_AT_GI, &
           NORMALISE, SUM2ONE, GET_GTHETA, QUAD_OVER_WHOLE_ELE

      character( len = option_path_len ) :: option_path, option_path2, path_temp, path_volf, &
           path_comp, path_spatial_discretisation


!      real, dimension(:), allocatable :: TUPWIND_MAT, TOLDUPWIND_MAT, DENUPWIND_MAT, &
!           DENOLDUPWIND_MAT, T2UPWIND_MAT, T2OLDUPWIND_MAT
      real, dimension(:,:), allocatable :: TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, &
           DENOLDUPWIND_MAT_ALL, T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL
      INTEGER :: IDUM(1)
      REAL :: RDUM(1),n1,n2,n3

      INTEGER :: I, NLEV, U_NLOC2, ILEV, IDIM, U_ILOC, U_INOD, ELE3
      INTEGER :: NFIELD, CV_KLOC, CV_NODK, CV_NODK_IPHA
      INTEGER :: IT, ITOLD, ID, IDOLD, IT2, IT2OLD, IFI
      INTEGER :: COUNT_IN, COUNT_OUT,CV_KLOC2,CV_NODK2,CV_NODK2_IPHA,CV_SKLOC, CV_SNODK, CV_SNODK_IPHA
      INTEGER :: IPT_IN, IPT_OUT
      INTEGER :: U_KLOC2,U_NODK2,U_NODK2_IPHA,U_SKLOC
      INTEGER :: IPT,ILOOP,IMID,JMID,JDIM,IJ
      LOGICAL :: STORE, integrate_other_side_and_not_boundary, prep_stop, GOT_VIS
      REAL :: R, NDOTQ_HAT
      REAL :: LIMT_keep(NPHASE ),  LIMTOLD_keep( NPHASE ), LIMD_keep( NPHASE ),   LIMDOLD_keep( NPHASE ), LIMT2_keep( NPHASE ),   LIMT2OLD_keep(NPHASE)
      REAL , DIMENSION( : ), ALLOCATABLE :: F_CV_NODI, F_CV_NODJ
      REAL , DIMENSION( :, : ), ALLOCATABLE :: NUOLDGI_ALL, NUGI_ALL, NU_LEV_GI
      REAL , DIMENSION( :, :, :, : ), ALLOCATABLE :: VECS_STRESS, VECS_GRAD_U
      REAL , DIMENSION( :, :, : ), ALLOCATABLE :: STRESS_IJ_THERM, STRESS_IJ_THERM_J
      REAL , DIMENSION( :, :, : ), ALLOCATABLE :: VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS, &
                                                  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
                                                  INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS
      REAL , DIMENSION( :, :, :, : ), ALLOCATABLE :: INV_V_OPT_VEL_UPWIND_COEFS
      REAL :: BCZERO(NPHASE),  T_ALL_J( NPHASE ), TOLD_ALL_J( NPHASE )
      INTEGER :: LOC_WIC_T_BC_ALL(NPHASE)


      REAL, DIMENSION( :, :, : ), ALLOCATABLE :: ABSORBT_ALL!U_ALL, NU_ALL, NUOLD_ALL,
      REAL, DIMENSION( :, : ), ALLOCATABLE :: &!X_ALL, T_ALL, TOLD_ALL,
            T2_ALL, T2OLD_ALL, &!FEMT_ALL, FEMTOLD_ALL, &
           FEMDEN_ALL, FEMDENOLD_ALL, FEMT2_ALL, FEMT2OLD_ALL, SOURCT_ALL
      LOGICAL, DIMENSION( : ), ALLOCATABLE :: DOWNWIND_EXTRAP_INDIVIDUAL
      LOGICAL, DIMENSION( :, : ), ALLOCATABLE :: IGOT_T_PACK, IGOT_T_CONST
      REAL, DIMENSION( :, : ), ALLOCATABLE :: IGOT_T_CONST_VALUE
      REAL, DIMENSION( :, :, : ), ALLOCATABLE :: TEN_XX_ONE


!      INTEGER, DIMENSION( :, : ), ALLOCATABLE :: WIC_T_BC_ALL, WIC_D_BC_ALL, WIC_T2_BC_ALL

!      TYPE( TENSOR_FIELD ), POINTER :: U_s, NU_s, NUOLD_s
!      REAL, DIMENSION( : ), ALLOCATABLE :: U, V, W, NU, NUOLD, NV, NVOLD, NW, NWOLD
!      type( vector_field ), pointer:: x
      !Working pointers
      real, dimension(:), allocatable :: T, TOLD !<= TEMPORARY, TO REMOVE WHEN CONVERSION IS FINISHED
      real, dimension(:), allocatable :: VOL_FRA_FLUID ! for solid coupling
      real, dimension(:, :), allocatable :: U_HAT_ALL ! for solid coupling
      real, dimension(:,:), allocatable, target :: T_ALL_TARGET, TOLD_ALL_TARGET, FEMT_ALL_TARGET, FEMTOLD_ALL_TARGET, T_TEMP, TOLD_TEMP
      real, dimension(:,:), pointer :: T_ALL, TOLD_ALL, X_ALL, FEMT_ALL, FEMTOLD_ALL
      real, dimension(:, :, :), pointer :: comp, comp_old, fecomp, fecomp_old, U_ALL, NU_ALL, NUOLD_ALL
      real, dimension(:,:), allocatable :: T_ALL_KEEP

!! boundary_condition fields
      type(tensor_field) :: velocity_BCs,tracer_BCs, density_BCs, saturation_BCs
      type(tensor_field) :: tracer_BCs_robin2, saturation_BCs_robin2

      INTEGER, DIMENSION( 1 , nphase , surface_element_count(tracer) ) :: WIC_T_BC_ALL
      INTEGER, DIMENSION( 1 , nphase , surface_element_count(tracer) ) ::&
           WIC_D_BC_ALL
      INTEGER, DIMENSION( 1 , nphase , igot_t2*surface_element_count(tracer) ) ::&
           WIC_T2_BC_ALL
      INTEGER, DIMENSION( ndim , nphase , surface_element_count(tracer) ) :: WIC_U_BC_ALL
      REAL, DIMENSION( :,:,: ), pointer :: SUF_T_BC_ALL,&
           SUF_T_BC_ROB1_ALL, SUF_T_BC_ROB2_ALL
      REAL, DIMENSION(:,:,: ), pointer :: SUF_D_BC_ALL,&
           SUF_T2_BC_ALL, SUF_T2_BC_ROB1_ALL, SUF_T2_BC_ROB2_ALL 
      REAL, DIMENSION(:,:,: ), pointer :: SUF_U_BC_ALL

!! femdem
      type( vector_field ), pointer :: delta_u_all, us_all
      type( scalar_field ), pointer :: solid_vol_fra
      real :: theta_cty_solid, VOL_FRA_FLUID_I, VOL_FRA_FLUID_J


      type( tensor_field_pointer ), dimension(4+2*IGOT_T2) :: psi,fempsi
      type( vector_field_pointer ), dimension(1) :: PSI_AVE,PSI_INT
      type(vector_field), pointer :: coord
      type( tensor_field ), pointer :: old_tracer, old_density, old_saturation
      integer :: FEM_IT

      integer, dimension(:), pointer :: neighbours
      integer :: nb
      logical :: skip




      !#################SET WORKING VARIABLES#################
      call get_var_from_packed_state(packed_state,PressureCoordinate = X_ALL,&
           OldNonlinearVelocity = NUOLD_ALL, NonlinearVelocity = NU_ALL)
      !For every Field_selector value but 3 (saturation) we need U_ALL to be NU_ALL
      U_ALL => NU_ALL
      
      allocate( t_all_keep( nphase, cv_nonods ) )

      old_tracer=>extract_tensor_field(packed_state,GetOldName(tracer))
      old_density=>extract_tensor_field(packed_state,GetOldName(density))
      if (present(saturation)) then
         old_saturation=>extract_tensor_field(packed_state,&
              GetOldName(saturation))
      end if


      if (.not.present(T_input)) then!<==TEMPORARY
         select case (Field_selector)
         case (1)!Temperature
            call get_var_from_packed_state(packed_state,Temperature = T_ALL,&
                 OldTemperature = TOLD_ALL, FETemperature =FEMT_ALL,OldFETemperature = FEMTOLD_ALL)
            T_ALL_KEEP = T_ALL
         case (2)!Component mass fraction
            if (present(icomp)) then
               call get_var_from_packed_state(packed_state,ComponentMassFraction = comp, &
                    OldComponentMassFraction = comp_old, FEComponentMassFraction = fecomp, &
                    OldFEComponentMassFraction = fecomp_old )
               T_ALL => comp(icomp,:,:)
               TOLD_ALL => comp_old(icomp,:,:)
               FEMT_ALL => fecomp(icomp,:,:)
               FEMTOLD_ALL => fecomp_old(icomp,:,:)
               T_ALL_KEEP = T_ALL
            else
               FLAbort('Component field require to introduce icomp')
            end if
         case (3)!Saturation
            call get_var_from_packed_state(packed_state,PhaseVolumeFraction = T_ALL,&
                 OldPhaseVolumeFraction = TOLD_ALL, Velocity = U_ALL,&
                 FEPhaseVolumeFraction = FEMT_ALL, OldFEPhaseVolumeFraction = FEMTOLD_ALL)
            
            T_ALL_KEEP = T_ALL
            
            IF( GETCT ) THEN
               IF( RETRIEVE_SOLID_CTY ) THEN
                  ALLOCATE(VOL_FRA_FLUID(CV_NONODS))
                  ALLOCATE(U_HAT_ALL(NDIM,U_NONODS))
                  
                  delta_u_all => extract_vector_field( packed_state, "delta_U" )
                  u_hat_all = delta_u_all%val + u_all( :, 1, :) ! ndim, u_nonods
                  
                  us_all => extract_vector_field( packed_state, "solid_U" )
                  
                  Solid_vol_fra => extract_scalar_field( packed_state, "SolidConcentration" )
                  VOL_FRA_FLUID = 1.0 - 1.0 * solid_vol_fra%val   ! cv_nonods
                  
                  
                  ALLOCATE(T_TEMP(NPHASE,CV_NONODS), TOLD_TEMP(NPHASE,CV_NONODS))
                  
                  IF(NPHASE==1) THEN
                     do cv_inod = 1, cv_nonods
                        do iphase = 1, nphase
                           ! Amend the saturations to produce the real voln fractions -only is we have just one phase.
                           T_TEMP(iphase, cv_inod) = VOL_FRA_FLUID(cv_inod)
                           TOLD_TEMP(iphase, cv_inod) = VOL_FRA_FLUID(cv_inod)
                        end do
                     end do
                  ELSE
                     T_TEMP= T_ALL
                     TOLD_TEMP=TOLD_ALL
                  ENDIF

                  ! switch off caching of CV face values as this will be wrong.
                  T_ALL=>T_TEMP
                  TOLD_ALL=>TOLD_TEMP
                  ! CONV = A*B ! conV is an allocatable target
                  ! T_ALL=>CONV ! conV is an allocatable target
                  
                  call get_option( '/blasting/theta_cty_solid', theta_cty_solid, default=1.  )
                  
               ENDIF
            ENDIF
         case default
            FLAbort('Invalid field_selector value')
         end select

      else
         ALLOCATE( T_ALL_TARGET( NPHASE, CV_NONODS ), TOLD_ALL_TARGET( NPHASE, CV_NONODS ), FEMT_ALL_TARGET(NPHASE, CV_NONODS) )
         do cv_inod = 1, cv_nonods
            do iphase = 1, nphase
               T_ALL_TARGET(iphase, cv_inod) = T_input(cv_inod+(iphase-1)*cv_nonods)
               TOLD_ALL_TARGET(iphase, cv_inod) = TOLD_input(cv_inod+(iphase-1)*cv_nonods)
               FEMT_ALL_TARGET(iphase, cv_inod) = FEMT_input(cv_inod+(iphase-1)*cv_nonods)
            end do
         end do
         T_ALL => T_ALL_TARGET
         TOLD_ALL => TOLD_ALL_TARGET
         FEMT_ALL => FEMT_ALL_TARGET
         allocate (FEMTOLD_ALL_TARGET(NPHASE, CV_NONODS))
         FEMTOLD_ALL_TARGET = 0.
         FEMTOLD_ALL => FEMTOLD_ALL_TARGET
         
         T_ALL_KEEP = T_ALL
      end if




      allocate (T (NPHASE* CV_NONODS), TOLD(NPHASE* CV_NONODS))!TEMPORARY




     !##################END OF SET VARIABLES##################


      !! Get boundary conditions from field

      call get_entire_boundary_condition(tracer,&
           ['weakdirichlet','robin        '],&
           tracer_BCs,WIC_T_BC_ALL,boundary_second_value=tracer_BCs_robin2)
      call get_entire_boundary_condition(density,&
           ['weakdirichlet'],&
           density_BCs,WIC_D_BC_ALL)
      if (present(saturation))&
           call get_entire_boundary_condition(saturation,&
           ['weakdirichlet','robin        '],&
           saturation_BCs,WIC_T2_BC_ALL,&
           boundary_second_value=saturation_BCs_robin2)
      call get_entire_boundary_condition(velocity,&
           ['weakdirichlet'],&
           velocity_BCs,WIC_U_BC_ALL)
      
      !! reassignments to old arrays, to be discussed

      SUF_T_BC_ALL=>tracer_BCs%val
      SUF_T_BC_ROB1_ALL=>tracer_BCs%val ! re-using memory from dirichlet bc.s for Robin bc
      SUF_T_BC_ROB2_ALL=>tracer_BCs_robin2%val
      SUF_D_BC_ALL=>density_BCs%val
      SUF_U_BC_ALL=>velocity_BCs%val
      if(present(saturation)) then
         SUF_T2_BC_ALL=>saturation_BCs%val
         SUF_T2_BC_ROB1_ALL=>saturation_BCs%val ! re-using memory from dirichlet bc.s for Robin bc
         SUF_T2_BC_ROB2_ALL=>saturation_BCs_robin2%val
      end if

!      x => extract_vector_field( packed_state, "PressureCoordinate" )
!      allocate( x_all( ndim, x_nonods ) ) ; x_all=0.0
!      x_all( 1, : ) = x % val( 1, : )
!      if (ndim >=2 ) x_all( 2, : ) = x % val( 2, : )
!      if (ndim >=3 ) x_all( 3, : ) = x % val( 3, : )




     ! LOCAL VARIABLE FOR INDIRECT ADDRESSING----------------------------------------------
!     REAL, DIMENSION ( :, :, : ), allocatable :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
!     ALLOCATE( UGI_COEF_ELE_ALL( NDIM, NPHASE, U_NLOC) )
!     ALLOCATE( UGI_COEF_ELE2_ALL( NDIM, NPHASE, U_NLOC) )
     !-------------------------------------------------------------------------------------

!         print *,'just entered sub'
      IDUM = 0
      RDUM = 0.

      ewrite(3,*) 'In CV_ASSEMB'


      GOT_VIS = .FALSE. 
      IF(IGOT_THERM_VIS==1) GOT_VIS = ( R2NORM( THERM_U_DIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0 )

      GOT_DIFFUS = ( R2NORM( TDIFFUSION, MAT_NONODS * NDIM * NDIM * NPHASE ) /= 0 )

     call get_option( "/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/viscosity_scheme/zero_or_two_thirds", zero_or_two_thirds, default=2./3. )

!      print *,'SECOND_THETA=',SECOND_THETA
!         stop 2821
      ewrite(3,*)'CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, SECOND_THETA, GOT_DIFFUS:', &
           CV_DISOPT, CV_DG_VEL_INT_OPT, DT, CV_THETA, CV_BETA, SECOND_THETA, GOT_DIFFUS
      ewrite(3,*)'GETCV_DISC, GETCT', GETCV_DISC, GETCT


!      ndotq = 0. ; ndotqold = 0.

!         print *,'just entered sub -1.1'
      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

    
      ! Allocate memory for the control volume surface shape functions, etc.

      IF(GETCT) THEN
         ALLOCATE( JCOUNT_KLOC( U_NLOC ))   
         ALLOCATE( JCOUNT_KLOC2( U_NLOC ))
         ALLOCATE( ICOUNT_KLOC( U_NLOC )) 
         ALLOCATE( ICOUNT_KLOC2( U_NLOC )) 
      ENDIF


      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( ndim, cv_ele_type, cv_nloc, u_nloc, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

  
!    DISTCONTINUOUS_METHOD=.false.  
!    IF(X_NONODS.NE.CV_NONODS) DISTCONTINUOUS_METHOD=.true.  ! *****need to change this...
    DISTCONTINUOUS_METHOD = ( CV_NONODS == TOTELE * CV_NLOC )
!     print *,'DISTCONTINUOUS_METHOD,X_NONODS,CV_NONODS:',DISTCONTINUOUS_METHOD,X_NONODS,CV_NONODS
!       stop 281
! Quadratic elements
    QUAD_ELEMENTS = ( ((NDIM==2).AND.(CV_NLOC==6)).or.((NDIM==3).AND.(CV_NLOC==10)) ) 


!         print *,'just entered sub -2'
      ALLOCATE( CVNORMX( SCVNGI ))
      ALLOCATE( CVNORMY( SCVNGI ))
      ALLOCATE( CVNORMZ( SCVNGI ))
      ALLOCATE( CVNORMX_ALL( NDIM, SCVNGI )) ; CVNORMX_ALL=0.0 

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

      ALLOCATE( SHAPE_CV_SNL( CV_NLOC ))

      ALLOCATE( SRA( SCVNGI ))

      ALLOCATE( DUMMY_ZERO_NDIM_NDIM(NDIM,NDIM))
      DUMMY_ZERO_NDIM_NDIM=0.0
      ALLOCATE( DUMMY_ZERO_NDIM_NDIM_NPHASE(NDIM,NDIM,NPHASE))
      DUMMY_ZERO_NDIM_NDIM_NPHASE=0.0

      ALLOCATE( CV_SLOC2LOC( CV_SNLOC ))
      ALLOCATE( U_SLOC2LOC( U_SNLOC ))

      ALLOCATE( UGI_COEF_ELE_ALL(NDIM,NPHASE,U_NLOC) )
      ALLOCATE( UGI_COEF_ELE2_ALL(NDIM,NPHASE,U_NLOC) )
      ! The procity mapped to the CV nodes
      ALLOCATE( SUM_CV( CV_NONODS ))

      ALLOCATE( TEN_XX_ONE( NDIM, NDIM, NPHASE )); TEN_XX_ONE=1.0

      ALLOCATE( ONE_PORE( TOTELE ))
      IF ( have_option( '/porous_media/actual_velocity' ) ) THEN
         ! solve for actual velocity
         ONE_PORE = VOLFRA_PORE
      ELSE
         ! solve for porosity * actual velocity
         ONE_PORE = 1.0
      END IF
!       PRINT *,'VOLFRA_PORE:',VOLFRA_PORE
!       PRINT *,'ONE_PORE:',ONE_PORE
!       STOP 22

      D1 = ( NDIM == 1 )
      D3 = ( NDIM == 3 )
      DCYL = ( NDIM == -2 )

      GETMAT = .TRUE.

      X_SHARE = .FALSE.

      ! If using the original limiting scheme, the first step is to estimate
      ! the upwind field value from the surrounding nodes


      !     ======= DEFINE THE SUB-CONTROL VOLUME & FEM SHAPE FUNCTIONS ========

      CALL cv_fem_shape_funs_plus_storage( &
                                ! Volume shape functions...
           NDIM, CV_ELE_TYPE,  &
           CV_NGI, CV_NGI_SHORT, CV_NLOC, U_NLOC, CVN, CVN_SHORT, &
           CVWEIGHT, CVFEN, CVFENLX_ALL, &
           CVWEIGHT_SHORT, CVFEN_SHORT, CVFENLX_SHORT_ALL, &
           UFEN, UFENLX_ALL, &
                                ! Surface of each CV shape functions...
           SCVNGI, CV_NEILOC, CV_ON_FACE, CVFEM_ON_FACE, &
           SCVFEN, SCVFENSLX, SCVFENSLY, SCVFEWEIGH, &
           SCVFENLX_ALL,  &
           SUFEN, SUFENSLX, SUFENSLY,  &
           SUFENLX_ALL,  &
                                ! Surface element shape funcs...
           U_ON_FACE, UFEM_ON_FACE, NFACE, &
           SBCVNGI, SBCVN, SBCVFEN,SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, SBCVFENLX_ALL, &
           SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX_ALL, &
           CV_SLOCLIST, U_SLOCLIST, CV_SNLOC, U_SNLOC, &
                                ! Define the gauss points that lie on the surface of the CV...
           FINDGPTS, COLGPTS, NCOLGPTS, &
           SELE_OVERLAP_SCALE, QUAD_OVER_WHOLE_ELE,&
           state, "cv-adv1" , StorageIndexes(40) )

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
      ! Also determine the CV mass matrix MASS_CV and centre of the CV's XC_CV_ALL
      ! This is for projecting to finite element basis functions...
      ALLOCATE( MASS_CV( CV_NONODS ))
      ALLOCATE( MASS_ELE( TOTELE ))
      ALLOCATE( XC_CV_ALL( NDIM, CV_NONODS ))

      ALLOCATE( DTX_ELE_ALL( NDIM, NPHASE, CV_NLOC, TOTELE ))
      ALLOCATE( DTOLDX_ELE_ALL( NDIM, NPHASE, CV_NLOC, TOTELE ))


      IGETCT = 0
      IF ( GETCT ) IGETCT = 1


! TEMP STUFF HERE

!      ALLOCATE( U_ALL( NDIM, NPHASE, U_NONODS ), NU_ALL( NDIM, NPHASE, U_NONODS ), NUOLD_ALL( NDIM, NPHASE, U_NONODS ) )
!      ALLOCATE( X_ALL( NDIM, X_NONODS ) )
!      ALLOCATE( T_ALL( NPHASE, CV_NONODS ), TOLD_ALL( NPHASE, CV_NONODS ) )

      ALLOCATE( T2_ALL( NPHASE, CV_NONODS ), T2OLD_ALL( NPHASE, CV_NONODS ) )

!      ALLOCATE( FEMT_ALL( NPHASE, CV_NONODS ), FEMTOLD_ALL( NPHASE, CV_NONODS ) )
      ALLOCATE( FEMDEN_ALL( NPHASE, CV_NONODS ), FEMDENOLD_ALL( NPHASE, CV_NONODS ) )
      ALLOCATE( FEMT2_ALL( NPHASE, CV_NONODS ), FEMT2OLD_ALL( NPHASE, CV_NONODS ) )

      ALLOCATE( SOURCT_ALL( NPHASE, CV_NONODS ), ABSORBT_ALL( NPHASE, NPHASE, CV_NONODS ) )


!      NU_ALL=0. ; NUOLD_ALL=0. ; X_ALL=0. ;  X_ALL=0.
!      T_ALL=0. ; TOLD_ALL=0. ; DEN_ALL=0. ; DENOLD_ALL=0. ; T2_ALL=0. ; T2OLD_ALL=0.
!      FEMT_ALL=0. ; FEMTOLD_ALL=0. ; FEMDEN_ALL=0. ; FEMDENOLD_ALL=0. ; FEMT2_ALL=0. ; FEMT2OLD_ALL=0.
!      SOURCT_ALL=0. ; ABSORBT_ALL=0.

!      NU_s => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedNonlinearVelocity" )
!      NU_ALL = NU_s % val
!
!      NUOLD_s => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedOldNonlinearVelocity" )
!      NUOLD_ALL = NUOLD_s % val
!
!      U_ALL = NU_ALL
!      U_s => EXTRACT_TENSOR_FIELD( PACKED_STATE, "PackedVelocity" )
!      if ( for_Sat ) U_ALL = U_s % val

!      ALLOCATE( U( NPHASE * U_NONODS ) ) ; U=0.0
!      ALLOCATE( V( NPHASE * U_NONODS ) ) ; V=0.0
!      ALLOCATE( W( NPHASE * U_NONODS ) ) ; W=0.0

!      ALLOCATE( NU( NPHASE * U_NONODS ), NUOLD( NPHASE * U_NONODS ) ) ; NU=0.0 ; NUOLD=0.0
!      ALLOCATE( NV( NPHASE * U_NONODS ), NVOLD( NPHASE * U_NONODS ) ) ; NV=0.0 ; NVOLD=0.0
!      ALLOCATE( NW( NPHASE * U_NONODS ), NWOLD( NPHASE * U_NONODS ) ) ; NW=0.0 ; NWOLD=0.0

      IF ( IS_OVERLAPPING ) THEN
         NLEV = CV_NLOC
         U_NLOC2 = MAX( 1, U_NLOC / CV_NLOC )
      ELSE
         NLEV = 1
         U_NLOC2 = U_NLOC
      END IF

!            DO U_INOD = 1, U_NONODS
!               DO IPHASE = 1, NPHASE
!                  DO IDIM = 1, NDIM
!                     IF ( IDIM==1 ) THEN
!                        !U_ALL( IDIM, IPHASE, U_INOD ) = U( U_INOD + (IPHASE-1)*U_NONODS )
!                        !NU_ALL( IDIM, IPHASE, U_INOD ) = NU( U_INOD + (IPHASE-1)*U_NONODS )
!                        !NUOLD_ALL( IDIM, IPHASE, U_INOD ) = NUOLD( U_INOD + (IPHASE-1)*U_NONODS )
!                        U( U_INOD + (IPHASE-1)*U_NONODS ) = U_ALL( IDIM, IPHASE, U_INOD )
!                        NU( U_INOD + (IPHASE-1)*U_NONODS ) = NU_ALL( IDIM, IPHASE, U_INOD )
!                        NUOLD( U_INOD + (IPHASE-1)*U_NONODS ) = NUOLD_ALL( IDIM, IPHASE, U_INOD )
!                     ELSE IF ( IDIM==2 ) THEN
!                        !U_ALL( IDIM, IPHASE, U_INOD ) = V( U_INOD + (IPHASE-1)*U_NONODS )
!                        !NU_ALL( IDIM, IPHASE, U_INOD ) = NV( U_INOD + (IPHASE-1)*U_NONODS )
!                        !NUOLD_ALL( IDIM, IPHASE, U_INOD ) = NVOLD( U_INOD + (IPHASE-1)*U_NONODS )
! 
!                        V( U_INOD + (IPHASE-1)*U_NONODS ) = U_ALL( IDIM, IPHASE, U_INOD )
!                        NV( U_INOD + (IPHASE-1)*U_NONODS ) = NU_ALL( IDIM, IPHASE, U_INOD )
!                        NVOLD( U_INOD + (IPHASE-1)*U_NONODS ) = NUOLD_ALL( IDIM, IPHASE, U_INOD )
!                     ELSE
!                        !U_ALL( IDIM, IPHASE, U_INOD ) = W( U_INOD + (IPHASE-1)*U_NONODS )
!                        !NU_ALL( IDIM, IPHASE, U_INOD ) = NW( U_INOD + (IPHASE-1)*U_NONODS )
!                        !NUOLD_ALL( IDIM, IPHASE, U_INOD ) = NWOLD( U_INOD + (IPHASE-1)*U_NONODS )
!                        W( U_INOD + (IPHASE-1)*U_NONODS ) = U_ALL( IDIM, IPHASE, U_INOD )
!                        NW( U_INOD + (IPHASE-1)*U_NONODS ) = NU_ALL( IDIM, IPHASE, U_INOD )
!                        NWOLD( U_INOD + (IPHASE-1)*U_NONODS ) = NUOLD_ALL( IDIM, IPHASE, U_INOD )
!                     END IF
!                  END DO
!               END DO
!            END DO

!print *, u_all(1,1,:)

      DO IPHASE = 1, NPHASE

!         DEN( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS ) = DEN_ALL( IPHASE, : )
!         DENOLD( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS ) = DENOLD_ALL( IPHASE, : )

!         T_ALL( IPHASE, : ) = T( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS )
!         TOLD_ALL( IPHASE, : ) = TOLD( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS )

!         T( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS ) = T_ALL( IPHASE, : )
!         TOLD( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS ) = TOLD_ALL( IPHASE, : )

         IF ( IGOT_T2 == 1 ) THEN
            T2_ALL( IPHASE, : ) = T2( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS )
            T2OLD_ALL( IPHASE, : ) = T2OLD( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS )
         END IF

         SOURCT_ALL( IPHASE, : ) = SOURCT( 1 + (IPHASE-1)*CV_NONODS : IPHASE*CV_NONODS )

         DO JPHASE = 1, NPHASE
            ABSORBT_ALL( JPHASE, IPHASE, : ) = ABSORBT( :, JPHASE, IPHASE )
         END DO

      END DO


!      print *,'before allocate'

! variables for get_int_tden********************
! Set up the fields...

      ALLOCATE( IGOT_T_PACK( NPHASE, 6 ), IGOT_T_CONST( NPHASE, 6 ), IGOT_T_CONST_VALUE( NPHASE, 6 ) )

! FOR packing as well as for detemining which variables to apply interface tracking**********
!          STORE=.TRUE.
          STORE=.FALSE.
          IF( GETCT .AND. RETRIEVE_SOLID_CTY) STORE=.FALSE. ! Avoid storing and retrieving solids voln frac. until we have sorted the code for this. 

          IGOT_T_PACK=.TRUE.
          IGOT_T_CONST      =.FALSE.
          IGOT_T_CONST_VALUE=0.0
! If we have any b.c then assume we ave a non-uniform field... 
          DO IPHASE=1,NPHASE
             IF( SUM(  WIC_T_BC_ALL( :, IPHASE, : ) ) == 0)  &
                          CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,1), IGOT_T_CONST_VALUE(IPHASE,1), T_ALL(IPHASE,:),CV_NONODS)
          END DO
          DO IPHASE=1,NPHASE
             IF( SUM(  WIC_T_BC_ALL( :, IPHASE, : ) ) == 0)  &
                          CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,2), IGOT_T_CONST_VALUE(IPHASE,2), TOLD_ALL(IPHASE,:),CV_NONODS)
          END DO
          DO IPHASE=1,NPHASE
             IF( SUM(  WIC_D_BC_ALL( :, IPHASE, : ) ) == 0)  &
                          CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,3), IGOT_T_CONST_VALUE(IPHASE,3), DEN_ALL(IPHASE,:),CV_NONODS)
          END DO
          DO IPHASE=1,NPHASE
             IF( SUM(  WIC_D_BC_ALL( :, IPHASE, : ) ) == 0)  &
                          CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,4), IGOT_T_CONST_VALUE(IPHASE,4), DENOLD_ALL(IPHASE,:),CV_NONODS)
          END DO

          DO IPHASE=1,NPHASE
             IF(IGOT_T2==1) THEN
                IF( SUM(  WIC_T2_BC_ALL(:,  IPHASE, : ) ) == 0)  &
                          CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,5), IGOT_T_CONST_VALUE(IPHASE,5), T2_ALL(IPHASE,:),CV_NONODS)
             ELSE
                IGOT_T_CONST(IPHASE,5)=.TRUE. 
                IGOT_T_CONST_VALUE(IPHASE,5)=1.0
             ENDIF
          END DO
          DO IPHASE=1,NPHASE
             IF(IGOT_T2==1) THEN
                IF( SUM(  WIC_T2_BC_ALL( :, IPHASE, : ) ) == 0)  &
                          CALL IS_FIELD_CONSTANT(IGOT_T_CONST(IPHASE,6), IGOT_T_CONST_VALUE(IPHASE,6), T2OLD_ALL(IPHASE,:),CV_NONODS)
             ELSE
                IGOT_T_CONST(IPHASE,6)=.TRUE. 
                IGOT_T_CONST_VALUE(IPHASE,6)=1.0
             ENDIF
          END DO

          NFIELD=0
          DO IFI=1,6
             DO IPHASE=1,NPHASE
                IF(.not.IGOT_T_CONST(IPHASE,IFI)) NFIELD=NFIELD+1
             END DO
          END DO

          ALLOCATE( DOWNWIND_EXTRAP_INDIVIDUAL( NFIELD ) )

          

! Determine IGOT_T_PACK(IPHASE,:): 
          IGOT_T_PACK=.FALSE.
          DO IPHASE=1,NPHASE
             DO ILOOP=1,6
                IF(.NOT.IGOT_T_CONST(IPHASE,ILOOP)) THEN
! here we might check to see if we have this in the local storage...
                   IGOT_T_PACK(IPHASE,ILOOP)=.TRUE.
                ENDIF
             END DO
          END DO

          DOWNWIND_EXTRAP_INDIVIDUAL=.FALSE.
          IPT=1
          IF( cv_disopt>=8 ) THEN
             IF(IGOT_T2==1) THEN
                DO IPHASE=1,NPHASE
                   IF(.NOT.IGOT_T_CONST(IPHASE,1)) THEN
                      DOWNWIND_EXTRAP_INDIVIDUAL(IPT)=.TRUE.
                      IPT=IPT+1
                   ENDIF
                END DO
             ENDIF
             IF(IGOT_T2==1) THEN
                DO IPHASE=1,NPHASE
                   IF(.NOT.IGOT_T_CONST(IPHASE,2)) THEN
                      DOWNWIND_EXTRAP_INDIVIDUAL(IPT)=.TRUE.
                      IPT=IPT+1
                   ENDIF
                END DO
             ENDIF
          ENDIF

!         print *,'DOWNWIND_EXTRAP_INDIVIDUAL:',DOWNWIND_EXTRAP_INDIVIDUAL
! FOR packing as well as for detemining which variables to apply interface tracking**********

           IF ( IS_OVERLAPPING .or. is_compact_overlapping ) THEN

               ALLOCATE( VI_LOC_OPT_VEL_UPWIND_COEFS(NDIM,NDIM,NPHASE),  GI_LOC_OPT_VEL_UPWIND_COEFS(NDIM,NDIM,NPHASE),  &
                         VJ_LOC_OPT_VEL_UPWIND_COEFS(NDIM,NDIM,NPHASE),  GJ_LOC_OPT_VEL_UPWIND_COEFS(NDIM,NDIM,NPHASE) )
               IF( is_compact_overlapping ) THEN ! The inverse of the sigma matrix...
                  ALLOCATE( INV_V_OPT_VEL_UPWIND_COEFS(NDIM,NDIM,NPHASE,MAT_NONODS) )
                  ALLOCATE( INV_VI_LOC_OPT_VEL_UPWIND_COEFS(NDIM,NDIM,NPHASE), INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(NDIM,NDIM,NPHASE) )

                  DO MAT_NODI=1,MAT_NONODS
                     DO IPHASE=1,NPHASE
                        DO IDIM=1,NDIM
                           DO JDIM=1,NDIM
                              IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODI-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                              INV_VI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) = OPT_VEL_UPWIND_COEFS(IJ)
                           END DO
                        END DO
                        call invert(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE))
                     END DO
                     INV_V_OPT_VEL_UPWIND_COEFS(:,:,:,MAT_NODI)=INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,:)
                  END DO
! endof IF( is_compact_overlapping ) THEN...
               ENDIF
            else
               ALLOCATE( VI_LOC_OPT_VEL_UPWIND_COEFS(0,0,0),  GI_LOC_OPT_VEL_UPWIND_COEFS(0,0,0),  &
                    VJ_LOC_OPT_VEL_UPWIND_COEFS(0,0,0),  GJ_LOC_OPT_VEL_UPWIND_COEFS(0,0,0) )
               ALLOCATE( INV_V_OPT_VEL_UPWIND_COEFS(0,0,0,0) )
               ALLOCATE( INV_VI_LOC_OPT_VEL_UPWIND_COEFS(0,0,0), INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(0,0,0) )

           ENDIF


           IF( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion equations
               IF(THERMAL) THEN
                     IF( RETRIEVE_SOLID_CTY ) THEN
                        ALLOCATE(VOL_FRA_FLUID(CV_NONODS))
!                        ALLOCATE(U_HAT_ALL(NDIM,U_NONODS))

!                        delta_u_all => extract_vector_field( packed_state, "delta_U" )
!                        u_hat_all = delta_u_all%val + u_all( :, 1, :) ! ndim, u_nonods

!                        us_all => extract_vector_field( packed_state, "solid_U" )

                        Solid_vol_fra => extract_scalar_field( packed_state, "SolidConcentration" )
                        VOL_FRA_FLUID = 1.0 - 1.0 * solid_vol_fra%val   ! cv_nonods
                      endif
               ENDIF
           ENDIF



! F and LOC_U:
      ALLOCATE(LOC_F(NFIELD,CV_NLOC)) 
      ALLOCATE(LOC_FEMF(NFIELD,CV_NLOC)) 

      ALLOCATE(SLOC_F(NFIELD,CV_SNLOC)) 
      ALLOCATE(SLOC_FEMF(NFIELD,CV_SNLOC)) 
      ALLOCATE(SLOC2_F(NFIELD,CV_SNLOC)) 
      ALLOCATE(SLOC2_FEMF(NFIELD,CV_SNLOC)) 

      ALLOCATE(LOC_UF(NDIM,NFIELD,U_NLOC))
      ALLOCATE(SLOC_UF(NDIM,NFIELD,U_SNLOC)) 
      ALLOCATE(SLOC2_UF(NDIM,NFIELD,U_SNLOC))
! Local variables used in GET_INT_VEL_NEW:
      ALLOCATE( LOC_U(NDIM, NPHASE, U_NLOC), LOC2_U(NDIM, NPHASE, U_NLOC) )
      ALLOCATE( LOC_NU(NDIM, NPHASE, U_NLOC), LOC2_NU(NDIM, NPHASE, U_NLOC) )
      ALLOCATE( SLOC_NU(NDIM, NPHASE, U_SNLOC) )
      ALLOCATE( LOC_NUOLD(NDIM, NPHASE, U_NLOC), LOC2_NUOLD(NDIM, NPHASE, U_NLOC) )
      IF(GETCT.AND.RETRIEVE_SOLID_CTY)  ALLOCATE( LOC_U_HAT(NDIM, U_NLOC), LOC2_U_HAT(NDIM, U_NLOC) )
      ALLOCATE( SLOC_NUOLD(NDIM, NPHASE, U_SNLOC) )
      ALLOCATE( LOC_FEMT(NPHASE, CV_NLOC), LOC2_FEMT(NPHASE, CV_NLOC) )
      ALLOCATE( LOC_FEMTOLD(NPHASE, CV_NLOC), LOC2_FEMTOLD(NPHASE, CV_NLOC) )
      ALLOCATE( LOC_FEMT2(NPHASE, CV_NLOC), LOC2_FEMT2(NPHASE, CV_NLOC) )
      ALLOCATE( LOC_FEMT2OLD(NPHASE, CV_NLOC), LOC2_FEMT2OLD(NPHASE, CV_NLOC) )
      ALLOCATE( LOC_T_I( NPHASE ), LOC_T_J( NPHASE ), LOC_TOLD_I( NPHASE ), LOC_TOLD_J( NPHASE ) )
      ALLOCATE( LOC_DEN_I( NPHASE ), LOC_DEN_J( NPHASE ), LOC_DENOLD_I( NPHASE ), LOC_DENOLD_J( NPHASE ) )
      ALLOCATE( LOC_T2_I( NPHASE ), LOC_T2_J( NPHASE ), LOC_T2OLD_I( NPHASE ), LOC_T2OLD_J( NPHASE ) )
! bc's:
      ALLOCATE( SELE_LOC_WIC_F_BC( NFIELD ) )
      ALLOCATE( SLOC_SUF_F_BC(NFIELD, CV_SNLOC) )  
! limiting values...
      ALLOCATE( FUPWIND_IN( NFIELD ) )
      ALLOCATE( FUPWIND_OUT( NFIELD ) )
! limiting and upwinding: 
      ALLOCATE( FVF( NFIELD ), LIMF( NFIELD ) )
      ALLOCATE( F_INCOME( NFIELD ), F_NDOTQ( NFIELD ) )
! 
      ALLOCATE( F_CV_NODI( NFIELD ), F_CV_NODJ( NFIELD ) )
! 
      ALLOCATE( NUOLDGI_ALL( NDIM, NPHASE ), NUGI_ALL( NDIM, NPHASE ),  NU_LEV_GI( NDIM, NPHASE ) )
! 
      IF(THERMAL.AND.GOT_VIS) THEN
         ALLOCATE( VECS_STRESS( NDIM, NDIM, NPHASE, CV_NONODS ), VECS_GRAD_U( NDIM, NDIM, NPHASE, CV_NONODS ) ) ; VECS_STRESS=0. ; VECS_GRAD_U=0.
         ALLOCATE( STRESS_IJ_THERM( NDIM, NDIM, NPHASE ), STRESS_IJ_THERM_J( NDIM, NDIM, NPHASE ) ) ; STRESS_IJ_THERM=0. ; STRESS_IJ_THERM_J=0.
      ENDIF

! NFIELD Variables: 
      ALLOCATE(DIFF_COEF_DIVDX(NPHASE), DIFF_COEFOLD_DIVDX(NPHASE), &
               NDOTQNEW(NPHASE),  NDOTQOLD(NPHASE), INCOMEOLD(NPHASE), LIMT2OLD(NPHASE), LIMDTOLD(NPHASE), &
               NDOTQ(NPHASE), INCOME(NPHASE), LIMT2(NPHASE), LIMTOLD(NPHASE), LIMT(NPHASE), LIMT_HAT(NPHASE), &
               FVTOLD(NPHASE), FVT2OLD(NPHASE), FVDOLD(NPHASE), &
               LIMDOLD(NPHASE), LIMDTT2OLD(NPHASE),&
               FVT(NPHASE), FVT2(NPHASE), FVD(NPHASE), LIMD(NPHASE),  &
               LIMDT(NPHASE), LIMDTT2(NPHASE), SUM_LIMT(NPHASE), SUM_LIMTOLD(NPHASE)  )
      LIMT_HAT=0.0
      ALLOCATE(INCOME_J(NPHASE),INCOMEold_J(NPHASE)) 

      ndotq = 0. ; ndotqold = 0.

! variables for get_int_tden********************


!      print *,'after allocate'

      ! END OF TEMP STUFF HERE

      psi(1)%ptr=>tracer
      psi(2)%ptr=>old_tracer
      FEM_IT=2
      if (.not. is_constant(density)) then
         psi(FEM_IT+1)%ptr=>density
         FEM_IT=FEM_IT+1
      end if
      if (.not. is_constant(old_density)) then
         psi(FEM_IT+1)%ptr=>old_density
         FEM_IT=FEM_IT+1
      end if
      if (present(saturation)) then
         if (.not. is_constant(saturation)) then
            psi(FEM_IT+1)%ptr=>saturation
            FEM_IT=FEM_IT+1
         end if
         if (.not. is_constant(old_saturation)) then
            psi(FEM_IT+1)%ptr=>old_saturation
            FEM_IT=FEM_IT+1
         end if
      end if
      
      do i=1,FEM_IT
         if (has_tensor_field(packed_state,&
              GetFEMName(psi(i)%ptr))) then
            fempsi(i)%ptr=>extract_tensor_field(packed_state,&
                 GetFEMName(psi(i)%ptr))
         else
            allocate(fempsi(i)%ptr)
            call allocate(fempsi(i)%ptr,psi(i)%ptr%mesh,"FEMPSI"//trim(psi(i)%ptr%name),&
                 dim=psi(i)%ptr%dim)
         end if
      end do

      allocate(psi_int(1)%ptr)
      call allocate(psi_int(1)%ptr,1,tracer%mesh,"CV_mass")
      call set(psi_int(1)%ptr,dim=1,val=1.0)

      coord=>extract_vector_field(packed_state,"PressureCoordinate")
      allocate(psi_ave(1)%ptr)
      call allocate(psi_ave(1)%ptr,ndim,tracer%mesh,"Barycentre")
      if (tracer%mesh%continuity<0) then
         psi_ave(1)%ptr%val=x_all(:,x_ndgln)
      else
         call set_all(psi_ave(1)%ptr,X_ALL)
      end if

      call PROJ_CV_TO_FEM_state( packed_state,FEMPSI(1:FEM_IT),&
           PSI(1:FEM_IT), NDIM, &
       PSI_AVE, PSI_INT, MASS_ELE, &
       CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
       CV_NGI_short, CV_NLOC, CVN_short, CVWEIGHT_short,&
       CVFEN_SHORT, CVFENLX_SHORT_ALL(1,:,:), CVFENLX_SHORT_ALL(2,:,:),&
       CVFENLX_SHORT_ALL(3,:,:), &
       X_NONODS, X_ALL, NCOLM, FINDM, COLM, MIDM, &
       IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC)

      XC_CV_ALL=0.0
      XC_CV_ALL(1:NDIM,:)=psi_ave(1)%ptr%val
      MASS_CV=psi_int(1)%ptr%val(1,:)

      FEMT_ALL(:,:)=FEMPSI(1)%ptr%val(1,:,:)
      FEMTOLD_ALL(:,:)=FEMPSI(2)%ptr%val(1,:,:)
      FEM_IT=3
      if (.not. is_constant(density)) then
         FEMDEN_ALL=psi(FEM_IT)%ptr%val(1,:,:)
         FEM_IT=FEM_IT+1
      else
         FEMDEN_ALL=density%val(1,:,:)
      end if
      if (.not. is_constant(old_density)) then
         FEMDENOLD_ALL=psi(FEM_IT)%ptr%val(1,:,:)
         FEM_IT=FEM_IT+1
      else
         FEMDENOLD_ALL=old_density%val(1,:,:)
      end if

      IF ( present(saturation) ) then
         if (.not. is_constant(saturation)) then
            FEMT2_ALL=psi(FEM_IT)%ptr%val(1,:,:)
            FEM_IT=FEM_IT+1
         else
            FEMT2_ALL=saturation%val(1,:,:)
         end if
         if (.not. is_constant(old_saturation)) then
            FEMt2OLD_ALL=psi(FEM_IT)%ptr%val(1,:,:)
            FEM_IT=FEM_IT+1
         else
            FEMt2OLD_ALL=old_saturation%val(1,:,:)
         end if
      end IF
         

      do i=1,FEM_IT-1
         if (fempsi(i)%ptr%name(1:6)=='FEMPSI') then
            call deallocate(fempsi(i)%ptr)
            deallocate(fempsi(i)%ptr)
         end if
      end do

      call deallocate(psi_int(1)%ptr)
      deallocate(psi_int(1)%ptr)
      call deallocate(psi_ave(1)%ptr)
      deallocate(psi_ave(1)%ptr)


       !###FEM VALUES###
! ***********LOOK AT T_FEMT,DEN_FEMT WITH A VIEW TO DELETING EVENTUALLY
!      DO CV_INOD = 1, CV_NONODS
!          DO IPHASE = 1, NPHASE
!              T_FEMT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMT_ALL( IPHASE, CV_INOD)
!              DEN_FEMT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMDEN_ALL( IPHASE, CV_INOD)
!          END DO
!      END DO
!      DO CV_INOD = 1, CV_NONODS
!          DO IPHASE = 1, NPHASE
!              FEMT( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMT_ALL( IPHASE, CV_INOD)
!              FEMTOLD( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMTOLD_ALL( IPHASE, CV_INOD)
!              FEMDEN( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMDEN_ALL( IPHASE, CV_INOD)
!              FEMDENOLD( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMDENOLD_ALL( IPHASE, CV_INOD)
!          END DO
!      END DO
!      if (IGOT_T2>0) then
!          DO CV_INOD = 1, CV_NONODS
!              DO IPHASE = 1, NPHASE
!                  FEMT2( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMT2_ALL( IPHASE, CV_INOD)
!                  FEMT2OLD( CV_INOD + ( IPHASE - 1 ) * CV_NONODS ) = FEMT2OLD_ALL( IPHASE, CV_INOD)
!              END DO
!          END DO
!      end if

!      XC_CV_ALL(1,:)=XC_CV(:)
!      IF(NDIM.GE.2) XC_CV_ALL(2,:)=YC_CV(:)
!      IF(NDIM.GE.3) XC_CV_ALL(3,:)=ZC_CV(:)




      MASS_ELE_TRANSP = MASS_ELE

      NORMALISE = .FALSE.
      IF ( NORMALISE ) THEN
         ! Make sure the FEM representation sums to unity so we don't get surprising results...
         DO CV_INOD = 1, CV_NONODS
            RSUM = SUM( FEMT_ALL( :, CV_INOD ) )
            FEMT_ALL( :, CV_INOD ) = FEMT_ALL( :, CV_INOD ) / RSUM
         END DO
      END IF

      ! Calculate MEAN_PORE_CV
      MEAN_PORE_CV = 0.0 ; SUM_CV = 0.0 
      DO ELE = 1, TOTELE
         DO CV_ILOC = 1, CV_NLOC
            CV_INOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
            SUM_CV( CV_INOD ) = SUM_CV( CV_INOD ) + MASS_ELE( ELE )
            MEAN_PORE_CV( CV_INOD ) = MEAN_PORE_CV( CV_INOD ) + &
                 MASS_ELE( ELE ) * VOLFRA_PORE( ELE )
         END DO
      END DO
      MEAN_PORE_CV = MEAN_PORE_CV / SUM_CV
      ewrite(3,*) 'MEAN_PORE_CV MIN/MAX:', MINVAL( MEAN_PORE_CV ), MAXVAL( MEAN_PORE_CV )

      IANISOLIM = 0
      IF ( CV_DISOPT >= 5 ) IANISOLIM = 1

      NSMALL_COLM = SIZE( SMALL_COLM )

      ALLOCATE( TUPWIND_MAT_ALL( NPHASE, NSMALL_COLM ), TOLDUPWIND_MAT_ALL( NPHASE, NSMALL_COLM ), &
              DENUPWIND_MAT_ALL( NPHASE, NSMALL_COLM ), DENOLDUPWIND_MAT_ALL( NPHASE, NSMALL_COLM ) )
      ALLOCATE( T2UPWIND_MAT_ALL( NPHASE*IGOT_T2, NSMALL_COLM* IGOT_T2), T2OLDUPWIND_MAT_ALL( NPHASE*IGOT_T2, NSMALL_COLM*IGOT_T2 ) )

    !#############CONVERSION FROM NEW VARIABLES TO OLD VARIABLES############
!         ALLOCATE( TUPWIND_MAT( NSMALL_COLM*NPHASE ), TOLDUPWIND_MAT( NSMALL_COLM*NPHASE ), &
!              DENUPWIND_MAT( NSMALL_COLM*NPHASE ), DENOLDUPWIND_MAT( NSMALL_COLM*NPHASE ) )
!         ALLOCATE( T2UPWIND_MAT( NSMALL_COLM*NPHASE*IGOT_T2 ), T2OLDUPWIND_MAT( NSMALL_COLM*NPHASE*IGOT_T2 ) )


      IF ( IANISOLIM == 0 ) THEN

! Isotropic limiting - calculate far field upwind maticies...
        CALL ISOTROPIC_LIMITER_ALL( &
! FOR SUB SURRO_CV_MINMAX:
           T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, size(small_colm), SMALL_CENTRM, SMALL_FINDRM, SMALL_COLM, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
           MASS_CV, &
! FOR SUB CALC_LIMIT_MATRIX_MAX_MIN:
           TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
           TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL )

      ELSE ! endof IF ( IANISOLIM == 0 ) THEN


         CALL CALC_ANISOTROP_LIM( &
              ! Caculate the upwind values stored in matrix form...
              T_ALL,TOLD_ALL,DEN_ALL,DENOLD_ALL,T2_ALL,T2OLD_ALL, &
              FEMT_ALL,FEMTOLD_ALL,FEMDEN_ALL,FEMDENOLD_ALL,FEMT2_ALL,FEMT2OLD_ALL, (CV_NONODS.NE.X_NONODS), &
              TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, &
              T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
              ! Store the upwind element for interpolation and its weights for
              ! faster results...
              IGOT_T2,NPHASE,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
              SMALL_FINDRM,SMALL_CENTRM,SMALL_COLM,NSMALL_COLM, &
              X_NDGLN,X_NONODS,NDIM, &
              X_ALL, XC_CV_ALL, IGOT_T_PACK, IGOT_T_CONST, IGOT_T_CONST_VALUE,&
              state, "anisotrop", storageindexes(42))

         
      END IF ! endof IF ( IANISOLIM == 0 ) THEN ELSE


!      DO COUNT = 1, NSMALL_COLM
!          DO IPHASE = 1, NPHASE
!              TUPWIND_MAT(COUNT + ( IPHASE - 1 ) * NSMALL_COLM) = TUPWIND_MAT_ALL(IPHASE, COUNT)
!              TOLDUPWIND_MAT(COUNT + ( IPHASE - 1 ) * NSMALL_COLM) = TOLDUPWIND_MAT_ALL(IPHASE, COUNT)
!              DENUPWIND_MAT(COUNT + ( IPHASE - 1 ) * NSMALL_COLM) = DENUPWIND_MAT_ALL(IPHASE, COUNT)
!              DENOLDUPWIND_MAT(COUNT + ( IPHASE - 1 ) * NSMALL_COLM) = DENOLDUPWIND_MAT_ALL(IPHASE, COUNT)
!              IF ( IGOT_T2 == 1 ) THEN
!                  T2UPWIND_MAT(COUNT + ( IPHASE - 1 ) * NSMALL_COLM) = T2UPWIND_MAT_ALL(IPHASE, COUNT)
!                  T2OLDUPWIND_MAT(COUNT + ( IPHASE - 1 ) * NSMALL_COLM) = T2OLDUPWIND_MAT_ALL(IPHASE, COUNT)
!              end if
!          end do
!      end do

      ALLOCATE( FACE_ELE( NFACE, TOTELE ) ) ; FACE_ELE = 0
      CALL CALC_FACE_ELE( FACE_ELE, TOTELE, STOTEL, NFACE, &
           NCOLELE, FINELE, COLELE, CV_NLOC, CV_SNLOC, CV_NONODS, CV_NDGLN, CV_SNDGLN, &
           CV_SLOCLIST, X_NLOC, X_NDGLN )

      IF ( GOT_DIFFUS ) THEN
         CALL DG_DERIVS_ALL2( FEMT_ALL, FEMTOLD_ALL, &
              DTX_ELE_ALL, DTOLDX_ELE_ALL, &
              NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
              X_NDGLN, X_NLOC, X_NDGLN,&
              CV_NGI, CV_NLOC, CVWEIGHT_SHORT, &
              CVFEN_SHORT, CVFENLX_SHORT_ALL(1,:,:), CVFENLX_SHORT_ALL(2,:,:), CVFENLX_SHORT_ALL(3,:,:), &
              CVFEN_SHORT, CVFENLX_SHORT_ALL(1,:,:), CVFENLX_SHORT_ALL(2,:,:), CVFENLX_SHORT_ALL(3,:,:), &
              X_NONODS, X_ALL(1,:),X_ALL(2,:),X_ALL(3,:), &
              NFACE, FACE_ELE, CV_SLOCLIST, CV_SLOCLIST, STOTEL, CV_SNLOC, CV_SNLOC, WIC_T_BC_ALL, SUF_T_BC_ALL, &
              SBCVNGI, SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFEWEIGH, &
              SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
              state, "CVN", StorageIndexes( 33 ) )

      END IF



      !     =============== DEFINE THETA FOR TIME-STEPPING ===================
      ! Define the type of time integration:
      ! Timopt is 0 if CV_DISOPT is even (theta specified);
      ! Timopt is 1 if CV_DISOPT is odd (non-linear theta scheme)
      TIMOPT = MOD( CV_DISOPT, 2 )

      VTHETA = 1.0
      FTHETA = CV_THETA

      IF ( GETCT ) THEN ! Obtain the CV discretised CT eqns plus RHS
         CT_RHS = 0.0
         CT = 0.0
      END IF

      IF ( GETCV_DISC ) THEN ! Obtain the CV discretised advection/diffusion eqns
         CV_RHS = 0.0
         IF ( GETMAT ) THEN
            CSR_ACV = 0.0
         END IF
      END IF


      GET_GTHETA = .FALSE.
      IF ( IGOT_THETA_FLUX == 1 ) THEN
         IF ( GET_THETA_FLUX ) THEN
            THETA_FLUX = 0.0
            ONE_M_THETA_FLUX = 0.0
            if(integrate_other_side) then
               THETA_FLUX_J = 0.0
               ONE_M_THETA_FLUX_J = 0.0
            endif
            IF ( IGOT_T2 == 1 ) THEN
               GET_GTHETA = .TRUE.
               THETA_GDIFF = 0.0
            END IF
         END IF
      END IF



      GLOBAL_FACE = 0

      Loop_Elements: DO ELE = 1, TOTELE

         if (IsParallel()) then
            if (.not. assemble_ele(tracer,ele)) then
               skip=.true.
               neighbours=>ele_neigh(tracer,ele)
               do nb=1,size(neighbours)
                  if (neighbours(nb)<=0) cycle 
                  if (assemble_ele(tracer,neighbours(nb))) then
                     skip=.false.
                     exit
                  end if
               end do
               if (skip) cycle
            end if
         end if
                  

         ! Calculate DETWEI, RA, NX, NY, NZ for element ELE
         CALL DETNLXR_INVJAC( ELE, X_ALL, X_NDGLN, TOTELE, X_NONODS, &
              CV_NLOC, SCVNGI, &
              SCVFEN, SCVFENLX_ALL, SCVFEWEIGH, SCVDETWEI, SCVRA, VOLUME, DCYL, &
              SCVFENX_ALL, &
              NDIM, INV_JAC, state, "CVI", StorageIndexes(29) )


! Generate some local F variables ***************

! loc_f
          DO CV_KLOC = 1, CV_NLOC
             CV_NODK = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_KLOC ) 
             IPT=1
             CALL PACK_LOC( LOC_F(:, CV_KLOC), T_ALL( :, CV_NODK ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) )
             CALL PACK_LOC( LOC_F(:, CV_KLOC), TOLD_ALL( :, CV_NODK ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) )
             CALL PACK_LOC( LOC_F(:, CV_KLOC), DEN_ALL( :, CV_NODK ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )
             CALL PACK_LOC( LOC_F(:, CV_KLOC), DENOLD_ALL( :, CV_NODK ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )
             CALL PACK_LOC( LOC_F(:, CV_KLOC), T2_ALL( :, CV_NODK ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )
             CALL PACK_LOC( LOC_F(:, CV_KLOC), T2OLD_ALL( :, CV_NODK ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )
! for FEM variables...
             IPT=1
             CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMT_ALL( :, CV_NODK ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) )
             CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMTOLD_ALL( :, CV_NODK ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) )
             CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMDEN_ALL( :, CV_NODK ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )
             CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMDENOLD_ALL( :, CV_NODK ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )
             CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMT2_ALL( :, CV_NODK ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )
             CALL PACK_LOC( LOC_FEMF(:, CV_KLOC), FEMT2OLD_ALL( :, CV_NODK ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )
          END DO
    
! LOC_UF
          DO U_KLOC = 1, U_NLOC
             U_NODK = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_KLOC )

             IPT=1
             IFI=1
             DO ILOOP=1,3
                DO IPHASE=1,NPHASE
                   IF(IGOT_T_PACK(IPHASE,IFI)) THEN ! T: 
                      LOC_UF(:, IPT, U_KLOC) =   U_ALL( :, IPHASE, U_NODK )
                      IPT=IPT+1
                   END IF 
                END DO
                IFI=IFI+1

                DO IPHASE=1,NPHASE
                   IF(IGOT_T_PACK(IPHASE,IFI)) THEN ! Told: 
         if(correct_method_petrov_method) then
                      LOC_UF(:, IPT, U_KLOC) =   NUOLD_ALL( :, IPHASE, U_NODK )
         else
                      LOC_UF(:, IPT, U_KLOC) =   U_ALL( :, IPHASE, U_NODK )
         endif
                      IPT=IPT+1
                   END IF 
                END DO
                IFI=IFI+1

             END DO

          END DO

! LOC_U, LOC_NU: 
         IF ( IS_OVERLAPPING ) THEN
            NLEV = CV_NLOC
            U_NLOC2 = MAX( 1, U_NLOC / CV_NLOC )
         ELSE
            NLEV = 1
            U_NLOC2 = U_NLOC
         END IF

         DO ILEV = 1, NLEV
            DO U_ILOC = 1 + (ILEV-1)*U_NLOC2, ILEV*U_NLOC2
               U_INOD = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_ILOC )

               LOC_U( :, :, U_ILOC)=U_ALL( :, :, U_INOD)
               LOC_NU( :, :, U_ILOC)=NU_ALL( :, :, U_INOD)
               LOC_NUOLD( :, :, U_ILOC)=NUOLD_ALL( :, :, U_INOD)
               IF(GETCT.AND.RETRIEVE_SOLID_CTY) LOC_U_HAT( :, U_ILOC)=U_HAT_ALL( :, U_INOD)
            END DO
         END DO 

         DO CV_KLOC=1,CV_NLOC
            CV_KNOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
 
            LOC_FEMT(:, CV_KLOC) = FEMT_ALL(:, CV_KNOD)
            LOC_FEMTOLD(:, CV_KLOC) = FEMTOLD_ALL(:, CV_KNOD)

            IF ( IGOT_T2 == 1 ) THEN
               LOC_FEMT2(:, CV_KLOC) = FEMT2_ALL(:, CV_KNOD)
               LOC_FEMT2OLD(:, CV_KLOC) = FEMT2OLD_ALL(:, CV_KNOD)
            END IF

          END DO

! Generate some local F variables ***************...
! 
! 


         Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC ! Loop over the nodes of the element

            ! Global node number of the local node
            CV_NODI = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
            X_NODI = X_NDGLN( ( ELE - 1 ) * X_NLOC  + CV_ILOC )
            MAT_NODI = MAT_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )
            IMID = SMALL_CENTRM(CV_NODI)

! Generate some local F variables ***************
            F_CV_NODI(:)= LOC_F(:, CV_ILOC)
            IF ( IS_OVERLAPPING .or. is_compact_overlapping ) THEN

               DO IPHASE=1,NPHASE
                  DO IDIM=1,NDIM
                     DO JDIM=1,NDIM
                        IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODI-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                        VI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) = OPT_VEL_UPWIND_COEFS(IJ) 
                        GI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) = OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM) 
                     END DO
                  END DO
               END DO
               IF( is_compact_overlapping ) THEN ! The inverse of the sigma matrix...
                  INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,:) = INV_V_OPT_VEL_UPWIND_COEFS(:,:,:,MAT_NODI)
               ENDIF

            ENDIF
! Generate some local F variables ***************

            ! Loop over quadrature (gauss) points in ELE neighbouring ILOC
            Loop_GCOUNT: DO GCOUNT = FINDGPTS( CV_ILOC ), FINDGPTS( CV_ILOC + 1 ) - 1
               ! COLGPTS stores the local Gauss-point number in the ELE
               GI = COLGPTS( GCOUNT )

               ! Get the neighbouring node for node ILOC and Gauss point GI
               CV_JLOC = CV_NEILOC( CV_ILOC, GI )

               ELE2 = 0
               SELE = 0
               CV_SILOC=0
               INTEGRAT_AT_GI = .TRUE.

               Conditional_CheckingNeighbourhood: IF ( CV_JLOC == -1 ) THEN

                  ! We are on the boundary or next to another element.  Determine CV_OTHER_LOC
                  ! CVFEM_ON_FACE(CV_KLOC,GI)=.TRUE. if CV_KLOC is on the face that GI is centred on.
                  ! Look for these nodes on the other elements.
                  CALL FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC, &
                       MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI, &
                       X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
                       CV_SNLOC, CVFEM_ON_FACE( :, GI ), X_SHARE, X_NONODS, ELE, ELE2,  &
                       FINELE, COLELE, NCOLELE, DISTCONTINUOUS_METHOD )

                  IF ( INTEGRAT_AT_GI ) THEN
                     CV_JLOC = CV_OTHER_LOC( CV_ILOC )
                     SELE = 0
                     ELE3=0

                     IF ( CV_JLOC == 0 ) THEN ! We are on the boundary of the domain or subdomain
                        CV_JLOC = CV_ILOC
                        ! Calculate SELE, CV_SILOC, U_SLOC2LOC, CV_SLOC2LOC
                        CALL CALC_SELE( ELE, ELE3, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
                             FACE_ELE, NFACE, CVFEM_ON_FACE( :, GI ), &
                             CV_NONODS, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
                             CV_NDGLN, U_NDGLN, CV_SNDGLN, U_SNDGLN )
                     END IF

                     INTEGRAT_AT_GI = .NOT.( (ELE==ELE2) .AND. (SELE==0) )
!                     INTEGRAT_AT_GI = .NOT.( (ELE3==0) .AND. (SELE==0) )
!                      if(
                  END IF

                  IF(INTEGRAT_AT_GI)  THEN 
! this is for DG and boundaries of the domain
                     IF(SELE.LE.0) THEN ! this is for DG
! Calculate U_SLOC2LOC, CV_SLOC2LOC: 
                        CV_SKLOC=0
                        DO CV_KLOC=1,CV_NLOC
                           CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                           IF(CV_KLOC2.NE.0) THEN
                              CV_SKLOC=CV_SKLOC+1
                              CV_SLOC2LOC(CV_SKLOC)=CV_KLOC
                              SHAPE_CV_SNL(CV_SKLOC) = SCVFEN(CV_KLOC,GI) 
                           ENDIF
                        END DO
                     ENDIF ! ENDOF IF(SELE.LE.0) THEN

                     DO CV_SKLOC=1,CV_SNLOC
                        CV_KLOC = CV_SLOC2LOC(CV_SKLOC)
                        SHAPE_CV_SNL(CV_SKLOC) = SCVFEN(CV_KLOC,GI) 
                     END DO

                  ENDIF

               END IF Conditional_CheckingNeighbourhood
               ! Avoid integrating across the middle of a CV on the boundaries of elements
               Conditional_integration: IF ( INTEGRAT_AT_GI ) THEN


                  ! Find its global node number
                  IF ( ELE2 == 0 ) THEN
                     CV_NODJ = CV_NDGLN( ( ELE - 1 )  * CV_NLOC + CV_JLOC )
                     MAT_NODJ = MAT_NDGLN( ( ELE - 1 )  * CV_NLOC + CV_JLOC )
                  ELSE
                     CV_NODJ = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_JLOC )
                     MAT_NODJ = MAT_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_JLOC )
                  END IF

           if((CV_NODJ.ge.CV_NODI).or.(.not.integrate_other_side)) then 

                  integrate_other_side_and_not_boundary = integrate_other_side.and.(SELE.LE.0)

                  GLOBAL_FACE = GLOBAL_FACE + 1
                  JMID = SMALL_CENTRM(CV_NODJ)


                  ! Calculate the control volume normals at the Gauss pts.
                  CALL SCVDETNX_new( ELE, GI, &
                       X_NLOC, SCVNGI, TOTELE, NDIM, &
                       X_NDGLN, X_NONODS, &
                       SCVDETWEI, CVNORMX, CVNORMY, &
                       CVNORMZ, SCVFEN, SCVFENSLX, &
                       SCVFENSLY, SCVFEWEIGH, XC_CV_ALL( 1:NDIM, CV_NODI ), &
                       X_ALL(1:NDIM,:),  &
                       D1, D3, DCYL )

                  CVNORMX_ALL(1,GI)=CVNORMX(GI)
                  IF(NDIM.GE.2) CVNORMX_ALL(2,GI)=CVNORMY(GI)
                  IF(NDIM.GE.3) CVNORMX_ALL(3,GI)=CVNORMZ(GI)
               !   CVNORMX(GI) = CVNORMX_ALL(1,GI)
               !   IF(NDIM.GE.2) CVNORMY(GI) =CVNORMX_ALL(2,GI)
               !   IF(NDIM.GE.3) CVNORMZ(GI) =CVNORMX_ALL(3,GI)

                  
                  IF( GETCT ) THEN
! could retrieve JCOUNT_KLOC and ICOUNT_KLOC from storage depending on quadrature point GLOBAL_FACE
                     DO U_KLOC = 1, U_NLOC
                        U_NODK = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_KLOC )
                        JCOUNT = 0
                        DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1
                           IF ( COLCT( COUNT ) == U_NODK ) THEN
                              JCOUNT = COUNT
                              EXIT
                           END IF
                        END DO
                        JCOUNT_KLOC( U_KLOC ) = JCOUNT
                if(integrate_other_side) then
! for integrating just on one side...
                        ICOUNT = 0
                        DO COUNT = FINDCT( CV_NODJ ), FINDCT( CV_NODJ + 1 ) - 1
                           IF ( COLCT( COUNT ) == U_NODK ) THEN
                              ICOUNT = COUNT
                              EXIT
                           END IF
                        END DO
                        ICOUNT_KLOC( U_KLOC ) = ICOUNT
                endif
                     END DO
                     IF ( ( ELE2 /= 0 ) .AND. ( ELE2 /= ELE ) ) THEN
                        DO U_KLOC =  1, U_NLOC
                           U_NODK = U_NDGLN( ( ELE2 - 1 ) * U_NLOC + U_KLOC )
                           JCOUNT = 0
                           DO COUNT = FINDCT( CV_NODI ), FINDCT( CV_NODI + 1 ) - 1
                              IF ( COLCT( COUNT ) == U_NODK ) THEN
                                 JCOUNT = COUNT
                                 EXIT
                              END IF
                           END DO
                           JCOUNT_KLOC2( U_KLOC ) = JCOUNT
                if(integrate_other_side) then
! for integrating just on one side...
                           ICOUNT = 0
                           DO COUNT = FINDCT( CV_NODJ ), FINDCT( CV_NODJ + 1 ) - 1
                              IF ( COLCT( COUNT ) == U_NODK ) THEN
                                 ICOUNT = COUNT
                                 EXIT
                              END IF
                           END DO
                           ICOUNT_KLOC2( U_KLOC ) = ICOUNT
                endif
                        END DO
                     END IF ! endof IF ( ( ELE2 /= 0 ) .AND. ( ELE2 /= ELE ) ) THEN
                  END IF ! endof IF( GETCT ) THEN

                  ! Compute the distance HDC between the nodes either side of the CV face
                  ! (this is needed to compute the local courant number and the non-linear theta)
                  
                  IF ( SELE == 0 ) THEN
                     HDC = SQRT( SUM( (XC_CV_ALL(1:NDIM,CV_NODI)-XC_CV_ALL(1:NDIM,CV_NODJ))**2) )
                  ELSE
                     HDC = SQRT( SUM( (XC_CV_ALL(1:NDIM,CV_NODI)-X_ALL(1:NDIM,X_NODI))**2) )
                  END IF


                  DO COUNT = SMALL_FINDRM( CV_NODI ), SMALL_FINDRM( CV_NODI + 1 ) - 1
                     IF ( SMALL_COLM( COUNT ) == CV_NODJ ) THEN 
                        JCOUNT_IPHA = COUNT
                        EXIT
                     END IF
                  END DO

                  DO COUNT = SMALL_FINDRM( CV_NODJ ), SMALL_FINDRM( CV_NODJ + 1 ) - 1
                     IF ( SMALL_COLM( COUNT ) == CV_NODI ) THEN 
                        ICOUNT_IPHA = COUNT
                        EXIT
                     END IF
                  END DO
 
!         if(.true.) then
! Generate some local F variables ***************
             IPT=1
             CALL PACK_LOC( F_CV_NODJ(:), T_ALL( :, CV_NODJ ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) )
             CALL PACK_LOC( F_CV_NODJ(:), TOLD_ALL( :, CV_NODJ ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) )
             CALL PACK_LOC( F_CV_NODJ(:), DEN_ALL( :, CV_NODJ ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )
             CALL PACK_LOC( F_CV_NODJ(:), DENOLD_ALL( :, CV_NODJ ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )
             CALL PACK_LOC( F_CV_NODJ(:), T2_ALL( :, CV_NODJ ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )
             CALL PACK_LOC( F_CV_NODJ(:), T2OLD_ALL( :, CV_NODJ ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )
! Generate some local F variables ***************

! local surface information***********

          IF( (ELE2 > 0) .OR. (SELE > 0) ) THEN
!          IF( (ELE2 > 0) .and. (SELE > 0) ) THEN
             DO CV_SKLOC = 1, CV_SNLOC
                CV_KLOC = CV_SLOC2LOC( CV_SKLOC )

                   SLOC_F(:, CV_SKLOC) = LOC_F(:, CV_KLOC) 
                   SLOC_FEMF(:, CV_SKLOC) = LOC_FEMF(:, CV_KLOC) 

                
                   IF(ELE2>0) THEN
                      CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
                      CV_NODK2 = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_KLOC2 ) 
             IPT=1
             CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), T_ALL( :, CV_NODK2 ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) )
             CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), TOLD_ALL( :, CV_NODK2 ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) )
             CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), DEN_ALL( :, CV_NODK2 ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )
             CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), DENOLD_ALL( :, CV_NODK2 ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )
             CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), T2_ALL( :, CV_NODK2 ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )
             CALL PACK_LOC( SLOC2_F(:, CV_SKLOC), T2OLD_ALL( :, CV_NODK2 ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )
! femf:
             IPT=1
             CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMT_ALL( :, CV_NODK2 ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) )
             CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMTOLD_ALL( :, CV_NODK2 ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) )
             CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMDEN_ALL( :, CV_NODK2 ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )
             CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMDENOLD_ALL( :, CV_NODK2 ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )
             CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMT2_ALL( :, CV_NODK2 ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )
             CALL PACK_LOC( SLOC2_FEMF(:, CV_SKLOC), FEMT2OLD_ALL( :, CV_NODK2 ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )

                   ELSE
                      SLOC2_F(:, CV_SKLOC)    = SLOC_F(:, CV_SKLOC)
                      SLOC2_FEMF(:, CV_SKLOC) = SLOC_FEMF(:, CV_SKLOC)
                   ENDIF
             END DO


          ENDIF ! ENDOF IF( (ELE2 > 0) .OR. (SELE > 0) ) THEN ELSE...


!         if(.true.) then
          IF( SELE > 0 ) THEN
! bcs:              
             ! Make allowances for no matrix stencil operating from outside the boundary.
             BCZERO=1.0-INCOME
! What type of b.c's -integer
             IPT=1
             CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),&
                  WIC_T_BC_ALL( : , : , SELE ),&
                  NPHASE, NFIELD, IPT, IGOT_T_PACK( :,1) )
             CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),&
                  WIC_T_BC_ALL( : , :, SELE ),&
                  NPHASE, NFIELD, IPT, IGOT_T_PACK( :,2) )
             CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),&
                  WIC_D_BC_ALL( :,:, SELE ),&
                  NPHASE, NFIELD, IPT, IGOT_T_PACK( :,3) )
             CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),&
                  WIC_D_BC_ALL( :,:, SELE ),&
                  NPHASE, NFIELD, IPT, IGOT_T_PACK( :,4) )

             IF(IGOT_T2==1) THEN
                CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),&
                     WIC_T2_BC_ALL( : , :, SELE ),&
                     NPHASE, NFIELD, IPT, IGOT_T_PACK( :,5) )
                CALL I_PACK_LOC( SELE_LOC_WIC_F_BC( : ),&
                     WIC_T2_BC_ALL( : , : , SELE ),&
                     NPHASE, NFIELD, IPT, IGOT_T_PACK( :,6) )
             ENDIF
! The b.c values: 
             DO CV_SKLOC=1,CV_SNLOC
                IPT=1
                CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T_BC_ALL( 1, :, CV_SKLOC + CV_SNLOC*( SELE- 1) ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) )
                CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T_BC_ALL( 1, :, CV_SKLOC+ CV_SNLOC*( SELE- 1 ) ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) )

                CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_D_BC_ALL( 1, :, CV_SKLOC+ CV_SNLOC*( SELE- 1) ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )
                CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_D_BC_ALL( 1, :, CV_SKLOC+ CV_SNLOC*( SELE- 1) ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )

               IF(IGOT_T2==1) THEN
                   CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T2_BC_ALL( 1, :, CV_SKLOC + CV_SNLOC*( SELE- 1) ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )
                   CALL PACK_LOC( SLOC_SUF_F_BC( :, CV_SKLOC ), SUF_T2_BC_ALL( 1, :, CV_SKLOC + CV_SNLOC*( SELE- 1) ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )
               ENDIF
            END DO

         ENDIF ! IF( SELE > 0 ) THEN
! local surface information***********

! limiting VALUES*************: 
          COUNT_OUT=JCOUNT_IPHA
          COUNT_IN =ICOUNT_IPHA
          
          IPT_IN =1
          IPT_OUT=1
          CALL PACK_LOC( FUPWIND_IN( : ),  TUPWIND_MAT_ALL( :, COUNT_IN),    NPHASE, NFIELD, IPT_IN, IGOT_T_PACK(:,1) )
          CALL PACK_LOC( FUPWIND_OUT( : ), TUPWIND_MAT_ALL( :, COUNT_OUT),    NPHASE, NFIELD, IPT_OUT, IGOT_T_PACK(:,1) )
          CALL PACK_LOC( FUPWIND_IN( : ),  TOLDUPWIND_MAT_ALL( :, COUNT_IN),    NPHASE, NFIELD, IPT_IN, IGOT_T_PACK(:,2) )
          CALL PACK_LOC( FUPWIND_OUT( : ), TOLDUPWIND_MAT_ALL( :, COUNT_OUT),    NPHASE, NFIELD, IPT_OUT, IGOT_T_PACK(:,2) )

          CALL PACK_LOC( FUPWIND_IN( : ),  DENUPWIND_MAT_ALL( :, COUNT_IN),    NPHASE, NFIELD, IPT_IN, IGOT_T_PACK(:,3) )
          CALL PACK_LOC( FUPWIND_OUT( : ), DENUPWIND_MAT_ALL( :, COUNT_OUT),    NPHASE, NFIELD, IPT_OUT, IGOT_T_PACK(:,3) )
          CALL PACK_LOC( FUPWIND_IN( : ),  DENOLDUPWIND_MAT_ALL( :, COUNT_IN),    NPHASE, NFIELD, IPT_IN, IGOT_T_PACK(:,4) )
          CALL PACK_LOC( FUPWIND_OUT( : ), DENOLDUPWIND_MAT_ALL( :, COUNT_OUT),    NPHASE, NFIELD, IPT_OUT, IGOT_T_PACK(:,4) )

          IF(IGOT_T2==1) THEN
             CALL PACK_LOC( FUPWIND_IN( : ),  T2UPWIND_MAT_ALL( :, COUNT_IN),    NPHASE, NFIELD, IPT_IN, IGOT_T_PACK(:,5) )
             CALL PACK_LOC( FUPWIND_OUT( : ), T2UPWIND_MAT_ALL( :, COUNT_OUT),    NPHASE, NFIELD, IPT_OUT, IGOT_T_PACK(:,5) )
             CALL PACK_LOC( FUPWIND_IN( : ),  T2OLDUPWIND_MAT_ALL( :, COUNT_IN),    NPHASE, NFIELD, IPT_IN, IGOT_T_PACK(:,6) )
             CALL PACK_LOC( FUPWIND_OUT( : ), T2OLDUPWIND_MAT_ALL( :, COUNT_OUT),    NPHASE, NFIELD, IPT_OUT, IGOT_T_PACK(:,6) )
          ENDIF
! limiting VALUES*************:
!     endif ! endof if(.false.) then
!         

! LOC2_U, LOC2_NU for GET_INT_VEL_NEW
       IF (ELE2/=0) THEN
          DO U_KLOC = 1, U_NLOC
             U_KLOC2 = U_OTHER_LOC( U_KLOC )
             IF ( U_KLOC2 /= 0 ) THEN
                U_NODK2 = U_NDGLN((ELE2-1)*U_NLOC+U_KLOC2)

                LOC2_U(:, :, U_KLOC) = U_ALL(:, :, U_NODK2)
                LOC2_NU(:, :, U_KLOC) = NU_ALL(:, :, U_NODK2)
                LOC2_NUOLD(:, :, U_KLOC) = NUOLD_ALL(:, :, U_NODK2)
                IF(GETCT.AND.RETRIEVE_SOLID_CTY) LOC2_U_HAT(:, U_KLOC) = U_HAT_ALL(:, U_NODK2)
             END IF
          END DO

         DO CV_KLOC=1,CV_NLOC
            CV_KNOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_KLOC)
 
            CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
            IF (CV_KLOC2 /= 0 ) THEN
               CV_KNOD2 = CV_NDGLN((ELE2-1)*CV_NLOC+CV_KLOC2)
               LOC2_FEMT(:, CV_KLOC) = FEMT_ALL(:, CV_KNOD2)
               LOC2_FEMTOLD(:, CV_KLOC) = FEMTOLD_ALL(:, CV_KNOD2)
               IF ( IGOT_T2 == 1 ) THEN
                  LOC2_FEMT2(:, CV_KLOC) = FEMT2_ALL(:, CV_KNOD2)
                  LOC2_FEMT2OLD(:, CV_KLOC) = FEMT2OLD_ALL(:, CV_KNOD2)
               END IF
            END IF
         END DO

       END IF

       IF ( SELE /= 0 ) THEN
          DO U_SKLOC = 1, U_SNLOC
             U_KLOC = U_SLOC2LOC( U_SKLOC )
             U_NODK = U_NDGLN(( ELE - 1 ) * U_NLOC + U_KLOC )
             U_SNODK = ( SELE - 1 ) * U_SNLOC + U_SKLOC 
         
             SLOC_NU(:, :, U_SKLOC) = NU_ALL(:, :, U_NODK)
             SLOC_NUOLD(:, :, U_SKLOC) = NUOLD_ALL(:, :, U_NODK)
    
          END DO
       END IF

       LOC_T_I( : ) = T_ALL(:, CV_NODI)
       LOC_T_J( : ) = T_ALL(:, CV_NODJ)
       LOC_TOLD_I( : ) = TOLD_ALL(:, CV_NODI)
       LOC_TOLD_J( : ) = TOLD_ALL(:, CV_NODJ)
       LOC_DEN_I( : ) = DEN_ALL(:, CV_NODI)
       LOC_DEN_J( : ) = DEN_ALL(:, CV_NODJ)
       LOC_DENOLD_I( : ) = DENOLD_ALL(:, CV_NODI)
       LOC_DENOLD_J( : ) = DENOLD_ALL(:, CV_NODJ)
       IF ( IGOT_T2 == 1 ) THEN
          LOC_T2_I( : ) = T2_ALL(:, CV_NODI)
          LOC_T2_J( : ) = T2_ALL(:, CV_NODJ)
          LOC_T2OLD_I( : ) = T2OLD_ALL(:, CV_NODI)
          LOC_T2OLD_J( : ) = T2OLD_ALL(:, CV_NODJ)
       END IF
!------------------

       If_GOT_DIFFUS2: IF ( GOT_DIFFUS ) THEN
          ! This sub caculates the effective diffusion
          ! coefficient DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
          T_ALL_J( : )   =LOC_T_J( : )
          TOLD_ALL_J( : )=LOC_TOLD_J( : )
          LOC_WIC_T_BC_ALL(:)=0
          IF(SELE.NE.0) THEN
             DO IPHASE=1,NPHASE
                LOC_WIC_T_BC_ALL(IPHASE)=WIC_T_BC_ALL(1, IPHASE, SELE)
                IF(LOC_WIC_T_BC_ALL(IPHASE)==WIC_T_BC_DIRICHLET) THEN
                   T_ALL_J( IPHASE ) = SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*( SELE- 1) )
                   TOLD_ALL_J( IPHASE )=SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*( SELE- 1) )
                ENDIF 
             END DO
          ENDIF
          CALL DIFFUS_CAL_COEFF( DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX,  &
               CV_NLOC, MAT_NLOC, CV_NONODS, NPHASE, TOTELE, MAT_NONODS, MAT_NDGLN, &
               SCVFEN, SCVFEN, SCVNGI, GI, NDIM, TDIFFUSION, DUMMY_ZERO_NDIM_NDIM_NPHASE, &
               HDC, &
               T_ALL_J( : ), LOC_T_I( : ), &
               TOLD_ALL_J( : ), LOC_TOLD_I( : ), &
               ELE, ELE2, CVNORMX_ALL( :, GI ), &
               DTX_ELE_ALL(:,:,:,ELE), DTOLDX_ELE_ALL(:,:,:,ELE),  DTX_ELE_ALL(:,:,:,MAX(1,ELE2)), DTOLDX_ELE_ALL(:,:,:,MAX(ELE2,1)), &
               SELE, STOTEL, LOC_WIC_T_BC_ALL, CV_OTHER_LOC, MAT_OTHER_LOC )
       ELSE
          DIFF_COEF_DIVDX = 0.0
          DIFF_COEFOLD_DIVDX = 0.0
       END IF If_GOT_DIFFUS2



                   IF ( IS_OVERLAPPING.or. is_compact_overlapping ) THEN

                      DO IPHASE=1,NPHASE
                         DO IDIM=1,NDIM
                            DO JDIM=1,NDIM
                               IJ=(IPHASE-1)*MAT_NONODS*NDIM*NDIM + (MAT_NODJ-1)*NDIM*NDIM + (IDIM-1)*NDIM +JDIM
                               VJ_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) = OPT_VEL_UPWIND_COEFS(IJ) 
                               GJ_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) = OPT_VEL_UPWIND_COEFS(IJ+NPHASE*MAT_NONODS*NDIM*NDIM) 
                            END DO
                         END DO
                      END DO
                      IF( is_compact_overlapping ) THEN ! The inverse of the sigma matrix...
                         INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,:) = INV_V_OPT_VEL_UPWIND_COEFS(:,:,:,MAT_NODJ)
                      ENDIF

                   ENDIF

                     NFACE_ITS = 1
                     FACE_ITS = 1
!                     DO FACE_ITS = 1, NFACE_ITS
                        ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI.
                        IF(IGOT_T2==1) THEN
                           CALL GET_INT_VEL_NEW( NPHASE, NDOTQNEW,  NDOTQOLD, INCOMEOLD, &
                                HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS,  &
                                LOC_T2OLD_I, LOC_T2OLD_J, LOC_FEMT2OLD,LOC2_FEMT2OLD, LOC_DENOLD_I, LOC_DENOLD_J, &
                                LOC_U, LOC2_U, LOC_NUOLD, LOC2_NUOLD, SLOC_NUOLD, &
                                CV_NODI, CV_NODJ, CVNORMX_ALL, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
                                SUF_SIG_DIAGTEN_BC, &
                                UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                ONE_PORE(ELE), ONE_PORE(MAX(1,ELE2)), CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS, &
       MASS_CV(CV_NODI), MASS_CV(CV_NODJ), NDIM, MAT_NLOC, MAT_NONODS, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                IANISOLIM,  &
                                NUOLDGI_ALL, T2OLDUPWIND_MAT_ALL( :, COUNT_IN), T2OLDUPWIND_MAT_ALL( :, COUNT_OUT) ) 
                            CALL GET_INT_VEL_NEW( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
                                LOC_T2_I, LOC_T2_J, LOC_FEMT2, LOC2_FEMT2, LOC_DEN_I, LOC_DEN_J, &
                                LOC_U, LOC2_U, LOC_NU, LOC2_NU, SLOC_NU, &
                                CV_NODI, CV_NODJ, CVNORMX_ALL, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
                                SUF_SIG_DIAGTEN_BC, &
                                UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                ONE_PORE(ELE), ONE_PORE(MAX(1,ELE2)), CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS, &
       MASS_CV(CV_NODI), MASS_CV(CV_NODJ), NDIM, MAT_NLOC, MAT_NONODS, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                IANISOLIM,  &
                                NUGI_ALL, T2UPWIND_MAT_ALL( :, COUNT_IN), T2UPWIND_MAT_ALL( :, COUNT_OUT) )
                        ELSE
                           CALL GET_INT_VEL_NEW( NPHASE, NDOTQNEW, NDOTQOLD, INCOMEOLD, &
                                HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
                                LOC_TOLD_I, LOC_TOLD_J, LOC_FEMTOLD, LOC2_FEMTOLD, LOC_DENOLD_I, LOC_DENOLD_J, &
                                LOC_U, LOC2_U, LOC_NUOLD, LOC2_NUOLD, SLOC_NUOLD, &
                                CV_NODI, CV_NODJ, CVNORMX_ALL, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
                                SUF_SIG_DIAGTEN_BC, &
                                UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                ONE_PORE(ELE), ONE_PORE(MAX(1,ELE2)), CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS, &
       MASS_CV(CV_NODI), MASS_CV(CV_NODJ), NDIM, MAT_NLOC, MAT_NONODS, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                IANISOLIM,  &
                                NUOLDGI_ALL, TOLDUPWIND_MAT_ALL( :, COUNT_IN), TOLDUPWIND_MAT_ALL( :, COUNT_OUT) ) 
                           CALL GET_INT_VEL_NEW( NPHASE, NDOTQNEW, NDOTQ, INCOME, &
                                HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
                                LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
                                LOC_U, LOC2_U, LOC_NU, LOC2_NU, SLOC_NU, &
                                CV_NODI, CV_NODJ, CVNORMX_ALL, &
                                CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
                                SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
                                SUF_SIG_DIAGTEN_BC, &
                                UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
                                ONE_PORE(ELE), ONE_PORE(MAX(1,ELE2)), CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS, &
       MASS_CV(CV_NODI), MASS_CV(CV_NODJ), NDIM, MAT_NLOC, MAT_NONODS, &
                                IN_ELE_UPWIND, DG_ELE_UPWIND, &
                                IANISOLIM, &
                                NUGI_ALL, TUPWIND_MAT_ALL( :, COUNT_IN), TUPWIND_MAT_ALL( :, COUNT_OUT) )
                        ENDIF

                     INCOME_J=1.-INCOME
                     INCOMEold_J=1.-INCOMEold


!          if(.true.) then
! Pack ndotq information: 
             IPT=1
             CALL PACK_LOC( F_INCOME(:), INCOME( : ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) ) ! t
             CALL PACK_LOC( F_INCOME(:), INCOMEOLD( : ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) ) ! TOLD
             CALL PACK_LOC( F_INCOME(:), INCOME( : ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )  ! d
             CALL PACK_LOC( F_INCOME(:), INCOMEOLD( : ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )  ! DOLD
             CALL PACK_LOC( F_INCOME(:), INCOME( : ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )  ! T2
             CALL PACK_LOC( F_INCOME(:), INCOMEOLD( : ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )  ! T2OLD
             IPT=1
             CALL PACK_LOC( F_NDOTQ(:), NDOTQ( : ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,1) ) ! t
             CALL PACK_LOC( F_NDOTQ(:), NDOTQOLD( : ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,2) ) ! TOLD
             CALL PACK_LOC( F_NDOTQ(:), NDOTQ( : ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,3) )  ! d
             CALL PACK_LOC( F_NDOTQ(:), NDOTQOLD( : ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,4) )  ! DOLD
             CALL PACK_LOC( F_NDOTQ(:), NDOTQ( : ),    NPHASE, NFIELD, IPT, IGOT_T_PACK(:,5) )  ! T2
             CALL PACK_LOC( F_NDOTQ(:), NDOTQOLD( : ), NPHASE, NFIELD, IPT, IGOT_T_PACK(:,6) )  ! T2OLD
!          endif

!             print *,'F_INCOME:',F_INCOME
!             print *,'F_NDOTQ:',F_NDOTQ



                        !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
                        ! Calculate T and DEN on the CV face at quadrature point GI.
                  IF(NFIELD.GT.0) THEN

                        CALL GET_INT_T_DEN_new( LIMF(:), & 
                             CV_DISOPT, CV_NONODS, NPHASE, NFIELD, CV_NODI, CV_NODJ, CV_ILOC, CV_JLOC, CV_SILOC, ELE, ELE2, GI,   &
                             CV_NLOC, TOTELE, CV_OTHER_LOC, SCVNGI, SCVFEN, F_INCOME, F_NDOTQ, &
                             LOC_F, LOC_FEMF, SLOC_F, SLOC_FEMF, SLOC2_F, SLOC2_FEMF, &
                             SELE, CV_SNLOC,  U_SNLOC,  STOTEL, CV_SLOC2LOC, SLOC_SUF_F_BC, &
                             U_SLOC2LOC, U_OTHER_LOC,  &
                             SELE_LOC_WIC_F_BC,   &
                             WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET,  & 
                             HDC, DT, &
                             SCVFENX_ALL(1:NDIM,:,:), CVNORMX_ALL, &
                             LOC_UF,  &
                             U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
                             FUPWIND_IN, FUPWIND_OUT, DISTCONTINUOUS_METHOD, QUAD_ELEMENTS, SHAPE_CV_SNL, DOWNWIND_EXTRAP_INDIVIDUAL, &
                             F_CV_NODI, F_CV_NODJ) 
                  ENDIF
! it does not matter about bcs for FVT below as its zero'ed out in the eqns:
                        FVT(:)=T_ALL(:,CV_NODI)*(1.0-INCOME(:)) + T_ALL(:,CV_NODJ)*INCOME(:) 
!                        FVD(:)=DEN_ALL(:,CV_NODI)*(1.0-INCOME(:)) + DEN_ALL(:,CV_NODJ)*INCOME(:) 


!           print *,'done temp'
! Generate some local F variables ***************
!       IF(.true.) THEN
! loc_f - Unpack into the limiting variables LIMT and may be store them in the cache.

             !###############TEMPORARY USAGE IN UNPACK_LOC #######################
             !Currently for UNPACK_LOC we are passing SCVNGI*TOTELE as the maximum Global_face value
             !However, that is an overstimate. THAT VALUE WILL NEED TO BE CHANGED!!
             !###############################################################
             IPT=1
             CALL UNPACK_LOC( LIMF(:), LIMT( : ),    NPHASE, NFIELD, IPT, STORE, IGOT_T_PACK(:,1), GLOBAL_FACE, IGOT_T_CONST(:,1), IGOT_T_CONST_VALUE(:,1),&
                SCVNGI*TOTELE,state, 'limf1', StorageIndexes(34) )
             CALL UNPACK_LOC( LIMF(:), LIMTOLD( : ), NPHASE, NFIELD, IPT, STORE, IGOT_T_PACK(:,2), GLOBAL_FACE, IGOT_T_CONST(:,2), IGOT_T_CONST_VALUE(:,2),&
                SCVNGI*TOTELE,state, 'limf2', StorageIndexes(35) )
             CALL UNPACK_LOC( LIMF(:), LIMD( : ),    NPHASE, NFIELD, IPT, STORE, IGOT_T_PACK(:,3), GLOBAL_FACE, IGOT_T_CONST(:,3), IGOT_T_CONST_VALUE(:,3),&
                SCVNGI*TOTELE,state, 'limf3', StorageIndexes(36) )
             CALL UNPACK_LOC( LIMF(:), LIMDOLD( : ), NPHASE, NFIELD, IPT, STORE, IGOT_T_PACK(:,4), GLOBAL_FACE, IGOT_T_CONST(:,4), IGOT_T_CONST_VALUE(:,4),&
                SCVNGI*TOTELE,state, 'limf4', StorageIndexes(37) )
             CALL UNPACK_LOC( LIMF(:), LIMT2( : ),    NPHASE, NFIELD, IPT, STORE, IGOT_T_PACK(:,5), GLOBAL_FACE, IGOT_T_CONST(:,5), IGOT_T_CONST_VALUE(:,5),&
                SCVNGI*TOTELE,state, 'limf5', StorageIndexes(38) )
             CALL UNPACK_LOC( LIMF(:), LIMT2OLD( : ), NPHASE, NFIELD, IPT, STORE, IGOT_T_PACK(:,6), GLOBAL_FACE, IGOT_T_CONST(:,6), IGOT_T_CONST_VALUE(:,6),&
                SCVNGI*TOTELE,state, 'limf6', StorageIndexes(39) )


                IF(GETCT.AND.RETRIEVE_SOLID_CTY) THEN
                   NDOTQ_HAT = 0.0
                   DO U_KLOC = 1, U_NLOC
                      IF (ELE2/=0) THEN ! Between elements...
                         NDOTQ_HAT =  NDOTQ_HAT + SUFEN( U_KLOC, GI ) * 0.5 * SUM( CVNORMX_ALL(:, GI) * (LOC_U_HAT( :, U_KLOC ) + LOC2_U_HAT( :, U_KLOC )) )
                      ELSE
                         NDOTQ_HAT =  NDOTQ_HAT + SUFEN( U_KLOC, GI ) * SUM( CVNORMX_ALL(:, GI) * LOC_U_HAT( :, U_KLOC ) )
                      ENDIF
                   END DO
                   
                   if(.false.) then ! use an average to make sure all is well in terms of propagation of information...
                      LIMT(:)=0.5*(LOC_T_I( : )+LOC_T_J( : )) 
                   endif
                      
                   DO IPHASE=1,NPHASE
                      LIMT_HAT(IPHASE) = MAX(1.E-7,LIMT(IPHASE))
                   END DO
                   R=SUM(LIMT_HAT(:))
                   LIMT_HAT(:)=LIMT_HAT(:)/R

                   if(sele.ne.0) then ! effectively apply the bcs to NDOTQ_HAT
                     NDOTQ_HAT =SUM(LIMT_HAT(:)*NDOTQNEW(:))
                   endif
                ENDIF

    ! Amend for porosity...
          IF ( ELE2 /= 0 ) THEN 
!             FVD   = 0.5 * ( ONE_PORE(ELE) + ONE_PORE(ELE2) ) * FVD
             LIMD   = 0.5 * ( ONE_PORE(ELE) + ONE_PORE(ELE2) ) * LIMD
             LIMDOLD   = 0.5 * ( ONE_PORE(ELE) + ONE_PORE(ELE2) ) * LIMDOLD
          ELSE
!             FVD   = ONE_PORE(ELE) * FVD
             LIMD   = ONE_PORE(ELE) * LIMD
             LIMDOLD   = ONE_PORE(ELE) * LIMDOLD
          END IF


          LIMDT=LIMD*LIMT
          LIMDTOLD=LIMDOLD*LIMTOLD

          LIMDTT2=LIMD*LIMT*LIMT2
          LIMDTT2OLD=LIMDOLD*LIMTOLD*LIMT2OLD

!      ENDIF
! Generate some local F variables ***************...

! Make allowances for no matrix stencil operating from outside the boundary.
          BCZERO=1.0
          IF( SELE > 0 ) BCZERO=1.0-INCOME

                  Loop_IPHASE: DO IPHASE = 1, NPHASE

                     CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
                     CV_NODJ_IPHA = CV_NODJ + ( IPHASE - 1 ) * CV_NONODS
                     RHS_NODI_IPHA = IPHASE +  (CV_NODI - 1 ) * NPHASE
                     RHS_NODJ_IPHA = IPHASE +  (CV_NODJ - 1 ) * NPHASE


                     ! Define face value of theta
                     IF ( IGOT_T2 == 1 ) THEN
                        FTHETA = FACE_THETA( DT, CV_THETA, ( CV_DISOPT>=8 ), HDC, NDOTQ(IPHASE), LIMDTT2(IPHASE), DIFF_COEF_DIVDX(IPHASE), &
                             T_ALL( IPHASE, CV_NODJ ) * DEN_ALL( IPHASE, CV_NODJ ) * T2_ALL( IPHASE, CV_NODJ ), &
                             T_ALL( IPHASE, CV_NODI ) * DEN_ALL( IPHASE, CV_NODI ) * T2_ALL( IPHASE, CV_NODI ), &
                             NDOTQOLD(IPHASE), LIMDTT2OLD(IPHASE), DIFF_COEFOLD_DIVDX(IPHASE), &
                             TOLD_ALL( IPHASE, CV_NODJ ) * DENOLD_ALL( IPHASE, CV_NODJ ) * T2OLD_ALL( IPHASE, CV_NODJ ), &
                             TOLD_ALL( IPHASE, CV_NODI ) * DENOLD_ALL( IPHASE, CV_NODI ) * T2OLD_ALL( IPHASE, CV_NODI ) )
                     ELSE
                        FTHETA = FACE_THETA( DT, CV_THETA, ( CV_DISOPT>=8 ), HDC, NDOTQ(IPHASE), LIMDTT2(IPHASE), DIFF_COEF_DIVDX(IPHASE), &
                             T_ALL( IPHASE, CV_NODJ ) * DEN_ALL( IPHASE, CV_NODJ ), &
                             T_ALL( IPHASE, CV_NODI ) * DEN_ALL( IPHASE, CV_NODI ), &
                             NDOTQOLD(IPHASE), LIMDTT2OLD(IPHASE), DIFF_COEFOLD_DIVDX(IPHASE), &
                             TOLD_ALL( IPHASE, CV_NODJ ) * DENOLD_ALL( IPHASE, CV_NODJ ), &
                             TOLD_ALL( IPHASE, CV_NODI ) * DENOLD_ALL( IPHASE, CV_NODI ) )
                     END IF

                     FTHETA_T2 = FTHETA * LIMT2(IPHASE)
                     ONE_M_FTHETA_T2OLD = (1.0-FTHETA) * LIMT2OLD(IPHASE)

                     FTHETA_T2_J = FTHETA * LIMT2(IPHASE)
                     ONE_M_FTHETA_T2OLD_J = (1.0-FTHETA) * LIMT2OLD(IPHASE)

                     IF(IGOT_THETA_FLUX == 1) THEN
                        IF ( GET_THETA_FLUX ) THEN
                           THETA_FLUX( IPHASE, GLOBAL_FACE ) = FTHETA * LIMDT(IPHASE) / DEN_ALL( IPHASE, CV_NODI )
                           ONE_M_THETA_FLUX( IPHASE, GLOBAL_FACE ) = (1.0-FTHETA) * LIMDTOLD(IPHASE) / DEN_ALL( IPHASE, CV_NODI )
                           if(integrate_other_side) then ! for the flux on the other side of the CV face...
                              THETA_FLUX_J( IPHASE, GLOBAL_FACE ) = FTHETA * LIMDT(IPHASE) / DEN_ALL( IPHASE, CV_NODJ )
                              ONE_M_THETA_FLUX_J( IPHASE, GLOBAL_FACE ) = (1.0-FTHETA) * LIMDTOLD(IPHASE) / DEN_ALL( IPHASE, CV_NODJ )
                           endif
                        END IF
                        IF ( USE_THETA_FLUX ) THEN
                           FTHETA_T2 = THETA_FLUX( IPHASE, GLOBAL_FACE )
                           ONE_M_FTHETA_T2OLD = ONE_M_THETA_FLUX( IPHASE, GLOBAL_FACE )
                           if(integrate_other_side) then ! for the flux on the other side of the CV face...
                              FTHETA_T2_J = THETA_FLUX_J( IPHASE, GLOBAL_FACE )
                              ONE_M_FTHETA_T2OLD_J = ONE_M_THETA_FLUX_J( IPHASE, GLOBAL_FACE )
!                              FTHETA_T2_J = THETA_FLUX( IPHASE, GLOBAL_FACE )
!                              ONE_M_FTHETA_T2OLD_J = ONE_M_THETA_FLUX( IPHASE, GLOBAL_FACE )
                           endif
                        END IF
                     END IF


                     ROBIN1=0.0
                     ROBIN2=0.0
                     IF ( SELE /= 0 ) THEN
                        IF ( WIC_T_BC_ALL(1,IPHASE,SELE) == WIC_T_BC_ROBIN ) THEN
! this needs to be corrected (its correct but misleading)...
                           ROBIN1 = SUF_T_BC_ROB1_ALL(1,iphase, CV_SILOC+CV_SNLOC*(sele-1))
                           ROBIN2 = SUF_T_BC_ROB2_ALL (1,iphase, CV_SILOC+CV_SNLOC*(sele-1))
                        END IF
                     END IF

                     !====================== ACV AND RHS ASSEMBLY ===================
                     Conditional_GETCT2 : IF ( GETCT ) THEN ! Obtain the CV discretised CT eqations plus RHS
                        CALL PUT_IN_CT_RHS( CT, CT_RHS, U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
                             CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE, ELE, ELE2, SELE, &
                             JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, U_OTHER_LOC, U_NDGLN, U_ALL, &
                             SUFEN, SCVDETWEI, CVNORMX_ALL, DEN_ALL, CV_NODI, CV_NODJ, &
                             UGI_COEF_ELE_ALL,  &
                             UGI_COEF_ELE2_ALL,  &
                             NDOTQNEW(IPHASE), NDOTQOLD(IPHASE), NDOTQ_HAT, LIMD(IPHASE), LIMT(IPHASE), LIMTOLD(IPHASE), LIMDT(IPHASE), LIMDTOLD(IPHASE), LIMT_HAT(IPHASE), &
                             FTHETA_T2, ONE_M_FTHETA_T2OLD, FTHETA_T2_J, ONE_M_FTHETA_T2OLD_J, integrate_other_side_and_not_boundary, &
                             RETRIEVE_SOLID_CTY,theta_cty_solid)
                     ENDIF Conditional_GETCT2


                     Conditional_GETCV_DISC: IF ( GETCV_DISC ) THEN
                        ! Obtain the CV discretised advection/diffusion equations
                        IF ( GETMAT ) THEN

                           ! - Calculate the integration of the limited, high-order flux over a face
                           ! Conservative discretisation. The matrix (PIVOT ON LOW ORDER SOLN)
                           IF ( ( CV_NODI /= CV_NODJ ) .AND. ( CV_NODJ /= 0 ) ) THEN
                              CSR_ACV( IPHASE+(JCOUNT_IPHA-1)*NPHASE ) =  CSR_ACV( IPHASE+(JCOUNT_IPHA-1)*NPHASE ) &
                                   + SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * INCOME(IPHASE) * LIMD(IPHASE) & ! Advection
                                   - FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) ! Diffusion contribution
! integrate the other CV side contribution (the sign is changed)...
                        if(integrate_other_side_and_not_boundary) then
                              CSR_ACV( IPHASE+(ICOUNT_IPHA-1)*NPHASE ) =  CSR_ACV( IPHASE+(ICOUNT_IPHA-1)*NPHASE ) &
                                   - SECOND_THETA * FTHETA_T2_J * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * INCOME_J(IPHASE) * LIMD(IPHASE) & ! Advection
                                   - FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) ! Diffusion contribution
                        endif

                              IF ( GET_GTHETA ) THEN
                                 THETA_GDIFF( IPHASE, CV_NODI ) =  THETA_GDIFF( IPHASE, CV_NODI ) &
                                      + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) * T_ALL( IPHASE, CV_NODJ ) ! Diffusion contribution
! integrate the other CV side contribution (the sign is changed)...
                        if(integrate_other_side_and_not_boundary) then
                                 THETA_GDIFF( IPHASE, CV_NODJ ) =  THETA_GDIFF( IPHASE, CV_NODJ ) &
                                      + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) * T_ALL( IPHASE, CV_NODI ) ! Diffusion contribution
                        endif
                              END IF
                           ELSE IF ( SELE /= 0 ) THEN
                              IF(WIC_T_BC_ALL(1,iphase,sele) == WIC_T_BC_DIRICHLET) THEN
                                 CV_RHS( RHS_NODI_IPHA ) =  CV_RHS( RHS_NODI_IPHA ) &
                                      + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) &
                                         * SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*( SELE- 1) )

                                 IF(GET_GTHETA) THEN
                                    THETA_GDIFF( IPHASE, CV_NODI ) =  THETA_GDIFF( IPHASE, CV_NODI ) &
                                         + FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) &
                                           * SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*( SELE- 1) )
                                 END IF
                              END IF
                           END IF


                           IMID_IPHA = IPHASE + (IMID-1)*NPHASE
                           JMID_IPHA = IPHASE + (JMID-1)*NPHASE

                           CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA ) &
                                +  SECOND_THETA * FTHETA_T2 * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * ( 1. - INCOME(IPHASE) ) * LIMD(IPHASE) & ! Advection
                                +  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE)  &  ! Diffusion contribution
                                +  SCVDETWEI( GI ) * ROBIN1  ! Robin bc
                        if(integrate_other_side_and_not_boundary) then
                           CSR_ACV( JMID_IPHA ) =  CSR_ACV( JMID_IPHA ) &
                                -  SECOND_THETA * FTHETA_T2_J * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * ( 1. - INCOME_J(IPHASE) ) * LIMD(IPHASE) & ! Advection
                                +  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE)    ! Diffusion contribution
                        endif

                           IF ( GET_GTHETA ) THEN
                              THETA_GDIFF( IPHASE, CV_NODI ) =  THETA_GDIFF( IPHASE, CV_NODI ) &
                                   -  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) * T_ALL( IPHASE, CV_NODI ) & ! Diffusion contribution
                                   -  SCVDETWEI( GI ) * ROBIN1 * T_ALL( IPHASE, CV_NODI )  ! Robin bc
                        if(integrate_other_side_and_not_boundary) then
                              THETA_GDIFF( IPHASE, CV_NODJ ) =  THETA_GDIFF( IPHASE, CV_NODJ ) &
                                   -  FTHETA * SCVDETWEI( GI ) * DIFF_COEF_DIVDX(IPHASE) * T_ALL( IPHASE, CV_NODJ ) ! Diffusion contribution
                        endif
                           END IF

                           ! CV_BETA=0 for Non-conservative discretisation (CV_BETA=1 for conservative disc)
                           CSR_ACV( IMID_IPHA ) = CSR_ACV( IMID_IPHA )  &
                                - SECOND_THETA * FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * LIMD(IPHASE)
                        if(integrate_other_side_and_not_boundary) then
                           CSR_ACV( JMID_IPHA ) = CSR_ACV( JMID_IPHA )  &
                                + SECOND_THETA * FTHETA_T2_J * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * LIMD(IPHASE)
                        endif
                        END IF

                        TMID = T_ALL( IPHASE, CV_NODI )
                        TOLDMID = TOLD_ALL( IPHASE, CV_NODI )
                        TMID_J = T_ALL( IPHASE, CV_NODJ )
                        TOLDMID_J = TOLD_ALL( IPHASE, CV_NODJ )

!                        ! Make allowances for no matrix stencil operating from outside the boundary.
!                        BCZERO(IPHASE) = 1.0
!                        IF ( (SELE /= 0) .AND. (INCOME(IPHASE) > 0.5) ) BCZERO(IPHASE)=0.0

                        ! Put results into the RHS vector
                        CV_RHS( RHS_NODI_IPHA ) =  CV_RHS( RHS_NODI_IPHA )  &
                                ! subtract 1st order adv. soln.
                             + SECOND_THETA * FTHETA_T2 * NDOTQNEW(IPHASE) * SCVDETWEI( GI ) * LIMD(IPHASE) * FVT(IPHASE) * BCZERO(IPHASE) &
                             -  SCVDETWEI( GI ) * ( FTHETA_T2 * NDOTQNEW(IPHASE) * LIMDT(IPHASE) &
                                + ONE_M_FTHETA_T2OLD * NDOTQOLD(IPHASE) * LIMDTOLD(IPHASE) ) ! hi order adv
                        if(integrate_other_side_and_not_boundary) then
                        CV_RHS( RHS_NODJ_IPHA ) =  CV_RHS( RHS_NODJ_IPHA )  &
                                ! subtract 1st order adv. soln.
                             - SECOND_THETA * FTHETA_T2_J * NDOTQNEW(IPHASE) * SCVDETWEI( GI ) * LIMD(IPHASE) * FVT(IPHASE) * BCZERO(IPHASE) &
                             +  SCVDETWEI( GI ) * ( FTHETA_T2_J * NDOTQNEW(IPHASE) * LIMDT(IPHASE) &
                                + ONE_M_FTHETA_T2OLD_J * NDOTQOLD(IPHASE) * LIMDTOLD(IPHASE) ) ! hi order adv
                        endif

                        ! Subtract out 1st order term non-conservative adv.
                        CV_RHS( RHS_NODI_IPHA ) =  CV_RHS( RHS_NODI_IPHA ) &
                             - FTHETA_T2 * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * LIMD(IPHASE) * TMID
                        if(integrate_other_side_and_not_boundary) then
                        CV_RHS( RHS_NODJ_IPHA ) =  CV_RHS( RHS_NODJ_IPHA ) &
                             + FTHETA_T2_J * ( 1. - CV_BETA ) * SCVDETWEI( GI ) * NDOTQNEW(IPHASE) * LIMD(IPHASE) * TMID_J
                        endif

                        ! High-order non-conservative advection contribution
                        CV_RHS( RHS_NODI_IPHA ) =  CV_RHS( RHS_NODI_IPHA ) &
                             + ( 1. - CV_BETA) * SCVDETWEI( GI ) &
                             * ( FTHETA_T2 * NDOTQNEW(IPHASE) * TMID * LIMD(IPHASE)  &
                                + ONE_M_FTHETA_T2OLD * NDOTQOLD(IPHASE) * LIMDOLD(IPHASE) * TOLDMID )
                        if(integrate_other_side_and_not_boundary) then
                        CV_RHS( RHS_NODJ_IPHA ) =  CV_RHS( RHS_NODJ_IPHA ) &
                             - ( 1. - CV_BETA) * SCVDETWEI( GI ) &
                             * ( FTHETA_T2_J * NDOTQNEW(IPHASE) * TMID_J * LIMD(IPHASE)  &
                                + ONE_M_FTHETA_T2OLD_J * NDOTQOLD(IPHASE) * LIMDOLD(IPHASE) * TOLDMID_J )
                        endif

                        ! Diffusion contribution
                        CV_RHS( RHS_NODI_IPHA ) =  CV_RHS( RHS_NODI_IPHA ) &
                             + (1.-FTHETA) * SCVDETWEI(GI) * DIFF_COEFOLD_DIVDX(IPHASE) &
                             * ( TOLD_ALL( IPHASE, CV_NODJ ) - TOLD_ALL( IPHASE, CV_NODI ) ) &
                                ! Robin bc
                             + SCVDETWEI( GI ) * ROBIN2
                        if(integrate_other_side_and_not_boundary) then
                        CV_RHS( RHS_NODJ_IPHA ) =  CV_RHS( RHS_NODJ_IPHA ) &
                             + (1.-FTHETA) * SCVDETWEI(GI) * DIFF_COEFOLD_DIVDX(IPHASE) &
                             * ( TOLD_ALL( IPHASE, CV_NODI ) - TOLD_ALL( IPHASE, CV_NODJ ) ) 
                        endif
                        IF ( GET_GTHETA ) THEN
                           THETA_GDIFF( IPHASE, CV_NODI ) =  THETA_GDIFF( IPHASE, CV_NODI ) &
                                + (1.-FTHETA) * SCVDETWEI(GI) * DIFF_COEFOLD_DIVDX(IPHASE) &
                                * ( TOLD_ALL( IPHASE, CV_NODJ ) - TOLD_ALL( IPHASE, CV_NODI ) ) &
                                ! Robin bc
                                + SCVDETWEI( GI ) * ROBIN2
                        if(integrate_other_side_and_not_boundary) then
                           THETA_GDIFF( IPHASE, CV_NODJ ) =  THETA_GDIFF( IPHASE, CV_NODJ ) &
                                + (1.-FTHETA) * SCVDETWEI(GI) * DIFF_COEFOLD_DIVDX(IPHASE) &
                                * ( TOLD_ALL( IPHASE, CV_NODI ) - TOLD_ALL( IPHASE, CV_NODJ ) ) 
                        endif
                        END IF

                        ! this is for the internal energy equation source term..
                        IF ( THERMAL ) THEN

                           THERM_FTHETA = 1.

                           IF( RETRIEVE_SOLID_CTY ) THEN
                              VOL_FRA_FLUID_I = VOL_FRA_FLUID(CV_NODI)
                              VOL_FRA_FLUID_J = VOL_FRA_FLUID(CV_NODJ)
                           ELSE
                              VOL_FRA_FLUID_I = 1.0
                              VOL_FRA_FLUID_J = 1.0
                           ENDIF


                           IF ( IGOT_T2 /= 0 ) THEN

                              CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                                   - CV_P( CV_NODI ) * SCVDETWEI( GI ) * ( &
                                   THERM_FTHETA * NDOTQNEW( IPHASE ) * LIMT2( IPHASE ) &
                                   + ( 1. - THERM_FTHETA ) * NDOTQOLD(IPHASE) * LIMT2OLD( IPHASE ) )*VOL_FRA_FLUID_I
                              if ( integrate_other_side_and_not_boundary ) then
                                 CV_RHS( RHS_NODJ_IPHA ) = CV_RHS( RHS_NODJ_IPHA ) &
                                      + CV_P( CV_NODJ ) * SCVDETWEI( GI ) * ( &
                                      THERM_FTHETA * NDOTQNEW(IPHASE) * LIMT2(IPHASE) &
                                      + ( 1. - THERM_FTHETA ) * NDOTQOLD(IPHASE) * LIMT2OLD(IPHASE) )*VOL_FRA_FLUID_J
                              end if

                           ELSE

                              CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                                   - CV_P( CV_NODI ) * SCVDETWEI( GI ) * ( &
                                   THERM_FTHETA * NDOTQNEW( IPHASE ) &
                                   + ( 1. - THERM_FTHETA ) * NDOTQOLD(IPHASE) )*VOL_FRA_FLUID_I
                              if ( integrate_other_side_and_not_boundary ) then
                                 CV_RHS( RHS_NODJ_IPHA ) = CV_RHS( RHS_NODJ_IPHA ) &
                                      + CV_P( CV_NODJ ) * SCVDETWEI( GI ) * ( &
                                      THERM_FTHETA * NDOTQNEW( IPHASE ) &
                                      + ( 1. - THERM_FTHETA ) * NDOTQOLD( IPHASE ) )*VOL_FRA_FLUID_J
                              end if
                           
                           END IF !IGOT_T2 

                           IF ( GOT_VIS ) THEN
                              ! stress form of viscosity...
                              NU_LEV_GI(:, IPHASE) = ( (1.-THERM_FTHETA) * NUOLDGI_ALL(:,IPHASE) + THERM_FTHETA * NUGI_ALL(:,IPHASE) )

                              STRESS_IJ_THERM(:,:,IPHASE) = 0.0 

                              CALL CALC_STRESS_TEN( STRESS_IJ_THERM(:,:,IPHASE), ZERO_OR_TWO_THIRDS, NDIM, &
                                   CVNORMX_ALL(:,GI), NU_LEV_GI(:,IPHASE) * SCVDETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,MAT_NODI) )
                                   !UFENX_ALL(1:NDIM,U_ILOC,GI), UFENX_ALL(1:NDIM,U_JLOC,GI) * DETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,CV_NODI) )

                              if ( integrate_other_side_and_not_boundary ) then
                                 STRESS_IJ_THERM_J(:,:,IPHASE) = 0.0 
                                 CALL CALC_STRESS_TEN( STRESS_IJ_THERM_J(:,:,IPHASE), ZERO_OR_TWO_THIRDS, NDIM, &
                                   CVNORMX_ALL(:,GI), NU_LEV_GI(:,IPHASE) * SCVDETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,MAT_NODJ) )
                                   !UFENX_ALL(1:NDIM,U_ILOC,GI), UFENX_ALL(1:NDIM,U_JLOC,GI) * DETWEI(GI), THERM_U_DIFFUSION(:,:,IPHASE,CV_NODJ) )
                              end if

                              DO IDIM = 1, NDIM
                                 DO JDIM = 1, NDIM
                                    VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODI) = VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODI) + STRESS_IJ_THERM(IDIM,JDIM,IPHASE)
                                    VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODI) = VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODI) + NU_LEV_GI(IDIM,IPHASE) * CVNORMX_ALL(JDIM,GI) * SCVDETWEI(GI)
                                    if ( integrate_other_side_and_not_boundary ) then
                                       VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODJ) = VECS_STRESS(IDIM,JDIM,IPHASE,CV_NODJ) - STRESS_IJ_THERM_J(IDIM,JDIM,IPHASE )
                                       VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODJ) = VECS_GRAD_U(IDIM,JDIM,IPHASE,CV_NODJ) - NU_LEV_GI(IDIM,IPHASE) * CVNORMX_ALL(JDIM,GI) * SCVDETWEI(GI)
                                    end if
                                 END DO
                              END DO
                           END IF ! GOT_VIS

                        END IF ! THERMAL

                     ENDIF Conditional_GETCV_DISC

                  END DO Loop_IPHASE

           endif ! if(CV_NODJ.ge.CV_NODI) then

               END IF Conditional_integration

            END DO Loop_GCOUNT

         END DO Loop_CV_ILOC

      END DO Loop_Elements




      IF(GET_GTHETA) THEN
         DO CV_NODI = 1, CV_NONODS
            DO IPHASE = 1, NPHASE
!               CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
               THETA_GDIFF(IPHASE, CV_NODI) = THETA_GDIFF(IPHASE, CV_NODI) / MASS_CV(CV_NODI)
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

            Loop_IPHASE2: DO IPHASE = 1, NPHASE
!               CV_NODI_IPHA = CV_NODI + ( IPHASE - 1 ) * CV_NONODS
               RHS_NODI_IPHA = IPHASE + ( CV_NODI -1 ) * NPHASE
               IMID_IPHA = IPHASE + (SMALL_CENTRM(CV_NODI)-1)*NPHASE

                R = MEAN_PORE_CV( CV_NODI ) * MASS_CV( CV_NODI ) / DT

               IF(THERMAL) THEN
                  IF(GOT_VIS) THEN
                     IF( RETRIEVE_SOLID_CTY ) THEN
                       CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                         + VOL_FRA_FLUID(cv_nodi)*SUM( VECS_STRESS(:,:,IPHASE,CV_NODI)*VECS_GRAD_U(:,:,IPHASE,CV_NODI)  )/MASS_CV(CV_NODI) 
                     else
                       CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                         + SUM( VECS_STRESS(:,:,IPHASE,CV_NODI)*VECS_GRAD_U(:,:,IPHASE,CV_NODI)  )/MASS_CV(CV_NODI) 
                     endif
                  ENDIF
               ENDIF

               IF ( IGOT_T2 == 1 ) THEN
                  CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                       + MASS_CV(CV_NODI) * SOURCT_ALL( IPHASE, CV_NODI )

                  CSR_ACV( IMID_IPHA ) = CSR_ACV( IMID_IPHA ) &
                       + (CV_BETA * DEN_ALL( IPHASE, CV_NODI ) * T2_ALL( IPHASE, CV_NODI ) &
                       + (1.-CV_BETA) * DEN_ALL( IPHASE, CV_NODI ) * T2_ALL( IPHASE, CV_NODI ) ) &
                       * R

                  CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                       + (CV_BETA * DENOLD_ALL( IPHASE, CV_NODI ) * T2OLD_ALL( IPHASE, CV_NODI ) &
                       + (1.-CV_BETA) * DEN_ALL( IPHASE, CV_NODI ) * T2_ALL( IPHASE, CV_NODI ) )  &
                       * R * TOLD_ALL( IPHASE, CV_NODI )
               ELSE

                  CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                       + MASS_CV( CV_NODI ) * SOURCT_ALL( IPHASE, CV_NODI )

                  CSR_ACV( IMID_IPHA ) =  CSR_ACV( IMID_IPHA ) &
                       + (CV_BETA * DEN_ALL( IPHASE, CV_NODI ) &
                       + (1.-CV_BETA) * DEN_ALL( IPHASE, CV_NODI ) )  &
                       * R

                  CV_RHS( RHS_NODI_IPHA ) = CV_RHS( RHS_NODI_IPHA ) &
                       + ( CV_BETA * DENOLD_ALL( IPHASE, CV_NODI ) &
                       + (1.-CV_BETA) * DEN_ALL( IPHASE, CV_NODI ) ) &
                       * R * TOLD_ALL( IPHASE, CV_NODI )
               END IF


               Conditional_GETMAT2: IF ( GETMAT ) THEN

                  DO JPHASE = 1, NPHASE
                     DENSE_ACV( JPHASE, IPHASE, CV_NODI )  = DENSE_ACV( JPHASE, IPHASE, CV_NODI ) &
                          + MASS_CV( CV_NODI ) * ABSORBT_ALL( IPHASE, JPHASE, CV_NODI )
                  END DO

               END IF Conditional_GETMAT2

            END DO Loop_IPHASE2

         END DO Loop_CVNODI2

      END IF Conditional_GETCV_DISC2

      IF ( GETCT ) THEN
      
      
         W_SUM_ONE1 = 1.0 !If == 1.0 applies constraint to T
         W_SUM_ONE2 = 0.0 !If == 1.0 applies constraint to TOLD


         DIAG_SCALE_PRES = 0.0

         DO CV_NODI = 1, CV_NONODS

            R = MASS_CV( CV_NODI ) * MEAN_PORE_CV( CV_NODI ) / DT

! Add constraint to force sum of volume fracts to be unity...
                  ! W_SUM_ONE==1 applies the constraint
                  ! W_SUM_ONE==0 does NOT apply the constraint
            CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - ( W_SUM_ONE1 - W_SUM_ONE2 ) * R

            IF(RETRIEVE_SOLID_CTY) THEN 
! VOL_FRA_FLUID is the old voln fraction of total fluid...
! multiply by solid-voln fraction: (1.-VOL_FRA_FLUID)
               CT_RHS( CV_NODI ) = CT_RHS( CV_NODI )  + (1.-VOL_FRA_FLUID( CV_NODI )) * R

            ENDIF

            DO IPHASE = 1, NPHASE

               CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) &
                    - R * ( &
                    + (1.0-W_SUM_ONE1) * T_ALL( IPHASE, CV_NODI ) - (1.0-W_SUM_ONE2) * TOLD_ALL( IPHASE, CV_NODI ) &
                    + ( TOLD_ALL( IPHASE, CV_NODI ) * ( DEN_ALL( IPHASE, CV_NODI ) - DENOLD_ALL( IPHASE, CV_NODI ) ) &
                    - DERIV( IPHASE, CV_NODI ) * CV_P( CV_NODI ) * T_ALL_KEEP( IPHASE, CV_NODI ) ) / DEN_ALL( IPHASE, CV_NODI ) )
                    !- DERIV( IPHASE, CV_NODI ) * CV_P( CV_NODI ) * max( 1., T_ALL( IPHASE, CV_NODI ) ) ) / DEN_ALL( IPHASE, CV_NODI ) )


               DIAG_SCALE_PRES( CV_NODI ) = DIAG_SCALE_PRES( CV_NODI ) + &
                    MEAN_PORE_CV( CV_NODI ) * T_ALL_KEEP( IPHASE, CV_NODI ) * DERIV( IPHASE, CV_NODI ) &
                    !MEAN_PORE_CV( CV_NODI ) * max( 1., T_ALL( IPHASE, CV_NODI ) ) * DERIV( IPHASE, CV_NODI ) &
                    / ( DT * DEN_ALL( IPHASE, CV_NODI ) )

               CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) + MASS_CV( CV_NODI ) * SOURCT_ALL( IPHASE, CV_NODI ) / DEN_ALL( IPHASE, CV_NODI )

               DO JPHASE = 1, NPHASE
                  CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) &
                       - MASS_CV( CV_NODI ) * ABSORBT_ALL( IPHASE, JPHASE, CV_NODI ) * T_ALL( JPHASE, CV_NODI ) / DEN_ALL( IPHASE, CV_NODI )
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
!      T_FEMT = FEMT
!      DEN_FEMT = FEMDEN

      ewrite(3,*) '----------sub cv_assemb--------'
      if( .false. .and. getct) then
         ewrite(3,*) 'ct_rhs:', ct_rhs
      end if
      if( .false. .and. GETCV_DISC ) then
         ewrite(3,*) 'cv_rhs:', cv_rhs
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
           T, X_NDGLN, CV_NDGLN, X_ALL(1,:))
      ewrite(3,*) 'just print out - in cv_assemb'

      ! Deallocating temporary working arrays

      IF(GETCT) THEN
         DEALLOCATE( JCOUNT_KLOC )
      ENDIF
      DEALLOCATE( CVNORMX )
      DEALLOCATE( CVNORMY )
      DEALLOCATE( CVNORMZ )
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

      DEALLOCATE( SRA )

      DEALLOCATE( CV_SLOC2LOC )
      DEALLOCATE( U_SLOC2LOC )

!      DEALLOCATE( FEMT )
!      DEALLOCATE( FEMTOLD )
!      DEALLOCATE( FEMDEN )
!      DEALLOCATE( FEMDENOLD )
      DEALLOCATE( MASS_CV )
      DEALLOCATE( XC_CV_ALL )
      DEALLOCATE( FACE_ELE )
      DEALLOCATE(TUPWIND_MAT_ALL)
      DEALLOCATE(TOLDUPWIND_MAT_ALL)
      DEALLOCATE(DENUPWIND_MAT_ALL)
      DEALLOCATE(DENOLDUPWIND_MAT_ALL)
      DEALLOCATE(T2UPWIND_MAT_ALL)
      DEALLOCATE(T2OLDUPWIND_MAT_ALL)

      call deallocate(tracer_BCs)
      call deallocate(tracer_BCs_robin2)
      call deallocate(density_BCs)
      call deallocate(velocity_BCs)
      if (present(saturation)) then
         call deallocate(saturation_BCs)
         call deallocate(saturation_BCs_robin2)
      end if

      deallocate (t, told)!TEMPORARY
if (present(T_input)) then!<==TEMPORARY
      deallocate( T_ALL_TARGET, TOLD_ALL_TARGET, FEMT_ALL_TARGET, FEMTOLD_ALL_TARGET)
end if


      ewrite(3,*) 'Leaving CV_ASSEMB'

      RETURN

    END SUBROUTINE CV_ASSEMB






             SUBROUTINE IS_FIELD_CONSTANT(IGOT_T_CONST, IGOT_T_CONST_VALUE, T_ALL, CV_NONODS)
             LOGICAL IGOT_T_CONST
             REAL IGOT_T_CONST_VALUE
             INTEGER CV_NONODS
             REAL T_ALL(CV_NONODS)
             real, parameter :: tolerance = 1.e-10, tolerance2 = 1.e-5
! Local variables...
             REAL RMAX,RMIN,PERCENT
             INTEGER I

             RMIN=+1.e+20
             RMAX=-1.e+20
             DO I=1,CV_NONODS
                RMAX=MAX(RMAX, T_ALL(I))
                RMIN=MIN(RMIN, T_ALL(I))
             END DO

             PERCENT = ABS(RMAX-RMIN)/MAX(ABS(RMIN),ABS(RMAX),tolerance)
             
             IGOT_T_CONST=PERCENT < tolerance2

             IGOT_T_CONST_VALUE = RMAX
             RETURN
             END SUBROUTINE IS_FIELD_CONSTANT
             




             SUBROUTINE PACK_LOC( LOC_F, T_ALL, NPHASE, NFIELD, IPT, IGOT_T_PACK ) 
! If PACK then pack T_ALL into LOC_F as long at IGOT_T==1 and STORE and not already in storage. 
             IMPLICIT NONE
             INTEGER, intent( in ) :: NPHASE, NFIELD
! GLOBAL_FACE is the quadrature point which helps point into the storage memory
             INTEGER, intent( inout ) :: IPT
             LOGICAL, DIMENSION(NPHASE), intent( in ) :: IGOT_T_PACK
             REAL, DIMENSION(NPHASE), intent( in ) :: T_ALL
             REAL, DIMENSION(NFIELD), intent( inout ) :: LOC_F
! local variables...
             INTEGER :: IPHASE

             DO IPHASE=1,NPHASE
                IF(IGOT_T_PACK(IPHASE)) THEN ! Put into packing vector LOC_F
!                   print *,'nphase,nfield,ipt:',nphase,nfield,ipt
                   LOC_F(IPT) = T_ALL(IPHASE)
                   IPT=IPT+1
                ENDIF
             END DO

             RETURN
             END SUBROUTINE PACK_LOC




             SUBROUTINE I_PACK_LOC( LOC_F, T_ALL, NPHASE, NFIELD, IPT, IGOT_T_PACK ) 
! If PACK then pack T_ALL into LOC_F as long at IGOT_T==1 and STORE and not already in storage. 
             IMPLICIT NONE
             INTEGER, intent( in ) :: NPHASE, NFIELD
! GLOBAL_FACE is the quadrature point which helps point into the storage memory
             INTEGER, intent( inout ) :: IPT
             LOGICAL, DIMENSION(NPHASE), intent( in ) :: IGOT_T_PACK
             INTEGER, DIMENSION(NPHASE), intent( in ) :: T_ALL
             INTEGER, DIMENSION(NFIELD), intent( inout ) :: LOC_F
! local variables...
             INTEGER :: IPHASE

             DO IPHASE=1,NPHASE
                IF(IGOT_T_PACK(IPHASE)) THEN ! Put into packing vector LOC_F
                   LOC_F(IPT) = T_ALL(IPHASE)
                   IPT=IPT+1
                ENDIF
             END DO

             RETURN
             END SUBROUTINE I_PACK_LOC


             SUBROUTINE UNPACK_LOC( LOC_F, T_ALL, NPHASE, NFIELD, IPT, STORE, IGOT_T_PACK, GLOBAL_FACE, IGOT_T_CONST, IGOT_T_CONST_VALUE,&
                TOTAL_GLOBAL_FACE, state, StorName, indx )
! If PACK then UNpack loc_f into T_ALL  as long at IGOT_T==1 and STORE and not already in storage.
             IMPLICIT NONE
             LOGICAL, intent( in ) :: STORE
             INTEGER, intent( in ) :: NPHASE,NFIELD, TOTAL_GLOBAL_FACE
             INTEGER, intent( in ) :: GLOBAL_FACE
! GLOBAL_FACE is the quadrature point which helps point into the storage memory
             INTEGER, intent( inout ) :: IPT
             LOGICAL, DIMENSION(NPHASE), intent( in ) :: IGOT_T_PACK, IGOT_T_CONST
             REAL, DIMENSION(NPHASE), intent( inout ) :: T_ALL
             REAL, DIMENSION(NPHASE), intent( in ) :: IGOT_T_CONST_VALUE
             REAL, DIMENSION(NFIELD), intent( inout ) :: LOC_F
             type( state_type ), intent( inout ), dimension(:) :: state
             character(len=*), intent(in) :: StorName
             integer, intent(inout) :: indx
             !Local variables
             !Variables to store things in state
             type(mesh_type), pointer :: fl_mesh
             type(mesh_type) :: Auxmesh
             type(scalar_field), target :: targ_fieldToStore
! local variables...
             INTEGER :: IPHASE

            !Initial check to know whether we have already stored all the values
            if (indx < 0) then
                if (GLOBAL_FACE<state(1)%scalar_fields(-indx)%ptr%val(size(state(1)%scalar_fields(-indx)%ptr%val,1))) then
                    !If the input index is smaller than the last stored index, then we have re-started and we should be extracting stored values
                    indx = abs(indx)
                end if
            end if

             DO IPHASE=1,NPHASE
                 IF(IGOT_T_PACK(IPHASE)) THEN
                     T_ALL(IPHASE) = LOC_F(IPT)
                      IPT=IPT+1
                     if (.not.IGOT_T_CONST(IPHASE)) then
                         IF(STORE) THEN ! Put in storage if not already in storage...
                             if (indx < 0) then!We may need to store a new value
                                 !Store in state, indx is an input
                                 state(1)%scalar_fields(-indx)%ptr%val(IPHASE+(GLOBAL_FACE-1)*NPHASE) = T_ALL(IPHASE)
                                 !Store input index to check when we start to get data instead of sting data
                                 state(1)%scalar_fields(-indx)%ptr%val(size(state(1)%scalar_fields(-indx)%ptr%val,1)) = GLOBAL_FACE
                             else if (GLOBAL_FACE==1) then !The first time we need to introduce the targets in state
                                 if (has_scalar_field(state(1), "Fld"//StorName)) then
                                     !If we are recalculating due to a mesh modification then
                                     !we return to the original situation
                                     call remove_scalar_field(state(1), "Fld"//StorName)
                                 end if
                                  !Get mesh file just to be able to allocate the fields we want to store
                                 fl_mesh => extract_mesh( state(1), "CoordinateMesh" )
                                 Auxmesh = fl_mesh
                                 !The number of nodes I want does not coincide
                                 !############################################################################
                                 !TOTAL_GLOBAL_FACE*NPHASE, SEEMS AN OVERSTIMATE, MAYBE TOTAL_GLOBAL_FACE*NFIELD IS ENOUGH???????
                                 !############################################################################
                                 Auxmesh%nodes = TOTAL_GLOBAL_FACE*NPHASE+1!(+1 to store the last input index)
                                 call allocate (targ_fieldToStore, Auxmesh)

                                 !Now we insert them in state and store the indexes
                                 call insert(state(1), targ_fieldToStore, "Fld"//StorName)
                                 !Store index with a negative value, because if the index is
                                 !zero or negative then we have to calculate stuff
                                 indx = -size(state(1)%scalar_fields)

                                  !Store in state
                                 state(1)%scalar_fields(-indx)%ptr%val(IPHASE+(GLOBAL_FACE-1)*NPHASE) = T_ALL(IPHASE)
                             end if
                         end if
                     ENDIF
                 ELSE IF(IGOT_T_CONST(IPHASE)) THEN
                     T_ALL(IPHASE) = IGOT_T_CONST_VALUE(IPHASE)
                 ELSE IF(STORE.and.indx>0) THEN ! Check storage and pull out of storage
                     T_ALL(IPHASE) = state(1)%scalar_fields(indx)%ptr%val(IPHASE+(GLOBAL_FACE-1)*NPHASE)
                 ELSE ! Set to 1 as last resort e.g. for T2, T2OLD
                     T_ALL(IPHASE) = 1.0
                 ENDIF
             END DO

             RETURN
             END SUBROUTINE UNPACK_LOC


             SUBROUTINE PACK_OR_UNPACK_LOC( LOC_F, T_ALL, NPHASE, NFIELD, IPT, PACK, STORE, IGOT_T , GLOBAL_FACE)
                 ! If PACK then pack T_ALL into LOC_F as long at IGOT_T==1 and STORE and not already in storage.
                 LOGICAL, intent( in ) :: STORE, PACK
                 INTEGER, intent( in ) :: NPHASE, IGOT_T
                 INTEGER, intent( in ) :: GLOBAL_FACE
                 ! GLOBAL_FACE is the quadrature point which helps point into the storage memory
                 INTEGER, intent( inout ) :: IPT
                 REAL, DIMENSION(NPHASE), intent( inout ) :: T_ALL
                 REAL, DIMENSION(NFIELD), intent( inout ) :: LOC_F
                 ! local variables...
                 LOGICAL :: IN_STORAGE

                 IN_STORAGE = .FALSE.
                 IF(STORE) THEN
                     ! IF STORE then look to see if in storage
                     IN_STORAGE = .FALSE.
                 ENDIF
             
                 IF(IGOT_T==1) THEN
                     IF(PACK) THEN
                         ! Pack solution into LOC_F
                         IF(.NOT.IN_STORAGE) THEN
                             LOC_F(IPT:IPT-1+NPHASE) = T_ALL(:)
                             IPT=IPT+NPHASE
                         ENDIF
                     ELSE
                         ! Unpack...
                         IF(STORE) THEN
                             IF(.NOT.IN_STORAGE) THEN ! See if we are already storing limited value
                             ENDIF
                             ! Put LOC_F(1:NPHASE, CV_KLOC) = DEN_ALL( 1:NPHASE, CV_NODK ) into storage...
                             !                      T_ALL = LOC_F ???
                             IPT=IPT+NPHASE
                         ELSE
                             ! Put LOC_F(1:NPHASE, CV_KLOC) = DEN_ALL( 1:NPHASE, CV_NODK ) into storage...
                             T_ALL(:) = LOC_F(IPT:IPT-1+NPHASE)
                             IPT=IPT+NPHASE
                         ENDIF
               
                     ENDIF ! ENDOF IF(PACK) THEN ELSE
                 ENDIF ! END OF IF(IGOT_T==1) THEN
                 RETURN
             END SUBROUTINE PACK_OR_UNPACK_LOC








    function CV_count_faces( packed_state,&
         CV_ELE_TYPE,  &
         STOTEL, CV_SNDGLN, U_SNDGLN, &
         face_sparsity) result(global_face)

      !  =====================================================================
      !     This subroutine counts then number of faces in the control volume space
      !

      ! Inputs/Outputs
      IMPLICIT NONE
      type(state_type), intent(inout) :: packed_state
      INTEGER, intent( in ) :: CV_ELE_TYPE, STOTEL
      INTEGER, DIMENSION( : ), pointer :: CV_NDGLN
      INTEGER, DIMENSION( : ), pointer ::  X_NDGLN
      INTEGER, DIMENSION( : ), pointer :: U_NDGLN
      INTEGER, DIMENSION( : ), pointer :: XU_NDGLN
      INTEGER, DIMENSION( : ), pointer :: MAT_NDGLN
      INTEGER, DIMENSION( : ) :: CV_SNDGLN
      INTEGER, DIMENSION(: )  :: U_SNDGLN
      ! Diagonal scaling of (distributed) pressure matrix (used to treat pressure implicitly)
      INTEGER, DIMENSION( : ), pointer :: FINDCMC, COLCMC

      INTEGER, DIMENSION( : ), pointer :: FINDM, COLM, MIDM
      INTEGER, DIMENSION( : ), pointer :: FINELE, COLELE
      integer, dimension(:), pointer :: SMALL_FINDRM, SMALL_COLM, SMALL_CENTRM
      !character( len = option_path_len ), intent( in ), optional :: option_path_spatial_discretisation


      ! Local variables
      integer :: CV_NONODS, U_NONODS, X_NONODS, TOTELE, &
           CV_NLOC, U_NLOC, X_NLOC, MAT_NLOC, &
           NDIM, XU_NLOC, cv_snloc, u_snloc
      LOGICAL, DIMENSION( : ), allocatable :: X_SHARE
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
           CVNORMX, &
           CVNORMY, CVNORMZ, SCVRA, SCVDETWEI, SRA, &
           SUM_CV, ONE_PORE, SELE_OVERLAP_SCALE, &
           UP_WIND_NOD, DU, DV, DW, PERM_ELE
      REAL, DIMENSION( : , : ), allocatable :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, &
           CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT,  &
           UFEN, UFENLX, UFENLY, UFENLZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
           SCVFENLX, SCVFENLY, SCVFENLZ, &
           SCVFENX, SCVFENY, SCVFENZ, &
           SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
           SBCVN,SBCVFEN, SBCVFENSLX, SBCVFENSLY, &
           SBCVFENLX, SBCVFENLY, SBCVFENLZ, SBUFEN, SBUFENSLX, SBUFENSLY, &
           SBUFENLX, SBUFENLY, SBUFENLZ

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

      !        ===>  LOGICALS  <===
      LOGICAL :: QUAD_OVER_WHOLE_ELE,integrat_at_gi

    INTEGER :: GLOBAL_FACE
    type(csr_sparsity), pointer :: connectivity
    type(mesh_type), pointer :: cv_mesh, x_mesh, xu_mesh, u_mesh, mat_mesh

    type(csr_sparsity), intent(out), optional :: face_sparsity
    type(ilist), dimension(:), allocatable :: face_list 

    GLOBAL_FACE=0

    connectivity=>extract_csr_sparsity(packed_state,"ElementConnectivity")
    cv_mesh=>extract_mesh(packed_state,"PressureMesh")
    cv_ndgln=>cv_mesh%ndglno
    totele=element_count(cv_mesh)
    ndim=mesh_dim(cv_mesh)
    cv_nonods=node_count(cv_mesh)
    cv_nloc=ele_loc(cv_mesh,1)
    cv_snloc=face_loc(cv_mesh,1)
    x_mesh=>extract_mesh(packed_state,"PressureMesh_Continuous")
    x_ndgln=>x_mesh%ndglno
    X_nonods=node_count(x_mesh)
    X_nloc=ele_loc(x_mesh,1)
    xu_mesh=>extract_mesh(packed_state,"VelocityMesh_Continuous")
    xu_ndgln=>xu_mesh%ndglno
    xu_nloc=ele_loc(xu_mesh,1)
    u_mesh=>extract_mesh(packed_state,"InternalVelocityMesh")
    u_ndgln=>u_mesh%ndglno
    u_nonods=node_count(u_mesh)
    if (is_overlapping) then
       u_nloc=ele_loc(u_mesh,1)*cv_nloc
       u_snloc=face_loc(u_mesh,1)*cv_nloc
    else
       u_nloc=ele_loc(u_mesh,1)
       u_snloc=face_loc(u_mesh,1)
    end if
    mat_mesh=>extract_mesh(packed_state,"PressureMesh_Discontinuous")
    mat_nloc=ele_loc(mat_mesh,1) 

    allocate( face_list( totele ) )

      ewrite(3,*) 'In CV_FACE_COUNT'

      QUAD_OVER_WHOLE_ELE=.FALSE.
      ! If QUAD_OVER_WHOLE_ELE=.true. then dont divide element into CV's to form quadrature.
      call retrieve_ngi( TOTELE, cv_ele_type, CV_NLOC,U_NLOC, &
           cv_ngi, cv_ngi_short, scvngi, sbcvngi, nface, QUAD_OVER_WHOLE_ELE )

                  global_face=totele * SCVNGI * 2
                  return


     
      return
    end function CV_COUNT_FACES





    SUBROUTINE FIND_OTHER_SIDE( CV_OTHER_LOC, CV_NLOC, CV_NODI, U_OTHER_LOC, U_NLOC,  &
         MAT_OTHER_LOC, MAT_NLOC, INTEGRAT_AT_GI, &
         X_NLOC, XU_NLOC, X_NDGLN, CV_NDGLN, XU_NDGLN, &
         CV_SNLOC, CVFEM_ON_FACE, X_SHARE, X_NONODS, ELE, ELE2,  &
         FINELE, COLELE, NCOLELE, DISTCONTINUOUS_METHOD ) 
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
      LOGICAL, intent( in ) :: DISTCONTINUOUS_METHOD
      ! Local variables
      INTEGER :: X_KLOC, X_NODK, X_NODK2, COUNT, ELE3, SUF_COUNT, CV_KLOC, CV_KLOC2, &
           &     U_KLOC, U_KLOC2, CV_NODK, XU_NODK, XU_NODK2, ILEV, JLEV
      LOGICAL :: INTEGRAT_AT_GI2

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
         IF( SUF_COUNT == CV_SNLOC ) THEN
            ELE3 = ELE2
            EXIT
         ENDIF
         !ewrite(3,*)'suf_count:', ele, ele2, suf_count, cv_snloc
      END DO
      ELE2 = ELE3

      DO X_KLOC = 1, X_NLOC
         X_NODK = X_NDGLN( ( ELE - 1 ) * X_NLOC + X_KLOC )
         X_SHARE( X_NODK ) = .FALSE.
      END DO


! Quite because there is no work to do here...
      IF(.NOT.DISTCONTINUOUS_METHOD) THEN
         IF(ELE2.NE.0) THEN ! this is not on the boundary of the domain.
            INTEGRAT_AT_GI=.FALSE.
            RETURN
         ENDIF
      ENDIF 
 

!      IF ( ELE2 /= 0 ) THEN 
!         ! Is CV_NODI in element ELE2 if yes set ELE2=0 as we don't want to integrate in 
!         ! the middle of a CV. 
!         DO CV_KLOC = 1, CV_NLOC
!            CV_NODK = CV_NDGLN( ( ELE2 - 1 ) * CV_NLOC + CV_KLOC )
!            !ewrite(3,*)'cv_nodi, cv_nodk:', cv_nodi, cv_nodk, ele, ele2
!            IF( CV_NODK == CV_NODI ) THEN
!               INTEGRAT_AT_GI = .FALSE.
!               stop 2911
!               EXIT
!            ENDIF
!         END DO
!      END IF

      IF ( ( ELE2 /= 0 ) .AND. INTEGRAT_AT_GI ) THEN ! Determine CV_OTHER_LOC(CV_KLOC)
         CV_OTHER_LOC = 0
         DO CV_KLOC = 1, CV_NLOC
            IF ( CVFEM_ON_FACE( CV_KLOC ) ) THEN ! Find opposite local node
               X_NODK = X_NDGLN( ( ELE - 1 ) * X_NLOC + CV_KLOC )
               DO CV_KLOC2 = 1, CV_NLOC
                  X_NODK2 = X_NDGLN( ( ELE2 - 1 ) * X_NLOC + CV_KLOC2 )
                  IF( X_NODK2 == X_NODK ) THEN
                     CV_OTHER_LOC( CV_KLOC ) = CV_KLOC2
                     EXIT
                  ENDIF
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






    subroutine proj_cv_to_fem_n( packed_state, cvn, n, nlx, nly, nlz, cvweight, &
         &                       cv_ngi, findm, colm, ncolm )

    implicit none

   
    ! DO NOT CALL THIS SUB...NOT FINISHED..

    integer, intent( in ) :: cv_ngi, ncolm
    type( state_type ), intent( inout ) :: packed_state
    real, dimension( :, : ), intent( in ) :: cvn, n, nlx, nly, nlz
    real, dimension( : ), intent( in ) :: cvweight
    integer, dimension( : ), intent( in ) :: findm, colm

    type( tensor_field ), pointer :: field
    type( vector_field ), pointer :: x
    type( scalar_field ), pointer :: pressure
    real, dimension( cv_ngi ) :: detwei
    integer, dimension( : ), pointer :: x_ndgln, cv_ndgln
    real, dimension( : ), allocatable :: rhs, mat, ra
    real, dimension( :, : ), allocatable :: xx
    real, dimension( :, :, : ), allocatable :: nx

    integer :: nfields, nphase, ncomp, cv_nonods, npsi, &
         ifield, iphase, icomp, ele, totele, x_nonods, &
         cv_iloc, cv_jloc, cv_nodi, cv_nodj, count, ipsi, &
         idx, idim, ndim, cv_nloc
    real :: nn, nm, volume
    character( len = option_path_len ) :: path
    character( len = field_name_len ), dimension( : ), allocatable :: cv_name
    character( len = field_name_len ) :: fe_name


    path = "/material_phase[0]/scalar_field::Pressure" 

    nfields = tensor_field_count( packed_state )

    pressure => extract_scalar_field( packed_state, "Pressure" )
    
    cv_nloc = ele_loc( pressure, 1 )


    x => extract_vector_field( packed_state, "PressureCoordinate" )
    x_nonods = node_count( x )


    ndim = x%dim

    allocate( nx( cv_nloc, cv_ngi, 3   ) ) ; xx = 0.0


    allocate( xx( x_nonods, 3 ) ) ; xx = 0.0
    do idim = 1, ndim
       xx( :, 1 ) = x%val( :, 1 )
    end do

    x_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh_Continuous" ) )
    cv_ndgln => get_ndglno( extract_mesh( packed_state, "PressureMesh" ) )

    cv_nonods = node_count( pressure )
    totele = ele_count( pressure )
    
    allocate( cv_name( nfields ) )


    npsi = 0
    do ifield = 1, nfields
       field => extract_tensor_field( packed_state, ifield )
       if ( trim( field%mesh%name ) == "PressureMesh" ) then
          nphase = size( field%val, 2 ) ; ncomp =  size( field%val, 1 )
          npsi = npsi + nphase * ncomp
          cv_name( ifield ) = trim( field%name )
       end if
    end do

    allocate( mat( ncolm ), ra( cv_ngi ), xx( x_nonods, 3 ) )
    allocate( rhs( npsi * cv_nonods ) ) ; rhs = 0.0

    do ele = 1, totele

       call detnlxr( ele, xx(:,1), xx(:,2), xx(:,3), x_ndgln, totele, x_nonods, cv_nloc, cv_ngi, &
            n, nlx, nly, nlz, cvweight, detwei, ra, volume, (ndim==1), (ndim==3), .false., &
            nx( :, :, 1 ), nx( :, :, 2 ), nx( :, :, 3 ) )

       do cv_iloc = 1, cv_nloc
          cv_nodi = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_iloc )

          do cv_jloc = 1, cv_nloc
             cv_nodj = cv_ndgln( ( ele - 1 ) * cv_nloc + cv_jloc )
       
             nn = dot_product( n( cv_iloc, : ), n(   cv_jloc,: ) * detwei )
             nm = dot_product( n( cv_iloc, : ), cvn( cv_jloc, : ) * detwei )

             call posinmat( count, cv_nodi, cv_nodj, cv_nonods, findm, colm, ncolm )
             mat( count ) = mat( count ) + nn

             ipsi = 1
             do ifield = 1, nfields
                field => extract_tensor_field( packed_state, ifield  )
                if ( trim( field%mesh%name ) == "PressureMesh" ) then
                   do iphase = 1, nphase
                      do icomp = 1, ncomp
                         idx = cv_nodi + (ipsi-1)*cv_nonods
                         rhs( idx ) = rhs( idx ) +  field%val( icomp, iphase, cv_nodj )
                         ipsi = ipsi + 1
                      end do
                   end do
                end if
             end do

          end do ! jloc
       end do ! iloc
    end do ! ele


    stop 666

    ! solve...
    idx = 1
    do ifield = 1, nfields

       fe_name = trim( cv_name( ifield ) ) ! come up with something clever here...


       ! check strings...

       fe_name = "PackedFE" // trim( cv_name( ifield )(5:) )



       field => extract_tensor_field( packed_state, trim( fe_name ) )

       !this is not good either..
       if ( trim( field%mesh%name ) == "PressureMesh" ) then
          do iphase = 1, nphase
             do icomp = 1, ncomp
                call solver( mat, field%val( icomp, iphase, :  ), &
                     rhs( 1 + ( idx - 1 ) * cv_nonods : idx * cv_nonods ), &
                     findm, colm, option_path = trim( path ) )
                idx = idx + 1
             end do
          end do
       end if
    end do

    deallocate( rhs, mat, ra, xx, nx )

    return

  end subroutine proj_cv_to_fem_n





  SUBROUTINE PROJ_CV_TO_FEM_4( state, &
       FEMT_ALL, FEMTOLD_ALL, FEMDEN_ALL, FEMDENOLD_ALL, T_ALL, TOLD_ALL, DEN_ALL, DENOLD_ALL, &
       IGOT_T2,T2_ALL,T2OLD_ALL, FEMT2_ALL,FEMT2OLD_ALL, &
       XC_CV_ALL, MASS_CV, MASS_ELE, &
       NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
       CV_NGI, CV_NLOC, CVN, CVWEIGHT, N, NLX, NLY, NLZ, &
       X_NONODS, X_ALL, NCOLM, FINDM, COLM, MIDM, &
       IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC )

    ! determine FEMT (finite element wise) etc from T (control volume wise)
    IMPLICIT NONE
    type( state_type ), dimension( : ), intent( in ) :: state
    INTEGER, intent( in ) :: NDIM, NPHASE, CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
         X_NONODS, NCOLM, IGOT_T2, IGETCT, NCOLCMC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    REAL, DIMENSION( :,: ), intent( inout ) :: FEMT_ALL, FEMTOLD_ALL, FEMDEN_ALL, FEMDENOLD_ALL
    REAL, DIMENSION( :,: ), intent( inout ) :: FEMT2_ALL, FEMT2OLD_ALL
    REAL, DIMENSION( :, : ), intent( in ) :: T_ALL, TOLD_ALL, DEN_ALL, DENOLD_ALL
    REAL, DIMENSION( :, : ), intent( in ) :: T2_ALL, T2OLD_ALL
    REAL, DIMENSION( : ), intent( inout ) :: MASS_CV
    REAL, DIMENSION( :, : ), intent( inout ) :: XC_CV_ALL ! (NDIM,X_NONODS)
    REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE
    REAL, DIMENSION( :, : ), intent( in ) :: CVN
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( :, : ), intent( in ) :: X_ALL
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
    INTEGER :: velocity_max_iterations, nstates, istate, iphase, k
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
!    PSI( 1 + 0 * NL : NL + 0 * NL )  =      T( 1 : NL )
!    PSI( 1 + 1 * NL : NL + 1 * NL )  =   TOLD( 1 : NL )
!    PSI( 1 + 2 * NL : NL + 2 * NL )  =    DEN( 1 : NL )
!    PSI( 1 + 3 * NL : NL + 3 * NL )  = DENOLD( 1 : NL )

!    IF(IGOT_T2==1) THEN
!       PSI( 1 + 4 * NL : NL + 4 * NL )  =      T2( 1 : NL )
!       PSI( 1 + 5 * NL : NL + 5 * NL )  =   T2OLD( 1 : NL )
!    ENDIF
         k = 1
         DO IPHASE = 1, NPHASE
             DO CV_INOD = 1, CV_NONODS
                 PSI( k ) = T_ALL( IPHASE, CV_INOD) 
                 PSI( k + NL ) = TOLD_ALL( IPHASE, CV_INOD) 
                 PSI( k + 2*NL ) = DEN_ALL( IPHASE, CV_INOD) 
                 PSI( k + 3*NL ) = DENOLD_ALL( IPHASE, CV_INOD) 
                 k = k + 1
             END DO
         END DO
         if (IGOT_T2>0) then
             k = 1
             DO IPHASE = 1, NPHASE
                 DO CV_INOD = 1, CV_NONODS
                     PSI( k + 4*NL ) = T2_ALL( IPHASE, CV_INOD) 
                     PSI( k + 5*NL ) = T2OLD_ALL( IPHASE, CV_INOD) 
                     k = k + 1
                 END DO
             END DO
         end if


    DO ELE=1,TOTELE
       DO CV_ILOC=1,CV_NLOC
          X_INOD = X_NDGLN((ELE-1)*X_NLOC +CV_ILOC)
          CV_INOD=CV_NDGLN((ELE-1)*CV_NLOC+CV_ILOC)
          !####TEMPORARY UNTIL WE VECTORIZE EVERYTHING####
          select case (ndim)
              case (1)
                  PSI_AVE(CV_INOD)            =X_ALL(1,X_INOD)
              case (2)
                  PSI_AVE(CV_INOD)            =X_ALL(1,X_INOD)
                  PSI_AVE(CV_INOD+CV_NONODS)  = X_ALL(2,X_INOD)
              case default
                  PSI_AVE(CV_INOD)            =X_ALL(1,X_INOD)
                  PSI_AVE(CV_INOD+CV_NONODS)  = X_ALL(2,X_INOD)
                  PSI_AVE(CV_INOD+2*CV_NONODS)= X_ALL(3,X_INOD)
          end select
          !####TEMPORARY UNTIL WE VECTORIZE EVERYTHING####
       END DO
    END DO
    PSI_INT=1.0

    FEMPSI = PSI

    CALL PROJ_CV_TO_FEM( FEMPSI, PSI, NTSOL, NDIM, &
         PSI_AVE,NTSOL_AVE, PSI_INT,NTSOL_INT, MASS_ELE, &
         CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
         CV_NGI, CV_NLOC, CVN, CVWEIGHT, N, NLX, NLY, NLZ, &
         X_NONODS, X_ALL(1,:), X_ALL(2,:), X_ALL(3,:), NCOLM, FINDM, COLM, MIDM, &
         IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC, PATH )

         NL=CV_NONODS*NPHASE


         k = 1
         DO IPHASE = 1, NPHASE
             DO CV_INOD = 1, CV_NONODS
                 FEMT_ALL( IPHASE, CV_INOD) = FEMPSI( k )
                 FEMTOLD_ALL( IPHASE, CV_INOD) = FEMPSI( k + NL )
                 FEMDEN_ALL( IPHASE, CV_INOD) = FEMPSI( k + 2*NL )
                 FEMDENOLD_ALL( IPHASE, CV_INOD) = FEMPSI( k + 3*NL )
                 k = k + 1
             END DO
         END DO
         if (IGOT_T2>0) then
             k = 1
             DO IPHASE = 1, NPHASE
                 DO CV_INOD = 1, CV_NONODS
                     FEMT2_ALL( IPHASE, CV_INOD) = FEMPSI( k + 4*NL )
                     FEMT2OLD_ALL( IPHASE, CV_INOD) = FEMPSI( k + 5*NL )
                     k = k + 1
                 END DO
             END DO
         end if



    XC_CV_ALL( 1, 1 : CV_NONODS ) = PSI_AVE( 1 : CV_NONODS )
    IF(NDIM.GE.2) XC_CV_ALL( 2, 1 : CV_NONODS ) = PSI_AVE( 1 +CV_NONODS:   2*CV_NONODS )
    IF(NDIM.GE.3) XC_CV_ALL( 3, 1 : CV_NONODS ) = PSI_AVE( 1 +2*CV_NONODS: 3*CV_NONODS )
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

!    FEMPSI=PSI

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


  SUBROUTINE PROJ_CV_TO_FEM_state( packed_state,FEMPSI, PSI, NDIM, &
       PSI_AVE, PSI_INT, MASS_ELE, &
       CV_NONODS, TOTELE, CV_NDGLN, X_NLOC, X_NDGLN, &
       CV_NGI, CV_NLOC, CVN, CVWEIGHT, N, NLX, NLY, NLZ, &
       X_NONODS, X, NCOLM, FINDM, COLM, MIDM, &
       IGETCT, MASS_MN_PRES, FINDCMC, COLCMC, NCOLCMC )

    ! Determine FEMT (finite element wise) etc from T (control volume wise)
    ! Also integrate PSI_INT over each CV and average PSI_AVE over each CV. 
    use shape_functions 
    use shape_functions_Linear_Quadratic
    use solvers_module
    use matrix_operations
    IMPLICIT NONE

    type(state_type) :: packed_state
    type(tensor_field_pointer), dimension(:) :: fempsi,psi
    type(vector_field_pointer), dimension(:) :: psi_int,  psi_ave

    INTEGER, intent( in ) :: NDIM, CV_NONODS, TOTELE, &
         X_NLOC, CV_NGI, CV_NLOC, X_NONODS, NCOLM, IGETCT, NCOLCMC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    REAL, DIMENSION( :, : ), intent( in ) :: CVN
    REAL, DIMENSION( : ), intent( inout ) :: CVWEIGHT
    REAL, DIMENSION( :, : ), intent( in ) :: N, NLX, NLY, NLZ
    REAL, DIMENSION( :,: ), intent( in ) :: X
    REAL, DIMENSION( : ), intent( inout ) :: MASS_ELE
    INTEGER, DIMENSION( : ), intent( in ) :: FINDM
    INTEGER, DIMENSION( : ), intent( in ) :: COLM
    INTEGER, DIMENSION( : ), intent( in ) :: MIDM

    REAL, DIMENSION( : ), intent( inout ) :: MASS_MN_PRES
    INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
    ! Local variables
    LOGICAL :: D1, D3, DCYL
    REAL, DIMENSION( : ), allocatable :: DETWEI, RA
    REAL, DIMENSION( :, : ), allocatable :: NX, NY, NZ
    REAL :: VOLUME, NN, NM, MN, MM
    INTEGER :: ELE, CV_ILOC, CV_JLOC, CV_NODI, CV_NODJ, CV_GI, COUNT, IT, idim,&
         max_iterations

    type(scalar_field) :: cv_mass
    type(tensor_field), dimension(size(psi)) :: fempsi_rhs
    type(vector_field), dimension(size(psi_ave)) :: psi_ave_temp
    type(vector_field), dimension(size(psi_int)) :: psi_int_temp
    type(tensor_field), pointer :: tfield
    type(csr_matrix) :: mat
    type(petsc_csr_matrix) :: pmat
    type(csr_sparsity), pointer :: sparsity

    ewrite(3,*) 'In PROJ_CV_TO_FEM_state'

    tfield=>psi(1)%ptr
    call allocate(cv_mass,psi(1)%ptr%mesh)
    call zero(cv_mass)


    ALLOCATE( DETWEI( CV_NGI )) 
    ALLOCATE( RA( CV_NGI ))
    ALLOCATE( NX( CV_NLOC, CV_NGI ))
    ALLOCATE( NY( CV_NLOC, CV_NGI ))
    ALLOCATE( NZ( CV_NLOC, CV_NGI ))

    do it=1,size(fempsi)
       call zero(fempsi(it)%ptr)
       call allocate(fempsi_rhs(it),psi(it)%ptr%mesh,"RHS",&
            dim=psi(it)%ptr%dim)
       call zero(fempsi_rhs(it))
       call halo_update(PSI(IT)%ptr)
    end do

    do it=1,size(psi_ave)
       call allocate(psi_ave_temp(it),psi_ave(it)%ptr%dim,&
            psi_ave(it)%ptr%mesh,&
            "PsiAveTemp")
       call set(psi_ave_temp(it),psi_ave(it)%ptr)
       call zero(psi_ave(it)%ptr)
    end do

    do it=1,size(psi_int)
       call allocate(psi_int_temp(it),psi_int(it)%ptr%dim,&
            psi_int(it)%ptr%mesh,&
            "PsiIntTemp")
       call set(psi_int_temp(it),psi_int(it)%ptr)
       call zero(psi_int(it)%ptr)
    end do

    sparsity=>extract_csr_sparsity(packed_state,"PressureMassMatrixSparsity")
    call allocate(mat,sparsity,name="ProjectionMatrix")
    call allocate(pmat,sparsity,[1,1],name="ProjectionMatrix")
    call zero(mat)

    
    IF(IGETCT.NE.0) MASS_MN_PRES=0.0

    D1 = ( NDIM == 1 )
    D3 = ( NDIM == 3 )
    DCYL = .FALSE. 

    Loop_Elements: DO ELE = 1, TOTELE
       if (isParallel()) then
          if (.not. assemble_ele(psi_int(1)%ptr,ele)) cycle
       end if

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNLXR( ELE, X(1,:), X(2,:), X(3,:), X_NDGLN, TOTELE, X_NONODS, CV_NLOC, CV_NGI, &
            N, NLX, NLY, NLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
            NX, NY, NZ ) 
       MASS_ELE( ELE ) = VOLUME

       Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC


          CV_NODI = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
          !ewrite(3,*)'ele,CV_NODI,CV_ILOC:',ele,CV_NODI,CV_ILOC, x(CV_NODI)

          if (.not.  node_owned(PSI(IT)%ptr,cv_nodi)) cycle

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

             call addto(mat,cv_nodi,cv_nodj,NN)

             call addto(CV_MASS, CV_NODI,  MM )
             


             DO IT = 1, size(fempsi_rhs)
                FEMPSI_RHS(it)%val(:,:,CV_NODI) = FEMPSI_RHS(it)%val(:,:,CV_NODI) &
                     + NM * PSI(IT)%ptr%val(:,:,CV_NODJ)
             END DO
             DO IT = 1, size(psi_ave)
                call addto( PSI_AVE(IT)%PTR, node_number=CV_NODI, val=MN * PSI_AVE_TEMP(IT)%val( : , CV_NODJ ) )
             END DO
             DO IT = 1, size(psi_int)
                call addto( PSI_INT(IT)%PTR, node_number=CV_NODI, val=MN * PSI_INT_TEMP(IT)%val( : , CV_NODJ ) )
             END DO

          END DO Loop_CV_JLOC

       END DO Loop_CV_ILOC

    END DO Loop_Elements

    call halo_update(cv_mass)
    call invert(cv_mass)

    ! Form average...
    DO IT = 1, size(psi_ave)
       call halo_update(PSI_AVE(it)%PTR)
       call scale(PSI_AVE(it)%PTR,cv_mass)
    END DO

    DO IT= 1, size(psi_int)
       call halo_update(PSI_int(it)%PTR)
    end do
    
    ! Solve...

    call get_option( trim(PSI(1)%ptr%option_path)//"/prognostic/solver/max_iterations", &
         max_iterations,  default =  500 )
    if (max_iterations ==0) THEN
       DO IT = 1, size(fempsi)
          call zero_non_owned(fempsi_rhs(it))
          CALL petsc_solve(  FEMPSI(IT)%ptr,  &
               MAT,  &
               FEMPSI_RHS(IT),  &
               option_path = "/material_phase[0]/scalar_field::Pressure/prognostic" )
!          FEMPSI(IT)%ptr%val(:,:,:)=PSI(IT)%ptr%val(:,:,:)
       END DO
    ELSE
       DO IT = 1, size(fempsi)
          call zero_non_owned(fempsi_rhs(it))
          CALL petsc_solve(  FEMPSI(IT)%ptr,  &
               MAT,  &
               FEMPSI_RHS(IT),  &
               option_path = trim(PSI(1)%ptr%option_path)//"/prognostic" )

!          FEMPSI(IT)%ptr%val(:,:,:)=PSI(IT)%ptr%val(:,:,:)
       END DO
     END IF

    

    call DEALLOCATE( cv_MASS )
    DO IT = 1, size(fempsi_RHS)
       CALL DEALLOCATE( FEMPSI_RHS(it) )
    END DO
    do it=1, size(psi_ave_temp)
       call DEALLOCATE( PSI_AVE_temp(it) )
    end do
    do it=1, size(psi_int_temp)
       call DEALLOCATE( PSI_INT_temp(it) )
    end do
    call DEALLOCATE( MAT )
    call DEALLOCATE( PMAT )
    DEALLOCATE( DETWEI )
    DEALLOCATE( RA )
    DEALLOCATE( NX )
    DEALLOCATE( NY )
    DEALLOCATE( NZ )

    ewrite(3,*) 'Leaving PROJ_CV_TO_FEM_state'

    RETURN

  END SUBROUTINE PROJ_CV_TO_FEM_state




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
      state, StorName, indx  )

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
    REAL, DIMENSION( :, :, : ), intent( in ) ::  SUF_T_BC
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
    integer, intent(inout) :: indx
    ! Local variables
    REAL, DIMENSION( :, :, : ), ALLOCATABLE :: MASELE
    REAL, DIMENSION( :, :, :, :, : ), ALLOCATABLE :: VTX_ELE, VTOLDX_ELE
    LOGICAL :: D1, D3, DCYL, APPLYBC( NCOMP, NPHASE )
    REAL, pointer, dimension( : ) :: DETWEI, RA
    REAL, pointer, DIMENSION( :,:,:):: NX_ALL
    REAL, pointer, DIMENSION( :, :, : ) :: X_NX_ALL
    REAL, pointer :: VOLUME
    REAL, DIMENSION( CV_NLOC, CV_NLOC )  :: MASS, INV_MASS
    REAL, DIMENSION( NDIM, X_SNLOC ) :: XSL( 3, X_SNLOC ), SNORMXN( 3, SBCVNGI ), SDETWE( SBCVNGI )
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
            , state,StorName , indx )

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
                VLM_NORX(:) = MATMUL( SNORMXN( 1:NDIM, : ), &
                     SDETWE(:) * SBCVFEN( CV_SILOC, : ) * SBCVFEN( CV_SJLOC, : ) )

                ! add diffusion term...
                DO IPHASE = 1, NPHASE
                   DO ICOMP = 1, NCOMP
                      IF ( APPLYBC( ICOMP, IPHASE ) ) THEN
                  
                         IF ( SELE2 /= 0 ) THEN
                            IF ( WIC_T_BC( ICOMP, IPHASE, SELE2 ) == WIC_T_BC_DIRICHLET ) THEN

                               RTBC = SUF_T_BC( ICOMP, IPHASE, CV_SJLOC2 + CV_SNLOC * ( SELE2 - 1 ) )

                               VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                              - VLM_NORX(:) * 0.5 * ( FEMT( ICOMP, IPHASE, CV_NODJ ) - RTBC  )

                               VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                              - VLM_NORX(:) * 0.5 * ( FEMTOLD( ICOMP, IPHASE, CV_NODJ ) - RTBC  )

                            END IF
                         ELSE
                            VTX_ELE(:, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                              VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                              - VLM_NORX(:) * 0.5 * ( FEMT( ICOMP, IPHASE, CV_NODJ ) - FEMT( ICOMP, IPHASE, CV_NODJ2 )  )
                            VTOLDX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) = &
                              VTX_ELE( :, ICOMP, IPHASE, CV_ILOC, ELE ) &
                              - VLM_NORX(:) * 0.5 * ( FEMTOLD( ICOMP, IPHASE, CV_NODJ ) - FEMTOLD( ICOMP, IPHASE, CV_NODJ2 )  )

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


  SUBROUTINE DG_DERIVS_ALL2( FEMT, FEMTOLD, &
       DTX_ELE, DTOLDX_ELE, &
       NDIM, NPHASE, CV_NONODS, TOTELE, CV_NDGLN, &
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

    INTEGER, intent( in ) :: NDIM, NPHASE, CV_NONODS, TOTELE, X_NLOC, CV_NGI, CV_NLOC, &
         &                   X_NONODS, STOTEL, CV_SNLOC, X_SNLOC, SBCVNGI, NFACE
    REAL, DIMENSION( :, : ), intent( in ) :: FEMT, FEMTOLD
    REAL, DIMENSION( :, :, :, : ), intent( inout ) :: DTX_ELE, DTOLDX_ELE
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) ::  XCV_NDGLN
    INTEGER, DIMENSION( :,  :, : ), intent( in ) ::  WIC_T_BC
    REAL, DIMENSION( :, :, : ), intent( in ) ::  SUF_T_BC
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
    !integer, dimension(:), intent(inout) :: StorageIndexes
    integer, intent(inout) :: StorageIndexes

    ! Local variables
    REAL, DIMENSION( :, :, : ), ALLOCATABLE :: MASELE
    REAL, DIMENSION( :, :, :, : ), ALLOCATABLE :: VTX_ELE, VTOLDX_ELE
    LOGICAL :: D1, D3, DCYL, APPLYBC( NPHASE )
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
         X_INOD, SGI, X_SILOC, X_ILOC, IDIM

    ewrite(3,*)'in DG_DERIVS'

    DTX_ELE = 0.0 ; DTOLDX_ELE = 0.0

    ALLOCATE( MASELE( CV_NLOC, CV_NLOC, TOTELE ) )
    ALLOCATE( VTX_ELE( NDIM, NPHASE, CV_NLOC, TOTELE ) )
    ALLOCATE( VTOLDX_ELE( NDIM, NPHASE, CV_NLOC, TOTELE ) )

    MASELE = 0.0 ; VTX_ELE = 0.0 ; VTOLDX_ELE = 0.0

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
            , state,StorName , StorageIndexes )

       Loop_CV_ILOC: DO CV_ILOC = 1, CV_NLOC

          CV_NODI = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC )

          Loop_CV_JLOC: DO CV_JLOC = 1, CV_NLOC

             CV_NODJ = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_JLOC )

             NN  = SUM( N( CV_ILOC, : ) * N(  CV_JLOC, : ) * DETWEI )
             NNX = MATMUL( NX_ALL( :, CV_JLOC, : ), N( CV_ILOC, : )  * DETWEI )

             MASELE( CV_ILOC, CV_JLOC, ELE) = MASELE( CV_ILOC, CV_JLOC, ELE ) + NN

             DO IPHASE = 1, NPHASE
                VTX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                     VTX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                     + NNX (:) * FEMT( IPHASE, CV_NODJ )

                VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                     VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                     + NNX (:) * FEMTOLD( IPHASE, CV_NODJ )
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
             APPLYBC = ( WIC_T_BC( 1, :, SELE2 ) == WIC_T_BC_DIRICHLET )
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
                   IF ( APPLYBC( IPHASE ) ) THEN

                      VTX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                           VTX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                           - VLM_NORX(:) * 0.5 * ( FEMT( IPHASE, CV_NODJ ) - FEMT( IPHASE, CV_NODJ2 ) * NRBC )
                      VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) = &
                           VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                           - VLM_NORX(:) * 0.5 * ( FEMTOLD( IPHASE, CV_NODJ ) - FEMTOLD( IPHASE, CV_NODJ2 ) * NRBC )

                      IF ( SELE2 /= 0 ) THEN
                         IF ( WIC_T_BC(1,  IPHASE, SELE2 ) == WIC_T_BC_DIRICHLET ) THEN

                            RTBC = SUF_T_BC( 1,IPHASE, CV_SJLOC2+ CV_SNLOC*( SELE2-1) )

                            VTX_ELE( :, IPHASE, CV_ILOC, ELE ) = VTX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                                 + VLM_NORX(:) * 0.5 * RTBC
                            VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) = VTOLDX_ELE( :, IPHASE, CV_ILOC, ELE ) &
                                 + VLM_NORX(:) * 0.5 * RTBC

                         END IF
                      END IF
                   END IF
                END DO
             END DO
          END DO

       END DO Between_Elements_And_Boundary

    END DO Loop_Elements2


    Loop_Elements3: DO ELE = 1, TOTELE

       MASS( :, : ) = MASELE( :, :, ELE )
       CALL MATDMATINV( MASS, INV_MASS, CV_NLOC )

       FORALL ( IDIM = 1:NDIM, IPHASE = 1:NPHASE )

          DTX_ELE( IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ), VTX_ELE( IDIM, IPHASE, :, ELE ) )
          DTOLDX_ELE( IDIM, IPHASE, :, ELE ) = MATMUL( INV_MASS( :, : ) , VTOLDX_ELE( IDIM, IPHASE, :, ELE ) )

          !DTX_ELE( IDIM, :, IPHASE, ELE ) = MATMUL( INV_MASS( :, : ), VTX_ELE( IDIM, IPHASE, :, ELE ) )
          !DTOLDX_ELE( IDIM, :, IPHASE, ELE ) = MATMUL( INV_MASS( :, : ) , VTOLDX_ELE( IDIM, IPHASE, :, ELE ) )
       END FORALL

    END DO Loop_Elements3

    DEALLOCATE( MASELE, VTX_ELE, VTOLDX_ELE )

    ewrite(3,*)'about to leave DG_DERIVS'

    RETURN

  END SUBROUTINE DG_DERIVS_ALL2

 SUBROUTINE ONVDLIM_ALL( TOTELE, &
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

        SUBROUTINE ONVDLIM_ANO( TOTELE, &
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

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD

!     if((PELE==433).and.(PELEOT==435)) then
!          print *,'***INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT:',INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT
!          print *,'***FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM:',FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM
!          print *,'***NOLIMI, FIRORD, TUPWIN2, TUPWI22:', NOLIMI, FIRORD, TUPWIN2, TUPWI22
!     endif

    RETURN

  END SUBROUTINE ONVDLIM_ANO


       SUBROUTINE ONVDLIM_2nd( TOTELE, &
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



 SUBROUTINE ONVDLIM( TOTELE, &
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






   SUBROUTINE ONVDLIMsqrt( TOTELE, &
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
         ! MAXUF = MAX( UC, UF )
          TILDEUF = MIN( 1.0, XI * UC, MAXUF )
          !TILDEUF = MIN( 1.0, 3.0 * UC, MAXUF )
       ENDIF

    ELSE ! Outside the region 0<UC<1 on the NVD, use first-order upwinding
       TILDEUF = UC
    ENDIF

    NVDFUNNEW = TILDEUF

  end function nvdfunnew





   REAL FUNCTION NVDFUNNEW_sqrt( UF, UC, COURAT )
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
         ! MAXUF = MAX( UC, UF )
          TILDEUF = MIN( 1.0, XI * UC, MAXUF )
          !TILDEUF = MIN( 1.0, 3.0 * UC, MAXUF )
       ENDIF

    ELSE ! Outside the region 0<UC<1 on the NVD, use first-order upwinding
       TILDEUF = UC
    ENDIF

    NVDFUNNEW_sqrt = TILDEUF

  end function nvdfunnew_sqrt






  FUNCTION NVDFUNNEW_MANY( UF, UC, XI_LIMIT ) result(nvd_limit) 
    implicit none
    ! The function computes NVDFUNNEW, the normalised value of the 
    ! advected variable on the face of the control volume, based on 
    ! the normalised value of the advected variable in the donor CV,
    ! UC, and the high-order estimate of the face value UF. 
    ! NVDFUNNEW is limited so that it is in the non-oscillatory 
    ! region of normalised variable diagram (NVD).
    !
    ! XI is the parameter in equation 38 of the Riemann paper. If XI is equal
    ! to 2 then this corresponds to a TVD condition in 1-D, a value of XI 
    ! equal to 3 has been recommended elsewhere
    ! 
    REAL, DIMENSION( : ), intent(in)  :: UC, UF, XI_LIMIT
    real, dimension(size(uc)) :: nvd_limit


    nvd_limit= MAX(  MIN(UF, XI_LIMIT*UC, 1.0), UC) 

  end function nvdfunnew_many




  FUNCTION NVDFUNNEW_MANY_sqrt( UF, UC, XI_LIMIT ) result(nvd_limit) 
    implicit none
    ! The function computes NVDFUNNEW, the normalised value of the 
    ! advected variable on the face of the control volume, based on 
    ! the normalised value of the advected variable in the donor CV,
    ! UC, and the high-order estimate of the face value UF. 
    ! NVDFUNNEW is limited so that it is in the non-oscillatory 
    ! region of normalised variable diagram (NVD).
    !
    ! XI is the parameter in equation 38 of the Riemann paper. If XI is equal
    ! to 2 then this corresponds to a TVD condition in 1-D, a value of XI 
    ! equal to 3 has been recommended elsewhere
    ! 
    REAL, DIMENSION( : ), intent(in)  :: UC, UF, XI_LIMIT
    real, dimension(size(uc)) :: nvd_limit


    nvd_limit= MAX(  MIN(UF, XI_LIMIT*UC, 1.0), UC) 

  end function nvdfunnew_many_sqrt





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




  SUBROUTINE SCVDETNX_new( ELE,      GI,        &
                                !     - INTEGERS
       NLOC,     SVNGI,   TOTELE, NDIM,  &
       XNDGLN,   XNONOD,&
                                !     - REALS
       CVDETWEI, CVNORMX, CVNORMy,&
       CVNORMz,SVN,     SVNLX,    &
       SVNLY,    SVWEIGH, XC_ALL,     &
       X_ALL,        &
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
    INTEGER , intent( in ) ::  SVNGI,  TOTELE, NDIM
    INTEGER, intent( in ) ::   XNONOD     
    REAL, DIMENSION( NDIM ), intent( in ) ::   XC_ALL
    INTEGER, DIMENSION( : ), intent( in ) :: XNDGLN
    REAL, DIMENSION( SVNGI ), intent( inout ) :: CVNORMX, CVNORMy, CVNORMz
    REAL, DIMENSION( : ), intent( inout ) :: CVDETWEI
    REAL, DIMENSION( :, : ), intent( in ) :: SVN, SVNLX, SVNLY
    REAL, DIMENSION( : ), intent( in ) :: SVWEIGH
    REAL, DIMENSION( :, : ), intent( in ) :: X_ALL ! dimension(NDIM,XNONOD)
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
    REAL :: RGI, RDUM

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

          DXDLX = DXDLX + SVNLX(JLOC,GI)*X_ALL(1,NODJ)
          DXDLY = DXDLY + SVNLY(JLOC,GI)*X_ALL(1,NODJ) 
          DYDLX = DYDLX + SVNLX(JLOC,GI)*X_ALL(2,NODJ) 
          DYDLY = DYDLY + SVNLY(JLOC,GI)*X_ALL(2,NODJ) 
          DZDLX = DZDLX + SVNLX(JLOC,GI)*X_ALL(3,NODJ) 
          DZDLY = DZDLY + SVNLY(JLOC,GI)*X_ALL(3,NODJ) 

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X_ALL(1,NODJ)
          POSVGIY = POSVGIY + SVN(JLOC,GI)*X_ALL(2,NODJ)
          POSVGIZ = POSVGIZ + SVN(JLOC,GI)*X_ALL(3,NODJ)
       end do

       !     - Note that POSVGIX,POSVGIY and POSVGIZ can be considered as the 
       !     - components of the Gauss pnt GI with the co-ordinate origin 
       !     - positioned at the current control volume NODI.

       POSVGIX = POSVGIX - XC_ALL(1)
       POSVGIY = POSVGIY - XC_ALL(2)
       POSVGIZ = POSVGIZ - XC_ALL(3)

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
!       CALL NORMGI( CVNORMX_ALL(1,GI), CVNORMX_ALL(2,GI), CVNORMX_ALL(3,GI),&
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

          DXDLX = DXDLX + SVNLX(JLOC,GI)*X_ALL(1,NODJ) 
          DYDLX = DYDLX + SVNLX(JLOC,GI)*X_ALL(2,NODJ) 

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X_ALL(1,NODJ)
          POSVGIY = POSVGIY + SVN(JLOC,GI)*X_ALL(2,NODJ)

          RGI = RGI + SVN(JLOC,GI)*X_ALL(2,NODJ)

       end do ! Was loop 300
       !     
       !     - Note that POSVGIX and POSVGIY can be considered as the components 
       !     - of the Gauss pnt GI with the co-ordinate origin positioned at the
       !     - current control volume NODI.
       !     

       POSVGIX = POSVGIX - XC_ALL(1)
       POSVGIY = POSVGIY - XC_ALL(2)

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
!       CALL NORMGI( CVNORMX_ALL(1,GI), CVNORMX_ALL(2,GI), RDUM,&
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

          POSVGIX = POSVGIX + SVN(JLOC,GI)*X_ALL(1,NODJ)

       end do ! Was loop 300
       !     
       !     - Note that POSVGIX and POSVGIY can be considered as the components 
       !     - of the Gauss pnt GI with the co-ordinate origin positioned at the
       !     - current control volume NODI.
       !     
       !          EWRITE(3,*)'POSVGIX, XC,POSVGIX - XC:',POSVGIX, XC,POSVGIX - XC
       POSVGIX = POSVGIX - XC_ALL(1)
! SIGN(A,B) sign of B times A. 
       CVNORMX(GI) = SIGN( 1.0, POSVGIX )
!       CVNORMX_ALL(1,GI) = SIGN( 1.0, POSVGIX )

!       IF(POSVGIX > 0 ) THEN
!          CVNORMX_ALL(1,GI) = +1.0
!       ELSE
!          CVNORMX_ALL(1,GI) = -1.0
!       ENDIF
!       CVNORMY(GI)=0.0 
!       CVNORMZ(GI)=0.0

       DETJ = 1.0
       CVDETWEI(GI)  = DETJ*SVWEIGH(GI)

       !       IF(GI.EQ.3) THEN
       !         EWRITE(3,*)'CVNORMX(GI),POSVGIX,XC:',CVNORMX(GI),POSVGIX,XC
       !         STOP 39344
       !       ENDIF

    ENDIF Conditional_Dimension

  END SUBROUTINE SCVDETNX_new





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



  real function ptolfun(value)
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




   real function tolfun(value)
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




  function tolfun_many(val) result(v_tolfun)

    implicit none
    real, dimension(:), intent(in) :: val
    real, dimension(size(val)) :: v_tolfun
    ! Local
    real, parameter :: tolerance = 1.e-10

    v_tolfun = sign( 1.0, val ) * max( tolerance, abs(val) )

    return

  end function tolfun_many







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
       SMATFEN, SCVFEN, SCVNGI, GI, NDIM, TDIFFUSION, DIFF_GI_ADDED, &
       HDC, &
       T_CV_NODJ, T_CV_NODI, &
       TOLD_CV_NODJ, TOLD_CV_NODI, &
       ELE, ELE2, CVNORMX_ALL,  &
       LOC_DTX_ELE_ALL, LOC_DTOLDX_ELE_ALL, LOC2_DTX_ELE_ALL, LOC2_DTOLDX_ELE_ALL, &
       SELE, STOTEL, LOC_WIC_T_BC, CV_OTHER_LOC, MAT_OTHER_LOC )
    ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
    ! based on a non-linear method and a non-oscillating scheme.
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NLOC, MAT_NLOC, CV_NONODS,NPHASE, TOTELE, MAT_NONODS, &
         &                   SCVNGI, GI, NDIM, ELE, ELE2, &
         &                   SELE, STOTEL
    REAL, intent( in ) :: HDC
    REAL, DIMENSION( NPHASE ), intent( in ) :: T_CV_NODJ, T_CV_NODI, &
         &                                     TOLD_CV_NODJ, TOLD_CV_NODI
    REAL, DIMENSION( NPHASE ), intent( inout ) :: DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
    INTEGER, DIMENSION( : ), intent( in ) :: MAT_NDGLN
    INTEGER, DIMENSION( NPHASE ), intent( in ) :: LOC_WIC_T_BC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: MAT_OTHER_LOC
    REAL, DIMENSION( :, : ), intent( in ) :: SMATFEN
    REAL, DIMENSION( :, : ), intent( in ) :: SCVFEN
    REAL, DIMENSION( :, :, :, : ), intent( in ) :: TDIFFUSION
    REAL, DIMENSION( :, :, : ), intent( in ) :: DIFF_GI_ADDED
    REAL, DIMENSION( NDIM, NPHASE, CV_NLOC ), intent( in ) :: LOC_DTX_ELE_ALL, LOC_DTOLDX_ELE_ALL, LOC2_DTX_ELE_ALL, LOC2_DTOLDX_ELE_ALL
    REAL, DIMENSION( : ), intent( in ) :: CVNORMX_ALL

    ! local variables

    ! DIFF_MIN_FRAC is the fraction of the standard diffusion coefficient to use 
    ! in the non-linear diffusion scheme. DIFF_MAX_FRAC is the maximum fraction. 
    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.05, DIFF_MAX_FRAC = 20.0

    REAL :: COEF
    INTEGER :: CV_KLOC, CV_KLOC2, MAT_KLOC, MAT_KLOC2, MAT_NODK, MAT_NODK2, CV_ILOC, IPHASE
    LOGICAL :: ZER_DIFF

    REAL, DIMENSION ( :, : ), allocatable :: DTDX_GI_ALL, DTOLDDX_GI_ALL, DTDX_GI2_ALL, DTOLDDX_GI2_ALL
    REAL, DIMENSION ( : ), allocatable :: N_DOT_DKDT_ALL, N_DOT_DKDTOLD_ALL, N_DOT_DKDT2_ALL, N_DOT_DKDTOLD2_ALL
    REAL, DIMENSION ( : ), allocatable :: DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL
    REAL, DIMENSION ( :, :, : ), allocatable :: DIFF_GI, DIFF_GI2

    ALLOCATE( DTDX_GI_ALL( NDIM, NPHASE ), DTOLDDX_GI_ALL( NDIM, NPHASE ) )
    ALLOCATE( DTDX_GI2_ALL( NDIM, NPHASE ), DTOLDDX_GI2_ALL( NDIM, NPHASE ) )
    ALLOCATE( N_DOT_DKDT_ALL( NPHASE ), N_DOT_DKDTOLD_ALL( NPHASE ) )
    ALLOCATE( N_DOT_DKDT2_ALL( NPHASE ), N_DOT_DKDTOLD2_ALL( NPHASE ) )
    ALLOCATE( DIFF_STAND_DIVDX_ALL( NPHASE ), DIFF_STAND_DIVDX2_ALL( NPHASE ) )

    ALLOCATE( DIFF_GI( NDIM, NDIM, NPHASE ) ) 
    ALLOCATE( DIFF_GI2( NDIM, NDIM, NPHASE ) ) 

    ZER_DIFF = .FALSE.
    IF ( SELE /= 0 ) ZER_DIFF = ANY ( LOC_WIC_T_BC( : ) /= WIC_T_BC_DIRICHLET )

    Cond_ZerDiff: IF ( ZER_DIFF ) THEN

       DIFF_COEF_DIVDX = 0.0 
       DIFF_COEFOLD_DIVDX = 0.0

    ELSE

       DTDX_GI_ALL = 0.0 ; DTOLDDX_GI_ALL = 0.0
       DO CV_KLOC = 1, CV_NLOC
          DTDX_GI_ALL( :, : ) = DTDX_GI_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) * LOC_DTX_ELE_ALL( :, :, CV_KLOC )
          DTOLDDX_GI_ALL( :, : ) = DTOLDDX_GI_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) * LOC_DTOLDX_ELE_ALL( :, :, CV_KLOC )
       END DO

       DIFF_GI = 0.0
       DO MAT_KLOC = 1, MAT_NLOC
          MAT_NODK = MAT_NDGLN( ( ELE - 1 ) * MAT_NLOC + MAT_KLOC )
          DO IPHASE = 1, NPHASE
             DIFF_GI( :, :, IPHASE ) = DIFF_GI( :, :, IPHASE ) &
                  + SMATFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK, :, :, IPHASE )
          END DO
       END DO
       DIFF_GI = DIFF_GI + DIFF_GI_ADDED
       DIFF_GI = MAX( 0.0, DIFF_GI ) 

       DO IPHASE = 1, NPHASE
          N_DOT_DKDT_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI( :, :, IPHASE ), DTDX_GI_ALL( :, IPHASE ) ) )
          N_DOT_DKDTOLD_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI( :, :, IPHASE ), DTOLDDX_GI_ALL( :, IPHASE ) ) )

          COEF = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI( :, :, IPHASE ), CVNORMX_ALL ) )
          DIFF_STAND_DIVDX_ALL( IPHASE ) = COEF / HDC
       END DO



       Conditional_MAT_DISOPT_ELE2: IF ( ( ELE2 /= 0 ) .AND. ( ELE2 /= ELE ) ) THEN


          DTDX_GI2_ALL = 0.0 ; DTOLDDX_GI2_ALL = 0.0
          DO CV_KLOC = 1, CV_NLOC
             CV_KLOC2 = CV_OTHER_LOC( CV_KLOC )
             IF ( CV_KLOC2 /= 0 )THEN
                DTDX_GI2_ALL( :, : ) = DTDX_GI2_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) *  LOC2_DTX_ELE_ALL( :, :, CV_KLOC )
                DTOLDDX_GI2_ALL( :, : ) = DTOLDDX_GI2_ALL( :, : ) + SCVFEN( CV_KLOC, GI ) * LOC2_DTOLDX_ELE_ALL( :, :, CV_KLOC )
             END IF
          END DO

          DIFF_GI2 = 0.0
          DO MAT_KLOC = 1, MAT_NLOC
             MAT_KLOC2 = MAT_OTHER_LOC( MAT_KLOC )
             IF ( MAT_KLOC2 /= 0 ) THEN
                MAT_NODK2 = MAT_NDGLN( ( ELE2 - 1 ) * MAT_NLOC + MAT_KLOC2 )
                DO IPHASE = 1, NPHASE
                   DIFF_GI2( :, :, IPHASE ) = DIFF_GI2( :, :, IPHASE ) &
                        + SMATFEN( MAT_KLOC, GI ) * TDIFFUSION( MAT_NODK2, :, :, IPHASE )
                END DO
             END IF
          END DO

          DO IPHASE = 1, NPHASE
             N_DOT_DKDT2_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), DTDX_GI2_ALL( :, IPHASE ) ) )
             N_DOT_DKDTOLD2_ALL( IPHASE ) = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), DTOLDDX_GI2_ALL( :, IPHASE ) ) )

             COEF = DOT_PRODUCT( CVNORMX_ALL, MATMUL( DIFF_GI2( :, :, IPHASE ), CVNORMX_ALL ) )
             DIFF_STAND_DIVDX2_ALL( IPHASE ) = COEF  /HDC
          END DO

          N_DOT_DKDT_ALL = 0.5 * ( N_DOT_DKDT_ALL + N_DOT_DKDT2_ALL )
          N_DOT_DKDTOLD_ALL = 0.5 * ( N_DOT_DKDTOLD_ALL + N_DOT_DKDTOLD2_ALL )

          DIFF_STAND_DIVDX_ALL = MIN( DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL ) 

       END IF Conditional_MAT_DISOPT_ELE2


       DIFF_COEF_DIVDX = MAX( DIFF_MIN_FRAC * DIFF_STAND_DIVDX_ALL, N_DOT_DKDT_ALL / TOLFUN_MANY( T_CV_NODJ - T_CV_NODI ) )
       DIFF_COEFOLD_DIVDX = MAX( DIFF_MIN_FRAC * DIFF_STAND_DIVDX_ALL, N_DOT_DKDTOLD_ALL / TOLFUN_MANY( TOLD_CV_NODJ - TOLD_CV_NODI ) )

       DIFF_COEF_DIVDX = MIN( DIFF_MAX_FRAC * DIFF_STAND_DIVDX_ALL, DIFF_COEF_DIVDX )
       DIFF_COEFOLD_DIVDX = MIN( DIFF_MAX_FRAC * DIFF_STAND_DIVDX_ALL, DIFF_COEFOLD_DIVDX )

    END IF Cond_ZerDiff

    IF ( SELE /= 0 ) THEN
       DO IPHASE=1,NPHASE
          IF( LOC_WIC_T_BC( IPHASE ) /= WIC_T_BC_DIRICHLET ) THEN
             DIFF_COEF_DIVDX( IPHASE ) = 0.0 
             DIFF_COEFOLD_DIVDX( IPHASE ) = 0.0
          ENDIF
       END DO
    ENDIF


    DEALLOCATE( DIFF_GI, DIFF_GI2, &
         DTDX_GI_ALL, DTOLDDX_GI_ALL, &
         DTDX_GI2_ALL, DTOLDDX_GI2_ALL, &
         N_DOT_DKDT_ALL, N_DOT_DKDTOLD_ALL, &
         N_DOT_DKDT2_ALL, N_DOT_DKDTOLD2_ALL, &
         DIFF_STAND_DIVDX_ALL, DIFF_STAND_DIVDX2_ALL )

    RETURN            

  END SUBROUTINE DIFFUS_CAL_COEFF




  SUBROUTINE DIFFUS_CAL_COEFF_STRESS_OR_TENSOR( DIFF_COEF_DIVDX, &
       DIFF_COEFOLD_DIVDX, STRESS_FORM, STRESS_FORM_STAB, ZERO_OR_TWO_THIRDS, &
       U_SNLOC, U_NLOC, CV_SNLOC, CV_NLOC, MAT_NLOC, NPHASE,  &
       SBUFEN,SBCVFEN,SBCVNGI, NDIM_VEL, NDIM, SLOC_UDIFFUSION, SLOC2_UDIFFUSION, DIFF_GI_ADDED, &
       HDC, &
       U_CV_NODJ_IPHA_ALL, U_CV_NODI_IPHA_ALL, &
       UOLD_CV_NODJ_IPHA_ALL, UOLD_CV_NODI_IPHA_ALL, &
       ELE, ELE2, SNORMXN_ALL, &
       SLOC_DUX_ELE_ALL, SLOC2_DUX_ELE_ALL,   SLOC_DUOLDX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL,  &
       SELE, STOTEL, WIC_U_BC, WIC_U_BC_DIRICHLET, SIMPLE_DIFF_CALC, DIFF_MIN_FRAC, DIFF_MAX_FRAC  )
    ! This sub calculates the effective diffusion coefficientd DIFF_COEF_DIVDX,DIFF_COEFOLD_DIVDX
    ! based on a non-linear method and a non-oscillating scheme.
! This implements the stress and tensor form of diffusion and calculates a jump conidition. 
! which is in DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
! The coefficient are in N_DOT_DKDU, N_DOT_DKDUOLD. 
! look at the manual DG treatment of viscocity. 
    IMPLICIT NONE
    LOGICAL, intent( in ) :: STRESS_FORM, STRESS_FORM_STAB, SIMPLE_DIFF_CALC
    INTEGER, intent( in ) :: U_SNLOC, U_NLOC, CV_SNLOC,CV_NLOC, MAT_NLOC, NPHASE,  &
         &                   SBCVNGI, NDIM_VEL, NDIM, ELE, ELE2, &
         &                   SELE, STOTEL, WIC_U_BC_DIRICHLET
    REAL, intent( in ) :: HDC, DIFF_MIN_FRAC, DIFF_MAX_FRAC
    REAL, DIMENSION(NDIM_VEL,NPHASE,SBCVNGI), intent( in ) :: U_CV_NODJ_IPHA_ALL, U_CV_NODI_IPHA_ALL, &
                                                          UOLD_CV_NODJ_IPHA_ALL, UOLD_CV_NODI_IPHA_ALL
    REAL, intent( in ) :: ZERO_OR_TWO_THIRDS
    REAL, DIMENSION( NDIM,NPHASE,SBCVNGI ), intent( inout ) :: DIFF_COEF_DIVDX, DIFF_COEFOLD_DIVDX
    INTEGER, DIMENSION( NDIM,NPHASE,STOTEL ), intent( in ) ::WIC_U_BC
    REAL, DIMENSION( CV_SNLOC, SBCVNGI  ), intent( in ) :: SBCVFEN
    REAL, DIMENSION( U_SNLOC, SBCVNGI  ), intent( in ) :: SBUFEN
    REAL, DIMENSION( NDIM,NDIM,NPHASE,CV_SNLOC ), intent( in ) :: SLOC_UDIFFUSION, SLOC2_UDIFFUSION
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
!    LOGICAL, PARAMETER :: SIMPLE_DIFF_CALC2 = .false.
    !REAL, PARAMETER :: DIFF_MIN_FRAC = 0.005, DIFF_MAX_FRAC = 200.0
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.1, DIFF_MAX_FRAC = 1000.0 ! works well but oscillations
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.5, DIFF_MAX_FRAC = 1000000.0 ! works well no oscillations
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.25, DIFF_MAX_FRAC = 100.0 ! works well no oscillations
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.01, DIFF_MAX_FRAC = 100.0 ! works well no oscillations
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.2, DIFF_MAX_FRAC = 100.0 ! works well no oscillations  ****recommended*****
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.05, DIFF_MAX_FRAC = 200.0 ! works well no oscillations
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.1, DIFF_MAX_FRAC = 200.0 ! works well no oscillations
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.25, DIFF_MAX_FRAC = 10000000.0 ! works well no oscillations
!    REAL, PARAMETER :: DIFF_MIN_FRAC = 0.25, DIFF_MAX_FRAC = 1000.0

    REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI, DIFF_GI2, DIFF_GI_BOTH

    REAL, DIMENSION( :, :, : ), allocatable :: N_DOT_DKDU, N_DOT_DKDUOLD, N_DOT_DKDU2, N_DOT_DKDUOLD2
    REAL, DIMENSION( :, :, : ), allocatable :: DIFF_STAND_DIVDX_U, DIFF_STAND_DIVDX2_U, &
             DIFF_COEF_DIVDX_U, DIFF_COEFOLD_DIVDX_U
    REAL, DIMENSION( :, : ), allocatable :: IDENT, RZER_DIFF_ALL
    REAL :: COEF
    INTEGER :: CV_KLOC,CV_KLOC2,MAT_KLOC,MAT_KLOC2,MAT_NODK,MAT_NODK2,IDIM,JDIM,CV_SKLOC
    INTEGER :: SGI,IPHASE
    LOGICAL :: ZER_DIFF

!    SIMPLE_DIFF_CALC=SIMPLE_DIFF_CALC2

    ALLOCATE( RZER_DIFF_ALL(NDIM,NPHASE) )


    ZER_DIFF=.FALSE.
    RZER_DIFF_ALL=1.0
    IF(SELE /= 0) THEN
       ZER_DIFF=.TRUE.
       RZER_DIFF_ALL=0.0
       DO IPHASE=1,NPHASE
       DO IDIM = 1, NDIM
          IF(WIC_U_BC( IDIM, IPHASE, SELE) == WIC_U_BC_DIRICHLET) THEN
             ZER_DIFF=.FALSE.
             RZER_DIFF_ALL(IDIM,IPHASE)=1.0
          ENDIF
       END DO
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
          ALLOCATE( DIFF_GI_BOTH(NDIM,NDIM,NPHASE,SBCVNGI) )

          ALLOCATE( IDENT(NDIM,NDIM) )
          IDENT=0.0
          DO IDIM=1,NDIM
            IDENT(IDIM,IDIM)=1.0
          END DO

          DIFF_GI = 0.0
          DO CV_SKLOC = 1, CV_SNLOC
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   DIFF_GI( 1:NDIM , 1:NDIM, IPHASE, SGI ) = DIFF_GI( 1:NDIM , 1:NDIM, IPHASE, SGI ) &
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
                      DIFF_GI2( 1:NDIM, 1:NDIM, IPHASE, SGI )= DIFF_GI2( 1:NDIM, 1:NDIM, IPHASE, SGI ) +SBCVFEN(CV_SKLOC,SGI) &
                        *SLOC2_UDIFFUSION(1:NDIM, 1:NDIM ,IPHASE, CV_SKLOC)
                   END DO
                END DO
             END DO
             DIFF_GI2=MAX(0.0, DIFF_GI2) 
             DIFF_GI=0.5*(DIFF_GI+DIFF_GI2)
          ENDIF Conditional_MAT_DISOPT_ELE2_2

          IF(STRESS_FORM) THEN
             
             IF(STRESS_FORM_STAB) THEN

                DIFF_GI_BOTH = DIFF_GI
                DO IDIM=1,NDIM
                   DO JDIM=1,NDIM
                      DIFF_GI_BOTH(IDIM, JDIM, :, :) = DIFF_GI_BOTH(IDIM, JDIM, :, :) &
                       + SQRT( DIFF_GI_ADDED(IDIM, 1,1, :, :) * DIFF_GI_ADDED(JDIM, 1,1, :, :) )
                   END DO
                END DO
                
                DO SGI=1,SBCVNGI
                   DO IPHASE=1, NPHASE
                      DO IDIM=1, NDIM_VEL
                         DIFF_COEF_DIVDX(IDIM,IPHASE,SGI)=8.* SUM( (1.+IDENT(IDIM,:))*SNORMXN_ALL(:,SGI)**2*DIFF_GI_BOTH(IDIM,:,IPHASE,SGI) ) /HDC
                      END DO
                   END DO
                END DO
             ELSE
                DO SGI=1,SBCVNGI
                   DO IPHASE=1, NPHASE
                      DO IDIM=1, NDIM_VEL
                         DIFF_COEF_DIVDX(IDIM,IPHASE,SGI)=8.*( SUM( (1.+IDENT(IDIM,:))*SNORMXN_ALL(:,SGI)**2*DIFF_GI(IDIM,:,IPHASE,SGI) ) &
                           +DIFF_GI_ADDED(IDIM, 1,1, IPHASE,SGI) ) /HDC
!                           +DIFF_GI_ADDED(IDIM, IPHASE,SGI,1,1) ) /HDC
                      END DO
                   END DO
                END DO
             ENDIF

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
               !  NDIM_VEL, NDIM, NPHASE, U_SNLOC, SBCVNGI, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, STRESS_FORM )
                 NDIM_VEL, NDIM, NPHASE, U_SNLOC, CV_SNLOC, SBCVNGI, SBUFEN, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, &
                 STRESS_FORM, STRESS_FORM_STAB )



          Conditional_MAT_DISOPT_ELE2: IF( ( ELE2 /= 0 ).AND.( ELE2 /= ELE) ) THEN



! Calculate DIFF_COEF_DIVDX, N_DOT_DKDU, N_DOT_DKDUOLD
             CALL FOR_TENS_DERIVS_NDOTS(DIFF_STAND_DIVDX2_U, N_DOT_DKDU2, N_DOT_DKDUOLD2,  &  
                    DIFF_GI_ADDED, SLOC2_DUX_ELE_ALL, SLOC2_DUOLDX_ELE_ALL, SLOC2_UDIFFUSION, &
                    NDIM_VEL, NDIM, NPHASE, U_SNLOC, CV_SNLOC, SBCVNGI, SBUFEN, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, &
                    STRESS_FORM, STRESS_FORM_STAB )




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
        DO IDIM=1, NDIM
           DIFF_COEF_DIVDX(IDIM,IPHASE,:)    =  RZER_DIFF_ALL(IDIM,IPHASE)*DIFF_COEF_DIVDX(IDIM,IPHASE,:)
           DIFF_COEFOLD_DIVDX(IDIM,IPHASE,:) =  RZER_DIFF_ALL(IDIM,IPHASE)*DIFF_COEFOLD_DIVDX(IDIM,IPHASE,:)
        END DO
     END DO


    RETURN            

  END SUBROUTINE DIFFUS_CAL_COEFF_STRESS_OR_TENSOR





        SUBROUTINE FOR_TENS_DERIVS_NDOTS( DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD,  &
                 DIFF_GI_ADDED, SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL, SLOC_UDIFFUSION, &
                 NDIM_VEL, NDIM, NPHASE, U_SNLOC, CV_SNLOC, SBCVNGI, SBUFEN, SBCVFEN, SNORMXN_ALL, HDC, ZERO_OR_TWO_THIRDS, &
                 STRESS_FORM, STRESS_FORM_STAB )

! Calculate DIFF_STAND_DIVDX_U, N_DOT_DKDU, N_DOT_DKDUOLD
! This implements the stress and tensor form of diffusion and calculates a jump conidition. 
! DIFF_STAND_DIVDX_U is the minimal amount of diffusion. 
! The coefficient are in N_DOT_DKDU, N_DOT_DKDUOLD. 
! look at the manual DG treatment of viscocity. 
    IMPLICIT NONE
      INTEGER, intent( in )  :: NDIM_VEL, NDIM, NPHASE, U_SNLOC, CV_SNLOC, SBCVNGI
      REAL, intent( in )  :: HDC, ZERO_OR_TWO_THIRDS
      LOGICAL, intent( in )  :: STRESS_FORM, STRESS_FORM_STAB
    REAL, DIMENSION( NDIM,NPHASE,SBCVNGI ), intent( inout ) :: DIFF_STAND_DIVDX_U
    REAL, DIMENSION( NDIM_VEL,NPHASE,SBCVNGI ), intent( inout ) :: N_DOT_DKDU, N_DOT_DKDUOLD
    ! DIFF_GI_ADDED( IDIM, :,:) is for dimension IDIM e.g IDIM=1 corresponds to U 
    ! the rest is for the diffusion tensor. 
    REAL, DIMENSION( NDIM_VEL, NDIM,NDIM, NPHASE, SBCVNGI), intent( in ) :: DIFF_GI_ADDED
    REAL, DIMENSION( NDIM_VEL, NDIM , NPHASE, U_SNLOC ), intent( in ) :: SLOC_DUX_ELE_ALL, SLOC_DUOLDX_ELE_ALL 
    REAL, DIMENSION( U_SNLOC, SBCVNGI  ), intent( in ) :: SBUFEN
    REAL, DIMENSION( CV_SNLOC, SBCVNGI  ), intent( in ) :: SBCVFEN
    REAL, DIMENSION( NDIM,NDIM,NPHASE,CV_SNLOC ), intent( in ) :: SLOC_UDIFFUSION
    REAL, DIMENSION( NDIM, SBCVNGI ), intent( in ) :: SNORMXN_ALL

    ! local variables
    REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI, STRESS_INDEX, STRESS_INDEXOLD
    REAL, DIMENSION( : , :, :, : ), allocatable :: DIFF_GI_BOTH
    REAL, DIMENSION( :, :, :, : ), allocatable :: DUDX_ALL_GI, DUOLDDX_ALL_GI
    REAL, DIMENSION( :, : ), allocatable :: IDENT
    REAL :: COEF, DIVU, DIVUOLD
    INTEGER :: U_KLOC,U_KLOC2,MAT_KLOC,MAT_KLOC2,IDIM,JDIM,IDIM_VEL,U_SKLOC,CV_SKLOC
    INTEGER :: SGI,IPHASE
    LOGICAL :: ZER_DIFF,SIMPLE_DIFF_CALC


       ALLOCATE( DIFF_GI(NDIM,NDIM,NPHASE,SBCVNGI) )
       ALLOCATE( STRESS_INDEX(NDIM,NDIM,NPHASE,SBCVNGI) )
       ALLOCATE( STRESS_INDEXOLD(NDIM,NDIM,NPHASE,SBCVNGI) )
       ALLOCATE( DIFF_GI_BOTH(NDIM,NDIM, NPHASE,SBCVNGI) )

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
                   DUDX_ALL_GI(:,:,:,SGI)    = DUDX_ALL_GI(:,:,:,SGI)    + SBUFEN(U_SKLOC,SGI) * SLOC_DUX_ELE_ALL(:,:,:,U_SKLOC)
                   DUOLDDX_ALL_GI(:,:,:,SGI) = DUOLDDX_ALL_GI(:,:,:,SGI) + SBUFEN(U_SKLOC,SGI) * SLOC_DUOLDX_ELE_ALL(:,:,:,U_SKLOC)
             END DO
          END DO

          DIFF_GI = 0.0
          DO CV_SKLOC = 1, CV_SNLOC
             DO SGI=1,SBCVNGI
                DO IPHASE=1, NPHASE
                   DIFF_GI( 1:NDIM , 1:NDIM, IPHASE,SGI ) = DIFF_GI( 1:NDIM , 1:NDIM, IPHASE,SGI ) &
                  + SBCVFEN(CV_SKLOC,SGI) * SLOC_UDIFFUSION( 1:NDIM , 1:NDIM , IPHASE, CV_SKLOC )
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

             DIFF_GI_BOTH = DIFF_GI
             IF(STRESS_FORM_STAB) THEN
                DO IDIM=1,NDIM
                   DO JDIM=1,NDIM
                      DIFF_GI_BOTH(IDIM, JDIM, :, :) = DIFF_GI_BOTH(IDIM, JDIM, :, :) &
                       + SQRT( DIFF_GI_ADDED(IDIM, 1,1, :, :) * DIFF_GI_ADDED(JDIM, 1,1, :, :) )
                   END DO
                END DO
             ELSE ! Tensor form
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
             ENDIF

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
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )  & 
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUDX_ALL_GI(:,IDIM_VEL,IPHASE,SGI) ) & 
! stress form addition...  
                                                           - ZERO_OR_TWO_THIRDS*SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI_BOTH(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*DIVU

! Stress form...
                         N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI)= N_DOT_DKDUOLD(IDIM_VEL,IPHASE,SGI) &
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUOLDDX_ALL_GI(IDIM_VEL,:,IPHASE,SGI) )  & 
                                                           + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*DUOLDDX_ALL_GI(:,IDIM_VEL,IPHASE,SGI) ) & 
! stress form addition...  
                                                           - ZERO_OR_TWO_THIRDS*SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI_BOTH(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*DIVUOLD
                                   
! This is for the minimum & max. diffusion...
                         DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI)   =  DIFF_STAND_DIVDX_U(IDIM_VEL,IPHASE,SGI) &
                                                           + (   SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*SNORMXN_ALL(:,SGI) )  & 
                                                               +  SNORMXN_ALL(IDIM_VEL,SGI)*DIFF_GI_BOTH(IDIM_VEL,IDIM_VEL,IPHASE,SGI)*SNORMXN_ALL(IDIM_VEL,SGI)     )/HDC
!                                                               + SUM( SNORMXN_ALL(:,SGI)*DIFF_GI_BOTH(IDIM_VEL,:,IPHASE,SGI)*SNORMXN_ALL(IDIM_VEL,SGI) )    )/HDC

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
!          DIFF_STAND_DIVDX_U=( 8.*DIFF_STAND_DIVDX_U )

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
!               STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) + FEN_TEN_XX(IDIM,JDIM) * UFENX_JLOC(JDIM) 
               STRESS_IJ( IDIM,JDIM ) = STRESS_IJ( IDIM,JDIM ) + FEN_TEN_XX(IDIM,JDIM) * UFENX_JLOC(IDIM) 
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
       TMIN_2ND_MC, TMAX_2ND_MC, &
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
    REAL, DIMENSION( : ), intent( inout ) :: TMAX_2ND_MC, TMIN_2ND_MC

    INTEGER, intent( in ) :: IANISOTROPIC
    INTEGER, intent( in ) :: NSMALL_COLM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT
    ! local variables
    INTEGER :: U_NLOC_LEV,U_KLOC_LEV,U_KLOC,U_NODK_IPHA, U_KLOC2, U_NODK2_IPHA

    logical, parameter :: LIMIT_USE_2ND=.false.


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


      SUBROUTINE GET_INT_VEL_OVERLAP( NPHASE, NDOTQ,INCOME, &
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




    SUBROUTINE GET_INT_VEL_ORIG( NPHASE, NDOTQ, INCOME, &
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
    UGI_COEF_ELE=0.0
    VGI_COEF_ELE=0.0
    WGI_COEF_ELE=0.0 

    UGI_COEF_ELE_ALL=0.0

    ! coefficients for this element ELE2
    UGI_COEF_ELE2=0.0
    VGI_COEF_ELE2=0.0
    WGI_COEF_ELE2=0.0 

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
    
    if ( sele /= 0 ) then
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
    end if

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




  SUBROUTINE GET_INT_T_DEN(FVT, FVT2, FVD, LIMD, LIMT, LIMT2, &
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
       TMAX_2ND_MC, T2MAX_2ND_MC, DENMAX_2ND_MC, &
       HDC, NDOTQ, DT, &
       SCVFENX, SCVFENY, SCVFENZ, CVNORMX, CVNORMY, CVNORMZ, &
       U,V,W, U_NDGLN,U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
       IANISOTROPIC, SMALL_FINDRM, SMALL_COLM, NSMALL_COLM, &
       TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT )
    !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
    IMPLICIT NONE
    ! Calculate T and DEN on the CV face at quadrature point GI.
    REAL, intent( inout ) :: FVT, FVT2, FVD, LIMD,LIMT,LIMT2, &
         LIMDT,LIMDTT2, &
         FEMDGI, FEMTGI, FEMT2GI
    REAL, intent( in ) :: INCOME,HDC,NDOTQ,DT
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

    LOGICAL, PARAMETER :: limit_use_2nd = .false.

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
!       if((cv_nodi_ipha==67).and.(cv_nodj_ipha==736)) then
!            print *,'--cv_kloc,FEMTGI, RSHAPE,FEMT( CV_NODK_IPHA ):',cv_kloc,FEMTGI, RSHAPE,FEMT( CV_NODK_IPHA )
!       endif
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

!        if( (CV_NODI_IPHA-(iphase-1)*cv_nonods==433) .and. (CV_NODJ_IPHA-(iphase-1)*cv_nonods==435) ) then
!            print *,'iphase, LIMT, FEMTGI:',iphase, LIMT, FEMTGI
!        endif

    END IF

!         print *,'before fvd,limd:',fvd,limd
    ! Amend for porosity...
    IF ( ELE2 /= 0 ) THEN 
       FVD    = 0.5 * ( VOLFRA_PORE(ELE) + VOLFRA_PORE(ELE2) ) * FVD
       LIMD   = 0.5 * ( VOLFRA_PORE(ELE) + VOLFRA_PORE(ELE2) ) * LIMD
    ELSE
       FVD    = VOLFRA_PORE(ELE) * FVD
       LIMD   = VOLFRA_PORE(ELE) * LIMD
    END IF
!         print *,'after fvd,limd,VOLFRA_PORE(ELE):',fvd,limd,VOLFRA_PORE(ELE)
!          stop  2911

    IF ( HI_ORDER_HALF ) THEN
       RELAX = MIN ( 2.*ABS(FEMTGI-0.5), 1.0 )
       LIMT = RELAX * LIMT + (1.-RELAX) * FEMTGI
    END IF

    LIMDT = LIMD * LIMT

    LIMDTT2 = LIMD * LIMT * LIMT2

!       if(sele=26) then
!           print *,'iphase, FEMDGI, LIMT, FIRSTORD, NOLIMI, CV_NODI_IPHA, CV_NODJ_IPHA
!       IF ( (CV_NODI_IPHA-(iphase-1)*cv_nonods==26).and.(CV_NODJ_IPHA-(iphase-1)*cv_nonods==26) )  then
!           print *,'iphase, FEMTGI, LIMT, FIRSTORD, NOLIMI, CV_NODI_IPHA, CV_NODJ_IPHA, income, FVT:', &
!                    iphase, FEMTGI, LIMT, FIRSTORD, NOLIMI, CV_NODI_IPHA, CV_NODJ_IPHA, income, FVT
!       endif

!    DEALLOCATE( SE_CV_NODK_IPHA, SE_CV_SNODK_IPHA )

    RETURN


    contains
      
      SUBROUTINE LIMITERS( &
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

!    if(IANISOTROPIC .ne. 1) stop 261

    IF ( IANISOTROPIC == 1 ) THEN
       IF( CV_NODJ /= CV_NODI ) THEN
          DO COUNT = SMALL_FINDRM(CV_NODJ), SMALL_FINDRM(CV_NODJ+1)-1
             IF ( SMALL_COLM(COUNT) == CV_NODI ) THEN
                TUPWIN = TUPWIND_MAT( COUNT + (IPHASE-1)*NSMALL_COLM)
                DENUPWIN = DENUPWIND_MAT( COUNT + (IPHASE-1)*NSMALL_COLM)
                IF ( IGOT_T2 == 1 ) T2UPWIN = T2UPWIND_MAT( COUNT + (IPHASE-1)*NSMALL_COLM)
                EXIT
             END IF
          END DO
          DO COUNT = SMALL_FINDRM(CV_NODI), SMALL_FINDRM(CV_NODI+1)-1
             IF ( SMALL_COLM(COUNT) == CV_NODJ ) THEN
                TUPWI2 = TUPWIND_MAT( COUNT + (IPHASE-1)*NSMALL_COLM)
                DENUPWI2 = DENUPWIND_MAT( COUNT + (IPHASE-1)*NSMALL_COLM)
                IF ( IGOT_T2 == 1 ) T2UPWI2 = T2UPWIND_MAT( COUNT + (IPHASE-1)*NSMALL_COLM)
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

!       if((cv_nodi==433).and.(cv_nodj==435)) then
!!       if((cv_nodi==67).and.(cv_nodj==736)) then
!       if((cv_nodi==26).and.(cv_nodj==26)) then
!          print *,'iphase,TUPWIN, TUPWI2,T(CV_NODI_IPHA), t(CV_NODj_IPHA), LIMT, FEMTGI,INCOME:', &
!                   iphase,TUPWIN, TUPWI2,T(CV_NODI_IPHA), t(CV_NODj_IPHA), LIMT, FEMTGI,INCOME
!       endif

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



 


  SUBROUTINE GET_INT_T_DEN_new(LIMF, &
       CV_DISOPT, CV_NONODS, NPHASE, NFIELD, CV_NODI, CV_NODJ, CV_ILOC, CV_JLOC, CV_SILOC, ELE, ELE2, GI, &
       CV_NLOC, TOTELE, CV_OTHER_LOC, SCVNGI, SCVFEN, F_INCOME, F_NDOTQ, &
       LOC_F, LOC_FEMF, SLOC_F, SLOC_FEMF, SLOC2_F, SLOC2_FEMF,  &
       SELE, CV_SNLOC,   U_SNLOC,     STOTEL, CV_SLOC2LOC, SLOC_SUF_F_BC, &
       U_SLOC2LOC, U_OTHER_LOC,  &
       SELE_LOC_WIC_F_BC,   &
       WIC_T_BC_DIRICHLET, WIC_D_BC_DIRICHLET, &
       HDC, DT, &
       SCVFENX_ALL, CVNORMX_ALL,  &
       LOC_UF,  &
       U_NLOC,U_NONODS,NDIM,SUFEN, INV_JAC, &
       FUPWIND_IN, FUPWIND_OUT, DISTCONTINUOUS_METHOD, QUAD_ELEMENTS, SHAPE_CV_SNL, DOWNWIND_EXTRAP_INDIVIDUAL, &
       F_CV_NODI, F_CV_NODJ) 
    !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
    IMPLICIT NONE
    ! Calculate T and DEN on the CV face at quadrature point GI.
    REAL, intent( in ) :: HDC,DT
    INTEGER, intent( in ) :: CV_DISOPT,CV_NONODS,NPHASE,NFIELD,CV_NODI,CV_NODJ,CV_ILOC,CV_JLOC,CV_SILOC,ELE,ELE2,  &
         CV_NLOC,TOTELE,SCVNGI,GI,SELE,CV_SNLOC,U_SNLOC,STOTEL, &
         WIC_T_BC_DIRICHLET,WIC_D_BC_DIRICHLET, U_NLOC,U_NONODS,NDIM
    LOGICAL, DIMENSION( NFIELD ), intent( in ) :: DOWNWIND_EXTRAP_INDIVIDUAL
    LOGICAL, intent( in ) :: DISTCONTINUOUS_METHOD, QUAD_ELEMENTS
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    INTEGER, DIMENSION( NFIELD ), intent( in ) :: SELE_LOC_WIC_F_BC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    REAL, DIMENSION(NDIM,CV_NLOC,SCVNGI), intent( in ) :: SCVFENX_ALL
    REAL, DIMENSION( CV_SNLOC ), intent( in ) :: SHAPE_CV_SNL
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFEN
!    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFENX, SCVFENY, SCVFENZ
    REAL, DIMENSION( :, :  ), intent( in ) :: CVNORMX_ALL
      REAL, DIMENSION ( NFIELD), intent( inout ) :: LIMF
      REAL, DIMENSION ( NFIELD), intent( in ) :: F_INCOME, F_NDOTQ
      REAL, DIMENSION ( NFIELD), intent( in ) :: F_CV_NODI, F_CV_NODJ

      REAL, DIMENSION ( NFIELD,CV_NLOC), intent( in ) :: LOC_F, LOC_FEMF
    REAL, DIMENSION(NFIELD, CV_SNLOC), intent( in ) :: SLOC_SUF_F_BC
    REAL, DIMENSION(NDIM,NFIELD,U_NLOC), intent( in ) ::  LOC_UF

    REAL, DIMENSION(NFIELD,CV_SNLOC), intent( in ) :: SLOC_F, SLOC_FEMF, SLOC2_F, SLOC2_FEMF

    REAL, DIMENSION( : , : , : ), intent( in ) :: INV_JAC
    REAL, DIMENSION( NFIELD  ), intent( in ) :: FUPWIND_IN, FUPWIND_OUT

!   REAL, DIMENSION( : ), intent( in ) :: TUPWIND_MAT, DENUPWIND_MAT
!   REAL, DIMENSION( : ), intent( in ) :: T2UPWIND_MAT
    ! Local variables
    ! If UPWIND then use upwind flux between elements else use central. 
    ! If HI_ORDER_HALF then use high order interpolation when around 
    ! a volume frac of 0.5 and gradually apply limiting near 0 and 1. 
    LOGICAL :: DOWNWIND_EXTRAP ! Extrapolate a downwind value for interface tracking.

    ! Scaling to reduce the downwind bias(=1downwind, =0central)
    LOGICAL, PARAMETER :: SCALE_DOWN_WIND = .true.
! The new simple limiter NEW_LIMITER can produce very slightly different results 
!    LOGICAL, PARAMETER :: NEW_LIMITER = .true.
    LOGICAL, PARAMETER :: NEW_LIMITER = .false.
    ! Non-linear Petrov-Galerkin option for interface tracking...
    ! =4 is anisotropic downwind diffusion based on a projected 1D system (1st recommend)
    ! =0 is anisotropic downwind diffusion based on a velocity projection like SUPG 
    ! (2nd recommend, most compressive)
    ! =2 is isotropic downwind diffusion  (3rd recommend,least compressive)
    ! =5 is isotropic downwind diffusion with magnitude of =0 option. 
    ! In tests they all produce similar results.
    !      INTEGER, PARAMETER :: NON_LIN_PETROV_INTERFACE = 5
    INTEGER, PARAMETER :: NON_LIN_PETROV_INTERFACE = 3
    real, parameter :: tolerance = 1.e-10

    LOGICAL :: FIRSTORD, NOLIMI, RESET_STORE, LIM_VOL_ADJUST
    REAL :: RELAX, TMIN_STORE, TMAX_STORE, TOLDMIN_STORE, &
         T2MIN_STORE, T2MAX_STORE, &
         DENMIN_STORE, DENMAX_STORE
    INTEGER :: CV_KLOC, CV_NODK, CV_NODK_IPHA, CV_KLOC2, CV_NODK2, CV_NODK2_IPHA, CV_STAR_IPHA, &
         CV_SKLOC, CV_SNODK, CV_SNODK_IPHA, U_KLOC,U_NODK,U_NODK_IPHA, IDIM, ELE_DOWN
    REAL :: T_BETWEEN_MIN, T_BETWEEN_MAX
    REAL :: T_AVE_EDGE, T_AVE_ELE
    REAL :: T_MIDVAL
    INTEGER :: U_KLOC2, U_NODK2, U_NODK2_IPHA
    INTEGER :: U_SKLOC, COUNT, COUNT_IN, COUNT_OUT, COUNT_IN_PHA, COUNT_OUT_PHA
    INTEGER :: IFIELD,IFI
 
! F:
      REAL :: FXGI_ALL(NDIM,NFIELD)

      REAL :: UDGI_ALL(NDIM,NFIELD)
      REAL :: A_STAR_X_ALL(NDIM,NFIELD)

      REAL :: courant_or_minus_one_new(NFIELD)
      REAL :: XI_LIMIT(NFIELD)

! Allocate quadrature pt variables...

      REAL :: FEMFGI(NFIELD), RGRAY(NFIELD), RSHAPE(NFIELD), DIFF_COEF(NFIELD), COEF(NFIELD) 
      REAL :: P_STAR(NFIELD), U_DOT_GRADF_GI(NFIELD), A_STAR_F(NFIELD) 
      REAL :: RESIDGI(NFIELD), ELE_LENGTH_SCALE(NFIELD), RSCALE(NFIELD), COEF2(NFIELD) 
      REAL :: FEMFGI_CENT(NFIELD), FEMFGI_UP(NFIELD) 
      REAL :: VEC_VEL2(NDIM,NFIELD) 


    ! By default do not use first-order upwinding
    FIRSTORD = .false.

    ! No limiting if CV_DISOPT is 6 or 7  (why not just define limt=femt and skip to assembly?)
    NOLIMI = ( INT( CV_DISOPT / 2 ) == 3 ) 

    ! Make a guess at the CV face value of advected field variable and density 
    ! (Depends on discetisation option, CV_DISOPT)                
    SELECT CASE( CV_DISOPT / 2 )

!    CASE( 0 ) ! First-order upwinding is achived through the limiting
!       FEMFGI(:)    = FVF(:)

    CASE( 1 ) ! Central differencing [Trapezoidal rule (2 OR 3)]
       FEMFGI(:)    = 0.5 * ( LOC_F( :, CV_ILOC ) + LOC_F( :, CV_JLOC ) )        

    CASE DEFAULT ! Finite element approximation (4 OR 5)(6 or 7)(8 or 9)
       FEMFGI(:)    = 0.0

       Conditional_CV_DISOPT_ELE2: IF ( SELE /= 0 ) THEN
          ! Is on boundary of the domain

          DO IFIELD=1,NFIELD
             IF ( SELE_LOC_WIC_F_BC(IFIELD) /= WIC_T_BC_DIRICHLET ) THEN ! Don't apply a Dirichlet bc
                LIMF(IFIELD)    = LOC_F( IFIELD, CV_ILOC )  
             ELSE
                LIMF(IFIELD)    = (1.-F_INCOME(IFIELD))*LOC_F( IFIELD, CV_ILOC )   + F_INCOME(IFIELD)*  SLOC_SUF_F_BC( IFIELD,  CV_SILOC)
             END IF
          END DO ! END OF DO IFIELD=1,NFIELD

       ELSE Conditional_CV_DISOPT_ELE2
    ! Extrapolate a downwind value for interface tracking.

          DOWNWIND_EXTRAP = ( cv_disopt>=8 )

          DO IFIELD=1,NFIELD
             IF( DOWNWIND_EXTRAP_INDIVIDUAL(IFIELD)  ) THEN 
                courant_or_minus_one_new(IFIELD) = abs ( dt * F_ndotq(IFIELD) / hdc )
                XI_LIMIT(IFIELD)=MAX(1./max(tolerance,(3.*courant_or_minus_one_new(IFIELD))),2.0) 
             else
                courant_or_minus_one_new(IFIELD) = -1.0
                XI_LIMIT(IFIELD) = 2.0
             end if
          END DO

         IF( ( ELE2 == 0 ) .OR. ( ELE2 == ELE ) ) THEN

          RSCALE(:) = 1.0 ! Scaling to reduce the downwind bias(=1downwind, =0central)
          IF ( SCALE_DOWN_WIND ) THEN

             IF ( DOWNWIND_EXTRAP  ) THEN

                DO IFIELD=1,MIN(2*NPHASE,NFIELD) 
                   DO IDIM=1,NDIM
                      FXGI_ALL(IDIM,IFIELD) = dot_product(SCVFENX_ALL(IDIM, : , GI ) , LOC_FEMF(IFIELD,:)) 
                   END DO
!                   FEMFGI(IFIELD) = dot_product( SCVFEN( : , GI ), LOC_FEMF(IFIELD,:)  )

                   DO IDIM=1,NDIM
                      UDGI_ALL(IDIM,IFIELD) = dot_product(SUFEN( : , GI ) , LOC_UF(IDIM,IFIELD, :) )
                   END DO
                END DO

                IF ( NON_LIN_PETROV_INTERFACE == 0 ) THEN ! NOT Petrov-Galerkin for interface capturing...  
                   ! no cosine rule :
                   DO IFIELD=1,MIN(2*NPHASE,NFIELD) ! It should be NPHASE because interface tracking only applied to the 1st set of fields.
                      RSCALE(IFIELD) = 1.0 / PTOLFUN( SQRT( SUM( UDGI_ALL(:,IFIELD)**2)   ) )

                      DO IDIM = 1, NDIM
                         VEC_VEL2(IDIM,IFIELD) = SUM( INV_JAC(IDIM, 1:NDIM, GI) * UDGI_ALL(1:NDIM,IFIELD) )
                      END DO
                   ! normalize the velocity in here: 
                      ELE_LENGTH_SCALE(IFIELD) = 0.5 * SQRT( SUM(UDGI_ALL(:,IFIELD)**2) ) / PTOLFUN( SUM( VEC_VEL2(:,IFIELD)**2 ) )
                   ! For discontinuous elements half the length scale...
                      IF(DISTCONTINUOUS_METHOD) ELE_LENGTH_SCALE(IFIELD)=0.5*ELE_LENGTH_SCALE(IFIELD)
                   ! For quadratic elements...
                      IF( QUAD_ELEMENTS ) ELE_LENGTH_SCALE(IFIELD)=0.5*ELE_LENGTH_SCALE(IFIELD)
                   END DO
                ELSE ! Interface capturing...

                   DO IFIELD=1,MIN(2*NPHASE,NFIELD) ! It should be NPHASE because interface tracking only applied to the 1st set of fields.

                      U_DOT_GRADF_GI(IFIELD) = SUM( UDGI_ALL(:,IFIELD)*FXGI_ALL(:,IFIELD)  ) 

                      IF ( NON_LIN_PETROV_INTERFACE == 5 ) THEN 
                         COEF(IFIELD) = 1.0 / PTOLFUN( SQRT( SUM(FXGI_ALL(IFIELD,:)**2)   ) )
                      ELSE
                         COEF(IFIELD) = U_DOT_GRADF_GI(IFIELD) / PTOLFUN( SUM( FXGI_ALL(:,IFIELD)**2 )  )
                      END IF

                      A_STAR_F(IFIELD) = 0.0
                      A_STAR_X_ALL(:,IFIELD) = COEF(IFIELD) * FXGI_ALL(:,IFIELD)

                      RESIDGI(IFIELD) = SQRT ( SUM( UDGI_ALL(:,IFIELD)**2 )  ) / HDC

                      VEC_VEL2(1:NDIM,IFIELD) = matmul( INV_JAC(:,:,GI), A_STAR_X_ALL(1:NDIM,IFIELD) )

                      P_STAR(IFIELD) = 0.5 * HDC / PTOLFUN( SQRT( SUM( A_STAR_X_ALL(:,IFIELD)**2 ))) 


                      select case (NON_LIN_PETROV_INTERFACE)
                      case ( 1 )     ! standard approach
                         DIFF_COEF(IFIELD) = COEF(IFIELD) * P_STAR(IFIELD) * RESIDGI(IFIELD)
                      case ( 2 )     ! standard approach making it +ve
                         DIFF_COEF(IFIELD) = MAX( 0.0, COEF(IFIELD) * P_STAR(IFIELD) * RESIDGI(IFIELD) )
                      case ( 3 )     ! residual squared approach
                         DIFF_COEF(IFIELD) = P_STAR(IFIELD) * RESIDGI(IFIELD)**2 / PTOLFUN( SUM(FXGI_ALL(:,IFIELD)**2)  ) 
                      case ( 4 )     ! anisotropic diffusion in the A* direction.
                         COEF2(IFIELD) =  SUM( CVNORMX_ALL(:,GI)*A_STAR_X_ALL(:,IFIELD) )  
                      case default   ! isotropic diffusion with u magnitide
                            DIFF_COEF(IFIELD) = SQRT( SUM( UDGI_ALL(:,IFIELD)**2 )) * P_STAR(IFIELD) 
                      END select
                      ! Make the diffusion coefficient negative (compressive)
                      DIFF_COEF(IFIELD) = -DIFF_COEF(IFIELD)
                      RSCALE(IFIELD) = 1. / TOLFUN( SUM(CVNORMX_ALL(:,GI)*UDGI_ALL(:,IFIELD))   )

                   END DO ! END OF DO IFIELD=1,NFIELD

 
                END IF ! Petrov-Galerkin end of IF(NON_LIN_PETROV_INTERFACE==0) THEN 

            
             END IF ! DOWNWIND_EXTRAP 
         
          END IF ! SCALE_DOWN_WIND

          FEMFGI=0.0
          DO CV_KLOC = 1, CV_NLOC
             DO IFIELD=1,NFIELD ! Only perform this loop for the 1st field which is the interface tracking field...

                IF( DOWNWIND_EXTRAP_INDIVIDUAL(IFIELD)  ) THEN ! Extrapolate to the downwind value...
             
                   IF ( NON_LIN_PETROV_INTERFACE.NE.0 ) THEN
                      IF ( NON_LIN_PETROV_INTERFACE == 4 ) THEN ! anisotropic diffusion...
                         RGRAY(IFIELD) = RSCALE(IFIELD) * COEF2(IFIELD) * P_STAR(IFIELD) * SUM( UDGI_ALL(:,IFIELD)*SCVFENX_ALL( :, CV_KLOC, GI ) )  
                      ELSE
                         RGRAY(IFIELD) = - DIFF_COEF(IFIELD) * RSCALE(IFIELD) * SUM( CVNORMX_ALL(:,GI)*SCVFENX_ALL( :, CV_KLOC, GI )  )
                      END IF
                   ELSE
                      RGRAY(IFIELD) = RSCALE(IFIELD) * ELE_LENGTH_SCALE(IFIELD) * SUM( UDGI_ALL(:,IFIELD)*SCVFENX_ALL( :, CV_KLOC, GI ) )
                   END IF
                  
                   RSHAPE(IFIELD)    = SCVFEN( CV_KLOC, GI ) + RGRAY(IFIELD)
                   FEMFGI(IFIELD)    = FEMFGI(IFIELD)     +  RSHAPE(IFIELD)     * LOC_FEMF( IFIELD, CV_KLOC) 
                ELSE

                   FEMFGI(IFIELD)    = FEMFGI(IFIELD)     +  SCVFEN( CV_KLOC, GI ) * LOC_FEMF( IFIELD, CV_KLOC)            
                END IF
             END DO ! ENDOF DO IFIELD=1,NPHASE



          END DO ! ENDOF DO CV_KLOC = 1, CV_NLOC

       ELSE  ! END OF IF( ( ELE2 == 0 ) .OR. ( ELE2 == ELE ) ) THEN  ---DG saturation across elements

             FEMFGI_CENT(:) = 0.0
             FEMFGI_UP(:)   = 0.0
             DO CV_SKLOC = 1, CV_SNLOC
! Central for DG...
                FEMFGI_CENT(:) = FEMFGI_CENT(:) +  SHAPE_CV_SNL( CV_SKLOC ) * 0.5 * ( SLOC_FEMF( :, CV_SKLOC ) & 
                              + SLOC2_FEMF( :, CV_SKLOC )    )
! Standard DG upwinding...
                FEMFGI_UP(:) = FEMFGI_UP(:) +  SHAPE_CV_SNL( CV_SKLOC ) * ( SLOC2_FEMF( :, CV_SKLOC)  & 
                     * F_INCOME(:) + SLOC_FEMF( :, CV_SKLOC) * ( 1. - F_INCOME(:) ) )
             END DO
             DO IFIELD=1,NFIELD
                IF( DOWNWIND_EXTRAP_INDIVIDUAL(IFIELD)  ) THEN ! Extrapolate to the downwind value...
                   FEMFGI(IFIELD) = FEMFGI_CENT(IFIELD)
                ELSE
                   FEMFGI(IFIELD) = FEMFGI_UP(IFIELD)
                ENDIF
             END DO
        ENDIF ! ENDOF IF( ( ELE2 == 0 ) .OR. ( ELE2 == ELE ) ) THEN ELSE


   if(NEW_LIMITER) then
    CALL ONVDLIM_ANO_MANY( NFIELD, &
       LIMF(:), FEMFGI(:), F_INCOME(:), & 
       F_CV_NODI(:), F_CV_NODJ(:),XI_LIMIT(:),  &
!       LOC_F(:,CV_ILOC), LOC_F(:,CV_JLOC),XI_LIMIT(:),  &
       FUPWIND_IN(:), FUPWIND_OUT(:) )
   else ! original limiter
     do ifield=1,nfield
       CALL ONVDLIM_ANO( cv_nonods, &
            LIMF(ifield), FEMFGI(ifield), f_INCOME(ifield), cv_nodi, cv_nodj, &
            F_CV_NODI(ifield), F_CV_NODJ(ifield),   FIRSTORD, .false., courant_or_minus_one_new(IFIELD), &
            FUPWIND_IN(ifield), FUPWIND_OUT(ifield) )
     end do
   endif

       ENDIF Conditional_CV_DISOPT_ELE2

    END SELECT

    RETURN

CONTAINS


        SUBROUTINE ONVDLIM_ANO_MANY( NFIELD, &
       TDLIM, TDCEN, INCOME, &
       ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT,  &
       TUPWIN, TUPWI2 )
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
    INTEGER, intent( in ) :: NFIELD
    REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM  
    REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2 
    REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
    ! Local variables   
    REAL :: UCIN(NFIELD), UCOU(NFIELD), TDELE(NFIELD), DENOIN(NFIELD), CTILIN(NFIELD), DENOOU(NFIELD), &
         CTILOU(NFIELD), FTILIN(NFIELD), FTILOU(NFIELD)


       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE 

       DENOIN = TOLFUN_MANY( TDELE - TUPWIN )

       UCIN = ETDNEW_PELEOT 
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT

       DENOOU = TOLFUN_MANY( TDELE - TUPWI2 )
       UCOU = ETDNEW_PELE 
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU



       FTILIN = ( TDCEN - TUPWIN ) / DENOIN
       FTILOU = ( TDCEN - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 )



    RETURN

  END SUBROUTINE ONVDLIM_ANO_MANY




        SUBROUTINE ONVDLIM_ANO( TOTELE, &
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
!       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, -1. ) * DENOIN ) &
!            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, -1. ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD
!       if(( PELE==433).and.( PELEOT==435))  print *,'++ firord,TDLIM, INCOME, TUPWI2, FTILOU, CTILOU, COURANT_OR_MINUS_ONE:', &
!                                                        firord,TDLIM, INCOME, TUPWI2, FTILOU, CTILOU, COURANT_OR_MINUS_ONE
!     if((PELE==433).and.(PELEOT==435)) then
!          print *,'---new INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT:',INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT
!          print *,'---new FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM:',FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM
!          print *,'---new NOLIMI, FIRORD, TUPWIN2, TUPWI22:', NOLIMI, FIRORD, TUPWIN2, TUPWI22
!     endif
    RETURN

  END SUBROUTINE ONVDLIM_ANO


  END SUBROUTINE GET_INT_T_DEN_NEW












  SUBROUTINE CAL_LIM_VOL_ADJUST( TMIN_STORE, TMIN, T, TMIN_NOD, RESET_STORE, MASS_CV, &
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





  REAL FUNCTION FACE_THETA( DT, CV_THETA, INTERFACE_TRACK, HDC, NDOTQ, LIMDT, DIFF_COEF_DIVDX, &
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
                IF(SELE2 == 0 .and. ele2>0) THEN ! is a volume element
                   CV_INOD=X_NDGLN((ELE-1)*X_NLOC+CV_ILOC)
                   DO CV_ILOC2=1,CV_NLOC
                      CV_INOD2=X_NDGLN((ELE2-1)*X_NLOC+CV_ILOC2) 
                      IF(CV_INOD == CV_INOD2) FOUND=.TRUE.
                   END DO
                ELSE if (sele2>0) then ! is a surface element
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






    SUBROUTINE CALC_LIMIT_MATRIX_MAX_MIN(TMAX_ALL, TMIN_ALL, DENMAX_ALL, DENMIN_ALL, &
       T2MAX_ALL, T2MIN_ALL, &
       T_ALL,  T2_ALL, DEN_ALL, IGOT_T2, NPHASE, CV_NONODS, &
       TMIN_NOD_ALL, TMAX_NOD_ALL,  &
       T2MIN_NOD_ALL, T2MAX_NOD_ALL, &
       DENMIN_NOD_ALL, DENMAX_NOD_ALL, &
       NSMALL_COLM, SMALL_FINDRM, SMALL_COLM, &
       TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL, MASS_CV)
! Populate  limiting matrix based on max and min values
    ! For each node, find the largest and smallest value of T and 
    ! DENSITY for both the current and previous timestep, out of 
    ! the node value and all its surrounding nodes including Dirichlet b.c's.
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, CV_NONODS, NSMALL_COLM, IGOT_T2
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( :, : ), intent( inout ) :: TMAX_ALL, TMIN_ALL, DENMAX_ALL, DENMIN_ALL
    REAL, DIMENSION( :, : ), intent( inout ) :: T2MAX_ALL, T2MIN_ALL
    REAL, DIMENSION( :, : ), intent( inout ) :: TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL

    REAL, DIMENSION( :, : ), intent( in ) :: T_ALL, DEN_ALL
    REAL, DIMENSION( :, :), intent( in ) :: T2_ALL
    REAL, DIMENSION( : ), intent( in ) :: MASS_CV
    INTEGER, DIMENSION( :, : ), intent( inout ) :: TMIN_NOD_ALL, TMAX_NOD_ALL, DENMIN_NOD_ALL, DENMAX_NOD_ALL
    INTEGER, DIMENSION( :, : ), intent( inout ) :: T2MIN_NOD_ALL, T2MAX_NOD_ALL
    ! Local variables
    INTEGER :: CV_NODI, CV_NODJ, IPHASE, COUNT
    INTEGER :: COUNT2, COUNT_IN, COUNT_OUT
    LOGICAL, PARAMETER :: LIM_VOL_ADJUST=.TRUE.
    REAL :: TMIN_STORE(NPHASE),TMAX_STORE(NPHASE),DENMIN_STORE(NPHASE),DENMAX_STORE(NPHASE)
    REAL :: T2MIN_STORE(NPHASE),T2MAX_STORE(NPHASE)
    LOGICAL :: RESET_STORE, NO_RESET_STORE
    REAL :: INCOME


    TUPWIND_MAT_ALL(1:NPHASE, 1: NSMALL_COLM)  =0.0
    DENUPWIND_MAT_ALL(1:NPHASE, 1: NSMALL_COLM)=0.0 
    IF(IGOT_T2==1) T2UPWIND_MAT_ALL(1:NPHASE, 1: NSMALL_COLM) =0.0


    DO CV_NODI=1,CV_NONODS
       DO COUNT=SMALL_FINDRM(CV_NODI), SMALL_FINDRM(CV_NODI+1)-1

          CV_NODJ=SMALL_COLM(COUNT)

! for outgoing information to CV_NODI ...

    INCOME=0.0

    DO IPHASE=1,NPHASE
!       CV_NODI_IPHA = CV_NODI + (IPHASE-1)*CV_NONODS
!       CV_NODJ_IPHA = CV_NODJ + (IPHASE-1)*CV_NONODS
       IF ( LIM_VOL_ADJUST ) THEN
          RESET_STORE = .FALSE.
          CALL CAL_LIM_VOL_ADJUST( TMIN_STORE(IPHASE), TMIN_ALL(IPHASE,:), T_ALL(IPHASE,:), TMIN_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
          CALL CAL_LIM_VOL_ADJUST( TMAX_STORE(IPHASE), TMAX_ALL(IPHASE,:), T_ALL(IPHASE,:), TMAX_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
          CALL CAL_LIM_VOL_ADJUST( DENMIN_STORE(IPHASE), DENMIN_ALL(IPHASE,:), DEN_ALL(IPHASE,:), DENMIN_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
          CALL CAL_LIM_VOL_ADJUST( DENMAX_STORE(IPHASE), DENMAX_ALL(IPHASE,:), DEN_ALL(IPHASE,:), DENMAX_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )

          IF(IGOT_T2==1) THEN
             CALL CAL_LIM_VOL_ADJUST( T2MIN_STORE(IPHASE), T2MIN_ALL(IPHASE,:), T2_ALL(IPHASE,:), T2MIN_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
               CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
             CALL CAL_LIM_VOL_ADJUST( T2MAX_STORE(IPHASE), T2MAX_ALL(IPHASE,:), T2_ALL(IPHASE,:), T2MAX_NOD_ALL(IPHASE,:), RESET_STORE, MASS_CV, &
               CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
          END IF

       END IF
    END DO



! ***PUT INTO MATRIX**************
! Populate  limiting matrix based on max and min values


       COUNT_OUT= COUNT

       DO IPHASE=1,NPHASE

           IF(T_ALL( IPHASE, CV_NODI ).GT.T_ALL( IPHASE, CV_NODJ )) THEN
              TUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = TMAX_ALL( IPHASE, CV_NODI ) 
              DENUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = DENMAX_ALL( IPHASE, CV_NODI ) 
              IF(IGOT_T2==1) THEN
                 T2UPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = T2MAX_ALL( IPHASE, CV_NODI )
              ENDIF
           ELSE
              TUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = TMIN_ALL( IPHASE, CV_NODI )
              DENUPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = DENMIN_ALL( IPHASE, CV_NODI )
              IF(IGOT_T2==1) THEN
                 T2UPWIND_MAT_ALL( IPHASE, COUNT_OUT ) = T2MIN_ALL( IPHASE, CV_NODI )
              ENDIF
           ENDIF
       END DO




    DO IPHASE=1,NPHASE
    IF ( LIM_VOL_ADJUST ) THEN
       RESET_STORE = .TRUE. 
       CALL CAL_LIM_VOL_ADJUST(TMIN_STORE(IPHASE),TMIN_ALL(IPHASE,:),T_ALL(IPHASE,:),TMIN_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
       CALL CAL_LIM_VOL_ADJUST(TMAX_STORE(IPHASE),TMAX_ALL(IPHASE,:),T_ALL(IPHASE,:),TMAX_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )

       CALL CAL_LIM_VOL_ADJUST(DENMIN_STORE(IPHASE),DENMIN_ALL(IPHASE,:),DEN_ALL(IPHASE,:),DENMIN_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
       CALL CAL_LIM_VOL_ADJUST(DENMAX_STORE(IPHASE),DENMAX_ALL(IPHASE,:),DEN_ALL(IPHASE,:),DENMAX_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
            CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )

       IF ( IGOT_T2 == 1 ) THEN
          CALL CAL_LIM_VOL_ADJUST(T2MIN_STORE(IPHASE),T2MIN_ALL(IPHASE,:),T2_ALL(IPHASE,:),T2MIN_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
               CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
          CALL CAL_LIM_VOL_ADJUST(T2MAX_STORE(IPHASE),T2MAX_ALL(IPHASE,:),T2_ALL(IPHASE,:),T2MAX_NOD_ALL(IPHASE,:),RESET_STORE,MASS_CV, &
               CV_NODI, CV_NODJ, 1, CV_NONODS, 1, INCOME )
       END IF

    END IF
    END DO



       END DO
    END DO



    RETURN

  END SUBROUTINE CALC_LIMIT_MATRIX_MAX_MIN






           SUBROUTINE ISOTROPIC_LIMITER_ALL( &
! FOR SUB SURRO_CV_MINMAX:
           T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, nsmall_colm, SMALL_CENTRM, SMALL_FINDRM, SMALL_COLM, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
           MASS_CV, &
! FOR SUB CALC_LIMIT_MATRIX_MAX_MIN:
           TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
           TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL )
      INTEGER, intent( in ) :: IGOT_T2, NPHASE, CV_NONODS, nsmall_colm, STOTEL, CV_SNLOC
      REAL, DIMENSION( :, : ), intent( in ) :: T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL
      REAL, DIMENSION( : ), intent( in ) :: MASS_CV
      REAL, DIMENSION( :, : ), intent( inout ) :: TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
                                                  TUPWIND_MAT_ALL,    DENUPWIND_MAT_ALL,    T2UPWIND_MAT_ALL
      INTEGER, DIMENSION( : ), intent( in ) :: SMALL_CENTRM, SMALL_FINDRM
      INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
      INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
      REAL, DIMENSION( :, :, : ), intent( in ), pointer :: SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL
      INTEGER, DIMENSION( :,:, : ), intent( in ) :: WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL

! Local variables...
      REAL, DIMENSION( :, : ), allocatable :: TMIN_ALL, TMAX_ALL, TOLDMIN_ALL, TOLDMAX_ALL, &
              T2MIN_ALL, T2MAX_ALL, T2OLDMIN_ALL, T2OLDMAX_ALL, DENMIN_ALL, DENMAX_ALL,  &
              DENOLDMIN_ALL, DENOLDMAX_ALL  
      INTEGER, DIMENSION( :, : ), allocatable ::TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
                   T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
                   DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL
      INTEGER :: CV_NODI, IMID, IPHASE

      ! Allocate memory for terms needed by GETGXYZ OR ONVDLIM
      ALLOCATE(      TMIN_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(      TMAX_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(   TOLDMIN_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(   TOLDMAX_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(      T2MIN_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(      T2MAX_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(   T2OLDMIN_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(   T2OLDMAX_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(    DENMIN_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(    DENMAX_ALL( NPHASE, CV_NONODS) )
      ALLOCATE( DENOLDMIN_ALL( NPHASE, CV_NONODS) )
      ALLOCATE( DENOLDMAX_ALL( NPHASE, CV_NONODS) )

      ALLOCATE(      TMIN_NOD_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(      TMAX_NOD_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(   TOLDMIN_NOD_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(   TOLDMAX_NOD_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(      T2MIN_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(      T2MAX_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(   T2OLDMIN_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(   T2OLDMAX_NOD_ALL( NPHASE, CV_NONODS* IGOT_T2) )
      ALLOCATE(    DENMIN_NOD_ALL( NPHASE, CV_NONODS) )
      ALLOCATE(    DENMAX_NOD_ALL( NPHASE, CV_NONODS) )
      ALLOCATE( DENOLDMIN_NOD_ALL( NPHASE, CV_NONODS) )
      ALLOCATE( DENOLDMAX_NOD_ALL( NPHASE, CV_NONODS) )

      ! For each node, find the largest and smallest value of T and
      ! DENSITY for both the current and previous timestep, out of
      ! the node value and all its surrounding nodes including Dirichlet b.c's.
      CALL SURRO_CV_MINMAX( TMAX_ALL, TMIN_ALL, TOLDMAX_ALL, TOLDMIN_ALL, DENMAX_ALL, DENMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL, &
           T2MAX_ALL, T2MIN_ALL, T2OLDMAX_ALL, T2OLDMIN_ALL, &
           T_ALL, TOLD_ALL, T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, size(small_colm), SMALL_FINDRM, SMALL_COLM, &
           STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC_ALL, SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
           MASS_CV, TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
           T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
           DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL )

! Populate  limiting matrix based on max and min values OLD:
         CALL CALC_LIMIT_MATRIX_MAX_MIN(TOLDMAX_ALL, TOLDMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL, &
       T2OLDMAX_ALL, T2OLDMIN_ALL,  &
       TOLD_ALL,  T2OLD_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, &
       TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL,  &
       T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL,  &
       DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL,  &
       NSMALL_COLM, SMALL_FINDRM, SMALL_COLM,   &
       TOLDUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, MASS_CV )

! Populate  limiting matrix based on max and min values
         CALL CALC_LIMIT_MATRIX_MAX_MIN(TMAX_ALL, TMIN_ALL, DENMAX_ALL, DENMIN_ALL, &
       T2MAX_ALL, T2MIN_ALL,  &
       T_ALL,  T2_ALL, DEN_ALL, IGOT_T2, NPHASE, CV_NONODS, &
       TMIN_NOD_ALL, TMAX_NOD_ALL,  &
       T2MIN_NOD_ALL, T2MAX_NOD_ALL,  &
       DENMIN_NOD_ALL, DENMAX_NOD_ALL,  &
       NSMALL_COLM, SMALL_FINDRM, SMALL_COLM,   &
       TUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, T2UPWIND_MAT_ALL, MASS_CV )

! make sure the diagonal is equal to the value:
         DO CV_NODI=1,CV_NONODS
            IMID=SMALL_CENTRM(CV_NODI)
            DO IPHASE=1,NPHASE
               TUPWIND_MAT_ALL( IPHASE, IMID)=T_ALL( IPHASE, CV_NODI)
               TOLDUPWIND_MAT_ALL(  IPHASE, IMID)=TOLD_ALL( IPHASE, CV_NODI)

               DENUPWIND_MAT_ALL(  IPHASE, IMID)=DEN_ALL( IPHASE, CV_NODI)
               DENOLDUPWIND_MAT_ALL(  IPHASE, IMID)=DENOLD_ALL( IPHASE, CV_NODI)

               IF( IGOT_T2 == 1 ) THEN
                  T2UPWIND_MAT_ALL(IPHASE, IMID)=T2_ALL( IPHASE, CV_NODI)
                  T2OLDUPWIND_MAT_ALL(IPHASE, IMID)=T2OLD_ALL( IPHASE, CV_NODI)
               ENDIF
            END DO
         END DO

      DEALLOCATE( TMIN_ALL, TMAX_ALL, TOLDMIN_ALL, TOLDMAX_ALL, &
              T2MIN_ALL, T2MAX_ALL, T2OLDMIN_ALL, T2OLDMAX_ALL, DENMIN_ALL, DENMAX_ALL,  &
              DENOLDMIN_ALL, DENOLDMAX_ALL )
      DEALLOCATE( TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
                   T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
                   DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL )


      RETURN
      END SUBROUTINE ISOTROPIC_LIMITER_ALL







  SUBROUTINE SURRO_CV_MINMAX( TMAX_ALL, TMIN_ALL, TOLDMAX_ALL, TOLDMIN_ALL, DENMAX_ALL, DENMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL, &
       T2MAX_ALL, T2MIN_ALL, T2OLDMAX_ALL, T2OLDMIN_ALL, &
       T_ALL, TOLD_ALL,  T2_ALL, T2OLD_ALL, DEN_ALL, DENOLD_ALL, IGOT_T2, NPHASE, CV_NONODS, NCOLACV, FINACV, COLACV, &
       STOTEL, CV_SNLOC, CV_SNDGLN, SUF_T_BC_ALL,  SUF_T2_BC_ALL, SUF_D_BC_ALL, WIC_T_BC_ALL, WIC_T2_BC_ALL, WIC_D_BC_ALL, &
      MASS_CV, TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, TOLDMAX_NOD_ALL, &
       T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, T2OLDMAX_NOD_ALL, &
       DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL )
    ! For each node, find the largest and smallest value of T and 
    ! DENSITY for both the current and previous timestep, out of 
    ! the node value and all its surrounding nodes including Dirichlet b.c's.
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE,CV_NONODS, NCOLACV,STOTEL,CV_SNLOC, &
         IGOT_T2
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SNDGLN
    REAL, DIMENSION( :, :, : ), intent( in ) :: SUF_T_BC_ALL, SUF_D_BC_ALL
    REAL, DIMENSION( :, :, : ), intent( in ), pointer :: SUF_T2_BC_ALL
    INTEGER, DIMENSION( : , : , : ), intent( in ) :: WIC_T_BC_ALL, WIC_D_BC_ALL
    INTEGER, DIMENSION( : , : , : ), intent( in ) :: WIC_T2_BC_ALL
    INTEGER, DIMENSION( : ), intent( in ) :: FINACV
    INTEGER, DIMENSION( : ), intent( in ), target :: COLACV
    REAL, DIMENSION( :, : ), intent( inout ) :: TMAX_ALL, TMIN_ALL, TOLDMAX_ALL, TOLDMIN_ALL,  &
         DENMAX_ALL, DENMIN_ALL, DENOLDMAX_ALL, DENOLDMIN_ALL
    REAL, DIMENSION( :, : ), intent( inout ) :: T2MAX_ALL, T2MIN_ALL, T2OLDMAX_ALL, T2OLDMIN_ALL

    REAL, DIMENSION( :, : ), intent( in ) :: T_ALL,TOLD_ALL,DEN_ALL,DENOLD_ALL
    REAL, DIMENSION( :, : ), intent( in ) :: T2_ALL,T2OLD_ALL
    REAL, DIMENSION( : ), intent( in ) :: MASS_CV
    INTEGER, DIMENSION( :, : ), intent( inout ) :: TMIN_NOD_ALL, TMAX_NOD_ALL, TOLDMIN_NOD_ALL, &
         TOLDMAX_NOD_ALL, DENMIN_NOD_ALL, DENMAX_NOD_ALL, DENOLDMIN_NOD_ALL, DENOLDMAX_NOD_ALL
    INTEGER, DIMENSION( :, : ), intent( inout ) :: T2MIN_NOD_ALL, T2MAX_NOD_ALL, T2OLDMIN_NOD_ALL, &
         T2OLDMAX_NOD_ALL
    ! Local variables
    INTEGER :: CV_NODI, CV_NODJ, IPHASE, COUNT, CV_SILOC, SELE, CV_INOD 
    integer, dimension(:), pointer :: cv_neigh_ptr

    Loop_CV_NODI: DO CV_NODI = 1, CV_NONODS

       cv_neigh_ptr=>colacv(finacv(cv_nodi):finacv(cv_nodi+1)-1)

       DO IPHASE = 1, NPHASE
          TMAX_ALL( IPHASE, CV_NODI ) = maxval(T_ALL( IPHASE, cv_neigh_ptr ))
          TMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(T_ALL( IPHASE, cv_neigh_ptr )))  ! COLN OF THE MAXIMUM VALUE
          TMIN_ALL( IPHASE, CV_NODI ) = minval(T_ALL( IPHASE, cv_neigh_ptr ))
          TMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(T_ALL( IPHASE, cv_neigh_ptr )))  
          TOLDMAX_ALL( IPHASE, CV_NODI ) = maxval(TOLD_ALL( IPHASE, cv_neigh_ptr ))
          TOLDMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(TOLD_ALL( IPHASE, cv_neigh_ptr )))  
          TOLDMIN_ALL( IPHASE, CV_NODI ) = minval(TOLD_ALL( IPHASE, cv_neigh_ptr ))
          TOLDMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(TOLD_ALL( IPHASE, cv_neigh_ptr )))  
          IF(IGOT_T2==1) THEN
             T2MAX_ALL( IPHASE, CV_NODI ) = maxval(T2_ALL( IPHASE, cv_neigh_ptr ))
             T2MAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) =cv_neigh_ptr(maxloc(T2_ALL( IPHASE, cv_neigh_ptr )))  
             T2MIN_ALL( IPHASE, CV_NODI ) = minval(T2_ALL( IPHASE, cv_neigh_ptr ))
             T2MIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(T2_ALL( IPHASE, cv_neigh_ptr )))  
             T2OLDMAX_ALL( IPHASE, CV_NODI ) = maxval(T2OLD_ALL( IPHASE, cv_neigh_ptr ))
             T2OLDMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(T2OLD_ALL( IPHASE, cv_neigh_ptr )))  
             T2OLDMIN_ALL( IPHASE, CV_NODI ) = minval(T2OLD_ALL( IPHASE, cv_neigh_ptr ))
             T2OLDMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(T2OLD_ALL( IPHASE, cv_neigh_ptr )))
          ENDIF
          DENMAX_ALL( IPHASE, CV_NODI ) = maxval(DEN_ALL( IPHASE, cv_neigh_ptr ))
          DENMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(DEN_ALL( IPHASE, cv_neigh_ptr )))  
          DENMIN_ALL( IPHASE, CV_NODI ) = minval(DEN_ALL( IPHASE, cv_neigh_ptr ))
          DENMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(DEN_ALL( IPHASE, cv_neigh_ptr )))  
          DENOLDMAX_ALL( IPHASE, CV_NODI ) = maxval(DENOLD_ALL( IPHASE, cv_neigh_ptr ))
          DENOLDMAX_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(maxloc(DENOLD_ALL( IPHASE, cv_neigh_ptr )))  
          DENOLDMIN_ALL( IPHASE, CV_NODI ) = minval(DENOLD_ALL( IPHASE, cv_neigh_ptr ))
          DENOLDMIN_NOD_ALL( IPHASE, CV_NODI:CV_NODI ) = cv_neigh_ptr(minloc(DENOLD_ALL( IPHASE, cv_neigh_ptr )))  
       END DO

    END DO Loop_CV_NODI

    ! Take into account the Dirichlet b.c's when working out max and min values.
    Loop_SELE: DO SELE= 1, STOTEL

       Loop_CV_SILOC: DO CV_SILOC = 1, CV_SNLOC

          CV_INOD=CV_SNDGLN((SELE-1)*CV_SNLOC+CV_SILOC)
!          SUF_CV_SI=(SELE-1)*CV_SNLOC+CV_SILOC

          DO IPHASE=1,NPHASE
!             SUF_CV_SI_IPHA = SUF_CV_SI + STOTEL * CV_SNLOC * ( IPHASE - 1 )
!             CV_INOD_IPHA=CV_INOD + CV_NONODS*(IPHASE-1)
             IF( (WIC_T_BC_ALL(1 ,IPHASE, SELE) == WIC_T_BC_DIRICHLET) &
                  .OR.(WIC_T_BC_ALL(1, IPHASE, SELE) == WIC_T_BC_DIRI_ADV_AND_ROBIN)) THEN
                IF(SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*(SELE-1) ) > TMAX_ALL( IPHASE, CV_INOD ) ) THEN
                   TMAX_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL( 1,  IPHASE, CV_SILOC+CV_SNLOC*(SELE-1) )
                   TMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF
                IF(SUF_T_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC*(SELE-1) ) < TMIN_ALL( IPHASE, CV_INOD ) ) THEN
                   TMIN_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL( 1, IPHASE, CV_SILOC+ CV_SNLOC*(SELE-1) )
                   TMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF

                IF(SUF_T_BC_ALL( 1, IPHASE, CV_SILOC+ CV_SNLOC*(SELE-1) ) > TOLDMAX_ALL( IPHASE, CV_INOD ) ) THEN
                   TOLDMAX_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL(1, IPHASE, CV_SILOC + CV_SNLOC*( SELE-1 ) )
                   TOLDMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF
                IF(SUF_T_BC_ALL( 1 , IPHASE, CV_SILOC + CV_SNLOC* ( SELE -1 ) ) < TOLDMIN_ALL( IPHASE, CV_INOD ) ) THEN
                   TOLDMIN_ALL( IPHASE, CV_INOD ) = SUF_T_BC_ALL( 1 , IPHASE, CV_SILOC + CV_SNLOC* ( SELE -1) )
                   TOLDMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF
             ENDIF

             ! T2:
             IF(IGOT_T2==1) THEN
                IF( (WIC_T2_BC_ALL(1, IPHASE, SELE) == WIC_T_BC_DIRICHLET) &
                     .OR.(WIC_T2_BC_ALL(1, IPHASE, SELE) == WIC_T_BC_DIRI_ADV_AND_ROBIN)) THEN
                   IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC* ( SELE-1 ) ) > T2MAX_ALL( IPHASE, CV_INOD ) ) THEN
                      T2MAX_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                      T2MAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                   ENDIF
                   IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1) ) < T2MIN_ALL( IPHASE, CV_INOD ) ) THEN
                      T2MIN_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOc * ( SELE - 1 ) )
                      T2MIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                   ENDIF

                   IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC* ( SELE - 1 ) ) > T2OLDMAX_ALL( IPHASE, CV_INOD ) ) THEN
                      T2OLDMAX_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                      T2OLDMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                   ENDIF
                   IF(SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) ) < T2OLDMIN_ALL( IPHASE, CV_INOD ) ) THEN
                      T2OLDMIN_ALL( IPHASE, CV_INOD ) = SUF_T2_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                      T2OLDMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                   ENDIF
                ENDIF
             ENDIF
             ! DEN: 
             IF( WIC_D_BC_ALL(1 , IPHASE, SELE) == WIC_D_BC_DIRICHLET ) THEN
                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) ) > DENMAX_ALL( IPHASE, CV_INOD ) ) THEN
                   DENMAX_ALL( IPHASE, CV_INOD ) = SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                   DENMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF
                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOc * ( SELE - 1 ) ) < DENMIN_ALL( IPHASE, CV_INOD ) ) THEN
                   DENMIN_ALL( IPHASE, CV_INOD ) = SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                   DENMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF

                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 )  ) > DENOLDMAX_ALL( IPHASE, CV_INOD ) ) THEN
                   DENOLDMAX_ALL( IPHASE, CV_INOD ) = SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                   DENOLDMAX_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF
                IF(SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) ) < DENOLDMIN_ALL( IPHASE, CV_INOD ) ) THEN
                   DENOLDMIN_ALL( IPHASE, CV_INOD )= SUF_D_BC_ALL( 1, IPHASE, CV_SILOC + CV_SNLOC * ( SELE - 1 ) )
                   DENOLDMIN_NOD_ALL( IPHASE, CV_INOD ) =  CV_INOD
                ENDIF
             ENDIF
          END DO

       END DO Loop_CV_SILOC

    END DO Loop_SELE

    RETURN
  END SUBROUTINE SURRO_CV_MINMAX
 




  SUBROUTINE CALC_SELE( ELE, ELE3, SELE, CV_SILOC, CV_ILOC, U_SLOC2LOC, CV_SLOC2LOC, &
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
    INTEGER, intent( inout ) :: SELE, ELE3, CV_SILOC 
    INTEGER, DIMENSION( : ), intent( inout ) :: U_SLOC2LOC
    INTEGER, DIMENSION( : ), intent( inout ) :: CV_SLOC2LOC
    ! local variables
    INTEGER :: IFACE, ELE2, SELE2, CV_JLOC, CV_JNOD, &
         U_JLOC, U_JNOD, CV_KLOC, CV_SKNOD, &
         U_KLOC, U_SKLOC, U_SKNOD, CV_SKLOC, CV_SKLOC2, I
    LOGICAL :: FOUND
    INTEGER, DIMENSION( CV_SNLOC ) :: LOG_ON_BOUND

    !ewrite(3,*)'In Calc_Sele'
    I = 1
    DO CV_JLOC = 1, CV_NLOC  
       CV_JNOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_JLOC )
       IF ( .NOT.CVFEM_ON_FACE( CV_JLOC ) ) THEN
          LOG_ON_BOUND( I ) = CV_JNOD
          I = I + 1
       END IF
    END DO

    SELE = 0
    ELE3 = 0
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
          IF( FOUND ) THEN
             SELE = SELE2
             ELE3 = ELE2
          ENDIF
       END IF
    END DO

!    IF(SELE==0) THEN
!       IF(ELE3==0) THEN
! dont integrate here as we are not between elements and or on the boundary of the domain - we are on the 
! boundary of a subdomain partition between subdomains. 
!          
!       ENDIF
!    ENDIF

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





  SUBROUTINE PUT_IN_CT_RHS( CT, CT_RHS, U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
       CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE, ELE, ELE2, SELE, &
       JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, U_OTHER_LOC, U_NDGLN,  NU_ALL,  &
       SUFEN, SCVDETWEI, CVNORMX_ALL, DEN_ALL, CV_NODI, CV_NODJ, &
       UGI_COEF_ELE_ALL,  &
       UGI_COEF_ELE2_ALL,  &
       NDOTQ, NDOTQOLD, NDOTQ_HAT, LIMD, LIMT, LIMTOLD, LIMDT, LIMDTOLD, LIMT_HAT, &
       FTHETA_T2, ONE_M_FTHETA_T2OLD, FTHETA_T2_J, ONE_M_FTHETA_T2OLD_J, integrate_other_side_and_not_boundary, &
       RETRIEVE_SOLID_CTY,theta_cty_solid)
    ! This subroutine caculates the discretised cty eqn acting on the velocities i.e. CT, CT_RHS
    IMPLICIT NONE
    INTEGER, intent( in ) :: U_NLOC, SCVNGI, GI, NCOLCT, NDIM, &
         CV_NONODS, U_NONODS, NPHASE, IPHASE, TOTELE,  ELE, ELE2, SELE, &
         CV_NODI, CV_NODJ
    LOGICAL, intent( in ) :: integrate_other_side_and_not_boundary, RETRIEVE_SOLID_CTY
    INTEGER, DIMENSION( : ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: JCOUNT_KLOC, JCOUNT_KLOC2, ICOUNT_KLOC, ICOUNT_KLOC2, U_OTHER_LOC
    REAL, DIMENSION( :, :, : ), intent( inout ) :: CT
    REAL, DIMENSION( : ), intent( inout ) :: CT_RHS
    REAL, DIMENSION( NDIM, NPHASE, U_NLOC ), intent( in ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
    REAL, DIMENSION( :, : ), intent( in ) :: SUFEN
    REAL, DIMENSION( : ), intent( in ) :: SCVDETWEI
    REAL, DIMENSION( NDIM, SCVNGI ), intent( in ) :: CVNORMX_ALL
    REAL, DIMENSION( NDIM, NPHASE, U_NONODS ), intent( in ) :: NU_ALL
    REAL, DIMENSION( NPHASE, CV_NONODS ), intent( in ) :: DEN_ALL
    REAL, intent( in ) :: NDOTQ, NDOTQOLD, NDOTQ_HAT, LIMT, LIMTOLD, LIMDT, LIMDTOLD, LIMT_HAT, LIMD
    ! LIMT_HAT is the normalised voln fraction
    REAL, intent( in ) :: FTHETA_T2, ONE_M_FTHETA_T2OLD, FTHETA_T2_J, ONE_M_FTHETA_T2OLD_J, theta_cty_solid 

    ! Local variables...
    INTEGER :: U_KLOC, U_KLOC2, JCOUNT_IPHA, IDIM, U_NODK, U_NODK_IPHA, JCOUNT2_IPHA, &
         U_KLOC_LEV, U_NLOC_LEV
    REAL :: RCON,RCON_J, UDGI_IMP_ALL(NDIM), NDOTQ_IMP
    REAL :: UDGI_ALL(NDIM), UOLDDGI_ALL(NDIM), UDGI_HAT_ALL(NDIM)

    IF ( RETRIEVE_SOLID_CTY ) THEN ! For solid modelling...
       ! Use backward Euler... (This is for the div uhat term - we subtract what we put in the CT matrix and add what we really want)
       CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) + THETA_CTY_SOLID * SCVDETWEI( GI ) * ( LIMT_HAT*NDOTQ - NDOTQ_HAT/REAL(NPHASE) ) 
! assume cty is satified for solids...
       CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) + (1.0-THETA_CTY_SOLID) * SCVDETWEI( GI ) * (LIMT_HAT - LIMT)*NDOTQ 
       ! flux from the other side (change of sign because normal is -ve)...
       if ( integrate_other_side_and_not_boundary ) then
! assume cty is satified for solids...
          CT_RHS( CV_NODJ ) = CT_RHS( CV_NODJ ) - THETA_CTY_SOLID * SCVDETWEI( GI ) * ( LIMT_HAT*NDOTQ - NDOTQ_HAT / REAL(NPHASE) )
! assume cty is satified for solids...
          CT_RHS( CV_NODJ ) = CT_RHS( CV_NODJ ) - (1.0-THETA_CTY_SOLID) * SCVDETWEI( GI ) * (LIMT_HAT - LIMT)*NDOTQ
       end if
    END IF ! For solid modelling...

    DO U_KLOC = 1, U_NLOC

       RCON = SCVDETWEI( GI ) * FTHETA_T2 * LIMDT &               
            * SUFEN( U_KLOC, GI ) / DEN_ALL( IPHASE, CV_NODI )

       IF ( RETRIEVE_SOLID_CTY ) THEN ! For solid modelling use backward Euler for this part...
          RCON = RCON + SCVDETWEI( GI ) * (LIMT_HAT - LIMT) &
               * SUFEN( U_KLOC, GI )
       END IF ! For solid modelling...

       DO IDIM = 1, NDIM
          CT( IDIM, IPHASE, JCOUNT_KLOC( U_KLOC ) ) &
               = CT( IDIM, IPHASE, JCOUNT_KLOC( U_KLOC ) ) &
               + RCON * UGI_COEF_ELE_ALL( IDIM, IPHASE, U_KLOC ) * CVNORMX_ALL( IDIM, GI )
       END DO
       ! flux from the other side (change of sign because normal is -ve)...
       if ( integrate_other_side_and_not_boundary ) then
          RCON_J = SCVDETWEI( GI ) * FTHETA_T2_J * LIMDT &
               * SUFEN( U_KLOC, GI ) / DEN_ALL( IPHASE, CV_NODJ )
          IF ( RETRIEVE_SOLID_CTY ) THEN ! For solid modelling...
             RCON_J = RCON_J  + SCVDETWEI( GI ) * (LIMT_HAT - LIMT) &
                  * SUFEN( U_KLOC, GI ) 
          END IF ! For solid modelling...

          DO IDIM = 1, NDIM
             CT( IDIM, IPHASE, ICOUNT_KLOC( U_KLOC ) ) &
                  = CT( IDIM, IPHASE, ICOUNT_KLOC( U_KLOC ) ) &
                  - RCON_J * UGI_COEF_ELE_ALL( IDIM, IPHASE, U_KLOC ) * CVNORMX_ALL( IDIM, GI )
          END DO
       end if

    END DO

    IF ( SELE /= 0 ) THEN
       UDGI_IMP_ALL=0.0
       DO U_KLOC = 1, U_NLOC
          U_NODK = U_NDGLN( ( ELE - 1 ) * U_NLOC + U_KLOC )
          UDGI_IMP_ALL(:) = UDGI_IMP_ALL(:) + SUFEN( U_KLOC, GI ) * &
               UGI_COEF_ELE_ALL( :, IPHASE, U_KLOC ) * NU_ALL( :, IPHASE, U_NODK ) 
       END DO

       NDOTQ_IMP= SUM( CVNORMX_ALL( :,GI ) * UDGI_IMP_ALL(:) ) 

       CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - SCVDETWEI( GI ) * ( &
            ONE_M_FTHETA_T2OLD * LIMDTOLD * NDOTQOLD &
            + FTHETA_T2  * LIMDT * (NDOTQ-NDOTQ_IMP) &
            ) / DEN_ALL( IPHASE, CV_NODI )

    ELSE

       CT_RHS( CV_NODI ) = CT_RHS( CV_NODI ) - SCVDETWEI( GI ) * ( &
            ONE_M_FTHETA_T2OLD * LIMDTOLD * NDOTQOLD &
            ) / DEN_ALL( IPHASE, CV_NODI ) 

       ! flux from the other side (change of sign because normal is -ve)...
       if ( integrate_other_side_and_not_boundary ) then
          CT_RHS( CV_NODJ ) = CT_RHS( CV_NODJ ) + SCVDETWEI( GI ) * ( &
               ONE_M_FTHETA_T2OLD_J * LIMDTOLD * NDOTQOLD &
               ) / DEN_ALL( IPHASE, CV_NODJ ) 
       end if
    END IF

    IF ( (ELE2 /= 0) .AND. (ELE2 /= ELE) ) THEN
       ! We have a discontinuity between elements so integrate along the face...
       DO U_KLOC = 1, U_NLOC
          U_KLOC2 = U_OTHER_LOC( U_KLOC)
          IF ( U_KLOC2 /= 0 ) THEN

             RCON = SCVDETWEI( GI ) * FTHETA_T2 * LIMDT  &
                  * SUFEN( U_KLOC, GI ) / DEN_ALL( IPHASE, CV_NODI )
             IF ( RETRIEVE_SOLID_CTY ) THEN ! For solid modelling use backward Euler for this part...
                RCON    = RCON    + SCVDETWEI( GI )  * (LIMT_HAT - LIMT)  &
                     * SUFEN( U_KLOC, GI ) 
             END IF ! For solid modelling...

             DO IDIM = 1, NDIM
                CT( IDIM, IPHASE, JCOUNT_KLOC2( U_KLOC2 ) ) &
                     = CT( IDIM, IPHASE, JCOUNT_KLOC2( U_KLOC2 ) ) &
                     + RCON * UGI_COEF_ELE2_ALL( IDIM, IPHASE, U_KLOC2 ) * CVNORMX_ALL( IDIM, GI )
             END DO
             ! flux from the other side (change of sign because normal is -ve)...
             if ( integrate_other_side_and_not_boundary ) then
                RCON_J = SCVDETWEI( GI ) * FTHETA_T2_J * LIMDT  &
                     * SUFEN( U_KLOC, GI ) / DEN_ALL( IPHASE, CV_NODJ )
                IF(RETRIEVE_SOLID_CTY) THEN ! For solid modelling...
                   RCON_J    = RCON_J  + SCVDETWEI( GI ) * (LIMT_HAT - LIMT)  &
                        * SUFEN( U_KLOC, GI ) 
                END IF ! For solid modelling...

                DO IDIM = 1, NDIM
                   CT( IDIM, IPHASE, ICOUNT_KLOC2( U_KLOC2 ) ) &
                        = CT( IDIM, IPHASE, ICOUNT_KLOC2( U_KLOC2 ) ) &
                        - RCON_J * UGI_COEF_ELE2_ALL( IDIM, IPHASE, U_KLOC2 ) * CVNORMX_ALL( IDIM, GI )
                END DO
             end if

          END IF
       END DO
    END IF

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
       T_ALL, TOLD_ALL, DEN_ALL, DENOLD_ALL, T2_ALL, T2OLD_ALL, &
       FEMT_ALL, FEMTOLD_ALL, FEMDEN_ALL, FEMDENOLD_ALL, FEMT2_ALL, FEMT2OLD_ALL, USE_FEMT, &
       TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, DENUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL, &
       T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL, &
       IGOT_T2, NPHASE, CV_NONODS,CV_NLOC, X_NLOC,TOTELE, CV_NDGLN, &
       SMALL_FINDRM, SMALL_CENTRM, SMALL_COLM,NSMALL_COLM, &
       X_NDGLN, X_NONODS, NDIM, &
       X_ALL, XC_CV_ALL, IGOT_T_PACK, IGOT_T_CONST, IGOT_T_CONST_VALUE,&
       state, storname, indx)
    ! For the anisotropic limiting scheme we find the upwind values
    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
    ! value for each node pair is stored in the matrices TUPWIND AND
    IMPLICIT NONE
    INTEGER, intent( in ) :: CV_NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, &
         NSMALL_COLM, NDIM,IGOT_T2,NPHASE
    REAL, DIMENSION( :, : ), intent( in ) :: T_ALL,TOLD_ALL,DEN_ALL,DENOLD_ALL
    REAL, DIMENSION( :,:), intent( in ) :: T2_ALL,T2OLD_ALL
    REAL, DIMENSION( :, :), intent( in ) :: FEMT_ALL,FEMTOLD_ALL,FEMDEN_ALL,FEMDENOLD_ALL
    REAL, DIMENSION( :, :), intent( in ) :: FEMT2_ALL,FEMT2OLD_ALL
    LOGICAL, intent( in ) :: USE_FEMT ! Use the FEM solns rather than CV's when interpolating soln
    REAL, DIMENSION( :, : ), intent( inout ) :: TUPWIND_MAT_ALL, TOLDUPWIND_MAT_ALL, &
         DENUPWIND_MAT_ALL, DENOLDUPWIND_MAT_ALL
    REAL, DIMENSION( :, : ), intent( inout ) :: T2UPWIND_MAT_ALL, T2OLDUPWIND_MAT_ALL
    REAL, DIMENSION( :,: ), intent( in ) :: X_ALL
    INTEGER, DIMENSION(: ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    INTEGER, DIMENSION( : ), intent( in) :: SMALL_CENTRM
    LOGICAL, DIMENSION( NPHASE, 6 ), intent( in) :: IGOT_T_PACK, IGOT_T_CONST
    REAL, DIMENSION( NPHASE, 6 ), intent( in) :: IGOT_T_CONST_VALUE
    REAL, DIMENSION( NDIM, CV_NONODS), intent( in ) :: XC_CV_ALL
    type( state_type ), intent( inout ), dimension(:) :: state
    character(len=*), intent(in) :: StorName
    integer, intent(inout) :: indx
! local variables...
    INTEGER :: NFIELD, IMID, NOD, IFIELD
    REAL, DIMENSION( :, : ), ALLOCATABLE :: F_ALL, FEMF_ALL, FUPWIND_MAT_ALL


! Pack all the variables in:
! Always pack in the T & TOLD & T2, T2OLD variables
    NFIELD=(4 + 2*IGOT_T2)*NPHASE
    ALLOCATE( FUPWIND_MAT_ALL(NFIELD,NSMALL_COLM) )
    ALLOCATE( F_ALL(NFIELD, CV_NONODS), FEMF_ALL(NFIELD, CV_NONODS)  )

    F_ALL(1:NPHASE,            :)=T_ALL(1:NPHASE, :)
    F_ALL(1+NPHASE:  2*NPHASE, :)=TOLD_ALL(1:NPHASE, :)
    F_ALL(1+2*NPHASE:3*NPHASE, :)=DEN_ALL(1:NPHASE, :)
    F_ALL(1+3*NPHASE:4*NPHASE, :)=DENOLD_ALL(1:NPHASE, :)

    IF(IGOT_T2.NE.0) THEN
       F_ALL(1+4*NPHASE:5*NPHASE, :)=T2_ALL(1:NPHASE, :)
       F_ALL(1+5*NPHASE:6*NPHASE, :)=T2OLD_ALL(1:NPHASE, :)
    ENDIF
! femf:
    FEMF_ALL(1:NPHASE,            :)=FEMT_ALL(1:NPHASE, :)
    FEMF_ALL(1+NPHASE:  2*NPHASE, :)=FEMTOLD_ALL(1:NPHASE, :)
    FEMF_ALL(1+2*NPHASE:3*NPHASE, :)=FEMDEN_ALL(1:NPHASE, :)
    FEMF_ALL(1+3*NPHASE:4*NPHASE, :)=FEMDENOLD_ALL(1:NPHASE, :)

    IF(IGOT_T2.NE.0) THEN
       FEMF_ALL(1+4*NPHASE:5*NPHASE, :)=FEMT2_ALL(1:NPHASE, :)
       FEMF_ALL(1+5*NPHASE:6*NPHASE, :)=FEMT2OLD_ALL(1:NPHASE, :)
    ENDIF


    !Find upwind field values for limiting
    CALL CALC_ANISOTROP_LIM_VALS( F_ALL, FEMF_ALL, USE_FEMT, FUPWIND_MAT_ALL,  &
    NFIELD,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, SMALL_FINDRM,&
    SMALL_COLM,NSMALL_COLM, X_NDGLN,X_NONODS,NDIM, X_ALL, XC_CV_ALL,&
    state, storname,indx )

! make sure the diagonal is equal to the value:
         DO NOD=1,CV_NONODS
            IMID=SMALL_CENTRM(NOD)
            DO IFIELD=1,NFIELD
               FUPWIND_MAT_ALL( IFIELD, IMID)=F_ALL( IFIELD, NOD)
            END DO
         END DO


    TUPWIND_MAT_ALL(1:NPHASE,      :)=FUPWIND_MAT_ALL(1:NPHASE, :)
    TOLDUPWIND_MAT_ALL(1:NPHASE,   :)=FUPWIND_MAT_ALL(1+NPHASE:2*NPHASE, :)
    DENUPWIND_MAT_ALL(1:NPHASE,    :)=FUPWIND_MAT_ALL(2*NPHASE+1:3*NPHASE, :)
    DENOLDUPWIND_MAT_ALL(1:NPHASE, :)=FUPWIND_MAT_ALL(3*NPHASE+1:4*NPHASE, :)

    IF(IGOT_T2.NE.0) THEN
       T2UPWIND_MAT_ALL(1:NPHASE,    :)=FUPWIND_MAT_ALL(4*NPHASE+1:5*NPHASE, :)
       T2OLDUPWIND_MAT_ALL(1:NPHASE, :)=FUPWIND_MAT_ALL(5*NPHASE+1:6*NPHASE, :)
    ENDIF

    RETURN
    END SUBROUTINE CALC_ANISOTROP_LIM





!##############COMMENTED SUBROUTINE SINCE IT IS NEVER CALLED####################
!  SUBROUTINE CALC_ANISOTROP_LIM_1time(&
!       ! Caculate the upwind values stored in matrix form...
!       T, DEN, T2, &
!       FEMT, FEMDEN, FEMT2, USE_FEMT, &
!       TUPWIND_MAT, DENUPWIND_MAT, &
!       T2UPWIND_MAT, &
!       IGOT_T2, NPHASE, CV_NONODS,CV_NLOC, X_NLOC,TOTELE, CV_NDGLN, &
!       SMALL_FINDRM, SMALL_CENTRM, SMALL_COLM,NSMALL_COLM, &
!       X_NDGLN, X_NONODS, NDIM, &
!       X, Y, Z, XC_CV, YC_CV, ZC_CV)
!    ! For the anisotropic limiting scheme we find the upwind values
!    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
!    ! value for each node pair is stored in the matrices TUPWIND AND
!    IMPLICIT NONE
!    INTEGER, intent( in ) :: CV_NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, &
!         NSMALL_COLM, NDIM,IGOT_T2,NPHASE
!    REAL, DIMENSION( : ), intent( in ) :: T,DEN
!    REAL, DIMENSION( :), intent( in ) :: T2
!    REAL, DIMENSION( :), intent( in ) :: FEMT,FEMDEN
!    REAL, DIMENSION( :), intent( in ) :: FEMT2
!    LOGICAL, intent( in ) :: USE_FEMT ! Use the FEM solns rather than CV's when interpolating soln
!    REAL, DIMENSION( : ), intent( inout ) :: TUPWIND_MAT, &
!         DENUPWIND_MAT
!    REAL, DIMENSION( : ), intent( inout ) :: T2UPWIND_MAT
!    INTEGER, DIMENSION(: ), intent( in ) :: X_NDGLN
!    INTEGER, DIMENSION( : ), intent( in ) :: CV_NDGLN
!    INTEGER, DIMENSION( : ), intent( in) :: SMALL_FINDRM
!    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
!    INTEGER, DIMENSION( : ), intent( in) :: SMALL_CENTRM
!    REAL, DIMENSION( : ), intent( in ) :: X,Y,Z
!    REAL, DIMENSION( : ), intent( in ) :: XC_CV, YC_CV, ZC_CV
!    REAL, DIMENSION(:), ALLOCATABLE :: SOL
!
!    ! Allocate memory
!    INTEGER :: COUNT, COUNT2, CV_NOD
!
!
!    ! Allocate memory and find upwind field values for limiting...
!    IF(IGOT_T2.NE.0) THEN
!
!       allocate( sol( 3*nsmall_colm*nphase) )
!       sol = (/TUPWIND_MAT, DENUPWIND_MAT, T2UPWIND_MAT/)
!
!       ! Obtain the weights
!       CALL CALC_ANISOTROP_LIM_VALS( &
!            ! Caculate the upwind values stored in matrix form...
!            (/T,DEN,T2/), &
!            (/FEMT,FEMDEN,FEMT2/), USE_FEMT, &
!            SOL,  &
!            NPHASE*3,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
!            SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
!            X_NDGLN,X_NONODS,NDIM, &
!            X,Y,Z, XC_CV, YC_CV, ZC_CV )
!
!       TUPWIND_MAT = sol( 1 : nsmall_colm*nphase )
!       DENUPWIND_MAT = sol( 1+nsmall_colm*nphase : 2*nsmall_colm*nphase )
!       T2UPWIND_MAT = sol( 1+2*nsmall_colm*nphase : 3*nsmall_colm*nphase )
!
!       deallocate( sol )
!
!    ELSE
!
!       allocate( sol( 2*nsmall_colm*nphase ) )
!       sol = (/TUPWIND_MAT, DENUPWIND_MAT/)
!
!       CALL CALC_ANISOTROP_LIM_VALS( &
!            ! Caculate the upwind values stored in matrix form...
!            (/T,DEN/),&
!            (/FEMT,FEMDEN/), USE_FEMT, &
!            SOL,  &
!            NPHASE*2,CV_NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
!            SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
!            X_NDGLN,X_NONODS,NDIM, &
!            X,Y,Z, XC_CV, YC_CV, ZC_CV )
!
!       TUPWIND_MAT = sol( 1 : nsmall_colm*nphase )
!       DENUPWIND_MAT = sol( 1+nsmall_colm*nphase : 2*nsmall_colm*nphase )
!
!       deallocate( sol )
!
!    ENDIF
!
!  END SUBROUTINE CALC_ANISOTROP_LIM_1time





  SUBROUTINE CALC_ANISOTROP_LIM_VALS( &
       ! Caculate the upwind values stored in matrix form...
       T_ALL, &
       FEMT_ALL, USE_FEMT, &
       TUPWIND_ALL, &
       NFIELD,NONODS,CV_NLOC,X_NLOC,TOTELE,CV_NDGLN, &
       SMALL_FINDRM,SMALL_COLM,NSMALL_COLM, &
       X_NDGLN,X_NONODS,NDIM, &
       X_ALL, XC_CV_ALL,&
       state, storname, indx )
    ! For the anisotropic limiting scheme we find the upwind values
    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
    ! value for each node pair is stored in the matrices TUPWIND AND
    IMPLICIT NONE
    INTEGER, intent(in) :: NONODS,X_NONODS,TOTELE,CV_NLOC, X_NLOC, NSMALL_COLM, NFIELD,NDIM
    REAL, DIMENSION( :, : ), intent( in ) :: T_ALL
    REAL, DIMENSION( :, : ), intent( in ) :: FEMT_ALL
    LOGICAL, intent( in ) :: USE_FEMT
    REAL, DIMENSION( :, : ), intent( inout ) :: TUPWIND_ALL
    INTEGER, DIMENSION( :  ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( :  ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_FINDRM
    INTEGER, DIMENSION( : ), intent( in ) :: SMALL_COLM
    REAL, DIMENSION( :, : ), intent( in ) :: X_ALL
    REAL, DIMENSION( NDIM, NONODS ), intent( in ) :: XC_CV_ALL
    type( state_type ), intent( inout ), dimension(:) :: state
    character(len=*), intent(in) :: StorName
    integer, intent(inout) :: indx
    ! the centre of each CV is: XC_CV, YC_CV, ZC_CV

    ! Allocate memory for the interpolated upwind values
    real, dimension( :, : ), allocatable :: N!, NLX, NLY, NLZ
    real, dimension (:, :, :), allocatable :: NLX_ALL
    real, dimension( : ), allocatable :: WEIGHT, L1, L2, L3, L4
    integer, dimension( : ), allocatable :: SUB_NDGLNO, SUB_XNDGLNO, ndgln_p2top1
    INTEGER :: COUNT, COUNT2, NOD, SUB_TOTELE, NGI,NLOC, ELE, IL_LOC, IQ_LOC, &
         LOC_ELE, SUB_ELE, SUB_LIN_TOTELE, IMID


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
    ALLOCATE( N(NLOC,NGI))
    ALLOCATE( WEIGHT(NGI) )
    ALLOCATE( L1(NGI), L2(NGI), L3(NGI), L4(NGI) )
    allocate(NLX_ALL(size(X_ALL,1), NLOC, NGI))
    !
    ! Shape functions for triangles and tets...
    CALL TRIQUAold( L1, L2, L3, L4, WEIGHT, ndim==3, NGI )
    ! Work out the shape functions and there derivatives...
    call SHATRInew(L1, L2, L3, L4, WEIGHT,  NLOC,NGI,  N,NLX_ALL)
!    CALL SHATRIold( L1, L2, L3, L4, WEIGHT, ndim==3, &
!         &          NLOC,NGI,&
!         &          N,NLX,NLY,NLZ)

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
         T_ALL, &
         FEMT_ALL, USE_FEMT, &
         TUPWIND_ALL,  &
         NFIELD, NONODS, NLOC, NGI, SUB_TOTELE, SUB_NDGLNO, &
         SMALL_FINDRM,SMALL_COLM, NSMALL_COLM, &
         SUB_XNDGLNO, X_NONODS, NDIM, &
         X_ALL, XC_CV_ALL, &
         N, NLX_ALL, WEIGHT,&
         state, storname, indx )


!    DEALLOCATE( N, NLX, NLY, NLZ, L1, L2, L3, L4, &
!         WEIGHT, SUB_NDGLNO, SUB_XNDGLNO )
    DEALLOCATE( N, NLX_ALL, L1, L2, L3, L4, &
         WEIGHT, SUB_NDGLNO, SUB_XNDGLNO )
    RETURN
  END SUBROUTINE CALC_ANISOTROP_LIM_VALS


  SUBROUTINE CALC_ANISOTROP_LIM_VALS2( &
       ! Caculate the upwind values stored in matrix form...
       T_ALL, &
       FEMT_ALL, USE_FEMT, &
       TUPWIND_ALL,  &
       NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
       FINDRM,COLM,NCOLM, &
       X_NDGLN,X_NONODS,NDIM, &
       X_ALL, XC_CV_ALL,  &
       N,NLX_ALL, WEIGHT,&
       state, storname, indx )
    ! For the anisotropic limiting scheme we find the upwind values
    ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
    ! value for each node pair is stored in the matrices TUPWIND AND
    IMPLICIT NONE
    INTEGER, intent(in) :: NONODS,X_NONODS,TOTELE,NLOC,NGI,NCOLM,NFIELD,NDIM
    REAL, DIMENSION( :,: ), intent( in ) :: T_ALL
    REAL, DIMENSION(  :,: ), intent( in ) :: FEMT_ALL
    LOGICAL, intent( in ) :: USE_FEMT
    REAL, DIMENSION( :,:  ), intent( inout ) :: TUPWIND_ALL
    INTEGER, DIMENSION( : ), INTENT(IN) :: NDGLNO,X_NDGLN
    INTEGER, DIMENSION( : ), INTENT(IN) :: FINDRM,COLM

    REAL, DIMENSION(:,:), intent( in ) :: X_ALL
    REAL, DIMENSION( NDIM, NONODS ), intent( in ) :: XC_CV_ALL
    REAL, DIMENSION(NLOC,NGI), INTENT(IN) :: N!,NLX,NLY,NLZ
    REAL, DIMENSION(:,:,:), INTENT(IN) :: NLX_ALL!DIMENSION(NDIM, NLOC,NGI)
    REAL, DIMENSION(NGI), INTENT(IN) :: WEIGHT
    type( state_type ), intent( inout ), dimension(:) :: state
    character(len=*), intent(in) :: StorName
    integer, intent(inout) :: indx
    !Local variables

    INTEGER, DIMENSION( : ), ALLOCATABLE, SAVE :: ELEMATPSI
    REAL, DIMENSION( :  ), ALLOCATABLE, SAVE :: ELEMATWEI
    LOGICAL, SAVE :: STORE_ELE=.TRUE., RET_STORE_ELE=.FALSE.

    LOGICAL :: D3,DCYL
    ! Allocate memory for the interpolated upwind values
    LOGICAL, PARAMETER :: BOUND  = .TRUE., REFLECT = .false. ! limiting options
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

       CALL FINPTSSTORE(T_ALL,FEMT_ALL,USE_FEMT,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
            TUPWIND_ALL,FINDRM,COLM,NCOLM,NDIM, &
            X_NDGLN,X_NONODS, &
            X_ALL, XC_CV_ALL, &
            N,NLX_ALL, WEIGHT, &
            NOD_FINDELE,NOD_COLELE,NCOLEL, &
            ELEMATPSI,ELEMATWEI,1, &
            BOUND, REFLECT, state, storName, indx)

    ELSE IF( RET_STORE_ELE ) THEN 
       
       ! Find the weights for the interpolation
       ! This does depend on the solns T when BOUND...
       CALL GETSTOREELEWEI(T_ALL,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
            TUPWIND_ALL,FINDRM,COLM,NCOLM,BOUND, &
            ELEMATPSI,ELEMATWEI)

    ELSE 

       ! Assume we have not stored anything (elements or weights)...
       ALLOCATE(DUMMYINT(NCOLM))
       ALLOCATE(DUMMYREAL(NCOLM*NLOC))

       CALL FINPTSSTORE(T_ALL,FEMT_ALL,USE_FEMT,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
            TUPWIND_ALL,FINDRM,COLM,NCOLM,NDIM, &
            X_NDGLN,X_NONODS, &
            X_ALL, XC_CV_ALL, &
            N,NLX_ALL, WEIGHT, &
            NOD_FINDELE,NOD_COLELE,NCOLEL, &
            DUMMYINT,DUMMYREAL,0, &
            BOUND, REFLECT, state, storname, indx)

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





  SUBROUTINE GETSTOREELEWEI(PSI_ALL,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
       &     MATPSI_ALL,FINDRM,COLM,NCOLM,BOUND,&
       &     ELEMATPSI,ELEMATWEI)
    ! use the stored interpolation coeffs to caclulate MATPSI.
    !     This sub finds the matrix values MATPSI for a given point on the 
    !     stencil 
    IMPLICIT NONE
    REAL FRALINE
    LOGICAL BOUND
    PARAMETER(FRALINE=0.001)
    INTEGER, intent(in) :: NFIELD,NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
    REAL, DIMENSION(:,:), INTENT(IN) :: PSI_ALL
    INTEGER, INTENT(IN) :: NCOLM
    INTEGER, DIMENSION(:), INTENT(IN) :: FINDRM
    INTEGER, DIMENSION(:), INTENT(IN) :: COLM
    REAL, DIMENSION(:,:), INTENT(INOUT) :: MATPSI_ALL
    INTEGER, DIMENSION(:), INTENT(IN) :: ELEMATPSI
    REAL, DIMENSION(NCOLM*NLOC),  INTENT(IN) ::  ELEMATWEI
    !  LOCAL VARIABLES...
    INTEGER NOD,COUNT,ELEWIC,ILOC,INOD,IFIELD
    INTEGER KNOD,COUNT2,JNOD
    REAL RMATPSI
    REAL, DIMENSION(NFIELD, TOTELE)::MINPSI
    REAL, DIMENSION(NFIELD, TOTELE)::MAXPSI

    if ( bound ) then

       ! find the max and min local to each element...
       CALL MINMAXELEWIC( PSI_ALL,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
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
                   RMATPSI = RMATPSI + ELEMATWEI( (COUNT-1)*NLOC+ILOC) * PSI_ALL(IFIELD,INOD)
                END DO

                RMATPSI   =PSI_ALL(IFIELD,NOD)   &
                     +(1./FRALINE)*(RMATPSI   -PSI_ALL(IFIELD,NOD))

                ! make locally bounded...
                if ( bound ) then
                   MATPSI_ALL(IFIELD, COUNT)   &
                        =MAX(MIN(RMATPSI,   MAXPSI(IFIELD, ELEWIC)),   &
                        &                            MINPSI(IFIELD, ELEWIC))
                else
                   MATPSI_ALL(IFIELD, COUNT)   =RMATPSI
                end if
             END DO
          END IF
       END DO
    END DO

!    if ( bound ) then
!       DEALLOCATE( MINPSI, MAXPSI )
!    end if

    RETURN

  end subroutine getstoreelewei

  SUBROUTINE MINMAXELEWIC(PSI_ALL,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
       &     FINDRM,COLM,NCOLM,&
       &     MINPSI,MAXPSI)
    ! This sub calculates the max and min values of PSI in local vacinity of 
    ! an element. 
    IMPLICIT NONE
    INTEGER, intent(in) :: NFIELD,NONODS,NLOC,TOTELE,NDGLNO(TOTELE*NLOC)
    REAL, DIMENSION(:,:), INTENT(IN) :: PSI_ALL
    INTEGER, INTENT(IN) :: NCOLM
    INTEGER, INTENT(IN) :: FINDRM(NONODS+1),COLM(NCOLM)
!    REAL, INTENT(INOUT) :: MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
    REAL, DIMENSION(:,:), INTENT(INOUT) :: MINPSI,MAXPSI
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
             !DO IFIELD = 1, NFIELD
                MINPSI( :, ELEWIC )  &
                     = MIN( PSI_ALL(:, JNOD), MINPSI(:, ELEWIC) )
                MAXPSI( :, ELEWIC )  &
                     = MAX( PSI_ALL(:, JNOD), MAXPSI(:, ELEWIC) )
!                     = MAX( PSI_ALL(JNOD+(IFIELD-1)*NONODS), MAXPSI(ELEWIC+(IFIELD-1)*TOTELE) )
             !END DO
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
  SUBROUTINE FINPTSSTORE(PSI_ALL,FEMPSI_ALL,USE_FEMPSI,NFIELD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
       &     MATPSI_ALL,FINDRM,COLM,NCOLM,NDIM, &
       &     X_NDGLN,X_NONODS, &
       &     X_ALL, XC_CV_ALL, &
       &     N,NLX_ALL, WEIGHT,&
       !     work space...
       &     FINDELE,COLELE,NCOLEL,&
       &     ELEMATPSI,ELEMATWEI,IGETSTOR,&
       &     BOUND, REFLECT,&
       state, storname,indx)
    !     This sub finds the matrix values MATPSI for a given point on the 
    !     stencil 
    ! IF IGETSTOR=1 then get ELEMATPSI,ELEMATWEI.
    IMPLICIT NONE
    LOGICAL BOUND,REFLECT
    ! IF REFLECT then use a reflection condition at boundary to 
    ! do limiting. 
    INTEGER, intent(in) :: NFIELD,NONODS,NLOC,NGI,TOTELE,NDIM,X_NONODS
    INTEGER, dimension(TOTELE*NLOC),intent(in) :: NDGLNO
    REAL, dimension(:,:), intent(in) :: PSI_ALL
    REAL, dimension(:,:), intent(in) :: FEMPSI_ALL
    LOGICAL, intent(in) :: USE_FEMPSI
    INTEGER, intent(in) :: NCOLM,NCOLEL
    INTEGER, dimension(NONODS+1), intent(in) :: FINDRM
    INTEGER, dimension(NCOLM),intent(in) :: COLM
    REAL, dimension(:,:), intent(inout) :: MATPSI_ALL
    INTEGER, dimension(TOTELE*NLOC),  intent(in) :: X_NDGLN
    REAL, dimension(:,:), intent(in) :: X_ALL
    REAL, DIMENSION( NDIM, NONODS ), intent( in ) :: XC_CV_ALL
    REAL, dimension(NLOC,NGI), intent(in) :: N!,NLX,NLY,NLZ
    REAL, dimension(:, :,:), intent(in) :: NLX_ALL!dimension(NDIM, NLOC,NGI)
    REAL, dimension(:), intent(in) :: WEIGHT!dimenson(NGI)
    type( state_type ), intent( inout ), dimension(:) :: state
    character(len=*), intent(in) :: StorName
    integer, intent(inout) :: indx
    !     work space...
    INTEGER, dimension(X_NONODS+1),intent(in) :: FINDELE
    INTEGER, dimension(NCOLEL),intent(in) :: COLELE
    INTEGER, intent(in) :: IGETSTOR
    INTEGER, dimension(NCOLM*IGETSTOR), intent(inout) :: ELEMATPSI
    REAL, dimension(NCOLM*NLOC*IGETSTOR), intent(inout) :: ELEMATWEI
    ! ELEWIC is the element to do interpolation from
    ! LOCCORDSK contains the weights. 
    !     Local variables...
    INTEGER NOD,COUNT,NODI,NODJ,ILOC,GI,ELE
    INTEGER ELEWIC,XNOD,XNODJ,IFIELD
    REAL LOCCORDSK(NLOC)
    REAL INVH,LENG
    !     work space...
    real, pointer :: volume
    REAL,  dimension(:), pointer :: DETWEI,RA!dimension(NGI)
    real, dimension (size(X_ALL,1)) :: NORMX1_ALL
    real, dimension (:, :, :), pointer :: NX_ALL ! dimension (size(X_ALL,1), NLOC, NGI)
    real, dimension (size(X_ALL,1), NONODS) :: NORMX_ALL
    REAL, ALLOCATABLE, DIMENSION(:)::MLUM
!    REAL, ALLOCATABLE, DIMENSION(:)::MINPSI
!    REAL, ALLOCATABLE, DIMENSION(:)::MAXPSI
    REAL, ALLOCATABLE, DIMENSION(:,:)::MINPSI
    REAL, ALLOCATABLE, DIMENSION(:,:)::MAXPSI
    INTEGER, ALLOCATABLE, DIMENSION(:)::NOD2XNOD
    real, dimension(size(X_ALL,1)) :: X1_ALL, X2_ALL


    NORMX1_ALL=0.0
!    NORMY1=0.0
!    NORMZ1=0.0
    IF(REFLECT) THEN
       !     calculate normals...********************
       ALLOCATE(MLUM(NONODS))
       NORMX_ALL = 0
       MLUM(1:NONODS) = 0.0
       DO ELE=1,TOTELE! Was loop 

           call DETNLXR_plus_storage( ELE, X_ALL, X_NDGLN, TOTELE, X_NONODS, NLOC, NGI, &
           N, NLX_ALL, WEIGHT, DETWEI, RA, VOLUME, .false., &
           NX_ALL, state, StorName , indx)

          !     
          DO ILOC=1,NLOC! Was loop 
             NODI=NDGLNO((ELE-1)*NLOC+ILOC)
             DO GI=1,NGI! Was loop 
                NORMX_ALL(:,NODI) = NORMX_ALL(:,NODI) + NX_ALL(:,ILOC,GI) * DETWEI(GI)
                MLUM(NODI) =MLUM(NODI) +N(ILOC,GI) *DETWEI(GI)
             END DO
          END DO
       END DO
       !     Renormalise
       DO NODI=1,NONODS! Was loop 
          INVH = SUM(ABS(NORMX_ALL(:,NODI)))/MLUM(NODI)

!          INVH=(ABS(NORMX_ALL(1,NODI))+ABS(NORMX_ALL(2,NODI))+ABS(NORMX_ALL(3,NODI)))&
!               &          /MLUM(NODI)
          IF(INVH.GT.1.E-5) THEN
             LENG = sqrt(dot_product(NORMX_ALL(:,NODI),NORMX_ALL(:,NODI)))
             NORMX_ALL(:,NODI) = NORMX_ALL(:,NODI) / LENG
          ELSE
             NORMX_ALL(:,NODI) = 0.0
          END IF
       END DO
    ENDIF

!    ALLOCATE(MINPSI(TOTELE*NFIELD))
!    ALLOCATE(MAXPSI(TOTELE*NFIELD))
    ALLOCATE(MINPSI(NFIELD, TOTELE))
    ALLOCATE(MAXPSI(NFIELD, TOTELE))

    IF(BOUND) THEN
       ! find the max and min local to each element...
       CALL MINMAXELEWIC(PSI_ALL,NFIELD,NONODS,NLOC,TOTELE,NDGLNO, &
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

    MATPSI_ALL=0.
    DO NOD=1,NONODS! Was loop 10

       XNOD=NOD2XNOD(NOD)
       !     
       DO COUNT=FINDRM(NOD ),FINDRM(NOD+1)-1! Was loop 20

          NODJ=COLM(COUNT)
          XNODJ=NOD2XNOD(NODJ)
          !     
          IF(NOD.NE.NODJ) THEN

             IF(REFLECT) THEN
                NORMX1_ALL = NORMX_ALL(:,NOD)
             ENDIF
             IF(NONODS.NE.X_NONODS) THEN ! Its a DG soln field...
               X1_ALL = XC_CV_ALL(:,NOD)

                X2_ALL = XC_CV_ALL(:, NODJ)
             ELSE
               X1_ALL = X_ALL(:,XNOD)
               X2_ALL = X_ALL(:,XNODJ)
             ENDIF
             CALL MATPTSSTORE(MATPSI_ALL,COUNT,NFIELD,NOD,XNOD,&
                  &              PSI_ALL,FEMPSI_ALL,USE_FEMPSI,NONODS,X_NONODS,&
                  &              NLOC,TOTELE,X_NDGLN,NDGLNO,&
                  &              NCOLM,&
                  &              X1_ALL,&
                  &              X2_ALL,&
                  &              NORMX1_ALL,&
                  &              X_ALL,&
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
      SUBROUTINE MATPTSSTORE(MATPSI_ALL,COUNT,NFIELD,NOD,XNOD,&
     &     PSI_ALL,FEMPSI_ALL,USE_FEMPSI,NONODS,X_NONODS,&
     &     NLOC,TOTELE,X_NDGLN,NDGLNO,&
     &     NCOLM,&
     &     X1_ALL,&
     &     X2_ALL,&
     &     NORMX1_ALL,&
     &     X_ALL,&
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
      REAL, dimension(:,:), intent(in) :: PSI_ALL
      REAL, dimension(:,:), intent(in) :: FEMPSI_ALL
      LOGICAL, intent(in) :: USE_FEMPSI
      REAL, dimension(:,:), intent(inout) :: MATPSI_ALL
      INTEGER, intent(in) :: X_NDGLN(NLOC*TOTELE),NDGLNO(NLOC*TOTELE)
      INTEGER, intent(in) :: NCOLM
!      REAL, intent(in) :: X1,Y1,Z1,X2,Y2,Z2,NORMX1,NORMY1,NORMZ1
      real, dimension(:) :: X1_ALL, X2_ALL, NORMX1_ALL!dimension(NDIM)
      REAL, dimension(:,:), intent(in) :: X_ALL
      INTEGER, intent(in) :: NCOLEL
      INTEGER, intent(in) :: FINDELE(X_NONODS+1),COLELE(NCOLEL)
!      REAL, intent(in) :: MINPSI(TOTELE*NFIELD),MAXPSI(TOTELE*NFIELD)
      REAL, DIMENSION(:, :), intent(in) :: MINPSI,MAXPSI
      INTEGER, intent(inout) :: ELEWIC
      REAL, intent(inout) :: LOCCORDSK(NLOC)
!     
!     Local variables...
      REAL, dimension(4):: LOCCORDS
      INTEGER , dimension(4) :: LOCNODS,LOCNODSK
      INTEGER, dimension(4) :: NLOCNODS,NLOCNODSK
      INTEGER :: ELE,ILOC,KNOD,JNOD,IFIELD, COUNT2
      REAL :: MINCOR,MINCORK,RSUM
      REAL :: DIST12,RN,RMATPSI
      REAL :: FRALINE
      LOGICAL IS_DG
      !The dimension of the variables below should be NDIM, however, due to cross products
      !we need three dimensions
      REAL, dimension(3) ::  XC_ALL, VX_ALL, REFX_ALL, REFX2_ALL, T2X_ALL, T1X_ALL, AUXNORMX1_ALL


      IS_DG=NONODS.NE.X_NONODS
!     
      FRALINE=FRALINE2
      IF(IS_DG) FRALINE=1.0
      XC_ALL = 0.
      XC_ALL(1:NDIM) = X1_ALL - FRALINE*(X2_ALL-X1_ALL)
!print *, "XC_ALL before reflect", XC_ALL

      IF(REFLECT) THEN
         IF(SUM(ABS(NORMX1_ALL)).NE.0.0) THEN
!  if (XC,YC,ZC) is outside the domain
!     The rotation matrix in 3-D is R=  
!     NX    NY    NZ
!     T1X   T1Y   T1Z
!     T2X   T2Y   T2Z
!
            VX_ALL = 0.
            VX_ALL(1:NDIM) = X1_ALL - X2_ALL
!     
            AUXNORMX1_ALL = 0.
            AUXNORMX1_ALL(1:NDIM) = NORMX1_ALL
            CALL XPROD(T2X_ALL, AUXNORMX1_ALL, VX_ALL)
!     
            !DIST12=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
            DIST12 = SQRT(DOT_PRODUCT(X1_ALL-X2_ALL,X1_ALL-X2_ALL))
            RN = SQRT(DOT_PRODUCT(T2X_ALL,T2X_ALL))
            IF(RN.LT.(1.E-5)*DIST12) THEN
!     Simply have VX,VY,VZ going in the opposite direction...
                XC_ALL(1:NDIM) = X1_ALL - VX_ALL(1:NDIM)*FRALINE
            ELSE
                T2X_ALL = T2X_ALL/RN
!     T1=Nx (-T2)
               CALL XPROD(T1X_ALL, AUXNORMX1_ALL, -T2X_ALL)
!     
               REFX2_ALL(1) = SUM(NORMX1_ALL(1:NDIM)*VX_ALL(1:NDIM))
               REFX2_ALL(2) = SUM(T1X_ALL(1:NDIM)*VX_ALL(1:NDIM))


!     Reflect...
                REFX2_ALL(1) = - REFX2_ALL(1)
!     MAP BACK USING R^T

!     (REFX,REFY,REFZ) is the reflected direction...
                REFX_ALL(1) =  NORMX1_ALL(1) * REFX2_ALL(1) + T1X_ALL(1) * REFX2_ALL(2)
                REFX_ALL(2) =  NORMX1_ALL(2) * REFX2_ALL(1) + T1X_ALL(2) * REFX2_ALL(2)
                IF (NDIM==3) THEN
                !SOME MORE THINGS NEED TO BE ADDED
                    REFX2_ALL(3) = SUM(T2X_ALL(:)*VX_ALL)

                    REFX_ALL(1) = REFX_ALL(1) + T2X_ALL(1) * REFX2_ALL(3)
                    REFX_ALL(2) = REFX_ALL(2) + T2X_ALL(2) * REFX2_ALL(3)
                    REFX_ALL(3) =  NORMX1_ALL(3) * REFX2_ALL(1) + T1X_ALL(3) * REFX2_ALL(2)+ T2X_ALL(3) * REFX2_ALL(3)
                END IF

                XC_ALL(1:NDIM) = X_ALL(:,1) + REFX_ALL(1:NDIM)*FRALINE
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
         
!     
         LOCNODS(1)=X_NDGLN((ELE-1)*NLOC+1)
         LOCNODS(2)=X_NDGLN((ELE-1)*NLOC+2)
         LOCNODS(3)=X_NDGLN((ELE-1)*NLOC+3)
!     
! Calculate the local coord but with 4th point replaced by INOD...
! Find local coords LOCCORDS of point INOD corresponding to these nodes LOCNODS...
 
         IF (NDIM==3) THEN
            !Two coordinates missing if 3D
            NLOCNODS(4)=NDGLNO((ELE-1)*NLOC+4)
            LOCNODS(4)=X_NDGLN((ELE-1)*NLOC+4)
         
            CALL TRILOCCORDS(XC_ALL(1),XC_ALL(2),XC_ALL(3), &
                 &        LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),LOCCORDS(4),&
                 !     The 4 corners of the tet...
                 &        X_ALL(1,LOCNODS(1)),X_ALL(2,LOCNODS(1)),X_ALL(3,LOCNODS(1)),&
                 &        X_ALL(1,LOCNODS(2)),X_ALL(2,LOCNODS(2)),X_ALL(3,LOCNODS(2)),&
                 &        X_ALL(1,LOCNODS(3)),X_ALL(2,LOCNODS(3)),X_ALL(3,LOCNODS(3)),&
                 &        X_ALL(1,LOCNODS(4)),X_ALL(2,LOCNODS(4)),X_ALL(3,LOCNODS(4)) )
         ELSE
            CALL TRILOCCORDS2D(XC_ALL(1),XC_ALL(2), &
                 &        LOCCORDS(1),LOCCORDS(2),LOCCORDS(3),&
                 !     The 3 corners of the tri...
                 &        X_ALL(1,LOCNODS(1)),X_ALL(2,LOCNODS(1)),&
                 &        X_ALL(1,LOCNODS(2)),X_ALL(2,LOCNODS(2)),&
                 &        X_ALL(1,LOCNODS(3)),X_ALL(2,LOCNODS(3)) )
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
               RMATPSI   =RMATPSI  +LOCCORDSK(ILOC)*FEMPSI_ALL(IFIELD,NLOCNODSK(ILOC))
            ELSE
               RMATPSI   =RMATPSI  +LOCCORDSK(ILOC)*PSI_ALL(IFIELD, NLOCNODSK(ILOC))
            ENDIF
!         XC=XC+LOCCORDSK(ILOC)*X(LOCNODSK(ILOC))
!         YC=YC+LOCCORDSK(ILOC)*Y(LOCNODSK(ILOC))
!         ZC=ZC+LOCCORDSK(ILOC)*Z(LOCNODSK(ILOC))
         END DO
!     Exaduate difference by a factor of 100.
         IF(USE_FEMPSI) THEN
            RMATPSI   = FEMPSI_ALL(IFIELD,  NOD )  &
         + (1./FRALINE) * ( RMATPSI - FEMPSI_ALL( IFIELD, NOD) )
         ELSE
            RMATPSI   = PSI_ALL( IFIELD, NOD )  &
         + (1./FRALINE) * ( RMATPSI - PSI_ALL(IFIELD,  NOD) )
         ENDIF

!     Now correct to make sure that we get a bounded soln...
         IF(BOUND) THEN
           RMATPSI   =MAX(MIN(RMATPSI,   MAXPSI(IFIELD, ELEWIC)),   MINPSI(IFIELD, ELEWIC))
         ENDIF
         MATPSI_ALL(IFIELD, COUNT)   =RMATPSI
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


   function tetvolume(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3)
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




  SUBROUTINE GET_INT_VEL_NEW( NPHASE, NDOTQNEW, NDOTQ,INCOME, &
       HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS,  &
       LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
       LOC_U, LOC2_U, LOC_NU, LOC2_NU, SLOC_NU, &
       CV_NODI, CV_NODJ, CVNORMX_ALL,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
       SUF_SIG_DIAGTEN_BC,  &
       UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
       VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS, &
       MASS_CV_I, MASS_CV_J, NDIM, MAT_NLOC, MAT_NONODS, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       IANISOTROPIC, &
       NUGI_ALL, TUPWIND_IN, TUPWIND_OUT)
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ, CV_NODI, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, CV_ELE_TYPE, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, &
         NDIM, MAT_NLOC, MAT_NONODS,  &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, MASS_CV_I, MASS_CV_J, VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2
    REAL, DIMENSION( : ), intent( inout ) :: NDOTQNEW,NDOTQ, INCOME
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( :, :, : ), intent( in ) :: WIC_U_BC_ALL
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :, : ), intent( in ) :: LOC_FEMT, LOC2_FEMT
    REAL, DIMENSION( : ), intent( in ) :: LOC_T_I, LOC_T_J, LOC_DEN_I, LOC_DEN_J
    REAL, DIMENSION( :, :, : ), intent( in ) :: LOC_U, LOC2_U, LOC_NU, LOC2_NU, SLOC_NU
    REAL, DIMENSION( :, : ), intent( in ) :: CVNORMX_ALL
    REAL, DIMENSION( :, :, : ), intent( in ) :: SUF_U_BC_ALL
    REAL, DIMENSION( :, : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
    REAL, DIMENSION( :, :, : ), intent( inout ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
    REAL, DIMENSION( :, : ), intent( inout ) :: NUGI_ALL
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    REAL, DIMENSION( :, :, : ), intent( in ) :: VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS
    REAL, DIMENSION( :, :, : ), intent( in ) :: INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS

    INTEGER, intent( in ) :: IANISOTROPIC
    REAL, DIMENSION( NPHASE ), intent( in ) :: TUPWIND_IN, TUPWIND_OUT
    ! local variables
!    LOGICAL, PARAMETER :: POROUS_VEL = .false. ! For reduced variable porous media treatment.
    INTEGER :: U_NLOC_LEV,U_KLOC_LEV,U_KLOC,U_NODK_IPHA, U_KLOC2, U_NODK2_IPHA, IPHASE
!    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS
    REAL, ALLOCATABLE, DIMENSION(:,:) :: NUGI_ALL_OTHER

    logical, parameter :: LIMIT_USE_2ND=.false.

    ! Local variable for indirect addressing
    INTEGER :: IDIM

!     print *,'is_overlapping,is_compact_overlapping:',is_overlapping,is_compact_overlapping
!      stop 21 

!          print *,'CV_DG_VEL_INT_OPT:',CV_DG_VEL_INT_OPT
!           stop 2197

    IF( is_overlapping ) THEN
       ! For overlapping basis function approach.
       CALL GET_INT_VEL_OVERLAP_NEW( NPHASE, NDOTQ, INCOME, &
            HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS,  &
            LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
            LOC_NU, LOC2_NU, SLOC_NU, &
            CV_NODI, CV_NODJ, CVNORMX_ALL,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
            SUF_SIG_DIAGTEN_BC, &
            UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
            VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
            VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
            MASS_CV_I, MASS_CV_J, NDIM, MAT_NLOC, MAT_NONODS, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            IANISOTROPIC, &
            TUPWIND_IN, TUPWIND_OUT)

    ELSE IF( is_compact_overlapping ) THEN
       ! For reduced variable porous media treatment.
       CALL GET_INT_VEL_POROUS_VEL( NPHASE, NDOTQ, INCOME, &
            HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS,  &
            LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
            LOC_NU, LOC2_NU, SLOC_NU, &
            CV_NODI, CV_NODJ, CVNORMX_ALL,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
            SUF_SIG_DIAGTEN_BC, &
            UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
            VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
            VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS, &
       NUGI_ALL, &
            MASS_CV_I, MASS_CV_J, NDIM, MAT_NLOC, MAT_NONODS, &
            IN_ELE_UPWIND, DG_ELE_UPWIND, &
            IANISOTROPIC, &
            TUPWIND_IN, TUPWIND_OUT)
    ELSE
       CALL GET_INT_VEL_ORIG_NEW( NPHASE, NDOTQ, INCOME, &
            HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS,  &
            LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
            LOC_NU, LOC2_NU, SLOC_NU, &
            CV_NODI, CV_NODJ, CVNORMX_ALL,  &
            CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
            SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
            UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
            VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
            NDIM, MAT_NLOC, MAT_NONODS, &
            IN_ELE_UPWIND, DG_ELE_UPWIND )
    END IF 

 
    ! Calculate NDOTQNEW from NDOTQ
    NDOTQNEW = NDOTQ ! initialize it like this so that it contains the b.c's

  IF(.NOT. is_compact_overlapping ) THEN
!  IF(.false.) THEN
    NUGI_ALL=0.0
    DO U_KLOC = 1, U_NLOC
       DO IDIM = 1, NDIM
          NDOTQNEW=NDOTQNEW + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE_ALL(IDIM, :,U_KLOC) & 
                           * ( LOC_U(IDIM,:, U_KLOC ) - LOC_NU(IDIM,:,U_KLOC ) ) * CVNORMX_ALL(IDIM, GI)
       END DO
       NUGI_ALL(:, :) = NUGI_ALL(:, :) + SUFEN( U_KLOC, GI )*LOC_NU(:, :, U_KLOC)
    END DO

    IF( (ELE2 /= 0) .AND. (ELE2 /= ELE) ) THEN
       ALLOCATE(NUGI_ALL_OTHER(NDIM,NPHASE))
       NUGI_ALL(:, :) = 0.5*NUGI_ALL(:, :) ! Reduce by half and take the other half from the other side of element...
       NUGI_ALL_OTHER(:, :) = 0.0
       ! We have a discontinuity between elements so integrate along the face...
       DO U_KLOC = 1, U_NLOC
          U_KLOC2 = U_OTHER_LOC( U_KLOC )
          IF( U_KLOC2 /= 0 ) THEN
             DO IDIM = 1, NDIM
                NDOTQNEW=NDOTQNEW + SUFEN( U_KLOC, GI ) * UGI_COEF_ELE2_ALL(IDIM, :,U_KLOC2) & 
                           * ( LOC2_U(IDIM,:, U_KLOC ) - LOC2_NU(IDIM,:,U_KLOC ) ) * CVNORMX_ALL(IDIM, GI)
             END DO
             NUGI_ALL_OTHER(:, :) = NUGI_ALL_OTHER(:, :) + 0.5*SUFEN( U_KLOC, GI )*LOC2_NU(:, :, U_KLOC)
          END IF
       END DO
       NUGI_ALL(:, :) = NUGI_ALL(:, :) +  NUGI_ALL_OTHER(:, :)
    END IF

! ENDOF IF(.NOT. is_compact_overlapping ) THEN
  ENDIF


    RETURN


contains





  SUBROUTINE GET_INT_VEL_ORIG_NEW( NPHASE, NDOTQ, INCOME, &
       HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS,  &
       LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
       LOC_NU, LOC2_NU, SLOC_NU, &
       CV_NODI, CV_NODJ, CVNORMX_ALL,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC,STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
       UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
       VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       NDIM, MAT_NLOC, MAT_NONODS, &
       IN_ELE_UPWIND, DG_ELE_UPWIND )

    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ, CV_NODI, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NDIM, MAT_NLOC, MAT_NONODS,  &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2
    REAL, DIMENSION( : ), intent( inout ) :: NDOTQ, INCOME
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( :, :, : ), intent( in ) :: WIC_U_BC_ALL
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :, : ), intent( in ) :: LOC_FEMT, LOC2_FEMT
    REAL, DIMENSION( : ), intent( in ) :: LOC_T_I, LOC_T_J, LOC_DEN_I, LOC_DEN_J
    REAL, DIMENSION( :, :, : ), intent( in ) :: LOC_NU, LOC2_NU, SLOC_NU
    REAL, DIMENSION( :, : ), intent( in ) :: CVNORMX_ALL
    REAL, DIMENSION( :, :, : ), intent( in ) :: SUF_U_BC_ALL
    REAL, DIMENSION( :, :, : ), intent( inout ) :: UGI_COEF_ELE_ALL, &
                                             UGI_COEF_ELE2_ALL
    REAL, DIMENSION( :, :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC

    ! Local variables
    REAL :: UDGI,VDGI,WDGI,  &
         UDGI2,VDGI2,WDGI2,  &
         UDGI_INT,VDGI_INT,WDGI_INT
    INTEGER :: U_KLOC,U_NODK,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_NODK2,U_SKLOC, &
         U_SNODK,U_SNODK_IPHA, II, IDIM
    integer, dimension(U_SNLOC) ::  U_NODK_IPHA_V, U_SNODK_IPHA_V
    !REAL, DIMENSION(NDIM) :: CVNORMX_ALL

    ! Local variable for indirect addressing
    REAL, DIMENSION ( NDIM, NPHASE ) :: UDGI_ALL, UDGI2_ALL, UDGI_INT_ALL
    INTEGER :: IPHASE, CV_NODI_IPHA, CV_NODJ_IPHA
    REAL, DIMENSION(NPHASE) :: DT_I,DT_J,NDOTQ_INT


!          stop 2928
    ! coefficients for this element ELE
    UGI_COEF_ELE_ALL = 0.0 

    ! coefficients for this element ELE2
    UGI_COEF_ELE2_ALL = 0.0


    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain. 
       DO IPHASE = 1, NPHASE
          IF( WIC_U_BC_ALL( 1, IPHASE, SELE) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
             UDGI_ALL(:, IPHASE) = 0.0
             DO U_KLOC = 1, U_NLOC
                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * LOC_NU( :, IPHASE, U_KLOC )
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) = 1.0
             END DO
          ELSE ! Specified vel bc.
             UDGI_ALL(:, IPHASE) = 0.0
             DO U_SKLOC = 1, U_SNLOC
                U_KLOC = U_SLOC2LOC( U_SKLOC )
                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * SUF_U_BC_ALL(:, IPHASE, U_SNLOC* (SELE-1) + U_SKLOC )
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) = 0.0  
             END DO
          END IF
       END DO

    ELSE ! Conditional_SELE. Not on the boundary of the domain.
       
       UDGI_ALL = 0.0
       DO U_KLOC = 1, U_NLOC
          UDGI_ALL = UDGI_ALL + SUFEN( U_KLOC, GI ) * LOC_NU( :, :, U_KLOC )
          UGI_COEF_ELE_ALL(:, :, U_KLOC) = 1.0
       END DO

       Conditional_ELE2: IF( ELE2 /= 0 ) THEN
          UDGI2_ALL = 0.0
          DO U_KLOC = 1, U_NLOC
             U_KLOC2 = U_OTHER_LOC( U_KLOC )
             IF ( U_KLOC2 /= 0 ) THEN ! MAKE SURE WE DONT NEED THIS...
                UDGI2_ALL = UDGI2_ALL + SUFEN( U_KLOC, GI ) * LOC2_NU(:, :, U_KLOC)
                UGI_COEF_ELE2_ALL( :, :, U_KLOC2) = 1.0
             ENDIF
          END DO

          IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
             DT_I=1.0
             DT_J=1.0
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
             DT_I=MAX(1.E-2,LOC_T_I)
             DT_J=MAX(1.E-2,LOC_T_J)
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
             DT_I=LOC_DEN_I*LOC_T_I
             DT_J=LOC_DEN_J*LOC_T_J
          ENDIF
          ! Amend weighting for porosity only across elements...
          IF(ABS(CV_DG_VEL_INT_OPT ) >= 2) THEN 
             IF(ELE /= ELE2) THEN 
                DT_I=VOLFRA_PORE_ELE *DT_I
                DT_J=VOLFRA_PORE_ELE2*DT_J
             ENDIF
          ENDIF

          DO IPHASE = 1, NPHASE
             UDGI_INT_ALL(:,IPHASE) = (DT_I(IPHASE) * UDGI_ALL(:,IPHASE) + DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)) &
                                      / (DT_I(IPHASE) + DT_J(IPHASE))
          

             IF( CV_DG_VEL_INT_OPT < 0 ) THEN

                NDOTQ_INT(IPHASE) = DOT_PRODUCT( CVNORMX_ALL(:, GI), UDGI_INT_ALL( :, IPHASE ) )

                IF( NDOTQ_INT(IPHASE) <= 0.0 ) THEN  !Incoming
                !   DT_I=1.0
                   DT_J(IPHASE)=DT_I(IPHASE)+DT_J(IPHASE)
                ELSE
                   DT_I(IPHASE)=DT_I(IPHASE)+DT_J(IPHASE)
                !   DT_J=1.0
                ENDIF

                UDGI_INT_ALL(:,IPHASE) = (DT_I(IPHASE) * UDGI_ALL(:,IPHASE) + &
                                          DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)) / (DT_I(IPHASE) + DT_J(IPHASE))

             ENDIF

             UDGI_ALL( :, IPHASE ) = UDGI_INT_ALL( :, IPHASE )

             UGI_COEF_ELE_ALL(:,IPHASE,:)=DT_I(IPHASE) * UGI_COEF_ELE_ALL(:,IPHASE,:) &
                                          /(DT_I(IPHASE) + DT_J(IPHASE))

             UGI_COEF_ELE2_ALL(:,IPHASE,:)=DT_J(IPHASE) * UGI_COEF_ELE2_ALL(:,IPHASE,:) &
                                          /(DT_I(IPHASE) + DT_J(IPHASE))

          END DO

       ENDIF Conditional_ELE2

    ENDIF Conditional_SELE


    NDOTQ =  MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL)

    ! Define whether flux is incoming or outgoing, depending on direction of flow
    INCOME = 0.5*( 1. + SIGN(1.0, -NDOTQ) )


    RETURN  

  END SUBROUTINE GET_INT_VEL_ORIG_NEW

end SUBROUTINE GET_INT_VEL_NEW




      SUBROUTINE GET_INT_VEL_OVERLAP_NEW( NPHASE, NDOTQ,INCOME, &
       HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
       LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
       LOC_NU, LOC2_NU, SLOC_NU, &
       CV_NODI, CV_NODJ, CVNORMX_ALL,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
       SUF_SIG_DIAGTEN_BC, &
       UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
       VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       MASS_CV_I, MASS_CV_J, NDIM, MAT_NLOC, MAT_NONODS, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       IANISOTROPIC,  &
       TUPWIND_IN, TUPWIND_OUT)
    !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    ! it assumes an overlapping decomposition approach for velocity. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, U_NLOC, CV_SNLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ, CV_NODI, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NDIM, MAT_NLOC, MAT_NONODS,  &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, MASS_CV_I, MASS_CV_J, VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2
    REAL, DIMENSION( : ), intent( inout ) :: NDOTQ, INCOME
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( :, :, : ), intent( in ) :: WIC_U_BC_ALL
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :, : ), intent( in ) :: LOC_FEMT, LOC2_FEMT
    REAL, DIMENSION( : ), intent( in ) :: LOC_T_I, LOC_T_J, LOC_DEN_I, LOC_DEN_J
    REAL, DIMENSION( :, :, : ), intent( in ) ::  LOC_NU, LOC2_NU, SLOC_NU
    REAL, DIMENSION( :, : ), intent( in ) :: CVNORMX_ALL
    REAL, DIMENSION( :, :, : ), intent( in ) :: SUF_U_BC_ALL
    REAL, DIMENSION( : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
    REAL, DIMENSION( :, :, : ), intent( inout ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
    REAL, DIMENSION( : , :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    REAL, DIMENSION( :, :, : ), intent( in ) :: VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS, &
                                                VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS

      INTEGER, intent( in ) :: IANISOTROPIC
      REAL, DIMENSION( NPHASE ), intent( in ) :: TUPWIND_IN, TUPWIND_OUT

    ! Local variables
    REAL :: UDGI,VDGI,WDGI,  &
         UDGI2,VDGI2,WDGI2,  &
         UDGI_INT,VDGI_INT,WDGI_INT, &
         FEMTOLDGI_IPHA, OVER_RELAX, V_NODI, G_NODI, V_NODJ, G_NODJ, &
         GEOMTOLDGI_IPHA, W_UPWIND, W_UPWINDOLD, INCOME4, &
         UPWIND_FRAC, TUPWIN, TUPWI2
    REAL, DIMENSION(NPHASE) :: NDOTQ_INT,NDOTQ2, FEMTGI_IPHA, GEOMTGI_IPHA, &
                               ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                               GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA, &
                               DT_I,DT_J, INCOME3
    REAL :: LIMT3(NPHASE) 
    REAL :: NVEC(NDIM),SUF_SIG_DIAGTEN_BC_GI(NDIM), UGI_TMP(NDIM)
    INTEGER :: U_KLOC,U_NODK,U_NODK2,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_SKLOC, &
         U_SNODK,U_SNODK_IPHA, II,  COUNT, &
         U_KLOC_LEV, U_NLOC_LEV, U_SKLOC_LEV, U_SNLOC_LEV, CV_KLOC, CV_KNOD, &
         CV_KLOC2, CV_KNOD2, IDIM, JDIM, IJ, MAT_NODI, MAT_NODJ, &
         CV_SKLOC, CV_SNODK, CV_SNODK_IPHA, CV_STAR_IPHA, CV_KNOD_IPHA, CV_KNOD2_IPHA, &
         U_NODK3,U_NODK3_IPHA, U_KLOC3
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
! The new simple limiter NEW_LIMITER can produce very slightly different results 
!    LOGICAL, PARAMETER :: NEW_LIMITER = .true.
    LOGICAL, PARAMETER :: NEW_LIMITER = .false.
    REAL :: TMIN_STORE, TMAX_STORE
    REAL, DIMENSION(NPHASE) :: PERM_TILDE,NDOTQ_TILDE, NDOTQ2_TILDE, NDOTQOLD_TILDE, NDOTQOLD2_TILDE, rden_ave, Q_UNDERLY
    REAL, DIMENSION(NPHASE) :: NDOTQ_KEEP_IN,  NDOTQ_KEEP, NDOTQ2_KEEP
    ! coefficients for this element ELE
    real :: gamma,  grad2nd
    real :: abs_tilde_2nd, abs_tildeold_2nd, max_nodtq_keep, min_nodtq_keep
    real :: w_weight, relax

    real, DIMENSION(NPHASE) :: DT_I_upwind, DT_J_upwind
    real :: dt_max, dt_min
    real, DIMENSION(NPHASE) :: abs_tilde1, abs_tilde2, abs_tilde, abs_max, abs_min, w_relax
    real, DIMENSION(NPHASE) :: wrelax, wrelax1, wrelax2
    real :: T_PELE, T_PELEOT, TMIN_PELE, TMAX_PELE, TMIN_PELEOT, TMAX_PELEOT

    ! Local variable for indirect addressing
    REAL, DIMENSION ( NDIM, NPHASE ) :: UDGI_ALL, UDGI2_ALL, UDGI_INT_ALL
    real :: courant_or_minus_one_new(Nphase),XI_LIMIT(Nphase)
    INTEGER :: IPHASE, CV_NODI_IPHA, CV_NODJ_IPHA


    ! print *,'IN_ELE_UPWIND,CV_DG_VEL_INT_OPT,  DG_ELE_UPWIND, IANISOTROPIC:',IN_ELE_UPWIND,CV_DG_VEL_INT_OPT,  DG_ELE_UPWIND, IANISOTROPIC
    !  stop 8721

    ! coefficients for this element ELE
    UGI_COEF_ELE_ALL=0.0

    ! coefficients for this element ELE2
    UGI_COEF_ELE2_ALL=0.0

    U_NLOC_LEV = U_NLOC / CV_NLOC
    U_SNLOC_LEV = U_SNLOC / CV_NLOC


    got_dt_ij=.false.

    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain. 
     DO IPHASE = 1, NPHASE
       IF( WIC_U_BC_ALL( 1, IPHASE, SELE) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
          UDGI_ALL(:, IPHASE) = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * LOC_NU( :, IPHASE, U_KLOC )
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
          UGI_COEF_ELE_ALL(:, IPHASE, :)=0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             IF (DOT_PRODUCT(UDGI_ALL(:, IPHASE), CVNORMX_ALL(:, GI)).LT.0.0) THEN ! Incomming...
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) &
                                                    +1.0*SUF_SIG_DIAGTEN_BC_GI(:)
             ELSE
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)+1.0
             ENDIF
          END DO

          IF(DOT_PRODUCT(UDGI_ALL(:, IPHASE), CVNORMX_ALL(:, GI)).LT.0.0) THEN ! Incomming...
             UGI_TMP = SUF_SIG_DIAGTEN_BC_GI(:) * UDGI_ALL(:, IPHASE)
             UDGI_ALL(:, IPHASE) = UGI_TMP
          ENDIF
       ELSE ! Specified vel bc.

          UDGI_ALL(:, IPHASE) = 0.0
          UGI_COEF_ELE_ALL(:, IPHASE, :) = 0.0 
          DO U_SKLOC_LEV = 1, U_SNLOC_LEV
             U_SKLOC = (CV_ILOC-1)*U_SNLOC_LEV + U_SKLOC_LEV
             U_KLOC = U_SLOC2LOC( U_SKLOC )

             IF(WIC_U_BC_ALL( 1, IPHASE, SELE ) == 10) THEN
                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * 0.5 * &
                                      (SLOC_NU(:, IPHASE, U_SKLOC) + SUF_U_BC_ALL(:, IPHASE, U_SNLOC_LEV* (SELE-1) +U_SKLOC_LEV))

                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) + 0.5
             ELSE

                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI )*SUF_U_BC_ALL(:, IPHASE, U_SNLOC_LEV* (SELE-1) +U_SKLOC_LEV)
             END IF
          END DO
       END IF
     END DO ! PHASE LOOP

    ELSE ! Conditional_SELE. Not on the boundary of the domain.
       Conditional_ELE2: IF(( ELE2 == 0 ).OR.( ELE2 == ELE)) THEN

          UDGI_ALL(:, :) = 0.0

          UDGI2_ALL(:, :) = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC =(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_KLOC2=(CV_JLOC-1)*U_NLOC_LEV + U_KLOC_LEV

             UDGI_ALL(:, :) = UDGI_ALL(:, :) + SUFEN( U_KLOC, GI ) * LOC_NU( :, :, U_KLOC )

             UDGI2_ALL(:, :) = UDGI2_ALL(:, :) + SUFEN( U_KLOC2, GI ) * LOC_NU( :, :, U_KLOC2 )

          END DO

          NDOTQ  = MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )
          NDOTQ2 = MATMUL( CVNORMX_ALL(:, GI), UDGI2_ALL )

          NDOTQ_TILDE = 0.5 * ( NDOTQ +  NDOTQ2 )
          NDOTQ2_TILDE =   NDOTQ_TILDE

          IF((IN_ELE_UPWIND==1).OR.FORCE_UPWIND_VEL) THEN

             WHERE (NDOTQ_TILDE < 0.0) 
                INCOME=1.0
             ELSE WHERE
                INCOME=0.0
             END WHERE

             WHERE (abs(NDOTQ2 - NDOTQ).lt. 1.e-7) 
                   INCOME=0.5
             END WHERE



          ELSE IF(IN_ELE_UPWIND==2) THEN ! the best

             UPWIND_FRAC=0.8
             WHERE (NDOTQ_TILDE < 0.0)
                INCOME=0.8
             ELSE WHERE
                INCOME=0.2 ! MAKE SURE IT IS RIGHT?
             END WHERE

          ELSE IF(IN_ELE_UPWIND==3) THEN ! the best optimal upwind frac.

             !Calculate vel at GI
             FEMTGI_IPHA = matmul(LOC_FEMT, SCVFEN(:,GI) )

             if(IANISOTROPIC==0) then ! this is the only needed for isotropic limiting for velocity...
                FEMTGI_IPHA = ( MASS_CV_J * LOC_T_I + &
                     MASS_CV_I * LOC_T_J ) / (MASS_CV_I+MASS_CV_J)
             endif
             GEOMTGI_IPHA = ( MASS_CV_J * LOC_T_I + &
                  MASS_CV_I * LOC_T_J ) / (MASS_CV_I+MASS_CV_J)

             !We perform: n' * sigma * n
             DO IPHASE = 1, NPHASE
                 ABS_CV_NODI_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
                 GRAD_ABS_CV_NODI_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(GI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
                 ABS_CV_NODJ_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
                 GRAD_ABS_CV_NODJ_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(GJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
             END DO
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
                NDOTQ_TILDE  = ( LOC_DEN_I * LOC_T_I * NDOTQ -  &
                     &           LOC_DEN_J * LOC_T_J * NDOTQ2 ) & 
                     / vtolfun( VOLFRA_PORE_ELE*LOC_DEN_I * LOC_T_I - VOLFRA_PORE_ELE*LOC_DEN_J * LOC_T_J )
                NDOTQ2_TILDE = NDOTQ_TILDE

                ! Make sure we have some sort of velocity (only needed between elements)...

             endif

             ! high order order (low order is an option above)

             WHERE (NDOTQ_TILDE+NDOTQ2_TILDE < 0.0)
                INCOME3=1.0
                ! Make sure the fem description is not biased to the downwind...         
             ELSE WHERE
                INCOME3=0.0
             END WHERE

! limiter *************

! ************NEW LIMITER**************************

             courant_or_minus_one_new(:) = -1.0
             XI_LIMIT(:) = 2.0

             if(NEW_LIMITER) then

                if(IANISOTROPIC==0) then ! this is the only good isotropic limiting for velocity...
! Use ONVDLIMsqrt limiter as its the only limiter that works well for isotropic limiting...

                   CALL ONVDLIM_ANO_MANY_SQRT( NPHASE, &
                      LIMT3(:), FEMTGI_IPHA(:), INCOME3(:), & 
                      LOC_T_I(:), LOC_T_J(:),XI_LIMIT(:),  &
                      TUPWIND_IN(:), TUPWIND_OUT(:) )

                ELSE ! ENDOF if(IANISOTROPIC==0) then

                   CALL ONVDLIM_ANO_MANY( NPHASE, &
                      LIMT3(:), FEMTGI_IPHA(:), INCOME3(:), & 
                      LOC_T_I(:), LOC_T_J(:),XI_LIMIT(:),  &
                      TUPWIND_IN(:), TUPWIND_OUT(:) )

                ENDIF

            else ! original limiter

               do iPHASE=1,nPHASE
                  if(IANISOTROPIC==0) then ! this is the only good isotropic limiting for velocity...
! Use ONVDLIMsqrt limiter as its the only limiter that works well for isotropic limiting...

                     CALL ONVDLIM_ANO_SQRT( cv_nonods, &
                        LIMT3(iPHASE), FEMTGI_IPHA(iPHASE), INCOME3(iPHASE), cv_nodi, cv_nodj, &
                        LOC_T_I(iPHASE), LOC_T_J(iPHASE),   .false., .false., courant_or_minus_one_new(IPHASE), &
                        TUPWIND_IN(iPHASE), TUPWIND_OUT(iPHASE) )

                  ELSE ! ENDOF if(IANISOTROPIC==0) then


                     CALL ONVDLIM_ANO( cv_nonods, &
                        LIMT3(iPHASE), FEMTGI_IPHA(iPHASE), INCOME3(iPHASE), cv_nodi, cv_nodj, &
                        LOC_T_I(iPHASE), LOC_T_J(iPHASE), .false., .false., courant_or_minus_one_new(IPHASE), &
                        TUPWIND_IN(iPHASE), TUPWIND_OUT(iPHASE) )

                  ENDIF ! ENDOF if(IANISOTROPIC==0) then ELSE
               end do

            endif


! limiter *************
             abs_tilde =  0.5*(  ABS_CV_NODI_IPHA  + ( LIMT3   -  LOC_T_I  ) * GRAD_ABS_CV_NODI_IPHA   +   &
                  ABS_CV_NODJ_IPHA  + ( LIMT3   -  LOC_T_J  ) * GRAD_ABS_CV_NODJ_IPHA )

             abs_max=max(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)
             abs_min=min(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)

             abs_tilde    = min(abs_max, max(abs_min,  abs_tilde ))

           if(.false.) then ! new rel perm function
!             NDOTQ_KEEP_IN(1)    =  0.5*( NDOTQ(1)*get_relperm_epsilon(LOC_T_I(1),1, 0.2, 0.3, 2., 1.) + NDOTQ2(1)*get_relperm_epsilon(LOC_T_J(1),1, 0.2, 0.3, 2., 1.) ) &
!               / get_relperm_epsilon(LIMT3(1),1, 0.2, 0.3, 2., 1.)
!             NDOTQ_KEEP_IN(2)    =  0.5*( NDOTQ(2)*get_relperm_epsilon(LOC_T_I(2),2, 0.3, 0.2, 2., 1.) + NDOTQ2(2)*get_relperm_epsilon(LOC_T_J(2),2, 0.3, 0.2, 2., 1.) ) &
!               / get_relperm_epsilon(LIMT3(2),2, 0.3, 0.2, 2., 1.)
                                                                        !VALUES FOR RELPERM ARE HARD CODED!!
             NDOTQ_KEEP_IN(1)    =  0.5*( NDOTQ(1)*inv_get_relperm_epsilon(LOC_T_I(1),1, 0.2, 0.3, 2., 1.) + NDOTQ2(1)*inv_get_relperm_epsilon(LOC_T_J(1),1, 0.2, 0.3, 2., 1.) ) &
               / inv_get_relperm_epsilon(LIMT3(1),1, 0.2, 0.3, 2., 1.)
             NDOTQ_KEEP_IN(2)    =  0.5*( NDOTQ(2)*inv_get_relperm_epsilon(LOC_T_I(2),2, 0.3, 0.2, 2., 1.) + NDOTQ2(2)*inv_get_relperm_epsilon(LOC_T_J(2),2, 0.3, 0.2, 2., 1.) ) &
               / inv_get_relperm_epsilon(LIMT3(2),2, 0.3, 0.2, 2., 1.)
           else

             NDOTQ_KEEP_IN    =  0.5*( NDOTQ*ABS_CV_NODI_IPHA    + NDOTQ2*ABS_CV_NODJ_IPHA )    /abs_tilde
           endif

             INCOME    = MIN(1.0, MAX(0.0,  (NDOTQ_KEEP_IN - NDOTQ)/VTOLFUN( NDOTQ2 - NDOTQ ) ))
             WHERE (abs(NDOTQ2 - NDOTQ).lt. 1.e-7) 
                   INCOME=0.5
             END WHERE



             ! based on switching to the 1st order scheme...
             if(LIMIT_SAT_BASED_UPWIND) then
                WHERE ( min(ABS_CV_NODI_IPHA,ABS_CV_NODJ_IPHA)/ (ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA).lt.1.e-3) ! resort to using the velocity from the highest absorption cell. 
                   !                 if( ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA.gt.1.0e+3 ) then ! resort to using the velocity from the highest absorption cell. 
                   !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) THEN 
                   WHERE (0.5*(NDOTQ+NDOTQ2) < 0.0) 
                      INCOME=1.0
                   ELSE WHERE
                      INCOME=0.0
                   END WHERE
                END WHERE
             endif

          ELSE
             UPWIND_FRAC=0.5
             WHERE (NDOTQ_TILDE < 0.0) 
                !INCOME=0.8
                INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV_J/(MASS_CV_I+MASS_CV_J)
             ELSE WHERE
                !INCOME=0.2
                INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV_I/(MASS_CV_I+MASS_CV_J)
             END WHERE
             !INCOME=0.5
          END IF

          DO IPHASE = 1, NPHASE

          UDGI_ALL(:, IPHASE) = 0.0

          UGI_COEF_ELE_ALL(:, IPHASE, :) = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC =(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_KLOC2=(CV_JLOC-1)*U_NLOC_LEV + U_KLOC_LEV

             UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) &
                                 + SUFEN( U_KLOC,  GI ) * LOC_NU( :, IPHASE, U_KLOC  ) * (1.0-INCOME(IPHASE)) &
                                 + SUFEN( U_KLOC2, GI ) * LOC_NU( :, IPHASE, U_KLOC2 ) * INCOME(IPHASE)


             UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) + 1.0-INCOME(IPHASE)

             UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC2)=UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC2) + INCOME(IPHASE)

          END DO
          END DO ! PHASE LOOP


       ELSE ! Conditional_ELE2: IF( ELE2 /= 0 ) THEN


          UDGI_ALL(:, :) = 0.0

          UGI_COEF_ELE_ALL(:, :, :) = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV

             UDGI_ALL(:, :) = UDGI_ALL(:, :) + SUFEN( U_KLOC,  GI ) * LOC_NU( :, :, U_KLOC  )

             UGI_COEF_ELE_ALL(:, :, U_KLOC)=UGI_COEF_ELE_ALL(:, :, U_KLOC) + 1.0
          END DO


          UDGI2_ALL(:, :) = 0.0
          UGI_COEF_ELE2_ALL(:, :, :) = 0.0
          DO U_KLOC_LEV = 1, U_NLOC_LEV
             U_KLOC=(CV_ILOC-1)*U_NLOC_LEV + U_KLOC_LEV
             U_KLOC2 = U_OTHER_LOC( U_KLOC )
             IF( U_KLOC2 /= 0 ) THEN
                UDGI2_ALL(:, :) = UDGI2_ALL(:, :) + SUFEN( U_KLOC,  GI ) * LOC2_NU( :, :, U_KLOC )

                UGI_COEF_ELE2_ALL(:, :, U_KLOC2)=UGI_COEF_ELE2_ALL(:, :, U_KLOC2) + 1.0

             END IF
          END DO

          NDOTQ  = MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )
          NDOTQ2 = MATMUL( CVNORMX_ALL(:, GI), UDGI2_ALL )

          IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
             DT_I=1.0
             DT_J=1.0
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
             DT_I=MAX(1.E-2,LOC_T_I)
             DT_J=MAX(1.E-2,LOC_T_J)
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
             DT_I=LOC_DEN_I*LOC_T_I
             DT_J=LOC_DEN_J*LOC_T_J
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 4) THEN


             ! DG_ELE_UPWIND==1: the upwind method. 
             ! DG_ELE_UPWIND==3: the best optimal upwind frac.


             !NVEC(1)=CVNORMX(GI)
             !NVEC(2)=CVNORMY(GI)
             !NVEC(3)=CVNORMZ(GI)
             NVEC=CVNORMX_ALL(:, GI)
             ABS_CV_NODI_IPHA      = 0.0
             GRAD_ABS_CV_NODI_IPHA = 0.0
             ABS_CV_NODJ_IPHA      = 0.0
             GRAD_ABS_CV_NODJ_IPHA = 0.0
             DO IPHASE = 1, NPHASE
             DO IDIM=1,NDIM
                V_NODI=0.0
                G_NODI=0.0
                V_NODJ=0.0
                G_NODJ=0.0
                DO JDIM=1,NDIM
                   V_NODI = V_NODI + VI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)
                   G_NODI = G_NODI + GI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)

                   V_NODJ = V_NODJ + VJ_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)
                   G_NODJ = G_NODJ + GJ_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)
                END DO
                ABS_CV_NODI_IPHA(IPHASE)      = ABS_CV_NODI_IPHA(IPHASE) + NVEC(IDIM)*V_NODI
                GRAD_ABS_CV_NODI_IPHA(IPHASE) = GRAD_ABS_CV_NODI_IPHA(IPHASE) + NVEC(IDIM)*G_NODI
                ABS_CV_NODJ_IPHA(IPHASE)      = ABS_CV_NODJ_IPHA(IPHASE) + NVEC(IDIM)*V_NODJ
                GRAD_ABS_CV_NODJ_IPHA(IPHASE) = GRAD_ABS_CV_NODJ_IPHA(IPHASE) + NVEC(IDIM)*G_NODJ
             END DO
             END DO ! PHASE LOOP


             ! Make sure we have some sort of velocity (only needed between elements)...
             ! take the mean of the underlying velocity...
             if(.false.) then
                Q_UNDERLY=( NDOTQ*ABS_CV_NODI_IPHA*MASS_CV_I+ NDOTQ2*ABS_CV_NODJ_IPHA*MASS_CV_J ) &
                     /(MASS_CV_I+MASS_CV_J)
                NDOTQ_TILDE =Q_UNDERLY/ABS_CV_NODI_IPHA
                NDOTQ2_TILDE=Q_UNDERLY/ABS_CV_NODJ_IPHA

             else ! better tested option...

                Q_UNDERLY=0.5*( NDOTQ*ABS_CV_NODI_IPHA + NDOTQ2*ABS_CV_NODJ_IPHA )
                NDOTQ_TILDE =Q_UNDERLY/ABS_CV_NODI_IPHA
                NDOTQ2_TILDE=Q_UNDERLY/ABS_CV_NODJ_IPHA

             endif

             ! These are the new limits of the velocities...
             NDOTQ_KEEP  = NDOTQ_TILDE   ! this is associated with saturation at NODI
             NDOTQ2_KEEP = NDOTQ2_TILDE  ! this is associated with saturation at NODJ

             ! between these limits work out which are associated with the flux limited saturation at this interface. 


             if(ROE_AVE) then

                ! do the Roe average of the rest of the velocity...
                NDOTQ_TILDE  = ( LOC_DEN_I * LOC_T_I * NDOTQ -  &
                     &           LOC_DEN_J * LOC_T_J * NDOTQ2 ) & 
                     / vtolfun( VOLFRA_PORE_ELE*LOC_DEN_I * LOC_T_I - VOLFRA_PORE_ELE2*LOC_DEN_J * LOC_T_J )
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

                         WHERE (0.5*(NDOTQ_TILDE+NDOTQ2_TILDE)>0.0) 
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 * LOC_FEMT(:, CV_KLOC)
                         ELSE WHERE
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 * LOC2_FEMT(:, CV_KLOC)
                         END WHERE

                   END IF
                END DO




                WHERE (0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) 
                   INCOME3=1.0
                ELSE WHERE
                   INCOME3=0.0
                END WHERE

                ! Between elements...




                ! redefine so that it detects oscillations...
                abs_tilde1 =  ABS_CV_NODI_IPHA  + 0.5*( LOC_T_J   -  LOC_T_I  ) * GRAD_ABS_CV_NODI_IPHA   
                abs_tilde2 =  ABS_CV_NODJ_IPHA  + 0.5*( LOC_T_I   -  LOC_T_J  ) * GRAD_ABS_CV_NODJ_IPHA 


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


                WHERE (income3.lt.0.5) ! flux limit
                   wrelax2=wrelax1
                ELSE WHERE
                   wrelax1=wrelax2
                END WHERE


                ! no 4 (more stable than 2 and 3):
                if(between_ele_dg_opt==1) then
                   WHERE ( income3 < 0.5 )  ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*0.5
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*0.5
                   END WHERE

                else if(between_ele_dg_opt==2) then
                   ! no 4 (more stable than 2 and 3) (also make sure we dont violate bounds):
                   WHERE ( income3 < 0.5 )  ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*0.5
                      !              IF(T(CV_NODI_IPHA).LT.0.2) THEN
                      WHERE (ABS_CV_NODI_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*0.5
                      !               IF(T(CV_NODJ_IPHA).LT.0.2) THEN
                      WHERE (ABS_CV_NODJ_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   END WHERE

                   ! no 2(again):
                else if(between_ele_dg_opt==3) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   END WHERE


                   ! no 2(again):
                else if(between_ele_dg_opt==4) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODI_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   ELSE WHERE
                      DT_I = ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODJ_IPHA.GT.1.0E+10) 
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   END WHERE

                   ! no 2(again):
                else if(between_ele_dg_opt==5) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODI_IPHA.GT.1.0E+10) 
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODJ_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   END WHERE

                   ! no 3(again): XXX failed gravity problem...
                else if(between_ele_dg_opt==6) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde1)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde2)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde2)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   END WHERE

                   ! no 2(again,again) XXXXfor the gravity problem:
                else if(between_ele_dg_opt==7) then
                   WHERE ( income3 < 0.5 )  ! upwind
                      DT_I = (1.-wrelax1)*(1.-min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1))  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde2)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*(1.-min(1.0,ABS_CV_NODI_IPHA/abs_tilde2))  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   END WHERE

                   ! no 3 soln...
                else if(between_ele_dg_opt==8) then
                   DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA  ) &
                        /(ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA)   + wrelax1*0.5
                   DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA  ) &
                        /(ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA)   + wrelax2*0.5

                   ! no 2 soln...
                   !        DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA*MASS_CV_I  ) &
                   !             /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)   + wrelax1*0.5
                   !        DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA*MASS_CV_J  ) &
                   !             /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)   + wrelax2*0.5

                   !        DTOLD_I = ( (1.-wrelaxold1)*ABS_CV_NODI_IPHA*MASS_CV_I  ) &
                   !                /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)  + wrelaxold1*0.5
                   !        DTOLD_J = ( (1.-wrelaxold2)*ABS_CV_NODJ_IPHA*MASS_CV_J  ) &
                   !                /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)  + wrelaxold2*0.5

                else if(between_ele_dg_opt==9) then
                   ! no 0 soln (more stable than 2 & 3)...
                   DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA*MASS_CV_I + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J ) &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA*MASS_CV_J + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I ) &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ! no 0 soln (more stable than 2 & 3)...
                else if(between_ele_dg_opt==10) then
                   DT_I = ABS_CV_NODJ_IPHA*MASS_CV_J  &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   DT_J = ABS_CV_NODI_IPHA*MASS_CV_I  &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                endif





                ! END OF IF ( DG_ELE_UPWIND==3 ) THEN ...
             ELSE  ! Low order...
                !if ( abs(T(CV_NODI_IPHA)-T(CV_NODJ_IPHA)).lt.1.e-4 ) then ! low order


                INCOME =0.5*ABS_CV_NODI_IPHA* MASS_CV_I /(0.5*(ABS_CV_NODI_IPHA*MASS_CV_I +ABS_CV_NODJ_IPHA*MASS_CV_J ))

                DO IPHASE = 1, NPHASE
                ! limit the volume fraction based on net letting flux go into certain elements...
                !     if(.false.) then
                if(LIMIT_SAT_BASED_INTERP) then
                   enforce_abs = .false.

                   if(iphase==1) then
                      !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) >= 0.0) THEN 
                      IF(0.5*(NDOTQ(IPHASE)+NDOTQ2(IPHASE)) >= 0.0) THEN 
                         IF(LOC_T_I(IPHASE) < 0.5) then
                            enforce_abs = .true.
                            relax=min( 10.*(0.5-LOC_T_I(IPHASE)), 1.0)
                         ENDIF
                      ELSE
                         IF(LOC_T_J(IPHASE) < 0.5) THEN
                            enforce_abs = .true. 
                            relax=min( 10.*(0.5-LOC_T_J(IPHASE)), 1.0)
                         ENDIF
                      ENDIF
                   endif

                   if(iphase==2) then
                      !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) >= 0.0) THEN 
                      IF(0.5*(NDOTQ(IPHASE)+NDOTQ2(IPHASE)) >= 0.0) THEN 
                         IF(LOC_T_I(IPHASE) < 0.5) then
                            enforce_abs = .true.
                            relax=min( 10.*(0.5-LOC_T_I(IPHASE)), 1.0)
                         ENDIF
                      ELSE
                         IF(LOC_T_J(IPHASE) < 0.5) THEN
                            enforce_abs = .true. 
                            relax=min( 10.*(0.5-LOC_T_J(IPHASE)), 1.0)
                         ENDIF
                      ENDIF
                   endif



                   if(enforce_abs) then
                      INCOME(IPHASE) =  (1.-RELAX)*INCOME(IPHASE)  + RELAX*( (1.-INCOME(IPHASE))**2 *(1. - 3.*INCOME(IPHASE)) +  INCOME(IPHASE)*(1.-INCOME(IPHASE))*(1.+3.*INCOME(IPHASE))  ) ! non-lumped
                   endif

                   income(IPHASE)=min(1.0, max(0.0, income(IPHASE)))
                endif

                END DO !PHASE LOOP



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
                DT_I=VOLFRA_PORE_ELE *DT_I
                DT_J=VOLFRA_PORE_ELE2*DT_J
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



          DO IPHASE = 1, NPHASE
             UDGI_INT_ALL(:, IPHASE) = DT_I(IPHASE) * UDGI_ALL(:, IPHASE) + DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)
          END DO

          IF( CV_DG_VEL_INT_OPT < 0 ) THEN
   !          FLAbort('2821') ! should not be going in here??
   !          if(got_dt_ij) FLAbort('2611') ! we can not use these CV_DG_VEL_INT_OPT < 0, got_dt_ij=.true. options together.

             !NDOTQ_INT = CVNORMX( GI ) * UDGI_INT + CVNORMY( GI ) * VDGI_INT + &
             !     CVNORMZ( GI ) * WDGI_INT
             NDOTQ_INT = MATMUL( CVNORMX_ALL(:, GI), UDGI_INT_ALL)

             WHERE ( NDOTQ_INT <= 0.0 )  ! Incoming
                !DT_I=1.0
                DT_J=DT_I+DT_J
             ELSE WHERE
                DT_I=DT_I+DT_J
                !DT_J=1.0
             END WHERE

             DO IPHASE = 1, NPHASE
                UDGI_INT_ALL(:, IPHASE) = (DT_I(IPHASE) * UDGI_ALL(:, IPHASE) + DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)) / (DT_I(IPHASE) + DT_J(IPHASE))
             END DO

          END IF

          UDGI_ALL = UDGI_INT_ALL

          DO IPHASE = 1, NPHASE
             UGI_COEF_ELE_ALL(:, IPHASE, :) = DT_I(IPHASE) * UGI_COEF_ELE_ALL(:, IPHASE, :) 

             UGI_COEF_ELE2_ALL(:, IPHASE, :) = DT_J(IPHASE) * UGI_COEF_ELE2_ALL(:, IPHASE, :) 
          END DO

       END IF Conditional_ELE2

    END IF Conditional_SELE

    NDOTQ =  MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )

    ! Define whether flux is incoming or outgoing, depending on direction of flow
    INCOME = 1.
    WHERE ( NDOTQ >= 0. ) 
          INCOME = 0.
    END WHERE
    !  IF( NDOTQ_tilde >= 0. ) INCOME = 0.
    !IF( NDOTQ > 0. ) INCOME = 0.

    RETURN  

CONTAINS


        SUBROUTINE ONVDLIM_ANO_MANY( NFIELD, &
       TDLIM, TDCEN, INCOME, &
       ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT,  &
       TUPWIN, TUPWI2 )
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
    INTEGER, intent( in ) :: NFIELD
    REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM  
    REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2 
    REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
    ! Local variables   
    REAL :: UCIN(NFIELD), UCOU(NFIELD), TDELE(NFIELD), DENOIN(NFIELD), CTILIN(NFIELD), DENOOU(NFIELD), &
         CTILOU(NFIELD), FTILIN(NFIELD), FTILOU(NFIELD)


       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE 

       DENOIN = TOLFUN_MANY( TDELE - TUPWIN )

       UCIN = ETDNEW_PELEOT 
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT

       DENOOU = TOLFUN_MANY( TDELE - TUPWI2 )
       UCOU = ETDNEW_PELE 
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU



       FTILIN = ( TDCEN - TUPWIN ) / DENOIN
       FTILOU = ( TDCEN - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 )



    RETURN

  END SUBROUTINE ONVDLIM_ANO_MANY



        SUBROUTINE ONVDLIM_ANO_MANY_SQRT( NFIELD, &
       TDLIM, TDCEN, INCOME, &
       ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT,  &
       TUPWIN, TUPWI2 )
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
    INTEGER, intent( in ) :: NFIELD
    REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM  
    REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2 
    REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
    ! Local variables   
    INTEGER, PARAMETER :: POWER=2
    REAL :: UCIN(NFIELD), UCOU(NFIELD), TDELE(NFIELD), DENOIN(NFIELD), CTILIN(NFIELD), DENOOU(NFIELD), &
         CTILOU(NFIELD), FTILIN(NFIELD), FTILOU(NFIELD)


       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE 

       DENOIN = TOLFUN_MANY( TDELE** POWER - TUPWIN** POWER )

       UCIN = ETDNEW_PELEOT 
       CTILIN = ( UCIN ** POWER- TUPWIN** POWER ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT

       DENOOU = TOLFUN_MANY( TDELE** POWER - TUPWI2** POWER )
       UCOU = ETDNEW_PELE 
       CTILOU = ( UCOU** POWER - TUPWI2** POWER ) / DENOOU



       FTILIN = ( TDCEN** POWER - TUPWIN** POWER ) / DENOIN
       FTILOU = ( TDCEN** POWER - TUPWI2 ** POWER) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN** POWER + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2** POWER + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)



    RETURN

  END SUBROUTINE ONVDLIM_ANO_MANY_SQRT




        SUBROUTINE ONVDLIM_ANO( TOTELE, &
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

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD

!     if((PELE==433).and.(PELEOT==435)) then
!          print *,'***INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT:',INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT
!          print *,'***FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM:',FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM
!          print *,'***NOLIMI, FIRORD, TUPWIN2, TUPWI22:', NOLIMI, FIRORD, TUPWIN2, TUPWI22
!     endif

    RETURN

  END SUBROUTINE ONVDLIM_ANO


        SUBROUTINE ONVDLIM_ANO_SQRT( TOTELE, &
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
    INTEGER, PARAMETER :: POWER=2
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
       UCIN = ETDNEW_PELE** POWER
       DENOIN = 1.
       CTILIN = 0.

       ! Calculate normalisation parameters for out going velocities 
       TUPWI2 = ETDNEW_PELE** POWER
       UCOU = ETDNEW_PELE ** POWER
       DENOOU = 1.
       CTILOU = 0.

    END IF Conditional_PELEOT

    Conditional_FIRORD: IF( FIRORD ) THEN ! Velocity is pointing into element

       ! Velocity is going out of element
       TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 
       TDLIM = MAX(TDLIM,0.0) ** (1.0/POWER)

    ELSE

       FTILIN = ( MAX(0.0,TDCEN) ** POWER - TUPWIN ) / DENOIN
       FTILOU = ( MAX(0.0,TDCEN) ** POWER - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, COURANT_OR_MINUS_ONE ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, COURANT_OR_MINUS_ONE ) * DENOOU )
!       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, -1. ) * DENOIN ) &
!            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, -1. ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD
!       if(( PELE==433).and.( PELEOT==435))  print *,'++ firord,TDLIM, INCOME, TUPWI2, FTILOU, CTILOU, COURANT_OR_MINUS_ONE:', &
!                                                        firord,TDLIM, INCOME, TUPWI2, FTILOU, CTILOU, COURANT_OR_MINUS_ONE
!     if((PELE==433).and.(PELEOT==435)) then
!          print *,'---new INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT:',INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT
!          print *,'---new FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM:',FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM
!          print *,'---new NOLIMI, FIRORD, TUPWIN2, TUPWI22:', NOLIMI, FIRORD, TUPWIN2, TUPWI22
!     endif
    RETURN

  END SUBROUTINE ONVDLIM_ANO_SQRT

  END SUBROUTINE GET_INT_VEL_OVERLAP_NEW








      SUBROUTINE GET_INT_VEL_POROUS_VEL( NPHASE, NDOTQ,INCOME, &
       HDC, GI, SUFEN, U_NLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
       LOC_T_I, LOC_T_J, LOC_FEMT, LOC2_FEMT, LOC_DEN_I, LOC_DEN_J, &
       LOC_NU, LOC2_NU, SLOC_NU, &
       CV_NODI, CV_NODJ, CVNORMX_ALL,  &
       CV_DG_VEL_INT_OPT, ELE, ELE2, U_OTHER_LOC, &
       SELE, U_SNLOC, STOTEL, U_SLOC2LOC, SUF_U_BC_ALL, WIC_U_BC_ALL, &
       SUF_SIG_DIAGTEN_BC, &
       UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL, &
       VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2, CV_ELE_TYPE, CV_SLOC2LOC, CV_NLOC, CV_SNLOC, CV_ILOC, CV_JLOC, SCVFEN, CV_OTHER_LOC, &
       VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS,  VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS, &
       INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS, &
       UDGI_ALL, &
       MASS_CV_I, MASS_CV_J, NDIM, MAT_NLOC, MAT_NONODS, &
       IN_ELE_UPWIND, DG_ELE_UPWIND, &
       IANISOTROPIC,  &
       TUPWIND_IN, TUPWIND_OUT)
    !================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===
    ! Calculate NDOTQ and INCOME on the CV boundary at quadrature pt GI. 
    ! it assumes an overlapping decomposition approach for velocity. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, GI, U_NLOC, CV_SNLOC, SCVNGI, TOTELE, U_NONODS, CV_NONODS, &
         CV_NODJ, CV_NODI, CV_DG_VEL_INT_OPT, ELE, ELE2, &
         SELE, U_SNLOC, STOTEL, CV_ELE_TYPE, CV_NLOC, CV_ILOC, CV_JLOC, &
         NDIM, MAT_NLOC, MAT_NONODS,  &
         IN_ELE_UPWIND, DG_ELE_UPWIND
    REAL, intent( in ) :: HDC, MASS_CV_I, MASS_CV_J, VOLFRA_PORE_ELE, VOLFRA_PORE_ELE2
    REAL, DIMENSION( : ), intent( inout ) :: NDOTQ, INCOME
    REAL, DIMENSION( :, :  ), intent( inout ) :: UDGI_ALL
    INTEGER, DIMENSION( : ), intent( in ) :: U_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: U_SLOC2LOC
    INTEGER, DIMENSION( :, :, : ), intent( in ) :: WIC_U_BC_ALL
    REAL, DIMENSION( :, :  ), intent( in ) :: SUFEN
    REAL, DIMENSION( :, : ), intent( in ) :: LOC_FEMT, LOC2_FEMT
    REAL, DIMENSION( : ), intent( in ) :: LOC_T_I, LOC_T_J, LOC_DEN_I, LOC_DEN_J
    REAL, DIMENSION( :, :, : ), intent( in ) ::  LOC_NU, LOC2_NU, SLOC_NU
    REAL, DIMENSION( :, : ), intent( in ) :: CVNORMX_ALL
    REAL, DIMENSION( :, :, : ), intent( in ) :: SUF_U_BC_ALL
    REAL, DIMENSION( : , : ), intent( in ) :: SUF_SIG_DIAGTEN_BC
    REAL, DIMENSION( :, :, : ), intent( inout ) :: UGI_COEF_ELE_ALL, UGI_COEF_ELE2_ALL
    REAL, DIMENSION( : , :  ), intent( in ) :: SCVFEN
    INTEGER, DIMENSION( : ), intent( in ) :: CV_OTHER_LOC
    INTEGER, DIMENSION( : ), intent( in ) :: CV_SLOC2LOC
    REAL, DIMENSION( :, :, : ), intent( in ) :: VI_LOC_OPT_VEL_UPWIND_COEFS, GI_LOC_OPT_VEL_UPWIND_COEFS, &
                                                VJ_LOC_OPT_VEL_UPWIND_COEFS, GJ_LOC_OPT_VEL_UPWIND_COEFS
    REAL, DIMENSION( :, :, : ), intent( in ) :: INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS

      INTEGER, intent( in ) :: IANISOTROPIC
      REAL, DIMENSION( NPHASE ), intent( in ) :: TUPWIND_IN, TUPWIND_OUT

    ! Local variables
    REAL :: UDGI,VDGI,WDGI,  &
         UDGI2,VDGI2,WDGI2,  &
         UDGI_INT,VDGI_INT,WDGI_INT, &
         FEMTOLDGI_IPHA, OVER_RELAX, V_NODI, G_NODI, V_NODJ, G_NODJ, & 
         GEOMTOLDGI_IPHA, W_UPWIND, W_UPWINDOLD, INCOME4, &
         UPWIND_FRAC, TUPWIN, TUPWI2
    REAL, DIMENSION(NPHASE) :: NDOTQ_INT,NDOTQ2, FEMTGI_IPHA, GEOMTGI_IPHA, &
                               ABS_CV_NODI_IPHA, ABS_CV_NODJ_IPHA, &
                               GRAD_ABS_CV_NODI_IPHA, GRAD_ABS_CV_NODJ_IPHA, &
                               DT_I,DT_J, INCOME3
    REAL :: LIMT3(NPHASE) 
    REAL :: NVEC(NDIM),SUF_SIG_DIAGTEN_BC_GI(NDIM), UGI_TMP(NDIM)
    INTEGER :: U_KLOC,U_NODK,U_NODK2,U_NODK2_IPHA,U_NODK_IPHA,U_KLOC2,U_SKLOC, &
         U_SNODK,U_SNODK_IPHA, II,  COUNT, &
         CV_KLOC, CV_KNOD, &
         CV_KLOC2, CV_KNOD2, IDIM, JDIM, IJ, MAT_NODI, MAT_NODJ, &
         CV_SKLOC, CV_SNODK, CV_SNODK_IPHA, CV_STAR_IPHA, CV_KNOD_IPHA, CV_KNOD2_IPHA, &
         U_NODK3,U_NODK3_IPHA, U_KLOC3
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
    LOGICAL, PARAMETER :: FORCE_UPWIND_VEL_DG_ELE = .false.
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
! The new simple limiter NEW_LIMITER can produce very slightly different results 
!    LOGICAL, PARAMETER :: NEW_LIMITER = .true.
    LOGICAL, PARAMETER :: NEW_LIMITER = .false.
    REAL :: TMIN_STORE, TMAX_STORE
    REAL, DIMENSION(NPHASE) :: PERM_TILDE,NDOTQ_TILDE, NDOTQ2_TILDE, NDOTQOLD_TILDE, NDOTQOLD2_TILDE, rden_ave, Q_UNDERLY
    REAL, DIMENSION(NPHASE) :: NDOTQ_KEEP_IN,  NDOTQ_KEEP, NDOTQ2_KEEP
    ! coefficients for this element ELE
    real :: gamma,  grad2nd
    real :: abs_tilde_2nd, abs_tildeold_2nd, max_nodtq_keep, min_nodtq_keep
    real :: w_weight, relax

    real, DIMENSION(NPHASE) :: DT_I_upwind, DT_J_upwind
    real :: dt_max, dt_min
    real, DIMENSION(NPHASE) :: abs_tilde1, abs_tilde2, abs_tilde, abs_max, abs_min, w_relax
    real, DIMENSION(NPHASE) :: wrelax, wrelax1, wrelax2
    real :: T_PELE, T_PELEOT, TMIN_PELE, TMAX_PELE, TMIN_PELEOT, TMAX_PELEOT

    ! Local variable for indirect addressing
    REAL, DIMENSION ( NDIM, NPHASE ) :: UDGI2_ALL, UDGI_INT_ALL, UDGI_ALL_FOR_INV, ROW_SUM_INV_VI, ROW_SUM_INV_VJ
!    REAL, DIMENSION ( NDIM, NPHASE ) :: UDGI_ALL, UDGI2_ALL, UDGI_INT_ALL, UDGI_ALL_FOR_INV, ROW_SUM_INV_VI, ROW_SUM_INV_VJ
!    REAL, DIMENSION ( NDIM, NDIM, NPHASE ) :: INV_VI_LOC_OPT_VEL_UPWIND_COEFS, INV_VJ_LOC_OPT_VEL_UPWIND_COEFS
    real :: courant_or_minus_one_new(Nphase),XI_LIMIT(Nphase), VEC_NDIM(NDIM), VEC2_NDIM(NDIm), UDGI_ALL_OTHER(NDIM, NPHASE)
    INTEGER :: IPHASE, CV_NODI_IPHA, CV_NODJ_IPHA


    ! coefficients for this element ELE
    UGI_COEF_ELE_ALL=0.0

    ! coefficients for this element ELE2
    UGI_COEF_ELE2_ALL=0.0


     DO IPHASE=1,NPHASE
        DO IDIM=1,NDIM
           ROW_SUM_INV_VI(IDIM,IPHASE)=SUM(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,:,IPHASE))
           ROW_SUM_INV_VJ(IDIM,IPHASE)=SUM(INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(IDIM,:,IPHASE))
        END DO
     END DO


    got_dt_ij=.false.


    Conditional_SELE: IF( SELE /= 0 ) THEN ! On the boundary of the domain. 
     DO IPHASE = 1, NPHASE
       IF( WIC_U_BC_ALL( 1, IPHASE, SELE) /= WIC_U_BC_DIRICHLET ) THEN ! velocity free boundary
          UDGI_ALL(:, IPHASE) = 0.0
          DO U_KLOC = 1, U_NLOC
             UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * LOC_NU( :, IPHASE, U_KLOC )
          END DO
          UDGI_ALL(:, IPHASE) = matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI_ALL(:, IPHASE))

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
!            print *,'sele, SUF_SIG_DIAGTEN_BC_GI:',sele, SUF_SIG_DIAGTEN_BC_GI

          ! Only modify boundary velocity for incoming velocity...
          DO U_KLOC = 1, U_NLOC
             IF (DOT_PRODUCT(UDGI_ALL(:, IPHASE), CVNORMX_ALL(:, GI)).LT.0.0) THEN ! Incomming...
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=SUF_SIG_DIAGTEN_BC_GI(:)
             ELSE
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=1.0
             ENDIF

             UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)= matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC))
          END DO

          IF(DOT_PRODUCT(UDGI_ALL(:, IPHASE), CVNORMX_ALL(:, GI)).LT.0.0) THEN ! Incomming...
             UGI_TMP = SUF_SIG_DIAGTEN_BC_GI(:) * UDGI_ALL(:, IPHASE)
             UDGI_ALL(:, IPHASE) = UGI_TMP
          ENDIF

       ELSE ! Specified vel bc.

          UDGI_ALL(:, IPHASE) = 0.0
          UDGI_ALL_FOR_INV(:, IPHASE) = 0.0
          UGI_COEF_ELE_ALL(:, IPHASE, :) = 0.0 
          DO U_SKLOC = 1, U_SNLOC
             U_KLOC = U_SLOC2LOC( U_SKLOC )

             IF(WIC_U_BC_ALL( 1, IPHASE, SELE ) == 10) THEN
                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI ) * 0.5 * &
                                      SUF_U_BC_ALL(:, IPHASE, U_SNLOC* (SELE-1) +U_SKLOC)

                UDGI_ALL_FOR_INV(:, IPHASE) = UDGI_ALL_FOR_INV(:, IPHASE) + SUFEN( U_KLOC, GI ) * 0.5 * &
                                      SLOC_NU(:, IPHASE, U_SKLOC)

!                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=0.5
!                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)= matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC))
                UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=0.5*ROW_SUM_INV_VI(:,IPHASE) 
             ELSE

                UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) + SUFEN( U_KLOC, GI )*SUF_U_BC_ALL(:, IPHASE, U_SNLOC* (SELE-1) +U_SKLOC)
             END IF
          END DO
          UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE)  + matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI_ALL_FOR_INV(:, IPHASE))


       END IF
     END DO ! PHASE LOOP

    ELSE ! Conditional_SELE. Not on the boundary of the domain.
       Conditional_ELE2: IF(( ELE2 == 0 ).OR.( ELE2 == ELE)) THEN!same element

          UDGI_ALL(:, :) = 0.0
          UDGI2_ALL(:, :) = 0.0

          DO U_KLOC = 1, U_NLOC
             UDGI_ALL(:, :) = UDGI_ALL(:, :) + SUFEN( U_KLOC, GI ) * LOC_NU( :, :, U_KLOC )
             UDGI2_ALL(:, :) = UDGI2_ALL(:, :) + SUFEN( U_KLOC, GI ) * LOC_NU( :, :, U_KLOC )
          END DO

          DO IPHASE=1,NPHASE
             UDGI_ALL(:, IPHASE)  = matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI_ALL(:, IPHASE))
             UDGI2_ALL(:, IPHASE) = matmul(INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI2_ALL(:, IPHASE))
          END DO 

          NDOTQ  = MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )
          NDOTQ2 = MATMUL( CVNORMX_ALL(:, GI), UDGI2_ALL )

          NDOTQ_TILDE = 0.5 * ( NDOTQ +  NDOTQ2 )
          NDOTQ2_TILDE =   NDOTQ_TILDE

          IF((IN_ELE_UPWIND==1).OR.FORCE_UPWIND_VEL) THEN

             WHERE (NDOTQ_TILDE < 0.0) 
                INCOME=1.0
             ELSE WHERE
                INCOME=0.0
             END WHERE

             WHERE (abs(NDOTQ2 - NDOTQ).lt. 1.e-7) 
                   INCOME=0.5
             END WHERE

             DO IPHASE=1,NPHASE
                UDGI_ALL(:, IPHASE)  =UDGI_ALL(:, IPHASE) * (1.-INCOME(IPHASE)) + UDGI2_ALL(:, IPHASE) *INCOME(IPHASE)
             END DO




          ELSE IF(IN_ELE_UPWIND==2) THEN ! the best

             UPWIND_FRAC=0.8
             WHERE (NDOTQ_TILDE < 0.0)
                INCOME=0.8
             ELSE WHERE
                INCOME=0.2 ! MAKE SURE IT IS RIGHT?
             END WHERE
             DO IPHASE=1,NPHASE
                UDGI_ALL(:, IPHASE)  =UDGI_ALL(:, IPHASE) * (1.-INCOME(IPHASE)) + UDGI2_ALL(:, IPHASE) *INCOME(IPHASE)
             END DO

          ELSE IF(IN_ELE_UPWIND==3) THEN ! the best optimal upwind frac.
             !Calculate vel at GI
             FEMTGI_IPHA = matmul(LOC_FEMT, SCVFEN(:,GI) )

             ! Central is fine as its within an element with equally spaced nodes. 
             !FEMTGI_IPHA = 0.5*( t(cv_nodi_ipha)+t(cv_nodj_ipha) )
             !FEMTOLDGI_IPHA = 0.5*( told(cv_nodi_ipha)+told(cv_nodj_ipha) )
             if(IANISOTROPIC==0) then ! this is the only needed for isotropic limiting for velocity...
                FEMTGI_IPHA = ( MASS_CV_J * LOC_T_I + &
                     MASS_CV_I * LOC_T_J ) / (MASS_CV_I+MASS_CV_J)
             endif
             GEOMTGI_IPHA = ( MASS_CV_J * LOC_T_I + &
                  MASS_CV_I * LOC_T_J ) / (MASS_CV_I+MASS_CV_J)

             !We perform: n' * sigma * n
             DO IPHASE = 1, NPHASE
                 ABS_CV_NODI_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
                 GRAD_ABS_CV_NODI_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(GI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
                 ABS_CV_NODJ_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
                 GRAD_ABS_CV_NODJ_IPHA(IPHASE) = dot_product(CVNORMX_ALL(:, GI),matmul(GJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE), CVNORMX_ALL(:, GI)))
             END DO


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
                NDOTQ_TILDE  = ( LOC_DEN_I * LOC_T_I * NDOTQ -  &
                     &           LOC_DEN_J * LOC_T_J * NDOTQ2 ) & 
                     / vtolfun( VOLFRA_PORE_ELE*LOC_DEN_I * LOC_T_I - VOLFRA_PORE_ELE*LOC_DEN_J * LOC_T_J )
                NDOTQ2_TILDE = NDOTQ_TILDE

                ! Make sure we have some sort of velocity (only needed between elements)...

             endif



             ! high order order (low order is an option above)



             WHERE (0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) 
                INCOME3=1.0
                ! Make sure the fem description is not biased to the downwind...         
             ELSE WHERE
                INCOME3=0.0
             END WHERE



! limiter *************

! ************NEW LIMITER**************************

             courant_or_minus_one_new(:) = -1.0
             XI_LIMIT(:) = 2.0

             if(NEW_LIMITER) then

                if(IANISOTROPIC==0) then ! this is the only good isotropic limiting for velocity...
! Use ONVDLIMsqrt limiter as its the only limiter that works well for isotropic limiting...

                   CALL ONVDLIM_ANO_MANY_SQRT( NPHASE, &
                      LIMT3(:), FEMTGI_IPHA(:), INCOME3(:), & 
                      LOC_T_I(:), LOC_T_J(:),XI_LIMIT(:),  &
                      TUPWIND_IN(:), TUPWIND_OUT(:) )

                ELSE ! ENDOF if(IANISOTROPIC==0) then

                   CALL ONVDLIM_ANO_MANY( NPHASE, &
                      LIMT3(:), FEMTGI_IPHA(:), INCOME3(:), & 
                      LOC_T_I(:), LOC_T_J(:),XI_LIMIT(:),  &
                      TUPWIND_IN(:), TUPWIND_OUT(:) )

                ENDIF

            else ! original limiter

               do iPHASE=1,nPHASE
                  if(IANISOTROPIC==0) then ! this is the only good isotropic limiting for velocity...
! Use ONVDLIMsqrt limiter as its the only limiter that works well for isotropic limiting...

                     CALL ONVDLIM_ANO_SQRT( cv_nonods, &
                        LIMT3(iPHASE), FEMTGI_IPHA(iPHASE), INCOME3(iPHASE), cv_nodi, cv_nodj, &
                        LOC_T_I(iPHASE), LOC_T_J(iPHASE),   .false., .false., courant_or_minus_one_new(IPHASE), &
                        TUPWIND_IN(iPHASE), TUPWIND_OUT(iPHASE) )

                  ELSE ! ENDOF if(IANISOTROPIC==0) then

!                     CALL ONVDLIM_ANO_sqrt( cv_nonods, &
                     CALL ONVDLIM_ANO( cv_nonods, &
                        LIMT3(iPHASE), FEMTGI_IPHA(iPHASE), INCOME3(iPHASE), cv_nodi, cv_nodj, &
                        LOC_T_I(iPHASE), LOC_T_J(iPHASE),   FORCE_UPWIND_VEL, .false., courant_or_minus_one_new(IPHASE), &
                        TUPWIND_IN(iPHASE), TUPWIND_OUT(iPHASE) )

                  ENDIF ! ENDOF if(IANISOTROPIC==0) then ELSE
               end do

            endif


! limiter *************



          DO IPHASE=1,NPHASE
             abs_tilde(IPHASE) =  0.5*(  ABS_CV_NODI_IPHA(IPHASE)  + ( LIMT3(IPHASE)   -  LOC_T_I(IPHASE)  ) * GRAD_ABS_CV_NODI_IPHA(IPHASE)   +   &
                  ABS_CV_NODJ_IPHA(IPHASE)  + ( LIMT3(IPHASE)   -  LOC_T_J(IPHASE)  ) * GRAD_ABS_CV_NODJ_IPHA(IPHASE) )
          END DO

             abs_max=max(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)
             abs_min=min(ABS_CV_NODI_IPHA,  ABS_CV_NODJ_IPHA)

             abs_tilde    = min(abs_max, max(abs_min,  abs_tilde ))

           if(.false.) then ! new rel perm function
!             NDOTQ_KEEP_IN(1)    =  0.5*( NDOTQ(1)*get_relperm_epsilon(LOC_T_I(1),1, 0.2, 0.3, 2., 1.) + NDOTQ2(1)*get_relperm_epsilon(LOC_T_J(1),1, 0.2, 0.3, 2., 1.) ) &
!               / get_relperm_epsilon(LIMT3(1),1, 0.2, 0.3, 2., 1.)
!             NDOTQ_KEEP_IN(2)    =  0.5*( NDOTQ(2)*get_relperm_epsilon(LOC_T_I(2),2, 0.3, 0.2, 2., 1.) + NDOTQ2(2)*get_relperm_epsilon(LOC_T_J(2),2, 0.3, 0.2, 2., 1.) ) &
!               / get_relperm_epsilon(LIMT3(2),2, 0.3, 0.2, 2., 1.)
                                                                        !VALUES FOR RELPERM ARE HARD CODED!!
             NDOTQ_KEEP_IN(1)    =  0.5*( NDOTQ(1)*inv_get_relperm_epsilon(LOC_T_I(1),1, 0.2, 0.3, 2., 1.) + NDOTQ2(1)*inv_get_relperm_epsilon(LOC_T_J(1),1, 0.2, 0.3, 2., 1.) ) &
               / inv_get_relperm_epsilon(LIMT3(1),1, 0.2, 0.3, 2., 1.)
             NDOTQ_KEEP_IN(2)    =  0.5*( NDOTQ(2)*inv_get_relperm_epsilon(LOC_T_I(2),2, 0.3, 0.2, 2., 1.) + NDOTQ2(2)*inv_get_relperm_epsilon(LOC_T_J(2),2, 0.3, 0.2, 2., 1.) ) &
               / inv_get_relperm_epsilon(LIMT3(2),2, 0.3, 0.2, 2., 1.)
           else

             NDOTQ_KEEP_IN    =  0.5*( NDOTQ*ABS_CV_NODI_IPHA    + NDOTQ2*ABS_CV_NODJ_IPHA )    /abs_tilde
           endif


             INCOME    = MIN(1.0, MAX(0.0,  (NDOTQ_KEEP_IN - NDOTQ)/VTOLFUN( NDOTQ2 - NDOTQ ) ))
             WHERE (abs(NDOTQ2 - NDOTQ).lt. 1.e-7) 
                   INCOME=0.5
             END WHERE



             ! based on switching to the 1st order scheme...
             if(LIMIT_SAT_BASED_UPWIND) then
                WHERE ( min(ABS_CV_NODI_IPHA,ABS_CV_NODJ_IPHA)/ (ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA).lt.1.e-3) ! resort to using the velocity from the highest absorption cell. 
                   !                 if( ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA.gt.1.0e+3 ) then ! resort to using the velocity from the highest absorption cell. 
                   !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) THEN 
                   WHERE (0.5*(NDOTQ+NDOTQ2) < 0.0) 
                      INCOME=1.0
                   ELSE WHERE
                      INCOME=0.0
                   END WHERE
                END WHERE
             endif
          ELSE

             UPWIND_FRAC=0.5
             WHERE (NDOTQ_TILDE < 0.0) 
                !INCOME=0.8
                INCOME= 1.0 - 2.*(1.-UPWIND_FRAC)*MASS_CV_J/(MASS_CV_I+MASS_CV_J)
             ELSE WHERE
                !INCOME=0.2
                INCOME= 2.*(1.-UPWIND_FRAC)*MASS_CV_I/(MASS_CV_I+MASS_CV_J)
             END WHERE
             !INCOME=0.5
          END IF

          DO IPHASE = 1, NPHASE

          UDGI_ALL(:, IPHASE) = 0.0
          UDGI_ALL_OTHER(:, IPHASE) = 0.0

          UGI_COEF_ELE_ALL(:, IPHASE, :) = 0.0
          DO U_KLOC = 1, U_NLOC

             UDGI_ALL(:, IPHASE) = UDGI_ALL(:, IPHASE) &
                                 + SUFEN( U_KLOC,  GI ) * LOC_NU( :, IPHASE, U_KLOC  ) * (1.0-INCOME(IPHASE))
             UDGI_ALL_OTHER(:, IPHASE) = UDGI_ALL_OTHER(:, IPHASE) &
                                 + SUFEN( U_KLOC, GI ) * LOC_NU( :, IPHASE, U_KLOC ) * INCOME(IPHASE)

!             VEC_NDIM(:) =1.0-INCOME(IPHASE)
!             VEC2_NDIM(:)=INCOME(IPHASE)

!             UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC) + matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),VEC_NDIM(:)) &
!                                                                                     + matmul(INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),VEC2_NDIM(:))

             UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)=ROW_SUM_INV_VI(:,IPHASE)* (1.0-INCOME(IPHASE)) &
                                                +ROW_SUM_INV_VJ(:,IPHASE)* INCOME(IPHASE)

          END DO

          UDGI_ALL(:, IPHASE) = matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI_ALL(:, IPHASE)) &
                              + matmul(INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI_ALL_OTHER(:, IPHASE))
          END DO ! PHASE LOOP


       ELSE ! Conditional_ELE2: IF( ELE2 /= 0 ) THEN

          UDGI_ALL(:, :) = 0.0
          UDGI2_ALL(:, :) = 0.0

          UGI_COEF_ELE_ALL(:, :, :) = 0.0
          UGI_COEF_ELE2_ALL(:, :, :) = 0.0

!          print *,'U_OTHER_LOC:', U_OTHER_LOC
!          stop 811
          DO U_KLOC = 1, U_NLOC
             U_KLOC2 = U_OTHER_LOC( U_KLOC )
             IF( U_KLOC2 /= 0 ) THEN
                UDGI_ALL(:, :)  = UDGI_ALL(:, :)  + SUFEN( U_KLOC,  GI ) * LOC_NU( :, :, U_KLOC  )
                UDGI2_ALL(:, :) = UDGI2_ALL(:, :) + SUFEN( U_KLOC,  GI ) * LOC2_NU( :, :, U_KLOC )

!                VEC2_NDIM(:) =1.0
                DO IPHASE=1,NPHASE
!                   UGI_COEF_ELE2_ALL(:, IPHASE, U_KLOC2)=UGI_COEF_ELE2_ALL(:, IPHASE, U_KLOC2) + matmul(INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),VEC2_NDIM(:))
                   UGI_COEF_ELE_ALL(:, IPHASE, U_KLOC)   = ROW_SUM_INV_VI(:,IPHASE)*1.0 
                   UGI_COEF_ELE2_ALL(:, IPHASE, U_KLOC2) = ROW_SUM_INV_VJ(:,IPHASE)*1.0 
                END DO

             END IF
          END DO


          DO IPHASE=1,NPHASE
             UDGI_ALL(:, IPHASE) = matmul(INV_VI_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI_ALL(:, IPHASE))
             UDGI2_ALL(:, IPHASE) = matmul(INV_VJ_LOC_OPT_VEL_UPWIND_COEFS(:,:,IPHASE),UDGI2_ALL(:, IPHASE))
          END DO

          NDOTQ  = MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )
          NDOTQ2 = MATMUL( CVNORMX_ALL(:, GI), UDGI2_ALL )


!          print *,'CV_DG_VEL_INT_OPT:',CV_DG_VEL_INT_OPT
!           stop 297

          IF( ABS( CV_DG_VEL_INT_OPT ) == 1 ) THEN
             DT_I=1.0
             DT_J=1.0
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 2) THEN
             DT_I=MAX(1.E-2,LOC_T_I)
             DT_J=MAX(1.E-2,LOC_T_J)
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 3) THEN
             DT_I=LOC_DEN_I*LOC_T_I
             DT_J=LOC_DEN_J*LOC_T_J
          ELSE IF( ABS(CV_DG_VEL_INT_OPT ) == 4) THEN


             ! DG_ELE_UPWIND==1: the upwind method. 
             ! DG_ELE_UPWIND==3: the best optimal upwind frac.


             !NVEC(1)=CVNORMX(GI)
             !NVEC(2)=CVNORMY(GI)
             !NVEC(3)=CVNORMZ(GI)
             NVEC=CVNORMX_ALL(:, GI)
             ABS_CV_NODI_IPHA      = 0.0
             GRAD_ABS_CV_NODI_IPHA = 0.0
             ABS_CV_NODJ_IPHA      = 0.0
             GRAD_ABS_CV_NODJ_IPHA = 0.0
             DO IPHASE = 1, NPHASE
             DO IDIM=1,NDIM
                V_NODI=0.0
                G_NODI=0.0
                V_NODJ=0.0
                G_NODJ=0.0
                DO JDIM=1,NDIM
                   V_NODI = V_NODI + VI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)
                   G_NODI = G_NODI + GI_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)

                   V_NODJ = V_NODJ + VJ_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)
                   G_NODJ = G_NODJ + GJ_LOC_OPT_VEL_UPWIND_COEFS(IDIM,JDIM,IPHASE) * NVEC(JDIM)
                END DO
                ABS_CV_NODI_IPHA(IPHASE)      = ABS_CV_NODI_IPHA(IPHASE) + NVEC(IDIM)*V_NODI
                GRAD_ABS_CV_NODI_IPHA(IPHASE) = GRAD_ABS_CV_NODI_IPHA(IPHASE) + NVEC(IDIM)*G_NODI
                ABS_CV_NODJ_IPHA(IPHASE)      = ABS_CV_NODJ_IPHA(IPHASE) + NVEC(IDIM)*V_NODJ
                GRAD_ABS_CV_NODJ_IPHA(IPHASE) = GRAD_ABS_CV_NODJ_IPHA(IPHASE) + NVEC(IDIM)*G_NODJ
             END DO
             END DO ! PHASE LOOP


             ! Make sure we have some sort of velocity (only needed between elements)...
             ! take the mean of the underlying velocity...
             if(.false.) then
                Q_UNDERLY=( NDOTQ*ABS_CV_NODI_IPHA*MASS_CV_I+ NDOTQ2*ABS_CV_NODJ_IPHA*MASS_CV_J ) &
                     /(MASS_CV_I+MASS_CV_J)
                NDOTQ_TILDE =Q_UNDERLY/ABS_CV_NODI_IPHA
                NDOTQ2_TILDE=Q_UNDERLY/ABS_CV_NODJ_IPHA

             else ! better tested option...

                Q_UNDERLY=0.5*( NDOTQ*ABS_CV_NODI_IPHA + NDOTQ2*ABS_CV_NODJ_IPHA )
                NDOTQ_TILDE =Q_UNDERLY/ABS_CV_NODI_IPHA
                NDOTQ2_TILDE=Q_UNDERLY/ABS_CV_NODJ_IPHA

             endif

             ! These are the new limits of the velocities...
             NDOTQ_KEEP  = NDOTQ_TILDE   ! this is associated with saturation at NODI
             NDOTQ2_KEEP = NDOTQ2_TILDE  ! this is associated with saturation at NODJ

             ! between these limits work out which are associated with the flux limited saturation at this interface. 

             if(ROE_AVE) then

                ! do the Roe average of the rest of the velocity...
                NDOTQ_TILDE  = ( LOC_DEN_I * LOC_T_I * NDOTQ -  &
                     &           LOC_DEN_J * LOC_T_J * NDOTQ2 ) & 
                     / vtolfun( VOLFRA_PORE_ELE*LOC_DEN_I * LOC_T_I - VOLFRA_PORE_ELE2*LOC_DEN_J * LOC_T_J )
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

                         WHERE (0.5*(NDOTQ_TILDE+NDOTQ2_TILDE)>0.0) 
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 * LOC_FEMT(:, CV_KLOC)
                         ELSE WHERE
                            FEMTGI_IPHA = FEMTGI_IPHA + SCVFEN(CV_KLOC,GI) &
                                 * LOC2_FEMT(:, CV_KLOC)
                         END WHERE

                   END IF
                END DO




                WHERE (0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) < 0.0) 
                   INCOME3=1.0
                ELSE WHERE
                   INCOME3=0.0
                END WHERE

                ! Between elements...


                ! redefine so that it detects oscillations...
                abs_tilde1 =  ABS_CV_NODI_IPHA  + 0.5*( LOC_T_J   -  LOC_T_I  ) * GRAD_ABS_CV_NODI_IPHA
                abs_tilde2 =  ABS_CV_NODJ_IPHA  + 0.5*( LOC_T_I   -  LOC_T_J  ) * GRAD_ABS_CV_NODJ_IPHA


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


                WHERE (income3.lt.0.5) ! flux limit
                   wrelax2=wrelax1
                ELSE WHERE
                   wrelax1=wrelax2
                END WHERE


                ! no 4 (more stable than 2 and 3):
                if(between_ele_dg_opt==1) then
                   WHERE ( income3 < 0.5 )  ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*0.5
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*0.5
                   END WHERE

                else if(between_ele_dg_opt==2) then
                   ! no 4 (more stable than 2 and 3) (also make sure we dont violate bounds):
                   WHERE ( income3 < 0.5 )  ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*0.5
                      !              IF(T(CV_NODI_IPHA).LT.0.2) THEN
                      WHERE (ABS_CV_NODI_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*0.5
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*0.5
                      !               IF(T(CV_NODJ_IPHA).LT.0.2) THEN
                      WHERE (ABS_CV_NODJ_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   END WHERE

                   ! no 2(again):
                else if(between_ele_dg_opt==3) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   END WHERE


                   ! no 2(again):
                else if(between_ele_dg_opt==4) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODI_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   ELSE WHERE
                      DT_I = ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODJ_IPHA.GT.1.0E+10) 
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   END WHERE

                   ! no 2(again):
                else if(between_ele_dg_opt==5) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = (1.-wrelax1)*1.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*0.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODI_IPHA.GT.1.0E+10) 
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*0.0  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*1.0  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      WHERE (ABS_CV_NODJ_IPHA.GT.1.0E+10)
                         DT_I = 0.0
                         DT_J = 0.0
                      END WHERE
                   END WHERE

                   ! no 3(again): XXX failed gravity problem...
                else if(between_ele_dg_opt==6) then
                   WHERE ( income3 < 0.5 ) ! upwind
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde1)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde2)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde2)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   END WHERE

                   ! no 2(again,again) XXXXfor the gravity problem:
                else if(between_ele_dg_opt==7) then
                   WHERE ( income3 < 0.5 )  ! upwind
                      DT_I = (1.-wrelax1)*(1.-min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1))  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*min(1.0,ABS_CV_NODJ_IPHA/abs_tilde1)  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ELSE WHERE
                      DT_I = (1.-wrelax1)*min(1.0,ABS_CV_NODI_IPHA/abs_tilde2)  + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                      DT_J = (1.-wrelax2)*(1.-min(1.0,ABS_CV_NODI_IPHA/abs_tilde2))  + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I  &
                           /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   END WHERE

                   ! no 3 soln...
                else if(between_ele_dg_opt==8) then
                   DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA  ) &
                        /(ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA)   + wrelax1*0.5
                   DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA  ) &
                        /(ABS_CV_NODI_IPHA+ABS_CV_NODJ_IPHA)   + wrelax2*0.5

                   ! no 2 soln...
                   !        DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA*MASS_CV_I  ) &
                   !             /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)   + wrelax1*0.5
                   !        DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA*MASS_CV_J  ) &
                   !             /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)   + wrelax2*0.5

                   !        DTOLD_I = ( (1.-wrelaxold1)*ABS_CV_NODI_IPHA*MASS_CV_I  ) &
                   !                /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)  + wrelaxold1*0.5
                   !        DTOLD_J = ( (1.-wrelaxold2)*ABS_CV_NODJ_IPHA*MASS_CV_J  ) &
                   !                /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)  + wrelaxold2*0.5

                else if(between_ele_dg_opt==9) then
                   ! no 0 soln (more stable than 2 & 3)...
                   DT_I = ( (1.-wrelax1)*ABS_CV_NODI_IPHA*MASS_CV_I + wrelax1*ABS_CV_NODJ_IPHA*MASS_CV_J ) &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   DT_J = ( (1.-wrelax2)*ABS_CV_NODJ_IPHA*MASS_CV_J + wrelax2*ABS_CV_NODI_IPHA*MASS_CV_I ) &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   ! no 0 soln (more stable than 2 & 3)...
                else if(between_ele_dg_opt==10) then
                   DT_I = ABS_CV_NODJ_IPHA*MASS_CV_J  &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                   DT_J = ABS_CV_NODI_IPHA*MASS_CV_I  &
                        /(ABS_CV_NODI_IPHA*MASS_CV_I+ABS_CV_NODJ_IPHA*MASS_CV_J)
                endif





                ! END OF IF ( DG_ELE_UPWIND==3 ) THEN ...
             ELSE  ! Low order...
                !if ( abs(T(CV_NODI_IPHA)-T(CV_NODJ_IPHA)).lt.1.e-4 ) then ! low order


                INCOME =0.5*ABS_CV_NODI_IPHA* MASS_CV_I /(0.5*(ABS_CV_NODI_IPHA*MASS_CV_I +ABS_CV_NODJ_IPHA*MASS_CV_J ))

                DO IPHASE = 1, NPHASE
                ! limit the volume fraction based on net letting flux go into certain elements...
                !     if(.false.) then
                if(LIMIT_SAT_BASED_INTERP) then
                   enforce_abs = .false.

                   if(iphase==1) then
                      !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) >= 0.0) THEN 
                      IF(0.5*(NDOTQ(IPHASE)+NDOTQ2(IPHASE)) >= 0.0) THEN 
                         IF(LOC_T_I(IPHASE) < 0.5) then
                            enforce_abs = .true.
                            relax=min( 10.*(0.5-LOC_T_I(IPHASE)), 1.0)
                         ENDIF
                      ELSE
                         IF(LOC_T_J(IPHASE) < 0.5) THEN
                            enforce_abs = .true. 
                            relax=min( 10.*(0.5-LOC_T_J(IPHASE)), 1.0)
                         ENDIF
                      ENDIF
                   endif

                   if(iphase==2) then
                      !  IF(0.5*(NDOTQ_TILDE+NDOTQ2_TILDE) >= 0.0) THEN 
                      IF(0.5*(NDOTQ(IPHASE)+NDOTQ2(IPHASE)) >= 0.0) THEN 
                         IF(LOC_T_I(IPHASE) < 0.5) then
                            enforce_abs = .true.
                            relax=min( 10.*(0.5-LOC_T_I(IPHASE)), 1.0)
                         ENDIF
                      ELSE
                         IF(LOC_T_J(IPHASE) < 0.5) THEN
                            enforce_abs = .true. 
                            relax=min( 10.*(0.5-LOC_T_J(IPHASE)), 1.0)
                         ENDIF
                      ENDIF
                   endif



                   if(enforce_abs) then
                      INCOME(IPHASE) =  (1.-RELAX)*INCOME(IPHASE)  + RELAX*( (1.-INCOME(IPHASE))**2 *(1. - 3.*INCOME(IPHASE)) +  INCOME(IPHASE)*(1.-INCOME(IPHASE))*(1.+3.*INCOME(IPHASE))  ) ! non-lumped
                   endif

                   income(IPHASE)=min(1.0, max(0.0, income(IPHASE)))
                endif

                END DO !PHASE LOOP



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
                DT_I=VOLFRA_PORE_ELE *DT_I
                DT_J=VOLFRA_PORE_ELE2*DT_J
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



          DO IPHASE = 1, NPHASE
             UDGI_INT_ALL(:, IPHASE) = DT_I(IPHASE) * UDGI_ALL(:, IPHASE) + DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)
          END DO

  !        print *,'CV_DG_VEL_INT_OPT:',CV_DG_VEL_INT_OPT
  !        stop 181

          IF( CV_DG_VEL_INT_OPT < 0 ) THEN
   !          FLAbort('2821') ! should not be going in here??
   !          if(got_dt_ij) FLAbort('2611') ! we can not use these CV_DG_VEL_INT_OPT < 0, got_dt_ij=.true. options together.

             !NDOTQ_INT = CVNORMX( GI ) * UDGI_INT + CVNORMY( GI ) * VDGI_INT + &
             !     CVNORMZ( GI ) * WDGI_INT
             NDOTQ_INT = MATMUL( CVNORMX_ALL(:, GI), UDGI_INT_ALL)

             WHERE ( NDOTQ_INT <= 0.0 )  ! Incoming
                !DT_I=1.0
                DT_J=DT_I+DT_J
             ELSE WHERE
                DT_I=DT_I+DT_J
                !DT_J=1.0
             END WHERE

             DO IPHASE = 1, NPHASE
                UDGI_INT_ALL(:, IPHASE) = (DT_I(IPHASE) * UDGI_ALL(:, IPHASE) + DT_J(IPHASE) * UDGI2_ALL(:, IPHASE)) / (DT_I(IPHASE) + DT_J(IPHASE))
             END DO

          END IF


          UDGI_ALL = UDGI_INT_ALL


          DO IPHASE = 1, NPHASE
             UGI_COEF_ELE_ALL(:, IPHASE, :) = DT_I(IPHASE) * UGI_COEF_ELE_ALL(:, IPHASE, :) 

             UGI_COEF_ELE2_ALL(:, IPHASE, :) = DT_J(IPHASE) * UGI_COEF_ELE2_ALL(:, IPHASE, :) 
          END DO


       END IF Conditional_ELE2

    END IF Conditional_SELE


    NDOTQ =  MATMUL( CVNORMX_ALL(:, GI), UDGI_ALL )

    ! Define whether flux is incoming or outgoing, depending on direction of flow
    INCOME = 1.
    WHERE ( NDOTQ >= 0. ) 
          INCOME = 0.
    END WHERE
    !  IF( NDOTQ_tilde >= 0. ) INCOME = 0.
    !IF( NDOTQ > 0. ) INCOME = 0.
    

    RETURN  

CONTAINS


        SUBROUTINE ONVDLIM_ANO_MANY( NFIELD, &
       TDLIM, TDCEN, INCOME, &
       ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT,  &
       TUPWIN, TUPWI2 )
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
    INTEGER, intent( in ) :: NFIELD
    REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM  
    REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2 
    REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
    ! Local variables   
    REAL :: UCIN(NFIELD), UCOU(NFIELD), TDELE(NFIELD), DENOIN(NFIELD), CTILIN(NFIELD), DENOOU(NFIELD), &
         CTILOU(NFIELD), FTILIN(NFIELD), FTILOU(NFIELD)


       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE 

       DENOIN = TOLFUN_MANY( TDELE - TUPWIN )

       UCIN = ETDNEW_PELEOT 
       CTILIN = ( UCIN - TUPWIN ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT

       DENOOU = TOLFUN_MANY( TDELE - TUPWI2 )
       UCOU = ETDNEW_PELE 
       CTILOU = ( UCOU - TUPWI2 ) / DENOOU



       FTILIN = ( TDCEN - TUPWIN ) / DENOIN
       FTILOU = ( TDCEN - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 )



    RETURN

  END SUBROUTINE ONVDLIM_ANO_MANY



        SUBROUTINE ONVDLIM_ANO_MANY_SQRT( NFIELD, &
       TDLIM, TDCEN, INCOME, &
       ETDNEW_PELE, ETDNEW_PELEOT, XI_LIMIT,  &
       TUPWIN, TUPWI2 )
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
    INTEGER, intent( in ) :: NFIELD
    REAL, DIMENSION( NFIELD ), intent( inout ) :: TDLIM  
    REAL, DIMENSION( NFIELD ), intent( in ) :: TDCEN, INCOME, XI_LIMIT, TUPWIN, TUPWI2 
    REAL, DIMENSION( NFIELD ), intent( in ) :: ETDNEW_PELE, ETDNEW_PELEOT
    ! Local variables   
    INTEGER, PARAMETER :: POWER=2
    REAL :: UCIN(NFIELD), UCOU(NFIELD), TDELE(NFIELD), DENOIN(NFIELD), CTILIN(NFIELD), DENOOU(NFIELD), &
         CTILOU(NFIELD), FTILIN(NFIELD), FTILOU(NFIELD)


       ! Calculate normalisation parameters for incomming velocities 
       TDELE = ETDNEW_PELE 

       DENOIN = TOLFUN_MANY( TDELE** POWER - TUPWIN** POWER )

       UCIN = ETDNEW_PELEOT 
       CTILIN = ( UCIN ** POWER- TUPWIN** POWER ) / DENOIN

       ! Calculate normalisation parameters for out going velocities 
       TDELE = ETDNEW_PELEOT

       DENOOU = TOLFUN_MANY( TDELE** POWER - TUPWI2** POWER )
       UCOU = ETDNEW_PELE 
       CTILOU = ( UCOU** POWER - TUPWI2** POWER ) / DENOOU



       FTILIN = ( TDCEN** POWER - TUPWIN** POWER ) / DENOIN
       FTILOU = ( TDCEN** POWER - TUPWI2 ** POWER) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN** POWER + NVDFUNNEW_MANY( FTILIN, CTILIN, XI_LIMIT ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2** POWER + NVDFUNNEW_MANY( FTILOU, CTILOU, XI_LIMIT ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)



    RETURN

  END SUBROUTINE ONVDLIM_ANO_MANY_SQRT




        SUBROUTINE ONVDLIM_ANO( TOTELE, &
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

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD

!     if((PELE==433).and.(PELEOT==435)) then
!          print *,'***INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT:',INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT
!          print *,'***FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM:',FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM
!          print *,'***NOLIMI, FIRORD, TUPWIN2, TUPWI22:', NOLIMI, FIRORD, TUPWIN2, TUPWI22
!     endif

    RETURN

  END SUBROUTINE ONVDLIM_ANO


        SUBROUTINE ONVDLIM_ANO_SQRT( TOTELE, &
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
    INTEGER, PARAMETER :: POWER=2
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
       UCIN = ETDNEW_PELE** POWER
       DENOIN = 1.
       CTILIN = 0.

       ! Calculate normalisation parameters for out going velocities 
       TUPWI2 = ETDNEW_PELE** POWER
       UCOU = ETDNEW_PELE ** POWER
       DENOOU = 1.
       CTILOU = 0.

    END IF Conditional_PELEOT

    Conditional_FIRORD: IF( FIRORD ) THEN ! Velocity is pointing into element

       ! Velocity is going out of element
       TDLIM = INCOME * UCIN + ( 1.0 - INCOME ) * UCOU 
       TDLIM = MAX(TDLIM,0.0) ** (1.0/POWER)

    ELSE

       FTILIN = ( MAX(0.0,TDCEN) ** POWER - TUPWIN ) / DENOIN
       FTILOU = ( MAX(0.0,TDCEN) ** POWER - TUPWI2 ) / DENOOU

       ! Velocity is going out of element
       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, COURANT_OR_MINUS_ONE ) * DENOIN ) &
            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, COURANT_OR_MINUS_ONE ) * DENOOU )
!       TDLIM =        INCOME   * ( TUPWIN + NVDFUNNEW( FTILIN, CTILIN, -1. ) * DENOIN ) &
!            + ( 1.0 - INCOME ) * ( TUPWI2 + NVDFUNNEW( FTILOU, CTILOU, -1. ) * DENOOU )

       TDLIM = MAX( TDLIM, 0.0 ) ** (1.0/POWER)

    ENDIF Conditional_FIRORD
!       if(( PELE==433).and.( PELEOT==435))  print *,'++ firord,TDLIM, INCOME, TUPWI2, FTILOU, CTILOU, COURANT_OR_MINUS_ONE:', &
!                                                        firord,TDLIM, INCOME, TUPWI2, FTILOU, CTILOU, COURANT_OR_MINUS_ONE
!     if((PELE==433).and.(PELEOT==435)) then
!          print *,'---new INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT:',INCOME,TUPWIN,TUPWI2,TDCEN,ETDNEW_PELE,ETDNEW_PELEOT
!          print *,'---new FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM:',FTILOU, CTILOU, COURANT_OR_MINUS_ONE, TDLIM
!          print *,'---new NOLIMI, FIRORD, TUPWIN2, TUPWI22:', NOLIMI, FIRORD, TUPWIN2, TUPWI22
!     endif
    RETURN

  END SUBROUTINE ONVDLIM_ANO_SQRT

  END SUBROUTINE GET_INT_VEL_POROUS_VEL



  function vtolfun(val) result(v_tolfun)

    implicit none
    real, dimension(:), intent(in) :: val
    real, dimension(size(val)) :: v_tolfun
    ! Local
    real, parameter :: tolerance = 1.e-10

    v_tolfun = sign( 1.0, val ) * max( tolerance, abs(val) )

    return

  end function vtolfun

!
! 
! 
! 
    real function get_relperm_epsilon(sat,iphase, Sr1, Sr2, kr_exp, krmax)
        Implicit none
        real, intent( in ) :: Sat, Sr1, Sr2, kr_exp, krmax
        integer, intent(in) :: iphase

            get_relperm_epsilon = min(max(1d-20, Krmax*( ( Sat - Sr1) / ( 1. - Sr1 - Sr2 ) ** kr_exp)), Krmax)
            !Make sure that the relperm is between bounds
            !Lower value just to make sure we do not divide by zero.

    end function get_relperm_epsilon
! 
! 
    real function inv_get_relperm_epsilon(sat,iphase, Sr1, Sr2, kr_exp, krmax)
        Implicit none
        real, intent( in ) :: Sat, Sr1, Sr2, kr_exp, krmax
        integer, intent(in) :: iphase

            inv_get_relperm_epsilon =  1.0/ min(max(1d-20, Krmax*( ( Sat - Sr1) / ( 1. - Sr1 - Sr2 ) ** kr_exp)), Krmax)
            !Make sure that the relperm is between bounds
            !Lower value just to make sure we do not divide by zero.

    end function inv_get_relperm_epsilon
! -----------------------------------------------------------------------------

! 
! -----------------------------------------------------------------------------

  end module cv_advection

