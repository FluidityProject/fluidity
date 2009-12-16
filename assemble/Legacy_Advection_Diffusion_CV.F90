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

module legacy_advection_diffusion_CV
  !!< This module contains the assembly subroutines for advection
  !!< using control volumes
  
  use allsorts
  use fields
  use field_derivatives
  use state_module
  use futils
  use fetools
  use spud
  use legacy_cv_shape_functions
  use legacy_cv_numbering
  use mesh_connections
  use shape_transformations
  use diff3d_module
  use hart3d_allsorts
  use surface_integration
  use flcomms_module
  use detnlxr_module
  implicit none

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public :: ASSEMBLE_MOMENTUM_ADVECTION_CV, ASSEMBLE_FIELD_ADVECTION_CV

contains


    SUBROUTINE ASSEMBLE_FIELD_ADVECTION_CV(&
        NONODS,VNONOD,XNONOD,TOTELE,NLOC,NGI, &
        SNLOC, STOTEL, &
        NDGLNO,XNDGLN,VONDGL, &
        SNDGLN, &
        FINDRM,COLM,NCOLM,CENTRM, & 
        BIGM,RHS,&
        X,Y,Z, &
        N,NLX,NLY,NLZ, WEIGHT, &
        NU,NV,NW, UG,VG,WG, T,TOLD,DENSTY,DENOLD, &
        DISOPT,DT,THETA,BETA, D3,DCYL, &
        PARA,halo_tag,NNODP, ULTRAC &
        , NOBCT, BCT1, BCT2 &
        , NITS, NDISOT2, NEWMES2, MATRIX &
        , NSUBNVLOC,NUSUB,NVSUB,NWSUB, &
        option_path &
        )
      !     ===================================================================
      !          !!!---THIS SUBROUTINE REPLACES THE OLD ADVOCEAN2---!!!
      !
      !     In this subroutine the advection terms in the advection-diffusion
      !     equation for a user-defined field variable (T) are calculated using 
      !     a control volume formulation.  They are assembled in MATRIX and RHS.  
      !
      !     They are then added to the rest of the advection-diffusion terms, 
      !     which were previously calculated in HART3D (or an alternative, 
      !     depending on geometry).
      !
      !     This routine is only called if the GEM option NDISOT(IT) is
      !     non zero; the last digit of NDISOT(IT) is passed in to this 
      !     routine as DISOPT.
      !
      !     The assembly of MATRIX and RHS is done in subroutine ADVOCEAN
      !     which is included in this file, below.  Please see this routine for
      !     more detailed comments on the advection routine.
      !
      !     This routine is called from ADVDIF
      !
      !     For more details see (enter manual/wiki page here...)
      !     
      !     Description                                   Programmer      Date
      !     ==================================================================
      !     Original version..................................CCP     Unknown!
      !     General clean up, improving efficiency and
      !     adding comments...................................GSC   2006-08-10
      !     Rename advoceanfld from advocean2 and rewritten 
      !     as free-form fortran 90 (including allocation)....GSC   2006-08-11
      !
      !***********************************************************************  
    
      ! Inputs/Outputs
      INTEGER :: NONODS,VNONOD,XNONOD,TOTELE,NLOC,NGI
      INTEGER, INTENT(IN) :: SNLOC, STOTEL
      INTEGER :: NDGLNO(TOTELE*NLOC),XNDGLN(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
      INTEGER, INTENT(IN) :: SNDGLN(STOTEL*SNLOC)
      INTEGER :: FINDRM(NONODS+1),NCOLM,COLM(NCOLM),CENTRM(NONODS)
      REAL    :: BIGM(NCOLM)
      REAL    :: RHS(NONODS)
      REAL    :: X(XNONOD),Y(XNONOD),Z(XNONOD)
      REAL    :: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL    :: WEIGHT(NGI) 
      INTEGER NSUBNVLOC
      REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
      REAL NWSUB(TOTELE*NSUBNVLOC)
      REAL    :: NU(VNONOD),NV(VNONOD),NW(VNONOD)
      REAL    :: UG(VNONOD),VG(VNONOD),WG(VNONOD)
      REAL    :: T(NONODS),TOLD(NONODS)
      REAL    :: DENSTY(NONODS),DENOLD(NONODS)
      INTEGER :: DISOPT
      REAL    :: DT,THETA,BETA 
      LOGICAL :: D3,DCYL
      INTEGER :: PARA,NNODP !Are these really inputs/outputs?... yes
      INTEGER :: halo_tag
      
      ! Element type informtation 
      INTEGER :: SCVNGI
      INTEGER :: ELETYP,OPTELM
      INTEGER :: CVNGI,CVNLOC,SELETYP          !Dummy integers output from SETELM
      LOGICAL :: CVD3,CVDCYL                                 !Dummy logicals output from SETELM
      INTEGER :: SNLOC2, SNGI
      
      INTEGER :: NOD,COUNT  !Loop indices 
      LOGICAL :: REDQUAD
      PARAMETER(REDQUAD=.FALSE.)
    
      INTEGER, INTENT(IN) :: NOBCT, BCT2(NOBCT)
      REAL, INTENT(IN) ::    BCT1(NOBCT)
      
      LOGICAL :: ULTRAC

      REAL :: CHANGE
    
      REAL, INTENT(INOUT) :: MATRIX(NCOLM)
      INTEGER, INTENT(IN) :: NITS, NDISOT2
      LOGICAL, INTENT(IN) :: NEWMES2
    
      character(len=*) :: option_path

      ewrite(2,*) 'NLOC, DISOPT, DT, THETA, BETA, D3 = ', &
                  NLOC, DISOPT, DT, THETA, BETA, D3
      CHANGE = 0.
      CALL COMPER(PARA, DENSTY, DENOLD, NONODS, CHANGE)
      ewrite(2,*)'DIFFERENCE CONSTRUCT_ADVECTION_DIFFUSION_CV (MATDENPT, MATDENOLD) = ', CHANGE
    !   DO NOD = 1, NONODS
    !     ewrite(2,*)DENSTY(NOD), DENOLD(NOD)
    !   ENDDO
      EWRITE_MINMAX(DENSTY)
      EWRITE_MINMAX(DENOLD)
    
      ! ============================= FIND SCVNGI ===============================
      ! Define the element type
      IF(D3) THEN                      !  element type                | optelm
          IF(NLOC.EQ.4)  OPTELM=5       ! -------------------------------------
          IF(NLOC.EQ.8)  OPTELM=7       ! - linear        triangles        1
          IF(NLOC.EQ.10) OPTELM=6       ! - quadratic     triangles        2
          IF(NLOC.EQ.27) OPTELM=8       ! - bi-linear     quadrilaterals   3
      ELSE                             ! - bi-quadratic  quadrilaterals   4
          IF(NLOC.EQ.3) OPTELM=1        ! - linear        tetrahedra       5
          IF(NLOC.EQ.4) OPTELM=3        ! - quadratic     tetrahedra       6
          IF(NLOC.EQ.6) OPTELM=2        ! - tri-linear    hexahedra        7
          IF(NLOC.EQ.9) OPTELM=4        ! - tri-quadratic hexahedra        8
      ENDIF
    
      ewrite(2,*)'OPTELM = ', OPTELM
      
      ! Use the SETELM routine to define the number of gauss points 
      ! on the control volume surface (SCVNGI), based on the element type
      ! SETELM takes the inputs OPTELM and REDQUAD, all other
      ! terms are outputs and are not used
      CALL SETELM( ELETYP, CVNGI, CVNLOC, &
      ! - INTEGERS
            OPTELM, SELETYP, &
            SNGI,   SNLOC2,   SCVNGI, &
      ! - LOGICALS
            CVD3,     CVDCYL,    REDQUAD   ) 
    
      ewrite(2,*)'SNLOC, SNLOC2, SNGI', SNLOC, SNLOC2, SNGI
    
      ! make sure we are using the correct shape function...
      IF(CVNLOC.NE.NLOC) STOP 'CVNLOC.NE.NLOC IN ASSEMBLE_FIELD_ADVECTION_CV'
      IF(SNLOC2.NE.SNLOC) STOP 'SNLOC2.NE.SNLOC IN ASSEMBLE_FIELD_ADVECTION_CV'
    
      ! ========================= COMPUTE MATRIX AND RHS ========================
      ! Allocate memory to the big matrix...
      
      ! Call the routine to assemble the matrix of advection terms for T
      ! Returns MATRIX and RHS...
      CALL CONSTRUCT_ADVECTION_DIFFUSION_CV(NONODS,VNONOD,XNONOD,TOTELE, &
            NLOC,NGI,SCVNGI,CVNGI,ELETYP,SELETYP, &
            SNLOC, SNGI, STOTEL, &
            NDGLNO,XNDGLN,VONDGL, &
            SNDGLN, &
            FINDRM,COLM,NCOLM,CENTRM, & 
            MATRIX,RHS, &
            X,Y,Z, &
            N,NLX,NLY,NLZ, WEIGHT, &
            NU,NV,NW, UG,VG,WG, T,TOLD,DENSTY,DENOLD, &
            DISOPT,DT,THETA,BETA, D3,DCYL, &
            PARA,halo_tag,NNODP &
            , ULTRAC &
            , NOBCT, BCT1, BCT2 &
            , NITS, NDISOT2, NEWMES2 &
            , NSUBNVLOC,NUSUB,NVSUB,NWSUB, option_path &
            )
    
      ! ===================== COMBINE WITH DIFFUSION TERMS, ETC ================
      ! Combine advection terms computed in CONSTRUCT_ADVECTION_DIFFUSION_CV with the
      ! diffusion terms already assembled in BIGM, and 
      ! subtract A*TOLD from RHS
      
      DO NOD=1,NONODS
          DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
            RHS(NOD)=RHS(NOD)-MATRIX(COUNT)*TOLD(COLM(COUNT)) !&
!                              -(1./DT)*BIGM(COUNT)*&
!                               (DENSTY(COLM(COUNT))-DENOLD(COLM(COUNT)))*&
!                               TOLD(COLM(COUNT))

!             BIGM(COUNT)=BIGM(COUNT)*DENSTY(COLM(COUNT))+DT*MATRIX(COUNT)
          END DO
      END DO
      
      DO COUNT=1,NCOLM
          BIGM(COUNT)=BIGM(COUNT)+DT*MATRIX(COUNT)
      END DO
    
      RETURN
    END SUBROUTINE ASSEMBLE_FIELD_ADVECTION_CV
    



    SUBROUTINE CONSTRUCT_ADVECTION_DIFFUSION_CV( &
        NONODS,VNONOD,XNONOD,TOTELE, &
        NLOC,NGI,SCVNGI,CVNGI,ELETYP,SELETYP, &
        SNLOC, SNGI, STOTEL, &
        NDGLNO,XNDGLN,VONDGL, &
        SNDGLN, &
        FINDRM,COLM,NCOLM,CENTRM, & 
        MATRIX,RHS,&
        X,Y,Z, &
        N,NLX,NLY,NLZ, WEIGHT, &
        NU,NV,NW, UG,VG,WG, T,TOLD,DENSTY,DENOLD, &
        DISOPT,DT,THETA,BETA, D3,DCYL, &
        PARA,halo_tag,NNODP &
        , ULTRAC &
        , NOBCT, BCT1, BCT2 &
        , NITS, NDISOT2, NEWMES2 &
        , NSUBNVLOC,NUSUB,NVSUB,NWSUB, option_path &
        )
      !  =====================================================================
      !     In this subroutine the advection terms in the advection-diffusion
      !     equation (in the matrix and RHS) are calculated as MATRIX and RHS.  
      !     They are then added to the rest of the advection-diffusion terms, 
      !     which were previously calculated in HART3D (or an alternative, 
      !     depending on geometry).
      !
      !     This routine is called to advect a user-defined field variable
      !     from ASSEMBLE_FIELD_ADVECTION_CV, or three times in succession to advect the 
      !     x/y/z momentum components from ASSEMBLE_MOMENTUM_ADVECTION_CV.  So, T and TOLD
      !     either represent the old and new field values, or the old and new
      !     U/V/W... 
      !
      !     !!! However, ASSEMBLE_MOMENTUM_ADVECTION_CV appears to be obsolete, so it is
      !     'safe' to assume that this routine is called from ASSEMBLE_FIELD_ADVECTION_CV,
      !     and that T always represents a field variable...            !!!
      !
      !     This routine is only called if the GEM option NDISOT(IT) is
      !     non zero; the last digit of NDISOT(IT) is passed in to this 
      !     routine as DISOPT.
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
      !        THETA value, or a non-linear THETA value that is determined 
      !        automatically.  
      !      - The face value may be limited using a univeral-limiter-type
      !        scheme, or a limited-downwind scheme that is ideal for INTERFACE
      !        TRACKING.  Alternatively no limiting can be applied.  
      !     
      !     These options are defined by the value of DISOPT, which corresponds 
      !     to the clast digit of the GEM option NDISOT for the field in question.
      !
      !     DISOPT=discretisation option in space and time
      !     ------------------------------------------------------------------
      !     DISOPT   Method for face-value est.    Time-stepping     Limiting   
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
      !     is selected (DISOPT=8/9).  ***In general ALOLIM appears to be set to 1 (GSC)
      !     
      !     NOTE: ANONVDLIM only works for TETS; for all other element types 
      !     ONVDLIM is used.
      !
      !
      !     IMPORTANT INPUTS:
      !     ----------------
      !     
      !     MATRIX   - Matrix for assembling the advection terms (empty on input)
      !     RHS      - Right-hand side vector for advection-diffusion terms
      !     X,Y,Z    - Node co-ordinates
      !     NU,NV,NZ - Nodal velocity components
      !     UG,VG,WG - Grid velocities
      !     T,TOLD   - New and old advected field values at nodes
      !     DENSTY,  - New and old "density" at nodes, which is actually a constant
      !     DENOLD     multiplying the advection diffusion equation for the field
      !     DISOPT   - The discretisation/limiting option (see above)
      !     DT       - The time step
      !     THETA    - The time-stepping discretisation parameter
      !     BETA     - Conservative(1.)/non-conservative(0.) flag
      !     D3,DCYL  - Geometry flags (D3 for 3D, DCYL for 2D cylindral)
      !     ELETYP   - Integer flag definining element type     
      !
      !     PARA,halo_tag,NNODP - Parallel parameters
      !
      !     IMPORTANT OUTPUTS:
      !     -----------------
      !
      !     MATRIX   - Matrix updated to include the advection terms
      !     RHS      - Right-hand side vector updated to include advection terms
      !
      !
      !     IMPORTANT LOCAL PARAMETERS:
      !     --------------------------
      !
      !     TIMOPT    - Temporal discretisation option, derived from DISOPT.
      !                (1 for non-linear theta; 0 for theta specified (THETA))
      !     ALOLIM    - Anisotropic limiting flag
      !                (1 to use best anisotropic limiter; 0 to turn off)
      !     ANISOLIM  - Logical flag true if ANONVDLIM is to be used
      !
      !     TUPWIND,  - Matrix containing the interpolated upwind value of
      !     TUPWINDOLD  the field variable (old and new) for universal limiter
      !
      !
      !
      !     Called from ADVDIF if GEM option NDISOT(IT) is non-zero
      !
      !     For more details see (enter manual/wiki page here...)
      !     
      !     Description                                   Programmer      Date
      !     ==================================================================
      !     Original version..................................CCP     Unknown!
      !     General clean up, improving efficiency and
      !     adding comments...................................GSC   2006-08-10
      !     Combined with advocean2 and rewritten as 
      !     free-form fortran 90 (including allocation).......GSC   2006-08-11
      !
      !***********************************************************************
    
      ! Inputs/Outputs
      INTEGER ISUBMOD,NQLOC
      LOGICAL QUAD_INTERP,PERT_APPROACH,COURANT_LIM
      LOGICAL ELE_GRAD
      LOGICAL LIMIT_PERT,CENT_HIGH_T,STAN_LIM,SEVER_LINK_GE2
! Adjustment for accurate non-conservative advectob ib straticied flows
! If QUAD_INTERP then use a quadratic interpolation method to obtain high order fluxes.
! If PERT_APPROACH then limit a pertabation from a particular plane defined 
! by the gradient at down wind node NODD and a point on the plane T(NODC). 
! If COURANT_LIM then include Courant no. in limiter for stratified flows (DISOPT=2 above)
! If ELE_GRA use the element gradeint to help caluclate the limiting.
! If LIMIT_PERT then limit the pertabation only
! If CENT_HIGH_T then use mid point rule to obtain the high order flux
! If STAN_LIM then use standard limiting approach.
! If SEVER_LINK_GE2 avoid applying limiting when there is more than 2 
! boundary faces on an element.
      PARAMETER(COURANT_LIM=.true.)
      PARAMETER(ELE_GRAD=.true.)
      PARAMETER(LIMIT_PERT=.true.,CENT_HIGH_T=.false.,STAN_LIM=.false.)
!      PARAMETER(LIMIT_PERT=.true.,CENT_HIGH_T=.true.,STAN_LIM=.true.)
!      PARAMETER(LIMIT_PERT=.false.,CENT_HIGH_T=.false.,STAN_LIM=.false.)
      PARAMETER(SEVER_LINK_GE2=.false.) 
      PARAMETER(ISUBMOD=1,NQLOC=10)
       INTEGER IFACE,NFACE
       PARAMETER(NFACE=4)
       INTEGER LOCLIST(NFACE,3)
! IF ISUBMOD.NE.0 then switch on submodel advection velocity if it is sent down.
      INTEGER :: NONODS,VNONOD,XNONOD,TOTELE,NLOC,NGI,SCVNGI,CVNGI,ELETYP,SELETYP
      INTEGER, INTENT(IN) :: SNLOC, SNGI, STOTEL
      INTEGER :: NDGLNO(TOTELE*NLOC),XNDGLN(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
      INTEGER, INTENT(IN) :: SNDGLN(STOTEL*SNLOC)
      INTEGER :: FINDRM(NONODS+1),NCOLM,COLM(NCOLM),CENTRM(NONODS)
      REAL    :: RHS(NONODS)
      REAL    :: X(XNONOD),Y(XNONOD),Z(XNONOD)
      REAL, INTENT(IN) :: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL :: SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      REAL, INTENT(IN) :: WEIGHT(NGI) 
      REAL :: SWEIGH(SNGI) 
      INTEGER NSUBNVLOC
      REAL NUSUB(TOTELE*NSUBNVLOC),NVSUB(TOTELE*NSUBNVLOC)
      REAL NWSUB(TOTELE*NSUBNVLOC)
      REAL    :: NU(VNONOD),NV(VNONOD),NW(VNONOD)
      REAL    :: UG(VNONOD),VG(VNONOD),WG(VNONOD)
      REAL    :: T(NONODS),TOLD(NONODS)
      REAL    :: DENSTY(NONODS),DENOLD(NONODS)
      INTEGER :: DISOPT
      REAL    :: DT,THETA,BETA 
      LOGICAL :: D3,DCYL
      INTEGER :: PARA,NNODP
      INTEGER :: halo_tag
      
      INTEGER, INTENT(IN) :: NOBCT, BCT2(NOBCT)
      REAL, INTENT(IN) ::    BCT1(NOBCT)

      ! Local allocatable arrays memory
      REAL, ALLOCATABLE :: SVN(:,:), SVNLX(:,:), SVNLY(:,:) 
      REAL, ALLOCATABLE :: SVN_QUAD(:,:)
      REAL, ALLOCATABLE :: M(:,:)
      REAL, ALLOCATABLE :: SVWEIGH(:)
      REAL, ALLOCATABLE :: TUPWIND(:),TUPWINDOLD(:)
      REAL, ALLOCATABLE :: DUPWIND(:),DUPWINDOLD(:)
      REAL, ALLOCATABLE :: TMAX(:),TMIN(:)
      REAL, ALLOCATABLE :: TOLDMAX(:),TOLDMIN(:)
      REAL, ALLOCATABLE :: DMAX(:),DMIN(:)
      REAL, ALLOCATABLE :: DOLDMAX(:),DOLDMIN(:)
      REAL, ALLOCATABLE :: TMAX2(:),TMIN2(:)
      REAL, ALLOCATABLE :: TOLDMAX2(:),TOLDMIN2(:)
      REAL, ALLOCATABLE :: DIFXT(:), DIFYT(:), DIFZT(:)
      REAL, ALLOCATABLE :: DIFXTOLD(:),DIFYTOLD(:),DIFZTOLD(:)
      REAL, ALLOCATABLE :: ML(:) 
      REAL, ALLOCATABLE :: CVDETWEI(:)
      REAL, ALLOCATABLE :: CVNORMX(:),CVNORMY(:),CVNORMZ(:)   
      REAL, ALLOCATABLE :: DETWEI(:),RA(:)
      REAL, ALLOCATABLE :: SDETWEI(:)
      REAL, ALLOCATABLE :: NX(:,:),NY(:,:),NZ(:,:) 
      INTEGER, ALLOCATABLE :: NEILOC(:,:),SNEILOC(:,:),COLGPTS(:),FINDGPTS(:)
      INTEGER, ALLOCATABLE :: FACLOC(:,:), FACNOD(:,:)
      INTEGER, ALLOCATABLE :: NEIELE(:,:)
    
      INTEGER, ALLOCATABLE :: FINDELE(:),COLELE(:),TEMPCOLELE(:)
      INTEGER, ALLOCATABLE :: NLIST(:),INLIST(:)
      
      REAL, ALLOCATABLE :: DM1(:), OLDDM1(:), DENDM1(:), DENOLDDM1(:)
      
      REAL, ALLOCATABLE :: SNORMXN(:), SNORMYN(:), SNORMZN(:)
      REAL :: SNORMX, SNORMY, SNORMZ, SAREA
      
      INTEGER :: NOD
      
      ! Local variables...
    
!       ! Loop indices and pointers
      INTEGER :: NODI,NODJ,COUNT,JCOUNT,KCOUNT,KKCOUNT, I
      INTEGER :: ELE,ELE2,ILOC,JLOC,GI,GCOUNT,II
      INTEGER :: SELE, SINOD, SILOC, SJNOD, SJLOC, SGI
      INTEGER :: NCOLGPTS,KLOC,NODK,VNODK,VNODI,XNODI,XNODJ
      INTEGER :: NODC,NODD,XNODC,XNODD
      
        REAL XD,YD,ZD,TDX,TDY,TDZ,TOLDDX,TOLDDY,TOLDDZ
        REAL DVECX,DVECY,DVECZ,DVEC,ULENG,UDNORM,VDNORM,WDNORM
        REAL HVECX,HVECY,HVECZ,TMID,TOLDMID,VOLUME,UDSIG
        REAL FULLLIM,FULLLIMOLD
        REAL OCVNORMX,OCVNORMY,OCVNORMZ,OFEMT,OFEMTOLD
        
        REAL GRAD_ELE_T_X,GRAD_ELE_T_Y,GRAD_ELE_T_Z  
        REAL GRAD_ELE_TOLD_X,GRAD_ELE_TOLD_Y,GRAD_ELE_TOLD_Z 
        REAL REM_GRAD_ELE_T_X,REM_GRAD_ELE_T_Y,REM_GRAD_ELE_T_Z   
        REAL REM_GRAD_ELE_TOLD_X,REM_GRAD_ELE_TOLD_Y,REM_GRAD_ELE_TOLD_Z
    
      !Options and variables...
      INTEGER :: TIMOPT,OPINTE
      REAL    :: UDGI,VDGI,WDGI,NDOTQ,INCOME,HDC
      REAL    :: SNDOTQ(SNGI)
      REAL    :: FVT,FVTOLD,FVD,FVDOLD
      REAL    :: FEMT,FEMTOLD,FEMD,FEMDOLD
      REAL    :: LIMT,LIMTOLD,LIMD,LIMDOLD
      REAL    :: GF,PINVTH,QINVTH,FTHETA
      REAL    :: PNTX,PNTY,PNTZ
      REAL    :: VTHETA,COURANT
      REAL    :: GPSIL,GPSIOLDL,FPSIL,FPSIOLDL,RN
      LOGICAL :: FIRSTORD,NOLIMI
      INTEGER :: MXNCOLEL, NCOLEL
      REAL MX_GRAD_PLAN_X_T,MX_GRAD_PLAN_Y_T,MX_GRAD_PLAN_Z_T
      REAL MX_GRAD_PLAN_X_TOLD,MX_GRAD_PLAN_Y_TOLD,MX_GRAD_PLAN_Z_TOLD
      REAL MIN_GRAD_PLAN_X_T,MIN_GRAD_PLAN_Y_T,MIN_GRAD_PLAN_Z_T
      REAL MIN_GRAD_PLAN_X_TOLD,MIN_GRAD_PLAN_Y_TOLD,MIN_GRAD_PLAN_Z_TOLD
    
      ! Local parameters...
      INTEGER :: ALOLIM,CROSTR
      LOGICAL :: GRADBO
      PARAMETER(CROSTR=0, GRADBO = .FALSE.)
    
      LOGICAL :: GOTPNTS,ANISOLIM,ULTRAC,ULTRACIN
      REAL :: TIN
      
      REAL XGI,YGI,ZGI,XDIR,YDIR,ZDIR
      REAL GRAD_PLAN_X_T,GRAD_PLAN_Y_T,GRAD_PLAN_Z_T
      REAL GRAD_PLAN_X_TOLD,GRAD_PLAN_Y_TOLD,GRAD_PLAN_Z_TOLD
      REAL XPT_PLAN,YPT_PLAN,ZPT_PLAN
      REAL EST_FEMT,PERT_FEMT,EST_FEMTOLD,PERT_FEMTOLD
      REAL EST_T_NOD,PERT_T_MXMIN,  EST_TOLD_NOD,PERT_TOLD_MXMIN
      REAL EST_T_NODD,EST_TOLD_NODD
      REAL PERT_LIMT,PERT_LIMTOLD    
      REAL GRAD_LINE_X_T,GRAD_LINE_Y_T,GRAD_LINE_Z_T
      REAL GRAD_LINE_X_TOLD,GRAD_LINE_Y_TOLD,GRAD_LINE_Z_TOLD
      REAL ALONG_GRAD_PLAN_X_T,ALONG_GRAD_PLAN_Y_T,ALONG_GRAD_PLAN_Z_T
      REAL ALONG_GRAD_PLAN_X_TOLD,ALONG_GRAD_PLAN_Y_TOLD,ALONG_GRAD_PLAN_Z_TOLD
      REAL REM_GRAD_PLAN_X_T,REM_GRAD_PLAN_Y_T,REM_GRAD_PLAN_Z_T
      REAL REM_GRAD_PLAN_X_TOLD,REM_GRAD_PLAN_Y_TOLD,REM_GRAD_PLAN_Z_TOLD
      REAL ALONG_GRAD_PLAN_X_T2,ALONG_GRAD_PLAN_Y_T2,ALONG_GRAD_PLAN_Z_T2
      REAL ALONG_GRAD_PLAN_X_TOLD2,ALONG_GRAD_PLAN_Y_TOLD2,ALONG_GRAD_PLAN_Z_TOLD2
      INTEGER ILOOP,NODCD
    
      REAL, INTENT(INOUT) :: MATRIX(NCOLM)
      INTEGER, INTENT(IN) :: NITS, NDISOT2
      LOGICAL, INTENT(IN) :: NEWMES2
    
      LOGICAL :: GETMAT, USE_INNER_STOR, STORE_ELE, RET_STORE_ELE
      LOGICAL :: BOUND, REFLECT
      REAL TQUAD(NQLOC),TOLDQUAD(NQLOC)
    
      LOGICAL, SAVE :: FIRSTVIS=.TRUE.
    
      INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: ELEMATPSI
      REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: ELEMATWEI
      INTEGER :: DUMMYINT(0),ICOUNT,ISUF
      REAL :: DUMMYREAL(0)
      
      REAL, ALLOCATABLE,DIMENSION(:) :: PERT_T,PERT_TOLD
      REAL, ALLOCATABLE,DIMENSION(:) :: PERT_TMAX,PERT_TMIN,PERT_TOLDMAX,PERT_TOLDMIN

      INTEGER, ALLOCATABLE,DIMENSION(:) :: FINELE,COLELE4
      LOGICAL, ALLOCATABLE,DIMENSION(:) :: ONBOU,SUFS_ONBOU_GE2
      LOGICAL, ALLOCATABLE,DIMENSION(:) :: SUF_LINK

      character(len=*) :: option_path

      real :: ptheta, btheta

      ewrite(2,*)'ENTERING CONSTRUCT_ADVECTION_DIFFUSION_CV()'

      ptheta = 1.0
      btheta = 1.0
    
      ! Allocate memory for the control volume surface shape functions, etc.
      ALLOCATE(   SVN(NLOC,SCVNGI) )
      ALLOCATE( SVNLX(NLOC,SCVNGI) )
      ALLOCATE( SVNLY(NLOC,SCVNGI) )
      ALLOCATE(     M(NLOC,CVNGI) )
      ALLOCATE(    SVWEIGH(SCVNGI) )
      ALLOCATE(   CVDETWEI(SCVNGI) )
      ALLOCATE(    CVNORMX(SCVNGI) )
      ALLOCATE(    CVNORMY(SCVNGI) )
      ALLOCATE(    CVNORMZ(SCVNGI) )
      ALLOCATE(  NEILOC(NLOC,SCVNGI) )
      ALLOCATE(  SNEILOC(SNLOC,SNGI) )
      ALLOCATE( COLGPTS(NLOC*SCVNGI) ) !The size of this vector is over-estimated
      ALLOCATE(    FINDGPTS(NLOC+1) )
      ALLOCATE(   SVN_QUAD(NQLOC,SCVNGI) )
      
        ALLOCATE(  DETWEI(NGI) )
        ALLOCATE(      RA(NGI) )
        ALLOCATE( NX(NLOC,NGI) )
        ALLOCATE( NY(NLOC,NGI) )
        ALLOCATE( NZ(NLOC,NGI) )
    
      ALOLIM = 1
      IF(NDISOT2.LT.0) ALOLIM=0
    
      IF((ALOLIM.EQ.0).AND.(DISOPT.EQ.2)) THEN
           PERT_APPROACH=.TRUE.
           QUAD_INTERP=.TRUE.
      ELSE
           PERT_APPROACH=.FALSE.
           QUAD_INTERP=.FALSE.
      ENDIF
        
      IF(PERT_APPROACH) THEN
        ALLOCATE(PERT_T(NONODS))
        ALLOCATE(PERT_TOLD(NONODS))
        ALLOCATE(PERT_TMAX(NONODS))
        ALLOCATE(PERT_TMIN(NONODS))
        ALLOCATE(PERT_TOLDMAX(NONODS))
        ALLOCATE(PERT_TOLDMIN(NONODS))
      ENDIF
    
      GETMAT = .FALSE.
      IF(NITS.EQ.1) GETMAT = .TRUE.
    
      IF(GETMAT) MATRIX = 0.0
    
      USE_INNER_STOR = .FALSE.
      IF(ABS(NDISOT2/1000).EQ.1) USE_INNER_STOR = .TRUE.
      
      ! Define whether original limiting scheme or anisotropic/interface tracking
      ! limiting scheme is being used...
      !    Use ANONVDLIM for TETS ONLY if:
      !    - ANISOTROPIC LIMITING is desired (parameter ALOLIM=1)
      !    - Or, interface tracking NDISOTT option is set
      ANISOLIM=.FALSE.
      IF( ((ALOLIM.EQ.1).or.(DISOPT.EQ.8).or.(DISOPT.EQ.9)).AND.(NLOC.EQ.4)) ANISOLIM=.TRUE.
    
      ! Allocate memory and find upwind field values for limiting...
      IF (ANISOLIM) THEN
    
        ! Over-estimate the size of the COLELE array
        MXNCOLEL=20*TOTELE+500
        
        ALLOCATE( FINDELE(NONODS+1) )
        ALLOCATE( TEMPCOLELE(MXNCOLEL) )
        
        ALLOCATE(   NLIST(NONODS) )
        ALLOCATE(  INLIST(NONODS) )
        
        ! Calculate node element list - moved from (I)FINPTS
        CALL PHILNODELE(NONODS,FINDELE,TEMPCOLELE, &
                        NCOLEL,MXNCOLEL, &
                        TOTELE,NLOC,NDGLNO, &
                        NLIST,INLIST)
        
        ALLOCATE( COLELE(NCOLEL) )
        
        DO I = 1, NCOLEL
        
          COLELE(I) = TEMPCOLELE(I)
        
        ENDDO
        
        DEALLOCATE( TEMPCOLELE, NLIST, INLIST )
        
        ! Allocate memory for the interpolated upwind values
        ALLOCATE(    TUPWIND(NCOLM) )
        ALLOCATE( TUPWINDOLD(NCOLM) )
        
        ! For the anisotropic limiting scheme we find the upwind values
        ! by interpolation using the subroutine FINPTS or IFINPTS; the upwind
        ! value for each node pair is stored in the matrices TUPWIND AND TUPWINDOLD...
    
        BOUND = .FALSE.
        REFLECT = .FALSE.
        IF((DISOPT.EQ.8).OR.(DISOPT.EQ.9)) THEN
          BOUND = .TRUE.
          REFLECT = .TRUE.
        ENDIF
    
        STORE_ELE = .FALSE.
        RET_STORE_ELE = .FALSE.
        IF(USE_INNER_STOR) THEN
          STORE_ELE=(FIRSTVIS.OR.NEWMES2)
          RET_STORE_ELE =(.NOT.STORE_ELE)
        ENDIF
        
        IF(STORE_ELE) THEN
          IF(.NOT.FIRSTVIS) THEN
            DEALLOCATE(ELEMATPSI)
            DEALLOCATE(ELEMATWEI)
          ENDIF
          
          FIRSTVIS=.FALSE.
          ALLOCATE(ELEMATPSI(NCOLM))
          ALLOCATE(ELEMATWEI(NCOLM*NLOC))
    
          CALL FINPTSSTORE(T,TOLD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                  TUPWIND,TUPWINDOLD,FINDRM,COLM,NCOLM, &
                  XNDGLN,XNONOD, &
                  X,Y,Z, &
                  N,NLX,NLY,NLZ, WEIGHT, &
                  FINDELE,COLELE,NCOLEL, &
                  ELEMATPSI,ELEMATWEI,1, &
                  BOUND, REFLECT)
    
        ELSEIF(RET_STORE_ELE) THEN
          
          CALL GETSTOREELEWEI(T,TOLD,NONODS,NLOC,TOTELE,NDGLNO, &
                    TUPWIND,TUPWINDOLD,FINDRM,COLM,NCOLM,BOUND, &
                    ELEMATPSI,ELEMATWEI)
    
        ELSE
    
          CALL FINPTSSTORE(T,TOLD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
                    TUPWIND,TUPWINDOLD,FINDRM,COLM,NCOLM, &
                    XNDGLN,XNONOD, &
                    X,Y,Z, &
                    N,NLX,NLY,NLZ, WEIGHT, &
                    FINDELE,COLELE,NCOLEL, &
                    DUMMYINT,DUMMYREAL,0, &
                    BOUND, REFLECT)
        ENDIF
    
      ELSE
    
        ! If using the original limiting scheme, the first step
        ! is to estimate the upwind field value from the surrounding nodes
    
        ! Allocate memory for terms needed by GETGXYZ OR ONVDLIM
        ALLOCATE(     TMIN(NONODS) )
        ALLOCATE(     TMAX(NONODS) )
        ALLOCATE(    TMAX2(NONODS) )
        ALLOCATE(    TMIN2(NONODS) )
        ALLOCATE( TOLDMIN2(NONODS) )
        ALLOCATE( TOLDMAX2(NONODS) )
        ALLOCATE(  TOLDMIN(NONODS) )
        ALLOCATE(  TOLDMAX(NONODS) )
        ALLOCATE(     DMIN(NONODS) )
        ALLOCATE(     DMAX(NONODS) )
        ALLOCATE(  DOLDMIN(NONODS) )
        ALLOCATE(  DOLDMAX(NONODS) )
    
        ! For each node, find the largest and smallest value
        ! of T and DENSITY for both the current and previous timestep,
        ! out of the node value and all its surrounding nodes...
        DO NODI=1,NONODS
          TMAX(NODI)=T(NODI)
          TMIN(NODI)=T(NODI)
          TOLDMAX(NODI)=TOLD(NODI)
          TOLDMIN(NODI)=TOLD(NODI)
          DMAX(NODI)=DENSTY(NODI)
          DMIN(NODI)=DENSTY(NODI)
          DOLDMAX(NODI)=DENOLD(NODI)
          DOLDMIN(NODI)=DENOLD(NODI)
          DO COUNT=FINDRM(NODI),FINDRM(NODI+1)-1
            NODJ=COLM(COUNT)
            TMAX(NODI)=MAX(TMAX(NODI),T(NODJ))
            TMIN(NODI)=MIN(TMIN(NODI),T(NODJ))
            TOLDMAX(NODI)=MAX(TOLDMAX(NODI),TOLD(NODJ))
            TOLDMIN(NODI)=MIN(TOLDMIN(NODI),TOLD(NODJ))
            DMAX(NODI)=MAX(DMAX(NODI),DENSTY(NODJ))
            DMIN(NODI)=MIN(DMIN(NODI),DENSTY(NODJ))
            DOLDMAX(NODI)=MAX(DOLDMAX(NODI),DENOLD(NODJ))
            DOLDMIN(NODI)=MIN(DOLDMIN(NODI),DENOLD(NODJ))
          END DO
        END DO
    
        ! Communicate halo values of TMAX,TMIN,TOLDMAX,TOLDMIN
        IF(PARA.EQ.1) THEN
          CALL HALGET(TMAX,   NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(TMIN,   NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(TOLDMAX,NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(TOLDMIN,NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(DMAX,   NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(DMIN,   NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(DOLDMAX,NONODS,NONODS,NNODP,halo_tag)
          CALL HALGET(DOLDMIN,NONODS,NONODS,NNODP,halo_tag)
        ENDIF
    
        IF(GRADBO.AND.(.NOT.(QUAD_INTERP.OR.PERT_APPROACH))) THEN
          ALLOCATE(    DIFXT(NONODS) )
          ALLOCATE(    DIFYT(NONODS) )
          ALLOCATE(    DIFZT(NONODS) )
          ALLOCATE( DIFXTOLD(NONODS) )
          ALLOCATE( DIFYTOLD(NONODS) )
          ALLOCATE( DIFZTOLD(NONODS) )
          ALLOCATE(       ML(NONODS) )
          ! Calculate the gradients of the field being advected for 
          ! this time step and the previous (DIFXT,DIFXTOLD,etc.)...
          CALL GETGXYZ(T,TOLD,DIFXT,DIFYT,DIFZT,DIFXTOLD,DIFYTOLD, &
                DIFZTOLD,ML,TOTELE,X,Y,Z,N,NLX,NLY,NLZ,WEIGHT,NLOC,NGI, &
                D3,DCYL,NDGLNO,NONODS,XNDGLN,XNONOD,DETWEI,RA,NX,NY,NZ)
      
          ! Communicate halo values of DIFXT, DIFYT, DIFZT, etc...
          IF(PARA.EQ.1) THEN
            CALL HALGET(DIFXT,   NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFYT,   NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFZT,   NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFXTOLD,NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFYTOLD,NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFZTOLD,NONODS,NONODS,NNODP,halo_tag)
          ENDIF
        ENDIF
    
        ! Store the maximum and minimum field values (don't know why yet).
          TMAX2(1:NONODS) = TMAX(1:NONODS)
          TMIN2(1:NONODS) = TMIN(1:NONODS)
          TOLDMAX2(1:NONODS) = TOLDMAX(1:NONODS)
          TOLDMIN2(1:NONODS) = TOLDMIN(1:NONODS)
    
      END IF
    

      ! weak boundary condition info is necessary for all meshes
      ALLOCATE(DM1(NONODS))
      ALLOCATE(OLDDM1(NONODS))
      
      DO NOD = 1, NONODS
      
        DM1(NOD) = T(NOD)
        OLDDM1(NOD) = TOLD(NOD)
      
      ENDDO
      
      DO NOD = 1, NOBCT
      
        COUNT = BCT2(NOD)
        
        DM1(COUNT) = BCT1(NOD)
        OLDDM1(COUNT) = BCT1(NOD) ! NEEDS TO CHANGE FOR TIME VARYING BOUNDARY CONDITIONS
      
      ENDDO
        
    !     ======= DEFINE THE SUB-CONTROL VOLUME SHAPE FUNCTIONS, ETC ========
    
      ! SVN(NLOC,SCVNGI)       - the shape function evaluated 
      !                         for each node at each surface gauss point
      ! SVNL[X/Y](NLOC,SCVNGI) - the surface derivatives of the shape 
      !                         function for each node at those same points, and
      ! SVWEIGH(SCVNGI)        - the Gauss weights to use when integrating around 
      !                         the control volume surface
      ! NEILOC(NLOC,SCVNGI)    - neighbour node for a given node/gauss-point pair 
      CALL SHAPESV( ELETYP,  NEILOC,  &  
                    CVNGI,     NLOC, &   
                    SCVNGI, &
                    SVN, &
                    SVNLX,   SVNLY,  &
                    SVWEIGH, &
                    M &
                    )

        IF(QUAD_INTERP.OR.PERT_APPROACH) THEN
! Get quadratic basis functions SVN_QUAD(KLOC,GI) at quadrature points 
!           CALL GET_SVN_QUAD2(SVN_QUAD,NQLOC,SCVNGI,SVN,NLOC, n,ngi) 
           CALL GET_SVN_QUAD(SVN_QUAD,NQLOC,SCVNGI,SVN,NLOC)
          
           ALLOCATE(    DIFXT(NONODS) )
           ALLOCATE(    DIFYT(NONODS) )
           ALLOCATE(    DIFZT(NONODS) )
           ALLOCATE( DIFXTOLD(NONODS) )
           ALLOCATE( DIFYTOLD(NONODS) )
           ALLOCATE( DIFZTOLD(NONODS) )
           ALLOCATE(       ML(NONODS) )
          ! Calculate the gradients of the field being advected for 
          ! this time step and the previous (DIFXT,DIFXTOLD,etc.)...
          CALL GETGXYZ(T,TOLD,DIFXT,DIFYT,DIFZT,DIFXTOLD,DIFYTOLD, &
                DIFZTOLD,ML,TOTELE,X,Y,Z,N,NLX,NLY,NLZ,WEIGHT,NLOC,NGI, &
                D3,DCYL,NDGLNO,NONODS,XNDGLN,XNONOD,DETWEI,RA,NX,NY,NZ)
      
          ! Communicate halo values of DIFXT, DIFYT, DIFZT, etc...
          IF(PARA.EQ.1) THEN
            CALL HALGET(DIFXT,   NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFYT,   NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFZT,   NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFXTOLD,NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFYTOLD,NONODS,NONODS,NNODP,halo_tag)
            CALL HALGET(DIFZTOLD,NONODS,NONODS,NNODP,halo_tag)
          ENDIF
        ENDIF
       IF(PERT_APPROACH) THEN
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
      
      ALLOCATE(ONBOU(NONODS))
      ALLOCATE(SUF_LINK(NCOLM))
      ALLOCATE(SUFS_ONBOU_GE2(TOTELE))
      ALLOCATE(FINELE(TOTELE+1))
      ALLOCATE(COLELE4(TOTELE*5))
      CALL GETFINELE4(TOTELE,NLOC,SNLOC,NDGLNO, COLELE4,NONODS)
      
       ONBOU=.FALSE.
       SUF_LINK=.FALSE.
       DO ELE=1,TOTELE
         ISUF=0
         DO IFACE=1,NFACE
! NB colele4 is ordered in terms of faces.
           ELE2=COLELE4((ELE-1)*5+IFACE)
           IF(ELE2.EQ.0) THEN
             ISUF=ISUF+1
! Calculate the nodes on the other side of the face:
             DO SILOC=1,SNLOC
               ILOC=LOCLIST(IFACE,SILOC)
               NODI=NDGLNO((ELE-1)*NLOC+ILOC)
               ONBOU(NODI)=.TRUE.
               DO SJLOC=1,SNLOC
                 JLOC=LOCLIST(IFACE,SJLOC)
                 NODJ=NDGLNO((ELE-1)*NLOC+JLOC)
                 DO KCOUNT=FINDRM(NODI),FINDRM(NODI+1)-1
                   IF(COLM(KCOUNT).EQ.NODJ) SUF_LINK(KCOUNT)=.TRUE.
                 END DO
               END DO
             END DO
           ENDIF
         END DO
         SUFS_ONBOU_GE2(ELE)=ISUF.GE.2
! ENDOF ELE=1,TOTELE...
        END DO
! Whenever there is more than 2 faces on an element make suf_link=.false.  
       IF(SEVER_LINK_GE2) THEN      
       DO ELE=1,TOTELE
           IF(SUFS_ONBOU_GE2(ELE)) THEN
! Calculate the nodes on the other side of the face:
             DO ILOC=1,NLOC
               NODI=NDGLNO((ELE-1)*NLOC+ILOC)
               DO JLOC=1,NLOC
                 NODJ=NDGLNO((ELE-1)*NLOC+JLOC)
                 DO KCOUNT=FINDRM(NODI),FINDRM(NODI+1)-1
                   IF(COLM(KCOUNT).EQ.NODJ) SUF_LINK(KCOUNT)=.FALSE.
                 END DO
               END DO
             END DO
           ENDIF
! ENDOF ELE=1,TOTELE...
        END DO
        ENDIF
        
        
       ENDIF
                    
      ! Define the gauss points that lie on the surface of the
      ! control volume surrounding a given local node (iloc)
      CALL GAUSSILOC( FINDGPTS, COLGPTS, NCOLGPTS, &
                      NEILOC,   NLOC,    SCVNGI     )
        
        ! get surface shape functions
      CALL SHAPESE( SELETYP,  SNEILOC, &    
                    SNGI,     SNLOC, &
                    SN, &
                    SNLX,   SNLY, &
                    SWEIGH &
                    )
      
      !     =============== DEFINE THETA FOR TIME-STEPPING ===================
      
      ! Define the type of time integration...
      ! Timopt is 0 if DISOPT is even (theta specified); 
      ! Timopt is 1 if DISOPT is odd (non-linear theta scheme) 
      TIMOPT=MOD(DISOPT,2) 
      IF(PERT_APPROACH) TIMOPT=1
      ! Always use the non-linear theta method with the pertabation approach.
      
      OPINTE=0
      GPSIL=0.
      FPSIL=0.
            
      VTHETA=1.0       !THIS SEEMS SUPERFLUOUS (SEE BELOW)    
            
      ! Now we begin the loop over elements to assemble the advection terms
      ! into the matrix (MATRIX) and the RHS...
      DO ELE = 1,TOTELE

        IF(QUAD_INTERP) THEN
! Calculate TQUAD,TOLDQUAD for quadratic interpolation: 
          CALL GET_QUAD_LOC(ELE,TQUAD,TOLDQUAD,NQLOC,T,TOLD, &
            DIFXT,DIFYT,DIFZT, DIFXTOLD,DIFYTOLD,DIFZTOLD, &
            NONODS,XNONOD,NLOC,TOTELE,NDGLNO,XNDGLN,X,Y,Z)
        ENDIF
        
        IF(PERT_APPROACH) THEN
! Calculate DETWEI,RA,NX,NY,NZ for element ELE
          CALL DETNLXR(ELE, X,Y,Z, XNDGLN, TOTELE,NONODS,NLOC,NGI, &
                     N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, D3,DCYL, &
                     NX,NY,NZ) 
        ENDIF
        ! Loop over the nodes of the element
        DO ILOC=1,NLOC 
    
          ! Global node number of the local node
          NODI = NDGLNO((ELE-1)*NLOC+ILOC)
    
          ! Loop over quadrature (gauss) points in ELE neighbouring ILOC
          DO GCOUNT = FINDGPTS(ILOC),FINDGPTS(ILOC+1)-1
            
            ! COLGPTS stores the local Gauss-point number in the ELE
            GI = COLGPTS(GCOUNT)
                    
            ! Get the neighbouring node for node ILOC and Gauss point GI
            JLOC = NEILOC(ILOC,GI)
              
            GOTPNTS = .FALSE.
            ! Calculate the control volume normals at the Gauss pts.
            CALL SCVDETNX( ELE,      GI,      NODI,    &
                          NLOC,     SCVNGI,   TOTELE,  &
                          XNDGLN,   XNONOD,           &
                          CVDETWEI, CVNORMX, CVNORMY, &  
                          CVNORMZ,  SVN,     SVNLX,   &   
                          SVNLY,    SVWEIGH, PNTX,    & 
                          PNTY,     PNTZ,    X,       & 
                          Y,        Z,                &
                          D3,       DCYL,    GOTPNTS )
                
      !     ================ COMPUTE THE FLUX ACROSS SUB-CV FACE ===============
      
                ! Compute the flux NDOTQ through the SUB-CV face at the Gauss point...
            UDGI    = 0.0
            VDGI    = 0.0
            WDGI    = 0.0
            DO KLOC = 1, NLOC
              VNODK = VONDGL((ELE-1)*NLOC+KLOC)
              UDGI  = UDGI + SVN(KLOC,GI)*(NU(VNODK)-UG(VNODK))
              VDGI  = VDGI + SVN(KLOC,GI)*(NV(VNODK)-VG(VNODK))
              WDGI  = WDGI + SVN(KLOC,GI)*(NW(VNODK)-WG(VNODK))
            END DO
            IF(ISUBMOD.NE.0) THEN
            DO KLOC = 1, NSUBNVLOC
              VNODK = (ELE-1)*NLOC+KLOC
              UDGI  = UDGI + SVN(KLOC,GI)*NUSUB(VNODK)
              VDGI  = VDGI + SVN(KLOC,GI)*NVSUB(VNODK)
              WDGI  = WDGI + SVN(KLOC,GI)*NWSUB(VNODK)
            END DO
            ENDIF
            XNODI=XNDGLN((ELE-1)*NLOC+ILOC)
            XNODJ=XNDGLN((ELE-1)*NLOC+JLOC)
            NDOTQ =  CVNORMX(GI)*UDGI + CVNORMY(GI)*VDGI + CVNORMZ(GI)*WDGI
                
            ! Define whether flux is incoming or outgoing, 
            ! depending on direction of flow...
            IF( NDOTQ .GT. 0.0 ) THEN
              INCOME = 0.0  !Outgoing
            ELSE
              INCOME = 1.0  !Incoming
            ENDIF
                
            ! Find its global node number
            NODJ = NDGLNO((ELE-1)*NLOC+JLOC)
    
            IF(GETMAT) THEN
              ! Calculate the the column in the MATRIX corresponding
              ! to this node (see below where matrix is assembled)
              DO COUNT = FINDRM(NODI),FINDRM(NODI+1)-1
                IF( NODJ .EQ. COLM(COUNT) )  JCOUNT = COUNT
              END DO
            ENDIF
    
            ! Compute the distance HDC between the nodes either
            ! side of the CV face (this is needed to compute the
            ! local courant number and the non-linear theta)...
            IF(D3) THEN
              HDC=SQRT((X(XNODI)-X(XNODJ))**2+(Y(XNODI)-Y(XNODJ))**2+(Z(XNODI)-Z(XNODJ))**2)
            ELSE
              HDC=SQRT((X(XNODI)-X(XNODJ))**2+(Y(XNODI)-Y(XNODJ))**2)
            ENDIF
                
      !     ================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
      
            FVT    = INCOME*T(NODJ)    + (1.-INCOME)*T(NODI)
            FVTOLD = INCOME*TOLD(NODJ) + (1.-INCOME)*TOLD(NODI)
            
            FVD    = INCOME*DENSTY(NODJ)    + (1.-INCOME)*DENSTY(NODI)
            FVDOLD = INCOME*DENOLD(NODJ) + (1.-INCOME)*DENOLD(NODI)
                
            ! By default do not use first-order upwinding
            FIRSTORD=.FALSE.
          
            ! No limiting if DISOPT is 6 or 7 
            ! (why not just define limt=femt and skip to assembly?)
            NOLIMI=(INT(DISOPT/2).EQ.3) 
    
            ! Make a guess at the CV face value of advected field 
            ! variable and density (Depends on discetisation option, DISOPT)
            ! First-order upwinding...
            IF((INT(DISOPT/2).EQ.0)) THEN
              FEMT    = FVT
              FEMTOLD = FVTOLD
              FEMD    = FVD
              FEMDOLD = FVDOLD
              FIRSTORD=.TRUE.    
      
            ! Central differencing [Trapezoidal rule (2 OR 3)]
            ELSE IF(INT(DISOPT/2).EQ.1) THEN
              FEMT    = 0.5*(T(NODI)+T(NODJ))   ! need to fix this!
              FEMTOLD = 0.5*(TOLD(NODI)+TOLD(NODJ)) 
              FEMD    = 0.5*(DENSTY(NODI)+DENSTY(NODJ))   ! need to fix this!
              FEMDOLD = 0.5*(DENOLD(NODI)+DENOLD(NODJ)) 
      
            ! Finite element approximation (4 OR 5)(6 or 7)(8 or 9)...
            ELSE
              FEMT    = 0.0
              FEMTOLD = 0.0
              FEMD    = 0.0
              FEMDOLD = 0.0
              DO KLOC = 1, NLOC
                NODK = NDGLNO((ELE-1)*NLOC+KLOC)
                FEMT    = FEMT    + SVN(KLOC,GI)*T(NODK)
                FEMTOLD = FEMTOLD + SVN(KLOC,GI)*TOLD(NODK)
                FEMD    = FEMD    + SVN(KLOC,GI)*DENSTY(NODK)
                FEMDOLD = FEMDOLD + SVN(KLOC,GI)*DENOLD(NODK)
              END DO
      
            ENDIF
              IF(QUAD_INTERP) THEN
                FEMT    = 0.0
                FEMTOLD = 0.0  
                DO KLOC = 1, NQLOC
                  FEMT    = FEMT    + SVN_QUAD(KLOC,GI)*TQUAD(KLOC)
                  FEMTOLD = FEMTOLD + SVN_QUAD(KLOC,GI)*TOLDQUAD(KLOC)
                END DO
              ENDIF
      
      !     ======================== THE LIMITING ============================
      
            ! Compute the limited face value of the advected field
            ! variable for both the old and current time...
            IF (ANISOLIM) THEN   !This anisotropic limiting scheme...
      
              ! Define the courant number for the limiting...
              IF((DISOPT.EQ.8).OR.(DISOPT.EQ.9)) THEN
      
                ! Interface tracking requires the local courant number...
                COURANT=1.7*DT*ABS(NDOTQ)/HDC
                ULTRACIN=ULTRAC
      
                ! -----------------------------------------------------
                ! ------------ EXPERIMENTAL CODE ADDED BY CCP ---------
                ! Attempt at better interface tracking by including 
                ! cross-stream term.  Turn on by setting CROSTR=1
                !
                ! IT IS SAFE TO COMPLETELY IGNORE THIS CODE
                !
                ! Note: Have tried this and it doesn't seem to make any improvement
                !       should check with Chris whether it can be removed...
                IF(CROSTR.EQ.1) THEN
                  GPSIL   =0.0
                  GPSIOLDL=0.0
                  FPSIL   =0.0
                  FPSIOLDL=0.0
                  RN=0.0
                  DO KLOC = 1, NLOC
                    NODK = NDGLNO((ELE-1)*NLOC+KLOC)
                    KCOUNT=0
                    IF(INCOME.GT.0.5) THEN
                      IF(NODK.NE.NODI) THEN
                        DO COUNT=FINDRM(NODK),FINDRM(NODK+1)-1
                          IF(COLM(COUNT).EQ.NODI) KCOUNT=COUNT
                        ENDDO
                      ENDIF
                    ELSE
                      IF(NODK.NE.NODJ) THEN
                        DO COUNT=FINDRM(NODK),FINDRM(NODK+1)-1
                          IF(COLM(COUNT).EQ.NODJ) KCOUNT=COUNT
                        ENDDO
                      ENDIF
                    ENDIF
                    IF(KCOUNT.NE.0) THEN
                      GPSIL    =GPSIL    +SVN(KLOC,GI)*T(NODK)
                      GPSIOLDL =GPSIOLDL +SVN(KLOC,GI)*TOLD(NODK)
                      FPSIL   =FPSIL   +SVN(KLOC,GI)*TUPWIND(KCOUNT)
                      FPSIOLDL=FPSIOLDL+SVN(KLOC,GI)*TUPWINDOLD(KCOUNT)
                      RN=RN+SVN(KLOC,GI)
                    ENDIF
                  END DO
                  GPSIL   =GPSIL   /RN
                  GPSIOLDL=GPSIOLDL/RN
                  FPSIL   =FPSIL   /RN
                  FPSIOLDL=FPSIOLDL/RN
                  OPINTE=1
                END IF
                ! -------------- END OF EXPERIMENTAL CODE ---------------
                ! ------------------------------------------------------- 
                      
              ELSE
      
                ! Anisotropic limiting method does not use courant number
                ! (set the courant number to -ve value)...
                COURANT=-1.0
                ULTRACIN=.FALSE.
      
              ENDIF
      
              IF((DISOPT.EQ.8).or.(DISOPT.EQ.9)) THEN
                ! Normalized variable diagram limiting (New T)
                CALL ANONVDLIM(NONODS,COURANT, &
                              LIMT,FEMT, INCOME, NODI,NODJ, &
                              T,TUPWIND,FIRSTORD,NOLIMI,  &
                              NCOLM,FINDRM,COLM, &
                              OPINTE,GPSIL,FPSIL,ULTRACIN) 
        
                ! Normalized variable diagram limiting (Old T)     
                CALL ANONVDLIM(NONODS,COURANT, &
                              LIMTOLD,FEMTOLD, INCOME, NODI,NODJ, &
                              TOLD,TUPWINDOLD,FIRSTORD,NOLIMI,  &
                              NCOLM,FINDRM,COLM, &
                              OPINTE,GPSIL,FPSIL,ULTRACIN)
              ELSE
    
                CALL BOTHANONVDLIM(NONODS,INCOME, NODI,NODJ, &
                              FIRSTORD,NOLIMI, &
                              NCOLM,FINDRM,COLM, &
                              OPINTE, &
                              LIMT,FEMT,FVT,T,TUPWIND, &
                              LIMTOLD,FEMTOLD,FVTOLD,TOLD,TUPWINDOLD)
              ENDIF                         
      !     ================= ESTIMATE THE FACE VALUE OF THE SUB-CV ===============
      
              LIMD    = FEMD
              LIMDOLD = FEMDOLD
      
                    
            ELSE  ! The original isotropic limiting method for all other cases...
      
              ! Define the donor and receptor nodes depending on the
              ! flux direction...
              NODC=NODI
              NODD=NODJ
              XNODC=XNODI
              XNODD=XNODJ
              IF(INCOME.GT.0.5) THEN
                NODC=NODJ
                NODD=NODI
                XNODC=XNODJ
                XNODD=XNODI
              ENDIF
      
              IF(GRADBO) THEN
                ! Not quite sure what this does.  Probably computes upwind values
                ! at the boundaries...
                CALL GRADBOUND(NODC,NODD,XNODC,XNODD,T,TOLD, &
                              X,Y,Z, &
                              TMIN,TMAX,TOLDMIN,TOLDMAX,  &
                              DIFXT, DIFYT, DIFZT,  &
                              DIFXTOLD, DIFYTOLD, DIFZTOLD, &
                              NONODS,XNONOD, &
                              TMAX2,TMIN2,TOLDMAX2,TOLDMIN2 )
              ENDIF
      
              ! Call the original NVD limiting routine to limit
              ! the face value of the old and new field variable and "density"
              
              
              
        if(PERT_APPROACH) then
! This seems the best COURANT no...  
            IF(COURANT_LIM) THEN          
              COURANT=1.7*dt*abs((X(XNODI)-X(XNODJ))*UDGI+(Y(XNODI)-Y(XNODJ))*VDGI &
              +(Z(XNODI)-Z(XNODJ))*WDGI)/(hdc**2)   
            ELSE 
              COURANT=100.
            ENDIF
            IF((INT(DISOPT/2).EQ.1).AND.(.NOT.QUAD_INTERP)) THEN
! trapizoidal rule...
              XGI    = 0.5*(X(XNODI)+X(XNODJ)) 
              YGI    = 0.5*(Y(XNODI)+Y(XNODJ)) 
              ZGI    = 0.5*(Z(XNODI)+Z(XNODJ)) 
            ELSE
! fem...
              XGI    = 0.0
              YGI    = 0.0
              ZGI    = 0.0
              DO KLOC = 1, NLOC
                NODK = XNDGLN((ELE-1)*NLOC+KLOC)
                XGI    = XGI    + SVN(KLOC,GI)*X(NODK)
                YGI    = YGI    + SVN(KLOC,GI)*Y(NODK)
                ZGI    = ZGI    + SVN(KLOC,GI)*Z(NODK)
              END DO
            ENDIF
            
            
! Calculate gradeint along unseen line interval beyond node D. 
              IF(CENT_HIGH_T) THEN
                XGI    = 0.5*(X(XNODI)+X(XNODJ))
                YGI    = 0.5*(Y(XNODI)+Y(XNODJ))
                ZGI    = 0.5*(Z(XNODI)+Z(XNODJ))
                FEMT   = 0.5*(T(NODI)+T(NODJ))
                FEMTOLD= 0.5*(TOLD(NODI)+TOLD(NODJ))
              ENDIF
              GRAD_ELE_T_X    = 0.0
              GRAD_ELE_T_Y    = 0.0
              GRAD_ELE_T_Z    = 0.0
              GRAD_ELE_TOLD_X    = 0.0
              GRAD_ELE_TOLD_Y    = 0.0
              GRAD_ELE_TOLD_Z    = 0.0
              DO KLOC = 1, NLOC
                NODK = NDGLNO((ELE-1)*NLOC+KLOC)
                GRAD_ELE_T_X    = GRAD_ELE_T_X    + NX(KLOC,1)*T(NODK)
                GRAD_ELE_T_Y    = GRAD_ELE_T_Y    + NY(KLOC,1)*T(NODK)
                GRAD_ELE_T_Z    = GRAD_ELE_T_Z    + NZ(KLOC,1)*T(NODK)
                GRAD_ELE_TOLD_X    = GRAD_ELE_TOLD_X    + NX(KLOC,1)*TOLD(NODK)
                GRAD_ELE_TOLD_Y    = GRAD_ELE_TOLD_Y    + NY(KLOC,1)*TOLD(NODK)
                GRAD_ELE_TOLD_Z    = GRAD_ELE_TOLD_Z    + NZ(KLOC,1)*TOLD(NODK)
              END DO            
                GRAD_PLAN_X_T=DIFXT(NODD)
                GRAD_PLAN_Y_T=DIFYT(NODD)
                GRAD_PLAN_Z_T=DIFZT(NODD)
                GRAD_PLAN_X_TOLD=DIFXTOLD(NODD)
                GRAD_PLAN_Y_TOLD=DIFYTOLD(NODD)
                GRAD_PLAN_Z_TOLD=DIFZTOLD(NODD)
                
                GRAD_LINE_X_T=(T(NODI)-T(NODJ))*(X(XNODI)-X(XNODJ))/(HDC**2)
                GRAD_LINE_Y_T=(T(NODI)-T(NODJ))*(Y(XNODI)-Y(XNODJ))/(HDC**2)
                GRAD_LINE_Z_T=(T(NODI)-T(NODJ))*(Z(XNODI)-Z(XNODJ))/(HDC**2)
                GRAD_LINE_X_TOLD=(TOLD(NODI)-TOLD(NODJ))*(X(XNODI)-X(XNODJ))/(HDC**2)
                GRAD_LINE_Y_TOLD=(TOLD(NODI)-TOLD(NODJ))*(Y(XNODI)-Y(XNODJ))/(HDC**2)
                GRAD_LINE_Z_TOLD=(TOLD(NODI)-TOLD(NODJ))*(Z(XNODI)-Z(XNODJ))/(HDC**2)
                
                ALONG_GRAD_PLAN_X_T=GRAD_PLAN_X_T*ABS(X(XNODI)-X(XNODJ))/HDC
                ALONG_GRAD_PLAN_Y_T=GRAD_PLAN_Y_T*ABS(Y(XNODI)-Y(XNODJ))/HDC
                ALONG_GRAD_PLAN_Z_T=GRAD_PLAN_Z_T*ABS(Z(XNODI)-Z(XNODJ))/HDC
                ALONG_GRAD_PLAN_X_TOLD=GRAD_PLAN_X_TOLD*ABS(X(XNODI)-X(XNODJ))/HDC
                ALONG_GRAD_PLAN_Y_TOLD=GRAD_PLAN_Y_TOLD*ABS(Y(XNODI)-Y(XNODJ))/HDC
                ALONG_GRAD_PLAN_Z_TOLD=GRAD_PLAN_Z_TOLD*ABS(Z(XNODI)-Z(XNODJ))/HDC
! remaining gradient after subtracting the line gradient...                
                REM_GRAD_PLAN_X_T=GRAD_PLAN_X_T -ALONG_GRAD_PLAN_X_T
                REM_GRAD_PLAN_Y_T=GRAD_PLAN_Y_T -ALONG_GRAD_PLAN_Y_T
                REM_GRAD_PLAN_Z_T=GRAD_PLAN_Z_T -ALONG_GRAD_PLAN_Z_T
                REM_GRAD_PLAN_X_TOLD=GRAD_PLAN_X_TOLD -ALONG_GRAD_PLAN_X_TOLD
                REM_GRAD_PLAN_Y_TOLD=GRAD_PLAN_Y_TOLD -ALONG_GRAD_PLAN_Y_TOLD
                REM_GRAD_PLAN_Z_TOLD=GRAD_PLAN_Z_TOLD -ALONG_GRAD_PLAN_Z_TOLD
! Calculate the remainder for the element gradient...
                REM_GRAD_ELE_T_X    =GRAD_ELE_T_X    *(1.0- ABS(X(XNODI)-X(XNODJ))/HDC)
                REM_GRAD_ELE_T_Y    =GRAD_ELE_T_Y    *(1.0- ABS(Y(XNODI)-Y(XNODJ))/HDC)
                REM_GRAD_ELE_T_Z    =GRAD_ELE_T_Z    *(1.0- ABS(Z(XNODI)-Z(XNODJ))/HDC)
                REM_GRAD_ELE_TOLD_X =GRAD_ELE_TOLD_X *(1.0- ABS(X(XNODI)-X(XNODJ))/HDC)
                REM_GRAD_ELE_TOLD_Y =GRAD_ELE_TOLD_Y *(1.0- ABS(Y(XNODI)-Y(XNODJ))/HDC)
                REM_GRAD_ELE_TOLD_Z =GRAD_ELE_TOLD_Z *(1.0- ABS(Z(XNODI)-Z(XNODJ))/HDC)
! calculate gradient along 2nd unseen interval...                
                ALONG_GRAD_PLAN_X_T2=2.0*ALONG_GRAD_PLAN_X_T-GRAD_LINE_X_T
                ALONG_GRAD_PLAN_Y_T2=2.0*ALONG_GRAD_PLAN_Y_T-GRAD_LINE_Y_T
                ALONG_GRAD_PLAN_Z_T2=2.0*ALONG_GRAD_PLAN_Z_T-GRAD_LINE_Z_T
                ALONG_GRAD_PLAN_X_TOLD2=2.0*ALONG_GRAD_PLAN_X_TOLD-GRAD_LINE_X_TOLD
                ALONG_GRAD_PLAN_Y_TOLD2=2.0*ALONG_GRAD_PLAN_Y_TOLD-GRAD_LINE_Y_TOLD
                ALONG_GRAD_PLAN_Z_TOLD2=2.0*ALONG_GRAD_PLAN_Z_TOLD-GRAD_LINE_Z_TOLD
! sum of remaining gradient plus the gradient along 2nd unseen interval... 
              IF(ELE_GRAD) THEN
                GRAD_PLAN_X_T=ALONG_GRAD_PLAN_X_T2+REM_GRAD_ELE_T_X                    
                GRAD_PLAN_Y_T=ALONG_GRAD_PLAN_Y_T2+REM_GRAD_ELE_T_Y             
                GRAD_PLAN_Z_T=ALONG_GRAD_PLAN_Z_T2+REM_GRAD_ELE_T_Z             
                GRAD_PLAN_X_TOLD=ALONG_GRAD_PLAN_X_TOLD2+REM_GRAD_ELE_TOLD_X             
                GRAD_PLAN_Y_TOLD=ALONG_GRAD_PLAN_Y_TOLD2+REM_GRAD_ELE_TOLD_Y    
                GRAD_PLAN_Z_TOLD=ALONG_GRAD_PLAN_Z_TOLD2+REM_GRAD_ELE_TOLD_Z
              ELSE               
                GRAD_PLAN_X_T=ALONG_GRAD_PLAN_X_T2+REM_GRAD_PLAN_X_T                      
                GRAD_PLAN_Y_T=ALONG_GRAD_PLAN_Y_T2+REM_GRAD_PLAN_Y_T              
                GRAD_PLAN_Z_T=ALONG_GRAD_PLAN_Z_T2+REM_GRAD_PLAN_Z_T              
                GRAD_PLAN_X_TOLD=ALONG_GRAD_PLAN_X_TOLD2+REM_GRAD_PLAN_X_TOLD              
                GRAD_PLAN_Y_TOLD=ALONG_GRAD_PLAN_Y_TOLD2+REM_GRAD_PLAN_Y_TOLD      
                GRAD_PLAN_Z_TOLD=ALONG_GRAD_PLAN_Z_TOLD2+REM_GRAD_PLAN_Z_TOLD
              ENDIF
! set gradient to zero if gradeint on unseen section 2 is greater than between nodes C and D...                
            IF(ALONG_GRAD_PLAN_X_T2**2+ALONG_GRAD_PLAN_Y_T2**2+ALONG_GRAD_PLAN_Z_T2**2 &
           .GT.ALONG_GRAD_PLAN_X_T**2+ALONG_GRAD_PLAN_Y_T**2+ALONG_GRAD_PLAN_Z_T**2) THEN               
                GRAD_PLAN_X_T=0.0                     
                GRAD_PLAN_Y_T=0.0             
                GRAD_PLAN_Z_T=0.0       
            ENDIF
            IF(ALONG_GRAD_PLAN_X_TOLD2**2+ALONG_GRAD_PLAN_Y_TOLD2**2+ALONG_GRAD_PLAN_Z_TOLD2**2 &
           .GT.ALONG_GRAD_PLAN_X_TOLD**2+ALONG_GRAD_PLAN_Y_TOLD**2+ALONG_GRAD_PLAN_Z_TOLD**2) THEN               
                GRAD_PLAN_X_TOLD=0.0                     
                GRAD_PLAN_Y_TOLD=0.0             
                GRAD_PLAN_Z_TOLD=0.0       
            ENDIF              
            
            IF(STAN_LIM) THEN
                GRAD_PLAN_X_T=0.0                     
                GRAD_PLAN_Y_T=0.0        
                GRAD_PLAN_Z_T=0.0                             
                GRAD_PLAN_X_TOLD=0.0                     
                GRAD_PLAN_Y_TOLD=0.0             
                GRAD_PLAN_Z_TOLD=0.0  
            ENDIF  
                
              XPT_PLAN=X(NODC)
              YPT_PLAN=Y(NODC)
              ZPT_PLAN=Z(NODC)
! Calculate the perutbed FEM value...
              XDIR=XGI-XPT_PLAN
              YDIR=YGI-YPT_PLAN
              ZDIR=ZGI-ZPT_PLAN
              EST_FEMT=T(NODC) + (XDIR*GRAD_PLAN_X_T + YDIR*GRAD_PLAN_Y_T +ZDIR*GRAD_PLAN_Z_T)
              PERT_FEMT=FEMT-EST_FEMT
              EST_FEMTOLD=TOLD(NODC) + (XDIR*GRAD_PLAN_X_TOLD + YDIR*GRAD_PLAN_Y_TOLD +ZDIR*GRAD_PLAN_Z_TOLD)
              PERT_FEMTOLD=FEMTOLD-EST_FEMTOLD
! Calculate the perturbed value at D node...
              XDIR=X(NODD)-XPT_PLAN
              YDIR=Y(NODD)-YPT_PLAN
              ZDIR=Z(NODD)-ZPT_PLAN
              EST_T_NODD=T(NODC) + (XDIR*GRAD_PLAN_X_T + YDIR*GRAD_PLAN_Y_T +ZDIR*GRAD_PLAN_Z_T)
              PERT_T(NODD)=T(NODD)-EST_T_NODD
              EST_TOLD_NODD=TOLD(NODC) + (XDIR*GRAD_PLAN_X_TOLD + YDIR*GRAD_PLAN_Y_TOLD +ZDIR*GRAD_PLAN_Z_TOLD)
              PERT_TOLD(NODD)=TOLD(NODD)-EST_TOLD_NODD
              PERT_T(NODC)=0.0
              PERT_TOLD(NODC)=0.0
! Calculate TMAX, TMIN at node NODC and NODD...
              DO ILOOP=0,1
                NODCD=NODC*ILOOP + NODD*(1-ILOOP)
                PERT_TMAX(NODCD)=-1.E+20
                PERT_TMIN(NODCD)=1.E+20
                PERT_TOLDMAX(NODCD)=-1.E+20
                PERT_TOLDMIN(NODCD)=1.E+20
                DO COUNT=FINDRM(NODCD),FINDRM(NODCD+1)-1
                  NOD=COLM(COUNT)
!                  IF(NOD.NE.NODCD) THEN
                    XDIR=X(NOD)-XPT_PLAN
                    YDIR=Y(NOD)-YPT_PLAN
                    ZDIR=Z(NOD)-ZPT_PLAN
                    EST_T_NOD=T(NODC) + (XDIR*GRAD_PLAN_X_T + YDIR*GRAD_PLAN_Y_T +ZDIR*GRAD_PLAN_Z_T)
                    PERT_T_MXMIN=T(NOD)-EST_T_NOD
                    PERT_TMAX(NODCD)=MAX(PERT_TMAX(NODCD),PERT_T_MXMIN)
                    PERT_TMIN(NODCD)=MIN(PERT_TMIN(NODCD),PERT_T_MXMIN)
                    EST_TOLD_NOD=TOLD(NODC) + (XDIR*GRAD_PLAN_X_TOLD + YDIR*GRAD_PLAN_Y_TOLD +ZDIR*GRAD_PLAN_Z_TOLD)
                    PERT_TOLD_MXMIN=TOLD(NOD)-EST_TOLD_NOD
                    PERT_TOLDMAX(NODCD)=MAX(PERT_TOLDMAX(NODCD),PERT_TOLD_MXMIN)
                    PERT_TOLDMIN(NODCD)=MIN(PERT_TOLDMIN(NODCD),PERT_TOLD_MXMIN)
!                  ENDIF
                END DO
              END DO
! This limiting method...               
              CALL ONVDLIM_C(COURANT,NONODS, &
                          PERT_LIMT,PERT_FEMT, INCOME, NODI,NODJ, &
                          PERT_T,PERT_TMIN,PERT_TMAX,FIRSTORD,NOLIMI)
      
              CALL ONVDLIM_C(COURANT,NONODS, &
                          PERT_LIMTOLD,PERT_FEMTOLD, INCOME, NODI,NODJ, &
                          PERT_TOLD,PERT_TOLDMIN,PERT_TOLDMAX,FIRSTORD,NOLIMI)
               
              IF(LIMIT_PERT) THEN
                LIMT   =PERT_LIMT    +EST_FEMT
                LIMTOLD=PERT_LIMTOLD +EST_FEMTOLD
              ELSE
                FULLLIM=(PERT_LIMT-0.0)/SIGN(max(1.E-20,abs(PERT_FEMT-0.0)),PERT_FEMT-0.0)
                FULLLIMOLD=(PERT_LIMTOLD-0.0) &
                   /SIGN(max(1.E-20,abs(PERT_FEMTOLD-0.0)),PERT_FEMTOLD-0.0)
                LIMT   =FULLLIM*FEMT    +(1.-FULLLIM)*FVT
                LIMTOLD=FULLLIMOLD*FEMTOLD +(1.-FULLLIMOLD)*FVTOLD
              ENDIF
!              LIMT   =FEMT
!              LIMTOLD=FEMTOLD
! Make sure LIMT lies between FEMT and FVT
              LIMT=MIN(MAX(FEMT,FVT),LIMT)
              LIMT=MAX(MIN(FEMT,FVT),LIMT)
              LIMTOLD=MIN(MAX(FEMTOLD,FVTOLD),LIMTOLD)
              LIMTOLD=MAX(MIN(FEMTOLD,FVTOLD),LIMTOLD)
                          
! ff:
!             if((abs(z(nodc)-500.0).lt.10.0).and.(abs(z(nodd)-500.0).gt.10.0)) then
             IF(ONBOU(NODC)) THEN
               DO KCOUNT=FINDRM(NODC),FINDRM(NODC+1)-1
                 IF(COLM(KCOUNT).EQ.NODD) KKCOUNT=KCOUNT
               END DO
!               IF(.NOT.SUF_LINK(KKCOUNT)) THEN
               IF((.NOT.SUF_LINK(KKCOUNT)).OR.(SUFS_ONBOU_GE2(ELE))) THEN
!               IF((.NOT.ONBOU(NODD)).OR.(SUF_LINK(KKCOUNT)) ) THEN            
!             if(.not.(ONBOU(NODC)).OR.(.NOT.ONBOU(NODD)))) THEN
!               print *,'nodc,nodd:',nodc,nodd
!               print *,'x(nodc),y(nodc),z(nodc):',x(nodc),y(nodc),z(nodc)
!               print *,'x(nodd),y(nodd),z(nodd):',x(nodd),y(nodd),z(nodd)
!               print *,'ONBOU(NODC),ONBOU(NODD):',ONBOU(NODC),ONBOU(NODD)
!             endif
! ff--:
!             if((abs(z(nodc)-500.0).lt.10.0).or.(abs(z(nodd)-500.0).lt.10.0)) then
! ff---:
!             if(((abs(z(nodc)-500.0).lt.10.0).and.(abs(z(nodd)-500.0).gt.10.0)) &
!            .or.((abs(z(nodd)-500.0).lt.10.0).and.(abs(z(nodc)-500.0).gt.10.0))) then
                 limt=femt
                 limtold=femtold
               endif
             endif
             
            
        else
              
              CALL ONVDLIM(NONODS, &
                          LIMT,FEMT, INCOME, NODI,NODJ, &
                          T,TMIN,TMAX,FIRSTORD,NOLIMI)
      
              CALL ONVDLIM(NONODS, &
                          LIMTOLD,FEMTOLD, INCOME, NODI,NODJ, &
                          TOLD,TOLDMIN,TOLDMAX,FIRSTORD,NOLIMI)
         endif
      
              CALL ONVDLIM(NONODS, &
                          LIMD,FEMD, INCOME, NODI,NODJ, &
                          DENSTY,DMIN,DMAX,FIRSTORD,NOLIMI)
      
              CALL ONVDLIM(NONODS, &
                          LIMDOLD,FEMDOLD, INCOME, NODI,NODJ, &
                          DENOLD,DOLDMIN,DOLDMAX,FIRSTORD,NOLIMI)
      
            ENDIF ! ENDOF IF((ALOLIM.EQ.1).AND.(NLOC.EQ.4)) THEN ELSE ...
      
            ! Define theta
            IF(TIMOPT.EQ.0) THEN !Specified
              FTHETA=THETA
            ELSE                 !Non-linear
              GF=TOLFUN( DT* NDOTQ*(LIMD*LIMT-LIMDOLD*LIMTOLD) )   ! IS INCLUDING LIMD APPROPRIATE?
              PINVTH=HDC*(T(NODI)-TOLD(NODI))/GF
              QINVTH=HDC*(T(NODJ)-TOLD(NODJ))/GF
              FTHETA=MAX(0.5, 1.-0.5*MIN(ABS(PINVTH),ABS(QINVTH))) 
            ENDIF
      
      !     ====================== MATRIX AND RHS ASSEMBLY ===================
      
            IF(GETMAT) THEN
              ! - Calculate the integration of the limited, 
              ! high-order flux over a face  
              ! Conservative discretisation...
              ! The matrix (PIVOT ON LOW ORDER SOLN)...
              MATRIX(JCOUNT) =  MATRIX(JCOUNT) &
                              + ptheta*CVDETWEI(GI)*NDOTQ*INCOME*LIMD
                  
              MATRIX(CENTRM(NODI)) =  MATRIX(CENTRM(NODI)) &
                                    + ptheta*CVDETWEI(GI)*NDOTQ*(1.-INCOME)*LIMD
                  
              ! BETA=0 for Non-conservative discretisation 
              ! (BETA=1 for conservative disc)...
              MATRIX(CENTRM(NODI)) =  MATRIX(CENTRM(NODI))  &
                                    - btheta*(1.-BETA)*CVDETWEI(GI)*NDOTQ*LIMD
            ENDIF
      
            TMID=T(NODI)
            TOLDMID=TOLD(NODI)
                
            ! - Put results into the RHS vector
            RHS(NODI) =  RHS(NODI) &
                      + ptheta*NDOTQ*CVDETWEI(GI)*FVT*LIMD &
                      - NDOTQ*CVDETWEI(GI) &
                        *(FTHETA*LIMD*LIMT+(1.-FTHETA)*LIMDOLD*LIMTOLD)
    ! 
    !         ! Subtract out 1st order...
             RHS(NODI) =  RHS(NODI) &
                        - btheta*(1.-BETA)*NDOTQ*CVDETWEI(GI)*T(NODI)*LIMD
    !   
    !         ! HIGH ORDER CONTRIBUTION...
             RHS(NODI) =  RHS(NODI) &
                        + (1.-BETA)*NDOTQ*CVDETWEI(GI) &
                          *( fTHETA*TMID*LIMD + (1.-fTHETA)*TOLDMID*LIMDOLD )
    !         ! 
    
      
          END DO              !Gauss-point loop
        END DO                 !Node loop
      END DO                    !Element loop

        
      ! Free up all dynamically allocated memory before exiting...
      DEALLOCATE( SVN,SVNLX,SVNLY,SVWEIGH )
      DEALLOCATE( CVDETWEI,CVNORMX,CVNORMY,CVNORMZ )
      DEALLOCATE( NEILOC,COLGPTS,FINDGPTS )
      IF (ANISOLIM) THEN
        DEALLOCATE( TUPWIND,TUPWINDOLD )
        DEALLOCATE( FINDELE, COLELE )
      ELSE
        DEALLOCATE( TMIN,TMAX,TOLDMIN,TOLDMAX,DMIN,DMAX,DOLDMIN,DOLDMAX )
        DEALLOCATE( TMIN2,TMAX2,TOLDMIN2,TOLDMAX2 )
!        DEALLOCATE( DIFXT,DIFYT,DIFZT,DIFXTOLD,DIFYTOLD,DIFZTOLD,ML )
        DEALLOCATE( DETWEI,RA,NX,NY,NZ )
      ENDIF
    
      RETURN
    END SUBROUTINE CONSTRUCT_ADVECTION_DIFFUSION_CV
    

!  
!
        SUBROUTINE GET_SVN_QUAD(SVN_QUAD,NQLOC,SCVNGI,SVN,NLOC)
! Get quadratic basis functions SVN_QUAD(KLOC,GI) at quadrature points 
        IMPLICIT NONE
        INTEGER NLOC,SCVNGI,NQLOC
        REAL SVN(NLOC,SCVNGI),SVN_QUAD(NQLOC,SCVNGI)
        REAL L1, L2, L3, L4
! Local variables...
        INTEGER GI
!
        DO GI=1,SCVNGI
! The local coords are calculated from the linear shape functions SVN.
         L1=SVN(1,GI)
         L2=SVN(2,GI)
         L3=SVN(3,GI)
         L4=SVN(4,GI)
         if(.false.) then
! original ordering...
           SVN_QUAD(1,GI)=(2.*L1-1.)*L1
           SVN_QUAD(3,GI)=(2.*L2-1.)*L2
           SVN_QUAD(5,GI)=(2.*L3-1.)*L3
           SVN_QUAD(10,GI)=(2.*L4-1.)*L4
!
           SVN_QUAD(2,GI)=4.*L1*L2
           SVN_QUAD(6,GI)=4.*L1*L3
           SVN_QUAD(7,GI)=4.*L1*L4
!
           SVN_QUAD(4,GI) =4.*L2*L3
           SVN_QUAD(9,GI) =4.*L3*L4
           SVN_QUAD(8,GI)=4.*L2*L4
         else
! modified ordering consistent with ILINK2...
           SVN_QUAD(1,GI)=(2.*L1-1.)*L1
           SVN_QUAD(3,GI)=(2.*L2-1.)*L2
           SVN_QUAD(6,GI)=(2.*L3-1.)*L3
           SVN_QUAD(10,GI)=(2.*L4-1.)*L4
!
           SVN_QUAD(2,GI)=4.*L1*L2
           SVN_QUAD(4,GI)=4.*L1*L3
           SVN_QUAD(7,GI)=4.*L1*L4
!
           SVN_QUAD(5,GI) =4.*L2*L3
           SVN_QUAD(9,GI) =4.*L3*L4
           SVN_QUAD(8,GI)=4.*L2*L4
         
         endif
          
        END DO
        RETURN
        END SUBROUTINE GET_SVN_QUAD
!
! 
    
    SUBROUTINE IFINPTS(PSI,PSIOLD,NONODS,NLOC,NGI,TOTELE,NDGLNO, &
        MATPSI,MATPSIOLD,FINDRM,COLM,NCOLM, &
        X,Y,Z, &
        N,NLX,NLY,NLZ, WEIGHT, &
        FINDELE, COLELE, NCOLEL ) 
      !======================================================================
      !  This subroutine takes the fields PSI and PSIOLD and computes, for
      !  each adjacent node-pair, a value of that field at a prescribed
      !  interpolation point.  These values are stored in the matrices,
      !  MATPSI and MATPSIOLD.  (THIS VERSION IS FOR INTERFACE TRACKING)
      !
      !  Currently, the interpolation point is prescribed to be...
      !
      !
      !  IMPORTANT INPUTS:
      !  ----------------
      !  PSI,PSIOLD - Field values used in interpolation
      !
      !
      !  IMPORTANT OUTPUTS:
      !  -----------------
      !  MATPSI,MATPSIOLD - Matrix of interpolation points
      !
      !
      !
      !  Called from CONSTRUCT_ADVECTION_DIFFUSION_CV for the interface tracking method; i.e.
      !  NDISOT(IT) ending in 8/9.
      !
      !  For more details see (enter manual/wiki page here...)
      !     
      !  Description                                   Programmer      Date
      !  ==================================================================
      !  Original version..................................CCP     Unknown!
      !  Rewritten as free-form fortran 90 (including 
      !  dynamic memoryallocation).........................GSC   2006-08-11
      !
      !***********************************************************************  
      
      !Inputs/Outputs
      INTEGER :: NONODS,NLOC,NGI,TOTELE
      INTEGER :: NDGLNO(TOTELE*NLOC)
      INTEGER :: NCOLM, NCOLEL
      REAL    :: PSI(NONODS),PSIOLD(NONODS)
      REAL    :: MATPSI(NCOLM),MATPSIOLD(NCOLM)
      INTEGER :: FINDRM(NONODS+1),COLM(NCOLM)
      REAL    :: X(NONODS),Y(NONODS),Z(NONODS)
      REAL    :: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL    :: WEIGHT(NGI)
      INTEGER :: FINDELE(NONODS+1),COLELE(NCOLEL)
    
      ! Local memory...
      REAL, ALLOCATABLE :: NORMX(:),NORMY(:),NORMZ(:)
      REAL, ALLOCATABLE :: MLUM(:) !NONODS
      REAL, ALLOCATABLE :: DETWEI(:),RA(:) !NGI
      REAL, ALLOCATABLE :: NX(:,:),NY(:,:),NZ(:,:) !NLOC,NGI
    
      !     Local variables...
      INTEGER NOD,COUNT,NODI,NODJ
      INTEGER ELE,ILOC,GI
      REAL LENG,VOLUME,INVH
    
      ! Dynamically allocate memory...
      ALLOCATE( NORMX(NONODS) )
      ALLOCATE( NORMY(NONODS) )
      ALLOCATE( NORMZ(NONODS) )
      ALLOCATE(  MLUM(NONODS) )
    
      ALLOCATE( DETWEI(NGI) )
      ALLOCATE(     RA(NGI) )
    
      ALLOCATE( NX(NLOC,NGI) )
      ALLOCATE( NY(NLOC,NGI) )
      ALLOCATE( NZ(NLOC,NGI) )
    
      ! calculate normals...************************    
      NORMX=0.
      NORMY=0.
      NORMZ=0.
      MLUM=0.
      DO ELE=1,TOTELE
        ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
        CALL DETNLXR(ELE, X,Y,Z, NDGLNO, TOTELE,NONODS,NLOC,NGI, &
              N,NLX,NLY,NLZ, WEIGHT, DETWEI,RA,VOLUME, .TRUE.,.FALSE., & 
              NX,NY,NZ) 
    
        DO ILOC=1,NLOC
            NODI=NDGLNO((ELE-1)*NLOC+ILOC)
            DO GI=1,NGI
              NORMX(NODI)=NORMX(NODI)+NX(ILOC,GI)*DETWEI(GI)
              NORMY(NODI)=NORMY(NODI)+NY(ILOC,GI)*DETWEI(GI)
              NORMZ(NODI)=NORMZ(NODI)+NZ(ILOC,GI)*DETWEI(GI)
              MLUM(NODI) =MLUM(NODI) +N(ILOC,GI) *DETWEI(GI)
            END DO
        END DO
    
      END DO
    
      ! Renormalise
      DO NODI=1,NONODS
        INVH=(ABS(NORMX(NODI))+ABS(NORMY(NODI))+ABS(NORMZ(NODI)))/MLUM(NODI)
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
      ! calculate normals...************************    
    
      DO NOD=1,NONODS
    
        DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
            NODJ=COLM(COUNT)
            MATPSI(COUNT)=0.
    
            IF(NOD.NE.NODJ) THEN
              CALL IMATPTS(MATPSI(COUNT),MATPSIOLD(COUNT),NOD, &
                    PSI,PSIOLD,NONODS, &
                    PSI(NODJ),PSIOLD(NODJ), &
                    NLOC,TOTELE,NDGLNO, &
                    FINDRM,COLM,NCOLM, &
                    X(NOD),Y(NOD),Z(NOD), &
                    X(NODJ),Y(NODJ),Z(NODJ), &
                    NORMX(NOD),NORMY(NOD),NORMZ(NOD), &
                    X,Y,Z, &
                    FINDELE,COLELE,NCOLEL) 
            ENDIF
    
        END DO
    
      END DO
    
      ! Free the dynamically allocated memory
      DEALLOCATE(NORMX,NORMY,NORMZ,MLUM)
      DEALLOCATE(DETWEI,RA,NX,NY,NZ)
    
      RETURN
    END SUBROUTINE IFINPTS
    
    SUBROUTINE FINPTS(PSI,PSIOLD,NONODS,NLOC,TOTELE,NDGLNO, &
        MATPSI,MATPSIOLD,FINDRM,COLM,NCOLM, &
        X,Y,Z, &
        FINDELE, COLELE, NCOLEL )
      !======================================================================
      !  This subroutine takes the fields PSI and PSIOLD and computes, for
      !  each adjacent node-pair, a value of that field at a prescribed
      !  interpolation point.  These values are stored in the matrices,
      !  MATPSI and MATPSIOLD.
      !
      !  Currently, the interpolation point is prescribed to be...
      !
      !
      !  IMPORTANT INPUTS:
      !  -----------------
      !  PSI,PSIOLD - Field values used in interpolation
      !
      !
      !  IMPORTANT OUTPUTS:
      !  ------------------
      !  MATPSI,MATPSIOLD - Matrix of interpolation points
      !
      !
      !  Called from CONSTRUCT_ADVECTION_DIFFUSION_CV for use with the anisotropic limiter; i.e.
      !  if the parameter ALOLIM=1 in CONSTRUCT_ADVECTION_DIFFUSION_CV...
      !
      !  For more details see (enter manual/wiki page here...)
      !     
      !  Description                                   Programmer      Date
      !  ==================================================================
      !  Original version..................................CCP     Unknown!
      !  Rewritten as free-form fortran 90 (including 
      !  dynamic memoryallocation).........................GSC   2006-08-11
      !
      !======================================================================
      
      !Inputs/Outputs
      INTEGER :: NONODS,NLOC,TOTELE,NCOLM
      INTEGER :: NDGLNO(TOTELE*NLOC)
      REAL    :: PSI(NONODS),PSIOLD(NONODS)
      REAL    :: MATPSI(NCOLM),MATPSIOLD(NCOLM)
      INTEGER :: FINDRM(NONODS+1),COLM(NCOLM)
      REAL    :: X(NONODS),Y(NONODS),Z(NONODS)
      
      INTEGER :: NCOLEL
      INTEGER :: FINDELE(NONODS+1),COLELE(NCOLEL)
    
      !     Local variables...
      INTEGER NOD,COUNT,NODJ
    
      ! Interpolate to find the value of...
      DO NOD=1,NONODS     
        DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
            NODJ=COLM(COUNT)
            MATPSI(COUNT)=0.
            IF(NOD.NE.NODJ) THEN
              CALL MATPTS(MATPSI(COUNT),MATPSIOLD(COUNT),NOD, &
                    PSI,PSIOLD,NONODS, &
                    NLOC,TOTELE,NDGLNO, &
                    FINDRM,COLM,NCOLM, &
                    X(NOD),Y(NOD),Z(NOD), &
                    X(NODJ),Y(NODJ),Z(NODJ), &
                    X,Y,Z, &
                    FINDELE,COLELE,NCOLEL) 
            ENDIF
        END DO
      END DO
      RETURN
      
      
    END SUBROUTINE FINPTS

    SUBROUTINE ASSEMBLE_MOMENTUM_ADVECTION_CV( &
        NONODS,VNONOD,XNONOD,TOTELE,NLOC,NGI, &
        SNLOC, STOTEL, &
        NDGLNO,XNDGLN,VONDGL, &
        SNDGLN, &
        FINDRM,COLM,NCOLM,CENTRM, & 
        NBIGM,BIGM, VECX,VECY,VECZ, & 
        X,Y,Z, &
        N,NLX,NLY,NLZ, WEIGHT, &
        NU,NV,NW, UG,VG,WG, &
        U,V,W, UOLD,VOLD,WOLD, &
        DENSTY,DENOLD, &
        DISOPT,DT,THETA,BETA, D3,DCYL, &
        PARA,halo_tag,NNODP &
        , NOBCU,BCU1,BCU2 &
        , NOBCV,BCV1,BCV2 &
        , NOBCW,BCW1,BCW2, option_path)
      !  =====================================================================
      !              !!! --- THIS ROUTINE MAY BE OBSOLETE --- !!!
      !     In this subroutine the advection terms in the momentum
      !     equation are calculated for each momentum component using 
      !     a control volume formulation. The terms for each momentum component 
      !     are calculated sequentially.  They are assembled in MATRIX and RHS.  
      !
      !     They are then added to the rest of the momentum equation terms, 
      !     which were previously calculated somewhere (NAVSTO?).
      !
      !     This routine is only called if the GEM option NDISOP is
      !     non zero; the last digit of NDISOP is passed in to this 
      !     routine as DISOPT. Note that at the time of writing this routine
      !     IS OBSOLETE as far as I can tell (GSC).
      !
      !     The assembly of MATRIX and RHS is done in subroutine CONSTRUCT_ADVECTION_DIFFUSION_CV
      !     which is included in this file, below.  Please see this routine for
      !     more detailed comments on the advection routine.
      !
      !     This routine is called from NAVSTO (but is probably obsolete).
      !
      !     For more details see (enter manual/wiki page here...)
      !     
      !     Description                                   Programmer      Date
      !     ==================================================================
      !     Original version..................................CCP     Unknown!
      !     General clean up, improving efficiency and
      !     adding comments...................................GSC   2006-08-10
      !     Rewritten as free-form fortran 90.................GSC   2006-08-11
      !
      !***********************************************************************
    
      ! Inputs/Outputs
      INTEGER :: NONODS,VNONOD,XNONOD,TOTELE,NLOC,NGI
      INTEGER, INTENT(IN) :: SNLOC, STOTEL
      INTEGER :: NDGLNO(TOTELE*NLOC),XNDGLN(TOTELE*NLOC),VONDGL(TOTELE*NLOC)
      INTEGER, INTENT(IN) :: SNDGLN(STOTEL*SNLOC)
      INTEGER :: FINDRM(NONODS+1),NCOLM,COLM(NCOLM),CENTRM(NONODS)
      INTEGER :: NBIGM
      REAL    :: BIGM(NCOLM)
      REAL    :: VECX(NONODS),VECY(NONODS),VECZ(NONODS)
      REAL    :: X(XNONOD),Y(XNONOD),Z(XNONOD)
      REAL    :: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL    :: WEIGHT(NGI) 
      REAL    :: NU(VNONOD),NV(VNONOD),NW(VNONOD)
      REAL    :: UG(VNONOD),VG(VNONOD),WG(VNONOD)
      REAL    :: U(NONODS),V(NONODS),W(NONODS)
      REAL    :: UOLD(NONODS),VOLD(NONODS),WOLD(NONODS)
      REAL    :: DENSTY(NONODS),DENOLD(NONODS)
      INTEGER :: DISOPT
      REAL    :: DT,THETA,BETA 
      LOGICAL :: D3,DCYL
      INTEGER :: PARA,NNODP !Are these really inputs/outputs?
      INTEGER :: halo_tag
    
      ! Element type informtation 
      INTEGER :: SCVNGI
      INTEGER :: ELETYP,OPTELM
      INTEGER :: CVNGI,CVNLOC,SELETYP          !Dummy integers output from SETELM
      LOGICAL :: CVD3,CVDCYL                                 !Dummy logicals output from SETELM
      INTEGER :: SNLOC2, SNGI
      INTEGER :: NOD,COUNT  !Loop indices 
      LOGICAL :: REDQUAD
      PARAMETER(REDQUAD=.FALSE.)
    
      INTEGER NOBCU,NOBCV,NOBCW
      INTEGER BCU2(NOBCU),BCV2(NOBCV),BCW2(NOBCW)
      REAL    BCU1(NOBCU),BCV1(NOBCV),BCW1(NOBCW)
        
      REAL, ALLOCATABLE :: MATRIX(:)
      REAL, ALLOCATABLE :: RHS(:)
      character(len=*) :: option_path
          
      ! ============================= FIND SCVNGI ===============================
      ! Define the element type
      IF(D3) THEN                      !  element type                | optelm
        IF(NLOC.EQ.4)  OPTELM=5       ! -------------------------------------
        IF(NLOC.EQ.8)  OPTELM=7       ! - linear        triangles        1
        IF(NLOC.EQ.10) OPTELM=6       ! - quadratic     triangles        2
        IF(NLOC.EQ.27) OPTELM=8       ! - bi-linear     quadrilaterals   3
      ELSE                             ! - bi-quadratic  quadrilaterals   4
        IF(NLOC.EQ.3) OPTELM=1        ! - linear        tetrahedra       5
        IF(NLOC.EQ.4) OPTELM=3        ! - quadratic     tetrahedra       6
        IF(NLOC.EQ.6) OPTELM=2        ! - tri-linear    hexahedra        7
        IF(NLOC.EQ.9) OPTELM=4        ! - tri-quadratic hexahedra        8
      ENDIF
    
      ! Use the SETELM routine to define the number of gauss points 
      ! on the control volume surface (SCVNGI), based on the element type
      ! SETELM takes the inputs ELETYP and REDQUAD, all other
      ! terms are outputs and are not used
      CALL SETELM( ELETYP, CVNGI, CVNLOC, &
      ! - INTEGERS
          OPTELM, SELETYP, &
          SNGI,   SNLOC2,   SCVNGI, &
      ! - LOGICALS
          CVD3,     CVDCYL,    REDQUAD   ) 
    
      ! make sure we are using the correct shape function...
      IF(CVNLOC.NE.NLOC) STOP 'CVNLOC.NE.NLOC IN ASSEMBLE_MOMENTUM_ADVECTION_CV'
      IF(SNLOC2.NE.SNLOC) STOP 'SNLOC2.NE.SNLOC IN ASSEMBLE_MOMENTUM_ADVECTION_CV'
    
      ! Allocate memory to the big matrix...
      IF(NBIGM.EQ.NCOLM) THEN
        !     ERROR("ERROR: DOES NOT WORK FOR NBIGM=NCOLM")
        STOP 'ERROR (ASSEMBLE_MOMENTUM_ADVECTION_CV): DOES NOT WORK FOR NBIGM=NCOLM' !8743
      ENDIF
      ALLOCATE( MATRIX(NCOLM) ) ! Do we need to zero this, or any of these for that matter?
      ALLOCATE(   RHS(NONODS) )
      
      ! ======================== MOMENTUM IN X-DIRECTION =========================
      ! Call the routine to assemble the matrix of advection terms for U
      ! Returns MATRIX and RHS...
      MATRIX=0.
      RHS=0.

      CALL CONSTRUCT_ADVECTION_DIFFUSION_CV(NONODS,VNONOD,XNONOD,TOTELE, &
          NLOC,NGI,SCVNGI,CVNGI,ELETYP,SELETYP, &
          SNLOC, SNGI, STOTEL, &
          NDGLNO,XNDGLN,VONDGL, &
          SNDGLN, &
          FINDRM,COLM,NCOLM,CENTRM, & 
          MATRIX,RHS, &
          X,Y,Z, &
          N,NLX,NLY,NLZ, WEIGHT, &
          NU,NV,NW, UG,VG,WG, U,UOLD,DENSTY,DENOLD, &
          DISOPT,DT,THETA,BETA, D3,DCYL, &
          PARA,halo_tag,NNODP &
          , .FALSE. &
          , NOBCU, BCU1, BCU2 &
            , 1, 0, .TRUE. &
              ,0,NU,NV,NW, option_path &
          )
    
      ! Combine advection terms computed in CONSTRUCT_ADVECTION_DIFFUSION_CV with the
      ! diffusion terms already assembled in BIGM, and 
      ! subtract A*TOLD from VECX
      DO NOD=1,NONODS
        VECX(NOD)=VECX(NOD)+RHS(NOD) 
      END DO
      DO NOD=1,NONODS
        DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
            VECX(NOD)=VECX(NOD)-MATRIX(COUNT)*UOLD(COLM(COUNT))
        END DO
      END DO
      DO COUNT=1,NCOLM
        BIGM(COUNT)=BIGM(COUNT)+DT*MATRIX(COUNT)
      END DO
    
      ! ======================== MOMENTUM IN Y-DIRECTION =========================
      ! Call the routine to assemble the matrix of advection terms for V
      ! Returns MATRIX and RHS...
      MATRIX=0.
      RHS=0.

      CALL CONSTRUCT_ADVECTION_DIFFUSION_CV(NONODS,VNONOD,XNONOD,TOTELE, &
          NLOC,NGI,SCVNGI,CVNGI,ELETYP,SELETYP, &
          SNLOC, SNGI, STOTEL, &
          NDGLNO,XNDGLN,VONDGL, &
          SNDGLN, &
          FINDRM,COLM,NCOLM,CENTRM, & 
          MATRIX,RHS, &
          X,Y,Z, &
          N,NLX,NLY,NLZ, WEIGHT, &
          NU,NV,NW, UG,VG,WG, V,VOLD,DENSTY,DENOLD, &
          DISOPT,DT,THETA,BETA, D3,DCYL, &
          PARA,halo_tag,NNODP &
          , .FALSE. &
          , NOBCV, BCV1, BCV2 &
            , 1, 0, .TRUE. &
              ,0,NU,NV,NW, option_path &
          )
    
      ! Combine advection terms computed in CONSTRUCT_ADVECTION_DIFFUSION_CV with the
      ! diffusion terms already assembled in BIGM, and 
      ! subtract A*TOLD from VECY
      DO NOD=1,NONODS
        VECY(NOD)=VECY(NOD)+RHS(NOD) 
      END DO
      DO NOD=1,NONODS
        DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
            VECY(NOD)=VECY(NOD)-MATRIX(COUNT)*VOLD(COLM(COUNT))
        END DO
      END DO
      IF(D3) THEN
        DO COUNT=1,NCOLM
            BIGM(COUNT+4*NCOLM)=BIGM(COUNT+4*NCOLM)+DT*MATRIX(COUNT)
        END DO
      ELSE
        DO COUNT=1,NCOLM
            BIGM(COUNT+3*NCOLM)=BIGM(COUNT+3*NCOLM)+DT*MATRIX(COUNT)
        END DO
      ENDIF
    
      ! ======================== MOMENTUM IN Z-DIRECTION =========================
      IF(D3) THEN
    
        ! Call the routine to assemble the matrix of advection terms for W
        ! Returns MATRIX and RHS...
        MATRIX=0.
        RHS=0.

        CALL CONSTRUCT_ADVECTION_DIFFUSION_CV(NONODS,VNONOD,XNONOD,TOTELE, &
              NLOC,NGI,SCVNGI,CVNGI,ELETYP,SELETYP, &
              SNLOC, SNGI, STOTEL, &
              NDGLNO,XNDGLN,VONDGL, &
              SNDGLN, &
              FINDRM,COLM,NCOLM,CENTRM, & 
              MATRIX,RHS, &
              X,Y,Z, &
              N,NLX,NLY,NLZ, WEIGHT, &
              NU,NV,NW, UG,VG,WG, W,WOLD,DENSTY,DENOLD, &
              DISOPT,DT,THETA,BETA, D3,DCYL, &
              PARA,halo_tag,NNODP &
              , .FALSE. &
              , NOBCW, BCU1, BCU2 &
              , 1, 0, .TRUE. &
              ,0,NU,NV,NW, option_path &
              )
        
        ! Combine advection terms computed in CONSTRUCT_ADVECTION_DIFFUSION_CV with the
        ! diffusion terms already assembled in BIGM, and 
        ! subtract A*TOLD from VECZ
        DO NOD=1,NONODS
            VECZ(NOD)=VECZ(NOD)+RHS(NOD) 
        END DO
        DO NOD=1,NONODS
            DO COUNT=FINDRM(NOD),FINDRM(NOD+1)-1
              VECZ(NOD)=VECZ(NOD)-MATRIX(COUNT)*WOLD(COLM(COUNT))
            END DO
        END DO
        DO COUNT=1,NCOLM
            BIGM(COUNT+8*NCOLM)=BIGM(COUNT+8*NCOLM)+DT*MATRIX(COUNT)
        END DO
        
      ENDIF
    
      ! Free the memory allocated to the matrix and right-hand side...
      DEALLOCATE(MATRIX,RHS)
    
      RETURN
    END SUBROUTINE ASSEMBLE_MOMENTUM_ADVECTION_CV
    

end module legacy_advection_diffusion_CV
