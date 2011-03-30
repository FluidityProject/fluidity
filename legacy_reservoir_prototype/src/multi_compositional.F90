
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

module Compositional_Terms
  use fldebug


contains


  SUBROUTINE CALC_COMP_ABSORB( ICOMP, NCOMP, DT, ALPHA_BETA, &
       NPHASE, CV_NONODS, &
       KCOMP_SIGMOID, K_COMP2, &
       DENOLD, SATURAOLD, &
       VOLFRA_PORE, COMP_ABSORB, &
       TOTELE, CV_NLOC, CV_NDGLN )

    ! Calculate compositional model linkage between the phase expressed in COMP_ABSORB. 
    ! Use values from the previous time step so its easier to converge. 
    ! ALPHA_BETA is the scaling coeff. of the compositional model e.g. =1.0

    IMPLICIT NONE
    INTEGER, intent( in ) :: ICOMP, NCOMP, NPHASE, CV_NONODS, TOTELE, CV_NLOC
    REAL, intent( in ) :: DT, ALPHA_BETA
    LOGICAL, intent( in ) :: KCOMP_SIGMOID
    REAL, DIMENSION( NCOMP, NPHASE, NPHASE ), intent( in ) :: K_COMP2
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: DENOLD, SATURAOLD
    REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( CV_NONODS, NPHASE, NPHASE  ), intent( inout ) :: COMP_ABSORB
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN

    ! Local Variables
    INTEGER :: IPHASE, JPHASE, ELE, CV_ILOC, CV_NOD, JCOMP
    REAL, DIMENSION( : ), allocatable :: ALPHA, SUM_NOD, VOLFRA_PORE_NOD
    REAL, DIMENSION( : , : , : , : ), allocatable :: K_COMP
    REAL :: MAX_K, MIN_K

    ALLOCATE( ALPHA( CV_NONODS ))
    ALLOCATE( SUM_NOD( CV_NONODS ))
    ALLOCATE( VOLFRA_PORE_NOD( CV_NONODS ))
    ALLOCATE( K_COMP( NCOMP, CV_NONODS, NPHASE, NPHASE ))

    ! Determine a node-wise representation of porosity VOLFRA_PORE_NOD.      
    SUM_NOD = 0.0
    VOLFRA_PORE_NOD = 0.0

    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
          CV_NOD = CV_NDGLN( ( ELE - 1 ) * CV_NLOC + CV_ILOC ) 
          SUM_NOD( CV_NOD ) = SUM_NOD( CV_NOD ) + 1.0
          VOLFRA_PORE_NOD( CV_NOD ) = VOLFRA_PORE_NOD( CV_NOD ) + VOLFRA_PORE( ELE )
       END DO
    END DO
    VOLFRA_PORE_NOD = VOLFRA_PORE_NOD / SUM_NOD

    COMP_ABSORB = 0.0
    MIN_K = max( 1.e-1, MINVAL( K_COMP2( ICOMP, : , : )))
    MAX_K = MAXVAL( K_COMP2( ICOMP, : , : ) )
    CALL Calc_KComp2( cv_nonods, nphase, ncomp, icomp, KComp_Sigmoid, &
         min( 1., max( 0., Saturaold )), K_Comp2, max_k, min_k, &
         K_Comp )

    ewrite(3,*)'min,max K:', ICOMP, MIN_K, MAX_K

    DO CV_NOD = 1, CV_NONODS
       DO IPHASE = 1, NPHASE
          DO JPHASE = IPHASE + 1, NPHASE, 1

             ALPHA( CV_NOD ) = ALPHA_BETA * VOLFRA_PORE_NOD( CV_NOD ) * &
                  ( max(0.0,SATURAOLD( ( IPHASE - 1) * CV_NONODS + CV_NOD ) * &
                  DENOLD( ( IPHASE - 1) * CV_NONODS + CV_NOD )) / &
                  K_COMP( ICOMP, CV_NOD, IPHASE, JPHASE ) + &
                  max(0.0,SATURAOLD( ( JPHASE - 1) * CV_NONODS + CV_NOD ) * &
                  DENOLD( ( JPHASE - 1) * CV_NONODS + CV_NOD ))) / DT

             COMP_ABSORB( CV_NOD, IPHASE, IPHASE ) = &
                  COMP_ABSORB( CV_NOD, IPHASE, IPHASE ) + &
                  ALPHA( CV_NOD ) * &
                  K_COMP( ICOMP, CV_NOD, IPHASE, JPHASE ) 

             COMP_ABSORB( CV_NOD, IPHASE, JPHASE ) = &
                  - ALPHA( CV_NOD )                  

          END DO
       END DO
    END DO

    DO CV_NOD = 1, CV_NONODS
       DO IPHASE = 1, NPHASE
          DO JPHASE = 1, IPHASE - 1

             ALPHA( CV_NOD ) = ALPHA_BETA * VOLFRA_PORE_NOD( CV_NOD ) * &
                  ( max(0.0,SATURAOLD( ( IPHASE - 1) * CV_NONODS + CV_NOD ) * &
                  DENOLD( ( IPHASE - 1) * CV_NONODS + CV_NOD )) + &
                  max(0.0,SATURAOLD( ( JPHASE - 1) * CV_NONODS + CV_NOD ) * &
                  DENOLD( ( JPHASE - 1) * CV_NONODS + CV_NOD )) / &
                  K_COMP( ICOMP, CV_NOD, JPHASE, IPHASE )) / DT
             COMP_ABSORB( CV_NOD, IPHASE, IPHASE ) = &
                  COMP_ABSORB( CV_NOD, IPHASE, IPHASE ) + ALPHA( CV_NOD )

             COMP_ABSORB( CV_NOD, IPHASE, JPHASE ) = - ALPHA( CV_NOD ) * &
                  K_COMP( ICOMP, CV_NOD, JPHASE, IPHASE ) 

          END DO
       END DO
    END DO


    ewrite(3,*)'comp_ABSORB:',((comp_ABSORB(1,iphase,jphase), iphase=1,nphase),jphase=1,nphase)

    DEALLOCATE( ALPHA )
    DEALLOCATE( SUM_NOD )
    DEALLOCATE( VOLFRA_PORE_NOD )
    DEALLOCATE( K_COMP )

    RETURN

  END SUBROUTINE CALC_COMP_ABSORB



  SUBROUTINE CALC_COMP_DIF( NDIM, NPHASE, COMP_DIFFUSION_OPT, MAT_NONODS, &
       TOTELE, MAT_NLOC, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
       COMP_DIFFUSION, NCOMP_DIFF_COEF, COMP_DIFF_COEF, &
       X_NONODS, X, Y, Z, NU, NV, NW, U_NONODS, CV_NONODS, &
       MAT_NDGLN, U_NDGLN, X_NDGLN, &
       U_ELE_TYPE, P_ELE_TYPE ) 
    ! Calculate the diffusion coefficient COMP_DIFFUSION for current composition...
    ! based on page 136 in Reservoir-Simulation-Mathematical-Techniques-In-Oil-Recovery-(2007).pdf
    ! COMP_DIFFUSION_OPT, integer option defining diffusion coeff
    ! NCOMP_DIFF_COEF,  integer defining how many coeff's are needed to define the diffusion
    ! COMP_DIFF_COEF( NCOMP,  NCOMP_DIFF_COEF, NPHASE  )

    implicit none
    INTEGER, intent( in ) :: NDIM, NPHASE, NCOMP_DIFF_COEF, &
         COMP_DIFFUSION_OPT, MAT_NONODS, TOTELE, MAT_NLOC, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
         X_NONODS, U_NONODS, CV_NONODS, U_ELE_TYPE, P_ELE_TYPE
    REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( inout ) :: COMP_DIFFUSION
    REAL, DIMENSION( NCOMP_DIFF_COEF, NPHASE ), intent( in ) :: COMP_DIFF_COEF
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW
    INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: X_NDGLN

    ! Local variables
    REAL, DIMENSION( : ), allocatable :: UD, MAT_U
    INTEGER :: ELE, CV_NOD, MAT_NOD, IPHASE, IDIM
    REAL :: DIFF_molecular, DIFF_longitudinal, DIFF_transverse

    IF( COMP_DIFFUSION_OPT == 0 ) THEN
       COMP_DIFFUSION = 0.0
       RETURN
    ENDIF

    ALLOCATE( MAT_U( NPHASE * NDIM * CV_NONODS ))

    CALL PROJ_U2MAT( NDIM, NPHASE, COMP_DIFFUSION_OPT, MAT_NONODS, &
         TOTELE, CV_NONODS, MAT_NLOC, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
         COMP_DIFFUSION, NCOMP_DIFF_COEF, COMP_DIFF_COEF, &
         X_NONODS, X, Y, Z, NU, NV, NW, U_NONODS, MAT_NDGLN, U_NDGLN, X_NDGLN, &
         U_ELE_TYPE, P_ELE_TYPE, &
         MAT_U ) 

    ! Determine the diffusion coeff tensor COMP_DIFFUSION from MAT_U and COMP_DIFF_COEF
    ALLOCATE( UD( NDIM ))

    DO MAT_NOD = 1, MAT_NONODS
       DO IPHASE = 1, NPHASE

          DO IDIM = 1,NDIM
             UD( IDIM ) = MAT_U( MAT_NOD + ( IDIM - 1 ) * MAT_NONODS + &
                  ( IPHASE - 1 ) * NDIM * MAT_NONODS )
          END DO

          DIFF_molecular    = COMP_DIFF_COEF( 1, IPHASE )
          DIFF_longitudinal = COMP_DIFF_COEF( 2, IPHASE )
          DIFF_transverse   = COMP_DIFF_COEF( 3, IPHASE )

          CALL CALC_COMP_DIF_TEN( NDIM, UD, DIFF_molecular, DIFF_longitudinal, DIFF_transverse, &
               COMP_DIFFUSION( MAT_NOD, : , : , IPHASE )) 
       END DO
    END DO

    DEALLOCATE( MAT_U )
    DEALLOCATE( UD )

    RETURN
  end subroutine CALC_COMP_DIF



  SUBROUTINE PROJ_U2MAT( NDIM, NPHASE, COMP_DIFFUSION_OPT, MAT_NONODS, &
       TOTELE, CV_NONODS, MAT_NLOC, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
       COMP_DIFFUSION, NCOMP_DIFF_COEF, COMP_DIFF_COEF, &
       X_NONODS, X, Y, Z, NU, NV, NW, U_NONODS, MAT_NDGLN, U_NDGLN, X_NDGLN, &
       U_ELE_TYPE, P_ELE_TYPE, &
       MAT_U ) 
    ! Determine MAT_U from NU,NV,NW which are these variables mapped to material mesh. 
    use shape_functions
    use matrix_operations

    implicit none
    INTEGER, intent( in ) :: NDIM, NPHASE, NCOMP_DIFF_COEF, &
         COMP_DIFFUSION_OPT, MAT_NONODS, TOTELE, MAT_NLOC, CV_NLOC, U_NLOC, CV_SNLOC, U_SNLOC, &
         X_NONODS, U_NONODS, CV_NONODS, U_ELE_TYPE, P_ELE_TYPE
    REAL, DIMENSION( MAT_NONODS, NDIM, NDIM, NPHASE ), intent( in ) :: COMP_DIFFUSION
    REAL, DIMENSION( NCOMP_DIFF_COEF, NPHASE ), intent( in ) :: COMP_DIFF_COEF
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X, Y, Z
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: NU, NV, NW
    INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) :: U_NDGLN
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: X_NDGLN
    REAL, DIMENSION( MAT_NONODS * NPHASE * NDIM), intent( inout ) :: MAT_U
    ! Determine MAT_U from NU,NV,NW which are these variables mapped to material mesh. 

    ! Local variables
    INTEGER, DIMENSION( :, : ), allocatable :: CV_SLOCLIST, U_SLOCLIST, CV_NEILOC, FACE_ELE
    INTEGER, DIMENSION( : ), allocatable :: FINDGPTS, COLGPTS
    REAL, DIMENSION( : ),    ALLOCATABLE :: CVWEIGHT, CVWEIGHT_SHORT, DETWEI,RA,  &
         SCVFEWEIGH, SBCVFEWEIGH, SELE_OVERLAP_SCALE
    REAL, DIMENSION( :, : ), ALLOCATABLE :: CVN, CVN_SHORT, CVFEN, CVFENLX, CVFENLY, CVFENLZ, & 
         CVFENX, CVFENY, CVFENZ, CVFEN_SHORT, CVFENLX_SHORT, CVFENLY_SHORT, CVFENLZ_SHORT, & 
         CVFENX_SHORT, CVFENY_SHORT, CVFENZ_SHORT, &
         UFEN, UFENLX, UFENLY, UFENLZ, UFENX, UFENY, UFENZ, SCVFEN, SCVFENSLX, SCVFENSLY, &
         SCVFENLX, SCVFENLY, SCVFENLZ, &
         SUFEN, SUFENSLX, SUFENSLY, SUFENLX, SUFENLY, SUFENLZ, &
         SBCVFEN, SBCVFENSLX, SBCVFENSLY, SBCVFENLX, SBCVFENLY, SBCVFENLZ, &
         SBUFEN, SBUFENSLX, SBUFENSLY, SBUFENLX, SBUFENLY, SBUFENLZ, &
         MASS, INV_MASS, MASS2U, INV_MASS_NM!, FEMU
    LOGICAL, DIMENSION( :, : ), allocatable :: CV_ON_FACE,U_ON_FACE
    INTEGER :: CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE, &
         ELE, MAT_ILOC, MAT_JLOC, CV_GI, U_JLOC, MAT_KLOC, MAT_NODI, MAT_NOD, &
         U_NODJ, IPHASE, U_NODJ_IP, IDIM, JDIM, MAT_NOD_ID_IP, CV_GI_SHORT, NCOLGPTS 
    REAL :: NN, NFEMU, VOLUME, MASELE
    LOGICAL :: D1, D3, DCYL


    CALL RETRIEVE_NGI( CV_NGI, CV_NGI_SHORT, SCVNGI, SBCVNGI, NFACE, &
         NDIM, U_ELE_TYPE, CV_NLOC, U_NLOC ) 

    ALLOCATE( DETWEI( CV_NGI ))
    ALLOCATE( RA( CV_NGI ))

    ALLOCATE( CVWEIGHT( CV_NGI ))
    ALLOCATE( CVN( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFEN( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENLX( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENLY( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENLZ( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENX( CV_NLOC, CV_NGI )) 
    ALLOCATE( CVFENY( CV_NLOC, CV_NGI ))
    ALLOCATE( CVFENZ( CV_NLOC, CV_NGI ))

    ALLOCATE( CVWEIGHT_SHORT( CV_NGI_SHORT ))
    ALLOCATE( CVN_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFEN_SHORT( CV_NLOC, CV_NGI_SHORT))
    ALLOCATE( CVFENLX_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFENLY_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFENLZ_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFENX_SHORT( CV_NLOC, CV_NGI_SHORT )) 
    ALLOCATE( CVFENY_SHORT( CV_NLOC, CV_NGI_SHORT ))
    ALLOCATE( CVFENZ_SHORT( CV_NLOC, CV_NGI_SHORT ))

    ALLOCATE( UFEN( U_NLOC, CV_NGI ))
    ALLOCATE( UFENLX( U_NLOC, CV_NGI ))
    ALLOCATE( UFENLY( U_NLOC, CV_NGI ))
    ALLOCATE( UFENLZ( U_NLOC, CV_NGI ))
    ALLOCATE( UFENX( U_NLOC, CV_NGI ))
    ALLOCATE( UFENY( U_NLOC, CV_NGI ))
    ALLOCATE( UFENZ( U_NLOC, CV_NGI ))

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

    ALLOCATE( CV_SLOCLIST( NFACE, CV_SNLOC ))
    ALLOCATE( U_SLOCLIST( NFACE, U_SNLOC ))
    ALLOCATE( CV_NEILOC( CV_NLOC, SCVNGI ))

    ALLOCATE( COLGPTS( CV_NLOC * SCVNGI )) !The size of this vector is over-estimated
    ALLOCATE( FINDGPTS( CV_NLOC + 1 ))
    NCOLGPTS = 0

    ALLOCATE( CV_ON_FACE( CV_NLOC, SCVNGI ))
    ALLOCATE( U_ON_FACE( U_NLOC, SCVNGI ))

    ALLOCATE( SELE_OVERLAP_SCALE( CV_NLOC ))

    !    ALLOCATE( FEMU( U_NLOC, CV_NGI ))

    CALL CV_FEM_SHAPE_FUNS( &
                                ! Volume shape functions...
         NDIM, P_ELE_TYPE,  & 
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
         SELE_OVERLAP_SCALE ) 


    ALLOCATE( MASS( MAT_NLOC, MAT_NLOC ))
    ALLOCATE( INV_MASS( MAT_NLOC, MAT_NLOC ))
    ALLOCATE( MASS2U( MAT_NLOC, U_NLOC ))
    ALLOCATE( INV_MASS_NM( MAT_NLOC, U_NLOC ))

    D1 = ( NDIM == 1 )
    D3 = ( NDIM == 3 )
    DCYL = .FALSE. 

    Loop_Elements1: DO ELE = 1, TOTELE

       ! Calculate DETWEI,RA,NX,NY,NZ for element ELE
       CALL DETNLXR_PLUS_U( ELE, X, Y, Z, X_NDGLN, TOTELE, X_NONODS, CV_NLOC, CV_NGI, &
            CVFEN, CVFENLX, CVFENLY, CVFENLZ, CVWEIGHT, DETWEI, RA, VOLUME, D1, D3, DCYL, &
            CVFENX, CVFENY, CVFENZ, &
            U_NLOC, UFENLX, UFENLY, UFENLZ, UFENX, UFENY, UFENZ ) 

       MASELE = 0.0
       Loop_MAT_ILOC: DO MAT_ILOC = 1, MAT_NLOC

          Loop_MAT_JLOC: DO MAT_JLOC = 1, MAT_NLOC

             NN = 0.0
             DO CV_GI_SHORT = 1, CV_NGI_SHORT
                NN = NN +  CVFEN_SHORT( MAT_ILOC, CV_GI_SHORT )  * CVFEN_SHORT(  MAT_JLOC, CV_GI_SHORT ) &
                     * DETWEI( CV_GI_SHORT )
             END DO

             MASS( MAT_ILOC,MAT_JLOC)  = MASS( MAT_ILOC,MAT_JLOC) + NN

          END DO Loop_MAT_JLOC

       END DO Loop_MAT_ILOC

       CALL MATDMATINV( MASS, INV_MASS, MAT_NLOC)

       MASS2U = 0.0
       Loop_MAT_ILOC2: DO MAT_ILOC = 1, MAT_NLOC

          Loop_U_JLOC: DO U_JLOC = 1, U_NLOC

             NFEMU = 0.0 !!!!!!!!!!!!HHHHHHHHEEEEERRRREEEEEEEE!!!
             DO CV_GI = 1, CV_NGI
                ! NFEMU = NFEMU +  CVFEN( MAT_ILOC, CV_GI ) * FEMU(  U_JLOC, CV_GI ) * DETWEI( CV_GI )
                NFEMU = NFEMU +  CVFEN( MAT_ILOC, CV_GI ) * UFEN(  U_JLOC, CV_GI ) * DETWEI( CV_GI )
             END DO

             MASS2U( MAT_ILOC,U_JLOC)  = MASS2U( MAT_ILOC,U_JLOC) + NFEMU

          END DO Loop_U_JLOC

       END DO Loop_MAT_ILOC2

       INV_MASS_NM = 0.0

       DO MAT_ILOC = 1, MAT_NLOC
          DO U_JLOC = 1, U_NLOC
             DO MAT_KLOC = 1, MAT_NLOC
                INV_MASS_NM( MAT_ILOC, U_JLOC ) = INV_MASS_NM( MAT_ILOC, U_JLOC ) &
                     + INV_MASS( MAT_ILOC, MAT_KLOC ) * MASS2U( MAT_KLOC, U_JLOC )

             END DO
          END DO
       END DO

       Loop_MAT_ILOC3: DO MAT_ILOC = 1, MAT_NLOC

          MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_ILOC )

          Loop_U_JLOC2: DO U_JLOC = 1, U_NLOC

             U_NODJ = U_NDGLN(( ELE - 1 ) * U_NLOC + U_JLOC )

             Loop_IPHASE: DO IPHASE = 1, NPHASE

                U_NODJ_IP = U_NODJ + ( IPHASE - 1 ) * U_NONODS

                IDIM = 1

                MAT_NOD_ID_IP = MAT_NOD + ( IDIM - 1 ) * MAT_NONODS  + &
                     ( IPHASE - 1 ) * NDIM * MAT_NONODS

                MAT_U( MAT_NOD_ID_IP ) = MAT_U( MAT_NOD_ID_IP ) + &
                     INV_MASS_NM( MAT_ILOC, U_JLOC ) * NU( U_NODJ_IP ) 

                IF( NDIM >= 2 ) THEN
                   IDIM = 2
                   MAT_NOD_ID_IP = MAT_NOD + ( IDIM - 1 ) * MAT_NONODS + &
                        ( IPHASE- 1 ) * NDIM * MAT_NONODS
                   MAT_U( MAT_NOD_ID_IP ) = MAT_U( MAT_NOD_ID_IP ) + &
                        INV_MASS_NM( MAT_ILOC, U_JLOC ) * NV( U_NODJ_IP ) 
                ENDIF

                IF( NDIM >= 3 ) THEN
                   IDIM = 3
                   MAT_NOD_ID_IP = MAT_NOD + ( IDIM - 1 ) * MAT_NONODS + &
                        ( IPHASE - 1 ) * NDIM * MAT_NONODS
                   MAT_U( MAT_NOD_ID_IP ) = MAT_U( MAT_NOD_ID_IP ) + &
                        INV_MASS_NM( MAT_ILOC, U_JLOC ) * NW( U_NODJ_IP ) 
                ENDIF

             END DO Loop_IPHASE

          END DO Loop_U_JLOC2

       END DO Loop_MAT_ILOC3

    END DO Loop_Elements1


    ! Deallocating temporary arrays

    DEALLOCATE( DETWEI )
    DEALLOCATE( RA )
    DEALLOCATE( CVWEIGHT )
    DEALLOCATE( CVN )
    DEALLOCATE( CVFEN )
    DEALLOCATE( CVFENLX )
    DEALLOCATE( CVFENLY )
    DEALLOCATE( CVFENLZ )
    DEALLOCATE( CVFENX ) 
    DEALLOCATE( CVFENY )
    DEALLOCATE( CVFENZ )

    DEALLOCATE( CVWEIGHT_SHORT )
    DEALLOCATE( CVN_SHORT )
    DEALLOCATE( CVFEN_SHORT )
    DEALLOCATE( CVFENLX_SHORT )
    DEALLOCATE( CVFENLY_SHORT )
    DEALLOCATE( CVFENLZ_SHORT )
    DEALLOCATE( CVFENX_SHORT ) 
    DEALLOCATE( CVFENY_SHORT )
    DEALLOCATE( CVFENZ_SHORT )

    DEALLOCATE( UFEN )
    DEALLOCATE( UFENLX )
    DEALLOCATE( UFENLY )
    DEALLOCATE( UFENLZ )
    DEALLOCATE( UFENX )
    DEALLOCATE( UFENY )
    DEALLOCATE( UFENZ )

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

    DEALLOCATE( CV_SLOCLIST )
    DEALLOCATE( U_SLOCLIST )
    DEALLOCATE( CV_NEILOC )

    DEALLOCATE( COLGPTS ) !The size of this vector is over-estimated
    DEALLOCATE( FINDGPTS )

    DEALLOCATE( CV_ON_FACE )
    DEALLOCATE( U_ON_FACE )

    DEALLOCATE( SELE_OVERLAP_SCALE )

    DEALLOCATE( MASS )
    DEALLOCATE( INV_MASS )
    DEALLOCATE( MASS2U )
    DEALLOCATE( INV_MASS_NM )

    RETURN
  end subroutine PROJ_U2MAT



  SUBROUTINE CALC_COMP_DIF_TEN( NDIM, UD, DIFF_molecular, DIFF_longitudinal, DIFF_transverse, &
       DIFF_TEN )         
    ! Calculate the diffusion coefficient COMP_DIFFUSION for current composition...
    ! based on page 136 in Reservoir-Simulation-Mathematical-Techniques-In-Oil-Recovery-(2007).pdf
    implicit none

    INTEGER, intent( in ) :: NDIM
    REAL, intent( in ) :: DIFF_molecular, DIFF_longitudinal, DIFF_transverse
    REAL, DIMENSION( NDIM ), intent( in ) :: UD
    REAL, DIMENSION( NDIM, NDIM ), intent( inout ) :: DIFF_TEN

    ! Local variables...
    REAL, PARAMETER :: TOLER = 1.0E-10
    REAL, DIMENSION( :, : ), allocatable :: E, E_ident, E_OTH
    REAL :: RN, RN2
    INTEGER :: IDIM, JDIM

    ALLOCATE( E( NDIM, NDIM ))
    ALLOCATE( E_ident( NDIM, NDIM ))
    ALLOCATE( E_OTH( NDIM, NDIM ))

    RN = 0.0
    DO IDIM = 1, NDIM
       RN = RN + UD( IDIM ) **2
    END DO

    RN2 = SQRT( RN )
    RN = MAX( RN2, TOLER )

    E_ident = 0.0
    DO IDIM = 1, NDIM 

       DO JDIM = 1, NDIM
          E( IDIM, JDIM ) = UD(IDIM) * UD( JDIM ) / RN
       END DO

       E_ident( IDIM, IDIM ) = 1.0
    END DO

    E_OTH = E_ident - E

    DIFF_TEN = DIFF_molecular * E_ident  &
         + RN2 * ( DIFF_longitudinal * E + DIFF_transverse * E_OTH )


    DEALLOCATE( E )
    DEALLOCATE( E_ident )
    DEALLOCATE( E_OTH )

    RETURN

  END SUBROUTINE CALC_COMP_DIF_TEN


  subroutine Calc_KComp2( cv_nonods, nphase, ncomp, icomp, KComp_Sigmoid, &
       Satura, K_Comp, max_k, min_k, &
       K_Comp2 )
    implicit none
    integer, intent( in ) :: cv_nonods, nphase, ncomp, icomp
    logical, intent( in ) :: KComp_Sigmoid
    real, dimension( cv_nonods * nphase ), intent( in ) :: Satura
    real, dimension( ncomp, nphase, nphase ), intent( in ) :: K_Comp
    real, intent( in ) :: max_k, min_k 
    real, dimension( ncomp, cv_nonods, nphase, nphase ), intent( inout ) :: K_Comp2
    ! Local variables
    integer :: iphase, jphase, cv_nod
    real, parameter :: Width = 0.1, Err = 1.e-6, Sat = 0.9
    real :: Sat0

    K_Comp2 = 0.

    if( .not. KComp_Sigmoid ) then

       do cv_nod = 1, cv_nonods
          do iphase = 1, nphase
             do jphase = iphase + 1, nphase, 1
                K_Comp2( icomp, cv_nod, iphase, jphase ) = &
                     1. / K_Comp( icomp, iphase, jphase )
             end do
          end do
       end do

    else

       Sat0 = Sat
       do cv_nod = 1, cv_nonods

          do iphase = 1, nphase

             do jphase = 1, nphase

                if ( jphase /= iphase ) then
                   K_Comp2( icomp, cv_nod, iphase, jphase ) = &
                        1. / sigmoid_function( satura( ( iphase - 1 ) * cv_nonods + cv_nod ), &
                        Sat0, Width, min_k, max_k )
                   !   Sat0, Width, max_k, min_k )
                endif
                !K_Comp2( icomp, cv_nod, iphase, jphase ) = K_Comp2( icomp, cv_nod, jphase, iphase )
             end do
             ewrite(3,*)'icomp, iphase, sat, kcomp:',icomp, cv_nod, iphase, &
                  satura( ( iphase - 1 ) * cv_nonods + cv_nod ), &
                  ( K_Comp2( icomp, cv_nod, iphase, jphase ), &
                  jphase = 1, nphase )

          end do
       end do


    end if

    return
  end subroutine Calc_KComp2


  real function sigmoid_function( Y, Y0, Width, LowMag, UpMag )
    implicit none
    real :: Y, Y0, Width, LowMag, UpMag
    ! Local Variables
    real :: alpha
    !
    ! Width: width of the sigmoid function.
    ! The sigmoid function, varies between ( LowMag, UpMag ).
    ! Y is the variable of the function and Y0 is the centre 
    ! of the function.
    ! The function looks like:
    !             -------------
    !           /
    !          /
    ! --------

    if( Y - Y0 < - 3. * Width ) then
       sigmoid_function = UpMag
    elseif( Y - Y0 > 3. * Width ) then
       sigmoid_function = LowMag
    else
       alpha = 10. / Width
       sigmoid_function = ( UpMag - LowMag ) / &
            ( 1. + exprep( alpha * ( Y - Y0 ))) + LowMag
    end if

    return
  end function sigmoid_function

  real function exprep( M )
    implicit none
    real :: M

    if( M > 174. ) then
       exprep = 3.69e+35
    elseif( M < -180. ) then
       exprep = 0.
    else
       exprep = exp( M )
    end if

    return
  end function exprep



  SUBROUTINE CAL_COMP_SUM2ONE_SOU( V_SOURCE_COMP, CV_NONODS, NPHASE, NCOMP2, DT, ITS, NITS, &  
       MEAN_PORE_CV, SATURA, SATURAOLD, DEN, DENOLD, COMP, COMPOLD ) 
    ! make sure the composition sums to 1.0 
    implicit none
    integer, intent( in ) :: cv_nonods, nphase, ncomp2, ITS, NITS
    REAL, INTENT( IN ) :: DT
    real, dimension( cv_nonods * nphase ), intent( inout ) :: V_SOURCE_COMP
    real, dimension( cv_nonods ), intent( in ) :: MEAN_PORE_CV
    real, dimension( cv_nonods * nphase ), intent( in ) :: SATURA, SATURAOLD, DEN, DENOLD
    real, dimension( cv_nonods * nphase * ncomp2 ), intent( in ) :: COMP, COMPOLD

    ! the relaxing E.G. 0.5 is to help convergence. 
    ! =1 is full adjustment to make sure we have sum to 1. 
    ! =0 is no adjustment. 
    ! Local variables...
    !    REAL, PARAMETER :: SUM2ONE_RELAX = 0.25
    REAL, PARAMETER :: SUM2ONE_RELAX = 0.1
    INTEGER :: IPHASE, CV_NODI, ICOMP
    REAL :: RSUM

    DO IPHASE = 1, NPHASE
       DO CV_NODI = 1, CV_NONODS

          IF( ITS == NITS ) THEN ! last iteration so prepare for next time step. 

             V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                  = V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) & 
                  - SUM2ONE_RELAX * MEAN_PORE_CV( CV_NODI ) * SATURA( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                  * DEN( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) / DT

             RSUM=0.0
             DO ICOMP = 1, NCOMP2
                V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                     = V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) & 
                     + SUM2ONE_RELAX * MEAN_PORE_CV( CV_NODI ) * SATURA( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                     * DEN( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                     * COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS &
                     + ( ICOMP - 1 ) * NPHASE * CV_NONODS ) / DT 
                RSUM=RSUM+COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS &
                     + ( ICOMP - 1 ) * NPHASE * CV_NONODS )
             END DO
             !          ewrite(3,*)'IPHASE,CV_NODI,S,COMP_SUM:',IPHASE,CV_NODI,SATURA( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ),RSUM

          ELSE

             V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                  = V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) & 
                  - SUM2ONE_RELAX * MEAN_PORE_CV( CV_NODI )* SATURAOLD( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                  * DENOLD( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) / DT 

             DO ICOMP = 1, NCOMP2

                V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                     = V_SOURCE_COMP( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) & 
                     + SUM2ONE_RELAX * MEAN_PORE_CV( CV_NODI ) * SATURAOLD( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                     * DENOLD( CV_NODI + ( IPHASE - 1 ) * CV_NONODS ) &
                     * COMPOLD( CV_NODI + ( IPHASE - 1 ) * CV_NONODS &
                     + ( ICOMP - 1 ) * NPHASE * CV_NONODS ) / DT 
             END DO

          ENDIF

       END DO

    END DO

    RETURN
  END SUBROUTINE CAL_COMP_SUM2ONE_SOU



end module Compositional_Terms
