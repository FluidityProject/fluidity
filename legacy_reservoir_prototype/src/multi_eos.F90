
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

module multiphase_EOS

  use fldebug

contains

  SUBROUTINE CAL_BULK_DENSITY( NPHASE, NCOMP, CV_NONODS, CV_PHA_NONODS, NCOEF, DEN, DERIV, &
       T, P, MASSFRAC, EOS_COEFS, EOS_OPTION )
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, NCOMP, CV_NONODS, CV_PHA_NONODS, NCOEF
    REAL, DIMENSION( CV_PHA_NONODS ), intent( inout ) :: DEN, DERIV
    REAL, DIMENSION( CV_PHA_NONODS ), intent( in ) :: T
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: P
    REAL, DIMENSION( CV_PHA_NONODS * NCOMP ), intent( in ) :: MASSFRAC
    REAL, DIMENSION( NPHASE, NCOEF ), intent( in ) :: EOS_COEFS
    INTEGER, DIMENSION( NPHASE ), intent( in ) ::  EOS_OPTION
    ! Local
    INTEGER :: ICOMP, CV_NOD, IPHASE, NOD, NCOMP2

    ewrite(3,*) 'In CAL_BULK_DENSITY'

    IF( NCOMP == 0 ) THEN
       NCOMP2 = 1
    ELSE
       NCOMP2 = NCOMP
    END IF

    Loop_Comp: DO ICOMP = 1, NCOMP2

       CALL CAL_DENSITY( NPHASE, ICOMP, CV_NONODS, CV_PHA_NONODS, NCOEF, DEN, DERIV, &
            T, P, EOS_COEFS, EOS_OPTION )

    END DO Loop_Comp

    ewrite(3,*) 'Leaving CAL_BULK_DENSITY'

    RETURN
  END SUBROUTINE CAL_BULK_DENSITY

  SUBROUTINE CAL_DENSITY( NPHASE, ICOMP, CV_NONODS, CV_PHA_NONODS, NCOEF, DEN, DERIV, &
       T, P, EOS_COEFS, EOS_OPTION ) 
    ! This sub calculates the DEN ,DERIV for multiphase flows i.e. the EoS 
    ! and the derivative. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, ICOMP, CV_NONODS, CV_PHA_NONODS, NCOEF
    REAL, DIMENSION( CV_PHA_NONODS ), intent( inout ) :: DEN, DERIV
    REAL, DIMENSION( CV_PHA_NONODS ), intent( in ) :: T
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: P
    REAL, DIMENSION( NPHASE, NCOEF ), intent( in ) :: EOS_COEFS
    INTEGER, DIMENSION( NPHASE ), intent( in ) ::  EOS_OPTION
    ! Local
    REAL, PARAMETER :: TOLER = 1.0E-10
    REAL :: DEN_P1, DEN_M1, PERT_P
    INTEGER :: IPHASE, CV_NOD, NODE

    ewrite(3,*) 'In CAL_DENSITY'

    Loop_Phase: DO IPHASE = 1, NPHASE

       Loop_CV: DO CV_NOD = 1, CV_NONODS

          NODE = CV_NOD + ( IPHASE - 1 ) * CV_NONODS

          CALL CAL_DENSITY_EOS( NCOEF, DEN( NODE ), T( NODE ), P( CV_NOD ), &
               EOS_COEFS( IPHASE, : ), EOS_OPTION( IPHASE ) )  

          IF( EOS_OPTION( IPHASE ) == 1) THEN
             PERT_P = 0.001 * ( ABS(P( CV_NOD )) + ABS(EOS_COEFS( IPHASE, 1 )) )
          ELSE
             PERT_P = 0.001 * ABS(P( CV_NOD ))
          ENDIF
          PERT_P =MAX(TOLER, PERT_P)

          CALL CAL_DENSITY_EOS( NCOEF, DEN_P1, T( NODE ), P( CV_NOD ) + PERT_P, &
               EOS_COEFS( IPHASE, : ), EOS_OPTION( IPHASE ) ) 
          CALL CAL_DENSITY_EOS( NCOEF, DEN_M1, T( NODE ), P( CV_NOD ) - PERT_P, &
               EOS_COEFS( IPHASE, : ), EOS_OPTION( IPHASE ) ) 

          ! Calculating d(DEN) / dP
          DERIV( NODE ) = ( DEN_P1 - DEN_M1 ) / ( 2. * PERT_P )

       END DO Loop_CV

    END DO Loop_Phase

    ewrite(3,*) 'Leaving CAL_DENSITY'

    RETURN

  END SUBROUTINE CAL_DENSITY



  SUBROUTINE CAL_DENSITY_EOS( NCOEF, DEN, T, P , &
       EOS_COEFS, EOS_OPTION )
    IMPLICIT NONE
    INTEGER, intent( in ) :: NCOEF
    REAL, intent( inout ) :: DEN
    REAL, intent( in ) :: T, P
    REAL, DIMENSION( NCOEF ), intent( in ) :: EOS_COEFS
    INTEGER, intent( in ) :: EOS_OPTION

    IF( EOS_OPTION == 1 ) THEN
       ! stiffened gas EoS...
       DEN = ( P + EOS_COEFS( 1 )) * EOS_COEFS( 2 ) / T
    ELSEIF( EOS_OPTION == 2 ) THEN
       ! Linear EoS in pressure and temperature...

       DEN =  EOS_COEFS( 1 )                 + EOS_COEFS( 2 ) * P             + EOS_COEFS( 3 ) * T           &
            + EOS_COEFS( 4 ) * P * T         + EOS_COEFS( 5 ) * ( P **2 )     + EOS_COEFS( 6 ) * ( T **2 )   &
            + EOS_COEFS( 7 ) * ( P **2 ) * T + EOS_COEFS( 8 ) * P * ( T **2 ) + EOS_COEFS( 9 ) * (( P * T ) **2 )

    ENDIF


    RETURN
  END SUBROUTINE CAL_DENSITY_EOS




  SUBROUTINE CAL_CPDEN( NPHASE, CV_NONODS, CV_PHA_NONODS, CPDEN, DEN, NCP_COEFS, CP_COEFS, CP_OPTION, STOTEL, CV_SNLOC, SUF_CPD_BCU, SUF_D_BCU ) 

    ! This sub calculates the CPDEN ie. CP*DEN
    IMPLICIT NONE
    INTEGER, intent( in ) :: NPHASE, CV_NONODS, CV_PHA_NONODS, NCP_COEFS, STOTEL, CV_SNLOC
    REAL, DIMENSION( CV_PHA_NONODS ), intent( in ) :: DEN
    REAL, DIMENSION( CV_PHA_NONODS ), intent( inout ) :: CPDEN
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( in ) :: SUF_D_BCU
    REAL, DIMENSION( STOTEL * CV_SNLOC * NPHASE ), intent( inout ) :: SUF_CPD_BCU
    REAL, DIMENSION( NPHASE, NCP_COEFS ), intent( in ) :: CP_COEFS
    INTEGER, DIMENSION( NPHASE ), intent( in ) ::  CP_OPTION
    ! Local
    INTEGER :: IPHASE, CV_NOD, CV_NOD_IPHA, IS_IPHA, II

    Loop_Phase1: DO IPHASE = 1, NPHASE

       Loop_CV: DO CV_NOD = 1, CV_NONODS

          CV_NOD_IPHA = CV_NOD + ( IPHASE - 1 ) * CV_NONODS

          IF( CP_OPTION( IPHASE ) == 1 ) THEN
             CPDEN( CV_NOD_IPHA ) = CP_COEFS( IPHASE, 1 ) * DEN( CV_NOD_IPHA )
          ELSE
             STOP 3931
          ENDIF

       END DO Loop_CV

    END DO Loop_Phase1

    ! For the b.c's
    Loop_Surfaces: DO II = 1, STOTEL * CV_SNLOC

       Loop_Phase2: DO IPHASE = 1, NPHASE
          IS_IPHA = II + ( IPHASE - 1 ) * NPHASE

          IF( CP_OPTION( IPHASE ) == 1 ) THEN
             SUF_CPD_BCU( IS_IPHA ) = CP_COEFS( IPHASE, 1 ) * SUF_D_BCU( IS_IPHA )
          ELSE
             STOP 3932
          ENDIF

       END DO Loop_Phase2
    END DO Loop_Surfaces

    RETURN

  END SUBROUTINE CAL_CPDEN



  SUBROUTINE CAL_U_ABSORB( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
       CV_NDGLN, MAT_NDGLN, &
       NUABS_COEFS, UABS_COEFS, UABS_OPTION, U_ABSORB, VOLFRA_PORE, PERM, MOBILITY, VISCOSITY, &
       X, X_NLOC, X_NONODS, X_NDGLN, &
       OPT_VEL_UPWIND_COEFS, NOPT_VEL_UPWIND_COEFS )
    ! Calculate absorption for momentum eqns
    use matrix_operations
    !    use cv_advection
    implicit none
    INTEGER, intent( in ) :: MAT_NONODS, CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC,MAT_NLOC, &
         NUABS_COEFS, X_NLOC, X_NONODS, NOPT_VEL_UPWIND_COEFS 
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SATURA
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( NPHASE, NUABS_COEFS ), intent( in ) :: UABS_COEFS
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X
    INTEGER, DIMENSION( NPHASE ), intent( in ) ::  UABS_OPTION
    REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( inout ) :: U_ABSORB
    REAL, DIMENSION( TOTELE ), intent( in ) :: VOLFRA_PORE
    REAL, DIMENSION( TOTELE, NDIM, NDIM ), intent( in ) :: PERM
    REAL, intent( in ) :: MOBILITY
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: VISCOSITY
    REAL, DIMENSION( NOPT_VEL_UPWIND_COEFS ), intent( inout ) :: OPT_VEL_UPWIND_COEFS
    ! local variable...
    REAL, DIMENSION( :, :, :), allocatable :: U_ABSORB2
    REAL, DIMENSION( : ), allocatable :: SATURA2
    REAL :: PERT
    INTEGER :: ELE, IMAT, ICV, IPHASE, CV_ILOC, IDIM, JDIM, IJ

    ewrite(3,*) 'In CAL_U_ABSORB'

    ALLOCATE( U_ABSORB2( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ))
    ALLOCATE( SATURA2( CV_NONODS * NPHASE ))
    U_ABSORB2 = 0.
    SATURA2 = 0.

    ewrite(3,*)'b4 in cal_u_absorb2, SATURA0',size(SATURA),SATURA
    CALL CAL_U_ABSORB2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
         CV_NDGLN, MAT_NDGLN, &
         NUABS_COEFS, UABS_COEFS, UABS_OPTION, U_ABSORB, PERM, MOBILITY, VISCOSITY, &
         X, X_NLOC, X_NONODS, X_NDGLN )


    PERT = 0.0001
    SATURA2( 1 : CV_NONODS ) = SATURA( 1 : CV_NONODS ) + PERT
    IF ( NPHASE > 1 ) SATURA2( 1 + CV_NONODS : 2 * CV_NONODS ) = SATURA( 1 + CV_NONODS : 2 * CV_NONODS ) - PERT

    !ewrite(3,*)'b4 in cal_u_absorb2, MAT_NONODS, CV_NONODS, NPHASE, NDIM:',MAT_NONODS, CV_NONODS, NPHASE, NDIM
    !ewrite(3,*)'b4 in cal_u_absorb2, TOTELE, CV_NLOC, MAT_NLOC, X_NLOC, X_NONODS:',TOTELE, CV_NLOC, MAT_NLOC, X_NLOC, X_NONODS
    !ewrite(3,*)'b4 in cal_u_absorb2, SATURA',size(SATURA),SATURA
    !ewrite(3,*)'b4 in cal_u_absorb2, SATURA2',size(SATURA2),SATURA2
    !ewrite(3,*)'b4 in cal_u_absorb2, CV_NDGLN',size(CV_NDGLN),CV_NDGLN
    !ewrite(3,*)'b4 in cal_u_absorb2, MAT_NDGLN',size(MAT_NDGLN),MAT_NDGLN
    !ewrite(3,*)'b4 in cal_u_absorb2, UABS_COEFS',size(UABS_COEFS),UABS_COEFS
    !ewrite(3,*)'b4 in cal_u_absorb2, UABS_OPTION',size(UABS_OPTION),UABS_OPTION
    !ewrite(3,*)'b4 in cal_u_absorb2, VOLFRA_PORE',size(VOLFRA_PORE),VOLFRA_PORE
    !ewrite(3,*)'b4 in cal_u_absorb2, PERM',size(PERM),PERM
    !ewrite(3,*)'b4 in cal_u_absorb2, X',size(X),X
    !ewrite(3,*)'b4 in cal_u_absorb2, X_NDGLN',size(X_NDGLN),X_NDGLN
    !ewrite(3,*)'b4 in cal_u_absorb2, U_ABSORB2:', size(U_ABSORB2),U_ABSORB2
    !ewrite(3,*)'b4 in cal_u_absorb2, OPT_VEL_UPWIND_COEFS:', size(OPT_VEL_UPWIND_COEFS),OPT_VEL_UPWIND_COEFS

    CALL CAL_U_ABSORB2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA2, TOTELE, CV_NLOC, MAT_NLOC, &
         CV_NDGLN, MAT_NDGLN, &
         NUABS_COEFS, UABS_COEFS, UABS_OPTION, U_ABSORB2, PERM, MOBILITY, VISCOSITY, &
         X, X_NLOC, X_NONODS, X_NDGLN )

    ewrite(3,*)'after in cal_u_absorb, U_ABSORB2:', size(U_ABSORB2),U_ABSORB2

    DO ELE = 1, TOTELE
       DO CV_ILOC = 1, CV_NLOC
          IMAT = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC +CV_ILOC )
          ICV = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC )
          DO IPHASE = 1, NPHASE
             DO IDIM = 1, NDIM
                DO JDIM = 1, NDIM
                   IJ = ( IPHASE - 1 ) * MAT_NONODS * NDIM * NDIM + ( IMAT - 1 ) * NDIM * NDIM + &
                        ( IDIM - 1 ) * NDIM + JDIM
                   OPT_VEL_UPWIND_COEFS( IJ ) &
                        = U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ) 
                   ! This is the gradient...
                   OPT_VEL_UPWIND_COEFS( IJ + NPHASE * MAT_NONODS * NDIM * NDIM ) &
                        = ( U_ABSORB2( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ) &
                        -U_ABSORB( IMAT, IDIM + ( IPHASE - 1 ) * NDIM, JDIM + ( IPHASE - 1 ) * NDIM ))  &
                        / (SATURA2( ICV + ( IPHASE - 1 ) * CV_NONODS ) - SATURA( ICV + ( IPHASE - 1 ) * CV_NONODS )) 
                END DO
             END DO
          END DO
       END DO
    END DO

    ewrite(3,*)'in cal_u_absorb, U_ABSORB:', size(U_ABSORB),U_ABSORB
    ewrite(3,*)'in cal_u_absorb, OPT_VEL_UPWIND_COEFS:', size(OPT_VEL_UPWIND_COEFS),OPT_VEL_UPWIND_COEFS
    ewrite(3,*) 'Leaving CAL_U_ABSORB'
    RETURN

  END SUBROUTINE CAL_U_ABSORB




  SUBROUTINE CAL_U_ABSORB2( MAT_NONODS, CV_NONODS, NPHASE, NDIM, SATURA, TOTELE, CV_NLOC, MAT_NLOC, &
       CV_NDGLN, MAT_NDGLN, &
       NUABS_COEFS, UABS_COEFS, UABS_OPTION, U_ABSORB, PERM, MOBILITY, VISCOSITY, &
       X, X_NLOC, X_NONODS, X_NDGLN ) 
    ! Calculate absorption for momentum eqns
    use matrix_operations
    !    use cv_advection
    implicit none
    INTEGER, intent( in ) :: MAT_NONODS, CV_NONODS, NPHASE, NDIM, TOTELE, CV_NLOC,MAT_NLOC, &
         NUABS_COEFS, X_NLOC, X_NONODS
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SATURA
    INTEGER, DIMENSION( TOTELE * CV_NLOC ), intent( in ) :: CV_NDGLN
    INTEGER, DIMENSION( TOTELE * X_NLOC ), intent( in ) :: X_NDGLN
    INTEGER, DIMENSION( TOTELE * MAT_NLOC ), intent( in ) :: MAT_NDGLN
    REAL, DIMENSION( NPHASE, NUABS_COEFS ), intent( in ) :: UABS_COEFS
    REAL, DIMENSION( X_NONODS ), intent( in ) :: X
    INTEGER, DIMENSION( NPHASE ), intent( in ) ::  UABS_OPTION
    REAL, DIMENSION( MAT_NONODS, NDIM * NPHASE, NDIM * NPHASE ), intent( inout ) :: U_ABSORB
    REAL, DIMENSION( TOTELE, NDIM, NDIM ), intent( in ) :: PERM
    REAL, intent( in ) :: MOBILITY
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: VISCOSITY
    ! Local variable
    REAL, PARAMETER :: TOLER = 1.E-10
    INTEGER :: ELE, CV_ILOC, CV_NOD, CV_PHA_NOD, MAT_NOD, JPHA_JDIM, &
         IPHA_IDIM, IDIM, JDIM, IPHASE, II, jphase, X_NODI, MAT_NOD1, MAT_NOD2, MAT_NOD3, &
         CV_NOD1, CV_NOD2, ELE2, CV_PHA_NOD1, CV_PHA_NOD2, &
         MAT_NOD_MID, CV_ILOC2
    REAL :: SATURATION, ABS_SUM, RSUM, &
         U_ABS1, U_ABS2, U_ABS3
    REAL, DIMENSION( :, :, :), allocatable :: INV_PERM
    LOGICAL :: MEAN_ELE

    ewrite(3,*) 'In CAL_U_ABSORB2'

    ALLOCATE( INV_PERM( TOTELE, NDIM, NDIM ))

    CALL PHA_BLOCK_INV( INV_PERM, PERM, TOTELE, NDIM )

    U_ABSORB = 0.0

    Loop_ELE: DO ELE = 1, TOTELE

       Loop_CVNLOC: DO CV_ILOC = 1, CV_NLOC

          MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
          CV_NOD = CV_NDGLN(( ELE - 1) * CV_NLOC + CV_ILOC )

          Loop_NPHASE: DO IPHASE = 1, NPHASE
             CV_PHA_NOD = CV_NOD + ( IPHASE - 1 ) * CV_NONODS

             Loop_DimensionsI: DO IDIM = 1, NDIM

                Loop_DimensionsJ: DO JDIM = 1, NDIM

                   SATURATION = SATURA( CV_PHA_NOD ) 
                   IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM 
                   JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM 
                   !JPHA_JDIM = IPHA_IDIM

                   Case_UABSOption: SELECT CASE( UABS_OPTION( IPHASE ) )
                   CASE( 0 ) 
                      ! no absorption option
                      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = 0.0

                   CASE( 1 )
                      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = INV_PERM( ELE, IDIM, JDIM )

                   CASE( 2 ) ! A standard polynomial representation of relative permeability             
                      ABS_SUM = 0.
                      DO II = 1, NUABS_COEFS
                         ABS_SUM = ABS_SUM + UABS_COEFS( IPHASE, II) * SATURATION** (II - 1 )
                      END DO
                      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = INV_PERM( ELE, IDIM, JDIM ) / &
                           MAX( TOLER, ABS_SUM )

                   CASE( 3 ) ! From Tara and Martin's notes for relative permeability

                      CALL ABS3( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), MOBILITY, &
                           INV_PERM( ELE, IDIM, JDIM ), min(1.0,max(0.0,SATURA(CV_NOD))), IPHASE)

                   CASE( 4 ) 
                      ABS_SUM = 0.0
                      DO II = 1, NUABS_COEFS
                         ABS_SUM = ABS_SUM + UABS_COEFS( IPHASE, II) * SATURATION** (II - 1 )
                      END DO
                      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = ABS_SUM

                   CASE( 5 ) ! same as option 3 but larger phase 2 abs when S^1=1.0

                      CALL ABS3_HP1( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), MOBILITY, &
                           INV_PERM( ELE, IDIM, JDIM ), min(1.0,max(0.0,SATURA(CV_NOD))), IPHASE )

                   CASE( 6 ) ! From Tara and Martin's notes for relative permeability

                      CALL ABS6( U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ), MOBILITY, &
                           INV_PERM( ELE, IDIM, JDIM ), min(1.0,max(0.0,SATURA(CV_NOD))), IPHASE)

                   CASE DEFAULT
                      ewrite(3,*) 'Unknown Option -- UABS_OPTION( IPHASE )'
                      ABS_SUM = 0.0
                      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = ABS_SUM
                      !                      stop 1023

                   END SELECT Case_UABSOption


                   !                   ewrite(3,*)'ELE, IDIM, JDIM, iphase, U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ):', &
                   !                        ELE, IDIM, JDIM, iphase, U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )
                   !U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) = INV_PERM( ELE, IDIM, JDIM ) / &
                   !     MAX( TOLER, UABS_COEFS( IPHASE, 1 ) + UABS_COEFS( IPHASE, 2) * SATURATION + &
                   !     UABS_COEFS( IPHASE, 3 ) * SATURATION **2 )
                   !ewrite(3,*) 'MAT_NOD, IPHA, JPHA, sat., U_ABSORB:', MAT_NOD, IPHA_IDIM, JPHA_JDIM, saturation, U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )

                END DO Loop_DimensionsJ

             END DO Loop_DimensionsI

          END DO Loop_NPHASE

       END DO Loop_CVNLOC

    END DO Loop_ELE
    !    stop 23

    Conditional_NotUsed1: IF(.false.) THEN 
       !    IF(.true.) THEN
       ! Use upwind absorptions. 
       !       IF((CV_NONODS==CV_NLOC*TOTELE).OR.(CV_NLOC==1)) THEN
       IF(.true.) THEN
          DO ELE = 2,TOTELE
             DO CV_ILOC=1,1
                MAT_NOD=MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                CV_NOD=CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC)

                DO IPHA_IDIM=1,NDIM*NPHASE
                   DO JPHA_JDIM=1,NDIM*NPHASE
                      if(IPHA_IDIM==JPHA_JDIM) then
                         CALL ABS3(rsum, MOBILITY, &
                              1.0, 0.5*(SATURA(CV_NOD)+SATURA(CV_NOD-1)), IPHA_IDIM)
                         U_ABSORB( MAT_NOD,   IPHA_IDIM, JPHA_JDIM )=RSUM
                         U_ABSORB( MAT_NOD-1, IPHA_IDIM, JPHA_JDIM )=RSUM
                      endif
                   END DO
                END DO
             END DO
          END DO
       ELSE
          DO ELE = TOTELE,1,-1
             DO CV_ILOC=CV_NLOC,1,-1
                MAT_NOD=MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                ELE2=ELE
                CV_ILOC2=CV_ILOC-1
                IF(CV_ILOC==1) THEN
                   CV_ILOC2=CV_NLOC-1
                   ELE2=ELE-1
                   IF(ELE2==0) THEN
                      CV_ILOC2=CV_ILOC
                      ELE2=ELE
                   ENDIF
                ENDIF
                MAT_NOD2=MAT_NDGLN(( ELE2 - 1 ) * MAT_NLOC + CV_ILOC2)
                DO IPHA_IDIM=1,NDIM*NPHASE
                   DO JPHA_JDIM=1,NDIM*NPHASE
                      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )=U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM )
                   END DO
                END DO
             END DO
          END DO
       ENDIF
    ENDIF Conditional_NotUsed1

    Conditional_NotUsed2: IF(.FALSE.) THEN
       DO ELE = 1, TOTELE
          DO IPHA_IDIM=1,NDIM*NPHASE
             DO JPHA_JDIM=1,NDIM*NPHASE
                RSUM=0.0
                !          RSUM=1.e+10
                DO CV_ILOC = 1, CV_NLOC
                   MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                   X_NODI = X_NDGLN(( ELE - 1 ) * X_NLOC + CV_ILOC)
                   U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )=X(X_NODI)
                END DO
             END DO
          END DO
       END DO
    ENDIF Conditional_NotUsed2

    MEAN_ELE=.false.
    IF( MEAN_ELE ) THEN
       ! Average absoption throughout element -using the harmonic average
       DO ELE = 1, TOTELE
          DO IPHA_IDIM=1,NDIM*NPHASE
             DO JPHA_JDIM=1,NDIM*NPHASE
                RSUM=0.0
                !          RSUM=1.e+10
                DO CV_ILOC = 1, CV_NLOC
                   MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                   !      RSUM=RSUM+ (1./max(1.e-10,U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )))/REAL(CV_NLOC+1)
                   !      IF(CV_ILOC==2) RSUM=RSUM+ (1./max(1.e-10,U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )))/REAL(CV_NLOC+1)
                   !      RSUM=RSUM+ U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )/REAL(CV_NLOC)
                   RSUM=RSUM+ U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )/REAL(CV_NLOC+1)
                   IF(CV_ILOC==2) RSUM=RSUM+ U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )/REAL(CV_NLOC+1)
                   !      RSUM=max(RSUM, U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ))
                   !      RSUM=min(RSUM, U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ))
                END DO
                MAT_NOD1 = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + 1)
                MAT_NOD2 = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + 2)
                MAT_NOD3 = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + 3)

                U_ABS1= &
                     0.66666* U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM ) &
                     +0.33333* U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM )

                U_ABS2= &
                     0.5*0.5* U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM ) &
                     +0.5* U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM ) &
                     +0.5*0.5* U_ABSORB( MAT_NOD3, IPHA_IDIM, JPHA_JDIM ) 

                U_ABS3= &
                     0.33333* U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM ) &
                     +0.66666* U_ABSORB( MAT_NOD3, IPHA_IDIM, JPHA_JDIM )

                U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM )=U_ABS1
                U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM )=U_ABS2
                U_ABSORB( MAT_NOD3, IPHA_IDIM, JPHA_JDIM )=U_ABS3

                DO CV_ILOC = 1, -CV_NLOC
                   MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                   MAT_NOD1 = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + 3)
                   !      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )=  &
                   !    0.25*U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM )  &
                   !    +0.25*U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )
                   U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )=RSUM
                   !      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )=0.65*1./RSUM &
                   !              + 0.35*U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )
                   !      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )=0.5*1./RSUM &
                   !              + 0.5*U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )
                   !      U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )=1./RSUM 
                END DO
             END DO
          END DO
       END DO

    ENDIF

    !    ewrite(3,*)'CV_NONODs,TOTELE*CV_NLOC,CV_NLOC:',CV_NONODs,TOTELE*CV_NLOC,CV_NLOC
    !     stop 383

    !    IF(CV_NONODS/=TOTELE*CV_NLOC) THEN
    Conditional_NotUsed3: IF(.false.) THEN
       IF(CV_NLOC==3) THEN
          !    IF(.false.) THEN
          DO ELE = 1, -TOTELE
             DO IPHASE = 1, NPHASE
                DO IDIM = 1, NDIM
                   DO JDIM = 1, NDIM
                      IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM 
                      JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM 
                      DO CV_ILOC=1,CV_NLOC
                         MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                         CV_NOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC)
                         SATURATION=SATURA(CV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                         U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )  &
                              = U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) &
                              / max(0.03,SATURATION)
                      END DO
                   END DO
                END DO
             END DO
          END DO
          DO ELE = 1, TOTELE
             DO IPHASE = 1, NPHASE
                DO IDIM = 1, NDIM
                   DO JDIM = 1, NDIM
                      IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM 
                      JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM 
                      MAT_NOD1 = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + MAT_NLOC)
                      MAT_NOD2 = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + 1)
                      MAT_NOD_MID = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + 2)

                      CV_NOD1 = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_NLOC)
                      CV_NOD2 = CV_NDGLN(( ELE - 1 ) * CV_NLOC + 1)

                      CV_PHA_NOD1 = CV_NOD1 + ( IPHASE - 1 ) * CV_NONODS
                      CV_PHA_NOD2 = CV_NOD2 + ( IPHASE - 1 ) * CV_NONODS
                      if(.true.) then
                         !        stop 56
                         ! max...
                         IF(U_ABSORB( MAT_NOD_MID, IPHA_IDIM, JPHA_JDIM )  &
                              >U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM ) ) THEN
                            IF(U_ABSORB( MAT_NOD_MID, IPHA_IDIM, JPHA_JDIM )  &
                                 >U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM ) ) THEN
                               U_ABSORB( MAT_NOD_MID, IPHA_IDIM, JPHA_JDIM )= &
                                    MAX(U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM ), &
                                    U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM ))
                            ENDIF
                         ENDIF
                         ! min...
                         IF(U_ABSORB( MAT_NOD_MID, IPHA_IDIM, JPHA_JDIM )  &
                              <U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM ) ) THEN
                            IF(U_ABSORB( MAT_NOD_MID, IPHA_IDIM, JPHA_JDIM )  &
                                 <U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM ) ) THEN
                               U_ABSORB( MAT_NOD_MID, IPHA_IDIM, JPHA_JDIM )= &
                                    MIN(U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM ), &
                                    U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM ))
                            ENDIF
                         ENDIF
                      ENDIF
                      if(.false.) then
                         U_ABSORB( MAT_NOD_MID, IPHA_IDIM, JPHA_JDIM )= &
                              0.5*U_ABSORB( MAT_NOD1, IPHA_IDIM, JPHA_JDIM ) &
                              +0.5*U_ABSORB( MAT_NOD2, IPHA_IDIM, JPHA_JDIM )
                      ENDIF
                   END DO
                END DO
             END DO
          END DO
          DO ELE = 1, -TOTELE
             DO IPHASE = 1, NPHASE
                DO IDIM = 1, NDIM
                   DO JDIM = 1, NDIM
                      IPHA_IDIM = ( IPHASE - 1 ) * NDIM + IDIM 
                      JPHA_JDIM = ( IPHASE - 1 ) * NDIM + JDIM 
                      DO CV_ILOC=1,CV_NLOC
                         MAT_NOD = MAT_NDGLN(( ELE - 1 ) * MAT_NLOC + CV_ILOC)
                         CV_NOD = CV_NDGLN(( ELE - 1 ) * CV_NLOC + CV_ILOC)
                         SATURATION=SATURA(CV_NOD + ( IPHASE - 1 ) * CV_NONODS )
                         U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM )  &
                              = U_ABSORB( MAT_NOD, IPHA_IDIM, JPHA_JDIM ) &
                              * max(0.03,SATURATION)
                      END DO
                   END DO
                END DO
             END DO
          END DO
       ENDIF
    ENDIF Conditional_NotUsed3

    ewrite(3,*)'NUABS_COEFS,UABS_COEFS:',NUABS_COEFS,UABS_COEFS
    ewrite(3,*)'U_ABSORB:'
    do mat_nod =1, -mat_nonods
       do iphase=1,nphase
          ewrite(3,*)mat_nod,iphase,(U_ABSORB( MAT_NOD, iphase, jphase ),jphase=1,2) 
       end do
    end do
    !    ewrite(3,*)'U_ABSORB:',U_ABSORB
    !    stop 4941
    DEALLOCATE( INV_PERM )

    ewrite(3,*) 'Leaving CAL_U_ABSORB2'

    RETURN

  END SUBROUTINE CAL_U_ABSORB2



  SUBROUTINE ABS3_HP1( ABS, MOBILITY, INV_PERM, SAT, IPHASE )
    IMPLICIT NONE
    REAL, intent( inout) :: ABS
    REAL, intent( in ) :: MOBILITY, INV_PERM, SAT
    INTEGER, intent( in ) :: IPHASE
    ! Local variables...
    REAL :: VISC1, VISC2, S_GC, S_OR, REF_MOBILITY, &
         KR1, KR2, KR, VISC, SATURATION, ABS_SUM, SAT2

    VISC1 = 1.0
    S_GC = 0.1
    S_OR = 0.3
    REF_MOBILITY = MOBILITY

    SATURATION = SAT
    IF( IPHASE == 2 ) SATURATION=1.-SAT

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
       !        ABS =  min( 1.e+6, ABS)
       if(saturation < 0.1) then
          ABS = (1. + max(100.*(0.1-saturation),0.0)) * ABS
          !                        abs=abs+100000.*exp(90.0*(0.1-sat))
          !          ABS = (0.001 + max(100.*(0.1-saturation),0.0)) * ABS
          !          ABS = saturation*ABS/10.
       endif
       ! new ************        
       if(sat>0.7) then
          ABS = max(0.01,0.7*((1.-sat)/0.3 ))
          !    abs = 0.7*exp(-40.*(sat-0.7))
       endif
    else
       !        ABS = min( 1.e+5, ABS)
       !        ABS = min( 2.013e+5, ABS)
       ABS = min( 4.0e+5, ABS)
       !        ABS = min( 1.e+7, ABS)
       if(saturation < 0.3) then
          ABS = (1. + max(100.*(0.3-saturation),0.0)) * ABS
          !   print *,'saturation,ABS :',saturation,ABS
          abs=abs+100000.*exp(30.0*(sat-0.7))
          !   print *,'sat,abs:',sat,abs
       endif
       ! new ************        
       !        if(sat<0.1) then
       !          ABS = max(0.01,9.0*(10.0*sat))
       !        endif
    endif
    RETURN   
  END SUBROUTINE ABS3_HP1





  SUBROUTINE ABS3( ABSP, MOBILITY, INV_PERM, SAT, IPHASE )
    IMPLICIT NONE
    REAL, intent( inout ) :: ABSP
    REAL, intent( in ) :: MOBILITY, SAT, INV_PERM
    INTEGER, intent( in ) :: IPHASE
    ! Local variables...
    REAL :: VISC1, VISC2, S_GC, S_OR, REF_MOBILITY, &
         KR1, KR2, KR, VISC, SATURATION, ABS_SUM, SAT2

    VISC1 = 1.0
    S_GC = 0.1
    S_OR = 0.3
    REF_MOBILITY = MOBILITY

    SATURATION = SAT
    IF( IPHASE == 2 ) SATURATION = 1. - SAT

    IF( SAT < S_GC ) THEN
       KR1 = 0.0
    ELSE IF( SAT > 1. -S_OR ) THEN
       KR1 = 1.0
    ELSE
       KR1 = ( ( SAT - S_GC) / ( 1. - S_GC - S_OR ))**2
    ENDIF

    SAT2 = 1.0 - SAT
    IF( SAT2 < S_OR ) THEN
       KR2 = 0.0
    ELSEIF( SAT2 > 1. - S_GC ) THEN
       KR2 = 1.0
    ELSE
       KR2 = ( ( SAT2 - S_OR ) / ( 1. - S_GC - S_OR ))**2
    ENDIF
    VISC2 = REF_MOBILITY

    IF( IPHASE == 1 ) THEN
       KR = KR1
       VISC = VISC1
    ELSE
       KR = KR2
       VISC = VISC2
    ENDIF

    ABS_SUM = KR / MAX( 1.e-6, VISC * max( 0.01, SATURATION ))

    ABSP = INV_PERM / MAX( 1.e-6, ABS_SUM )

    if( iphase == 1 ) then
       ABSP =  min( 1.e+4, ABSP )
       if( saturation < 0.1 ) then
          ABSP = ( 1. + max( 100. * ( 0.1 - saturation ), 0.0 )) * ABSP
       endif
    else
       ABSP = min( 1.e+5, ABSP )
       if( saturation < 0.3 ) then
          ABSP = ( 1. + max( 100. * ( 0.3 - saturation ), 0.0 )) * ABSP
       endif
    endif

    RETURN
  END SUBROUTINE ABS3




  SUBROUTINE ABS6( ABSP, MOBILITY, INV_PERM, SAT, IPHASE )
    IMPLICIT NONE
    REAL, intent( inout ) :: ABSP
    REAL, intent( in ) :: MOBILITY, SAT, INV_PERM
    INTEGER, intent( in ) :: IPHASE
    ! Local variables...
    REAL :: VISC1, VISC2, S_GC, S_OR, REF_MOBILITY, &
         KR1, KR2, KR, VISC, SATURATION, ABS_SUM, SAT2

    VISC1 = 1.0
    S_GC = 0.1
    S_OR = 0.3
    !REF_MOBILITY = 6.0
    REF_MOBILITY = MOBILITY

    SATURATION = SAT
    IF( IPHASE == 2 ) SATURATION = 1. - SAT

    IF( SAT < S_GC ) THEN
       KR1 = 0.0
    ELSE IF( SAT > 1. -S_OR ) THEN
       KR1 = 1.0
    ELSE
       KR1 = ( ( SAT - S_GC) / ( 1. - S_GC - S_OR ))**2
    ENDIF

    SAT2 = 1.0 - SAT
    IF( SAT2 < S_OR ) THEN
       KR2 = 0.0
    ELSEIF( SAT2 > 1. - S_GC ) THEN
       KR2 = 1.0
    ELSE
       KR2 = ( ( SAT2 - S_OR ) / ( 1. - S_GC - S_OR ))**2
    ENDIF
    VISC2 = REF_MOBILITY

    IF( IPHASE == 1 ) THEN
       KR = KR1
       VISC = VISC1
    ELSE
       KR = KR2
       VISC = VISC2
    ENDIF

    ABS_SUM = KR / MAX( 1.e-6, VISC * max( 0.01, SATURATION ))

    ABSP = INV_PERM / MAX( 1.e-6, ABS_SUM )

    if( iphase == 1 ) then
       ABSP =  min( 1.e+4, ABSP )
       if( saturation < 0.1 ) then
          ABSP = ( 1. + max( 100. * ( 0.1 - saturation ), 0.0 )) * ABSP
       endif
    else
       ABSP = min( 1.e+5, ABSP )
       if( saturation < 0.3 ) then
          ABSP = ( 1. + max( 100. * ( 0.3 - saturation ), 0.0 )) * ABSP
       endif
    endif

    RETURN
  END SUBROUTINE ABS6



  SUBROUTINE CALC_CAPIL_PRES( CAPIL_PRES_OPT, CV_NONODS, NPHASE, NCAPIL_PRES_COEF, &
       CAPIL_PRES_COEF, IPLIKE_GRAD_SOU, PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD, SATURA )

    ! CAPIL_PRES_OPT is the capillary pressure option for deciding what form it might take. 
    ! CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE, NPHASE ) are the coefficients 
    ! Capillary pressure coefs have the dims CAPIL_PRES_COEF( NCAPIL_PRES_COEF, NPHASE,NPHASE )
    ! used to calulcate the capillary pressure. 

    IMPLICIT NONE
    INTEGER, intent( in ) :: IPLIKE_GRAD_SOU, CAPIL_PRES_OPT, CV_NONODS, NPHASE, NCAPIL_PRES_COEF
    REAL, DIMENSION( NCAPIL_PRES_COEF, NPHASE, NPHASE ), intent( in ) :: CAPIL_PRES_COEF
    REAL, DIMENSION( IPLIKE_GRAD_SOU*CV_NONODS * NPHASE ), intent( inout ) :: &
         PLIKE_GRAD_SOU_COEF, PLIKE_GRAD_SOU_GRAD
    REAL, DIMENSION( CV_NONODS * NPHASE ), intent( in ) :: SATURA
    ! Local Variables
    INTEGER :: IPHASE, JPHASE


    IF( IPLIKE_GRAD_SOU /= 1 ) THEN
       EWRITE(3,*),'PROBLEM WITH THE CAPILLARY OPTION SET UP' 
       STOP 3831
    ENDIF


    Case_CapillaryPressure: SELECT CASE( CAPIL_PRES_OPT )

    CASE DEFAULT

       EWRITE(3,*),'PROBLEM WITH THE CAPILLARY OPTION SET UP - OPTION NOT AVAILABLE' 
       STOP 3832

    CASE( 1 )
       ! Capillary pressure coefs 
       ! In this form PLIKE_GRAD_SOU_GRAD is the actual capillary pressure. 
       PLIKE_GRAD_SOU_COEF = 1.0
       PLIKE_GRAD_SOU_GRAD = 0.0

       DO IPHASE = 1, NPHASE

          DO JPHASE = 1, NPHASE

             PLIKE_GRAD_SOU_GRAD( 1 + ( IPHASE - 1 ) * CV_NONODS : IPHASE * CV_NONODS ) = &
                  PLIKE_GRAD_SOU_GRAD( 1 + ( IPHASE - 1 ) * CV_NONODS : IPHASE * CV_NONODS ) + &
                  CAPIL_PRES_COEF( 1, IPHASE, JPHASE ) * &
                  MAX( SATURA( 1 + ( JPHASE - 1 ) * CV_NONODS : JPHASE * CV_NONODS ), 0.0 ) &
                  ** CAPIL_PRES_COEF( 2, IPHASE, JPHASE )

          END DO

       END DO

    end SELECT Case_CapillaryPressure

    RETURN
  END SUBROUTINE CALC_CAPIL_PRES


end module multiphase_EOS
