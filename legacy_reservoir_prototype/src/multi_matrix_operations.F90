
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

!!!!========================================!!!!
!!!!           MATRIX OPERATIONS            !!!!
!!!!========================================!!!!
module matrix_operations

  use fldebug

contains

  SUBROUTINE MULMAT( VEC, CT, U, FREDOP, NONODS, NCOLCT, FINDCT, COLCT ) 
    IMPLICIT NONE
    ! perform VEC = C^T * P, matrix vector multiplication

    INTEGER, intent( in ) :: FREDOP, NONODS, NCOLCT
    REAL, DIMENSION( FREDOP ), intent( inout ) :: VEC 
    REAL, DIMENSION( NCOLCT ), intent( inout ) :: CT
    REAL, DIMENSION( NONODS ), intent( in ) :: U
    INTEGER, DIMENSION( NONODS + 1 ), intent( in ) :: FINDCT 
    INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
    ! Local variables
    INTEGER :: PNOD, COUNT, COL

    VEC( 1: FREDOP ) = 0.0

    DO PNOD = 1, FREDOP

       DO COUNT = FINDCT( PNOD ), FINDCT( PNOD + 1 ) - 1
          COL = COLCT( COUNT )
          VEC( PNOD ) = VEC( PNOD ) + CT( COUNT ) * U( COL )
       END DO

    END DO

    RETURN
  END SUBROUTINE MULMAT


  SUBROUTINE MULMATtransp( VEC, CT, P, FREDOP, NONODS, NCOLCT, FINDCT, COLCT ) 
    IMPLICIT NONE
    ! perform VEC = C * P, matrix vector multiplication

    INTEGER, intent( in ) :: FREDOP, NONODS, NCOLCT
    REAL, DIMENSION( NONODS ), intent( inout ) :: VEC 
    REAL, DIMENSION( NCOLCT ), intent( inout ) :: CT
    REAL, DIMENSION( FREDOP ), intent( in ) :: P
    INTEGER, DIMENSION( NONODS + 1 ), intent( in ) :: FINDCT 
    INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
    ! Local variables
    INTEGER :: PNOD, COUNT, COL

    VEC( 1: NONODS ) = 0.0

    DO PNOD = 1, FREDOP
       DO COUNT = FINDCT( PNOD ), FINDCT( PNOD + 1 ) - 1
          COL = COLCT( COUNT )
          VEC( COL ) = VEC( COL ) +CT( COUNT ) * P( PNOD )
       END DO
    END DO

    RETURN
  END SUBROUTINE MULMATtransp


  SUBROUTINE MATDMATINV( DMAT, DMATINV, NLOC )
    ! calculate DMATINV
    IMPLICIT NONE
    INTEGER, intent( in ) :: NLOC
    REAL, DIMENSION( NLOC, NLOC ), intent( inout ) :: DMAT, DMATINV
    ! Local variables
    REAL, DIMENSION( : , : ), allocatable :: MAT, MAT2
    REAL, DIMENSION( : ), allocatable :: X, B

    ALLOCATE( MAT( NLOC, NLOC ))
    ALLOCATE( MAT2( NLOC, NLOC ))
    ALLOCATE( X( NLOC ))
    ALLOCATE( B( NLOC ))

    DMATINV = DMAT
    CALL MATINV( DMATINV, NLOC, NLOC, MAT, MAT2, X, B )

    DEALLOCATE( MAT )
    DEALLOCATE( MAT2 )
    DEALLOCATE( X )
    DEALLOCATE( B )

    RETURN
  END SUBROUTINE MATDMATINV
  !


  SUBROUTINE MATINV( A, N, NMAX, MAT, MAT2, X, B)
    ! This sub finds the inverse of the matrix A and puts it back in A. 
    ! MAT, MAT2, X and B are working vectors. 
    IMPLICIT NONE
    INTEGER, intent( in ) :: N, NMAX
    REAL, DIMENSION( NMAX, NMAX ), intent( inout ) ::  A
    REAL, DIMENSION( N, N ), intent( inout ) :: MAT, MAT2
    REAL, DIMENSION( N ), intent( inout ) :: X, B
    ! Local variables
    INTEGER :: ICOL, IM, JM

    Loop_Boolean: IF( .true. ) THEN
       ! Solve MAT XL=BX (NB BDIAG is overwritten)
       ICOL = 1
       B( 1 : N ) = 0.
       B( ICOL ) = 1.

       DO IM = 1, N

          DO JM = 1, N

             MAT( IM, JM ) = A( IM, JM )

          END DO

       END DO

       CALL SMLINNGOT( MAT, X, B, N, N, .FALSE. ) ! X contains the column ICOL of inverse

       DO IM = 1, N
          MAT2( IM, ICOL ) = X( IM )
       END DO

       DO ICOL = 2, N ! Form column ICOL of the inverse. 
          B = 0.
          B( ICOL ) = 1.0 ! Solve MAT X=B (NB MAT is overwritten).  

          CALL SMLINNGOT( MAT, X, B, N, N, .TRUE. ) ! X contains the column ICOL of inverse

          DO IM = 1, N
             MAT2( IM, ICOL ) = X( IM )
          END DO
       END DO

    ELSE
       !
       Loop_Col: DO ICOL = 1, N ! Form column ICOL of the inverse. 

          Loop_IM1: DO IM = 1, N

             B( IM ) = 0.
             Loop_JM1: DO JM = 1, N
                MAT( IM, JM ) = A( IM , JM )
             END DO Loop_JM1
          END DO Loop_IM1

          B( ICOL ) = 1.0

          ! Solve MAT X=B (NB MAT is overwritten).  

          CALL SMLINN( MAT, X, B, N, N )

          ! X contains the column ICOL of inverse
          DO IM = 1, N
             MAT2( IM, ICOL ) = X( IM )
          END DO
          !
       END DO Loop_Col

    ENDIF Loop_Boolean
    !
    ! Set A to MAT2
    Loop_IM2: DO IM = 1, N
       Loop_JM2: DO JM = 1, N
          A( IM, JM ) = MAT2( IM, JM )
       END DO Loop_JM2
    END DO Loop_IM2

    RETURN
  END SUBROUTINE MATINV
  !
  !

  SUBROUTINE MATMASSINV( MASINV, MMAT, NONODS, NLOC, TOTELE )
    IMPLICIT NONE
    INTEGER, intent( in ) :: NONODS, NLOC, TOTELE
    REAL, DIMENSION( NONODS, NONODS ), intent( inout ) ::  MASINV, MMAT
    ! matrix is   AGI   BGI
    !             CGI   DGI
    ! Local variables
    REAL :: AGI, BGI, CGI, DGI, DETJ
    REAL :: AI11, AI12, AI21, AI22
    INTEGER :: ELE, GLOBI1, GLOBI2

    MASINV( 1 : NONODS, 1 : NONODS ) = 0.0

    Loop_ELE: DO ELE = 1, TOTELE

       GLOBI1 = ( ELE - 1 ) * NLOC + 1
       GLOBI2 = ( ELE - 1 ) * NLOC + 2

       AGI = MMAT( GLOBI1, GLOBI1 )
       BGI = MMAT( GLOBI1, GLOBI2 )
       CGI = MMAT( GLOBI2, GLOBI1 )
       DGI = MMAT( GLOBI2, GLOBI2 )
       DETJ = AGI * DGI - BGI * CGI
       AI11 = DGI / DETJ
       AI12 = -BGI / DETJ
       AI21 = -CGI / DETJ
       AI22 = AGI / DETJ

       MASINV( GLOBI1, GLOBI1 ) = AI11
       MASINV( GLOBI1, GLOBI2 ) = AI12
       MASINV( GLOBI2, GLOBI1 ) = AI21
       MASINV( GLOBI2, GLOBI2 ) = AI22

    END DO Loop_ELE

    RETURN
  END SUBROUTINE MATMASSINV

  !

  SUBROUTINE SMLINN( A, X, B, NMX, N )
    IMPLICIT NONE
    INTEGER, intent( in ) :: NMX, N
    REAL, DIMENSION( NMX, NMX ), intent( inout ) :: A
    REAL, DIMENSION( NMX ), intent( inout ) :: X
    REAL, DIMENSION( NMX ), intent( in ) :: B
    ! Local
    REAL :: R
    INTEGER ::  K, I, J
    !     Form X = A^{-1} B
    !     Useful subroutine for inverse. This sub overwrites matrix A. 
    Loop_K: DO K = 1, N - 1
       Loop_I: DO I = K + 1, N
          A( I, K ) = A( I, K ) / A( K, K )
       END DO Loop_I

       Loop_J: DO J = K + 1, N
          Loop_I1: DO I = K + 1, N
             A( I, J ) = A( I, J ) - A( I, K ) * A( K, J )
          END DO Loop_I1
       END DO Loop_J

    END DO Loop_K
    !     
    !     Solve L_1 x=b
    Loop_I2: DO I = 1, N
       R = 0.
       Loop_J2: DO J = 1, I - 1
          R = R + A( I, J ) * X( J )
       END DO Loop_J2
       X( I ) = B( I ) - R
    END DO Loop_I2
    !     
    !     Solve U x=y
    Loop_I3: DO I = N, 1, -1
       R = 0.
       Loop_J3: DO J = I + 1, N
          R = R + A( I, J) * X( J )
       END DO Loop_J3
       X( I ) = ( X( I ) - R ) / A( I, I )
    END DO Loop_I3

    RETURN

  END SUBROUTINE SMLINN
  !     

  SUBROUTINE SMLINNGOT( A, X, B, NMX, N, GOTDEC )
    IMPLICIT NONE
    INTEGER :: NMX, N
    REAL, DIMENSION( NMX, NMX ), intent( inout ) :: A
    real, DIMENSION( NMX ), intent( inout ) :: X ! inout as n might be < nmx
    real, DIMENSION( NMX ), intent( in ) ::  B
    LOGICAL, intent( in ) :: GOTDEC
    ! Local     
    REAL :: R
    INTEGER :: K, I, J

    ! IF GOTDEC then assume we have already got the LU decomposition in A
    ! Form X = A^{-1} B ;  Useful subroutine for inverse
    ! This sub overwrites the matrix A. 

    Cond_LUDecomp: IF( .NOT. GOTDEC ) THEN

       DO K = 1, N - 1

          DO I = K + 1, N, 1
             A( I, K ) = A( I, K ) / A( K, K )
          END DO

          DO J = K + 1, N, 1

             DO I = K + 1, N, 1
                A( I, J ) = A( I, J ) - A( I, K ) * A( K, J )
             END DO

          END DO

       END DO

    ENDIF Cond_LUDecomp
    !     
    ! Solve L_1 x=b
    DO I = 1, N
       R = 0.

       DO J = 1, I - 1
          R = R + A( I, J ) * X( J )
       END DO

       X( I ) = B( I ) - R

    END DO
    !     
    ! Solve U x=y
    DO I = N, 1, -1
       R = 0.

       DO J = I + 1, N, 1
          R = R + A( I, J ) * X( J )
       END DO

       X( I ) = ( X( I ) - R ) / A( I, I )

    END DO

    RETURN
  END SUBROUTINE SMLINNGOT
  !     


  SUBROUTINE ABMATRIXMUL( AB, A, NONODS1, NONODS2, &
       B, NONODS3, NONODS4 )
    !
    ! Perform matrix matrix multiplication: AB = A * B
    IMPLICIT NONE
    INTEGER, intent( in ) :: NONODS1, NONODS2, NONODS3, NONODS4
    REAL, DIMENSION( NONODS1, NONODS4 ), intent( inout ) :: AB
    REAL, DIMENSION( NONODS1, NONODS2 ), intent( in )    :: A
    REAL, DIMENSION( NONODS3, NONODS4 ), intent( in )    :: B
    ! Local
    INTEGER :: I, J, II
    !          IF(NONODS2.NE.NONODS3) STOP 8329

    Loop_I: DO I = 1, NONODS1

       Loop_J: DO J = 1, NONODS4
          AB( I, J ) = 0.0

          Loop_II: DO II = 1, NONODS2
             AB( I, J ) = AB( I, J ) + A( I, II ) * B( II, J )
          END DO Loop_II

       END DO Loop_J

    END DO Loop_I

    RETURN
  END SUBROUTINE ABMATRIXMUL




  SUBROUTINE COLOR_GET_CMC_PHA( CV_NONODS, U_NONODS, NDIM, NPHASE, &
       NCOLC, FINDC, COLC, &
       INV_PIVIT_MAT,  &
       TOTELE, U_NLOC, U_NDGLN, &
       NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
       CMC, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
       C, CT, NDPSET )
    !use multiphase_1D_engine

    implicit none
    ! form pressure matrix CMC using a colouring approach
    INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, &
         TOTELE, U_NLOC, NCOLCT, NCOLCMC
    INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) ::FINDC
    INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
    REAL, DIMENSION( TOTELE, U_NLOC * NPHASE * NDIM, U_NLOC * NPHASE * NDIM ), intent( inout ) :: INV_PIVIT_MAT
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) ::  U_NDGLN
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
    INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
    REAL, DIMENSION( CV_NONODS ), intent( in ) :: DIAG_SCALE_PRES
    REAL, DIMENSION( NCOLCMC ), intent( inout ) :: CMC
    REAL, DIMENSION( NCOLCMC ), intent( in ) :: MASS_MN_PRES
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCMC
    INTEGER, DIMENSION( NCOLCMC ), intent( in ) :: COLCMC
    REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( in ) :: C 
    REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( inout ) :: CT
    INTEGER, intent( in ) :: NDPSET
    ! Local variables
    INTEGER, PARAMETER :: MX_NCOLOR = 1000
    REAL, PARAMETER :: INFINY = 1.0E+10
    LOGICAL :: DONE
    REAL, DIMENSION( :), allocatable :: NEED_COLOR, COLOR_VEC, CMC_COLOR_VEC
    REAL, DIMENSION( :), allocatable :: CDP, DU, DV, DW, DU_LONG
    INTEGER :: NCOLOR, CV_NOD, CV_JNOD, COUNT, COUNT2, IDIM, IPHASE, CV_COLJ, U_JNOD, CV_JNOD2
    REAL :: SUM

    ALLOCATE( NEED_COLOR( CV_NONODS ))
    ALLOCATE( COLOR_VEC( CV_NONODS ))
    ALLOCATE( CDP( U_NONODS * NDIM * NPHASE ))
    ALLOCATE( DU_LONG( U_NONODS * NDIM * NPHASE ))
    ALLOCATE( DU( U_NONODS * NPHASE ))
    ALLOCATE( DV( U_NONODS * NPHASE ))
    ALLOCATE( DW( U_NONODS * NPHASE ))
    ALLOCATE( CMC_COLOR_VEC( CV_NONODS )) 


    if(.false.) then
       DO CV_NOD = 1, CV_NONODS
          DO COUNT = FINDCT( CV_NOD ), FINDCT( CV_NOD + 1 ) - 1
             U_JNOD = COLCT( COUNT )
             ! Find CT for row CV_NOD and coln U_JNOD
             DO COUNT2=FINDC(U_JNOD),FINDC(U_JNOD+1)-1
                CV_COLJ=COLC(COUNT2)
                IF(CV_COLJ == CV_NOD) THEN
                   DO IDIM=1,NDIM
                      DO IPHASE=1,NPHASE
                         CT(COUNT2+(IDIM-1)*NCOLCT+NCOLCT*NDIM*(IPHASE-1))=  &
                              C(COUNT+(IDIM-1)*NCOLC+NCOLC*NDIM*(IPHASE-1))
                      END DO
                   END DO
                ENDIF
             END DO
          END DO
       END DO

    endif


    CMC = 0.0

    NEED_COLOR = 1.0
    NCOLOR = 0
    DONE = .FALSE. 

    Loop_while: DO WHILE (.NOT.DONE) 

       NCOLOR=NCOLOR+1  ! Determine what nodes can be coloured with the new color       
       COLOR_VEC = 0.0

       Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
          IF( ABS(NEED_COLOR( CV_NOD )) > 0.5 ) THEN 

             SUM = 0.0
             ! use a distance-2 colouring...
             Loop_Row: DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                CV_JNOD = COLCMC( COUNT )
                Loop_Row2: DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                   CV_JNOD2 = COLCMC( COUNT2 )
                   SUM = SUM + COLOR_VEC( CV_JNOD2 ) * NEED_COLOR( CV_JNOD2 )
                END DO Loop_Row2
             END DO Loop_Row

             IF( ABS(SUM) < 0.5 ) THEN
                COLOR_VEC( CV_NOD ) = 1.0
                !!                NEED_COLOR( CV_NOD ) = 0.0
             ENDIF

          ENDIF
          !   ewrite(3,*)'NCOLOR,CV_NOD,COLOR_VEC( CV_NOD ),NEED_COLOR( CV_NOD ):', &
          !               NCOLOR,CV_NOD,COLOR_VEC( CV_NOD ),NEED_COLOR( CV_NOD )
       END DO Loop_CVNOD

       !       COLOR_VEC=0
       !       COLOR_VEC(NCOLOR)=1. 

       NEED_COLOR = NEED_COLOR -  COLOR_VEC

       !       EWRITE(3,*)'NEED_COLOR:',NEED_COLOR
       !       EWRITE(3,*)'COLOR_VEC:',COLOR_VEC
       !       STOP 7756

       !       COLOR_VEC=0
       !       COLOR_VEC(2)=1.

       !       stop 393

       ! CMC_COLOR_VEC=CMC*COLOR_VEC
       ! CDP=C*COLOR_VEC
       CALL C_MULT( CDP, COLOR_VEC, CV_NONODS, U_NONODS, NDIM, NPHASE, &
            C, NCOLC, FINDC, COLC )
       !!DU_LONG = BLOCK_MAT * CDP
       if(.true.) then
          !         do ele=1,totele
          !     ewrite(3,*)'ele=',ele
          !     do i=1,U_nloc*nphase
          !         ewrite(3,*)INV_PIVIT_MAT(ele,i,:)
          !     end do
          !  end do
          !       INV_PIVIT_MAT=0.0
          !       do i=1,u_nloc*ndim*nphase
          !         INV_PIVIT_MAT(:,i,i)=1.0
          !       end do
          CALL PHA_BLOCK_MAT_VEC( DU_LONG, INV_PIVIT_MAT, CDP, U_NONODS, NDIM, NPHASE, &
               TOTELE, U_NLOC, U_NDGLN )
          ! NB. P_RHS=CT*U + CV_RHS 

          CALL ULONG_2_UVW( DU, DV, DW, DU_LONG, U_NONODS, NDIM, NPHASE )

       else
          DO IPHASE=1,NPHASE
             IDIM=1
             DU(1+(IPHASE-1)*U_NONODS:U_NONODS+(IPHASE-1)*U_NONODS)  &
                  =CDP( 1+U_NONODS*(IDIM-1)+U_NONODS*NDIM*(IPHASE-1): &
                  U_NONODS+U_NONODS*(IDIM-1)+U_NONODS*NDIM*(IPHASE-1) )
          END DO
       endif
       CALL CT_MULT( CMC_COLOR_VEC, DU, DV, DW, CV_NONODS, U_NONODS, NDIM, NPHASE, &
            CT, NCOLCT, FINDCT, COLCT )

       ! Matrix vector involving the mass diagonal term
       DO CV_NOD = 1, CV_NONODS
          DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
             CV_JNOD = COLCMC( COUNT )
             CMC_COLOR_VEC(CV_NOD) = CMC_COLOR_VEC(CV_NOD) &
                  +  DIAG_SCALE_PRES(CV_NOD)*MASS_MN_PRES(COUNT)*COLOR_VEC(CV_JNOD)
          END DO
       END DO
       IF(NDPSET /= 0) THEN
          CV_NOD=NDPSET
          CMC_COLOR_VEC(CV_NOD) = CMC_COLOR_VEC(CV_NOD) &
               +  INFINY*COLOR_VEC(CV_JNOD)
       ENDIF

       !Put into matrix CMC
       DO CV_NOD = 1, CV_NONODS 
          !         ewrite(3,*)CV_NOD,CMC_COLOR_VEC( CV_NOD )

          DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
             CV_JNOD = COLCMC( COUNT )
             CMC( COUNT ) = CMC( COUNT ) + CMC_COLOR_VEC( CV_NOD ) * COLOR_VEC( CV_JNOD )
          END DO

       END DO
       !       stop 383

       SUM = 0.0
       DO CV_NOD = 1, CV_NONODS
          SUM = SUM + NEED_COLOR( CV_NOD )
       END DO
       DONE = ( ABS(SUM) < 0.5 )

       ewrite(3,*)'************ sum,done,NCOLOR=',sum,done,NCOLOR

    END DO Loop_while


    if(.false.) then
       DO CV_NOD = 1, CV_NONODS

          DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
             CV_JNOD = COLCMC( COUNT )

          END DO

       END DO
       !       cmc(FINDCMC( CV_NONODS + 1 ) - 1)=50000.
    endif

    DEALLOCATE( NEED_COLOR )
    DEALLOCATE( COLOR_VEC )
    DEALLOCATE( CDP )
    DEALLOCATE( DU_LONG)
    DEALLOCATE( DU)
    DEALLOCATE( DV)
    DEALLOCATE( DW)
    DEALLOCATE( CMC_COLOR_VEC )

    RETURN
  END SUBROUTINE COLOR_GET_CMC_PHA




  SUBROUTINE PHA_BLOCK_INV( INV_PIVIT_MAT, PIVIT_MAT, TOTELE, NBLOCK )
    implicit none
    INTEGER, intent( in ) :: TOTELE, NBLOCK
    REAL, DIMENSION( TOTELE, NBLOCK, NBLOCK), intent( inout ) ::  INV_PIVIT_MAT 
    REAL, DIMENSION( TOTELE, NBLOCK, NBLOCK), intent( in ) ::  PIVIT_MAT
    ! Local variables
    REAL, DIMENSION( :, :), allocatable :: DMAT, DMATINV
    INTEGER :: ELE

    ALLOCATE( DMAT( NBLOCK, NBLOCK ))
    ALLOCATE( DMATINV( NBLOCK, NBLOCK ))
    !    

    DO ELE = 1, TOTELE

       DMAT( 1 : NBLOCK, 1 : NBLOCK ) = PIVIT_MAT( ELE, 1 : NBLOCK, 1 : NBLOCK )

       CALL MATDMATINV( DMAT, DMATINV, NBLOCK )

       INV_PIVIT_MAT( ELE, 1 : NBLOCK, 1 : NBLOCK ) = DMATINV( 1 : NBLOCK, 1 : NBLOCK )

    END DO

    DEALLOCATE( DMAT )
    DEALLOCATE( DMATINV )

    RETURN 

  END SUBROUTINE PHA_BLOCK_INV



  SUBROUTINE PHA_BLOCK_MAT_VEC( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, &
       TOTELE, U_NLOC, U_NDGLN ) 
    implicit none
    ! U = BLOCK_MAT * CDP
    INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC
    INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) ::  U_NDGLN
    REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( inout ) :: U
    REAL, DIMENSION( TOTELE, U_NLOC * NDIM * NPHASE, U_NLOC * NDIM * NPHASE ), intent( in ) :: BLOCK_MAT
    REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( in ) :: CDP
    ! Local 
    INTEGER :: ELE, U_ILOC, U_INOD, IDIM, IPHASE, I, U_JLOC, U_JNOD, JDIM, JPHASE, J, II, JJ

    U = 0.0 

    Loop_Elements: DO ELE = 1, TOTELE

       Loop_VelocNodsI: DO U_ILOC = 1, U_NLOC
          U_INOD = U_NDGLN(( ELE - 1 ) *U_NLOC + U_ILOC )

          Loop_DimensionsI: DO IDIM = 1, NDIM

             Loop_PhasesI: DO IPHASE = 1, NPHASE
                I = U_INOD + ( IDIM - 1 ) * U_NONODS +( IPHASE - 1 ) * NDIM * U_NONODS

                Loop_VelocNodsJ: DO U_JLOC = 1, U_NLOC
                   U_JNOD = U_NDGLN(( ELE - 1 ) * U_NLOC + U_JLOC )

                   Loop_DimensionsJ: DO JDIM = 1, NDIM

                      Loop_PhasesJ: DO JPHASE = 1, NPHASE

                         J = U_JNOD + ( JDIM - 1 ) * U_NONODS + ( JPHASE - 1 ) * NDIM * U_NONODS
                         II = U_ILOC + ( IDIM - 1 ) * U_NLOC + ( IPHASE - 1 ) * NDIM * U_NLOC
                         JJ = U_JLOC + ( JDIM - 1 ) * U_NLOC + ( JPHASE - 1 ) * NDIM * U_NLOC
                         U( I ) = U( I ) + BLOCK_MAT( ELE, II, JJ ) * CDP( J )

                      END DO Loop_PhasesJ

                   END DO Loop_DimensionsJ

                END DO Loop_VelocNodsJ

             END DO Loop_PhasesI

          END DO Loop_DimensionsI

       END DO Loop_VelocNodsI

    END DO Loop_Elements

    RETURN

  END SUBROUTINE PHA_BLOCK_MAT_VEC




  SUBROUTINE CT_MULT( CV_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
       CT, NCOLCT, FINDCT, COLCT ) 
    ! CV_RHS=CT*U
    implicit none
    INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT
    REAL, DIMENSION( CV_NONODS ), intent( inout) :: CV_RHS
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W  
    INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
    INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
    REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( in ) :: CT

    ! Local variables
    INTEGER :: CV_INOD, COUNT, U_JNOD, IPHASE, J

    CV_RHS = 0.0

    DO CV_INOD = 1, CV_NONODS

       DO COUNT = FINDCT( CV_INOD ), FINDCT( CV_INOD + 1 ) - 1, 1
          U_JNOD = COLCT( COUNT )

          DO IPHASE = 1, NPHASE
             J = U_JNOD + ( IPHASE - 1 ) * U_NONODS

             CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( COUNT + NCOLCT * NDIM * ( IPHASE - 1 )) * U( J )
             IF( NDIM >= 2 ) CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( COUNT + &
                  NCOLCT + NCOLCT     * NDIM * ( IPHASE - 1 )) * V( J )
             IF( NDIM >= 3 ) CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( COUNT + &
                  2 * NCOLCT + NCOLCT * NDIM * ( IPHASE-  1 )) * W( J )

          END DO

       END DO

    END DO

    RETURN

  END SUBROUTINE CT_MULT




  SUBROUTINE C_MULT( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, &
       C, NCOLC, FINDC, COLC ) 
    implicit none
    ! CDP=C*DP
    INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC
    REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( inout ) :: CDP
    REAL, DIMENSION( CV_NONODS ), intent( in )  :: DP
    REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( in ) :: C
    INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) ::FINDC
    INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
    ! Local variables
    INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA

    CDP = 0.0

    Loop_VelNodes: DO U_INOD = 1, U_NONODS

       Loop_Crow: DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1, 1
          P_JNOD = COLC( COUNT )

          Loop_Phase: DO IPHASE = 1, NPHASE
             Loop_Dim: DO IDIM = 1, NDIM
                COUNT_DIM_PHA = COUNT + NCOLC*(IDIM-1) + NCOLC*NDIM*(IPHASE-1)
                I1 = U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS
                CDP( I1 ) = CDP( I1 ) + C( COUNT_DIM_PHA ) * DP( P_JNOD )
             END DO Loop_Dim
          END DO Loop_Phase

       END DO Loop_Crow

    END DO Loop_VelNodes

    RETURN

  END SUBROUTINE C_MULT




  SUBROUTINE ULONG_2_UVW( U, V, W, UP, U_NONODS, NDIM, NPHASE)
    implicit none
    INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE
    REAL, DIMENSION( U_NONODS * NPHASE ), intent( inout ) :: U,V,W
    REAL, DIMENSION( U_NONODS * NDIM * NPHASE ), intent( in ) :: UP
    ! local variables...
    INTEGER :: IPHASE
    DO IPHASE = 1, NPHASE
       U( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
            UP( 1 + ( IPHASE - 1 ) * NDIM * U_NONODS : U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
       IF( NDIM >= 2 ) &
            V( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
            UP( 1 + U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS : 2 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
       IF( NDIM >= 3 ) &
            W( 1 + ( IPHASE - 1 ) * U_NONODS : U_NONODS + ( IPHASE - 1 ) * U_NONODS ) = &
            UP( 1 + 2 * U_NONODS + ( IPHASE - 1) * NDIM * U_NONODS : 3 * U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS )
    END DO
    RETURN
  END SUBROUTINE ULONG_2_UVW






  SUBROUTINE INITIAL_INT_ARRAY(N, VEC, IVALUE )
    implicit none

    INTEGER, intent( in ) :: N
    INTEGER, DIMENSION( N ), intent( inout ) :: VEC
    INTEGER, intent( in ) :: IVALUE

    ! Local
    INTEGER :: I

    DO I = 1, N
       VEC( I ) = IVALUE
    END DO

    RETURN
  END SUBROUTINE INITIAL_INT_ARRAY



  SUBROUTINE INITIAL_REAL_ARRAY(N, VEC, RVALUE )
    implicit none

    INTEGER, intent( in ) :: N
    REAL, DIMENSION( N ), intent( inout ) :: VEC
    REAL, intent( in ) :: RVALUE

    ! Local
    INTEGER :: I

    DO I = 1, N
       VEC( I ) = RVALUE
    END DO

    RETURN
  END SUBROUTINE INITIAL_REAL_ARRAY


end module matrix_operations


