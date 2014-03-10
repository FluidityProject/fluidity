
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
    use spud
    implicit none

  contains

    SUBROUTINE MULMAT( VEC, CT, U, FREDOP, NONODS, NCOLCT, FINDCT, COLCT ) 
      IMPLICIT NONE
      ! perform VEC = C^T * P, matrix vector multiplication

      INTEGER, intent( in ) :: FREDOP, NONODS, NCOLCT
      REAL, DIMENSION( : ), intent( inout ) :: VEC
      REAL, DIMENSION( : ), intent( inout ) :: CT
      REAL, DIMENSION( : ), intent( in ) :: U
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
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
      REAL, DIMENSION( : ), intent( inout ) :: VEC
      REAL, DIMENSION( : ), intent( inout ) :: CT
      REAL, DIMENSION( : ), intent( in ) :: P
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
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
      REAL, DIMENSION( :, : ), intent( in ) :: DMAT
      REAL, DIMENSION( :, : ), intent( inout ) ::  DMATINV
      ! Local variables
    !  REAL, DIMENSION( NLOC , NLOC ) :: MAT, MAT2
    !  REAL, DIMENSION( NLOC ) :: X, B

      DMATINV = DMAT
      CALL MATINV( DMATINV, NLOC, NLOC)!, MAT, MAT2, X, B )

      RETURN
    END SUBROUTINE MATDMATINV
    !


    SUBROUTINE MATINV( A, N, NMAX)!, MAT, MAT2, X, B)
      ! This sub finds the inverse of the matrix A and puts it back in A. 
      ! MAT, MAT2, X and B are working vectors. 
      IMPLICIT NONE
      INTEGER, intent( in ) :: N, NMAX
      REAL, DIMENSION( :, : ), intent( inout ) ::  A
!      REAL, DIMENSION( N, N ), intent( inout ) :: MAT, MAT2
!      REAL, DIMENSION( N ), intent( inout ) :: X, B
      ! Local variables
      INTEGER :: ICOL, IM, JM
!      INTEGER , DIMENSION(N,N) :: IPIV
      INTEGER , DIMENSION(N ) :: IPIV
      REAL, DIMENSION( N, N ) :: MAT, MAT2
      REAL, DIMENSION( N ) :: X, B

      real, dimension(max(1,N*N)) :: WORK
      integer :: LWORK

      integer info

      interface
         subroutine dgetrf(N,M,A,NMAX,IPIV,INFO)
           integer :: N, M, NMAX, info
           real, dimension(NMAX,M) :: A
           integer, dimension(N) :: IPIV
         end subroutine dgetrf
      end interface

      interface
         subroutine dgetri(N,A,NMAX,IPIV,WORK,LWORK,INFO)
           integer :: N, NMAX, info, LWORK
           real, dimension(NMAX,N) :: A
           integer, dimension(N) :: IPIV
           real, dimension(LWORK) :: WORK
         end subroutine dgetri
      end interface

      LWORK=N*N

      call dgetrf(N,N,A,NMAX,IPIV,INFO)
      call dgetri(N,A,NMAX,IPIV,WORK,max(1,LWORK),INFO)
      if (info==0) then
         return
      else
         FLAbort("PIVIT Matrix block inversion failed")
      end if

    END SUBROUTINE MATINV

    SUBROUTINE MATMASSINV( MASINV, MMAT, NONODS, NLOC, TOTELE )
      IMPLICIT NONE
      INTEGER, intent( in ) :: NONODS, NLOC, TOTELE
      REAL, DIMENSION( :, : ), intent( inout ) ::  MASINV, MMAT
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
      REAL, DIMENSION( :, : ), intent( inout ) :: A
      REAL, DIMENSION( : ), intent( inout ) :: X
      REAL, DIMENSION( : ), intent( in ) :: B
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

    SUBROUTINE SMLINNGOT( A, X, B, NMX, N, IPIV, GOTDEC )
      IMPLICIT NONE
      INTEGER :: NMX, N
      REAL, DIMENSION( :, : ), intent( inout ) :: A
      real, DIMENSION( : ), intent( inout ) :: X ! inout as n might be < nmx
      real, DIMENSION( : ), intent( in ) ::  B
      LOGICAL, intent( in ) :: GOTDEC
      ! Local     
      REAL :: R
      INTEGER :: K, I, J
      INTEGER , DIMENSION(NMX,NMX), intent(inout) :: IPIV
      INTEGER :: INFO

      real, DIMENSION( NMX,1 ) :: Bloc

      external dgetrf, dgetrs

      ! IF GOTDEC then assume we have already got the LU decomposition in A
      ! Form X = A^{-1} B ;  Useful subroutine for inverse
      ! This sub overwrites the matrix A. 

     Cond_LUDecomp: IF( .NOT. GOTDEC ) THEN

         call dgetrf(NMX,NMX,A,NMX,IPIV,INFO)

      ENDIF Cond_LUDecomp

      Bloc(:,1)=B
      call dgetrs('N',NMX,1,A, NMX,IPIV, Bloc, NMX, INFO)
      X=Bloc(:,1)

      RETURN
    END SUBROUTINE SMLINNGOT
    !     


    SUBROUTINE ABMATRIXMUL( AB, A, NONODS1, NONODS2, &
         B, NONODS3, NONODS4 )
      !
      ! Perform matrix matrix multiplication: AB = A * B
      IMPLICIT NONE
      INTEGER, intent( in ) :: NONODS1, NONODS2, NONODS3, NONODS4
      REAL, DIMENSION( :, : ), intent( inout ) :: AB
      REAL, DIMENSION( :, : ), intent( in )    :: A
      REAL, DIMENSION( :, : ), intent( in )    :: B
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
         CMC, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
         C, CT )
      !use multiphase_1D_engine

      implicit none
      ! form pressure matrix CMC using a colouring approach
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, &
           TOTELE, U_NLOC, NCOLCT, NCOLCMC, IGOT_CMC_PRECON
      INTEGER, DIMENSION( : ), intent( in ) ::FINDC
      INTEGER, DIMENSION( : ), intent( in ) :: COLC
      REAL, DIMENSION( :, :, : ), intent( in ) :: INV_PIVIT_MAT
      INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      REAL, DIMENSION( : ), intent( in ) :: DIAG_SCALE_PRES
      REAL, DIMENSION( : ), intent( inout ) :: CMC
      REAL, DIMENSION( : ), intent( inout ) :: CMC_PRECON
      REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      REAL, DIMENSION( : ), intent( in ) :: C
      REAL, DIMENSION( : ), intent( inout ) :: CT

      ! Local variables
      LOGICAL, PARAMETER :: FAST = .true.
      REAL, PARAMETER :: INFINY = 1.0E+10
      LOGICAL :: UNDONE, LCOL
      logical, DIMENSION( : ), allocatable :: NEED_COLOR
      logical, DIMENSION( : ), allocatable :: to_color
      logical, DIMENSION( :, : ), allocatable :: COLOR_VEC_MANY_LOGICAL
      REAL, DIMENSION( : ), allocatable :: COLOR_VEC, CMC_COLOR_VEC, CMC_COLOR_VEC2
      REAL, DIMENSION( : ), allocatable :: CDP, DU, DV, DW, DU_LONG
      REAL, DIMENSION( :, : ), allocatable :: COLOR_VEC_MANY, CMC_COLOR_VECC_MANY, CMC_COLOR_VEC2C_MANY
      REAL, DIMENSION( :, : ), allocatable :: CDPC_MANY, DUC_MANY, DVC_MANY, DWC_MANY, DU_LONGC_MANY
      INTEGER :: NCOLOR, CV_NOD, CV_JNOD, COUNT, COUNT2, IDIM, IPHASE, CV_COLJ, U_JNOD, CV_JNOD2
      INTEGER :: I, ELE,u_inod,u_nod
      integer, save :: ndpset=-1
      REAL :: RSUM


      if (ndpset<0) call get_option( '/material_phase[0]/scalar_field::Pressure/' // &
           'prognostic/reference_node', ndpset, default = 0 )


      IF ( FAST ) THEN
         ! Fast but memory intensive...
         CALL COLOR_GET_CMC_PHA_FAST( CV_NONODS, U_NONODS, NDIM, NPHASE, &
              NCOLC, FINDC, COLC, &
              INV_PIVIT_MAT,  &
              TOTELE, U_NLOC, U_NDGLN, &
              NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
              CMC, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
              C, CT, ndpset)
      ELSE
         ! Slow but memory efficient...
         CALL COLOR_GET_CMC_PHA_SLOW( CV_NONODS, U_NONODS, NDIM, NPHASE, &
              NCOLC, FINDC, COLC, &
              INV_PIVIT_MAT,  &
              TOTELE, U_NLOC, U_NDGLN, &
              NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
              CMC, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
              C, CT, ndpset)
      END IF

    END SUBROUTINE COLOR_GET_CMC_PHA





       SUBROUTINE COLOR_GET_CMC_PHA_FAST( CV_NONODS, U_NONODS, NDIM, NPHASE, &
            NCOLC, FINDC, COLC, &
            INV_PIVIT_MAT,  &
            TOTELE, U_NLOC, U_NDGLN, &
            NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
            CMC, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
            C, CT, ndpset )

         implicit none
         ! form pressure matrix CMC using a colouring approach
         INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, &
              TOTELE, U_NLOC, NCOLCT, NCOLCMC, IGOT_CMC_PRECON
         INTEGER, intent( in ) :: ndpset 
         INTEGER, DIMENSION( : ), intent( in ) ::FINDC
         INTEGER, DIMENSION( : ), intent( in ) :: COLC
         REAL, DIMENSION( :, :, : ), intent( in ) :: INV_PIVIT_MAT
         INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
         INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
         INTEGER, DIMENSION( : ), intent( in ) :: COLCT
         REAL, DIMENSION( : ), intent( in ) :: DIAG_SCALE_PRES
         REAL, DIMENSION( : ), intent( inout ) :: CMC
         REAL, DIMENSION( : ), intent( inout ) :: CMC_PRECON
         REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
         INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
         INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
         REAL, DIMENSION( : ), intent( in ) :: C
         REAL, DIMENSION( : ), intent( inout ) :: CT

         ! Local variables
         INTEGER, PARAMETER :: MX_NCOLOR = 1000
         REAL, PARAMETER :: INFINY = 1.0E+10
         LOGICAL, PARAMETER :: SAVED_CMC_COLOR = .FALSE.
         LOGICAL :: LCOL
         INTEGER, DIMENSION( : ), allocatable :: ICOLOR
         INTEGER, DIMENSION( : ), allocatable :: COLOR_IN_ROW, COLOR_IN_ROW2
         REAL, DIMENSION( : ), allocatable :: COLOR_VEC, CMC_COLOR_VEC, CMC_COLOR_VEC2
         REAL, DIMENSION( : ), allocatable :: CDP, DU, DV, DW, DU_LONG
         REAL, DIMENSION( :, : ), allocatable :: COLOR_VEC_MANY, CMC_COLOR_VECC_MANY, CMC_COLOR_VEC2C_MANY
         REAL, DIMENSION( :, : ), allocatable :: CDP_MANY, CDPC_MANY
         REAL, DIMENSION( :, : ), allocatable :: DU_MANY, DV_MANY, DW_MANY, DU_LONG_MANY 
         REAL, DIMENSION( :, : ), allocatable :: CMC_COLOR_VEC_MANY, CMC_COLOR_VEC2_MANY
         INTEGER :: NCOLOR, CV_NOD, CV_JNOD, COUNT, COUNT2, COUNT3, IDIM, IPHASE, CV_COLJ, U_JNOD, CV_JNOD2
         INTEGER :: MAX_COLOR_IN_ROW, ICHOOSE
         INTEGER :: KVEC
         INTEGER :: I, ELE,u_inod,u_nod
         REAL :: RSUM

         ALLOCATE( ICOLOR( CV_NONODS ))
         ICOLOR = 0

         ALLOCATE( COLOR_VEC( CV_NONODS ))
         ALLOCATE( CDP( U_NONODS * NDIM * NPHASE ))
         ALLOCATE( DU_LONG( U_NONODS * NDIM * NPHASE ))
         ALLOCATE( DU( U_NONODS * NPHASE ))
         ALLOCATE( DV( U_NONODS * NPHASE ))
         ALLOCATE( DW( U_NONODS * NPHASE ))
         ALLOCATE( CMC_COLOR_VEC( CV_NONODS )) 
         ALLOCATE( CMC_COLOR_VEC2( CV_NONODS )) 


         CMC = 0.0
         IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON = 0.0

         MAX_COLOR_IN_ROW = 0
         DO CV_NOD = 1, CV_NONODS
            MAX_COLOR_IN_ROW = MAX( MAX_COLOR_IN_ROW, FINDCMC( CV_NOD + 1 ) - FINDCMC( CV_NOD )  )
         END DO

         ALLOCATE( COLOR_IN_ROW( MAX_COLOR_IN_ROW**2 ) ) 
         ALLOCATE( COLOR_IN_ROW2( MAX_COLOR_IN_ROW**2 ) ) 

         IF ( SAVED_CMC_COLOR ) THEN
            !         NCOLOR=
            !         ICOLOR=
         ELSE

            NCOLOR=0
            Loop_CVNOD7: DO CV_NOD = 1, CV_NONODS

               ! Color this node CV_NOD
               COUNT3=0
               ! use a distance-2 colouring...
               Loop_Row7: DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                  CV_JNOD = COLCMC( COUNT )
                  Loop_Row2: DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                     CV_JNOD2 = COLCMC( COUNT2 )
                     IF(CV_NOD.NE.CV_JNOD2) THEN
                        IF(ICOLOR(CV_JNOD2).NE.0) THEN
                           COUNT3=COUNT3+1
                           COLOR_IN_ROW(COUNT3)= ICOLOR(CV_JNOD2)
                        ENDIF
                     ENDIF
                  END DO Loop_Row2
               END DO Loop_Row7


               ! Perform a bubble sort...
               CALL IBUBLE(COLOR_IN_ROW(1:COUNT3),COUNT3) 
               IF(COUNT3==0) THEN
                  COUNT2=0
               ELSE
                  COUNT2=1
                  I=1
                  COLOR_IN_ROW2(COUNT2)=COLOR_IN_ROW(I)
                  DO I=2,COUNT3
                     IF(COLOR_IN_ROW(I).NE.COLOR_IN_ROW(I-1)) THEN
                        COUNT2=COUNT2+1
                        COLOR_IN_ROW2(COUNT2)=COLOR_IN_ROW(I)
                     END IF
                  END DO
               END IF

               ICHOOSE=0
               DO I = 1, COUNT2
                  IF ( COLOR_IN_ROW2(I).NE.I ) THEN
                     ICHOOSE=I
                     EXIT
                  END IF
               END DO
               IF ( COUNT2==0 ) THEN
                  ICHOOSE=1
               ELSE
                  IF ( ICHOOSE==0 ) ICHOOSE = COUNT2 + 1
               END IF
               ICOLOR(CV_NOD)=ICHOOSE 
               NCOLOR = MAX( NCOLOR, ICHOOSE)
            END DO Loop_CVNOD7



            IF ( NCOLOR > MX_NCOLOR ) THEN
               PRINT *,'NOT ENOUGH COLOURS STOPPING - NEED TO MAKE MX_COLOR BIGGER'
               STOP 281
            END IF

         END IF ! ENDOF IF(SAVED_CMC_COLOR) THEN ELSE


         ALLOCATE( COLOR_VEC_MANY( NCOLOR, U_NONODS * NDIM * NPHASE ) ) 
         COLOR_VEC_MANY = 0.0

         Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
            COLOR_VEC_MANY( ICOLOR( CV_NOD ), CV_NOD ) = 1.0
         END DO Loop_CVNOD


         ALLOCATE( CDP_MANY( NCOLOR, NDIM*NPHASE*U_NONODS ) )

         CALL C_MULT_MANY( CDP_MANY, COLOR_VEC_MANY, CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLOR, &
              C, NCOLC, FINDC, COLC)
      
         ! DU_LONG = INV_PIVIT_MAT * CDP
         ALLOCATE( DU_LONG_MANY( NCOLOR, NDIM*NPHASE*U_NONODS ) ) 

         CALL PHA_BLOCK_MAT_VEC_MANY( DU_LONG_MANY, INV_PIVIT_MAT, CDP_MANY, U_NONODS, NDIM, NPHASE, NCOLOR, &
              TOTELE, U_NLOC, U_NDGLN )
  

         ! NB. P_RHS = CT * U + CV_RHS 
         ! DU_LONG = CDP

         ALLOCATE( DU_MANY( NCOLOR, NPHASE*U_NONODS ), DV_MANY( NCOLOR, NPHASE*U_NONODS ), DW_MANY( NCOLOR, NPHASE*U_NONODS ) ) 
         DO KVEC = 1, NCOLOR
            CALL ULONG_2_UVW( DU_MANY( KVEC, : ), DV_MANY( KVEC, : ), DW_MANY( KVEC, : ), DU_LONG_MANY( KVEC, : ), U_NONODS, NDIM, NPHASE )
         END DO

         ALLOCATE( CMC_COLOR_VEC_MANY( NCOLOR, CV_NONODS ) ) 

         CALL CT_MULT_MANY( CMC_COLOR_VEC_MANY, DU_MANY, DV_MANY, DW_MANY, CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLOR, &
              CT, NCOLCT, FINDCT, COLCT )

         IF ( IGOT_CMC_PRECON /= 0 ) THEN
            ALLOCATE( CMC_COLOR_VEC2_MANY( NCOLOR, CV_NONODS ) ) 
            CALL CT_MULT_WITH_C_MANY( CMC_COLOR_VEC2_MANY, DU_LONG_MANY, CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLOR, &
                 C, NCOLC, FINDC, COLC )
         END IF
         DEALLOCATE( DU_LONG_MANY )
         DEALLOCATE( DU_MANY, DV_MANY, DW_MANY ) 

         ! Matrix vector involving the mass diagonal term
         DO CV_NOD = 1, CV_NONODS
            RSUM = 0.0
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               CV_JNOD = COLCMC( COUNT )
               DO KVEC = 1, NCOLOR
                  CMC_COLOR_VEC_MANY( KVEC, CV_NOD ) = CMC_COLOR_VEC_MANY( KVEC, CV_NOD ) &
                       +  DIAG_SCALE_PRES( CV_NOD ) * MASS_MN_PRES( COUNT ) * COLOR_VEC_MANY( KVEC, CV_JNOD )
               END DO
               RSUM = RSUM + MASS_MN_PRES( COUNT )
            END DO
            IF ( IGOT_CMC_PRECON /= 0 ) THEN ! Use lumping of MASS_MN_PRES...
               DO KVEC = 1, NCOLOR
                  CMC_COLOR_VEC2_MANY( KVEC, CV_NOD ) = CMC_COLOR_VEC2_MANY( KVEC, CV_NOD ) &
                       +  DIAG_SCALE_PRES( CV_NOD ) * RSUM * COLOR_VEC_MANY( KVEC, CV_NOD )
               END DO
            END IF
         END DO

         !Put into matrix CMC

         DO CV_NOD = 1, CV_NONODS 
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               CV_JNOD = COLCMC( COUNT )
               DO KVEC = 1, NCOLOR
                  CMC( COUNT ) = CMC( COUNT ) + CMC_COLOR_VEC_MANY( KVEC, CV_NOD ) * COLOR_VEC_MANY( KVEC, CV_JNOD )
                  IF ( IGOT_CMC_PRECON /= 0 ) THEN 
                     CMC_PRECON( COUNT ) = CMC_PRECON( COUNT ) + CMC_COLOR_VEC2_MANY( KVEC, CV_NOD ) * COLOR_VEC_MANY( KVEC, CV_JNOD )
                  END IF
               END DO
            END DO
         END DO



         IF ( NDPSET /= 0 ) THEN
            CV_NOD = NDPSET
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               CV_JNOD = COLCMC( COUNT )
               IF ( CV_JNOD /= CV_NOD ) THEN
                  CMC( COUNT ) = 0.0 ! not the diagonal
                  IF ( IGOT_CMC_PRECON /= 0 ) CMC_PRECON( COUNT ) = 0.0
                  DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                     CV_JNOD2 = COLCMC( COUNT2 )
                     IF ( CV_JNOD2 == CV_NOD ) CMC( COUNT2 ) = 0.0 ! not the diagonal
                     IF ( IGOT_CMC_PRECON/=0 ) THEN
                        IF ( CV_JNOD2 == CV_NOD ) CMC_PRECON( COUNT2 ) = 0.0
                     END IF
                  END DO
               END IF
            END DO
         END IF

         !DEALLOCATE( NEED_COLOR )
         !DEALLOCATE( COLOR_VEC )
         DEALLOCATE( CDP )
         DEALLOCATE( DU_LONG )
         DEALLOCATE( DU, DV, DW )
         DEALLOCATE( CMC_COLOR_VEC ) 
         DEALLOCATE( CMC_COLOR_VEC2 )

         RETURN
       END SUBROUTINE COLOR_GET_CMC_PHA_FAST




    SUBROUTINE IBUBLE(LIST,NLIST)

      INTEGER NLIST,LIST(NLIST)
      INTEGER I,J,II
      do I=1,NLIST
         do J=2,NLIST
            IF(LIST(J-1).GT.LIST(J)) THEN
               !     SWOP
               II=LIST(J-1)
               LIST(J-1)=LIST(J)
               LIST(J)=II
            ENDIF
         END DO
      END DO
    END SUBROUTINE IBUBLE




       SUBROUTINE COLOR_GET_CMC_PHA_SLOW( CV_NONODS, U_NONODS, NDIM, NPHASE, &
         NCOLC, FINDC, COLC, &
         INV_PIVIT_MAT,  &
         TOTELE, U_NLOC, U_NDGLN, &
         NCOLCT, FINDCT, COLCT, DIAG_SCALE_PRES, &
         CMC, CMC_PRECON, IGOT_CMC_PRECON, NCOLCMC, FINDCMC, COLCMC, MASS_MN_PRES, &
         C, CT, ndpset )
      !use multiphase_1D_engine

      implicit none
      ! form pressure matrix CMC using a colouring approach
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, &
           TOTELE, U_NLOC, NCOLCT, NCOLCMC, IGOT_CMC_PRECON
      INTEGER, intent( in ) :: ndpset 
      INTEGER, DIMENSION( : ), intent( in ) ::FINDC
      INTEGER, DIMENSION( : ), intent( in ) :: COLC
      REAL, DIMENSION( :, :, : ), intent( in ) :: INV_PIVIT_MAT
      INTEGER, DIMENSION( : ), intent( in ) ::  U_NDGLN
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      REAL, DIMENSION( : ), intent( in ) :: DIAG_SCALE_PRES
      REAL, DIMENSION( : ), intent( inout ) :: CMC
      REAL, DIMENSION( : ), intent( inout ) :: CMC_PRECON
      REAL, DIMENSION( : ), intent( in ) :: MASS_MN_PRES
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCMC
      INTEGER, DIMENSION( : ), intent( in ) :: COLCMC
      REAL, DIMENSION( : ), intent( in ) :: C
      REAL, DIMENSION( : ), intent( inout ) :: CT

      ! Local variables
      INTEGER, PARAMETER :: MX_NCOLOR = 1000
      REAL, PARAMETER :: INFINY = 1.0E+10
      LOGICAL :: UNDONE, LCOL
      logical, DIMENSION( : ), allocatable :: NEED_COLOR
      logical, DIMENSION( CV_NONODS ) :: to_color
      REAL, DIMENSION( : ), allocatable :: COLOR_VEC, CMC_COLOR_VEC, CMC_COLOR_VEC2
      REAL, DIMENSION( : ), allocatable :: CDP, DU, DV, DW, DU_LONG
      INTEGER :: NCOLOR, CV_NOD, CV_JNOD, COUNT, COUNT2, IDIM, IPHASE, CV_COLJ, U_JNOD, CV_JNOD2
      INTEGER :: I, ELE,u_inod,u_nod
      REAL :: RSUM


      ALLOCATE( NEED_COLOR( CV_NONODS ))
      ALLOCATE( COLOR_VEC( CV_NONODS ))
      ALLOCATE( CDP( U_NONODS * NDIM * NPHASE ))
      ALLOCATE( DU_LONG( U_NONODS * NDIM * NPHASE ))
      ALLOCATE( DU( U_NONODS * NPHASE ))
      ALLOCATE( DV( U_NONODS * NPHASE ))
      ALLOCATE( DW( U_NONODS * NPHASE ))
      ALLOCATE( CMC_COLOR_VEC( CV_NONODS )) 
      ALLOCATE( CMC_COLOR_VEC2( CV_NONODS )) 


      CMC = 0.0
      IF(IGOT_CMC_PRECON.NE.0) CMC_PRECON = 0.0

      NEED_COLOR = .true.
      NCOLOR = 0
      UNDONE = .true. 

      Loop_while: DO WHILE (UNDONE) 

         NCOLOR=NCOLOR+1  ! Determine what nodes can be coloured with the new color       
         to_color=.false.

         Loop_CVNOD: DO CV_NOD = 1, CV_NONODS
            IF(NEED_COLOR( CV_NOD )) THEN 

               LCOL= .FALSE.
               ! use a distance-2 colouring...
               Loop_Row: DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
                  CV_JNOD = COLCMC( COUNT )
                  if(TO_COLOR( CV_JNOD )) then
                     LCOL=.true.
                     exit
                  end if
                  Loop_Row2: DO COUNT2 = FINDCMC( CV_JNOD ), FINDCMC( CV_JNOD + 1 ) - 1
                     if(TO_COLOR( COLCMC( COUNT2 ))) then
                        LCOL=.true.
                        exit
                     end if
                  END DO Loop_Row2
                  if (LCOL) exit 
               END DO Loop_Row

               IF(.not. LCOL) THEN
                  TO_COLOR( CV_NOD ) = .true.
               ENDIF

            ENDIF
         END DO Loop_CVNOD

         NEED_COLOR = NEED_COLOR .and. .not. TO_COLOR
         COLOR_VEC=merge(1.0,0.0,TO_COLOR)

         CALL C_MULT( CDP, COLOR_VEC, CV_NONODS, U_NONODS, NDIM, NPHASE, &
              C, NCOLC, FINDC, COLC)
         !!DU_LONG = BLOCK_MAT * CDP
         if(.true.) then
            !            do ele=1,totele
            !               ewrite(3,*)'ele=',ele
            !               do i=1,U_nloc*nphase*ndim
            !                  ewrite(3,*)i, INV_PIVIT_MAT(ele,i,:)
            !               end do
            !            end do
            !                   INV_PIVIT_MAT=0.0
            !                   do i=1,u_nloc*ndim*nphase
            !                     INV_PIVIT_MAT(:,i,i)=1.0
            !                   end do
            CALL PHA_BLOCK_MAT_VEC( DU_LONG, INV_PIVIT_MAT, CDP, U_NONODS, NDIM, NPHASE, &
                 TOTELE, U_NLOC, U_NDGLN )
            ! NB. P_RHS=CT*U + CV_RHS 
            !               DU_LONG=CDP

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
              CT, NCOLCT, FINDCT, COLCT)
      IF(IGOT_CMC_PRECON.NE.0) THEN
         CALL CT_MULT_WITH_C( CMC_COLOR_VEC2, DU_LONG, CV_NONODS, U_NONODS, NDIM, NPHASE, &
                C, NCOLC, FINDC, COLC )
!         CALL CT_MULT( CMC_COLOR_VEC2, DU, DV, DW, CV_NONODS, U_NONODS, NDIM, NPHASE, &
!              *****CT, NCOLCT, FINDCT, COLCT )
      ENDIF

         ! Matrix vector involving the mass diagonal term
         DO CV_NOD = 1, CV_NONODS
            !ewrite(3,*) 'cv_nod=',cv_nod
            RSUM=0.0
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               CV_JNOD = COLCMC( COUNT )
               !ewrite(3,*) 'CV_JNOD, diag:', CV_JNOD, &
               !     DIAG_SCALE_PRES(CV_NOD), MASS_MN_PRES(COUNT), COLOR_VEC(CV_JNOD)
               CMC_COLOR_VEC(CV_NOD) = CMC_COLOR_VEC(CV_NOD) &
                    +  DIAG_SCALE_PRES(CV_NOD) * MASS_MN_PRES(COUNT) * COLOR_VEC(CV_JNOD)
               RSUM=RSUM+MASS_MN_PRES(COUNT)
            END DO
            IF(IGOT_CMC_PRECON.NE.0) THEN ! Use lumping of MASS_MN_PRES...
               CMC_COLOR_VEC2(CV_NOD) = CMC_COLOR_VEC2(CV_NOD) &
                    +  DIAG_SCALE_PRES(CV_NOD) * RSUM * COLOR_VEC(CV_NOD)
            ENDIF
         END DO
         !IF(NDPSET /= 0) THEN
         !   CV_NOD=NDPSET
         !   CMC_COLOR_VEC(CV_NOD) = CMC_COLOR_VEC(CV_NOD) &
         !        +  INFINY * COLOR_VEC(CV_JNOD)
         !ENDIF

         !Put into matrix CMC

         DO CV_NOD = 1, CV_NONODS 
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               CV_JNOD = COLCMC( COUNT )
               CMC( COUNT ) = CMC( COUNT ) + CMC_COLOR_VEC( CV_NOD ) * COLOR_VEC( CV_JNOD )
               IF(IGOT_CMC_PRECON.NE.0) THEN 
                  CMC_PRECON( COUNT ) = CMC_PRECON( COUNT ) + CMC_COLOR_VEC2( CV_NOD ) * COLOR_VEC( CV_JNOD )
               ENDIF
               !ewrite(3,*) cv_nod, count, CV_JNOD, CMC_COLOR_VEC( CV_NOD ), COLOR_VEC( CV_JNOD )
            END DO
         END DO
         !stop 383

         UNDONE = any(NEED_COLOR)

         ewrite(3,*)'************ rsum,undone,NCOLOR=',rsum,undone,NCOLOR

      END DO Loop_while





      if (.false.) then

         do ncolor=1,-cv_nonods

            COLOR_VEC=0.0
            COLOR_VEC(ncolor)=1.0

            du_long=0.0
            cdp=0.0

            CALL C_MULT( CDP, COLOR_VEC, CV_NONODS, U_NONODS, NDIM, NPHASE, &
                 C, NCOLC, FINDC, COLC )

            if(.true.) then
               du_long=cdp

               CALL ULONG_2_UVW( DU, DV, DW, DU_LONG, U_NONODS, NDIM, NPHASE )

               CALL CT_MULT( CMC_COLOR_VEC, DU, DV, DW, CV_NONODS, U_NONODS, NDIM, NPHASE, &
                    CT, NCOLCT, FINDCT, COLCT )
            else
               ! CMC_COLOR_VEC=c^T CDP
               CMC_COLOR_VEC=0.0
               do u_nod=1,u_nonods
                  do count=findc(u_nod),findc(u_nod+1)-1
                     cv_jnod=colc(count)
                     CMC_COLOR_VEC(cv_jnod)=CMC_COLOR_VEC(cv_jnod)+c(count)*CDP(u_nod)
                  end do
               end do
            endif

         end do

      end if

      if(.false.) then
         DO CV_NOD = 1, CV_NONODS
            ewrite(3,*) 'cv_nod=',cv_nod
            DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
               CV_JNOD = COLCMC( COUNT )
               ewrite(3,*) 'CV_JNOD,cmc(count), diag:', CV_JNOD, cmc(count), &
                    DIAG_SCALE_PRES(CV_NOD), MASS_MN_PRES(COUNT)

            END DO
         END DO
         stop 740
      endif

      IF(NDPSET /= 0) THEN
         CV_NOD=NDPSET
         DO COUNT = FINDCMC( CV_NOD ), FINDCMC( CV_NOD + 1 ) - 1
            CV_JNOD = COLCMC( COUNT )
            if(cv_jnod/=cv_nod) then
               cmc(count)=0.0 ! not the diagonal
               IF(IGOT_CMC_PRECON.NE.0) cmc_PRECON(count)=0.0
               DO COUNT2 = FINDCMC( CV_jNOD ), FINDCMC( CV_jNOD + 1 ) - 1
                  CV_JNOD2 = COLCMC( COUNT2 )
                  if(cv_jnod2==cv_nod) cmc(count2)=0.0 ! not the diagonal
                  IF(IGOT_CMC_PRECON.NE.0) THEN
                     if(cv_jnod2==cv_nod) cmc_PRECON(count2)=0.0
                  ENDIF
               END DO
            endif
         END DO
      ENDIF

      DEALLOCATE( NEED_COLOR )
      DEALLOCATE( COLOR_VEC )
      DEALLOCATE( CDP )
      DEALLOCATE( DU_LONG )
      DEALLOCATE( DU,DV,DW )
      DEALLOCATE( CMC_COLOR_VEC ) 
      DEALLOCATE( CMC_COLOR_VEC2 )

      RETURN
    END SUBROUTINE COLOR_GET_CMC_PHA_SLOW




    SUBROUTINE PHA_BLOCK_INV( PIVIT_MAT, TOTELE, NBLOCK )
      implicit none
      INTEGER, intent( in ) :: TOTELE, NBLOCK
      REAL, DIMENSION( : , : , : ), intent( inout ), CONTIGUOUS ::  PIVIT_MAT
      ! Local variables
      INTEGER :: ELE


      DO ELE = 1, TOTELE

!         INV_PIVIT_MAT(:,:,ele)=PIVIT_MAT(:,:,ele)

!         CALL MATDMATINV( PIVIT_MAT(:,:,ele), INV_PIVIT_MAT(:,:,ele), NBLOCK )
         CALL MATINV( PIVIT_MAT(:,:,ele), NBLOCK, nblock )

      END DO


      RETURN 

    END SUBROUTINE PHA_BLOCK_INV





     SUBROUTINE PHA_BLOCK_MAT_VEC( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, &
         TOTELE, U_NLOC, U_NDGLN ) 
      implicit none
      ! U = BLOCK_MAT * CDP
      INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC
      INTEGER, DIMENSION( : ), intent( in ), target ::  U_NDGLN
      REAL, DIMENSION( : ), intent( inout ) :: U
      REAL, DIMENSION( :, :,: ), intent( in ), target :: BLOCK_MAT
      REAL, DIMENSION( : ), intent( in ) :: CDP
      ! Local 
      INTEGER :: ELE, U_ILOC, U_INOD, IDIM, IPHASE, I, U_JLOC, U_JNOD, JDIM, JPHASE, J, II, JJ

      integer, dimension(:), pointer :: U_NOD

      real, dimension(U_NLOC*NDIM*NPHASE) :: lcdp, lu
      integer, dimension(U_NLOC*NDIM*NPHASE) :: u_nodi
      integer :: N
      
      interface 
         subroutine dgemv(T,M,N,alpha,MAT,NMAX,X,Xinc,beta,Y,yinc)
           implicit none
           character(len=1) :: T
           integer :: m,n,nmax,xinc,yinc
           real ::  alpha, beta
           real, dimension(nmax,n) :: MAT
           real, dimension(N) :: X
           real, dimension(M) :: Y
         end subroutine dgemv
      end interface
           

      N=U_NLOC * NDIM * NPHASE

      Loop_Elements: DO ELE = 1, TOTELE

         U_NOD => U_NDGLN(( ELE - 1 ) * U_NLOC +1: ELE * U_NLOC)            


         Loop_PhasesJ: DO JPHASE = 1, NPHASE
            Loop_DimensionsJ: DO JDIM = 1, NDIM
                     
!               J = ( JDIM - 1 ) * U_NONODS + ( JPHASE - 1 ) * NDIM * U_NONODS
               JJ = ( JDIM - 1 ) * U_NLOC + ( JPHASE - 1 ) * NDIM * U_NLOC

               J=JDIM+(JPHASE-1)*NDIM

               lcdp([(J+(i-1)*ndim*nphase,i=1,u_NLOC)])=CDP(U_NOD+(J-1)*U_NONODS)
               U_NODI([(J+(i-1)*ndim*nphase,i=1,u_NLOC)])=U_NOD+(J-1)*U_NONODS
            end do Loop_DimensionsJ
         end do Loop_PhasesJ
                           

!         LU=U(U_NODI)
         call dgemv('N',N,N,1.0d0,BLOCK_MAT( : , : ,ele),N,LCDP,1,0.0d0,LU,1)
         U(U_NODI)=LU
                           
!         U( U_NODI) = U( U_NODI ) + matmul(LOC_BLOCK_MAT( : , : ), LCDP( : ))


      END DO Loop_Elements

      RETURN


    END SUBROUTINE PHA_BLOCK_MAT_VEC







    SUBROUTINE PHA_BLOCK_MAT_VEC_MANY2( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, NBLOCK, &
         TOTELE, U_NLOC, U_NDGLN ) 
      implicit none
      ! U = BLOCK_MAT * CDP
      INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC, NBLOCK
      INTEGER, DIMENSION( : ), intent( in ), target ::  U_NDGLN
      REAL, DIMENSION( :, : ), intent( inout ) :: U
      REAL, DIMENSION( :, :,: ), intent( in ), target :: BLOCK_MAT
      REAL, DIMENSION( :, : ), intent( in ) :: CDP
      ! Local 
      INTEGER :: ELE, U_ILOC, U_INOD, IDIM, IPHASE, I, U_JLOC, U_JNOD, JDIM, JPHASE, J, II, JJ, IVEC

      integer, dimension(:), pointer :: U_NOD

      real, dimension(U_NLOC*NDIM*NPHASE) :: lcdp, lu
      integer, dimension(U_NLOC*NDIM*NPHASE) :: u_nodi
      integer :: N
      
      interface 
         subroutine dgemv(T,M,N,alpha,MAT,NMAX,X,Xinc,beta,Y,yinc)
           implicit none
           character(len=1) :: T
           integer :: m,n,nmax,xinc,yinc
           real ::  alpha, beta
           real, dimension(nmax,n) :: MAT
           real, dimension(N) :: X
           real, dimension(M) :: Y
         end subroutine dgemv
      end interface

      N = U_NLOC * NDIM * NPHASE

      Loop_Elements: DO ELE = 1, TOTELE

         U_NOD => U_NDGLN( ( ELE - 1 ) * U_NLOC +1 : ELE * U_NLOC )            

         DO IVEC = 1, NBLOCK

            Loop_PhasesJ: DO JPHASE = 1, NPHASE
               Loop_DimensionsJ: DO JDIM = 1, NDIM

                  J = JDIM + ( JPHASE - 1 ) * NDIM
                  JJ = ( JDIM - 1 ) * U_NLOC + ( JPHASE - 1 ) * NDIM * U_NLOC

                  lcdp([(J+(i-1)*ndim*nphase, i = 1, u_NLOC)]) = CDP(IVEC, U_NOD+(J-1)*U_NONODS)
                  U_NODI([(J+(i-1)*ndim*nphase, i = 1, u_NLOC)]) = U_NOD+(J-1)*U_NONODS
               end do Loop_DimensionsJ
            end do Loop_PhasesJ

            call dgemv('N', N, N, 1.0d0, BLOCK_MAT( : , : , ELE ), N, LCDP, 1, 0.0d0, LU, 1 )
            U( IVEC, U_NODI ) = LU

         END DO

      END DO Loop_Elements

      RETURN


    END SUBROUTINE PHA_BLOCK_MAT_VEC_MANY2


 

    SUBROUTINE PHA_BLOCK_MAT_VEC_MANY( U, BLOCK_MAT, CDP, U_NONODS, NDIM, NPHASE, NBLOCK, &
         TOTELE, U_NLOC, U_NDGLN ) 
      implicit none
      ! U = BLOCK_MAT * CDP
      INTEGER, intent( in )  :: U_NONODS, NDIM, NPHASE, TOTELE, U_NLOC, NBLOCK
      INTEGER, DIMENSION( TOTELE * U_NLOC ), intent( in ) ::  U_NDGLN
      REAL, DIMENSION( NBLOCK, U_NONODS * NDIM * NPHASE ), intent( inout ) :: U
      REAL, DIMENSION( U_NLOC * NDIM * NPHASE, U_NLOC * NDIM * NPHASE, TOTELE ), intent( in ), target :: BLOCK_MAT
      REAL, DIMENSION( NBLOCK, U_NONODS * NDIM * NPHASE ), intent( in ) :: CDP
      ! Local 
      real, dimension(:,:), allocatable :: CDP_CONV, U_CONV

      INTEGER :: ELE, U_ILOC, U_INOD, IDIM, IPHASE, I, J, U_JLOC, U_JNOD, JDIM, JPHASE, II, JJ, IORIG

! convert 
      ALLOCATE( CDP_CONV(NBLOCK, NDIM*NPHASE*U_NONODS) )
      ALLOCATE( U_CONV(NBLOCK, NDIM*NPHASE*U_NONODS) )

!      PRINT *,'HERE1'
      
      DO ELE = 1, TOTELE
         DO U_ILOC = 1, U_NLOC
            U_INOD = U_NDGLN( (ELE-1)*U_NLOC+ U_ILOC)
!               PRINT *,'U_INOD=',U_INOD

            DO IPHASE = 1, NPHASE
               DO IDIM = 1, NDIM
                  I = IDIM+(IPHASE-1)*NDIM+(U_ILOC-1)*NDIM*NPHASE
                  II= I + (ELE-1)*NDIM*NPHASE*U_NLOC

                  IORIG=U_INOD + (IPHASE-1)*NDIM*U_NONODS + (IDIM-1)*U_NONODS
                  CDP_CONV(:,II)=CDP(:,IORIG)
               END DO
            END DO

         END DO
      END DO

!      PRINT *,'HERE2'

      U_CONV = 0.0 


      Loop_Elements: DO ELE = 1, TOTELE


                  Loop_VelocNodsJ: DO U_JLOC = 1, U_NLOC

                        Loop_PhasesJ: DO JPHASE = 1, NPHASE

                     Loop_DimensionsJ: DO JDIM = 1, NDIM
                              J = JDIM+(JPHASE-1)*NDIM+(U_JLOC-1)*NDIM*NPHASE
                              JJ= J + (ELE-1)*NDIM*NPHASE*U_NLOC


         Loop_VelocNodsI: DO U_ILOC = 1, U_NLOC

               Loop_PhasesI: DO IPHASE = 1, NPHASE

            Loop_DimensionsI: DO IDIM = 1, NDIM

                              I = IDIM+(IPHASE-1)*NDIM+(U_ILOC-1)*NDIM*NPHASE
                              II= I + (ELE-1)*NDIM*NPHASE*U_NLOC

                           U_CONV( :, II ) = U_CONV( :, II ) + BLOCK_MAT( I, J, ELE ) * CDP_CONV( :, JJ )

            END DO Loop_DimensionsI

               END DO Loop_PhasesI

         END DO Loop_VelocNodsI

                     END DO Loop_DimensionsJ

                        END DO Loop_PhasesJ

                  END DO Loop_VelocNodsJ

      END DO Loop_Elements

!      PRINT *,'HERE3'

! convert back
      
      DO ELE = 1, TOTELE
         DO U_ILOC = 1, U_NLOC
            U_INOD = U_NDGLN( (ELE-1)*U_NLOC+ U_ILOC)
!               PRINT *,'U_INOD=',U_INOD

            DO IPHASE = 1, NPHASE
               DO IDIM = 1, NDIM
                  I = IDIM+(IPHASE-1)*NDIM+(U_ILOC-1)*NDIM*NPHASE
                  II= I + (ELE-1)*NDIM*NPHASE*U_NLOC

                  IORIG=U_INOD + (IPHASE-1)*NDIM*U_NONODS + (IDIM-1)*U_NONODS
                  U(:,IORIG)=U_CONV(:,II)
               END DO
            END DO

         END DO
      END DO
!      PRINT *,'HERE4'


      RETURN

    END SUBROUTINE PHA_BLOCK_MAT_VEC_MANY



!!$    SUBROUTINE CT_MULT( CV_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
!!$         CT, NCOLCT, FINDCT, COLCT ) 
!!$      ! CV_RHS=CT*U
!!$      implicit none
!!$      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT
!!$      REAL, DIMENSION( CV_NONODS ), intent( inout) :: CV_RHS
!!$      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W  
!!$      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
!!$      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
!!$      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( in ) :: CT
!!$
!!$      ! Local variables
!!$      INTEGER :: CV_INOD, COUNT, U_JNOD, IPHASE, J
!!$
!!$      CV_RHS = 0.0
!!$
!!$      DO CV_INOD = 1, CV_NONODS
!!$
!!$         DO COUNT = FINDCT( CV_INOD ), FINDCT( CV_INOD + 1 ) - 1, 1
!!$            U_JNOD = COLCT( COUNT )
!!$
!!$            DO IPHASE = 1, NPHASE
!!$               J = U_JNOD + ( IPHASE - 1 ) * U_NONODS
!!$
!!$               CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( COUNT + NCOLCT * NDIM * ( IPHASE - 1 )) * U( J )
!!$               IF( NDIM >= 2 ) CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( COUNT + &
!!$                    NCOLCT + NCOLCT     * NDIM * ( IPHASE - 1 )) * V( J )
!!$               IF( NDIM >= 3 ) CV_RHS( CV_INOD ) = CV_RHS( CV_INOD ) + CT( COUNT + &
!!$                    2 * NCOLCT + NCOLCT * NDIM * ( IPHASE-  1 )) * W( J )
!!$            END DO
!!$
!!$         END DO
!!$
!!$      END DO
!!$
!!$      RETURN
!!$
!!$    END SUBROUTINE CT_MULT





    SUBROUTINE CT_MULT( CV_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
         CT, NCOLCT, FINDCT, COLCT ) 
      ! CV_RHS=CT*U
      implicit none
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT
      REAL, DIMENSION( : ), intent( out) :: CV_RHS
      REAL, DIMENSION( : ), intent( in ) :: U, V, W
      INTEGER, DIMENSION( : ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( : ), intent( in ) :: COLCT
      REAL, DIMENSION( : ), intent( in ) :: CT


      ! Local variables
      real, dimension( ncolct) :: CTU
      INTEGER :: CV_INOD, COUNT, U_JNOD, IPHASE, i,J,k
      integer, pointer :: countp

    ! when code is realigned, it will be possible to form the dot product
    !              sum(CT(1:ndim*nphase,CV_NOD,U(1:ndim*nphase,cv_nod) 
    ! directly, unrolling one loop and reducing indirection.


      DO CV_INOD = 1, CV_NONODS
         CTU( FINDCT( CV_INOD ) :FINDCT( CV_INOD + 1 )-1)=CT( FINDCT( CV_INOD ) :FINDCT( CV_INOD + 1 ) -1)&
              *u(colct( FINDCT( CV_INOD ) :FINDCT( CV_INOD + 1 )-1))
      end do
      if (ndim>=2) CTU=CTU+CT(ncolct+1:2*ncolct)*V(colct)
      if (ndim>=3) CTU=CTU+CT(2*ncolct+1:3*ncolct)*W(colct)

      DO IPHASE = 2, NPHASE
         CTU=CTU+CT((iphase-1)*ncolct*ndim+1:(iphase-1)*ncolct*ndim+ncolct)*U(colct + ( IPHASE - 1 ) * U_NONODS)
         if (ndim>=2) &
              CTU=CTU+CT((iphase-1)*ncolct*ndim+ncolct+1:(iphase-1)*ncolct*ndim+2*ncolct)&
              *V(colct + ( IPHASE - 1 ) * U_NONODS)
         if (ndim>=3) CTU=CTU+CT((iphase-1)*ncolct*ndim+2*ncolct+1:(iphase-1)*ncolct*ndim+3*ncolct)*W(colct + ( IPHASE - 1 ) * U_NONODS)
      end DO

      DO CV_INOD = 1, CV_NONODS
         CV_RHS( CV_INOD ) = sum(CTU( FINDCT( CV_INOD ) :FINDCT( CV_INOD + 1 ) - 1))
      END DO

      RETURN

    END SUBROUTINE CT_MULT




    SUBROUTINE CT_MULT_MANY( CV_RHS, U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, NBLOCK, &
         CT, NCOLCT, FINDCT, COLCT ) 
      ! CV_RHS = CT * U
      implicit none
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT, NBLOCK
      REAL, DIMENSION( NBLOCK, CV_NONODS ), intent( inout) :: CV_RHS
      REAL, DIMENSION(NBLOCK,  U_NONODS * NPHASE ), intent( in ) :: U, V, W  
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( in ) :: CT

      ! Local variables
      INTEGER :: CV_INOD, COUNT, U_JNOD, IPHASE, J, IVEC

      CV_RHS = 0.0

      DO CV_INOD = 1, CV_NONODS

         DO COUNT = FINDCT( CV_INOD ), FINDCT( CV_INOD + 1 ) - 1
            U_JNOD = COLCT( COUNT )

            DO IPHASE = 1, NPHASE
               J = U_JNOD + ( IPHASE - 1 ) * U_NONODS

               DO IVEC = 1, NBLOCK
                  CV_RHS( IVEC, CV_INOD ) = CV_RHS( IVEC, CV_INOD ) + CT( COUNT + NCOLCT * NDIM * ( IPHASE - 1 ) ) * U( IVEC, J )
                  IF( NDIM >= 2 ) CV_RHS( IVEC, CV_INOD ) = CV_RHS( IVEC, CV_INOD ) &
                       + CT( COUNT + NCOLCT + NCOLCT * NDIM * ( IPHASE - 1 ) ) * V( IVEC, J )
                  IF( NDIM >= 3 ) CV_RHS( IVEC, CV_INOD ) = CV_RHS( IVEC, CV_INOD ) &
                       + CT( COUNT + 2 * NCOLCT + NCOLCT * NDIM * ( IPHASE-  1 ) ) * W( IVEC, J )
               END DO
            END DO

         END DO

      END DO

      RETURN

    END SUBROUTINE CT_MULT_MANY





    SUBROUTINE C_MULT( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, &
         C, NCOLC, FINDC, COLC ) 
      implicit none
      ! CDP=C*DP
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC
      REAL, DIMENSION( : ), intent( inout ) :: CDP
      REAL, DIMENSION( : ), intent( in )  :: DP
      REAL, DIMENSION( : ), intent( in ), target :: C
      INTEGER, DIMENSION( : ), intent( in ) ::FINDC
      INTEGER, DIMENSION( : ), intent( in ), target :: COLC
      ! Local variables
      INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA,j,dim_pha
      real, dimension(NCOLC) :: ldp

      !interface
      !   real function ddot( N, dx,incx,dy, incy)
      !     integer :: n, incx, incy
      !     real, dimension(N) :: dx,dy
      !   end function ddot
      !end interface

      LDP = DP( COLC )

      Loop_Phase: DO IPHASE = 0, NPHASE-1
         Loop_Dim: DO IDIM = 1, NDIM
            DIM_PHA = (IDIM-1) + NDIM*IPHASE
            Loop_VelNodes: DO U_INOD = 1, U_NONODS 
               CDP( U_INOD + DIM_PHA * U_NONODS ) = &
                    dot_product( &
                    C( FINDC( U_INOD ) + NCOLC * DIM_PHA : FINDC( U_INOD + 1 ) - 1 + NCOLC * DIM_PHA ), &
                    LDP( FINDC( U_INOD ) : FINDC( U_INOD + 1 ) - 1 ) )
            END DO Loop_VelNodes
         END DO Loop_Dim
      END DO Loop_Phase


      return
    end SUBROUTINE C_MULT




    SUBROUTINE C_MULT_MANY( CDP, DP, CV_NONODS, U_NONODS, NDIM, NPHASE, NBLOCK, &
         C, NCOLC, FINDC, COLC ) 
      implicit none
      ! CDP=C*DP
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, NBLOCK
      REAL, DIMENSION( NBLOCK, U_NONODS * NDIM * NPHASE ), intent( inout ) :: CDP
      REAL, DIMENSION( NBLOCK, CV_NONODS ), intent( in )  :: DP
      REAL, DIMENSION( NCOLC * NDIM * NPHASE ), intent( in ) :: C
      INTEGER, DIMENSION( U_NONODS + 1 ), intent( in ) ::FINDC
      INTEGER, DIMENSION( NCOLC ), intent( in ) :: COLC
      ! Local variables
      INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA

      CDP = 0.0

      Loop_VelNodes: DO U_INOD = 1, U_NONODS

         Loop_Crow: DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1
            P_JNOD = COLC( COUNT )

            Loop_Phase: DO IPHASE = 1, NPHASE
               Loop_Dim: DO IDIM = 1, NDIM
                  COUNT_DIM_PHA = COUNT + NCOLC*(IDIM-1) + NCOLC*NDIM*(IPHASE-1)
                  I1 = U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS
                  CDP( :, I1 ) = CDP( :, I1 ) + C( COUNT_DIM_PHA ) * DP( :, P_JNOD )
               END DO Loop_Dim
            END DO Loop_Phase

         END DO Loop_Crow

      END DO Loop_VelNodes

      RETURN

    END SUBROUTINE C_MULT_MANY


   


    SUBROUTINE CT_MULT_WITH_C( DP, U_LONG, CV_NONODS, U_NONODS, NDIM, NPHASE, &
         C, NCOLC, FINDC, COLC ) 
      implicit none
      ! DP = (C)^T U_LONG
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC
      REAL, DIMENSION( : ), intent( in ) :: U_LONG
      REAL, DIMENSION( : ), intent( inout )  :: DP
      REAL, DIMENSION( : ), intent( in ) :: C
      INTEGER, DIMENSION( : ), intent( in ) ::FINDC
      INTEGER, DIMENSION( : ), intent( in ) :: COLC
      ! Local variables
      INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA

!      CDP = 0.0
      DP = 0.0

      Loop_VelNodes: DO U_INOD = 1, U_NONODS

         Loop_Crow: DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1, 1
            P_JNOD = COLC( COUNT )

            Loop_Phase: DO IPHASE = 1, NPHASE
               Loop_Dim: DO IDIM = 1, NDIM
                  COUNT_DIM_PHA = COUNT + NCOLC*(IDIM-1) + NCOLC*NDIM*(IPHASE-1)
                  I1 = U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS
!                  CDP( I1 ) = CDP( I1 ) + C( COUNT_DIM_PHA ) * DP( P_JNOD )
              DP( P_JNOD ) = DP( P_JNOD ) + C( COUNT_DIM_PHA ) * U_LONG( I1 )
               END DO Loop_Dim
            END DO Loop_Phase

         END DO Loop_Crow

      END DO Loop_VelNodes

      RETURN

    END SUBROUTINE CT_MULT_WITH_C







    SUBROUTINE CT_MULT_WITH_C_MANY( DP, U_LONG, CV_NONODS, U_NONODS, NDIM, NPHASE, NBLOCK, &
         C, NCOLC, FINDC, COLC ) 
      implicit none
      ! DP = (C)^T U_LONG
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLC, NBLOCK
      REAL, DIMENSION( :, : ), intent( in ) :: U_LONG
      REAL, DIMENSION( :, : ), intent( inout )  :: DP
      REAL, DIMENSION( : ), intent( in ) :: C
      INTEGER, DIMENSION( : ), intent( in ) ::FINDC
      INTEGER, DIMENSION( : ), intent( in ) :: COLC
      ! Local variables
      INTEGER :: U_INOD, COUNT, P_JNOD, IPHASE, I1, IDIM, COUNT_DIM_PHA

      !CDP = 0.0
      DP = 0.0

      Loop_VelNodes: DO U_INOD = 1, U_NONODS

         Loop_Crow: DO COUNT = FINDC( U_INOD ), FINDC( U_INOD + 1 ) - 1
            P_JNOD = COLC( COUNT )

            Loop_Phase: DO IPHASE = 1, NPHASE
               Loop_Dim: DO IDIM = 1, NDIM
                  COUNT_DIM_PHA = COUNT + NCOLC*(IDIM-1) + NCOLC*NDIM*(IPHASE-1)
                  I1 = U_INOD + (IDIM-1)*U_NONODS + ( IPHASE - 1 ) * NDIM * U_NONODS
                  !CDP( I1 ) = CDP( I1 ) + C( COUNT_DIM_PHA ) * DP( P_JNOD )
                  DP( :, P_JNOD ) = DP( :, P_JNOD ) + C( COUNT_DIM_PHA ) * U_LONG( :, I1 )
               END DO Loop_Dim
            END DO Loop_Phase

         END DO Loop_Crow

      END DO Loop_VelNodes

      RETURN

    END SUBROUTINE CT_MULT_WITH_C_MANY







    SUBROUTINE ULONG_2_UVW( U, V, W, UP, U_NONODS, NDIM, NPHASE)
      implicit none
      INTEGER, intent( in ) :: U_NONODS, NDIM, NPHASE
      REAL, DIMENSION( : ), intent( inout ) :: U,V,W
      REAL, DIMENSION( : ), intent( in ) :: UP
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



    subroutine posinmat( posmat, globi, globj, &
         nonods, findrm, colm, ncolm )
      ! Find position in matrix POSMAT which has column GLOBJ
      implicit none
      integer, intent( inout ) :: posmat
      integer, intent( in ) :: globi, globj, nonods, ncolm
      integer, dimension( : ), intent( in ) :: findrm
      integer, dimension( : ), intent( in ) :: colm
      ! Local variables
      integer, parameter :: nmax = 1000000
      integer :: inum, lower, upper, count

      lower = findrm( globi )
      upper = findrm( globi + 1 ) - 1

      count = 1
      Loop_While: do while ( count <= nmax )
         inum = lower + ( upper - lower + 1 ) / 2

         if( globj >= colm( inum ) ) then
            lower = inum
         else
            upper = inum
         end if

         if( ( upper - lower ) <= 1 ) then
            if( globj == colm( lower ) ) then
               posmat = lower
            else  
               posmat = upper
            end if
            return
         end if

         count = count + 1

      end do Loop_While

      return
    end subroutine posinmat

    subroutine assemble_global_multiphase_csr(global_csr,&
         block_csr,dense_block_matrix,block_to_global,global_dense_block)

      real, dimension(:), intent(out)    ::  global_csr
      real, dimension(:), intent(in)     :: block_csr
      real, dimension(:,:,:), intent(in) :: dense_block_matrix
      integer, dimension(:), intent(in)  :: block_to_global
      integer, dimension(:,:), intent(in)  :: global_dense_block

      integer :: node, jphase, count, node_count, nphase


      ewrite(3,*), "In  assemble_global_multiphase_csr"

      ! copy the block_csr to global using the assigned map

      ewrite(3,*) 'size(global_csr), size(block_csr)', size(global_csr), size(block_csr)
      global_csr=0.0
      global_csr(block_to_global)=block_csr(:)

      ! now for the dense block
      node_count=size(dense_block_matrix,3)
      nphase=size(dense_block_matrix,2)

      do node=1,node_count
         do jphase=1,nphase
            global_csr(global_dense_block(jphase,node):&
                 global_dense_block(jphase,node)+nphase-1)=&
                 global_csr(global_dense_block(jphase,node):&
                 global_dense_block(jphase,node)+nphase-1)+dense_block_matrix(:,jphase,node)
         end do
      end do


      ewrite(3,*), "Leaving assemble_global_multiphase_csr"


    end subroutine assemble_global_multiphase_csr

  end module matrix_operations


