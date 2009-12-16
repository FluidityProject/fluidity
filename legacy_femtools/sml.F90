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

module sml
contains

  SUBROUTINE SMLIN2(&
       AGI, BGI, &
       CGI, DGI, &
       VEC1,VEC2,&
                                ! =
       VECX,VECY)
    ! Form vec1 = A^{-1} vecx
    ! Useful subroutine for inverse
    REAL &
         AGI, BGI, &
         CGI, DGI, &
         VEC1,VEC2,&
                                ! =
         VECX,VECY,&
                                !
         A11, A12, &
         A21, A22, &
         DETJ
    !
    DETJ=AGI*DGI-BGI*CGI 
    A11= DGI/DETJ
    A12=-BGI /DETJ
    !
    A21=-CGI /DETJ
    A22= AGI /DETJ
    !
    VEC1=A11*VECX+A12*VECY
    VEC2=A21*VECX+A22*VECY
  END SUBROUTINE SMLIN2

  SUBROUTINE SMLIN3(&
       AGI, BGI, CGI,&
       DGI, EGI, FGI,&
       GGI, HGI, KGI,&
       VEC1,VEC2,VEC3,&
                                ! =
       VECX,VECY,VECZ  )
    use FLDebug
    IMPLICIT NONE
    !       Form vec1 = A^{-1} vecx
    !       Useful subroutine for inverse
    REAL &
         AGI, BGI, CGI,&
         DGI, EGI, FGI,&
         GGI, HGI, KGI,&
         VEC1,VEC2,VEC3,&
                                ! =
         VECX,VECY,VECZ
    REAL&
         A11, A12, A13,&
         A21, A22, A23,&
         A31, A32, A33,&
         DETJ

    DETJ=AGI*(EGI*KGI-FGI*HGI)&
         -BGI*(DGI*KGI-FGI*GGI)&
         +CGI*(DGI*HGI-EGI*GGI)

    IF(DETJ .EQ. 0.0) THEN
       ewrite(-1, *) '*** ERROR: Found zero determinant in smlin3()'
       ewrite(3, *) '*** Cannot invert matrix'
       ewrite(3, *) AGI, BGI, CGI
       ewrite(3, *) DGI, EGI, FGI
       ewrite(3, *) GGI, HGI, KGI
       stop 4398
       !       Subsequently let it throw a SIGFPE in the normal fashion
    END IF

    A11= (EGI*KGI-FGI*HGI) /DETJ
    A21=-(DGI*KGI-FGI*GGI) /DETJ
    A31= (DGI*HGI-EGI*GGI) /DETJ
    !       
    A12=-(BGI*KGI-CGI*HGI) /DETJ
    A22= (AGI*KGI-CGI*GGI) /DETJ
    A32=-(AGI*HGI-BGI*GGI) /DETJ
    !       
    A13= (BGI*FGI-CGI*EGI) /DETJ
    A23=-(AGI*FGI-CGI*DGI) /DETJ
    A33= (AGI*EGI-BGI*DGI) /DETJ
    !       
    VEC1=A11*VECX+A12*VECY+A13*VECZ
    VEC2=A21*VECX+A22*VECY+A23*VECZ
    VEC3=A31*VECX+A32*VECY+A33*VECZ

  END SUBROUTINE SMLIN3

  SUBROUTINE SMLIN4(&
       A11, A12, A13, A14,&
       A21, A22, A23, A24,&
       A31, A32, A33, A34,&
       A41, A42, A43, A44,&
       VEC1,VEC2,VEC3,VEC4,&
                                ! =
       VECX,VECY,VECZ,VECZZ  )
    IMPLICIT NONE
    !       Form vec1 = A^{-1} vecx
    !       Useful subroutine for inverse
    REAL &
         A11, A12, A13, A14,&
         A21, A22, A23, A24,&
         A31, A32, A33, A34,&
         A41, A42, A43, A44,&
         VEC1,VEC2,VEC3,VEC4,&
                                ! =
         VECX,VECY,VECZ, VECZZ
    REAL A(4,4),X(4),B(4),R
    INTEGER N,K,I,J
    !
    A(1,1)=A11
    A(1,2)=A12
    A(1,3)=A13
    A(1,4)=A14
    !
    A(2,1)=A21
    A(2,2)=A22
    A(2,3)=A23
    A(2,4)=A24
    !
    A(3,1)=A31
    A(3,2)=A32
    A(3,3)=A33
    A(3,4)=A34
    !
    A(4,1)=A41
    A(4,2)=A42
    A(4,3)=A43
    A(4,4)=A44
    !
    B(1)=VECX
    B(2)=VECY
    B(3)=VECZ
    B(4)=VECZZ
    !
    N=4
    do  K=1,N-1! Was loop 10
       do  I=K+1,N! Was loop 20
          A(I,K)=A(I,K)/A(K,K)
       end do ! Was loop 20
       do  J=K+1,N! Was loop 30
          do  I=K+1,N! Was loop 30
             A(I,J)=A(I,J) - A(I,K)*A(K,J)
          end do ! Was loop 30
       end do ! Was loop 30
    end do ! Was loop 10
    !
    ! Solve L x=b
    do  I=1,N! Was loop 40
       R=0.
       do  J=1,I-1! Was loop 50
          R=R+A(I,J)*X(J)
       end do ! Was loop 50
       X(I)=B(I)-R
    end do ! Was loop 40
    !
    ! Solve U x=y
    do  I=N,1,-1! Was loop 140
       R=0.
       do  J=I+1,N! Was loop 150
          R=R+A(I,J)*X(J)
       end do ! Was loop 150
       X(I)=(X(I)-R)/A(I,I)
    end do ! Was loop 140
    !
    VEC1=X(1)
    VEC2=X(2)
    VEC3=X(3)
    VEC4=X(4)
  END SUBROUTINE SMLIN4

  SUBROUTINE SMLINN(A,X,B,NMX,N)
    IMPLICIT NONE
    INTEGER NMX,N
    REAL A(NMX,NMX),X(NMX),B(NMX)
    REAL R
    INTEGER K,I,J
    !     Form X = A^{-1} B
    !     Useful subroutine for inverse
    !     This sub overwrites the matrix A. 
    do K=1,N-1
       do I=K+1,N
          A(I,K)=A(I,K)/A(K,K)
       END DO
       do J=K+1,N
          do I=K+1,N
             A(I,J)=A(I,J) - A(I,K)*A(K,J)
          END DO
       END DO
    END DO
    !     
    !     Solve L_1 x=b
    do I=1,N
       R=0.
       do J=1,I-1
          R=R+A(I,J)*X(J)
       END DO
       X(I)=B(I)-R
    END DO
    !     
    !     Solve U x=y
    do I=N,1,-1
       R=0.
       do J=I+1,N
          R=R+A(I,J)*X(J)
       END DO
       X(I)=(X(I)-R)/A(I,I)
    END DO
    RETURN
  END SUBROUTINE SMLINN
  !     
  !     

  SUBROUTINE SMLINNGOT(A,X,B,NMX,N,GOTDEC)
    IMPLICIT NONE
    INTEGER NMX,N
    REAL, intent(inout):: A(NMX,NMX)
    real, intent(inout):: X(NMX) ! inout as n might be < nmx
    real, intent(in)::  B(NMX)
    LOGICAL, intent(in):: GOTDEC
    !      
    REAL R
    INTEGER K,I,J
    !     IF GOTDEC then assume we have already got the LU decomposition in A
    !     Form X = A^{-1} B
    !     Useful subroutine for inverse
    !     This sub overwrites the matrix A. 
    IF(.NOT.GOTDEC) THEN
       do K=1,N-1
          do I=K+1,N
             A(I,K)=A(I,K)/A(K,K)
          END DO
          do J=K+1,N
             do I=K+1,N
                A(I,J)=A(I,J) - A(I,K)*A(K,J)
             END DO
          END DO
       END DO
    ENDIF
    !     
    !     Solve L_1 x=b
    do I=1,N
       R=0.
       do J=1,I-1
          R=R+A(I,J)*X(J)
       END DO
       X(I)=B(I)-R
    END DO
    !     
    !     Solve U x=y
    do I=N,1,-1
       R=0.
       do J=I+1,N
          R=R+A(I,J)*X(J)
       END DO
       X(I)=(X(I)-R)/A(I,I)
    END DO
    RETURN
  END SUBROUTINE SMLINNGOT
  !     
  !     
  !     
  !     
  SUBROUTINE SMLMA2(&
       R11, R12,     &
       R21, R22,      &
                                ! =
       A11, A12, &
       A21, A22, &
                                ! * 
       B11, B12, &
       B21, B22                         )
    ! ****************
    ! THIS SUB PERFORMS MATRIX MATRIX MULT OF ORDER 2
    REAL  R11, R12,     &
         R21, R22,      &
                                ! =
         A11, A12, &
         A21, A22, &
                                ! * 
         B11, B12, &
         B21, B22       
    !
    R11=A11*B11+A12*B21
    R12=A11*B12+A12*B22
    !
    R21=A21*B11+A22*B21
    R22=A21*B12+A22*B22
  END SUBROUTINE SMLMA2

  SUBROUTINE SMLMA3(&
       R11, R12, R13,        &
       R21, R22, R23,      &
       R31, R32, R33,        &
       ! =
       A11, A12, A13,&
       A21, A22, A23,&
       A31, A32, A33,&
       ! * 
       B11, B12, B13,&
       B21, B22, B23,&
       B31, B32, B33             )
    ! ****************
    ! THIS SUB PERFORMS MATRIX MATRIX MULT OF ORDER 3
    REAL  R11, R12, R13,        &
         R21, R22, R23,      &
         R31, R32, R33,        &
         ! =
         A11, A12, A13,&
         A21, A22, A23,&
         A31, A32, A33,&
         ! * 
         B11, B12, B13,&
         B21, B22, B23,&
         B31, B32, B33   
    !
    R11=A11*B11+A12*B21+A13*B31
    R12=A11*B12+A12*B22+A13*B32
    R13=A11*B13+A12*B23+A13*B33
    !
    R21=A21*B11+A22*B21+A23*B31
    R22=A21*B12+A22*B22+A23*B32
    R23=A21*B13+A22*B23+A23*B33
    !
    R31=A31*B11+A32*B21+A33*B31
    R32=A31*B12+A32*B22+A33*B32
    R33=A31*B13+A32*B23+A33*B33
  END SUBROUTINE SMLMA3

end module sml
