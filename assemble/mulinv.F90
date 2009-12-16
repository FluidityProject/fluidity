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

module mulinv_module

  use fldebug
  use sml

  implicit none

  private

  public :: mulinv

contains

      SUBROUTINE MULINV(VEC1,VEC2,VEC3,&
     &            ML1,ML2,ML3,&
     &            VECX,VECY,VECZ,NONODS,D3,&
! The following is for rotations only ROTAT=.TRUE.
     &         ROTAT,NNODRO,  NODROT,&
     &         NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &         DM1,DM2,DM3 )
!       this sub Forms   vec1=(R M_L R^T+D)^{-1} vecx
      INTEGER NONODS
      REAL VEC1(NONODS),VEC2(NONODS),VEC3(NONODS),&
     &            ML1(NONODS),ML2(NONODS),ML3(NONODS),&
     &            VECX(NONODS),VECY(NONODS),VECZ(NONODS)
!       NB The normal is pointing out of the domain. 
!       The rotation matrix in 3-D is R=  
!       T1X   T1Y   T1Z
!       T2X   T2Y   T2Z
!       NX    NY    NZ
!       The rotation matrix in 2-D is R=  
!       T1X   T1Y   
!       NX    NY    
!       
      INTEGER NNODRO
      REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
      REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
      REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
      INTEGER NODROT(NNODRO)
      REAL DM1(NONODS),DM2(NONODS),DM3(NONODS)
      LOGICAL ROTAT
!       
      INTEGER D3
      REAL &
     &        A11, A12, A13,&
     &        A21, A22, A23,&
     &        A31, A32, A33
!
      INTEGER IDIM, I, II, NOD
      
!       
!       Form simple inversion part
      IF(D3.EQ.0) THEN
!       For 2-D ...
         IDIM=2
!       LD3=.FALSE.
      do I=1,NONODS
            VEC1(I)=VECX(I)/(ML1(I)+DM1(I))
            VEC2(I)=VECY(I)/(ML2(I)+DM2(I))
         END DO
!       
      do II=1,NNODRO
            NOD=NODROT(II)
!       Form R*(M_L R^T)
            CALL SMLMA2(&
     &        A11, A12, &
     &        A21, A22, &
! =
     &       T1X (II),  T1Y(II) ,&
     &        NX(II) ,   NY(II),    &
! *
     &  T1X(II)*ML1(NOD),NX(II)*ML1(NOD),&
     &  T1Y(II)*ML2(NOD), NY(II)*ML2(NOD)          )
! add in D
        A11=A11+DM1(NOD)
        A22=A22+DM2(NOD)
!
! Form vec1 = A^{-1} vecx
! Useful subroutine for inverse
         CALL SMLIN2(&
     &        A11, A12, &
     &        A21, A22, &
     &      VEC1(NOD),VEC2(NOD),&
! =
     &      VECX(NOD),VECY(NOD)    )
      END DO
!
! IF(D3.EQ.0) ELSE
      ELSE
         IDIM=3
!       ewrite(3,*) 'vecz:',vecz
      do I=1,NONODS
            VEC1(I)=VECX(I)/(ML1(I)+DM1(I))
            VEC2(I)=VECY(I)/(ML2(I)+DM2(I))
            VEC3(I)=VECZ(I)/(ML3(I)+DM3(I))
         END DO
!       
      do II=1,NNODRO
            NOD=NODROT(II)
!       Form R *(M_L R^T)
            CALL SMLMA3(&
     &        A11, A12, A13,&
     &        A21, A22, A23,&
     &        A31, A32, A33,&
! =
     &        T1X(II),   T1Y(II),   T1Z(II),&
     &        T2X(II),   T2Y(II),   T2Z(II),&
     &         NX(II),    NY(II),    NZ(II),&
! *
     &  T1X(II)*ML1(NOD), T2X(II)*ML1(NOD), NX(II)*ML1(NOD),&
     &  T1Y(II)*ML2(NOD), T2Y(II)*ML2(NOD), NY(II)*ML2(NOD),&
     &  T1Z(II)*ML3(NOD), T2Z(II)*ML3(NOD), NZ(II)*ML3(NOD)     )
              IF(NOD.EQ.20) THEN
             ewrite(3,*)'A11, A12, A13:'
             ewrite(3,*)'A21, A22, A23:'
             ewrite(3,*)'A31, A32, A33:'
             ewrite(3,*) A11, A12, A13
             ewrite(3,*) A21, A22, A23
             ewrite(3,*) A31, A32, A33
             ewrite(3,*)'T1X(II)*ML1(),T2X(II)*ML1(N),NX(II)*ML1(N)'
             ewrite(3,*)'T1Y(II)*ML2(NOD),T2Y(II)*ML2(),NY(II)*ML2()'
             ewrite(3,*)'T1Z(II)*ML3(NOD),T2Z(II)*ML3(),NZ(II)*ML3()'
             ewrite(3,*)T1X(II)*ML1(NOD),T2X(II)*ML1(NOD),NX(II)*ML1(NOD)
             ewrite(3,*)T1Y(II)*ML2(NOD),T2Y(II)*ML2(NOD),NY(II)*ML2(NOD)
             ewrite(3,*)T1Z(II)*ML3(NOD),T2Z(II)*ML3(NOD),NZ(II)*ML3(NOD)
            endif
!       add in D
            A11=A11+DM1(NOD)
            A22=A22+DM2(NOD)
            A33=A33+DM3(NOD)
!       
!       Form vec1 = A^{-1} vecx
!       Useful subroutine for inverse
            CALL SMLIN3(&
     &        A11, A12, A13,&
     &        A21, A22, A23,&
     &        A31, A32, A33,&
     &      VEC1(NOD),VEC2(NOD),VEC3(NOD), &
! =
     &      VECX(NOD),VECY(NOD),VECZ(NOD)  )
              IF(NOD.EQ.20) THEN
             ewrite(3,*)'A11, A12, A13:'
             ewrite(3,*)'A21, A22, A23:'
             ewrite(3,*)'A31, A32, A33:'
             ewrite(3,*) A11, A12, A13
             ewrite(3,*) A21, A22, A23
             ewrite(3,*) A31, A32, A33
             ewrite(3,*)'VEC1(NOD),VEC2(NOD),VEC3(NOD):',&
     &            VEC1(NOD),VEC2(NOD),VEC3(NOD)
             ewrite(3,*)'VECX(NOD),VECY(NOD),VECZ(NOD):',&
     &            VECX(NOD),VECY(NOD),VECZ(NOD)
             ewrite(3,*)'T1X(II),T1Y(II),T1Z(II):'
             ewrite(3,*)'T2X(II),T2Y(II),T2Z(II):'
             ewrite(3,*)'NX(II), NY(II), NZ(II) :'
             ewrite(3,*) T1X(II),T1Y(II),T1Z(II)
             ewrite(3,*) T2X(II),T2Y(II),T2Z(II)
             ewrite(3,*) NX(II), NY(II), NZ(II)
             ewrite(3,*)'ML1(NOD), ML2(NOD), ML3(NOD) :',&
     &                 ML1(NOD), ML2(NOD), ML3(NOD)
             ewrite(3,*)'DM1(NOD), DM2(NOD), DM3(NOD) :',&
     &                 DM1(NOD), DM2(NOD), DM3(NOD)
             ewrite(3,*)'II,NOD,NNODRO:',II,NOD,NNODRO
              ENDIF
         END DO
!       
!       ELSE IF(D3.EQ.0)
      ENDIF

  end subroutine mulinv

end module mulinv_module
