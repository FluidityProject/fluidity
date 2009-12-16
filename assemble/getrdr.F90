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

module getrdr_module

  use sml

  implicit none

  private

  public :: getrdr

contains

         SUBROUTINE GETRDR(   NONODS, &
     &         ROTAT,NNODRO,NRTDR,   NODROT,&
     &         NX,NY,NZ, T1X,T1Y,T1Z, T2X,T2Y,T2Z,  &
     &         DM1,DM2,DM3, RTDR)
! This sub forms the matrix R^T D R in which D is DM1 etc.
! the commponents of the normal direction to the boundary(NX,NY,NZ)
! NODROT contains the nodes to rotate. 
! NB The normal is pointing out of the domain. 
!  The rotation matrix in 3-D is R=  
!   T1X   T1Y   T1Z
!   T2X   T2Y   T2Z
!   NX    NY    NZ
!  The rotation matrix in 2-D is R=  
!   T1X   T1Y   
!   NX    NY    
! 
         INTEGER NRTDR
       INTEGER NONODS,NNODRO,NODROT(NNODRO)
       REAL DM1(NONODS),DM2(NONODS),DM3(NONODS)
       REAL NX(NNODRO),NY(NNODRO),NZ(NNODRO)
       REAL T1X(NNODRO),T1Y(NNODRO),T1Z(NNODRO)
       REAL T2X(NNODRO),T2Y(NNODRO),T2Z(NNODRO)
!
       LOGICAL ROTAT
       LOGICAL D3
       REAL RTDR(NRTDR)
       REAL &
     &        A11, A12, A13,&
     &        A21, A22, A23,&
     &        A31, A32, A33
       INTEGER IDIM,I,II,NOD
!
          
!          ewrite(3,*)'NONODS,NNODRO,NRTDR,ROTAT:',
!     &                  NONODS,NNODRO,NRTDR,ROTAT
!
         IF(ROTAT) THEN
!
         D3=.FALSE.
         IF(NRTDR.EQ.6*NONODS) D3=.TRUE.
!          ewrite(3,*)'D3=',D3
!
! Form (R^T D R) = RTDR
! NB in 2-D it has block numbering ...
! 1 *
! 2 3
! and in 3-D ...
! 1 * *
! 2 3 *
! 4 5 6       *-not stored
!
           IF(.NOT.D3) THEN
! For 2-D ...
            IDIM=2
!            LD3=.FALSE.
      do  I=1,3*NONODS! Was loop 19
            RTDR(I) =0.
      end do ! Was loop 19
      do  I=1,NONODS! Was loop 21
             RTDR(I)=DM1(I)
             RTDR(I+2*NONODS)=DM2(I)
      end do ! Was loop 21
!
      do  II=1,-NNODRO! Was loop 31
             NOD=NODROT(II)
! Form R^T*(D R)
           CALL SMLMA2(&
     &        A11, A12, &
     &        A21, A22, &
! =
     &       T1X (II),   NX(II) ,&
     &        T1Y(II) ,  NY(II),    &
! *
     &       T1X (II)*DM1(NOD),  T1Y(II)*DM1(NOD) ,&
     &        NX(II)*DM2(NOD) ,   NY(II)*DM2(NOD)    ) 
! set values of (R^T D R)
             RTDR(NOD)         =A11
             RTDR(NOD+NONODS  )=A21
             RTDR(NOD+2*NONODS)=A22
      end do ! Was loop 31
!
! IF(D3)
           ELSE
! For 3-D ...
            IDIM=3
!            LD3=.TRUE.
      do  I=1,6*NONODS! Was loop 191
            RTDR(I) =0.
      end do ! Was loop 191
      do  I=1,NONODS! Was loop 22
             RTDR(I)=DM1(I)
             RTDR(I+2*NONODS)=DM2(I)
             RTDR(I+5*NONODS)=DM3(I)
      end do ! Was loop 22
!
      do  II=1,-NNODRO! Was loop 30
!                 ewrite(3,*)'II:',II
             NOD=NODROT(II)
!                 ewrite(3,*)'NOD,II:',NOD,II
! Form R ^T*(D R)
             CALL SMLMA3(&
     &        A11, A12, A13,&
     &        A21, A22, A23,&
     &        A31, A32, A33,&
! =
     &       T1X(II), T2X(II), NX(II),&
     &       T1Y(II), T2Y(II), NY(II),&
     &       T1Z(II), T2Z(II), NZ(II),&
! *
     &  T1X(II)*DM1(NOD),   T1Y(II)*DM1(NOD),   T1Z(II)*DM1(NOD),&
     &  T2X(II)*DM2(NOD),   T2Y(II)*DM2(NOD),   T2Z(II)*DM2(NOD),&
     &   NX(II)*DM3(NOD),    NY(II)*DM3(NOD),    NZ(II)*DM3(NOD) )
! set values of (RMR^T +D)^{-1}
             RTDR(NOD)         =A11
             RTDR(NOD+NONODS  )=A21
             RTDR(NOD+2*NONODS)=A22
             RTDR(NOD+3*NONODS)=A31
             RTDR(NOD+4*NONODS)=A32
             RTDR(NOD+5*NONODS)=A33
      end do ! Was loop 30
!
! ELSE IF(D3)
           ENDIF
!
! ENOF IF(ROTAT)
           ENDIF

  end subroutine getrdr

end module getrdr_module
