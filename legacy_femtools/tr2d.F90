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

module tr2d_module
   ! module full of random crap that we should get rid of
   use FLDebug
   use Quadrature
   use Shape_Functions
   use FETools
   use lagrot_module

   implicit none
   
   private
   
   public tr2d, tr3d, triqua, shatri, quadtetshapes
   
   contains
      
   SUBROUTINE TR2D(LOWQUA,NGI,NLOC,MLOC,&
     &     M,WEIGHT,N,NLX,NLY, &
     &     SNGI,SNLOC,SWEIGH,SN,SNLX )
!     This subroutine defines the shape functions M and N and their
!     derivatives at the Gauss points
!     For 3-D FLOW. 
      INTEGER NGI,NLOC,MLOC
      REAL M(MLOC,NGI),WEIGHT(NGI)
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI)
      INTEGER SNGI,SNLOC
      REAL SWEIGH(SNGI)
      REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI)
      REAL POSI,TLY
      INTEGER P,Q,CORN,GPOI,ILOC,JLOC,GI
      LOGICAL LOWQUA,GETNDP
      REAL WEIT(20),LX(20),LXP(2)
      INTEGER I
! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I
       
      ewrite(3,*) 'HERE 1 MLOC,NLOC,NGI=',MLOC,NLOC,NGI
      ewrite(3,*) 'HERE 2'

      IF((NLOC.NE.3).OR.(NGI.NE.3)) THEN
         FLAbort('PROBLEM IN TR2D')
      ENDIF

! This is for one point quadrature. 
      do  I=1,3
         NLX(I,1)=1.
         NLY(I,1)=1.
      end do

      WEIGHT(1)=1.

      N(1,1)=1./3.
      N(2,1)=1./3.
      N(3,1)=1./3.

! This is for 3 point quadrature. 
      do  GI=1,3
         NLX(1,GI)=1.
         NLY(1,GI)=0.

         NLX(2,GI)=0.
         NLY(2,GI)=1.

         NLX(3,GI)=-1.
         NLY(3,GI)=-1.
      end do


      N(1,1)=0.5
      N(2,1)=0.5
      N(3,1)=0.

      N(1,2)=0.
      N(2,2)=0.5
      N(3,2)=0.5

      N(1,3)=0.5
      N(2,3)=0.
      N(3,3)=0.5

      WEIGHT(1)=1./3.
      WEIGHT(2)=1./3.
      WEIGHT(3)=1./3.

      IF(MLOC.EQ.1) THEN
         M(1,1)=1.0
         M(1,2)=1.0
         M(1,3)=1.0
      ENDIF
      IF(MLOC.EQ.NLOC) M = N

      IF(SNGI.GT.0) THEN
         LXP(1)=-1
         LXP(2)= 1
         GETNDP=.FALSE.
         CALL LAGROT(WEIT,LX,SNGI,GETNDP)
         do  P=1,SNGI
            do  CORN=1,2! Was loop 327
               GPOI=P
               SN(CORN,GPOI)=0.5*(1.+LXP(CORN)*LX(P))
               SNLX(CORN,GPOI)=0.5*LXP(CORN)
               SWEIGH(GPOI)=WEIT(P)
            end do
         end do 
      ENDIF
      
   END subroutine tr2d

   SUBROUTINE TR3D(LOWQUA,NGI,NLOC,MLOC,&
     &     M,WEIGHT,N,NLX,NLY,NLZ,&
     &     SNGI,SNLOC,SWEIGH,SN,SNLX,SNLY)
!     This subroutine defines the shape functions M and N and their
!     derivatives at the Gauss points
!     For 3-D FLOW.
      INTEGER NGI,NLOC,MLOC
      REAL ALPHA,BETA
      PARAMETER(ALPHA=0.58541020,BETA=0.13819660)
      INTEGER SNGI,SNLOC
      REAL SWEIGH(SNGI)
      REAL SN(SNLOC,SNGI),SNLX(SNLOC,SNGI),SNLY(SNLOC,SNGI)
      REAL M(MLOC,NGI),WEIGHT(NGI)
      REAL N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
      REAL RUB(20)
      INTEGER P,Q,CORN,GPOI,ILOC,JLOC,GI
      LOGICAL LOWQUA
      INTEGER I
! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I

      IF((NLOC.NE.4).OR.(NGI.NE.4)) THEN
         ewrite(3,*) 'PROBLEM IN TR3D'
         FLAbort("STOP 201")
      ENDIF

! This is for one point. 

! This is for 4 point quadrature. 
      do  GI=1,NGI
         NLX(1,GI)=1.
         NLY(1,GI)=0.
         NLZ(1,GI)=0.

         NLX(2,GI)=0.
         NLY(2,GI)=1.
         NLZ(2,GI)=0.

         NLX(3,GI)=0.
         NLY(3,GI)=0.
         NLZ(3,GI)=1.

         NLX(4,GI)=-1.
         NLY(4,GI)=-1.
         NLZ(4,GI)=-1.
      end do 

      do I=1,4
         do GI=1,4
            N(I,GI)=BETA
         END DO
      END DO

      do I=1,4
         N(I,I)=ALPHA
         WEIGHT(I)=0.25
         IF(MLOC.EQ.1) M(1,I)=1.0
      END DO
      IF(MLOC.EQ.NLOC) M = N

      IF(SNGI.GT.0) THEN
         CALL TR2D(.FALSE.,SNGI,SNLOC,SNLOC,&
     &      RUB,SWEIGH,SN,SNLX,SNLY, &
     &      0,0,RUB,RUB,RUB )
      ENDIF 

      do I=1,NGI
         WEIGHT(I)=WEIGHT(I)/6.
      END DO
      do I=1,SNGI
         SWEIGH(I)=SWEIGH(I)/2.
      END DO
      
   end subroutine tr3d
      
   SUBROUTINE TRIQUA(L1, L2, L3, L4, WEIGHT, D3,NGI)
! This sub calculates the local corrds L1, L2, L3, L4 and 
! weights at the quadrature points. 
! If D3 it does this for 3Dtetrahedra elements else 
! triangular elements.
      INTEGER, intent(in) :: NGI
      LOGICAL, intent(in) :: D3
      REAL, intent(out) ::  L1(NGI), L2(NGI), L3(NGI), L4(NGI), WEIGHT(NGI)
! Local variables...
      type(quadrature_type) :: quad
      INTEGER I

      IF(D3) THEN
!     3D Tetrahedral elements.

         quad=make_quadrature(vertices=4, dim=3, ngi=ngi)
           
         weight=quad%weight
         L1=quad%l(:,1)
         L2=quad%l(:,2)
         L3=quad%l(:,3)
         L4=quad%l(:,4)

      else
!     2D Triangular elements.

         quad=make_quadrature(vertices=3, dim=2, ngi=ngi)

         weight=quad%weight
         L1=quad%l(:,1)
         L2=quad%l(:,2)
         L3=quad%l(:,3)

      end if
        
      call deallocate(quad)
      
   END subroutine triqua
      
   SUBROUTINE SHATRI(L1, L2, L3, L4, WEIGHT, D3,&
     &               NLOC,NGI,&
     &               N,NLX,NLY,NLZ) 
! Work out the shape functions and there derivatives...
      INTEGER, intent(in):: NLOC,NGI
      LOGICAL, intent(in):: D3
      REAL, intent(in):: L1(NGI), L2(NGI), L3(NGI), L4(NGI)
      REAL, intent(in):: WEIGHT(NGI)
      REAL, intent(out):: N(NLOC,NGI),NLX(NLOC,NGI),NLY(NLOC,NGI),NLZ(NLOC,NGI)
! Local variables...
      INTEGER GI
!
      IF(.NOT.D3) THEN
! Assume a triangle...
         IF((NLOC.EQ.6).OR.(NLOC.EQ.7)) THEN
            do  GI=1,NGI! Was loop 10
               N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
               N(2,GI)=(2.*L2(GI)-1.)*L2(GI)
               N(3,GI)=(2.*L3(GI)-1.)*L3(GI)

               N(4,GI)=4.*L1(GI)*L2(GI)
               N(5,GI)=4.*L2(GI)*L3(GI)
               N(6,GI)=4.*L1(GI)*L3(GI)
         
         ! nb L1+L2+L3+L4=1
         ! x-derivative...
               NLX(1,GI)=4.*L1(GI)-1.
               NLX(2,GI)=0.
               NLX(3,GI)=-4.*(1.-L2(GI))+4.*L1(GI) + 1.
               
               NLX(4,GI)=4.*L2(GI)
               NLX(5,GI)=-4.*L2(GI)
               NLX(6,GI)=4.*(1.-L2(GI))-8.*L1(GI)

         ! y-derivative...
               NLY(1,GI)=0.
               NLY(2,GI)=4.*L2(GI)-1.0
               NLY(3,GI)=-4.*(1.-L1(GI))+4.*L2(GI) + 1.

               NLY(4,GI)=4.*L1(GI)
               NLY(5,GI)=4.*(1.-L1(GI))-8.*L2(GI)
               NLY(6,GI)=-4.*L1(GI)
               IF(NLOC.EQ.7) THEN
         ! Bubble function...
                  N(7,GI)  =L1(GI)*L2(GI)*L3(GI)
                  NLX(7,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
                  NLY(7,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
               ENDIF
            end do
         ENDIF

         IF((NLOC.EQ.3).OR.(NLOC.EQ.4)) THEN
            do  GI=1,NGI! Was loop 20
               N(1,GI)=L1(GI)
               N(2,GI)=L2(GI)
               N(3,GI)=L3(GI)
               !
               NLX(1,GI)=1.0
               NLX(2,GI)=0.0
               NLX(3,GI)=-1.0
               !
               NLY(1,GI)=0.0
               NLY(2,GI)=1.0
               NLY(3,GI)=-1.0
               IF(NLOC.EQ.4) THEN
               ! Bubble function...
                  N(4,GI)  =L1(GI)*L2(GI)*L3(GI)
                  NLX(4,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
                  NLY(4,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
               ENDIF
            end do ! Was loop 20
         ENDIF

         IF(NLOC.EQ.1) THEN
            do  GI=1,NGI
                N(1,GI)=1.0
                NLX(1,GI)=0.0
                NLY(1,GI)=0.0
            end do
         ENDIF

      ENDIF

      IF(D3) THEN
      ! Assume a tet...
      ! This is for 5 point quadrature. 
         IF((NLOC.EQ.10).OR.(NLOC.EQ.11)) THEN
            call QuadTetShapes(n, nlx, nly, nlz, nloc, ngi)
            IF(NLOC.EQ.11) THEN
         ! Bubble function...
               N(11,:)  =L1*L2*L3*L4
               NLX(11,:)=L2*L3*(-L2-L3+1.)-2.*L1*L2*L3
               NLY(11,:)=L1*L3*(-L1-L3+1.)-2.*L1*L2*L3
               NLZ(11,:)=L1*L2*(-L1-L2+1.)-2.*L1*L2*L3
            ENDIF
         ENDIF 

         IF((NLOC.EQ.4).OR.(NLOC.EQ.5)) THEN
            do  GI=1,NGI
               N(1,GI)=L1(GI)
               N(2,GI)=L2(GI)
               N(3,GI)=L3(GI)
               N(4,GI)=L4(GI)

               NLX(1,GI)=1.0
               NLX(2,GI)=0
               NLX(3,GI)=0
               NLX(4,GI)=-1.0

               NLY(1,GI)=0.0
               NLY(2,GI)=1.0
               NLY(3,GI)=0.0
               NLY(4,GI)=-1.0

               NLZ(1,GI)=0.0
               NLZ(2,GI)=0.0
               NLZ(3,GI)=1.0
               NLZ(4,GI)=-1.0
               IF(NLOC.EQ.5) THEN
               ! Bubble function...
                  N(5,GI)  =L1(GI)*L2(GI)*L3(GI)*L4(GI)
                  NLX(5,GI)=L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                  NLY(5,GI)=L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
                  NLZ(5,GI)=L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
               ENDIF
            end do
         ENDIF

         IF(NLOC.EQ.1) THEN
            do  GI=1,NGI
              N(1,GI)=1.0
              NLX(1,GI)=0.0
              NLY(1,GI)=0.0
              NLZ(1,GI)=0.0
            end do
         ENDIF
      END if
   
   end subroutine shatri

   subroutine QuadTetShapes(n, nlx, nly, nlz, nloc, ngi)
   ! Work out the shape functions and their derivatives...
      integer, intent(in):: NLOC,NGI
      real n(nloc,ngi),nlx(nloc,ngi),nly(nloc,ngi),nlz(nloc,ngi)

      type(quadrature_type)::quad
      type(element_type)::shape

      ASSERT(nloc==10 .or. nloc==11)

      ! 4 vertices, 3 dimensions
      quad=make_quadrature(4, 3, ngi=ngi)

      ! 4 vertices, 3 dimensions, 2nd degree
      shape=make_element_shape(4, 3, 2, quad)

      n(1:10,:) = shape%n
      nlx(1:10,:)=shape%dn(:,:,X_)
      nly(1:10,:)=shape%dn(:,:,Y_)
      nlz(1:10,:)=shape%dn(:,:,Z_)

      call deallocate(shape)
      call deallocate(quad)

   end subroutine QuadTetShapes

end module tr2d_module
