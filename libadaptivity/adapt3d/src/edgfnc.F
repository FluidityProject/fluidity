C Copyright (C) 2006 Imperial College London and others.
C 
C Please see the AUTHORS file in the main source directory for a full list
C of copyright holders.
C 
C Adrian Umpleby
C Applied Modelling and Computation Group
C Department of Earth Science and Engineering
C Imperial College London
C 
C adrian@Imperial.ac.uk
C 
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public
C License as published by the Free Software Foundation; either
C version 2.1 of the License.
C 
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C 
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C USA
      REAL FUNCTION EDGFNC( BIGLST, NODLST, NOD1, NOD2, UNDLEN )
C-----------------------------------------------------------------------
C
C - Work out the functional for the edge joining nodes pointed to by NOD1, NOD2
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NOD1, NOD2
C
      REAL UNDLEN
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER J
C
      REAL X(2), Y(2), Z(2), XD, YD, ZD, M(9),
     :     V1(3), V2(3), V3(3), D1, D2, D3, D, drel
C
c      print*,'Into EDGFNC: ',nod1,nod2
      EDGFNC = 1E+30
      UNDLEN = 0.0
C
      IF( IERR .NE. 0 ) RETURN
C
cc      IF( NOD1 .LE. 0 ) THEN
cc         PRINT*,'*** EDGFNC: INVALID POINTER FOR NODE 1'
cc         IERR = -1
cc         RETURN
cc      END IF
C
cc      IF( NOD2 .LE. 0 ) THEN
cc         PRINT*,'*** EDGFNC: INVALID POINTER FOR NODE 2'
cc         IERR = -1
cc         RETURN
cc      END IF
C
cc      IF( NODLST(6,NOD1) .EQ. 0.0 ) THEN
cc         PRINT*,'*** EDGFNC: NODE 1 BLOCK IS EMPTY'
cc         IERR = -2
cc         RETURN
cc      END IF
C
cc      IF( NODLST(6,NOD2) .EQ. 0.0 ) THEN
cc         PRINT*,'*** EDGFNC: NODE 2 BLOCK IS EMPTY'
cc         IERR = -2
cc         RETURN
cc      END IF
C
C - find the vector joining the nodes...
C
      XD = NODLST(1,NOD1) - NODLST(1,NOD2)
      YD = NODLST(2,NOD1) - NODLST(2,NOD2)
      ZD = NODLST(3,NOD1) - NODLST(3,NOD2)
C
C - and the metric for the edge (the average for the two node metrics)
C
      DO J = 1, 9
         M(J) = ( NODLST(J+6,NOD1) + NODLST(J+6,NOD2) ) / 2
      END DO
C
C - work out the edge vector in the metric space...
C
      D1 = M(1)*XD + M(4)*YD + M(7)*ZD
      D2 = M(2)*XD + M(5)*YD + M(8)*ZD
      D3 = M(3)*XD + M(6)*YD + M(9)*ZD
C
C - and the length (squared) of the edge in metric space (D) and in non-metric space (DREL)
C
      D = D1*XD + D2*YD + D3*ZD
      DREL = XD*XD + YD*YD + ZD*ZD
C
C      X(1) = NODLST(1,NOD1)
C      Y(1) = NODLST(2,NOD1)
C      Z(1) = NODLST(3,NOD1)
C
C      X(2) = NODLST(1,NOD2)
C      Y(2) = NODLST(2,NOD2)
C      Z(2) = NODLST(3,NOD2)
C
c      CALL EDGEIG( BIGLST, NODLST, NOD1, NOD2, V1, V2, V3, D1, D2, D3 )
c      CALL UNDSTT( XD, YD, ZD, 1, V1, V2, V3, D1, D2, D3 )
C
c      IF( IERR .NE. 0 ) RETURN
C
C      XD = X(2) - X(1)
C      YD = Y(2) - Y(1)
C      ZD = Z(2) - Z(1)
C
c      D = XD*XD + YD*YD + ZD*ZD
C
C - Now find the length of the edge in metric space (UNDLEN) - provided
C - the metric is not too distorted!
C
      IF( D/DREL .GT. 1E-14 .or. d .gt. 1E-4 ) THEN
         UNDLEN = SQRT(D)
C
      else if( d .lt. 0.0 ) then
         print*,'+++ EDGFNC: Warning, got negative edge: ',d,drel
         print*,m(1),m(4),m(7)
         print*,m(2),m(5),m(8)
         print*,m(3),m(6),m(9)
         print*,'node 1: ',(nodlst(j,nod1),j=1,3)
         print*,nodlst(7,nod1),nodlst(10,nod1),nodlst(13,nod1)
         print*,nodlst(8,nod1),nodlst(11,nod1),nodlst(14,nod1)
         print*,nodlst(9,nod1),nodlst(12,nod1),nodlst(15,nod1)
         print*,'node 2: ',(nodlst(j,nod2),j=1,3)
         print*,nodlst(7,nod2),nodlst(10,nod2),nodlst(13,nod2)
         print*,nodlst(8,nod2),nodlst(11,nod2),nodlst(14,nod2)
         print*,nodlst(9,nod2),nodlst(12,nod2),nodlst(15,nod2)
         undlen = 0.0
      else
         print*,'+++ EDGFNC: Warning, got very small edge: ',d,drel
         print*,m(1),m(4),m(7)
         print*,m(2),m(5),m(8)
         print*,m(3),m(6),m(9)
         print*,'node 1: ',(nodlst(j,nod1),j=1,3)
         print*,nodlst(7,nod1),nodlst(10,nod1),nodlst(13,nod1)
         print*,nodlst(8,nod1),nodlst(11,nod1),nodlst(14,nod1)
         print*,nodlst(9,nod1),nodlst(12,nod1),nodlst(15,nod1)
         print*,'node 2: ',(nodlst(j,nod2),j=1,3)
         print*,nodlst(7,nod2),nodlst(10,nod2),nodlst(13,nod2)
         print*,nodlst(8,nod2),nodlst(11,nod2),nodlst(14,nod2)
         print*,nodlst(9,nod2),nodlst(12,nod2),nodlst(15,nod2)
         undlen = sqrt(d)
      end if
C
C - Work out the edge functional
C
      D = 1.0 - UNDLEN
C
      IF( EDGON ) THEN
         EDGFNC = D*D
      ELSE
         EDGFNC = 0.0
      END IF
C
c      print*,'Leaving EDGFNC: ',d*d,undlen
      RETURN
      END
C
