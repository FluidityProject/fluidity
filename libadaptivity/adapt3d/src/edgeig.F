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
      SUBROUTINE EDGEIG( BIGLST, NODLST, NOD1, NOD2, V1, V2, V3,
     :                   D1, D2, D3 )
C-----------------------------------------------------------------------
C
C - This subroutine takes the hessian(s) at the two nodes supplied and
C - forms an 'average' hessian for the edge, which is then made +ve def.
C - and 'solved' for the eigenvectors (V1,V2,V3) and the eigenvalues,
C - which are then turned into length-scales (D1,D2,D3).
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NOD1, NOD2
C
      REAL V1(3), V2(3), V3(3), D1, D2, D3, a, b, c
C
      REAL X(2), Y(2), Z(2), XX, YY, ZZ, D1N1, D2N1, D3N1, D1N2, D2N2,
     :     D3N2, V1N1(3), V2N1(3), V3N1(3), V1N2(3), V2N2(3), V3N2(3)
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      X(1) = NODLST( 1, NOD1 )
      Y(1) = NODLST( 2, NOD1 )
      Z(1) = NODLST( 3, NOD1 )
C
      X(2) = NODLST( 1, NOD2 )
      Y(2) = NODLST( 2, NOD2 )
      Z(2) = NODLST( 3, NOD2 )
C
      XX = ( X(1) + X(2) ) / 2
      YY = ( Y(1) + Y(2) ) / 2
      ZZ = ( Z(1) + Z(2) ) / 2
C
c      CALL NODEIG( NOD1, V1N1, V2N1, V3N1, D1N1, D2N1, D3N1 )
c      CALL NODEIG( NOD2, V1N2, V2N2, V3N2, D1N2, D2N2, D3N2 )
C
C - just use a set length scale over set directions for now
C
      INCLUDE 'fixeig.i'
C
      return
c
      IF( ABS(D1N1/D1) .LT. 0.4 .OR.
     :    ABS(D2N1/D2) .LT. 0.4 .OR.
     :    ABS(D3N1/D3) .LT. 0.4 ) THEN
         D1 = D1N1
         D2 = D2N1
         D3 = D3N1
         V1(1) = V1N1(1)
         V1(2) = V1N1(2)
         V1(3) = V1N1(3)
         V2(1) = V2N1(1)
         V2(2) = V2N1(2)
         V2(3) = V2N1(3)
         V3(1) = V3N1(1)
         V3(2) = V3N1(2)
         V3(3) = V3N1(3)
      ELSE IF( ABS(D1N2/D1) .LT. 0.4 .OR.
     :         ABS(D2N2/D2) .LT. 0.4 .OR.
     :         ABS(D3N2/D3) .LT. 0.4 ) THEN
         D1 = D1N2
         D2 = D2N2
         D3 = D3N2
         V1(1) = V1N2(1)
         V1(2) = V1N2(2)
         V1(3) = V1N2(3)
         V2(1) = V2N2(1)
         V2(2) = V2N2(2)
         V2(3) = V2N2(3)
         V3(1) = V3N2(1)
         V3(2) = V3N2(2)
         V3(3) = V3N2(3)
      END IF
C
      RETURN
      END
C
