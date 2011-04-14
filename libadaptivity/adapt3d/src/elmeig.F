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
      SUBROUTINE ELMEIG( BIGLST, NODLST, NOD1, NOD2, NOD3, NOD4,
     :                   V1, V2, V3, D1, D2, D3 )
C-----------------------------------------------------------------------
C
C - This subroutine takes the Hessian(s) at the four nodes supplied and
C - forms an 'average' Hessian for the element, which is then made +ve def.
C - and 'solved' for the eigenvectors (V1,V2,V3) and the eigenvalues,
C - which are then turned into length-scales (D1,D2,D3).
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NOD1, NOD2, NOD3, NOD4
C
      REAL V1(3), V2(3), V3(3), D1, D2, D3, a, b, c
C
      REAL X(4), Y(4), Z(4), XX, YY, ZZ, D1N1, D2N1, D3N1, D1N2, D2N2,
     :     D3N2, D1N3, D2N3, D3N3, D1N4, D2N4, D3N4,
     :     V1N1(3), V2N1(3), V3N1(3), V1N2(3), V2N2(3), V3N2(3),
     :     V1N3(3), V2N3(3), V3N3(3), V1N4(3), V2N4(3), V3N4(3),
     :     a1234, a2134, a3214, a4231, a1324, a1432, a1243
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
      X(3) = NODLST( 1, NOD3 )
      Y(3) = NODLST( 2, NOD3 )
      Z(3) = NODLST( 3, NOD3 )
C
      X(4) = NODLST( 1, NOD4 )
      Y(4) = NODLST( 2, NOD4 )
      Z(4) = NODLST( 3, NOD4 )
C
      XX = ( X(1) + X(2) + X(3) + X(4) ) / 4
      YY = ( Y(1) + Y(2) + Y(3) + Y(4) ) / 4
      ZZ = ( Z(1) + Z(2) + Z(3) + Z(4) ) / 4
c
c      if( abs(xx) .gt. 1e-5 ) then
c      a1234 = ( ( x(1) + x(2) ) + x(3) ) + x(4)
c      a2134 = ( ( x(2) + x(1) ) + x(3) ) + x(4)
c      a3214 = ( ( x(3) + x(2) ) + x(1) ) + x(4)
c      a4231 = ( ( x(4) + x(2) ) + x(3) ) + x(1)
c      a1324 = ( ( x(1) + x(3) ) + x(2) ) + x(4)
c      a1432 = ( ( x(1) + x(4) ) + x(3) ) + x(2)
c      a1243 = ( ( x(1) + x(2) ) + x(4) ) + x(3)
c      if( abs(a1234/a2134)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 2 on X: ',a1234,a2134
c      if( abs(a1234/a3214)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 3 on X: ',a1234,a3214
c      if( abs(a1234/a4231)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 4 on X: ',a1234,a4231
c      if( abs(a1234/a1324)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 2 & 3 on X: ',a1234,a1324
c      if( abs(a1234/a1432)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 2 & 4 on X: ',a1234,a1432
c      if( abs(a1234/a1243)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 3 & 4 on X: ',a1234,a1243
c      end if
cc
c      if( abs(yy) .gt. 1e-5 ) then
c      a1234 = ( ( y(1) + y(2) ) + y(3) ) + y(4)
c      a2134 = ( ( y(2) + y(1) ) + y(3) ) + y(4)
c      a3214 = ( ( y(3) + y(2) ) + y(1) ) + y(4)
c      a4231 = ( ( y(4) + y(2) ) + y(3) ) + y(1)
c      a1324 = ( ( y(1) + y(3) ) + y(2) ) + y(4)
c      a1432 = ( ( y(1) + y(4) ) + y(3) ) + y(2)
c      a1243 = ( ( y(1) + y(2) ) + y(4) ) + y(3)
c      if( abs(a1234/a2134)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 2 on Y: ',a1234,a2134
c      if( abs(a1234/a3214)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 3 on Y: ',a1234,a3214
c      if( abs(a1234/a4231)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 4 on Y: ',a1234,a4231
c      if( abs(a1234/a1324)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 2 & 3 on Y: ',a1234,a1324
c      if( abs(a1234/a1432)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 2 & 4 on Y: ',a1234,a1432
c      if( abs(a1234/a1243)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 3 & 4 on Y: ',a1234,a1243
c      end if
cc
c      if( abs(zz) .gt. 1e-5 ) then
c      a1234 = ( ( z(1) + z(2) ) + z(3) ) + z(4)
c      a2134 = ( ( z(2) + z(1) ) + z(3) ) + z(4)
c      a3214 = ( ( z(3) + z(2) ) + z(1) ) + z(4)
c      a4231 = ( ( z(4) + z(2) ) + z(3) ) + z(1)
c      a1324 = ( ( z(1) + z(3) ) + z(2) ) + z(4)
c      a1432 = ( ( z(1) + z(4) ) + z(3) ) + z(2)
c      a1243 = ( ( z(1) + z(2) ) + z(4) ) + z(3)
c      if( abs(a1234/a2134)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 2 on Z: ',a1234,a2134
c      if( abs(a1234/a3214)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 3 on Z: ',a1234,a3214
c      if( abs(a1234/a4231)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 1 & 4 on Z: ',a1234,a4231
c      if( abs(a1234/a1324)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 2 & 3 on Z: ',a1234,a1324
c      if( abs(a1234/a1432)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 2 & 4 on Z: ',a1234,a1432
c      if( abs(a1234/a1243)-1 .gt. 1e-4 )
c     :         print*,'ELMEIG: SWAPPED 3 & 4 on Z: ',a1234,a1243
c      end if
C
c      CALL NODEIG( NOD1, V1N1, V2N1, V3N1, D1N1, D2N1, D3N1 )
c      CALL NODEIG( NOD2, V1N2, V2N2, V3N2, D1N2, D2N2, D3N2 )
c      CALL NODEIG( NOD3, V1N3, V2N3, V3N3, D1N3, D2N3, D3N3 )
c      CALL NODEIG( NOD4, V1N4, V2N4, V3N4, D1N4, D2N4, D3N4 )
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
      ELSE IF( ABS(D1N3/D1) .LT. 0.4 .OR.
     :         ABS(D2N3/D2) .LT. 0.4 .OR.
     :         ABS(D3N3/D3) .LT. 0.4 ) THEN
         D1 = D1N3
         D2 = D2N3
         D3 = D3N3
         V1(1) = V1N3(1)
         V1(2) = V1N3(2)
         V1(3) = V1N3(3)
         V2(1) = V2N3(1)
         V2(2) = V2N3(2)
         V2(3) = V2N3(3)
         V3(1) = V3N3(1)
         V3(2) = V3N3(2)
         V3(3) = V3N3(3)
      ELSE IF( ABS(D1N4/D1) .LT. 0.4 .OR.
     :         ABS(D2N4/D2) .LT. 0.4 .OR.
     :         ABS(D3N4/D3) .LT. 0.4 ) THEN
         D1 = D1N4
         D2 = D2N4
         D3 = D3N4
         V1(1) = V1N4(1)
         V1(2) = V1N4(2)
         V1(3) = V1N4(3)
         V2(1) = V2N4(1)
         V2(2) = V2N4(2)
         V2(3) = V2N4(3)
         V3(1) = V3N4(1)
         V3(2) = V3N4(2)
         V3(3) = V3N4(3)
      END IF
C
      RETURN
      END
C
