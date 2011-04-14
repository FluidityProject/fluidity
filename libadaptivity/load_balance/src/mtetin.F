C Copyright (C) 2006 Imperial College London and others.
C
C Please see the AUTHORS file in the main source directory for a full list
C of copyright holders.
C
C Mr Adrian Umpleby
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

      SUBROUTINE MTETIN( X, Y, Z, M, VOL, AREAS, L, RADIUS, QUALTY )
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      REAL X(4), Y(4), Z(4), VOL, AREAS(4), RADIUS, M(3,3), QUALTY,
     :     L(6)
C
      REAL A, S, X12, X13, X14, X23, X24, X34, DET,
     :     Y12, Y13, Y14, Y23, Y24, Y34, Z12, Z13, Z14, Z23, Z24, Z34,
     :     V1, V2, V3, L12, L13, L14, L23, L24, L34, ALPHA, ALPH2
C
c      PARAMETER( ALPHA = SQRT(6.0)/12, ALPH2 = 1.0/12/SQRT(6.0) )
c
      PARAMETER( ALPHA = 0.204124145, ALPH2 = 0.03402069 )
C
      X12 = X(2) - X(1)
      X13 = X(3) - X(1)
      X14 = X(4) - X(1)
      Y12 = Y(2) - Y(1)
      Y13 = Y(3) - Y(1)
      Y14 = Y(4) - Y(1)
      Z12 = Z(2) - Z(1)
      Z13 = Z(3) - Z(1)
      Z14 = Z(4) - Z(1)
C
      X23 = X(3) - X(2)
      X24 = X(4) - X(2)
      X34 = X(4) - X(3)
      Y23 = Y(3) - Y(2)
      Y24 = Y(4) - Y(2)
      Y34 = Y(4) - Y(3)
      Z23 = Z(3) - Z(2)
      Z24 = Z(4) - Z(2)
      Z34 = Z(4) - Z(3)
C
      VOL = X12*( Y13*Z14 - Y14*Z13 )
     :    + X13*( Y14*Z12 - Y12*Z14 )
     :    + X14*( Y12*Z13 - Y13*Z12 )
C
      DET = M(1,1)*( M(2,2)*M(3,3) - M(2,3)*M(3,2) )
     :    + M(2,1)*( M(3,2)*M(1,3) - M(1,2)*M(3,3) )
     :    + M(3,1)*( M(1,2)*M(2,3) - M(2,2)*M(1,3) )
C
      if( det .gt. 1e-20 ) then
         VOL = SQRT(DET)*VOL
      else if( det .lt. 0.0 ) then
         write(0, *) "WARNING: negative determinant = ", det
         write(0, *) m(1,1), m(1,2), m(1,3)
         write(0, *) m(2,1), m(2,2), m(2,3)
         write(0, *) m(3,1), m(3,2), m(3,3)
         vol = 0.0
         det = 0.0
      else
         vol = sqrt(det)*vol
         write(0, *) "WARNING: got small det", det
         write(0, *) m(1,1), m(1,2), m(1,3)
         write(0, *) m(2,1), m(2,2), m(2,3)
         write(0, *) m(3,1), m(3,2), m(3,3)
      end if
C
      S = 0.0
C
      V1 = M(1,1)*X12 + M(1,2)*Y12 + M(1,3)*Z12
      V2 = M(2,1)*X12 + M(2,2)*Y12 + M(2,3)*Z12
      V3 = M(3,1)*X12 + M(3,2)*Y12 + M(3,3)*Z12
C
      L12 = MAX( X12*V1 + Y12*V2 + Z12*V3, 0.0 )
      L(1) = SQRT(L12)
C
      V1 = M(1,1)*X13 + M(1,2)*Y13 + M(1,3)*Z13
      V2 = M(2,1)*X13 + M(2,2)*Y13 + M(2,3)*Z13
      V3 = M(3,1)*X13 + M(3,2)*Y13 + M(3,3)*Z13
C
      L13 = MAX( X13*V1 + Y13*V2 + Z13*V3, 0.0 )
      L(2) = SQRT(L13)
c      QUALTY = MAX( L12, L13 )
C
      V1 = M(1,1)*X14 + M(1,2)*Y14 + M(1,3)*Z14
      V2 = M(2,1)*X14 + M(2,2)*Y14 + M(2,3)*Z14
      V3 = M(3,1)*X14 + M(3,2)*Y14 + M(3,3)*Z14
C
      L14 = MAX( X14*V1 + Y14*V2 + Z14*V3, 0.0 )
      L(3) = SQRT(L14)
c      QUALTY = MAX( QUALTY, L14 )
C
      V1 = M(1,1)*X23 + M(1,2)*Y23 + M(1,3)*Z23
      V2 = M(2,1)*X23 + M(2,2)*Y23 + M(2,3)*Z23
      V3 = M(3,1)*X23 + M(3,2)*Y23 + M(3,3)*Z23
C
      L23 = MAX( X23*V1 + Y23*V2 + Z23*V3, 0.0 )
      L(4) = SQRT(L23)
c      QUALTY = MAX( QUALTY, L23 )
C
      V1 = M(1,1)*X24 + M(1,2)*Y24 + M(1,3)*Z24
      V2 = M(2,1)*X24 + M(2,2)*Y24 + M(2,3)*Z24
      V3 = M(3,1)*X24 + M(3,2)*Y24 + M(3,3)*Z24
C
      L24 = MAX( X24*V1 + Y24*V2 + Z24*V3, 0.0 )
      L(5) = SQRT(L24)
c      QUALTY = MAX( QUALTY, L24 )
C
      V1 = M(1,1)*X34 + M(1,2)*Y34 + M(1,3)*Z34
      V2 = M(2,1)*X34 + M(2,2)*Y34 + M(2,3)*Z34
      V3 = M(3,1)*X34 + M(3,2)*Y34 + M(3,3)*Z34
C
      L34 = MAX( X34*V1 + Y34*V2 + Z34*V3, 0.0 )
      L(6) = SQRT(L34)
c      QUALTY = MAX( QUALTY, L34 )
      qualty = l(1) + l(2) + l(3) + l(4) + l(5) + l(6)
C
      A = L12 + L23 - L13
      A = L12*L23*4 - A*A
      if( a .gt. 0.0 ) then
         A = sqrt(A)
      else
         a = 0.0
      end if
      AREAS(1) = A/4
      S = S + A
C
      A = L12 + L24 - L14
      A = L12*L24*4 - A*A
      if( a .gt. 0.0 ) then
         A = sqrt(A)
      else
         a = 0.0
      end if
      AREAS(2) = A/4
      S = S + A
C
      A = L23 + L34 - L24
      A = L23*L34*4 - A*A
      if( a .gt. 0.0 ) then
         A = sqrt(A)
      else
         a = 0.0
      end if
      AREAS(4) = A/4
      S = S + A
C
      A = L13 + L14 - L34
      A = L13*L14*4 - A*A
      if( a .gt. 0.0 ) then
         A = sqrt(A)
      else
         a = 0.0
      end if
      AREAS(3) = A/4
      S = S + A
C
      RADIUS = VOL*2/S
C
      VOL = VOL/6
C
      IF( ABS(RADIUS) .LT. 1E-20 ) THEN
         IF( RADIUS .LT. 0.0 ) THEN
            QUALTY = -1E+20
         ELSE
            QUALTY = 1E+20
         END IF
      ELSE
c         QUALTY = ALPHA*SQRT(QUALTY)/ABS(RADIUS)
         QUALTY = ALPH2*QUALTY/RADIUS
      END IF
C
      RETURN
      END
