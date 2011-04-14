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
      SUBROUTINE UNDSTT( X, Y, Z, NPTS, V1, V2, V3, D1, D2, D3 )
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NPTS
C
      REAL X(NPTS), Y(NPTS), Z(NPTS), V1(3), V2(3), V3(3), D1, D2, D3
C
      INTEGER I
C
      REAL XX, YY, ZZ, A, DOT, A1, A2, A3
C
      A1 = V1(1)*V1(1) + V1(2)*V1(2) + V1(3)*V1(3)
      A2 = V2(1)*V2(1) + V2(2)*V2(2) + V2(3)*V2(3)
      A3 = V3(1)*V3(1) + V3(2)*V3(2) + V3(3)*V3(3)
C
      A1 = (1.0/D1 - 1.0)/SQRT(A1)
      A2 = (1.0/D2 - 1.0)/SQRT(A2)
      A3 = (1.0/D3 - 1.0)/SQRT(A3)
C
      DO I = 1, NPTS
C
         DOT = X(I)*V1(1) + Y(I)*V1(2) + Z(I)*V1(3)
         A = DOT*A1
C
         X(I) = X(I) + A*V1(1)
         Y(I) = Y(I) + A*V1(2)
         Z(I) = Z(I) + A*V1(3)
C
         DOT = X(I)*V2(1) + Y(I)*V2(2) + Z(I)*V2(3)
         A = DOT*A2
C
         X(I) = X(I) + A*V2(1)
         Y(I) = Y(I) + A*V2(2)
         Z(I) = Z(I) + A*V2(3)
C
         DOT = X(I)*V3(1) + Y(I)*V3(2) + Z(I)*V3(3)
         A = DOT*A3
C
         X(I) = X(I) + A*V3(1)
         Y(I) = Y(I) + A*V3(2)
         Z(I) = Z(I) + A*V3(3)
C
      END DO
C
      RETURN
      END
C
