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
      SUBROUTINE FIXEIG( XYZ, M, ANISOT )
C-----------------------------------------------------------------------
C
C - This subroutine sets up the metric with fixed eigvals and eigvects 
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER ANISOT
C
      REAL XYZ(3), M(3,3)
C
      REAL V1(3), V2(3), V3(3), D1, D2, D3, A, B, C, XX, YY, ZZ
C
      XX = XYZ(1)
      YY = XYZ(2)
      ZZ = XYZ(3)
C
      INCLUDE 'fixeig.i'
C
      V1(1) = V1(1)/D1
      V2(1) = V2(1)/D1
      V3(1) = V3(1)/D1
C
      V1(2) = V1(2)/D2
      V2(2) = V2(2)/D2
      V3(2) = V3(2)/D2
C
      V1(3) = V1(3)/D3
      V2(3) = V2(3)/D3
      V3(3) = V3(3)/D3
C
      M(1,1) = V1(1)*V1(1) + V1(2)*V1(2) + V1(3)*V1(3)
      M(2,1) = V2(1)*V1(1) + V2(2)*V1(2) + V2(3)*V1(3)
      M(3,1) = V3(1)*V1(1) + V3(2)*V1(2) + V3(3)*V1(3)
      M(1,2) = M(2,1)
      M(2,2) = V2(1)*V2(1) + V2(2)*V2(2) + V2(3)*V2(3)
      M(3,2) = V3(1)*V2(1) + V3(2)*V2(2) + V3(3)*V2(3)
      M(1,3) = M(3,1)
      M(2,3) = M(3,2)
      M(3,3) = V3(1)*V3(1) + V3(2)*V3(2) + V3(3)*V3(3)
C
      RETURN
      END
C
