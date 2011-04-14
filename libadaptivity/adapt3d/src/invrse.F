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
      SUBROUTINE INVRS2( A, DET, IFAIL )
C
      IMPLICIT NONE
C
      REAL A(4), DET
      INTEGER IFAIL
C
      REAL*8 DET8, V
C
      DET8  = 1.0D+00*A(1)*A(4) - 1.0D+00*A(2)*A(3)
      IFAIL = -1
      IF( ABS(DET8) .LT. 1E-30 ) RETURN
C
      DET8  = 1.0/DET8
      IFAIL = 0
      DET   = DET8
C
      V = A(1)
      A(1) = DET8*A(4)
      A(4) = DET8*V
C
      A(2) = -DET8*A(2)
      A(3) = -DET8*A(3)
C
      RETURN
      END
C
C
C
      SUBROUTINE INVRS3( A, DET, IFAIL )
C
      IMPLICIT NONE
C
      REAL A(9), DET
      INTEGER IFAIL
C
      REAL*8 V9, V1, V2, V3, V4, V5, V6, V7, V8, DET8
C
C - start with A(1) to hopefully get rest of A into cache...
C
      V9 = 1.0D+00*A(1)*A(5) - 1.0D+00*A(2)*A(4)
      V1 = 1.0D+00*A(5)*A(9) - 1.0D+00*A(6)*A(8)
      V2 = 1.0D+00*A(8)*A(3) - 1.0D+00*A(2)*A(9)
      V3 = 1.0D+00*A(2)*A(6) - 1.0D+00*A(3)*A(5)
      V4 = 1.0D+00*A(6)*A(7) - 1.0D+00*A(4)*A(9)
      V5 = 1.0D+00*A(1)*A(9) - 1.0D+00*A(3)*A(7)
      V6 = 1.0D+00*A(3)*A(4) - 1.0D+00*A(1)*A(6)
      V7 = 1.0D+00*A(4)*A(8) - 1.0D+00*A(5)*A(7)
      V8 = 1.0D+00*A(7)*A(2) - 1.0D+00*A(1)*A(8)
C
      DET8 = V1*A(1) + V4*A(2) + V7*A(3)
      IFAIL = -1
      IF( ABS(DET8) .LT. 1E-30 ) RETURN
C
      DET8  = 1.0/DET8
      IFAIL = 0
      DET   = DET8
C
      A(1) = V1*DET8
      A(2) = V2*DET8
      A(3) = V3*DET8
      A(4) = V4*DET8
      A(5) = V5*DET8
      A(6) = V6*DET8
      A(7) = V7*DET8
      A(8) = V8*DET8
      A(9) = V9*DET8
C
      RETURN
      END
