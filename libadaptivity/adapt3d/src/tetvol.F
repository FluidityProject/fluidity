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
      REAL FUNCTION TETVOL( X, Y, Z )
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      REAL X(4), Y(4), Z(4)
C
      REAL VOL, X12, X13, X14, Y12, Y13, Y14, Z12, Z13, Z14
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
c      X12 = X(2)
c      X13 = X(3)
c      X14 = X(4)
c      Y12 = Y(2)
c      Y13 = Y(3)
c      Y14 = Y(4)
c      Z12 = Z(2)
c      Z13 = Z(3)
c      Z14 = Z(4)
C
c      X23 = X(3) - X(2)
c      X24 = X(4) - X(2)
c      X34 = X(4) - X(3)
c      Y23 = Y(3) - Y(2)
c      Y24 = Y(4) - Y(2)
c      Y34 = Y(4) - Y(3)
c      Z23 = Z(3) - Z(2)
c      Z24 = Z(4) - Z(2)
c      Z34 = Z(4) - Z(3)
C
C      VOL = X12*Y13*Z14 + X13*Y14*Z12 + X14*Y12*Z13
C     :    - X14*Y13*Z12 - X13*Y12*Z14 - X12*Y14*Z13
C
      VOL = X12*( Y13*Z14 - Y14*Z13 )
     :    + X13*( Y14*Z12 - Y12*Z14 )
     :    + X14*( Y12*Z13 - Y13*Z12 )
C
      TETVOL = VOL/6
C
      RETURN
      END
C
