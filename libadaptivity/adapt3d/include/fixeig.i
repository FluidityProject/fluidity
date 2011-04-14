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
C
C - this code is for fixing the eigenvalues and eigenvectors for the
C - subroutines ELMEIG, EDGEIG and NODEIG.
C
c      a = 1.0/sqrt(3.0)
c
c      V1(1) = a
c      V1(2) = -a
c      V1(3) = -a
C
c      V2(1) = a
c      V2(2) = a
c      V2(3) = -a
C
c      V3(1) = a
c      V3(2) = a
c      V3(3) = a
c
      D1 = 0.5
      D2 = 0.5
      D3 = 0.5
c
      a = 1.0/sqrt(2.0)
c
      v1(1) = a
      v1(2) = a
      v1(3) = 0.0
c
      v2(1) = -a
      v2(2) = a
      v2(3) = 0.0
C
      V3(1) = 0.0
      V3(2) = 0.0
      V3(3) = 1.0
c
c      b = -abs(1.0-xx-yy)
c
c      D1 = 0.13 - 0.124*exp(b*3)
c      D2 = 0.13 - 0.124*exp(a*3)
c      D1 = 0.1 - 0.094*exp(b*3)
c      D2 = 0.1 - 0.094*exp(a*3)
c      D3 = 0.35 - 0.34*exp(c*3)
c      D3 = 0.35 - zz/3
c      D1 = 0.2 - 0.19*exp(b)
c      D2 = 0.2 - 0.19*exp(a)
      if( anisot .lt. 0 ) then
c
         d1 = float(-anisot)/100
         d2 = d1
         d3 = d1
c
      else if( anisot .eq. 2 ) then
         a = -abs(xx-yy)
         b = -abs(1.0-xx-yy)
         c = -abs(10.5-zz)
c         d2 = 0.2
c         D1 = 0.2 - 0.19*exp(a)
c         D3 = 0.3
         d1 = 0.1  - 0.092*exp(a*3)
         D2 = 0.1  - 0.092*exp(b*3)
         d3 = 0.25 - 0.237*exp(c*2)
      else if( anisot .eq. 1 ) then
         a = 1.0/sqrt(3.0)
c
c         V1(1) = a
c         V1(2) = -a
c         V1(3) = -a
C
c         V2(1) = a
c         V2(2) = a
c         V2(3) = -a
C
c         V3(1) = a
c         V3(2) = a
c         V3(3) = a
c
c         d1 = 0.03
c         d2 = 0.15
c         D3 = 0.4
         d1 = 0.05
         d2 = 0.05
         d3 = 0.05
      else if( anisot .eq. 3 ) then
         a = -abs(0.0707107+xx-yy)
c         b = -abs(1.0-xx-yy)
         b = -abs(1.0707107-xx-yy)
         c = -abs(0.55-zz)
c         D2 = 0.2 - 0.19*exp(b)
c         d1 = 0.2
c         D3 = 0.3
         d1 = 0.1  - 0.092*exp(a*3)
         D2 = 0.1  - 0.092*exp(b*3)
         d3 = 0.25 - 0.237*exp(c*2)
      end if
c      D2 = 0.2
c      D1 = 0.2
c      D2 = 0.15
c      D1 = 0.05
c      D2 = 0.1
c      D1 = 0.1
c      D1 = 0.05
c      D2 = 0.04 - 0.039*exp(a*10)
c      D3 = 0.05
c      d2 = 0.05
c
c      V1(1) = 1.0
c      V1(2) = 0.0
c      V1(3) = 0.0
C
c      V2(1) = 0.0
c      V2(2) = 1.0
c      V2(3) = 0.0
C
c      V3(1) = 0.0
c      V3(2) = 0.0
c      V3(3) = 1.0
C
c      D1 = 0.125
c      D2 = 0.25
c      D3 = 0.5
C
