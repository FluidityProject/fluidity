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
      SUBROUTINE NODEIG( BIGLST, NODLST, NOD1, V1, V2, V3, D1, D2, D3 )
C-----------------------------------------------------------------------
C
C - This subroutine takes the Hessian at the node supplied and
C - 'solves' it for the eigenvectors (V1,V2,V3) and the eigenvalues,
C - which are then turned into length-scales (D1,D2,D3).
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NOD1
C
      REAL V1(3), V2(3), V3(3), D1, D2, D3, a, b, c
C
      REAL XX, YY, ZZ
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      XX = NODLST( 1, NOD1 )
      YY = NODLST( 2, NOD1 )
      ZZ = NODLST( 3, NOD1 )
C
C - just use a set length scale over set directions for now
C
      INCLUDE 'fixeig.i'
C
      RETURN
      END
C
