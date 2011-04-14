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
      SUBROUTINE FACNRM( BIGLST, NODLST, N1, N2, N3, NX, NY, NZ )
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER N1, N2, N3
C
      REAL NX, NY, NZ
C
C -local
C
      integer i
c
      REAL*8 X12, Y12, Z12, X13, Y13, Z13, DOT, CROSS, NXD, NYD, NZD
c      real*8 x12err, y12err, z12err, x13err, y13err, z13err,
c     :       nxerr, nyerr, nzerr
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( N1. LE. 0 ) THEN
         IERR = -55
         PRINT*,'*** FACNRM: INVALID NODE POINTER 1'
         print*,n1,n2,n3
         RETURN
      END IF
C
      IF( N2. LE. 0 ) THEN
         IERR = -55
         PRINT*,'*** FACNRM: INVALID NODE POINTER 2'
         print*,n1,n2,n3
         RETURN
      END IF
C
      IF( N3. LE. 0 ) THEN
         IERR = -55
         PRINT*,'*** FACNRM: INVALID NODE POINTER 3'
         print*,n1,n2,n3
         RETURN
      END IF
C
      IF( N2 .EQ. N3 ) THEN
         IERR = -56
         PRINT*,'*** FACNRM: NODES 2 AND 3 ARE IDENTICAL'
         print*,n1,n2,n3
         RETURN
      END IF
C
      IF( N1 .EQ. N3 ) THEN
         IERR = -56
         PRINT*,'*** FACNRM: NODES 1 AND 3 ARE IDENTICAL'
         print*,n1,n2,n3
         RETURN
      END IF
C
      IF( N2 .EQ. N1 ) THEN
         IERR = -56
         PRINT*,'*** FACNRM: NODES 1 AND 2 ARE IDENTICAL'
         print*,n1,n2,n3
         RETURN
      END IF
C
      X12 = 1.0D+0*NODLST(1,N2) - 1.0D+0*NODLST(1,N1)
      Y12 = 1.0D+0*NODLST(2,N2) - 1.0D+0*NODLST(2,N1)
      Z12 = 1.0D+0*NODLST(3,N2) - 1.0D+0*NODLST(3,N1)
C
      X13 = 1.0D+0*NODLST(1,N3) - 1.0D+0*NODLST(1,N1)
      Y13 = 1.0D+0*NODLST(2,N3) - 1.0D+0*NODLST(2,N1)
      Z13 = 1.0D+0*NODLST(3,N3) - 1.0D+0*NODLST(3,N1)
c
c      x12err = max(abs(NODLST(1,N2)),abs(NODLST(1,N1))))*1e-14
c      y12err = max(abs(NODLST(2,N2)),abs(NODLST(2,N1))))*1e-14
c      z12err = max(abs(NODLST(3,N2)),abs(NODLST(3,N1))))*1e-14
c
c      x13err = max(abs(NODLST(1,N3)),abs(NODLST(1,N1))))*1e-14
c      y13err = max(abs(NODLST(2,N3)),abs(NODLST(2,N1))))*1e-14
c      z13err = max(abs(NODLST(3,N3)),abs(NODLST(3,N1))))*1e-14
C
      NXD = Y12*Z13 - Z12*Y13
c      nxerr = max(abs(y12*z13),abs(z12*y13))*1e-14
      NYD = Z12*X13 - X12*Z13
c      nyerr = max(abs(z12*x13),abs(x12*z13))*1e-14
      NZD = X12*Y13 - Y12*X13
c      nzerr = max(abs(x12*y13),abs(y12*x13))*1e-14
C
      DOT = X12*X12 + Y12*Y12 + Z12*Z12
      DOT = DOT*( X13*X13 + Y13*Y13 + Z13*Z13 )
C
      CROSS = NXD*NXD + NYD*NYD + NZD*NZD
      if( CROSS/DOT .le. 1e-12 ) then
         print*,'*** FACNRM: GOT COLINEAR POINTS IN FACE'
         IERR = -10000
         print*,'   nodes: ',n1,n2,n3,dot
         print*,'   node 1: ',(nodlst(i,n1),i=1,3)
         print*,'   1 to 2: ',(nodlst(i,n2)-nodlst(i,n1),i=1,3)
         print*,'   1 to 3: ',(nodlst(i,n3)-nodlst(i,n1),i=1,3)
         print*,'   normal: ',nx,ny,nz
         print*,'   1-2&3: ',sqrt(X12*X12 + Y12*Y12 + Z12*Z12),
     :                       sqrt(X13*X13 + Y13*Y13 + Z13*Z13)
         print*,'   dot: ',x12*x13+y12*y13+z12*z13
         NX = 0.0
         NY = 0.0
         NZ = 0.0
         RETURN
      else
         CROSS = 1.0/SQRT(CROSS)
      end if
C
      NX = NXD*CROSS
      NY = NYD*CROSS
      NZ = NZD*CROSS
C
      RETURN
      END
C
