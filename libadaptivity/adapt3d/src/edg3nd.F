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
      SUBROUTINE EDG3ND( BIGLST, NODLST, EDG1, EDG2, EDG3,
     :                   N123, N1, N2, N3 )
C-----------------------------------------------------------------------
C
C - This subroutine takes three edges with a common node and returns
C - the common node (N123) and the other three hanging nodes for
C - each edge (N1 for EDG1, N2 for EDG2 & N3 for EDG3).
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
      INTEGER EDG1, EDG2, EDG3, N123, N1, N2, N3
C
      INTEGER NOD1(2), NOD2(2), NOD3(2), IFLAG, NXT, LST, INEWED
C
c      print*,' '
c      print*,'Into EDG3ND...',ierr
c      print*,'   EDG1: ',EDG1
c      print*,'   EDG2: ',EDG2
c      print*,'   EDG3: ',EDG3
C
      IF( IERR .NE. 0 ) RETURN
C
      CALL EDGINF( BIGLST, NODLST, EDG1, NOD1, INEWED, IFLAG, NXT, LST )
      CALL EDGINF( BIGLST, NODLST, EDG2, NOD2, INEWED, IFLAG, NXT, LST )
      CALL EDGINF( BIGLST, NODLST, EDG3, NOD3, INEWED, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( NOD1(1) .EQ. NOD2(1) ) THEN
         N123 = NOD1(1)
         N1   = NOD1(2)
         N2   = NOD2(2)
      ELSE IF( NOD1(1) .EQ. NOD2(2) ) THEN
         N123 = NOD1(1)
         N1   = NOD1(2)
         N2   = NOD2(1)
      ELSE IF( NOD1(2) .EQ. NOD2(1) ) THEN
         N123 = NOD1(2)
         N1   = NOD1(1)
         N2   = NOD2(2)
      ELSE IF( NOD1(2) .EQ. NOD2(2) ) THEN
         N123 = NOD1(2)
         N1   = NOD1(1)
         N2   = NOD2(1)
      ELSE
         IERR = -50
         PRINT*,'*** EDG3ND: EDGES 1 AND 2 DO NOT HAVE A COMMON NODE'
         RETURN
      END IF
C
      IF( NOD3(1) .EQ. N123 ) THEN
         N3 = NOD3(2)
      ELSE IF( NOD3(2) .EQ. N123 ) THEN
         N3 = NOD3(1)
      ELSE
         IERR = -51
         PRINT*,'*** EDG3ND: EDGE 3 DOES NOT HAVE THE COMMON NODE'
         RETURN
      END IF
C
c      print*,'Leaving EDG3ND...'
c      print*,'     N1: ',N1
c      print*,'     N2: ',N2
c      print*,'     N3: ',N3
c      print*,'   N123: ',N123
C
      RETURN
      END
C
