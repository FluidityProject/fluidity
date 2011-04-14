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
      SUBROUTINE EDG4ND( BIGLST, NODLST, EDG1, EDG2, EDG3, EDG4,
     :                   N12, N34, N14, N23 )
C-----------------------------------------------------------------------
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
      INTEGER EDG1, EDG2, EDG3, EDG4, N12, N34, N14, N23
C
      INTEGER NOD1(2), NOD2(2), NOD3(2), IFLAG, NXT, LST, INEWED,nod4(2)
C
      CALL EDGINF( BIGLST, NODLST, EDG1, NOD1, INEWED, IFLAG, NXT, LST )
      CALL EDGINF( BIGLST, NODLST, EDG2, NOD2, INEWED, IFLAG, NXT, LST )
      CALL EDGINF( BIGLST, NODLST, EDG3, NOD3, INEWED, IFLAG, NXT, LST )
      CALL EDGINF( BIGLST, NODLST, EDG4, NOD4, INEWED, IFLAG, NXT, LST )
C
c      debug = .true.
      if( debug ) then
       print*,' '
       print*,'Into EDG4ND...',ierr
       print*,'   EDG1: ',EDG1,nod1
       print*,'   EDG2: ',EDG2,nod2
       print*,'   EDG3: ',EDG3,nod3
       print*,'   EDG4: ',EDG4,nod4
      end if
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( NOD1(1) .EQ. NOD2(1) ) THEN
         N12 = NOD1(1)
         N14 = NOD1(2)
         N23 = NOD2(2)
      ELSE IF( NOD1(2) .EQ. NOD2(1) ) THEN
         N12 = NOD1(2)
         N14 = NOD1(1)
         N23 = NOD2(2)
      ELSE IF( NOD1(1) .EQ. NOD2(2) ) THEN
         N12 = NOD1(1)
         N14 = NOD1(2)
         N23 = NOD2(1)
      ELSE IF( NOD1(2) .EQ. NOD2(2) ) THEN
         N12 = NOD1(2)
         N14 = NOD1(1)
         N23 = NOD2(1)
      ELSE
         IERR = -50
         PRINT*,'*** EDG4ND: EDGES 1 AND 2 DO NOT HAVE A COMMON NODE'
         RETURN
      END IF
C
      IF( NOD3(1) .EQ. N23 ) THEN
         N34 = NOD3(2)
      ELSE IF( NOD3(2) .EQ. N23 ) THEN
         N34 = NOD3(1)
      ELSE
         IERR = -50
         PRINT*,'*** EDG4ND: EDGES 2 AND 3 DO NOT HAVE A COMMON NODE'
         RETURN
      END IF
C
C - these are just consistency checks and can be removed...
C
      IF( NOD3(1) .NE. N23 .AND. NOD3(2) .NE. N23 ) THEN
         IERR = -51
         PRINT*,'*** EDG4ND: EDGE 3 DOES NOT CONTAIN N23'
         RETURN
      END IF
C
      CALL EDGINF( BIGLST, NODLST, EDG4, NOD1, INEWED, IFLAG, NXT, LST )
C
      IF( NOD1(1) .NE. N34 .AND. NOD1(2) .NE. N34 ) THEN
         IERR = -51
         PRINT*,'*** EDG4ND: EDGE 4 DOES NOT CONTAIN N34'
         RETURN
      END IF
C
      IF( NOD1(1) .NE. N14 .AND. NOD1(2) .NE. N14 ) THEN
         IERR = -51
         PRINT*,'*** EDG4ND: EDGE 4 DOES NOT CONTAIN N14'
         RETURN
      END IF
C
C - end of consistency checks
C
c      print*,'Leaving EDG4ND...'
c      print*,'   N12:',N12
c      print*,'   N34:',N34
c      print*,'   N14:',N14
c      print*,'   N23:',N23
c
      RETURN
      END
C
