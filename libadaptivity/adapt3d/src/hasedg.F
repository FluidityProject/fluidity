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
      INTEGER FUNCTION HASEDG( BIGLST, NODLST, IELM, IEDGE )
C-----------------------------------------------------------------------
C
C - This function checks to see if element IELM has edge IEDGE.
C - If so it returns .TRUE. otherwise it returns .FALSE.
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
      INTEGER IELM, IEDGE
C
      INTEGER NXT, LST
C
      HASEDG = 0
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** HASEDG: GOT ERROR ON ENTRY'
         RETURN
      END IF
C
      IF( IELM .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** HASEDG: INVALID ELEMENT POINTER GIVEN'
         print*,ielm,iedge
         RETURN
      END IF
C
      IF( IEDGE .LE. 0 ) THEN
         IERR = -2
         PRINT*,'*** HASEDG: INVALID EDGE POINTER GIVEN'
         print*,ielm,iedge
         RETURN
      END IF
C
      IF( BIGLST(4,IELM) .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** HASEDG: FIRST BLOCK IS NOT AN ELEMENT'
         print*,ielm,iedge,(biglst(nxt,ielm),nxt=1,4)
         RETURN
      END IF
C
      NXT = BIGLST(NXTBIG,IELM)
C
      IF( NXT .LE. 0 ) THEN
         IERR = -4
         PRINT*,'*** HASEDG: INVALID POINTER TO SECOND BLOCK'
         print*,ielm,iedge,nxt
         RETURN
      END IF
C
      IF( BIGLST(4,NXT) .LT. 0 ) THEN
         IERR = -5
         PRINT*,'*** HASEDG: SECOND BLOCK IS NOT AN ELEMENT'
         print*,ielm,iedge,nxt,biglst(4,nxt)
         RETURN
      END IF
C
      HASEDG = 1
      IF( BIGLST(1,NXT) .EQ. IEDGE ) GOTO 100
      HASEDG = 2
      IF( BIGLST(2,NXT) .EQ. IEDGE ) GOTO 100
      HASEDG = 3
      IF( BIGLST(3,NXT) .EQ. IEDGE ) GOTO 100
      HASEDG = 4
      IF( BIGLST(4,NXT) .EQ. IEDGE ) GOTO 100
C
      HASEDG = 0
C
      NXT = BIGLST(NXTBIG,NXT)
C
      IF( NXT .LE. 0 ) THEN
         IERR = -6
         PRINT*,'*** HASEDG: INVALID POINTER TO THIRD BLOCK'
         print*,ielm,iedge,nxt
         RETURN
      END IF
C
      IF( BIGLST(4,NXT) .LT. 0 ) THEN
         IERR = -7
         PRINT*,'*** HASEDG: THIRD BLOCK IS NOT AN ELEMENT'
         print*,ielm,iedge,nxt,biglst(4,nxt)
         RETURN
      END IF
C
      HASEDG = 5
      IF( BIGLST(1,NXT) .EQ. IEDGE ) GOTO 100
      HASEDG = 6
      IF( BIGLST(2,NXT) .EQ. IEDGE ) GOTO 100
C
      HASEDG = 0
C
      RETURN
C
 100  CONTINUE
C
      RETURN
      END
C
