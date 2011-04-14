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
      SUBROUTINE ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, ORDER )
C-----------------------------------------------------------------------
C
C - This subroutine finds the nodes of an element by looking at the
C - first edge and finding the opposite edge.
C - The nodes of these two edges should define the nodes of the element.
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
      INTEGER IPTR, NODS(4), NXT, LST
C
      LOGICAL ORDER
C
      INTEGER ED1, ED2, ED5, NXT1, NXT2, J
C
      NODS(1) = 0
      NODS(2) = 0
      NODS(3) = 0
      NODS(4) = 0
C
      LST  = BIGLST(LSTBIG,IPTR)
      NXT1 = BIGLST(NXTBIG,IPTR)
      ED1  = BIGLST(1,NXT1)
      ED2  = BIGLST(2,NXT1)
      NXT2 = BIGLST(NXTBIG,NXT1)
      ED5  = BIGLST(1,NXT2)
      NXT  = BIGLST(NXTBIG,NXT2)
C
      NODS(1) = BIGLST(1,ED1)
      NODS(2) = BIGLST(2,ED1)
      NODS(3) = BIGLST(1,ED5)
      NODS(4) = BIGLST(2,ED5)
C
      IF( .NOT. ORDER ) RETURN
C
      IF( BIGLST(1,ED2) .EQ. NODS(1) .OR.
     :    BIGLST(2,ED2) .EQ. NODS(1) ) THEN
         J = NODS(1)
         NODS(1) = NODS(2)
         NODS(2) = J
      ELSE IF( BIGLST(1,ED2) .NE. NODS(2) .AND.
     :         BIGLST(2,ED2) .NE. NODS(2) ) THEN
         PRINT*,'*** ELNODS: FAILED TO FIND COMMON NODE (1)'
         print*,'  elem (lst,nxt): ',iptr,lst,nxt
         print*,'  edges 1,2,5: ',ed1,ed2,ed5
         print*,'  nodes of ed1: ',(biglst(j,ed1),j=1,2)
         print*,'  nodes of ed2: ',(biglst(j,ed2),j=1,2)
         print*,'  nodes of ed5: ',(biglst(j,ed5),j=1,2)
         ed1 = biglst(3,nxt1)
         ed2 = biglst(4,nxt1)
         ed5 = biglst(2,nxt2)
         print*,'  edges 3,4,6: ',ed1,ed2,ed5
         print*,'  nodes of ed3: ',(biglst(j,ed1),j=1,2)
         print*,'  nodes of ed4: ',(biglst(j,ed2),j=1,2)
         print*,'  nodes of ed6: ',(biglst(j,ed5),j=1,2)
         IERR = -3895
         RETURN
      END IF
C
      IF( BIGLST(1,ED2) .EQ. NODS(4) .OR.
     :    BIGLST(2,ED2) .EQ. NODS(4) ) THEN
         J = NODS(4)
         NODS(4) = NODS(3)
         NODS(3) = J
      ELSE IF( BIGLST(1,ED2) .NE. NODS(3) .AND.
     :         BIGLST(2,ED2) .NE. NODS(3) ) THEN
         PRINT*,'*** ELNODS: FAILED TO FIND COMMON NODE (2)'
         print*,'  elem (lst,nxt): ',iptr,lst,nxt
         print*,'  edges 1,2,5: ',ed1,ed2,ed5
         print*,'  nodes of ed1: ',(biglst(j,ed1),j=1,2)
         print*,'  nodes of ed2: ',(biglst(j,ed2),j=1,2)
         print*,'  nodes of ed5: ',(biglst(j,ed5),j=1,2)
         ed1 = biglst(3,nxt1)
         ed2 = biglst(4,nxt1)
         ed5 = biglst(2,nxt2)
         print*,'  edges 3,4,6: ',ed1,ed2,ed5
         print*,'  nodes of ed3: ',(biglst(j,ed1),j=1,2)
         print*,'  nodes of ed4: ',(biglst(j,ed2),j=1,2)
         print*,'  nodes of ed6: ',(biglst(j,ed5),j=1,2)
         IERR = -3895
         RETURN
      END IF
C
      RETURN
      END
C
