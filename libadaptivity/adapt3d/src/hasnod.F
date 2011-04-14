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
      INTEGER FUNCTION HASNOD( BIGLST, NODLST, IELM, INOD )
C-----------------------------------------------------------------------
C
C - This function checks to see if element IELM has node INOD.
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
      INTEGER IELM, INOD
C
      INTEGER NXT, LST, EDGLST(6), CONELM(4), INEW, IFLAG,
     :        NODS1(2), NODS2(2), I, OPEDG, IREG
C
      HASNOD = 0
C
      CALL ELMINF( BIGLST, NODLST, IELM, CONELM, EDGLST,
     :             INEW, IFLAG, IREG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      CALL EDGINF( BIGLST, NODLST, EDGLST(1), NODS1,
     :             INEW, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** HASNOD: GOT ERROR FROM EDGINF'
         print*,'   element,node: ',IELM,INOD
         print*,'   conelm: ',conelm
         print*,'   edglst: ',edglst
         RETURN
      END IF
C
      HASNOD = 1
      IF( NODS1(1) .EQ. INOD ) GOTO 100
      HASNOD = 2
      IF( NODS1(2) .EQ. INOD ) GOTO 100
C
      I = 1
C
  10  IF( I .LT. 6 ) THEN
         I = I + 1
         OPEDG = EDGLST(I)
         IF( BIGLST(1,OPEDG) .EQ. NODS1(1) ) GOTO 10
         IF( BIGLST(2,OPEDG) .EQ. NODS1(1) ) GOTO 10
         IF( BIGLST(1,OPEDG) .EQ. NODS1(2) ) GOTO 10
         IF( BIGLST(2,OPEDG) .EQ. NODS1(2) ) GOTO 10
      ELSE
         OPEDG = 0
         PRINT*,'*** HASNOD: COULD NOT FIND OPPOSITE EDGE'
         PRINT*,EDGLST
         IERR = -150
         RETURN
      END IF
C
      CALL EDGINF( BIGLST, NODLST, OPEDG, NODS2,
     :             INEW, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** HASNOD: GOT ERROR FROM 2nd EDGINF'
         print*,'   element,node,edge: ',IELM,INOD,OPEDG
         print*,'   conelm: ',conelm
         print*,'   edglst: ',edglst
         RETURN
      END IF
C
      HASNOD = 3
      IF( NODS2(1) .EQ. INOD ) GOTO 100
      HASNOD = 4
      IF( NODS2(2) .EQ. INOD ) GOTO 100
C
      HASNOD = 0
C
 100  CONTINUE
C
      RETURN
      END
C
