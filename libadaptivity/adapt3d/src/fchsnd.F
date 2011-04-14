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
      INTEGER FUNCTION FCHSND( BIGLST, NODLST, IFAC, INOD )
C-----------------------------------------------------------------------
C
C - This function checks to see if face IFAC has node INOD.
C - If so it returns .TRUE. otherwise it returns .FALSE.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkfrt.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IFAC, INOD
C
      INTEGER NXT, LST, EDGLST(3), CONFAC(3), IELM, INEW, IFLAG,
     :        NODS1(2), NODS2(2), NODS3(2)
C
      FCHSND = 0
C
      CALL FACINF( IFAC, CONFAC, EDGLST, IELM, INEW, IFLAG, NXT, LST )
C
      IF( IERR .LT. 0 ) RETURN
C
      IF( INOD .LE. 0 ) THEN
         IERR = -2
         PRINT*,'*** FCHSND: INVALID NODE POINTER GIVEN'
         RETURN
      END IF
C
      CALL EDGINF( BIGLST, NODLST, EDGLST(1), NODS1,
     :             INEW, IFLAG, NXT, LST )
      CALL EDGINF( BIGLST, NODLST, EDGLST(2), NODS2,
     :             INEW, IFLAG, NXT, LST )
C      CALL EDGINF( BIGLST, NODLST, EDGLST(3), NODS3, INEW, IFLAG, NXT, LST )
C
      FCHSND = 1
      IF( NODS1(1) .EQ. INOD ) GOTO 100
      FCHSND = 2
      IF( NODS1(2) .EQ. INOD ) GOTO 100
      FCHSND = 3
      IF( NODS2(1) .EQ. INOD ) GOTO 100
      IF( NODS2(2) .EQ. INOD ) GOTO 100
C
      FCHSND = 0
C
      RETURN
C
 100  CONTINUE
C
      RETURN
      END
C
