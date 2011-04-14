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
      SUBROUTINE CNCTFC( FACE1, FACE2 )
C-----------------------------------------------------------------------
C
C - This subroutine adds FACE2 to the list of connected faces for FACE1,
C - and vice versa. If they already exist in the lists then it does nothing.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blkfrt.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER FACE1, FACE2
c
c      print*,'Into CNCTFC: ',face1,face2,ierr
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( FACE1 .LE. 0 ) THEN
         PRINT*,'*** CNCTFC: FACE1 IS NOT A VALID POINTER'
         IERR = -1
         RETURN
      END IF
C
      IF( FACE2 .LE. 0 ) THEN
         PRINT*,'*** CNCTFC: FACE2 IS NOT A VALID POINTER'
         IERR = -1
         RETURN
      END IF
C
      IF( FRTLST(1,FACE1) .EQ. FACE2 ) GOTO 100
      IF( FRTLST(2,FACE1) .EQ. FACE2 ) GOTO 100
      IF( FRTLST(3,FACE1) .EQ. FACE2 ) GOTO 100
C
      IF( FRTLST(1,FACE1) .EQ. 0 ) THEN
         FRTLST(1,FACE1) = FACE2
      ELSE IF( FRTLST(2,FACE1) .EQ. 0 ) THEN
         FRTLST(2,FACE1) = FACE2
      ELSE IF( FRTLST(3,FACE1) .EQ. 0 ) THEN
         FRTLST(3,FACE1) = FACE2
      ELSE
         PRINT*,'*** CNCTFC: FACE1 IS FULLY CONNECTED'
         IERR = -26
      END IF
C
 100  CONTINUE
C
      IF( FRTLST(1,FACE2) .EQ. FACE1 ) GOTO 200
      IF( FRTLST(2,FACE2) .EQ. FACE1 ) GOTO 200
      IF( FRTLST(3,FACE2) .EQ. FACE1 ) GOTO 200
C
      IF( FRTLST(1,FACE2) .EQ. 0 ) THEN
         FRTLST(1,FACE2) = FACE1
      ELSE IF( FRTLST(2,FACE2) .EQ. 0 ) THEN
         FRTLST(2,FACE2) = FACE1
      ELSE IF( FRTLST(3,FACE2) .EQ. 0 ) THEN
         FRTLST(3,FACE2) = FACE1
      ELSE
         PRINT*,'*** CNCTFC: FACE2 IS FULLY CONNECTED'
         IERR = -26
      END IF
C
 200  CONTINUE
c
c      print*,'Leaving CNCTFC'
C
      RETURN
      END
C
