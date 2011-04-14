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
      SUBROUTINE RMVFAC( IPTR )
C-----------------------------------------------------------------------
C
C - This subroutine removes a face from the front list - and also removes
C - that face from connected faces' lists of connected faces (to be filled
C - in by new face(s) which should be created immediately).
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blkfrt.i'
C
      INCLUDE 'blknwf.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IPTR
C
      INTEGER FACE1, FACE2, FACE3, CONFAC(3), EDGLST(3), IELM, INEWFC,
     :        NXT, LST, IFLAG
C
c      print*,' '
c      print*,'Into RMVFAC...',ierr
c      print*,'     IPTR:',iptr
c      print*,'     sttfrt,emtfrt,endfrt',sttfrt,emtfrt,endfrt
c      print*,'     nfaces,nnewfc:',nfaces,nnewfc
C
      IF( IERR .NE. 0 ) RETURN
C
      CALL FACINF( IPTR, CONFAC, EDGLST, IELM, INEWFC, IFLAG,
     :             NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      FACE1 = CONFAC(1)
      FACE2 = CONFAC(2)
      FACE3 = CONFAC(3)
C
      IF( FACE1 .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** RMVFAC: FOUND INVALID CONNECTED FACE 1'
         RETURN
      END IF
C
      IF( FACE1 .NE. 0 ) THEN
c         print*,'   removing connection of face ',face1
         IF( FRTLST(1,FACE1) .EQ. IPTR ) THEN
            FRTLST(1,FACE1) = 0
         ELSE IF( FRTLST(2,FACE1) .EQ. IPTR ) THEN
            FRTLST(2,FACE1) = 0
         ELSE IF( FRTLST(3,FACE1) .EQ. IPTR ) THEN
            FRTLST(3,FACE1) = 0
         END IF
      END IF
C
      IF( FACE2 .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** RMVFAC: FOUND INVALID CONNECTED FACE 2'
         RETURN
      END IF
C
      IF( FACE2 .NE. 0 ) THEN
c         print*,'   removing connection of face ',face2
         IF( FRTLST(1,FACE2) .EQ. IPTR ) THEN
            FRTLST(1,FACE2) = 0
         ELSE IF( FRTLST(2,FACE2) .EQ. IPTR ) THEN
            FRTLST(2,FACE2) = 0
         ELSE IF( FRTLST(3,FACE2) .EQ. IPTR ) THEN
            FRTLST(3,FACE2) = 0
         END IF
      END IF
C
      IF( FACE3 .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** RMVFAC: FOUND INVALID CONNECTED FACE 3'
         RETURN
      END IF
C
      IF( FACE3 .NE. 0 ) THEN
c         print*,'   removing connection of face ',face3
         IF( FRTLST(1,FACE3) .EQ. IPTR ) THEN
            FRTLST(1,FACE3) = 0
         ELSE IF( FRTLST(2,FACE3) .EQ. IPTR ) THEN
            FRTLST(2,FACE3) = 0
         ELSE IF( FRTLST(3,FACE3) .EQ. IPTR ) THEN
            FRTLST(3,FACE3) = 0
         END IF
      END IF
C
      FRTLST( 1, IPTR ) = 0
      FRTLST( 2, IPTR ) = 0
      FRTLST( 3, IPTR ) = 0
      FRTLST( 4, IPTR ) = 0
      FRTLST( 5, IPTR ) = 0
      FRTLST( 6, IPTR ) = 0
      FRTLST( 7, IPTR ) = 0
      FRTLST( 8, IPTR ) = 0
C
c      IF( INEWFC .GT. 0 ) THEN
c         NEWFAC(INEWFC) = NEWFAC(NNEWFC)
c         IFLAG = FRTLST(8,NEWFAC(INEWFAC)) - NNEWFC
c         FRTLST(8,NEWFAC(INEWFAC)) = IFLAG + INEWFAC
c         NNEWFC = NNEWFC - 1
c      END IF
C
      IF( IPTR .EQ. STTFRT ) STTFRT = NXT
C
      IF( LST .GT. 0 ) FRTLST( NXTFRT, LST ) = NXT
      IF( NXT .GT. 0 ) FRTLST( LSTFRT, NXT ) = LST
C
      IF( ENDFRT .EQ. IPTR ) ENDFRT = LST
C
      FRTLST( LSTFRT, IPTR   ) = ENDFRT
      FRTLST( NXTFRT, IPTR   ) = EMTFRT
C
      IF( EMTFRT .GT. 0 ) FRTLST( LSTFRT, EMTFRT ) = IPTR
      IF( ENDFRT .GT. 0 ) FRTLST( NXTFRT, ENDFRT ) = IPTR
C
      EMTFRT = IPTR
C
      NFACES = NFACES - 1
C
      IF( FRTLST(LSTFRT,EMTFRT) .NE. ENDFRT ) THEN
         IERR = -110
         PRINT*,'*** RMVFAC: INCONSISTENCY AT END OF FRTLST'
         print*,'lst,endfrt: ',FRTLST(LSTFRT,EMTFRT),ENDFRT
         RETURN
      END IF
C
c      print*,'Leaving RMVFAC'
c      print*,'     sttfrt,emtfrt,endfrt',sttfrt,emtfrt,endfrt
c      print*,'     nfaces,nnewfc:',nfaces,nnewfc
C
      RETURN
      END
C
