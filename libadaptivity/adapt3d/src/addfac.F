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
      SUBROUTINE ADDFAC( CONFAC, CNCT, EDGLST, IELM, IPTR, INEWFC )
C-----------------------------------------------------------------------
C     
C     
C     
C     
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
      INTEGER CONFAC(3), EDGLST(3), IPTR, INEWFC, IELM
C     
      LOGICAL CNCT
C     
      INTEGER NXT, LST
c     
c      PRINT*,' '
c      print*,'Into ADDFAC.. ',ierr
c      print*,'    CONFAC:',CONFAC
c      print*,'    EDGLST:',EDGLST
c      print*,' IELM,INEW:',IELM,INEWFC
c      print*,'     sttfrt,emtfrt,endfrt',sttfrt,emtfrt,endfrt
c      print*,'     nfaces,nnewfc:',nfaces,nnewfc
C
      IF( IERR .NE. 0 ) then
         RETURN
      END IF
C
      IPTR = EMTFRT
      IF( IPTR .LE. 0 ) THEN
         IERR = -4
         PRINT*,'*** FAILED TO ADD FACE TO FRTLST - FULL'
         RETURN
      END IF
C
      IF( FRTLST(LSTFRT,IPTR) .NE. ENDFRT ) THEN
         IERR = -110
         PRINT*,'*** ADDFAC: INCONSISTENCY AT END OF FRTLST'
         print*,'lst,endfrt: ',FRTLST(LSTFRT,IPTR),ENDFRT
         RETURN
      END IF
C
      IF( ENDFRT .GT. 0 ) THEN
         IF( FRTLST(NXTFRT,ENDFRT) .NE. IPTR ) THEN
            IERR = -110
            PRINT*,'*** ADDFAC: INCONSISTENCY NEAR END OF FRTLST'
            RETURN
         END IF
      END IF
C
      FRTLST(1,IPTR) = CONFAC(1)
      FRTLST(2,IPTR) = CONFAC(2)
      FRTLST(3,IPTR) = CONFAC(3)
C
C - Need to tell connected faces this new face is connected to them.
C - Since rmvfac should have been called already (to remove the faces
C - this one will replace on the front) there should be a zero in each
C - of the connected face's list of connected faces ready to be filled
C - by this new face.
C
      IF( CNCT ) THEN
         IF( CONFAC(1) .GT. 0 )
     :       CALL CNCTFC( CONFAC(1), IPTR )
         IF( CONFAC(2) .GT. 0 )
     :       CALL CNCTFC( CONFAC(2), IPTR )
         IF( CONFAC(3) .GT. 0 )
     :       CALL CNCTFC( CONFAC(3), IPTR )
      END IF
C
      IF( IERR .NE. 0 ) RETURN
C
      FRTLST(4,IPTR) = EDGLST(1)
      FRTLST(5,IPTR) = EDGLST(2)
      FRTLST(6,IPTR) = EDGLST(3)
C
c     IF( INEWFC .GE. 0 ) THEN
c         NNEWFC = NNEWFC + 1
c         INEWFC = NNEWFC
c         NEWFAC(INEWFC) = IPTR
c      ELSE
         INEWFC = 0
c      END IF
C
      FRTLST(7,IPTR) = IELM
      FRTLST(8,IPTR) = INEWFC + MAXNWF
C
C      FRTLST(LSTFRT,IPTR) = ENDFRT
C      IF( ENDFRT .GT. 0 ) FRTLST(NXTFRT,ENDFRT) = IPTR
      ENDFRT = IPTR
      EMTFRT = FRTLST(NXTFRT,IPTR)
C
      NFACES = NFACES + 1
c     
c     print*,'Leaving ADDFAC...'
c     print*,'   IPTR,INEW:',IPTR,INEWFC
c     print*,'     sttfrt,emtfrt,endfrt',sttfrt,emtfrt,endfrt
c     print*,'     nfaces,nnewfc:',nfaces,nnewfc
C     
      RETURN
      END SUBROUTINE ADDFAC
