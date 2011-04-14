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
      SUBROUTINE STCHFR( BIGLST, NODLST, IREG, MINDOT, MINREL, 
     :                   AWSDOT, AWSREL, CHNGED )
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
      INCLUDE 'blkfrt.i'
C
      INCLUDE 'blknwf.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IREG
C
      REAL MINDOT, MINREL, AWSDOT, AWSREL
C
      LOGICAL CHNGED, DUNONE
C
      INTEGER I, J, K, IPTR, NXT, LST, INEW, IFLAG, CONFC(3),
     :        EDGLST(3), IELM, BSTFC1, BSTFC2, BSTFC3
C
      REAL BSTDOT, BSTREL, THSDOT, THSREL
C
      LOGICAL GTFRFL
c
      print*,' '
      print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      print*,' '
      print*,'Into STCHFR: ',ierr,nnewfc
C
      IF( IERR .NE. 0 ) RETURN
C
  10  CONTINUE
C
      BSTDOT = MINDOT
      BSTREL = MINREL
      BSTFC1 = 0
      BSTFC2 = 0
C
      DUNONE = .FALSE.
C
      IPTR = STTFRT
C
      DO I = 1, NFACES
         CALL STFRFL( IPTR, 2, .FALSE. )
         IPTR = FRTLST(NXTFRT,IPTR)
      END DO
C
      IF( IPTR .NE. EMTFRT ) THEN
         IERR = -90
         PRINT*,'*** STCHFR: FRTLST DOES NOT END CORRECTLY'
         RETURN
      END IF
C
      IPTR = STTFRT
C
  20  CONTINUE
C
      IF( IS3DMS ) THEN
         IF( NFACES .EQ. 4 ) THEN
            GOTO 900
         ELSE IF( NFACES .LT. 4 ) THEN
            IERR = -100
            PRINT*,'*** STCHFR: GONE BELOW FOUR FACES'
            RETURN
         END IF
      ELSE IF( NFACES .EQ. 3 ) THEN
         GOTO 900
      ELSE IF( NFACES .LT. 3 ) THEN
         IERR = -100
         PRINT*,'*** STCHFR: GONE BELOW THREE EDGES'
         RETURN
      END IF
C
      CALL FACINF( IPTR, CONFC, EDGLST, IELM, INEW, IFLAG,
     :             NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( GTFRFL(CONFC(1),2) ) THEN
         CALL DOFACS( BIGLST, NODLST, IPTR, CONFC, EDGLST, IELM, NXT, 1,
     :                THSDOT, THSREL, MINDOT, MINREL, .FALSE., IREG )
         IF( IERR .NE. 0 ) RETURN
         IF( (THSDOT .GT. BSTDOT) .AND. (THSREL .GT. BSTREL) ) THEN
            BSTDOT = THSDOT
            BSTREL = THSREL
            BSTFC1 = IPTR
            BSTFC2 = 1
         END IF
C         IF( IPTR .EQ. 0 ) THEN
C            CHNGED = .TRUE.
C            GOTO 30
C         END IF
      END IF
C
      IF( GTFRFL(CONFC(2),2) ) THEN
         CALL DOFACS( BIGLST, NODLST, IPTR, CONFC, EDGLST, IELM, NXT, 2,
     :                THSDOT, THSREL, MINDOT, MINREL, .FALSE., IREG )
         IF( IERR .NE. 0 ) RETURN
         IF( (THSDOT .GT. BSTDOT) .AND. (THSREL .GT. BSTREL) ) THEN
            BSTDOT = THSDOT
            BSTREL = THSREL
            BSTFC1 = IPTR
            BSTFC2 = 2
         END IF
C         IF( IPTR .EQ. 0 ) THEN
C            CHNGED = .TRUE.
C            GOTO 30
C         END IF
      END IF
C
      IF( GTFRFL(CONFC(3),2) ) THEN
         CALL DOFACS( BIGLST, NODLST, IPTR, CONFC, EDGLST, IELM, NXT, 3,
     :                THSDOT, THSREL, MINDOT, MINREL, .FALSE., IREG )
         IF( IERR .NE. 0 ) RETURN
         IF( (THSDOT .GT. BSTDOT) .AND. (THSREL .GT. BSTREL) ) THEN
            BSTDOT = THSDOT
            BSTREL = THSREL
            BSTFC1 = IPTR
            BSTFC2 = 3
         END IF
C         IF( IPTR .EQ. 0 ) CHNGED = .TRUE.
      END IF
C
  30  IF( IPTR .NE. 0 ) CALL STFRFL( IPTR, 2, .TRUE. )
C
      IF( BSTDOT .GT. AWSDOT .AND. BSTREL .GT. AWSREL ) GOTO 100
C
      IPTR = NXT
C
      IF( GTFRFL(IPTR,1) ) GOTO 20
C
 100  IF( BSTFC1 .GT. 0 ) THEN
C
         IPTR = BSTFC1
C
         CALL FACINF( IPTR, CONFC, EDGLST, IELM, INEW, IFLAG,
     :                NXT, LST )
C
         PRINT*,' '
         PRINT*,'* FOUND BEST FACES: ',BSTFC1,CONFC(BSTFC2)
         PRINT*,'   BSTDOT,BSTREL: ',BSTDOT,BSTREL
c         PRINT*,' '
C
c         PAUSE
C
         CALL DOFACS( BIGLST, NODLST, IPTR, CONFC, EDGLST, IELM, NXT,
     :                BSTFC2, THSDOT, THSREL, BSTDOT/2, BSTREL/2, 
     :               .TRUE., IREG )
C
         IF( IPTR .NE. 0 ) THEN
            PRINT*,'*** STCHFR: DID NOT CHANGE FRONT'
            IERR = -101
            RETURN
         END IF
C
         DUNONE = .TRUE.
C
         IPTR   = BSTFC1
         BSTFC1 = 0
         BSTFC2 = 0
         BSTDOT = MINDOT
         BSTREL = MINREL
C
      END IF
C
      IPTR = NXT
C
c      if( nelems .lt. 7 ) then
      IF( GTFRFL(IPTR,1) ) GOTO 20
c      end if
C
      if( dunone ) then
         PRINT*,'CHECKING FRONT AGAIN...'
c         GOTO 10
      end if
C
      PRINT*,'Reduced to ',NFACES,' faces'
C
      GOTO 1000
C
 900  CONTINUE
C
      PRINT*,'Successfully reduced to minimum number of faces'
c
 1000 CONTINUE
C
      IPTR = STTFRT
C
      DO I = 1, NFACES
         CALL STFRFL( IPTR, 2, .FALSE. )
         IPTR = FRTLST(NXTFRT,IPTR)
      END DO
C
      IF( IPTR .NE. EMTFRT ) THEN
         IERR = -90
         PRINT*,'*** STCHFR: FRTLST DOES NOT END CORRECTLY'
         RETURN
      END IF
C
      print*,'Leaving STCHFR'
      print*,' '
      print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      print*,' '
C
      RETURN
      END
C
