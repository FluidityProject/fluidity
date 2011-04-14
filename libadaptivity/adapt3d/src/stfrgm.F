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
      SUBROUTINE STFRGM( BIGLST, NODLST )
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
      INCLUDE 'blkfrt.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IPTR, CONFC(3), EDGLST(3), IELM, INEW, IFLAG, NXT, LST,
     :        N1, N2, N3, CONFC2(3), EDGLS2(3), NXT2, LST2, I
C
      REAL NX1, NY1, NZ1, NX2, NY2, NZ2, DOT

C
      LOGICAL GTFRFL
C
c      PRINT*,' '
c      PRINT*,'Into STFRGM... ',ierr
C
      IF( IERR .NE. 0 ) RETURN
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
         PRINT*,'*** STFRGM: FRTLST DOES NOT END CORRECTLY'
         RETURN
      END IF
C
      IPTR = STTFRT
C
  10  CONTINUE
C
      CALL FACINF( IPTR, CONFC, EDGLST, IELM, INEW, IFLAG,
     :             NXT, LST )
C
      CALL CMNNOD( BIGLST, NODLST, EDGLST(1), EDGLST(2), N1, N2, N3 )
      CALL FACNRM( BIGLST, NODLST, N1, N2, N3, NX1, NY1, NZ1 )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( GTFRFL(CONFC(1),2) ) THEN
C
         CALL FACINF( CONFC(1), CONFC2, EDGLS2, IELM, INEW, IFLAG,
     :                NXT2, LST2 )
         CALL CMNNOD( BIGLST, NODLST, EDGLS2(1), EDGLS2(2), N1, N2, N3 )
         CALL FACNRM( BIGLST, NODLST, N1, N2, N3, NX2, NY2, NZ2 )
C
         IF( IERR .NE. 0 ) RETURN
C
         DOT = NX1*NX2 + NY1*NY2 + NZ1*NZ2
C
c         print*,'   edge,dot ',edglst(1),dot
         CALL SETFLG( BIGLST, NODLST, EDGLST(1), 2, DOT.LT.0.99999 )
         CALL SETFLG( BIGLST, NODLST, EDGLST(1), 3, .FALSE. )
C
      END IF
C
      IF( GTFRFL(CONFC(2),2) ) THEN
C
         CALL FACINF( CONFC(2), CONFC2, EDGLS2, IELM, INEW, IFLAG,
     :                NXT2, LST2 )
         CALL CMNNOD( BIGLST, NODLST, EDGLS2(1), EDGLS2(2), N1, N2, N3 )
         CALL FACNRM( BIGLST, NODLST, N1, N2, N3, NX2, NY2, NZ2 )
C
         IF( IERR .NE. 0 ) RETURN
C
         DOT = NX1*NX2 + NY1*NY2 + NZ1*NZ2
C
c         print*,'   edge,dot ',edglst(1),dot
         CALL SETFLG( BIGLST, NODLST, EDGLST(2), 2, DOT.LT.0.99999 )
         CALL SETFLG( BIGLST, NODLST, EDGLST(2), 3, .FALSE. )
C
      END IF
C
      IF( GTFRFL(CONFC(3),2) ) THEN
C
         CALL FACINF( CONFC(3), CONFC2, EDGLS2, IELM, INEW, IFLAG,
     :                NXT2, LST2 )
         CALL CMNNOD( BIGLST, NODLST, EDGLS2(1), EDGLS2(2), N1, N2, N3 )
         CALL FACNRM( BIGLST, NODLST, N1, N2, N3, NX2, NY2, NZ2 )
C
         IF( IERR .NE. 0 ) RETURN
C
         DOT = NX1*NX2 + NY1*NY2 + NZ1*NZ2
C
c         print*,'   edge,dot ',edglst(1),dot
         CALL SETFLG( BIGLST, NODLST, EDGLST(3), 2, DOT.LT.0.99999 )
         CALL SETFLG( BIGLST, NODLST, EDGLST(3), 3, .FALSE. )
C
      END IF
C
  30  CALL STFRFL( IPTR, 2, .TRUE. )
C
      IPTR = NXT
C
      IF( GTFRFL(IPTR,1) ) GOTO 10
C
c      PRINT*,'Leaving STFRGM'
C
      RETURN
      END
C
