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
      SUBROUTINE NDCNCT( BIGLST, NODLST, INOD, IELM,
     :                   ELPTRS, EDPTRS, NDPTRS, NUMEDG, NUMELM )
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
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      INTEGER ELPTRS(MAXEDE), EDPTRS(MAXEDE), NDPTRS(MAXEDE), NUMELM,
     :        NUMEDG, INOD, IELM
C
      INTEGER ELMADD, EDGADD, CONELM(4), EDGLST(6), INEW, IFLAG, NXT,
     :        LST, HASNOD, I, IREG
C
      LOGICAL GETFLG
C
      CALL ADAPT_SECOND( TIMSTT )
C
      NUMELM = 0
      NUMEDG = 0
      ELMADD = 0
C
      IF( IERR .NE. 0 ) RETURN
C
      NUMELM = NUMELM + 1
      ELPTRS(NUMELM) = IELM
      CALL SETFLG( BIGLST, NODLST, ELPTRS(NUMELM), 4, .TRUE. )
      ELMADD = ELMADD + 1
C
 10   CALL ELMINF( BIGLST, NODLST, ELPTRS(NUMELM), CONELM, EDGLST,
     :             INEW, IFLAG, IREG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      DO I = 1, 6
C
         IF( .NOT. GETFLG(BIGLST,NODLST,EDGLST(I),5) ) THEN
            IF( BIGLST(1,EDGLST(I)) .EQ. INOD ) THEN
               NUMEDG = NUMEDG + 1
               EDPTRS(NUMEDG) = EDGLST(I)
               NDPTRS(NUMEDG) = BIGLST(2,EDGLST(I))
               CALL SETFLG( BIGLST, NODLST, EDGLST(I), 5, .TRUE. )
            ELSE IF( BIGLST(2,EDGLST(I)) .EQ. INOD ) THEN
               NUMEDG = NUMEDG + 1
               EDPTRS(NUMEDG) = EDGLST(I)
               NDPTRS(NUMEDG) = BIGLST(1,EDGLST(I))
               CALL SETFLG( BIGLST, NODLST, EDGLST(I), 5, .TRUE. )
            END IF
         END IF
C
      END DO
C
      IF( IERR .NE. 0 ) RETURN
C
      DO I = 1, 4
C
         IF( CONELM(I) .GT. 0 ) THEN
            IF( .NOT. GETFLG(BIGLST,NODLST,CONELM(I),4) ) THEN
               IF( HASNOD(BIGLST,NODLST,CONELM(I),INOD) .NE. 0 ) THEN
                  ELMADD = ELMADD + 1
                  ELPTRS(ELMADD) = CONELM(I)
                  CALL SETFLG( BIGLST, NODLST, CONELM(I), 4, .TRUE. )
               END IF
            END IF
         END IF
C
      END DO
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( NUMELM .LT. ELMADD ) THEN
         NUMELM = NUMELM + 1
         GOTO 10
      END IF
C
      DO I = 1, NUMEDG
         CALL SETFLG( BIGLST, NODLST, EDPTRS(I), 5, .FALSE. )
      END DO
C
      DO I = 1, NUMELM
         CALL SETFLG( BIGLST, NODLST, ELPTRS(I), 4, .FALSE. )
      END DO
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMNCN = TIMNCN + TIMEND - TIMSTT
      TOTNCN = TOTNCN + 1
C
      RETURN
      END
C
