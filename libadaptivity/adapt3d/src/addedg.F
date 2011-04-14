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
      SUBROUTINE ADDEDG( BIGLST, NODLST, NODS, IPTR, INEWED )
C-----------------------------------------------------------------------
C
C - This subroutine adds an edge into biglst.
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
      INTEGER NODS(2), IPTR, INEWED
C
      INTEGER NXT, LST, JFLAG
c
c      PRINT*,' '
c      print*,'Into ADDEDG.. ',ierr
c      print*,'     NODS:',NODS
c      PRINT*,'   INEWED:',INEWED
C
      IF( IERR .NE. 0 ) RETURN
C
      IPTR = EMTBIG
      IF( IPTR .LE. 0 ) THEN
         IERR = -4
         PRINT*,'*** FAILED TO ADD EDGE TO BIGLST - FULL'
         RETURN
      END IF
C
      JFLAG = BIGLST(4,IPTR)
C
      IF( JFLAG .NE. 0 ) THEN
         PRINT*,'*** FAILED TO ADD EDGE TO BIGLST ',
     :          '- BLOCK NOT EMPTY!'
         IERR = -2
         RETURN
      END IF
C
      IF( BIGLST(LSTBIG,IPTR) .NE. ENDBIG ) THEN
         IERR = -110
         PRINT*,'*** ADDEDG: INCONSISTENCY AT END OF BIGLST'
         RETURN
      END IF
C
      IF( ENDBIG .GT. 0 ) THEN
         IF( BIGLST(NXTBIG,ENDBIG) .NE. IPTR ) THEN
            IERR = -110
            PRINT*,'*** ADDEDG: INCONSISTENCY NEAR END OF BIGLST'
            RETURN
         END IF
      END IF
C
      BIGLST(1,IPTR) = NODS(1)
      BIGLST(2,IPTR) = NODS(2)
C
c      IF( INEWED .GE. 0 ) THEN
c         NNEWED = NNEWED + 1
c         INEWED = NNEWED
c         NEWEDG(INEWED) = IPTR
c      ELSE
         INEWED = 0
c      END IF
C
      BIGLST(3,IPTR) = INEWED
      BIGLST(4,IPTR) = -1
C
C      BIGLST(LSTBIG,IPTR) = ENDBIG
C      IF( ENDBIG .GT. 0 ) BIGLST(NXTBIG,ENDBIG) = IPTR
      ENDBIG = IPTR
      EMTBIG = BIGLST(NXTBIG,IPTR)
C
      BIGLST(7,IPTR) = 0
C
      NEDGES = NEDGES + 1
c
c      if( iptr .eq. 32454 ) then
c       PRINT*,' '
c       print*,'Into ADDEDG.. ',ierr
c       print*,'     NODS:',NODS
c       print*,'Leaving ADDEDG... '
c       print*,'  IPTR,INEW:',IPTR,INEWED
c       print*,' '
c       pause
c      end if
C
      RETURN
      END
C
