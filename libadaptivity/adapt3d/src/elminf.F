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
      SUBROUTINE ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                   INEWEL, IFLAG, IREG, NXT, LST )
C-----------------------------------------------------------------------
C
C - This subroutine returns the element connectivity list and edge list for
C - the element in BIGLST pointed to by IPTR, as well as its position in
C - NEWELE (if it is in there) and pointers to the next and last data in BIGLST.
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
      INTEGER CONELM(4), EDGLST(6), IPTR, INEWEL, IFLAG, IREG, NXT, LST
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** ELMINF: INVALID POINTER'
         RETURN
      END IF
C
      IF( BIGLST(4,IPTR) .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** ELMINF: FIRST BLOCK IS AN EDGE BLOCK'
         print*,'    element: ',iptr,biglst(4,iptr)
         GOTO 900
      ELSE IF( BIGLST(4,IPTR) .EQ. 0 ) THEN
         IERR = -2
         PRINT*,'*** ELMINF: FIRST BLOCK IS EMPTY'
         print*,'    element: ',iptr,biglst(4,iptr)
         GOTO 900
      END IF
C
      CONELM(4) = BIGLST(4,IPTR)/2
      CONELM(1) = BIGLST(1,IPTR)
      CONELM(2) = BIGLST(2,IPTR)
      CONELM(3) = BIGLST(3,IPTR)
C
      NXT = BIGLST(NXTBIG,IPTR)
C
      IF( NXT .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** ELMINF: INVALID POINTER IN FIRST BLOCK'
         RETURN
      END IF
C
      EDGLST(1) = BIGLST(1,NXT)
      EDGLST(2) = BIGLST(2,NXT)
      EDGLST(3) = BIGLST(3,NXT)
      EDGLST(4) = BIGLST(4,NXT)
      IF( BIGLST(4,NXT) .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** ELMINF: SECOND BLOCK IS AN EDGE BLOCK'
         GOTO 900
      ELSE IF( BIGLST(4,NXT) .EQ. 0 ) THEN
         IERR = -2
         PRINT*,'*** ELMINF: SECOND BLOCK IS EMPTY'
         GOTO 900
      END IF
C
      NXT = BIGLST(NXTBIG,NXT)
C
      IF( NXT .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** ELMINF: INVALID POINTER IN SECOND BLOCK'
         RETURN
      END IF
C
      EDGLST(5) = BIGLST(1,NXT)
      EDGLST(6) = BIGLST(2,NXT)
      INEWEL    = BIGLST(3,NXT)
      IFLAG     = BIGLST(4,NXT)
      IREG      = BIGLST(7,NXT)
C
      IF( IFLAG .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** ELMINF: THIRD BLOCK IS AN EDGE BLOCK'
      ELSE IF( IFLAG .EQ. 0 ) THEN
         IERR = -2
         PRINT*,'*** ELMINF: THIRD BLOCK IS EMPTY'
      END IF
C
 900  LST = BIGLST(LSTBIG,IPTR)
      NXT = BIGLST(NXTBIG,NXT)
C
      RETURN
      END
C


