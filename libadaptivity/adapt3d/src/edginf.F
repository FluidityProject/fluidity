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
      SUBROUTINE EDGINF( BIGLST, NODLST, IPTR, NODS,
     :                   INEWED, IFLAG, NXT, LST )
C-----------------------------------------------------------------------
C
C - This subroutine returns the nodes list for the edge in BIGLST pointed
C - to by IPTR, as well as its position in NEWEDG (if it is in there) and
C - pointers to the next and last data in BIGLST.
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
      INTEGER NODS(2), IPTR, INEWED, IFLAG, NXT, LST
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** EDGINF: INVALID POINTER'
         print*,iptr
         RETURN
      END IF
C
      NODS(1) = BIGLST(1,IPTR)
      NODS(2) = BIGLST(2,IPTR)
      INEWED  = BIGLST(3,IPTR)
      IFLAG   = BIGLST(4,IPTR)
C
      IF( IFLAG .LT. 0 ) THEN
         IFLAG = -IFLAG
      ELSE IF( IFLAG .GT. 0 ) THEN
         IERR = -2
         PRINT*,'*** EDGINF: POINTER POINTS TO AN ELEMENT BLOCK'
         print*,iptr,biglst(4,iptr)
      ELSE
         IERR = -3
         PRINT*,'*** EDGINF: BLOCK IS EMPTY'
         print*,iptr,biglst(4,iptr)
      END IF
C
      LST = BIGLST(LSTBIG,IPTR)
      NXT = BIGLST(NXTBIG,IPTR)
C
      RETURN
      END
C


