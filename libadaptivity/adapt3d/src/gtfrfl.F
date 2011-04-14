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
      LOGICAL FUNCTION GTFRFL( IPTR, IBIT )
C-----------------------------------------------------------------------
C
C - This subroutine checks the flag for the element or edge in block IPTR
C - to see if bit IBIT (from 1 to 9) is set (.TRUE.) or unset (.FALSE.).
C
C - The bits are used in the following ways:
C - Bit 1: set if block is not empty (cannot be set using STFRFL)
C - Bit 2: used temporarily (by STCHFR) for looping through pairs of faces
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
      INCLUDE 'blkpwr.i'
C
      INTEGER IPTR, IBIT
C
      INTEGER IFLAG, IVAL, NXT, LST
C
      GTFRFL = .FALSE.
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IBIT .LE. 0 .OR. IBIT .GE. 10 ) THEN
         IERR = -10
         PRINT*,'*** GTFRFL: INVALID BIT NUMBER'
         RETURN
      END IF
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** GTFRFL: INVALID POINTER'
         RETURN
      END IF
C
      IFLAG = FRTLST(8,IPTR) / MAXNWF
C
      IF( IFLAG .EQ. 0 ) THEN
         IF( IBIT .EQ. 1 ) THEN
            GTFRFL = .FALSE.
         ELSE
            IERR = -2
            PRINT*,'*** GTFRFL: FACE BLOCK IS EMPTY'
         END IF
         RETURN
      ELSE IF( IFLAG .LT. 0 ) THEN
         IERR = -22
         PRINT*,'*** GTFRFL: FACE BLOCK IS INVALID'
         RETURN
      END IF
C
c      IVAL = 2**(IBIT-1)
      ival = powers(ibit)
C
      GTFRFL = ( MOD(IFLAG/IVAL,2) .EQ. 1 )
C
      RETURN
      END
C
