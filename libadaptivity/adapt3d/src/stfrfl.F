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
      SUBROUTINE STFRFL( IPTR, IBIT, VALUE )
C-----------------------------------------------------------------------
C
C - This subroutine sets/unsets bit IBIT (from 2 to 9) of the flag for
C - the face in block IPTR to VALUE (.TRUE. or .FALSE.)
C
C - The bits are used in the following ways:
C - Bit 1: set if block is not empty - not allowed by this subroutine
C - Bit 2: used temporarily (by STCHFR) for looping through pairs of faces
C - Bit 3: used temporarily (by NODFCS) for find faces round a node
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
      LOGICAL VALUE
C
      INTEGER IFLAG, IVAL, NXT, INEWFC
C
      LOGICAL GETFLG
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IBIT .LE. 1 .OR. IBIT .GE. 10 ) THEN
         IERR = -10
         PRINT*,'*** STFRFL: INVALID BIT NUMBER'
         RETURN
      END IF
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** STFRFL: INVALID POINTER'
         RETURN
      END IF
C
      IFLAG = FRTLST(8,IPTR) / MAXNWF
C
      IF( IFLAG .EQ. 0 ) THEN
         IERR = -2
         PRINT*,'*** STFRFL: FACE BLOCK IS EMPTY'
         RETURN
      END IF
C
      INEWFC = MOD(FRTLST(8,IPTR),MAXNWF)
C
c      IVAL = 2**(IBIT-1)
      ival = powers(ibit)
C
      GETFLG = ( MOD(IFLAG/IVAL,2) .EQ. 1 )
C
      IF( GETFLG .NEQV. VALUE ) THEN
C
         IF( VALUE ) THEN
            IFLAG = IFLAG + IVAL
         ELSE
            IFLAG = IFLAG - IVAL
         END IF
C
         FRTLST(8,IPTR) = IFLAG*MAXNWF + INEWFC
C
      END IF
C
      RETURN
      END
C
