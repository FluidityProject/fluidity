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
      SUBROUTINE STNDFL( BIGLST, NODLST, IPTR, IBIT, VALUE )
C-----------------------------------------------------------------------
C
C - This subroutine sets/unsets bit IBIT (from 2 to 17) of the flag for
C - the node in block IPTR to VALUE (.TRUE. or .FALSE.)
C
C - Bits 1 to 10 are reserved for internal use...
C
C - These bits are used in the following ways:
C - Bit 1: set if block is not empty - not allowed by this subroutine
C - Bit 2: set if the node is an essential geometry node
C - Bit 3: set if the node is internal (should be overridden by bit 2)
C - Bit 4: set if the node splits a previous geometry edge
C - Bit 5: used temporarily (by ??????) to loop through nodes using edges
C - Bit 6: used temporarily (by ??????) to loop through nodes using elements
C - Bit 7: used temporarily (by EDGTST) during edge shrink evaluation
C - Bit 8: used to flag nodes not to be touched during adapt
C - Bit 9: used temporarily (by MKADPT) to flag halo nodes in a parallel problem
C
C-----------------------------------------------------------------------
C
      use flag_handling
C
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blkpwr.i'
C
      INTEGER IPTR, IBIT
C
      LOGICAL VALUE
C
      call set_flag(nodlst, iptr, ibit, value)
C
      end subroutine stndfl
C
#if 0
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blkpwr.i'
C
      INTEGER IPTR, IBIT
C
      LOGICAL VALUE
C
      INTEGER IFLAG, IVAL, NXT
C
      LOGICAL GETFLG
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IBIT .LE. 1 .OR. IBIT .GE. 18 ) THEN
         IERR = -10
         PRINT*,'*** STNDFL: INVALID BIT NUMBER'
         print*,iptr,ibit
         RETURN
      END IF
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** STNDFL: INVALID POINTER'
         print*,iptr,ibit
         RETURN
      END IF
C
      IFLAG = INT(NODLST(6,IPTR))
C
      IF( IFLAG .EQ. 0 ) THEN
         IERR = -2
         PRINT*,'*** STNDFL: NODE BLOCK IS EMPTY'
         print*,iptr,ibit
         RETURN
      END IF
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
         NODLST(6,IPTR) = FLOAT(IFLAG)
C
      END IF
C
      RETURN
      END
C
#endif
