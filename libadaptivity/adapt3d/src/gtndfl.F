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
      LOGICAL FUNCTION GTNDFL( BIGLST, NODLST, IPTR, IBIT )
C-----------------------------------------------------------------------
C
C - This subroutine checks the flag for the node in block IPTR
C - to see if bit IBIT (from 1 to 17) is set (.TRUE.) or unset (.FALSE.).
C
C - Bits 1 to 10 are reserved for internal use...
C
C - These bits are used in the following ways:
C - Bit 1: set if block is not empty (cannot be set using STNDFL)
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
      INTEGER IPTR, IBIT
C
      gtndfl = get_flag(nodlst, iptr, ibit)
C
      end function gtndfl
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
      INTEGER IFLAG, IVAL, NXT, LST
C
      GTNDFL = .FALSE.
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IBIT .LE. 0 .OR. IBIT .GE. 18 ) THEN
         IERR = -10
         PRINT*,'*** GTNDFL: INVALID BIT NUMBER'
         print*,iptr,ibit
         print*,(nodlst(ival,iptr),ival=1,3)
         print*,(int(nodlst(ival,iptr)),ival=4,5)
         RETURN
      END IF
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** GTNDFL: INVALID POINTER'
         print*,iptr,ibit
         print*,(nodlst(ival,iptr),ival=1,3)
         print*,(int(nodlst(ival,iptr)),ival=4,5)
         RETURN
      END IF
C
      IFLAG = INT(NODLST(6,IPTR))
C
      IF( IFLAG .EQ. 0 ) THEN
         IF( IBIT .EQ. 1 ) THEN
            GTNDFL = .FALSE.
         ELSE
            IERR = -2
            PRINT*,'*** GTNDFL: NODE BLOCK IS EMPTY'
            print*,iptr,ibit
            print*,(nodlst(ival,iptr),ival=1,3)
            print*,(int(nodlst(ival,iptr)),ival=4,5)
         END IF
         RETURN
      ELSE IF( IFLAG .LT. 0 ) THEN
         IERR = -22
         PRINT*,'*** GTNDFL: NODE BLOCK IS INVALID'
         print*,iptr,ibit,iflag
         print*,(nodlst(ival,iptr),ival=1,3)
         print*,(int(nodlst(ival,iptr)),ival=4,5)
         RETURN
      END IF
C
c      IVAL = 2**(IBIT-1)
      ival = powers(ibit)
C
      GTNDFL = ( IAND(IFLAG,IVAL) .NE. 0 )
C
      RETURN
      END
C
#endif
