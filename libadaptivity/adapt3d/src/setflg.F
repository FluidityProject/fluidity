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
      SUBROUTINE SETFLG( BIGLST, NODLST, IPTR, IBIT, VALUE )
C-----------------------------------------------------------------------
C
C - This subroutine sets/unsets bit IBIT (from 2 to 29) of the flag for the
C - element or edge in block IPTR to VALUE (.TRUE. or .FALSE.)
C
C - Bits 1 to 10 are reserved for internal use...
C
C - Edges make use of bits 1 to 10 as follows:
C - Bit 1: set if block is not empty - not allowed using this subroutine
C - Bit 2: used to show geometry edges
C - Bit 3: used to show internal edges (should be overidden by bit 2)
C - Bit 4: used temporarily (by CHKFAC) for finding edges through faces
C - Bit 4: used temporarily (by ADAPT)  for finding edges through elements
C - Bit 5: used temporarily (by NDCNCT) for finding edges round a node
C - Bit 6: used to flag edges not to be touched
C - Bit 7: used (by ADAPT1) to keep a tally of elements added/removed/checked
C - Bit 8: used temporarily (by SETGMY) for calculating geometry information
C
C - Elements make use of bits 1 to 10 as follows:
C - Bit 1: set if block is not empty (cannot be set using SETFLG)
C - Bit 2: used temporarily (by EDGELS) for finding elements round an edge
C - Bit 3: used temporarily (by ADAPT1)  for finding pairs of elements
C - Bit 4: used temporarily (by NDCNCT) for finding elements round a node
C - Bit 5: used to flag elements not to be touched
C - Bit 6: used temporarily (by EDGTST) to flag surrounding elements during edge collapse
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
      LOGICAL VALUE
C
      call set_flag(biglst, iptr, ibit, value)
C
      end subroutine setflg
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
      INTEGER IFLAG, IVAL, NXT, LST
C
      LOGICAL GETFLG, ISEDGE
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IBIT .LE. 1 .OR. IBIT .GE. 30 ) THEN
         IERR = -10
         PRINT*,'*** SETFLG: INVALID BIT NUMBER'
         print*,iptr,ibit,value
         RETURN
      END IF
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** SETFLG: INVALID POINTER'
         print*,iptr,ibit,value
         RETURN
      END IF
C
      IFLAG = BIGLST(4,IPTR)
C
      IF( IFLAG .EQ. 0 ) THEN
         IERR = -2
         PRINT*,'*** SETFLG: BLOCK IS EMPTY'
         print*,iptr,ibit,value
         RETURN
      END IF
C
      IF( IFLAG .LT. 0 ) THEN
C
         IFLAG = -IFLAG
         NXT = IPTR
         ISEDGE = .TRUE.
C
      ELSE
C
         ISEDGE = .FALSE.
         NXT = BIGLST(NXTBIG,IPTR)
         IFLAG = BIGLST(4,NXT)
C
         IF( IFLAG .LE. 0 ) THEN
            IERR = -3
            PRINT*,'*** SETFLG: SECOND BLOCK IS NOT AN ELEMENT'
            RETURN
         END IF
C
         NXT = BIGLST(NXTBIG,NXT)
         IFLAG = BIGLST(4,NXT)
C
         IF( IFLAG .LE. 0 ) THEN
            IERR = -4
            PRINT*,'*** SETFLG: THIRD BLOCK IS NOT AN ELEMENT'
            RETURN
         END IF
C
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
         IF( ISEDGE ) IFLAG = -IFLAG
C
         BIGLST(4,NXT) = IFLAG
C
      END IF
C
      RETURN
      END
C
#endif
