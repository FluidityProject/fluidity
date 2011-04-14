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
      LOGICAL FUNCTION GETFLG( BIGLST, NODLST, IPTR, IBIT )
C-----------------------------------------------------------------------
C
C - This subroutine checks the flag for the element or edge in block IPTR
C - to see if bit IBIT (from 1 to 29) is set (.TRUE.) or unset (.FALSE.).
C
C - Bits 1 to 10 are reserved for internal use...
C
C - Edges make use of bits 1 to 10 as follows:
C - Bit 1: set if block is not empty (cannot be set using SETFLG)
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
      getflg = get_flag(biglst, iptr, ibit)
C
      end function getflg
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
      GETFLG = .FALSE.
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IBIT .LE. 0 .OR. IBIT .GE. 30 ) THEN
         IERR = -10
         PRINT*,'*** GETFLG: INVALID BIT NUMBER'
         print*,'   iptr,bit: ',iptr,ibit
         RETURN
      END IF
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** GETFLG: INVALID POINTER'
         print*,'   iptr,bit: ',iptr,ibit
         RETURN
      END IF
C
      IFLAG = BIGLST(4,IPTR)
C
      IF( IFLAG .EQ. 0 ) THEN
C
         IF( IBIT .EQ. 1 ) THEN
            GETFLG = .FALSE.
         ELSE
            IERR = -2
            PRINT*,'*** GETFLG: BLOCK IS EMPTY'
            print*,'   iptr,bit: ',iptr,ibit
         END IF
         RETURN
C
      ELSE IF( IFLAG .LT. 0 ) THEN
C
         IFLAG = -IFLAG
C
      ELSE
C
         NXT = BIGLST(NXTBIG,IPTR)
         IFLAG = BIGLST(4,NXT)
C
         IF( IFLAG .LE. 0 ) THEN
            IERR = -3
            PRINT*,'*** GETFLG: SECOND BLOCK IS NOT AN ELEMENT'
            print*,'   1st block: ',iptr,ibit
            print*,'   2nd block: ',nxt,iflag
            RETURN
         END IF
C
         NXT = BIGLST(NXTBIG,NXT)
         IFLAG = BIGLST(4,NXT)
C
         IF( IFLAG .LE. 0 ) THEN
            IERR = -4
            PRINT*,'*** GETFLG: THIRD BLOCK IS NOT AN ELEMENT'
            print*,'   1st block: ',iptr,ibit
            print*,'   2nd block: ',biglst(nxtbig,iptr)
            print*,'   3rd block: ',nxt,iflag
            RETURN
         END IF
C
      END IF
C
c      IVAL = 2**(IBIT-1)
      ival = powers(ibit)
C
      GETFLG = ( IAND(IFLAG,IVAL) .NE. 0 )
C
      RETURN
      END
C
#endif
