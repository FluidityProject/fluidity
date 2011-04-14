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
      SUBROUTINE FACINF( IPTR, CONFC, EDGLST, ELMPTR, INEWFC, IFLAG,
     :                   NXT, LST )
C-----------------------------------------------------------------------
C
C - This subroutine returns the connected face list and the edge list for
C - the face in FRTLST pointed to by IPTR, as well as its position in NEWFAC
C - (if it has one) and pointers to the next and last data in BIGFAC.
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
      INTEGER CONFC(3), EDGLST(3), IPTR, INEWFC, ELMPTR, IFLAG,
     :        NXT, LST
C
c      PRINT*,' '
c      print*,'Into FACINF...',ierr
c      print*,'    IPTR:',IPTR
c      PRINT*,'   CONFC:',CONFC
c      PRINT*,'  EDGLST:',EDGLST
      IF( IERR .NE. 0 ) RETURN
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** FACINF: INVALID POINTER'
         RETURN
      END IF
C
      IFLAG  = FRTLST(8,IPTR) / MAXNWF
C
      IF( IFLAG .LE. 0 ) THEN
         IERR = -2
         PRINT*,'*** FACINF: FACE BLOCK IS EMPTY'
         RETURN
      END IF
C
      CONFC(1) = FRTLST(1,IPTR)
      CONFC(2) = FRTLST(2,IPTR)
      CONFC(3) = FRTLST(3,IPTR)
C
      EDGLST(1) = FRTLST(4,IPTR)
      EDGLST(2) = FRTLST(5,IPTR)
      EDGLST(3) = FRTLST(6,IPTR)
C
      ELMPTR = FRTLST(7,IPTR)
      INEWFC = MOD(FRTLST(8,IPTR),MAXNWF)
C
      LST = FRTLST(LSTFRT,IPTR)
      NXT = FRTLST(NXTFRT,IPTR)
C
c      print*,'Leaving FACINF...'
c      print*,'      IELM:',ELMPTR
c      PRINT*,' INEW,FLAG:',INEWFC,IFLAG
c      PRINT*,'   NXT,LST:',NXT,LST
C
      RETURN
      END
C


