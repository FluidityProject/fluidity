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
      SUBROUTINE CMNEDS( BIGLST, NODLST, ELM1, ELM2, ED1, ED2, ED3 )
C-----------------------------------------------------------------------
C
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
      INTEGER ELM1, ELM2, ED1, ED2, ED3
C
      INTEGER CONELM(4), EDGLST(6), INEW, IFLAG, NXT, LST, HASEDG, I
C
      ED1 = 0
      ED2 = 0
      ED3 = 0
C
      CALL ELMINF( BIGLST, NODLST, ELM1, CONELM, EDGLST,
     :             INEW, IFLAG, I, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      I = 0
C
 10   IF( I .LT. 6 ) THEN
         I = I + 1
         ED1 = HASEDG( BIGLST, NODLST, ELM2, EDGLST(I) )
         IF( ED1 .EQ. 0 ) GOTO 10
         ED1 = EDGLST(I)
      END IF
C
 20   IF( I .LT. 6 ) THEN
         I = I + 1
         ED2 = HASEDG( BIGLST, NODLST, ELM2, EDGLST(I) )
         IF( ED2 .EQ. 0 ) GOTO 20
         ED2 = EDGLST(I)
      END IF
C
 30   IF( I .LT. 6 ) THEN
         I = I + 1
         ED3 = HASEDG( BIGLST, NODLST, ELM2, EDGLST(I) )
         IF( ED3 .EQ. 0 ) GOTO 30
         ED3 = EDGLST(I)
      END IF
C
      RETURN
      END
C
