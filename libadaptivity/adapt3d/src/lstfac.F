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
      SUBROUTINE LSTFAC
C-----------------------------------------------------------------------
C
C
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
      INTEGER IPTR, NXT, LST, INEW, I, IFLAG, IELM,
     :        CONFAC(3), EDGLST(3)
C
      LOGICAL GTFRFL
C
      IPTR = STTFRT
C
      print*,'-------------------------------------------------------'
      print*,' '
      PRINT*,'Number of faces: ',NFACES
C
      DO I = 1, NFACES
C
         IF( .NOT. GTFRFL(IPTR,1) ) THEN
            NXT = FRTLST(NXTFRT,IPTR)
            LST = FRTLST(LSTFRT,IPTR)
            PRINT*,'*** PREMATURE EMPTY BLOCK ',IPTR
            PRINT*,' INEW&IFLAG: ',FRTLST(8,IPTR)
            PRINT*,'    NXT,LST: ',NXT,LST
         ELSE
            CALL FACINF( IPTR, CONFAC, EDGLST, IELM,
     :                   INEW, IFLAG, NXT, LST )
            PRINT*,'FACE ',IPTR
            PRINT*,'     CONFAC: ',CONFAC
            PRINT*,'     EDGLST: ',EDGLST
            PRINT*,'       IELM: ',IELM
            PRINT*,'  INEW,FLAG: ',INEW,IFLAG
            PRINT*,'    NXT,LST: ',NXT,LST
         END IF
C
         IPTR = NXT
C
         IF( NXT .LE. 0 ) THEN
            IERR = -1
            PRINT*,'*** FATAL ERROR: NXT =',NXT
            RETURN
         END IF
C
      END DO
C
      RETURN
      END
C
