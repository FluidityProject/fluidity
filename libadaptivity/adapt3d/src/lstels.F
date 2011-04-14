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
      SUBROUTINE LSTELS( BIGLST, NODLST, outfil )
C-----------------------------------------------------------------------
C
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
      integer outfil
C
      INTEGER IPTR, NXT, LST, INEW, I, IFLAG, IREG,
     :        NODS(4), CONELM(4), EDGLST(6), ICNT
C
      LOGICAL GETFLG
C
c      print*,'-------------------------------------------------------'
c      print*,' '
c      PRINT*,'Number of nodes: ',NUMNDS,STTNOD,ENDNOD
C
      IPTR = STTBIG
      ICNT = 0
C
      DO I = 1, NEDGES+NELEMS
C
         IF( .NOT. GETFLG(BIGLST,NODLST,IPTR,1) ) THEN
            NXT = BIGLST(NXTBIG,IPTR)
            LST = BIGLST(LSTBIG,IPTR)
            PRINT*,'*** PREMATURE EMPTY BLOCK ',IPTR
            PRINT*,'      IFLAG: ',BIGLST(4,IPTR)
            PRINT*,'    NXT,LST: ',NXT,LST
         ELSE IF( BIGLST(4,IPTR) .LT. 0 ) THEN
            CALL EDGINF( BIGLST, NODLST, IPTR, NODS,
     :                   INEW, IFLAG, NXT, LST )
C            PRINT*,'EDGE ',IPTR
C            PRINT*,'      NODES: ',NODS
C            PRINT*,'  INEW,FLAG: ',INEW,IFLAG
C            PRINT*,'    NXT,LST: ',NXT,LST
         ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
            CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                   INEW, IFLAG, IREG, NXT, LST )
            CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, .FALSE. )
C            PRINT*,'ELEMENT ',IPTR
C            PRINT*,'     CONELM: ',CONELM
C            PRINT*,'     EDGLST: ',EDGLST
C            PRINT*,'  INEW,FLAG: ',INEW,IFLAG
C            PRINT*,'    NXT,LST: ',NXT,LST
            ICNT = ICNT + 1
            NODS(1) = ABS(INT(NODLST(6,NODS(1))))
            NODS(2) = ABS(INT(NODLST(6,NODS(2))))
            NODS(3) = ABS(INT(NODLST(6,NODS(3))))
            NODS(4) = ABS(INT(NODLST(6,NODS(4))))
            WRITE( outfil, 999 ) ICNT,1,1,IREG,1,NODS
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
 999  FORMAT( 'ELEM ',9I5 )
      RETURN
      END
C
