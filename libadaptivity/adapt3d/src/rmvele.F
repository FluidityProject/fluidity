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
#include "ewrite.h"
      SUBROUTINE RMVELE( BIGLST, NODLST, IPTR )
C-----------------------------------------------------------------------
C
C - This subroutine removes an element from BIGLST - and also removes
C - that element from connected elements' lists of connected elements
C - (to be filled in by new element(s) which should be created immediately).
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blknxt.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      INTEGER IPTR
C
      INTEGER ELM1, ELM2, ELM3, ELM4, CONELM(4), EDGLST(6), IELM,
     :        INEWEL, NXT, LST, IFLAG, IPTR2, IPTR3
C
      logical getflg
c
      CALL ADAPT_SECOND( TIMSTT )
c      debug = .false.
c      if( iptr .eq. 34117 .or. iptr .eq. 57095 ) debug = .true.
C
C      IF( IERR .NE. 0 ) RETURN
C
C      CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
C     :             INEWEL, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
cc      IF( IPTR .LE. 0 ) THEN
cc         IERR = -1
cc         PRINT*,'*** RMVELE: INVALID POINTER'
cc         RETURN
cc      ELSE IF( BIGLST(4,IPTR) .LT. 0 ) THEN
cc         IERR = -3
cc         PRINT*,'*** RMVELE: FIRST BLOCK IS AN EDGE BLOCK'
cc         print*,'    element: ',iptr,biglst(4,iptr)
cc         RETURN
cc      ELSE IF( BIGLST(4,IPTR) .EQ. 0 ) THEN
cc         IERR = -2
cc         PRINT*,'*** RMVELE: FIRST BLOCK IS EMPTY'
cc         print*,'    element: ',iptr,biglst(4,iptr)
cc         RETURN
cc      END IF
C
      ELM1 = BIGLST(1,IPTR)
      ELM2 = BIGLST(2,IPTR)
      ELM3 = BIGLST(3,IPTR)
      ELM4 = BIGLST(4,IPTR)/2
C
      LST  = BIGLST(LSTBIG,IPTR)
c
c      if( debug ) then
c         print*,'>>>>>>>>> removing element ',iptr,' <<<<<<<<<'
c      print*,nextls
c         print*,'   conelm: ',conelm
c      end if
C
      IPTR2 = BIGLST(NXTBIG,IPTR)
C
cc      IF( IPTR2 .LE. 0 ) THEN
cc         IERR = -1
cc         PRINT*,'*** RMVELE: INVALID POINTER IN FIRST BLOCK'
cc         RETURN
cc      ELSE IF( BIGLST(4,IPTR2) .LT. 0 ) THEN
cc         IERR = -3
cc         PRINT*,'*** RMVELE: SECOND BLOCK IS AN EDGE BLOCK'
cc         print*,'    element: ',iptr,biglst(4,iptr)
cc         RETURN
cc      ELSE IF( BIGLST(4,IPTR2) .EQ. 0 ) THEN
cc         IERR = -2
cc         PRINT*,'*** RMVELE: SECOND BLOCK IS EMPTY'
cc         print*,'    element: ',iptr,biglst(4,iptr)
cc         RETURN
cc      END IF
C
      IPTR3 = BIGLST(NXTBIG,IPTR2)
C
cc      IF( IPTR3 .LE. 0 ) THEN
cc         IERR = -1
cc         PRINT*,'*** RMVELE: INVALID POINTER IN SECOND BLOCK'
cc         RETURN
cc      ELSE IF( BIGLST(4,IPTR3) .LT. 0 ) THEN
cc         IERR = -3
cc         PRINT*,'*** RMVELE: THIRD BLOCK IS AN EDGE BLOCK'
cc         print*,'    element: ',iptr,biglst(4,iptr)
cc         RETURN
cc      ELSE IF( BIGLST(4,IPTR3) .EQ. 0 ) THEN
cc         IERR = -2
cc         PRINT*,'*** RMVELE: THIRD BLOCK IS EMPTY'
cc         print*,'    element: ',iptr,biglst(4,iptr)
cc         RETURN
cc      END IF
C
      NXT = BIGLST(NXTBIG,IPTR3)
C
      IF( ADPTNG ) THEN
C
         NELSUB = NELSUB + 1
C
         IF( .NOT. GETFLG(BIGLST,NODLST,IPTR,7) ) THEN
C
            ELSLFT = ELSLFT - 1
C
            IF( ABS(ELSLFT-LSTLFT) .EQ. 4000 ) THEN
               LSTLFT = ELSLFT
               ewrite(1,3) 
     :              'LEFT/TOTAL:',ELSLFT,NELEMS,
     :              '      ADD/SUB/CHK:',NELADD,NELSUB,NELCHK,
     :              '      NODES:',NUMNDS
  3            format( a, 2i9, a, 3i8, a, i9 )
               NELADD = 0
               NELSUB = 0
               NELCHK = 0
               mvcnt  = 0
            END IF
C
         END IF
C
      END IF
C
C      IF( ELM1 .LT. 0 ) THEN
C         IERR = -2
C         PRINT*,'*** RMVELE: FOUND INVALID CONNECTED ELEMENT 1'
C         RETURN
C      END IF
C
      IF( ELM1 .GT. 0 ) THEN
         IF( BIGLST(1,ELM1) .EQ. IPTR ) THEN
            BIGLST(1,ELM1) = 0
         ELSE IF( BIGLST(2,ELM1) .EQ. IPTR ) THEN
            BIGLST(2,ELM1) = 0
         ELSE IF( BIGLST(3,ELM1) .EQ. IPTR ) THEN
            BIGLST(3,ELM1) = 0
         ELSE IF( BIGLST(4,ELM1)/2 .EQ. IPTR ) THEN
            BIGLST(4,ELM1) = 1
         END IF
      END IF
C
C      IF( ELM2 .LT. 0 ) THEN
C         IERR = -2
C         PRINT*,'*** RMVELE: FOUND INVALID CONNECTED ELEMENT 2'
C         RETURN
C      END IF
C
      IF( ELM2 .GT. 0 ) THEN
         IF( BIGLST(1,ELM2) .EQ. IPTR ) THEN
            BIGLST(1,ELM2) = 0
         ELSE IF( BIGLST(2,ELM2) .EQ. IPTR ) THEN
            BIGLST(2,ELM2) = 0
         ELSE IF( BIGLST(3,ELM2) .EQ. IPTR ) THEN
            BIGLST(3,ELM2) = 0
         ELSE IF( BIGLST(4,ELM2)/2 .EQ. IPTR ) THEN
            BIGLST(4,ELM2) = 1
         END IF
      END IF
C
C      IF( ELM3 .LT. 0 ) THEN
C         IERR = -2
C         PRINT*,'*** RMVELE: FOUND INVALID CONNECTED ELEMENT 3'
C         RETURN
C      END IF
C
      IF( ELM3 .GT. 0 ) THEN
         IF( BIGLST(1,ELM3) .EQ. IPTR ) THEN
            BIGLST(1,ELM3) = 0
         ELSE IF( BIGLST(2,ELM3) .EQ. IPTR ) THEN
            BIGLST(2,ELM3) = 0
         ELSE IF( BIGLST(3,ELM3) .EQ. IPTR ) THEN
            BIGLST(3,ELM3) = 0
         ELSE IF( BIGLST(4,ELM3)/2 .EQ. IPTR ) THEN
            BIGLST(4,ELM3) = 1
         END IF
      END IF
C
C      IF( ELM4 .LT. 0 ) THEN
C         IERR = -2
C         PRINT*,'*** RMVELE: FOUND INVALID CONNECTED ELEMENT 4'
C         RETURN
C      END IF
C
      IF( ELM4 .GT. 0 ) THEN
         IF( BIGLST(1,ELM4) .EQ. IPTR ) THEN
            BIGLST(1,ELM4) = 0
         ELSE IF( BIGLST(2,ELM4) .EQ. IPTR ) THEN
            BIGLST(2,ELM4) = 0
         ELSE IF( BIGLST(3,ELM4) .EQ. IPTR ) THEN
            BIGLST(3,ELM4) = 0
         ELSE IF( BIGLST(4,ELM4)/2 .EQ. IPTR ) THEN
            BIGLST(4,ELM4) = 1
         END IF
      END IF
C
      BIGLST( 1, IPTR )  = 0
      BIGLST( 2, IPTR )  = 0
      BIGLST( 3, IPTR )  = 0
      BIGLST( 4, IPTR )  = 0
c
      biglst( 7, iptr )  = 0
C
      BIGLST( 1, IPTR2 ) = 0
      BIGLST( 2, IPTR2 ) = 0
      BIGLST( 3, IPTR2 ) = 0
      BIGLST( 4, IPTR2 ) = 0
c
      biglst( 7, iptr2 ) = 0
C
      BIGLST( 1, IPTR3 ) = 0
      BIGLST( 2, IPTR3 ) = 0
      BIGLST( 3, IPTR3 ) = 0
      BIGLST( 4, IPTR3 ) = 0
c
      biglst( 7, iptr3 ) = 0
C
c      IF( INEWEL .GT. 0 ) THEN
c         NEWELE(INEWEL) = NEWELE(NNEWEL)
c         NNEWEL = NNEWEL - 1
c      END IF
C
      DO IFLAG = 1, MXNEXT
         IF( NEXTLS(IFLAG) .EQ. IPTR .AND. ENDBIG .NE. IPTR3 ) THEN
c            print*,'>>>> RMVELE: REMOVING NEXT ELEMENT: ',
c     :              IPTR,NXT,NELEMS,NEDGES
            NEXTLS(IFLAG) = NXT
         ELSE IF( NEXTLS(IFLAG) .EQ. EMTBIG ) THEN
c            print*,'>>>> RMVELE: REMOVING LAST ELEMENT: ',
c     :              IPTR,EMTBIG,NELEMS,NEDGES
            NEXTLS(IFLAG) = IPTR
         END IF
      END DO
C
      IF( IPTR .EQ. STTBIG ) STTBIG = NXT
C
      IF( LST .GT. 0 ) BIGLST( NXTBIG, LST ) = NXT
      IF( NXT .GT. 0 ) BIGLST( LSTBIG, NXT ) = LST
C
      IF( ENDBIG .EQ. IPTR3 ) ENDBIG = LST
C
      BIGLST( LSTBIG, IPTR   ) = ENDBIG
      BIGLST( NXTBIG, IPTR3  ) = EMTBIG
C
      IF( EMTBIG .GT. 0 ) BIGLST( LSTBIG, EMTBIG ) = IPTR3
      IF( ENDBIG .GT. 0 ) BIGLST( NXTBIG, ENDBIG ) = IPTR
C
      EMTBIG = IPTR
C
      NELEMS = NELEMS - 1
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMSUB = TIMSUB + TIMEND - TIMSTT
      TOTSUB = TOTSUB + 1
c
c      print*,'>>>>>>>>> removed element ',iptr,' <<<<<<<<<'
c      print*,nextls
C
      RETURN
      END
C
