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
      SUBROUTINE ADDELE( BIGLST, NODLST, CONELM, EDGLST, IPTR,
     :                   INEWEL, IREG )
C-----------------------------------------------------------------------
C
C - This subroutine adds a new element into biglst, and also adds the new
C - element to connected elements' lists of connected elements.
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      LOGICAL BLDING
      COMMON / BLKBLD / BLDING
C
      INTEGER CONELM(4), EDGLST(6), IPTR, INEWEL, IREG
C
      INTEGER NXT, LST, JFLAG, i, nods(4), nod1, nod2, nod3, nod4,
     :        i1, i2, dumnxt, dumlst, opp1, opp2, opp3,
     :        IFAC, IED1, IED2, IED3
C
c      logical debug
c      debug = .false.
c      if( emtbig .eq. 1677 ) debug = .true.
      CALL ADAPT_SECOND( TIMSTT )
c
c      print*,'Entered ADDELE...  IREG: ',ireg
c      print*,'   CONELM: ',conelm
c      print*,'   EDGLST: ',edglst
C
      IF( IERR .NE. 0 ) RETURN
C
      IPTR = EMTBIG
      IF( IPTR .LE. 0 ) THEN
         IERR = -4
         WRITE(0,*) '*** FAILED TO ADD ELEMENT TO BIGLST - FULL 1'
         RETURN
      END IF
C
      JFLAG = BIGLST(4,IPTR)
C
      IF( JFLAG .NE. 0 ) THEN
         WRITE(0,*) '*** FAILED TO ADD ELEMENT TO BIGLST ',
     :          '- BLOCK NOT EMPTY!'
         IERR = -2
         RETURN
      END IF
C
      IF( BIGLST(LSTBIG,IPTR) .NE. ENDBIG ) THEN
         IERR = -110
         WRITE(0, *)'*** ADDELE: INCONSISTENCY AT END OF BIGLST'
         RETURN
      END IF
C
      IF( ENDBIG .GT. 0 ) THEN
         IF( BIGLST(NXTBIG,ENDBIG) .NE. IPTR ) THEN
            IERR = -110
            WRITE(0,*) '*** ADDELE: INCONSISTENCY NEAR END OF BIGLST'
            RETURN
         END IF
      END IF
C
C - We must have the fourth connected element non-negative so that the
C - element block cannot be mistaken for an edge. A negative value for
C - a connected element means that face lies on a surface.
C
c      if(iptr.eq.1315 .or. iptr.eq.551) then
c         print*,'elm: ',iptr,conelm
c         print*,' edges: ',edglst
c      endif
      IF( CONELM(4) .LT. 0 ) THEN
         IF( CONELM(1) .GE. 0 ) THEN
C - swap faces 1 and 4
            JFLAG = CONELM(4)
            CONELM(4) = CONELM(1)
            CONELM(1) = JFLAG
C - also swap edges 1/3 with 4/5 to keep it consistent
            JFLAG = EDGLST(1)
            EDGLST(1) = EDGLST(4)
            EDGLST(4) = JFLAG
            JFLAG = EDGLST(3)
            EDGLST(3) = EDGLST(5)
            EDGLST(5) = JFLAG
         ELSE IF( CONELM(2) .GE. 0 ) THEN
C - swap faces 2 and 4
            JFLAG = CONELM(4)
            CONELM(4) = CONELM(2)
            CONELM(2) = JFLAG
C - also swap edges 1/5 with 2/6 to keep it consistent
            JFLAG = EDGLST(1)
            EDGLST(1) = EDGLST(2)
            EDGLST(2) = JFLAG
            JFLAG = EDGLST(5)
            EDGLST(5) = EDGLST(6)
            EDGLST(6) = JFLAG
         ELSE IF( CONELM(3) .GE. 0 ) THEN
C - swap faces 3 and 4
            JFLAG = CONELM(4)
            CONELM(4) = CONELM(3)
            CONELM(3) = JFLAG
C - also swap edges 2/4 with 3/6 to keep it consistent
            JFLAG = EDGLST(2)
            EDGLST(2) = EDGLST(3)
            EDGLST(3) = JFLAG
            JFLAG = EDGLST(4)
            EDGLST(4) = EDGLST(6)
            EDGLST(6) = JFLAG
         ELSE
            WRITE(0,*) '*** ADDELE: CANNOT HAVE ALL SIDES ON SURFACES'
            WRITE(0,*) CONELM
            IERR = -98
            RETURN
         END IF
      END IF
C
c      IF( .NOT. BLDING ) THEN
c        IF( CONELM(4) .EQ. 0 .OR.
c     :    CONELM(3) .EQ. 0 .OR.
c     :    CONELM(2) .EQ. 0 .OR.
c     :    CONELM(1) .EQ. 0 ) THEN
c           PRINT*,'*** ADDELE: FACE MISSING A CONNECTION VALUE'
c           PRINT*,'   elem: ',iptr
c           PRINT*,'   conelm: ',conelm
c           PRINT*,'   edglst: ',edglst
c           IERR = -20202
c           RETURN
c        END IF
c      END IF
C
      BIGLST(1,IPTR) = CONELM(1)
      BIGLST(2,IPTR) = CONELM(2)
      BIGLST(3,IPTR) = CONELM(3)
      BIGLST(4,IPTR) = CONELM(4)*2 + 1
c
      biglst(7,iptr) = 0
C
      BIGLST(LSTBIG,IPTR) = ENDBIG
      ENDBIG = IPTR
      NXT = BIGLST(NXTBIG,IPTR)
      IF( NXT .LE. 0 ) THEN
         IERR = -4
         WRITE(0,*)'*** FAILED TO ADD ELEMENT TO BIGLST - FULL 2'
         RETURN
      END IF
C
      BIGLST(1,NXT) = EDGLST(1)
      BIGLST(2,NXT) = EDGLST(2)
      BIGLST(3,NXT) = EDGLST(3)
      BIGLST(4,NXT) = EDGLST(4)
c
      biglst(7,nxt) = 0
C
      BIGLST(LSTBIG,NXT) = ENDBIG
      ENDBIG = NXT
      NXT = BIGLST(NXTBIG,NXT)
      IF( NXT .LE. 0 ) THEN
         IERR = -4
         WRITE(0,*)'*** FAILED TO ADD ELEMENT TO BIGLST - FULL 3'
         RETURN
      END IF
C
      BIGLST(1,NXT) = EDGLST(5)
      BIGLST(2,NXT) = EDGLST(6)
C
      CALL SETFLG( BIGLST, NODLST, EDGLST(1), 4, .FALSE. )
      CALL SETFLG( BIGLST, NODLST, EDGLST(2), 4, .FALSE. )
      CALL SETFLG( BIGLST, NODLST, EDGLST(3), 4, .FALSE. )
      CALL SETFLG( BIGLST, NODLST, EDGLST(4), 4, .FALSE. )
      CALL SETFLG( BIGLST, NODLST, EDGLST(5), 4, .FALSE. )
      CALL SETFLG( BIGLST, NODLST, EDGLST(6), 4, .FALSE. )
c
      NOD1 = BIGLST( 1, EDGLST(1) )
      NOD2 = BIGLST( 2, EDGLST(1) )
      nod3 = biglst( 1, edglst(5) )
      nod4 = biglst( 2, edglst(5) )
C
c      I = 1
cC
c 100  IF( I .LT. 6 ) THEN
cC
c         I = I + 1
cC
c         NOD3 = BIGLST( 1, EDGLST(I) )
c         NOD4 = BIGLST( 2, EDGLST(I) )
cC
c         IF( NOD3 .EQ. NOD1 ) GOTO 100
c         IF( NOD3 .EQ. NOD2 ) GOTO 100
c         IF( NOD4 .EQ. NOD1 ) GOTO 100
c         IF( NOD4 .EQ. NOD2 ) GOTO 100
cC
c      ELSE
cC
c         PRINT*,'*** ADDELE: COULD NOT FIND OPPOSITE EDGE'
c         PRINT*,'  newel: ',iptr
c         print*,'    con: ',conelm
c         do i = 1, 6
c            print*,edglst(i),biglst(1,edglst(i)),biglst(2,edglst(i))
c         end do
c         IERR = -150
c         RETURN
cC
c      END IF
cC
c      opp1 = i
c      
      do i = 1, 6
c
         i1 = biglst(1,edglst(i))
         i2 = biglst(2,edglst(i))
c
         if( i1 .ne. nod1 .and.
     :       i1 .ne. nod2 .and.
     :       i1 .ne. nod3 .and.
     :       i1 .ne. nod4 ) then
            write(0,*)'*** ADDELE: EDGE HAS NODE NOT IN ELEMENT!'
            write(0,*)'   newel: ',iptr,ireg,i
            write(0,*)'   nodes: ',nod1,nod2,nod3,nod4
            write(0,*)'   conel: ',conelm
            write(0,*)'   edges: ',edglst
            write(0,*)'   ednds: ',i1,i2,edglst(i)
            ierr = -11112
            return
         end if
c
         if( i2 .ne. nod1 .and.
     :       i2 .ne. nod2 .and.
     :       i2 .ne. nod3 .and.
     :       i2 .ne. nod4 ) then
            write(0,*)'*** ADDELE: EDGE HAS NODE NOT IN ELEMENT!'
            write(0,*)'   newel: ',iptr,ireg,i
            write(0,*)'   nodes: ',nod1,nod2,nod3,nod4
            write(0,*)'   conel: ',conelm
            write(0,*)'   edges: ',edglst
            write(0,*)'   ednds: ',i1,i2,edglst(i)
            ierr = -11112
            return
         end if
c
      end do
C
      CALL STNDFL( BIGLST, NODLST, NOD1, 6, .FALSE. )
      CALL STNDFL( BIGLST, NODLST, NOD2, 6, .FALSE. )
      CALL STNDFL( BIGLST, NODLST, NOD3, 6, .FALSE. )
      CALL STNDFL( BIGLST, NODLST, NOD4, 6, .FALSE. )
C
c      IF( INEWEL .GE. 0 ) THEN
c         NNEWEL = NNEWEL + 1
c         INEWEL = NNEWEL
c         NEWELE(INEWEL) = IPTR
c      ELSE
         INEWEL = 0
c      END IF
C
      BIGLST(3,NXT) = INEWEL
      BIGLST(4,NXT) = 1
      BIGLST(7,NXT) = IREG
C
c      BIGLST(LSTBIG,NXT) = ENDBIG
c      IF( ENDBIG .GT. 0 ) BIGLST(NXTBIG,ENDBIG) = IPTR
      ENDBIG = NXT
      EMTBIG = BIGLST(NXTBIG,NXT)
C
      NELEMS = NELEMS + 1
C
C - Need to tell connected elements this new element is connected to them.
C - Since rmvele should have been called already (to remove elements this
C - one is overwriting - or else the mesh is just being formed so no
C - element exists there anyway) there should be zeroes in the connected
C - elements' lists of connected elements ready to be filled by this new
C - element.
C
      IF( CONELM(1) .GT. 0 ) THEN
         IFAC = 1
         IED1 = 0
         CALL NWCNEL( BIGLST, NODLST, IPTR, CONELM(1),
     :                IFAC, IED1, IED2, IED3 )
C         CALL CNCTEL( BIGLST, NODLST, IPTR, CONELM(1) )
         CALL SETFLG( BIGLST, NODLST, CONELM(1), 3, .FALSE. )
      END IF
      IF( CONELM(2) .GT. 0 ) THEN
         IFAC = 2
         IED1 = 0
         CALL NWCNEL( BIGLST, NODLST, IPTR, CONELM(2),
     :                IFAC, IED1, IED2, IED3 )
C         CALL CNCTEL( BIGLST, NODLST, IPTR, CONELM(2) )
         CALL SETFLG( BIGLST, NODLST, CONELM(2), 3, .FALSE. )
      END IF
      IF( CONELM(3) .GT. 0 ) THEN
         IFAC = 3
         IED1 = 0
         CALL NWCNEL( BIGLST, NODLST, IPTR, CONELM(3),
     :                IFAC, IED1, IED2, IED3 )
C         CALL CNCTEL( BIGLST, NODLST, IPTR, CONELM(3) )
         CALL SETFLG( BIGLST, NODLST, CONELM(3), 3, .FALSE. )
      END IF
      IF( CONELM(4) .GT. 0 ) THEN
         IFAC = 4
         IED1 = 0
         CALL NWCNEL( BIGLST, NODLST, IPTR, CONELM(4),
     :                IFAC, IED1, IED2, IED3 )
C         CALL CNCTEL( BIGLST, NODLST, IPTR, CONELM(4) )
         CALL SETFLG( BIGLST, NODLST, CONELM(4), 3, .FALSE. )
      END IF
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,*)'*** ADDELE: Got error from NWCNEL: ',IERR
         RETURN
      END IF
C
      IF( ADPTNG ) THEN
C
         ELSLFT = ELSLFT + 1
C
         NELADD = NELADD + 1
C
         IF( ABS(ELSLFT-LSTLFT) .EQ. 4000 ) THEN
            LSTLFT = ELSLFT
            ewrite(1,3)
     :           'LEFT/TOTAL:',ELSLFT,NELEMS,
     :           '      ADD/SUB/CHK:',NELADD,NELSUB,NELCHK,
     :           '      NODES:',NUMNDS
  3         format( a, 2i8, a, 3i8, a, i8 )
            NELADD = 0
            NELSUB = 0
            NELCHK = 0
            mvcnt  = 0
         END IF
C
      END IF
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMADD = TIMADD + TIMEND - TIMSTT
      TOTADD = TOTADD + 1
c
c      if( iptr .eq. 5845656 ) then
c        WRITE(0,*)' '
c        ewrite(3,*)'>>>>>>>>>>>> ADDING ELEMENT AT ',iptr,' <<<<<<<<<<<<'
c        WRITE(0,*)'   CONELM: ',conelm
c        WRITE(0,*)'   EDGLST: ',edglst
c        do i = 1, 6
c          ewrite(3,*)'   edg: ',edglst(i),
c     :          biglst(1,edglst(i)),biglst(2,edglst(i))
c        end do
c        call elnods( biglst, nodlst, iptr, nods, nxt, lst, .true. )
c        ewrite(3,*)'   NODS: ',nods
c        ewrite(3,*)'   node 1: ',(nodlst(i,nods(1)),i=1,3)
c        ewrite(3,*)'   node 2: ',(nodlst(i,nods(2)),i=1,3)
c        ewrite(3,*)'   node 3: ',(nodlst(i,nods(3)),i=1,3)
c        ewrite(3,*)'   node 4: ',(nodlst(i,nods(4)),i=1,3)
c        WRITE(0,*)'   INEWEL: ',INEWEL,nxt,lst,ireg
c        if( ierr .ne. 0 ) stop
cc        pause
c      end if
C
c      call OPPEDG(BIGLST,NODLST,edglst(2),iptr,opp2,dumnxt,dumlst)
c      call OPPEDG(BIGLST,NODLST,edglst(3),iptr,opp3,dumnxt,dumlst)
cc
c      if(opp1.ne.5) then
c         ewrite(3,*)"+++ ADDELE INFO: inconsistent oppedg ",1,opp1
c         ierr = 12345
c      end if
cc
c      if(opp2.ne.edglst(6)) then
c         ewrite(3,*)"+++ ADDELE INFO: inconsistent oppedg ",2,opp2,edglst
c         ierr = 12345
c      end if
cc
c      if(opp3.ne.edglst(4)) then
c         ewrite(3,*)"+++ ADDELE INFO: inconsistent oppedg ",3,opp3,edglst
c         ierr = 12345
c      end if
cc
      RETURN
      END SUBROUTINE
C
