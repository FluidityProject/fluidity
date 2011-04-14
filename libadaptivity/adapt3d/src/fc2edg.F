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
      SUBROUTINE FC2EDG( BIGLST, NODLST, IELM1, IELM2 )
C-----------------------------------------------------------------------
C
C
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
      INTEGER IELM1, IELM2
C
      INTEGER CONEL1(4), CONEL2(4), EDGLS1(6), EDGLS2(6), INEW, IFLAG,
     :        NXT, LST, NODS1(4), NODS2(4), I, J, CONELS(6), N1, N2,
     :        N3, NE1, NE2, CONELM(4), EDGLST(6), IR1, IR2, INFO(4),
     :        FAC1, FAC2, ED1, ED2, ED3, OPED1, OPED2, OPED3, OPED4,
     :        OPED5, OPED6, E1, E2, E3, NODS(4), ED4, ED5, ED6, ED7,
     :        ED8, ED9, FAC3, FAC4, FAC5, OPPLST(6), IEDG
C
      REAL DX, DY, DZ, DOT1, DOT2, DOT3, elmfnc
C
      LOGICAL goterr, addit, gtndfl, getflg
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,*) '*** FC2EDG: GOT ERROR ON ENTRY'
         RETURN
      END IF
C
      CALL ADAPT_SECOND( TIMSTT )
C
      goterr = debug
c      if( ielm1 .eq. 38431 .and. ielm2 .eq. 34490 ) goterr = .true.
c      if( ielm1 .eq. 34490 .and. ielm2 .eq. 44665 ) stop
c
      if( goterr ) then
        WRITE(0,*) ' '
        WRITE(0,*) 'Into FC2EDG... ',ierr
        WRITE(0,*) '   IELM1,IELM2: ',IELM1,IELM2
      end if
C
      CALL ELMINF( BIGLST, NODLST, IELM1, CONEL1, EDGLS1,
     :             INEW, IFLAG, IR1, NXT, LST )
      CALL ELNODS( BIGLST, NODLST, IELM1, NODS1, NXT, LST, .TRUE. )
      CALL ELMINF( BIGLST, NODLST, IELM2, CONEL2, EDGLS2,
     :             INEW, IFLAG, IR2, NXT, LST )
      CALL ELNODS( BIGLST, NODLST, IELM2, NODS2, NXT, LST, .TRUE. )
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,*) '*** FC2EDG: GOT ERROR FROM ELMINF/ELNODS'
         RETURN
      END IF
C
      IF( IR1 .NE. IR2 ) THEN
         ewrite(0, *)
     :        '+++ Warning: attempted to remove separator face'
         RETURN
      END IF
C
      INFO(1) = 340382
      INFO(2) = 418254
      INFO(3) = 975565
      INFO(4) = 678342
C
      OPPLST(1) = 5
      OPPLST(2) = 6
      OPPLST(3) = 4
      OPPLST(4) = 3
      OPPLST(5) = 1
      OPPLST(6) = 2
C
      ED1 = 0
      ED2 = 0
      ED3 = 0
C
      FAC1 = 0
 10   IF( FAC1 .LT. 4 ) THEN
         FAC1 = FAC1 + 1
         IF( CONEL1(FAC1) .NE. IELM2 ) GOTO 10
      END IF
C
      IF( FAC1 .GT. 4 ) THEN
         WRITE(0,*) '*** FC2EDG: FAILED TO FIND CONNECTED ELEMENT'
         WRITE(0,*) ielm2,ielm2
         WRITE(0,*) conel1
         WRITE(0,*) conel2
         WRITE(0,*) edgls1
         WRITE(0,*) edgls2
         IERR = -160
         RETURN
      END IF
C
      FAC2 = 0
 20   IF( FAC2 .LT. 4 ) THEN
         FAC2 = FAC2 + 1
         IF( CONEL2(FAC2) .NE. IELM1 ) GOTO 20
      END IF
C
      IF( FAC2 .GT. 4 ) THEN
         WRITE(0,*) '*** FC2EDG: FAILED TO FIND CONNECTED ELEMENT'
         WRITE(0,*) ielm2,ielm2
         WRITE(0,*) conel1
         WRITE(0,*) conel2
         WRITE(0,*) edgls1
         WRITE(0,*) edgls2
         IERR = -160
         RETURN
      END IF
C
      NE1 = NODS1(5 - FAC1)
      NE2 = NODS2(5 - FAC2)
      I = INFO(FAC1)
      CONELS(1) = CONEL1(MOD(I,4) + 1)
      I = I / 4
      CONELS(2) = CONEL1(MOD(I,4) + 1)
      I = I / 4
      CONELS(3) = CONEL1(MOD(I,4) + 1)
      I = I / 4
      N1 = NODS1(MOD(I,4) + 1)
      I = I / 4
      N2 = NODS1(MOD(I,4) + 1)
      I = I / 4
      N3 = NODS1(MOD(I,4) + 1)
      I = I / 4
      J = MOD(I,8)
      ED1 = EDGLS1(J)
      OPED1 = EDGLS1(OPPLST(J))
      I = I / 8
      J = MOD(I,8)
      ED2 = EDGLS1(J)
      OPED2 = EDGLS1(OPPLST(J))
      I = I / 8
      J = MOD(I,8)
      ED3 = EDGLS1(J)
      OPED3 = EDGLS1(OPPLST(J))
C
      I = INFO(FAC2)
      FAC3 = CONEL2(MOD(I,4) + 1)
      I = I / 4
      FAC4 = CONEL2(MOD(I,4) + 1)
      I = I / 4
      FAC5 = CONEL2(MOD(I,4) + 1)
      I = I / 256
      J = MOD(I,8)
      ED4 = EDGLS2(J)
      ED7 = EDGLS2(OPPLST(J))
      I = I / 8
      J = MOD(I,8)
      ED5 = EDGLS2(J)
      ED8 = EDGLS2(OPPLST(J))
      I = I / 8
      J = MOD(I,8)
      ED6 = EDGLS2(J)
      ED9 = EDGLS2(OPPLST(J))
C
      CONELS(4) = 0
      IF( ED4 .EQ. ED1 ) THEN
         CONELS(4) = FAC3
         OPED4 = ED7
      ELSE IF( ED5 .EQ. ED1 ) THEN
         CONELS(4) = FAC4
         OPED4 = ED8
      ELSE IF( ED6 .EQ. ED1 ) THEN
         CONELS(4) = FAC5
         OPED4 = ED9
      ELSE
         WRITE(0,*) '*** FC2EDG: FAILED TO FIND COMMON EDGE 1'
         IERR = -161
      END IF
C
      IF( CONELS(1) .LT. 0 .AND. CONELS(4) .LT. 0 ) THEN
         IF( .NOT. GETFLG( BIGLST, NODLST, ED1, 2 ) ) THEN
            ewrite(0,*)
     :         '+++ Aborted face to edge: no geom edge between surfaces'
            RETURN
         END IF
      END IF
C
      CONELS(5) = 0
      IF( ED4 .EQ. ED2 ) THEN
         CONELS(5) = FAC3
         OPED5 = ED7
      ELSE IF( ED5 .EQ. ED2 ) THEN
         CONELS(5) = FAC4
         OPED5 = ED8
      ELSE IF( ED6 .EQ. ED2 ) THEN
         CONELS(5) = FAC5
         OPED5 = ED9
      ELSE
         WRITE(0,*) '*** FC2EDG: FAILED TO FIND COMMON EDGE 2'
         IERR = -161
      END IF
C
      IF( CONELS(2) .LT. 0 .AND. CONELS(5) .LT. 0 ) THEN
         IF( .NOT. GETFLG( BIGLST, NODLST, ED2, 2 ) ) THEN
            ewrite(0,*)
     :         '+++ Aborted face to edge: no geom edge between surfaces'
            RETURN
         END IF
      END IF
C
      CONELS(6) = 0
      IF( ED4 .EQ. ED3 ) THEN
         CONELS(6) = FAC3
         OPED6 = ED7
      ELSE IF( ED5 .EQ. ED3 ) THEN
         CONELS(6) = FAC4
         OPED6 = ED8
      ELSE IF( ED6 .EQ. ED3 ) THEN
         CONELS(6) = FAC5
         OPED6 = ED9
      ELSE
         WRITE(0,*) '*** FC2EDG: FAILED TO FIND COMMON EDGE 3'
         IERR = -161
      END IF
C
      IF( CONELS(3) .LT. 0 .AND. CONELS(6) .LT. 0 ) THEN
         IF( .NOT. GETFLG( BIGLST, NODLST, ED3, 2 ) ) THEN
            ewrite(0,*)
     :         '+++ Aborted face to edge: no geom edge between surfaces'
            RETURN
         END IF
      END IF
C    
      IF( IERR .NE. 0 ) THEN
         WRITE(0,*) ielm1,ielm2
         WRITE(0,*) fac1,fac2
         WRITE(0,*) ed1,ed2,ed3
         WRITE(0,*) oped1,oped2,oped3
         WRITE(0,*) oped4,oped5,oped6
         WRITE(0,*) ed4,ed5,ed6
         WRITE(0,*) ed7,ed8,ed9
         WRITE(0,*) conels
         WRITE(0,*) conel1
         WRITE(0,*) conel2
         WRITE(0,*) edgls1
         WRITE(0,*) edgls2
         return
      END IF
C      
      NODS(1) = NE1
      NODS(2) = NE2
C
      INEW = -1
      CALL ADDEDG( BIGLST, NODLST, NODS, IEDG, INEW )
      if( goterr ) WRITE(0,*) '   Added edge ',iedg,nods(1),nods(2)
C
      CALL SETFLG( BIGLST, NODLST, IEDG, 3, .TRUE. )
C
      CALL RMVELE( BIGLST, NODLST, IELM1 )
      CALL RMVELE( BIGLST, NODLST, IELM2 )
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,*) '*** FC2EDG: GOT ERROR FROM RMVELE'
         RETURN
      END IF
C
      CONELM(1) = CONELS(1)
      CONELM(2) = CONELS(4)
      CONELM(3) = 0
      CONELM(4) = 0
C
      EDGLST(1) = ED1
      EDGLST(2) = OPED2
      EDGLST(3) = OPED3
      EDGLST(4) = OPED5
      EDGLST(6) = OPED6
      EDGLST(5) = IEDG
C
      INEW = -1
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, E1, INEW, IR1 )
      if( goterr ) write(0,*) '   added elm ',e1
      if( ierr.eq.12345 ) then
         ewrite(0,*)
     :        '^^^ inconsistent oppedg from addele in fc2edg (1) ^^^'
         ierr = 0
      else if( ierr .lt. 0 ) then
         ewrite(-1,*) 
     :        '*** FC2EDG: Got error from ADDELE (1): ',IERR
         return
      end if
C
      CONELM(1) = CONELS(2)
      CONELM(2) = CONELS(5)
      CONELM(4) = E1
      CONELM(3) = 0
C
      EDGLST(1) = ED2
      EDGLST(2) = OPED3
      EDGLST(3) = OPED1
      EDGLST(4) = OPED6
      EDGLST(6) = OPED4
      EDGLST(5) = IEDG
C
      INEW = -1
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, E2, INEW, IR1 )
      if( goterr ) write(0,*) '   added elm ',e2
      if( ierr.eq.12345 ) then
         ewrite(0,*)
     :        '^^^ inconsistent oppedg from addele in fc2edg (2) ^^^'
         ierr = 0
      else if( ierr .lt. 0 ) then
         ewrite(-1,*) 
     :        '*** FC2EDG: Got error from ADDELE (2): ',IERR
         return
      end if
C
      CONELM(1) = CONELS(3)
      CONELM(2) = CONELS(6)
      CONELM(4) = E2
      CONELM(3) = E1
C
      EDGLST(1) = ED3
      EDGLST(2) = OPED1
      EDGLST(3) = OPED2
      EDGLST(4) = OPED4
      EDGLST(6) = OPED5
      EDGLST(5) = IEDG
C
      INEW = -1
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, E3, INEW, IR1 )
      if( goterr ) write(0,*) '   added elm ',e3
      if( ierr.eq.12345 ) then
         ewrite(0,*)
     :     '^^^ inconsistent oppedg from addele in fc2edg (3) ^^^'
         ierr = 0
      else if( ierr .lt. 0 ) then
         ewrite(-1,*) 
     :        '*** FC2EDG: Got error from ADDELE (3): ',IERR
         return
      end if
C
c      ewrite(3,*) 'Leaving FC2EDG... ',ierr
c      ewrite(3,*) '   E1,E2,E3: ',E1,E2,E3
C
  950 CONTINUE
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMCHG = TIMCHG + TIMEND - TIMSTT
      TOTCHG = TOTCHG + 1
c
      return
C
C - this is a final temporary consistency check
C
      I = STTBIG
      goterr = .false.
C
 999  addit = .false.
c
      IF( BIGLST(4,I) .LT. 0 ) THEN
         NXT = BIGLST(NXTBIG,I)
      ELSE
         CALL ELMINF( BIGLST, NODLST, I, CONEL1, EDGLS1,
     :                INEW, IFLAG, IR1, NXT, LST )
         DO J = 1, 4
            IF( CONEL1(J) .EQ. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** FC2EDG: INCONSISTENT CONNECTED ELEMENTS'
                  print*,'   element: ',i
                  print*,'   conelm: ',conel1
                  print*,'   edglst: ',edgls1
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,conel1(j)) .le. 0 ) then
               if( .not. addit ) then
                  print*,'*** FC2EDG: INVALID CONNECTED ELEMENT'
                  print*,'   element: ',i
                  print*,'   conelm: ',conel1
                  print*,'   edglst: ',edgls1
                  addit  = .true.
                  goterr = .true.
               end if
               PRINT*,'   invalid element: ',
     :                      conel1(j),biglst(4,conel1(j))
            end if
         end do
         call elnods( BIGLST,NODLST, i, nods, nxt, lst, .false. )
         dx = elmfnc( BIGLST,NODLST, 0,
     :                nods(1), nods(2), nods(3), nods(4), dy )
         if( abs(dy) .lt. 1e-8 ) then
            if( .not. addit ) then
               PRINT*,'*** FC2EDG: SMALL ELEMENT'
               print*,'   element: ',i
               print*,'   conelm: ',conelm
               print*,'   edglst: ',edglst
               addit  = .true.
               goterr = .true.
            end if
            PRINT*,'   vol,fnc: ',dy,dx
         end if
         do j = 1, 6
            nods(1) = biglst(1,edgls1(j))
            nods(2) = biglst(2,edgls1(j))
            IF( edgls1(J) .LE. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** FC2EDG: INCONSISTENT EDGE LIST'
                  print*,'   element: ',i
                  print*,'   conelm: ',conel1
                  print*,'   edglst: ',edgls1
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,edgls1(j)) .ge. 0 ) then
               if( .not. addit ) then
                  PRINT*,'*** FC2EDG: INVALID EDGE'
                  print*,'   element: ',i
                  print*,'   conelm: ',conel1
                  print*,'   edglst: ',edgls1
                  addit  = .true.
                  goterr = .true.
               end if
               PRINT*,'   invalid edge: ',edgls1(j),biglst(4,edgls1(j))
            else if( getflg(BIGLST,NODLST,edgls1(j),2) ) then
               if( .not. gtndfl(BIGLST,NODLST,nods(1),2) .and.
     :             .not. gtndfl(BIGLST,NODLST,nods(1),4) ) then
                  if( .not. addit ) then
                    PRINT*,'*** FC2EDG: INVALID GEOM-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conel1
                    print*,'   edglst: ',edgls1
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edgls1(j),nods(1),nods(2)
               else if(.not. gtndfl(BIGLST,NODLST,nods(2),2) .and.
     :                 .not. gtndfl(BIGLST,NODLST,nods(2),4) ) then
                  if( .not. addit ) then
                    PRINT*,'*** FC2EDG: INVALID GEOM-EDGE NODE 2 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conel1
                    print*,'   edglst: ',edgls1
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edgls1(j),nods(1),nods(2)
               end if
            else if( .not. getflg(BIGLST,NODLST,edgls1(j),3) ) then
               if( gtndfl(BIGLST,NODLST,nods(1),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** FC2EDG: INVALID SURF-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conel1
                    print*,'   edglst: ',edgls1
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edgls1(j),nods(1),nods(2)
               else if( gtndfl(BIGLST,NODLST,nods(2),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** FC2EDG: INVALID SURF-EDGE NODE 2 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conel1
                    print*,'   edglst: ',edgls1
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edgls1(j),nods(1),nods(2)
               end if
            end if
         end do
      end if
c
      if( goterr ) then
         ierr = -10000
         return
      end if
c
      i = nxt
c
      if( biglst(4,i) .ne. 0 ) goto 999
c
      if( i.ne.emtbig ) then
         print*,'*** FC2EDG: INCONSISTENCY AT END OF BIGLST'
         print*,i,emtbig
         ierr = -10000
      end if
C
      RETURN
      END
C
