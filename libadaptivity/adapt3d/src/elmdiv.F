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
      SUBROUTINE ELMDIV( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   IELM )
C-----------------------------------------------------------------------
C
C - This subroutine splits an element into four elements.
C
C *** NOT CURRENTLY USED BY ADAPTIVITY *** - see lines 257-261 adptvy.F
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
      INCLUDE 'blkorg.i'
C
      INTEGER IELM
C
      INTEGER I, J, K, E1, E2, E3, N1, N2, NCMN, P0, P1, P2, NDS(2),
     :        NEWEDS(4), NEWELS(4), EDGLS1(6), CONEL1(4), OPED1, OPED2,
     :        CMNED(12), EDGFLG(6), INEW, INOD, NODS(4), EDGLST(6),
     :        CONELM(4), NXT, LST, IFLAG, INELM
C
      REAL X, Y, Z, d1, d2, dx, dy, elmfnc, edgfnc
C
      LOGICAL GETFLG, goterr, addit, gtndfl
c
      goterr = .false.
      goterr = .true.
C
      CALL ELNODS( BIGLST, NODLST, IELM, NODS, NXT, LST, .FALSE. )
C
      CALL ELMINF( BIGLST, NODLST, IELM, CONELM, EDGLST,
     :             INEW, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      X = 0.0
      Y = 0.0
      Z = 0.0
C
      DO I = 1, 4
         X = X + NODLST(1,NODS(I))
         Y = Y + NODLST(2,NODS(I))
         Z = Z + NODLST(3,NODS(I))
      END DO
C
      X = X / 4
      Y = Y / 4
      Z = Z / 4
C
      INEW = -1
      INELM = INT(NODLST(16,NODS(1)))  ! start search from this elem
      CALL ADDND( BIGLST, NODLST,
     :            ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :            SZENLS, SZNELS, NNOD,   NELM,
     :            NODX,   NODY,   NODZ,   ORGMTX,
     :            X, Y, Z, INOD,  INEW,   INELM,  .TRUE. )
C
      CALL STNDFL( BIGLST, NODLST, INOD, 3, .TRUE. )
C
      IF( IERR .NE. 0 ) RETURN
C
      DO I = 1, 6
         EDGFLG(I) = 0
      END DO
C
      DO I = 1, 4
C
         NDS(1) = NODS(I)
         NDS(2) = INOD
         INEW = -1
         CALL ADDEDG( BIGLST, NODLST, NDS, NEWEDS(I), INEW )
         CALL SETFLG( BIGLST, NODLST, NEWEDS(I), 3, .TRUE. )
         if( goterr ) print*,'   added edge ',neweds(i),nds
C
         IF( CONELM(I) .GT. 0 ) THEN
C
            CALL CMLCED( BIGLST, NODLST, IELM, CONELM(I), EDGLST,
     :                   E1, E2, E3 )
C
            IF( E1 .EQ. 0 ) THEN
               PRINT*,'*** ELMDIV: LOST COMMON EDGES (1)'
               PRINT*,'   elements: ',ielm,conelm(i)
               print*,'   edglst: ',edglst
               print*,i,e1,e2,e3
               print*,biglst(1,e1),biglst(2,e1)
               print*,biglst(1,e2),biglst(2,e2)
               print*,biglst(1,e3),biglst(2,e3)
               ierr = -1010
               RETURN
            END IF
C
            IF( E2 .EQ. 0 ) THEN
               PRINT*,'*** ELMDIV: LOST COMMON EDGES (2)'
               PRINT*,'   elements: ',ielm,conelm(i)
               print*,'   edglst: ',edglst
               print*,'   edge 1: ',e1
               ierr = -1010
               RETURN
            END IF
C
            IF( E3 .EQ. 0 ) THEN
               PRINT*,'*** ELMDIV: LOST COMMON EDGES (3)'
               PRINT*,'   elements: ',ielm,conelm(i)
               print*,'   edglst: ',edglst
               print*,'   edge12: ',e1,e2
               ierr = -1010
               RETURN
            END IF
C
            IF( EDGFLG(E1) .GT. 0 ) THEN
               EDGFLG(E1) = -I
            ELSE IF( EDGFLG(E1) .EQ. 0 ) THEN
               EDGFLG(E1) = I
            ELSE
               PRINT*,'*** ELMDIV: EDGE COUNTED THREE TIMES'
               IERR = -1005
               RETURN
            END IF
C
            IF( EDGFLG(E2) .GT. 0 ) THEN
               EDGFLG(E2) = -I
            ELSE IF( EDGFLG(E2) .EQ. 0 ) THEN
               EDGFLG(E2) = I
            ELSE
               PRINT*,'*** ELMDIV: EDGE COUNTED THREE TIMES'
               IERR = -1005
               RETURN
            END IF
C
            IF( EDGFLG(E3) .GT. 0 ) THEN
               EDGFLG(E3) = -I
            ELSE IF( EDGFLG(E3) .EQ. 0 ) THEN
               EDGFLG(E3) = I
            ELSE
               PRINT*,'*** ELMDIV: EDGE COUNTED THREE TIMES'
               IERR = -1005
               RETURN
            END IF
C
            CMNED(I*3-2) = EDGLST(E1)
            CMNED(I*3-1) = EDGLST(E2)
            CMNED(I*3)   = EDGLST(E3)
C
         END IF
C
      END DO
C
      IF( IERR .NE. 0 ) RETURN
C
      DO I = 1, 4
C
         IF( CONELM(I) .LT. 0 ) THEN
C
            IF( EDGFLG(1) .GE. 0 ) THEN
               E1 = 1
            ELSE IF( EDGFLG(2) .GE. 0 ) THEN
               E1 = 2
            ELSE IF( EDGFLG(3) .GE. 0 ) THEN
               E1 = 3
            ELSE IF( EDGFLG(4) .GE. 0 ) THEN
               E1 = 4
            ELSE
               PRINT*,'*** ELMDIV: LOST EDGE'
               IERR = -1002
               RETURN
            END IF
C
            IF( EDGFLG(E1) .EQ. 0 ) EDGFLG(E1) = I
C
            E2 = E1
C
 100        IF( E2 .LT. 5 ) THEN
               E2 = E2 + 1
               IF( EDGFLG(E2) .LT. 0 ) GOTO 100
               IF( EDGFLG(E2) .EQ. EDGFLG(E1) ) GOTO 100
            ELSE
               PRINT*,'*** ELMDIV: LOST EDGE'
               IERR = -1002
               RETURN
            END IF
C
            IF( EDGFLG(E2) .EQ. 0 ) EDGFLG(E2) = I
C
            E3 = E2
C
 200        IF( E3 .LT. 6 ) THEN
               E3 = E3 + 1
               IF( EDGFLG(E3) .LT. 0 ) GOTO 200
               IF( EDGFLG(E3) .EQ. EDGFLG(E1) ) GOTO 100
               IF( EDGFLG(E3) .EQ. EDGFLG(E2) ) GOTO 100
            ELSE
               PRINT*,'*** ELMDIV: LOST EDGE'
               IERR = -1002
               RETURN
            END IF
C
            IF( EDGFLG(E1) .NE. I ) EDGFLG(E1) = -I
            IF( EDGFLG(E2) .NE. I ) EDGFLG(E2) = -I
C
            IF( EDGFLG(E3) .GT. 0 ) THEN
               EDGFLG(E3) = -I
            ELSE IF( EDGFLG(E3) .EQ. 0 ) THEN
               EDGFLG(E3) = I
            ELSE
               PRINT*,'*** ELMDIV: EDGE COUNTED THREE TIMES'
               IERR = -1005
               RETURN
            END IF
C
            CMNED(I*3-2) = EDGLST(E1)
            CMNED(I*3-1) = EDGLST(E2)
            CMNED(I*3)   = EDGLST(E3)
C
         END IF
C
      END DO
c
      if( goterr ) then
         print*,'Testing split of element ',ielm
         print*,'   edglst: ',edglst
         print*,'    nodes: ',nods
         print*,'   cmned1: ',(cmned(i),i=1,3)
         print*,'   cmned2: ',(cmned(i),i=4,6)
         print*,'   cmned3: ',(cmned(i),i=7,9)
         print*,'   cmned4: ',(cmned(i),i=10,12)
         print*,'   edgflg: ',edgflg
      end if
C
      CALL RMVELE( BIGLST, NODLST, IELM )
C
      IF( IERR .NE. 0 ) RETURN
C
      DO I = 1, 4
C
         E1 = CMNED(I*3-2)
         E2 = CMNED(I*3-1)
         E3 = CMNED(I*3)
C
         CALL CMNNOD( BIGLST, NODLST, E1, E2, NCMN, N1, N2 )
c
         if( ncmn .eq. 0 ) then
            print*,'*** ELMDIV: LOST COMMON NODES (cmn)'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1011
            RETURN
         END IF
c
         if( n1 .eq. 0 ) then
            print*,'*** ELMDIV: LOST COMMON NODES (1)'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1011
            RETURN
         END IF
c
         if( n2 .eq. 0 ) then
            print*,'*** ELMDIV: LOST COMMON NODES (2)'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1011
            RETURN
         END IF
c
         if( biglst(1,e3) .ne. n1 .and. biglst(2,e3) .ne. n1 ) then
            print*,'*** ELMDIV: NODE 1 NOT ON THIRD EDGE'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1021
            RETURN
         END IF
c
         if( biglst(1,e3) .ne. n2 .and. biglst(2,e3) .ne. n2 ) then
            print*,'*** ELMDIV: NODE 2 NOT ON THIRD EDGE'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1021
            RETURN
         END IF
C
         P0 = 0
         P1 = 0
         P2 = 0
C
         DO J = 1, 4
            IF( NCMN .EQ. NODS(J) ) THEN
               P0 = J
            ELSE IF( N1 .EQ. NODS(J) ) THEN
               P1 = J
            ELSE IF( N2 .EQ. NODS(J) ) THEN
               P2 = J
            END IF
         END DO
C
         IF( P0 .EQ. 0 ) THEN
            PRINT*,'*** ELMDIV: LOST FIRST NODE'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1032
            RETURN
         END IF
C
         IF( P1 .EQ. 0 ) THEN
            PRINT*,'*** ELMDIV: LOST SECOND NODE'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1032
            RETURN
         END IF
C
         IF( P2 .EQ. 0 ) THEN
            PRINT*,'*** ELMDIV: LOST THIRD NODE'
            print*,i,e1,e2,e3
            print*,biglst(1,e1),biglst(2,e1)
            print*,biglst(1,e2),biglst(2,e2)
            print*,biglst(1,e3),biglst(2,e3)
            IERR = -1032
            RETURN
         END IF
C
         CONEL1(1) = CONELM(I)
         CONEL1(2) = 0
         CONEL1(3) = 0
         CONEL1(4) = 0
C
         EDGLS1(1) = E1
         EDGLS1(2) = E2
         EDGLS1(3) = E3
         EDGLS1(4) = NEWEDS(P0)
         EDGLS1(5) = NEWEDS(P1)
         EDGLS1(6) = NEWEDS(P2)
C
         INEW = -1
         CALL ADDELE( BIGLST, NODLST, CONEL1, EDGLS1, NEWELS(I), INEW )
C
         IF( IERR .NE. 0 ) RETURN
C
      END DO
C
      CALL CNCTEL( BIGLST, NODLST, NEWELS(1), NEWELS(2) )
      CALL CNCTEL( BIGLST, NODLST, NEWELS(1), NEWELS(3) )
      CALL CNCTEL( BIGLST, NODLST, NEWELS(1), NEWELS(4) )
      CALL CNCTEL( BIGLST, NODLST, NEWELS(2), NEWELS(3) )
      CALL CNCTEL( BIGLST, NODLST, NEWELS(2), NEWELS(4) )
      CALL CNCTEL( BIGLST, NODLST, NEWELS(3), NEWELS(4) )
C
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
         CALL ELMINF( BIGLST, NODLST, I, CONELM, EDGLST,
     :                INEW, IFLAG, NXT, LST )
         DO J = 1, 4
            IF( CONELM(J) .EQ. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** ELMDIV: INCONSISTENT CONNECTED ELEMENTS'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,conelm(j)) .le. 0 ) then
               if( .not. addit ) then
                  print*,'*** ELMDIV: INVALID CONNECTED ELEMENT'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
               PRINT*,'   invalid element: ',
     :                      conelm(j),biglst(4,conelm(j))
            end if
         end do
         call elnods( BIGLST, NODLST, i, nods, nxt, lst, .false. )
         d1 = elmfnc( BIGLST, NODLST, i,
     :                nods(1), nods(2), nods(3), nods(4), d2 )
         if( abs(d2) .lt. 1e-8 ) then
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
            nods(1) = biglst(1,edglst(j))
            nods(2) = biglst(2,edglst(j))
            IF( edglst(J) .LE. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** ELMDIV: INCONSISTENT EDGE LIST'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,edglst(j)) .ge. 0 ) then
               if( .not. addit ) then
                  PRINT*,'*** ELMDIV: INVALID EDGE'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
               PRINT*,'   invalid edge: ',edglst(j),biglst(4,edglst(j))
            else if( getflg(BIGLST,NODLST,edglst(j),2) ) then
               if( .not. gtndfl(BIGLST,NODLST,nods(1),2) .and.
     :             .not. gtndfl(BIGLST,NODLST,nods(1),4) ) then
                  if( .not. addit ) then
                    PRINT*,'*** ELMDIV: INVALID GEOM-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
               else if(.not. gtndfl(BIGLST,NODLST,nods(2),2) .and.
     :                 .not. gtndfl(BIGLST,NODLST,nods(2),4) ) then
                  if( .not. addit ) then
                    PRINT*,'*** ELMDIV: INVALID GEOM-EDGE NODE 2 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
               end if
            else if( .not. getflg(BIGLST,NODLST,edglst(j),3) ) then
               if( gtndfl(BIGLST,NODLST,nods(1),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** ELMDIV: INVALID SURF-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
               else if( gtndfl(BIGLST,NODLST,nods(2),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** ELMDIV: INVALID SURF-EDGE NODE 2 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
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
         print*,'*** ELMDIV: INCONSISTENCY AT END OF BIGLST'
         print*,i,emtbig
         ierr = -10000
      end if
C
      RETURN
      END
C
