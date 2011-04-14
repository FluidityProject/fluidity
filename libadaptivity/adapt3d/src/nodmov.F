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
      SUBROUTINE NODMOV( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX, ONCE,
     :                   INOD, IELM, XN, YN, ZN, INELM )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INTEGER INOD, IELM, INELM
C
      REAL FNCDIF, FNCORG, XN, YN, ZN
C
      LOGICAL ONCE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkorg.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      INTEGER ELPTRS(MAXEDE), EDPTRS(MAXEDE), NDPTRS(MAXEDE), I, J, K,
     :        NUMEDE, NUMELM, NODS(4), INEW, JNOD, CONELM(4), EDGLST(6),
     :        IFLAG, NXT, LST, IREG, nodt(4)
C
      REAL F0, F1, F2, F3, X(MAXEDE), Y(MAXEDE), Z(MAXEDE), D1, D2, D3,
     :     V1(3), V2(3), V3(3), XO, YO, ZO, XX(2), YY(2), ZZ(2), VOL1,
     :     VOL2, dx, dy, elmfnc, edgfnc, t, v(maxede), rad, xt, yt, zt
C
      LOGICAL ISINSD, ISSPLT, GTNDFL, GETFLG, goterr, addit
C
      CALL ADAPT_SECOND( TIMSTT )
C
      FNCORG = 1E+30
      FNCDIF = -1E+30
C
      IF( GTNDFL(BIGLST,NODLST,INOD,2) ) THEN
         ewrite(-1,*) 
     :        '*** NODMOV: attempt to relocate geometry node'
         ierr = -10000
         return
      end if
C
      CALL NDCNCT( BIGLST, NODLST, INOD, IELM,
     :             ELPTRS, EDPTRS, NDPTRS, NUMEDE, NUMELM )
C
      IF( IERR .NE. 0 ) RETURN
C
c      CALL NODEIG( BIGLST, NODLST, INOD, V1, V2, V3, D1, D2, D3 )
C
      XO = NODLST( 1, INOD )
      YO = NODLST( 2, INOD )
      ZO = NODLST( 3, INOD )
C
c      INELM = INT(NODLST(16,INOD))
C
c      XX(1) = 0.0
c      YY(1) = 0.0
c      ZZ(1) = 0.0
C
c      ISSPLT = GTNDFL( BIGLST, NODLST, INOD, 4 )
c      ISINSD = GTNDFL( BIGLST, NODLST, INOD, 3 ) .AND. .NOT. ISSPLT
C
c      I = 0
C
c      IF( ISINSD ) THEN
C
c         DO J = 1, NUMEDE
c            XX(2) = NODLST(1,NDPTRS(J)) - XO
c            YY(2) = NODLST(2,NDPTRS(J)) - YO
c            ZZ(2) = NODLST(3,NDPTRS(J)) - ZO
c            CALL EDGEIG( INOD, NDPTRS(J), V1, V2, V3, D1, D2, D3 )
cc            CALL UNDSTT( XX, YY, ZZ, 2, V1, V2, V3, D1, D2, D3 )
c            if( tstdst ) t = edgfnc( ndptrs(j), inod, v(j) )
c            X(J) = XX(2)
c            Y(J) = YY(2)
c            Z(J) = ZZ(2)
c         END DO
C
c         I = NUMEDE
C
c      ELSE IF( ISSPLT ) THEN
C
c         DO J = 1, NUMEDE
c            IF( GETFLG(EDPTRS(J),2) ) THEN
c               I = I + 1
c               XX(2) = NODLST(1,NDPTRS(J)) - XO
c               YY(2) = NODLST(2,NDPTRS(J)) - YO
c               ZZ(2) = NODLST(3,NDPTRS(J)) - ZO
c               CALL EDGEIG( INOD, NDPTRS(J), V1, V2, V3, D1, D2, D3 )
cc               CALL UNDSTT( XX, YY, ZZ, 2, V1, V2, V3, D1, D2, D3 )
c               if( tstdst ) t = edgfnc( ndptrs(j), inod, v(i) )
c               X(I) = XX(2)
c               Y(I) = YY(2)
c               Z(I) = ZZ(2)
c            END IF
c         END DO
C
c         IF( I .NE. 2 ) THEN
c            PRINT*,'*** NODMOV: GEOM-SPLIT NODE NOT ON TWO EDGES'
c            IERR = 601
c            RETURN
c         END IF
C
c      ELSE
C
c         DO J = 1, NUMEDE
c            IF( .NOT. GETFLG(EDPTRS(J),3) ) THEN
c               I = I + 1
c               XX(2) = NODLST(1,NDPTRS(J)) - XO
c               YY(2) = NODLST(2,NDPTRS(J)) - YO
c               ZZ(2) = NODLST(3,NDPTRS(J)) - ZO
c               CALL EDGEIG( INOD, NDPTRS(J), V1, V2, V3, D1, D2, D3 )
cc               CALL UNDSTT( XX, YY, ZZ, 2, V1, V2, V3, D1, D2, D3 )
c               if( tstdst ) t = edgfnc( ndptrs(j), inod, v(i) )
c               X(I) = XX(2)
c               Y(I) = YY(2)
c               Z(I) = ZZ(2)
c            END IF
c         END DO
C
c         IF( I .LT. 3 ) THEN
c            PRINT*,'*** NODMOV: SURFACE NODE ON LESS THAN 3 EDGES'
c            IERR = 601
c            RETURN
c         END IF
C
c      END IF
C
c      IF( I .LE. 0 ) THEN
c         PRINT*,'*** NODMOV: BAD NUMBER OF CONNECTIONS'
c         PRINT*,INOD,I,NUMEDE,ISSPLT,ISINSD
c         IERR = -2000
c         RETURN
c      END IF
C
c      if( .not. tstdst ) then
c         do j = 1, i
c            v(j) = 1.0
c         end do
c      end if
c
c      T = 0.0
C
c      DO J = 1, I
C
c         XX(1) = XX(1) + X(J)*V(J)
c         YY(1) = YY(1) + Y(J)*V(J)
c         ZZ(1) = ZZ(1) + Z(J)*V(J)
C
c         T = T + V(J)
C
c      END DO
C
c      XX(1) = XX(1)/T + XO
c      YY(1) = YY(1)/T + YO
c      ZZ(1) = ZZ(1)/T + ZO
C
c      IFLAG = INT(NODLST(6,INOD))
c      CALL RMVND( BIGLST, NODLST, INOD )
c      INEW = -1
c      CALL ADDND( BIGLST, NODLST, XN, YN, ZN, JNOD, INEW )
c      NODLST(6,JNOD) = FLOAT(IFLAG)
C
      xt = xn
      yt = yn
      zt = zn
c
      CALL FNDELM( ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :             SZENLS, SZNELS, NNOD,   NELM,
     :             NODX,   NODY,   NODZ,   ORGMTX,
     :             XN, YN, ZN,     INELM,  .FALSE., NODLST(7,INOD) )
c
      if( inelm .lt. 0 ) then
         ewrite(1,*) 
     :        '--+++ NODMOV: Moved node outside mesh!! +++--'
         ewrite(1,*) xo,yo,zo
         ewrite(1,*) xt,yt,zt
         ewrite(1,*) xn,yn,zn
      end if
c
      NODLST(1,INOD) = XN
      NODLST(2,INOD) = YN
      NODLST(3,INOD) = ZN
c
      NODLST(16,INOD) = FLOAT(ABS(INELM))
C      CALL FIXEIG( NODLST(1,INOD), NODLST(7,INOD), ANISOT )
C
c      if( inod .ne. jnod ) then
c         print*,'*** NODMOV: GOT DIFFERENT NODE NUMBERS'
c         ierr = -610
c      end if
C
c      IF( IERR .NE. 0 ) RETURN
C
      if( .not. once ) then
       DO I = 1, NUMELM
C
         CALL ELMINF( BIGLST, NODLST, ELPTRS(I), CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
C
         if( debug ) then
            debug = .false.
            call elnods( BIGLST, NODLST, ELPTRS(I), nodt, nxt, lst,
     :                  .FALSE. )
            ewrite(3,'(a,i8,1PE13.5)')
     :           '   removing elm ',ELPTRS(I), elmfnc( BIGLST, NODLST,
     :           ELPTRS(I),nodt(1),nodt(2),nodt(3),nodt(4),rad)
            ewrite(3,'(a,4i8)')
     :           '     nodes: ',(nodt(ierr),ierr=1,4)
            ierr = 0
            debug = .true.
         end if
         CALL RMVELE( BIGLST, NODLST, ELPTRS(I) )
C
         INEW = -1
         CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, ELPTRS(I),
     :                INEW, IREG )
c
         if( ierr.eq.12345 ) then
           ewrite(0,*)
     :            '--- inconsistent oppedg from addele in nodmov (1)'
           ierr = 0
         end if
c
         if( debug ) then
            debug = .false.
            call elnods( BIGLST, NODLST, ELPTRS(I), nodt, nxt, lst,
     :                  .FALSE. )
            ewrite(3,*)
     :        '   added elm ',ELPTRS(I), elmfnc( BIGLST, NODLST,
     :                ELPTRS(I),nodt(1),nodt(2),nodt(3),nodt(4),rad)
            ewrite(3,*)
     :              '     nodes: ',(nodt(ierr),ierr=1,4)
            ierr = 0
            biglst(7,ELPTRS(I)) = 0
            nxt = biglst(nxtbig,ELPTRS(I))
            biglst(7,nxt) = 0
            debug = .true.
         end if
C
         IF( IERR .NE. 0 ) RETURN
C
       END DO
      end if
c
c      do i = 1, numede
c         call stndfl( BIGLST, NODLST, ndptrs(i), 6, .false. )
c      end do
c
      CALL ADAPT_SECOND( TIMEND )
C
      TIMCHG = TIMCHG + TIMEND - TIMSTT
      TOTCHG = TOTCHG + 1
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
     :                INEW, IFLAG, IREG, NXT, LST )
         DO J = 1, 4
            IF( CONELM(J) .EQ. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** NODMOV: INCONSISTENT CONNECTED ELEMENTS'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,conelm(j)) .le. 0 ) then
               if( .not. addit ) then
                  print*,'*** NODMOV: INVALID CONNECTED ELEMENT'
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
                  PRINT*,'*** NODMOV: INCONSISTENT EDGE LIST'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,edglst(j)) .ge. 0 ) then
               if( .not. addit ) then
                  PRINT*,'*** NODMOV: INVALID EDGE'
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
                    PRINT*,'*** NODMOV: INVALID GEOM-EDGE NODE 1 FLAGS'
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
                    PRINT*,'*** NODMOV: INVALID GEOM-EDGE NODE 2 FLAGS'
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
                    PRINT*,'*** NODMOV: INVALID SURF-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
               else if( gtndfl(BIGLST,NODLST,nods(2),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** NODMOV: INVALID SURF-EDGE NODE 2 FLAGS'
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
         print*,'*** NODMOV: INCONSISTENCY AT END OF BIGLST'
         print*,i,emtbig
         ierr = -10000
      end if
C
      RETURN
      END
C
