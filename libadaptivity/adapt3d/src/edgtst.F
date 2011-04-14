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
      SUBROUTINE EDGTST( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   SPLTED, CLPSED, ED2FC, MVNODE,
     :                   IEDG, IELM, ELPTRS, EDPTRS, NDPTRS,
     :                   NUMEDE, SURFAC, INTRNL,
     :                   ELPTS2, EDPTS2, NDPTS2,
     :                   NUMEL1, NUMEL2, NUMED1, NUMED2,
     :                   FNCDIF, FNCORG, IWHAT,
     :                   AVEDIF, AVEORG, AWHAT,
     :                   BSTRAT, BSTX,  BSTY,  BSTZ,  BSTINE,
     :                   BSTRAA, BSTXA, BSTYA, BSTZA, BSTINA,
     :                   ALWAYS, DOTOP, MINCHG )
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
      INCLUDE 'blkorg.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      INCLUDE 'badnod.i'
C
      INTEGER IEDG, IELM, IWHAT, ELPTRS(3*MAXEDE), EDPTRS(3*MAXEDE+3),
     :        NDPTRS(MAXEDE+1), NUMEDE, SURFAC, INTRNL, AWHAT, BSTINE,
     :        BSTINA
C
      REAL BSTX, BSTY, BSTZ, FNCDIF, FNCORG, AVEDIF, AVEORG, BSTRAT,
     :     ALWAYS, DOTOP, MINCHG, BSTXA, BSTYA, BSTZA, BSTRAA
C
      LOGICAL MVNODE, ED2FC, SPLTED, CLPSED
C
      INTEGER ELPTS2(2*MAXEDE), EDPTS2(2*MAXEDE), NDPTS2(2*MAXEDE),
     :        NUMEL1, NUMEL2, NUMED1, NUMED2
C
      INTEGER NOD1, NOD2, NODS(4), I, J, INEW, STRTND, LASTND, NUMEDG,
     :        LSTED1, LSTED2, NXT, LST,
     :        INOD, NUMELS, NODO1, NODO2, CNTLFT, INELM, CNT1, CNT2,
     :        CNT3
C
      REAL XX, YY, ZZ, F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, FA, FB,
     :     SPLFNC, SHRFNC, E2FFNC, EDGFNC, ELMFNC,
     :     VOL1, VOL2, UNDLEN, ORGLEN, ORGRDS, U1, U2,
     :     FUP, FLW, OX, OY, OZ, a,b,c,d,e, vol3, vol4,
     :     SPLAVE, SHRAVE, E2FAVE, E2EAVE, RATIO
C
      LOGICAL ISGEOM, ISINSD, GTNDFL, GETFLG, didadd
c
c      debug = .false.
c      if( iedg .eq. 7825 ) debug = .true.
C
      IWHAT  = 0
      AWHAT  = 0
      RATIO  = BSTRAT
C
      FNCDIF = 1E+31
      FNCORG = 1E+32
      AVEDIF = 1E+31
      AVEORG = 1E+32
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT ERROR ON ENTRY'
         RETURN
      END IF
C
      IF( GETFLG( BIGLST, NODLST, IEDG, 6 ) ) RETURN
C
cc      IF( IEDG .LE. 0 ) THEN
cc         PRINT*,'*** EDGTST: INVALID EDGE POINTER'
cc         IERR = -2
cc         RETURN
cc      END IF
C
cc      IF( BIGLST(4,IEDG) .EQ. 0 ) THEN
cc         PRINT*,'*** EDGTST: BLOCK IS EMPTY'
cc         IERR = -3
cc         RETURN
cc      ELSE IF( BIGLST(4,IEDG) .GT. 0 ) THEN
cc         PRINT*,'*** EDGTST: BLOCK IS AN ELEMENT'
cc         IERR = -3
cc         RETURN
cc      END IF
C
      ISGEOM = GETFLG( BIGLST, NODLST, IEDG, 2 )
      ISINSD = GETFLG( BIGLST, NODLST, IEDG, 3 ) .AND. .NOT. ISGEOM
C
      IF( SURFAC .GT. 0 ) THEN
         IF( ISINSD ) THEN
            PRINT*,'*** EDGTST: SURFACE EDGE INCONSISTENCY'
            print*,'   els: ',(elptrs(i*3-2),i=1,numede)
            print*,'   eds: ',(edptrs(i*3),i=1,numede)
            print*,'   nds: ',(ndptrs(i),i=1,numede+1)
            print*,'   up els: ',(elptrs(i*3-1),i=1,numede)
            print*,'   lw els: ',(elptrs(i*3),i=1,numede)
            print*,'   up eds: ',(edptrs(i*3-2),i=1,numede+1)
            print*,'   lw eds: ',(edptrs(i*3-1),i=1,numede+1)
            IERR = -111
            RETURN
         END IF
         NUMEDG = NUMEDE+1
      ELSE
         NUMEDG = NUMEDE
      END IF
C
      STRTND = NDPTRS(NUMEDG)
C
      NOD1 = BIGLST(1,IEDG)
      NOD2 = BIGLST(2,IEDG)
C
      IF( NOD1 .LE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT INVALID NODE FROM EDGE:'
         PRINT*,NOD1,NOD2,IEDG,BIGLST(4,IEDG)
         IERR = -10102
         RETURN
      ELSE IF( NOD2 .LE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT INVALID NODE FROM EDGE:'
         PRINT*,NOD1,NOD2,IEDG,BIGLST(4,IEDG)
         IERR = -10102
         RETURN
      END IF
c
      xx = nodlst(1,nod1)
      yy = nodlst(2,nod1)
      zz = nodlst(3,nod1)
      if( abs(xx-badndx).lt.1e-5.and.abs(yy-badndy).lt.1e-5 ) then
         print*,'*** EDGTST: Found bad node (1) !!'
         print*,nod1,xx,yy,zz
         print*,nod2,nodlst(1,nod2),nodlst(2,nod2),nodlst(3,nod2)
      end if
c
      xx = nodlst(1,nod2)
      yy = nodlst(2,nod2)
      zz = nodlst(3,nod2)
      if( abs(xx-badndx).lt.1e-5.and.abs(yy-badndy).lt.1e-5 ) then
         print*,'*** EDGTST: Found bad node (2) !!'
         print*,nod2,xx,yy,zz
         print*,nod1,nodlst(1,nod1),nodlst(2,nod1),nodlst(3,nod1)
      end if
c
      didadd = .false.
c
c      if( debug ) print*,'Testing edge: ',iedg,nod1,nod2
C
C - this is the functional for the elements surrounding the edge
C
      F0 = 0.0
      F1 = 0.0
      AVEORG = 0.0
c      ORGRDS = 1E+30
C
      LASTND = STRTND
C
      DO I = 1, NUMEDE
         F0 = ELMFNC( BIGLST, NODLST, ELPTRS(I*3-2),
     :                NOD1, NOD2, NDPTRS(I), LASTND, VOL1 )
         AVEORG = AVEORG + F0
         F1 = MAX( F1, F0 )
c         if( debug ) print*,'   old in-sphere: ',vol1
c         ORGRDS = MIN(ORGRDS,ABS(VOL1))
         LASTND = NDPTRS(I)
         IF( GETFLG(BIGLST,NODLST,ELPTRS(I*3-2),5) ) THEN
            AVEORG = 1E+32
            GOTO 501
         END IF
      END DO
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT ERROR FROM SURROUNDING FNCLS'
         print*,'   ',iedg,nod1,nod2,ielm
         print*,'   ',numede,numedg,surfac
         print*,'   els: ',(elptrs(i*3-2),i=1,numede)
         print*,'   eds: ',(edptrs(i*3),i=1,numede)
         print*,'   nds: ',(ndptrs(i),i=1,numedg)
         print*,'   up els: ',(elptrs(i*3-1),i=1,numede)
         print*,'   lw els: ',(elptrs(i*3),i=1,numede)
         print*,'   up eds: ',(edptrs(i*3-2),i=1,numedg)
         print*,'   lw eds: ',(edptrs(i*3-1),i=1,numedg)
         RETURN
      END IF
C
      AVEORG = AVEORG / NUMEDE
c
      if( f1 .le. dotop ) goto 401
C
c      IF( ORGRDS .GT. 0.1 ) THEN
c         ORGRDS = ORGRDS*0.5
c      ELSE
c         ORGRDS = ORGRDS*1.1
c      END IF
c      ORGRDS = 0.01
C
 100  XX = NODLST(1,NOD2)*RATIO + NODLST(1,NOD1)*(1.0-RATIO)
      YY = NODLST(2,NOD2)*RATIO + NODLST(2,NOD1)*(1.0-RATIO)
      ZZ = NODLST(3,NOD2)*RATIO + NODLST(3,NOD1)*(1.0-RATIO)
C
      INEW = -1
      INELM = INT(NODLST(16,NOD1))
      CALL ADDND( BIGLST, NODLST,
     :            ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :            SZENLS, SZNELS, NNOD,   NELM,
     :            NODX,   NODY,   NODZ,   ORGMTX,
     :            XX, YY, ZZ, INOD, INEW, INELM, .TRUE. )
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT ERROR FROM ADDND'
         RETURN
      END IF
c
      if( abs(xx-badndx).lt.1e-5.and.abs(yy-badndy).lt.1e-5) then
         print*,'*** EDGTST: ABOUT TO ADD BAD NODE (1)'
         print*, xx, yy, zz, ratio
         print*,NODLST(1,NOD1),NODLST(2,NOD1),NODLST(3,NOD1)
         print*,NODLST(1,NOD2),NODLST(2,NOD2),NODLST(3,NOD2)
         print*,'   els: ',(elptrs(i*3-2),i=1,numede)
         print*,'   eds: ',(edptrs(i*3),i=1,numede)
         print*,'   nds: ',(ndptrs(i),i=1,numede+1)
c         print*,'   up els: ',(elptrs(i*3-1),i=1,numede)
c         print*,'   lw els: ',(elptrs(i*3),i=1,numede)
c         print*,'   up eds: ',(edptrs(i*3-2),i=1,numede+1)
c         print*,'   lw eds: ',(edptrs(i*3-1),i=1,numede+1)
         do i = 1, numede+1
         print*,'   nd',i,' : ',(nodlst(j,ndptrs(i)),j=1,3)
         end do
         do i = 1, numede
            j = elptrs(i*3-2)
            call elnods(biglst,nodlst,j,nods,nxt,lst,.false.)
            print*,'   el',j,nods
         end do
      end if
c
      didadd = .true.
c
      if( inelm .lt. 0 ) then
         print*,'---+++ EDGTST(1): Got node outside mesh +++---'
         inelm = -inelm
         if( isgeom ) then
            print*,'  A geometry edge: ',ratio
         else if( isinsd ) then
            print*,'  An internal edge: ',ratio
         else
            print*,'  A plane surface edge: ',ratio
         end if
         print*,'node 1: ',(nodlst(i,nod1),i=1,3)
         print*,'node 2: ',(nodlst(i,nod2),i=1,3)
      end if
C
      BSTX   = XX
      BSTY   = YY
      BSTZ   = ZZ
      BSTRAT = RATIO
      BSTXA  = XX
      BSTYA  = YY
      BSTZA  = ZZ
      BSTRAA = RATIO
      BSTINE = INELM
      BSTINA = INELM
C
c      IF( .NOT. SPLTED .OR. ORGLEN .LT. 0.2 ) GOTO 301
      IF( .NOT. SPLTED ) GOTO 301
C
      CALL ADAPT_SECOND( TIMSTT )
C
      LASTND = STRTND
C
C - work out the functional of the new elements around the new edges
C
      F0 = 0.0
      F3 = 0.0
      F4 = 0.0
C
      DO I = 1, NUMEDE
         F0 = ELMFNC( BIGLST, NODLST, 0,
     :                NOD1, INOD, NDPTRS(I), LASTND, VOL1 )
         F4 = F4 + F0
         F3 = MAX( F0, F3 )
         F0 = ELMFNC( BIGLST, NODLST, 0,
     :                NOD2, INOD, NDPTRS(I), LASTND, VOL2 )
         F4 = F4 + F0
         F3 = MAX( F0, F3 )
         LASTND = NDPTRS(I)
c        if( debug ) then
c           print*,'  new in-spheres: ',vol1,vol2
c        end if
        IF( VOL1/VOL2 .GT. -1E-4 ) GOTO 300
c         IF( ABS(VOL1) .LT. ORGRDS ) GOTO 300
c         IF( ABS(VOL2) .LT. ORGRDS ) GOTO 300
      END DO
C
C - find the difference in the functional for the edge split case
C
      FNCORG = F1
      SPLFNC = F3 - FNCORG
      SPLAVE = F4 / NUMEDE / 2 - AVEORG
      FNCDIF = SPLFNC
      IF( FNCDIF .LT. 0.0 ) AVEDIF = SPLAVE
C
      IWHAT = 1
      AWHAT = 1
C
 300  CONTINUE
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMSPL = TIMSPL + TIMEND - TIMSTT
      TOTSPL = TOTSPL + 1
C
C - this is an attempt to cut down the time by accepting an edge split
C - straight away if it makes a significant improvement; thus avoiding
C - checking edge-to-face, edge-to-edge and edge collapsing.
C
      IF( FNCDIF .LT. ALWAYS ) GOTO 501
      IF( AVEDIF .LT. ALWAYS ) GOTO 501
      if( fncdif+fncorg .le. dotop ) goto 501
C
 301  IF( .NOT. ED2FC ) GOTO 401
C
C - The case of edge to face swapping can only be checked if there
C - are three or four elements round an edge and there is no surface.
C - In theory we could also do an edge swap for two elements around
C - a non-geometry surface edge - later, perhaps!
C
      IF( SURFAC .GT. 0 ) GOTO 401
C
      IF( INTRNL .LT. 0 ) GOTO 401
C
      IF( NUMEDE .NE. 3 ) GOTO 351
C
      CALL ADAPT_SECOND( TIMSTT )
C
C - we already have the functional of the surrounding elements from
C - the previous check (stored in F1). Now find the functional for the
C - two elements that would replace them.
C
      F0 = ELMFNC( BIGLST, NODLST, 0,
     :             NOD1, NDPTRS(1), NDPTRS(2), NDPTRS(3), VOL1 )
      F4 = ELMFNC( BIGLST, NODLST, 0,
     :             NOD2, NDPTRS(1), NDPTRS(2), NDPTRS(3), VOL2 )
      F3 = MAX( F0, F4 )
      F4 = ( F4 + F0 ) / 2
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT ERROR FROM ELMFNC'
         RETURN
      END IF
C
C - find the difference in the functional for the edge to face case
C
      IF( VOL1*VOL2 .GE. 0.0 ) THEN
         GOTO 400
c      ELSE IF( ABS(VOL1) .LT. ORGRDS ) THEN
c         GOTO 400
c      ELSE IF( ABS(VOL2) .LT. ORGRDS ) THEN
c         GOTO 400
      ELSE
         E2FFNC = F3 - F1
         E2FAVE = F4 - AVEORG
      END IF
C
      IF( E2FFNC .LT. FNCDIF ) THEN
         IWHAT  = 3
         FNCDIF = E2FFNC
         FNCORG = F1
      END IF
C
      IF( E2FAVE .LT. AVEDIF .AND. E2FFNC .LT. 0.0 ) THEN
         AWHAT  = 3
         AVEDIF = E2FAVE
      END IF
c
c      if( debug ) then
c        print*,'   edge to face functionals: ',F1,F3
c        print*,nod1,nod2,ndptrs(1),ndptrs(2),ndptrs(3)
c      end if
C
 350  CONTINUE
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIME2F = TIME2F + TIMEND - TIMSTT
      TOTE2F = TOTE2F + 1
C
C - this is an attempt to cut down the time by accepting edge-to-face
C - straight away if it makes a significant improvement; thus avoiding
C - checking edge-to-edge and edge collapsing.
C
      IF( FNCDIF .LT. ALWAYS ) GOTO 501
      IF( AVEDIF .LT. ALWAYS ) GOTO 501
      if( fncdif+fncorg .le. dotop ) goto 501
C
 351  IF( NUMEDE .NE. 4 ) GOTO 401
c 351  goto 401
C
      CALL ADAPT_SECOND( TIMSTT )
C
      F0 = ELMFNC( BIGLST, NODLST, 0,
     :             NOD1, NDPTRS(1), NDPTRS(2), NDPTRS(3), VOL1 )
      F2 = ELMFNC( BIGLST, NODLST, 0,
     :             NOD1, NDPTRS(3), NDPTRS(4), NDPTRS(1), VOL2 )
      F3 = MAX( F0, F2 )
      F2 = F2 + F0
C
      IF( VOL1*VOL2 .LE. 0.0 ) THEN
         F3 = -1.0
      ELSE
C
         F0 = ELMFNC( BIGLST, NODLST, 0,
     :                NOD2, NDPTRS(1), NDPTRS(2), NDPTRS(3), VOL3 )
         F2 = F2 + F0
         F3 = MAX( F3, F0 )
C
         IF( VOL1*VOL3 .GE. 0.0 ) THEN
            F3 = -1.0
         ELSE
C
            F0 = ELMFNC( BIGLST, NODLST, 0,
     :                   NOD2, NDPTRS(3), NDPTRS(4), NDPTRS(1), VOL4 )
            F2 = ( F2 + F0 ) / 4
            F3 = MAX( F3, F0 )
C
            IF( VOL3*VOL4 .LE. 0.0 ) THEN
               F3 = -1.0
            ELSE IF( VOL2*VOL4 .GE. 0.0 ) THEN
               F3 = -1.0
            END IF
C
         END IF
C
      END IF
C
      F0 = ELMFNC( BIGLST, NODLST, 0,
     :             NOD1, NDPTRS(4), NDPTRS(1), NDPTRS(2), VOL1 )
      F4 = ELMFNC( BIGLST, NODLST, 0,
     :             NOD1, NDPTRS(2), NDPTRS(3), NDPTRS(4), VOL2 )
      F5 = MAX( F0, F4 )
      F4 = F4 + F0
C
      IF( VOL1*VOL2 .LE. 0.0 ) THEN
         F5 = -1.0
      ELSE
C
         F0 = ELMFNC( BIGLST, NODLST, 0,
     :                NOD2, NDPTRS(4), NDPTRS(1), NDPTRS(2), VOL3 )
         F4 = F4 + F0
         F5 = MAX( F5, F0 )
C
         IF( VOL1*VOL3 .GE. 0.0 ) THEN
            F5 = -1.0
         ELSE
C
            F0 = ELMFNC( BIGLST, NODLST, 0,
     :                   NOD2, NDPTRS(2), NDPTRS(3), NDPTRS(4), VOL4 )
            F4 = ( F4 + F0 ) / 4
            F5 = MAX( F5, F0 )
C
            IF( VOL3*VOL4 .LE. 0.0 ) THEN
               F5 = -1.0
            ELSE IF( VOL2*VOL4 .GE. 0.0 ) THEN
               F5 = -1.0
            END IF
C
         END IF
C
      END IF
C
      IF( F3 .GE. 0.0 ) THEN
C
         IF( F3 - F1 .LT. FNCDIF ) THEN
            FNCORG = F1
            FNCDIF = F3 - F1
            IWHAT = 4
         END IF
C
         IF( F2 - AVEORG .LT. AVEDIF .AND. F3-F1 .LT. 0.0 ) THEN
            AVEDIF = F2 - AVEORG
            AWHAT = 4
         END IF
C
      END IF
C
      IF( F5 .GE. 0.0 ) THEN
C
         IF( F5 - F1 .LT. FNCDIF ) THEN
            FNCORG = F1
            FNCDIF = F5 - F1
            IWHAT = -4
         END IF
C
         IF( F4 - AVEORG .LT. AVEDIF .AND. F5-F1 .LT. 0.0 ) THEN
            AVEDIF = F4 - AVEORG
            AWHAT = -4
         END IF
C
      END IF
c      if( debug )
c     :      print*,'   edge to edge fnls: ',f1,f3,f5
C
 400  CONTINUE
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIME2E = TIME2E + TIMEND - TIMSTT
      TOTE2E = TOTE2E + 1
C
C - this is an attempt to cut down the time by accepting edge-to-edge
C - straight away if it makes a significant improvement; thus avoiding
C - checking edge collapsing.
C
      IF( FNCDIF .LT. ALWAYS ) GOTO 501
      IF( AVEDIF .LT. ALWAYS ) GOTO 501
      if( fncdif+fncorg .le. dotop ) goto 501
C
 401  CONTINUE
C
c      IF( .NOT. CLPSED .OR. ORGLEN .GT. 5.0 ) GOTO 501
      IF( .NOT. CLPSED ) GOTO 501
C
      CALL ADAPT_SECOND( TIMSTT )
C
C - Finally, check the case of edge shrinking...
C - Before shrinking, need to see if we can, or if it has to be to one
C - of the edge's nodes rather than to the middle of the edge.
C
      NODO1 = NOD1
      NODO2 = NOD2
C
      IF( ISINSD ) THEN
         IF( .NOT. GTNDFL(BIGLST,NODLST,NOD1,3) .AND.
     :       .NOT. GTNDFL(BIGLST,NODLST,NOD2,3) ) THEN
c            print*,'+++ Warning: attempted to shrink internal edge +++'
c            print*,iedg,nod1,nod2
            GOTO 500
         ELSE IF( .NOT. GTNDFL(BIGLST,NODLST,NOD1,3) ) THEN
c            print*,'--- Info: shrinking internal edge to node 1 ---'
            if( gtndfl(BIGLST,NODLST,nod2,2) .or.
     :          gtndfl(BIGLST,NODLST,nod2,4) ) then
               print*,'*** EDGTST: INCONSISTENT NODE 2 FLAGS'
               print*,iedg,nod1,nod2
               print*,gtndfl(BIGLST,NODLST,nod1,2),
     :                gtndfl(BIGLST,NODLST,nod1,3),
     :                gtndfl(BIGLST,NODLST,nod1,4)
               print*,gtndfl(BIGLST,NODLST,nod2,2),
     :                gtndfl(BIGLST,NODLST,nod2,3),
     :                gtndfl(BIGLST,NODLST,nod2,4)
               ierr = -1000
               return
            end if
            nod1 = 0
         ELSE IF( .NOT. GTNDFL(BIGLST,NODLST,NOD2,3) ) THEN
c            print*,'--- Info: shrinking internal edge to node 2 ---'
            if( gtndfl(BIGLST,NODLST,nod1,2) .or.
     :          gtndfl(BIGLST,NODLST,nod1,4) ) then
               print*,'*** EDGTST: INCONSISTENT NODE 1 FLAGS'
               print*,iedg,nod1,nod2
               print*,gtndfl(BIGLST,NODLST,nod1,2),
     :                gtndfl(BIGLST,NODLST,nod1,3),
     :                gtndfl(BIGLST,NODLST,nod1,4)
               print*,gtndfl(BIGLST,NODLST,nod2,2),
     :                gtndfl(BIGLST,NODLST,nod2,3),
     :                gtndfl(BIGLST,NODLST,nod2,4)
               ierr = -1000
               return
            end if
            nod2 = 0
         END IF
      ELSE IF( ISGEOM ) THEN
         IF( GTNDFL(BIGLST,NODLST,NOD1,2) .AND.
     :       GTNDFL(BIGLST,NODLST,NOD2,2) ) THEN
c            print*,'+++ Warning: attempted to shrink geometry edge +++'
c            print*,iedg,nod1,nod2
            GOTO 500
         ELSE IF( GTNDFL(BIGLST,NODLST,NOD1,2) .AND.
     :            GTNDFL(BIGLST,NODLST,NOD2,4) ) THEN
            nod1 = 0
c            print*,'--- Info: shrinking geometry edge to node 1 ---'
         ELSE IF( GTNDFL(BIGLST,NODLST,NOD2,2) .AND.
     :            GTNDFL(BIGLST,NODLST,NOD1,4) ) THEN
            nod2 = 0
c            print*,'--- Info: shrinking geometry edge to node 2 ---'
         ELSE IF( GTNDFL(BIGLST,NODLST,NOD1,4) .AND.
     :            GTNDFL(BIGLST,NODLST,NOD2,4) ) THEN
c            print*,'--- INFO: shrinking a pre-split geometry edge ---'
         ELSE
            print*,'*** EDGTST: INCONSISTENT NODE FLAGS ON GEOM EDGE'
            print*,'   edge: ',iedg
            print*,'   nodes: ',nod1,nod2
            print*,'   edge flags: ',biglst(4,iedg)
            print*,'   node flags: ',nodlst(6,nod1),nodlst(6,nod2)
            print*,'   node 1 co-ords: ',(nodlst(i,nod1),i=1,3)
            print*,'   node 2 co-ords: ',(nodlst(i,nod2),i=1,3)
            print*,'   els: ',(elptrs(i*3-2),i=1,numede)
            print*,'   eds: ',(edptrs(i*3),i=1,numede)
            print*,'   nds: ',(ndptrs(i),i=1,numede+1)
            print*,'   up els: ',(elptrs(i*3-1),i=1,numede)
            print*,'   lw els: ',(elptrs(i*3),i=1,numede)
            print*,'   up eds: ',(edptrs(i*3-2),i=1,numede+1)
            print*,'   lw eds: ',(edptrs(i*3-1),i=1,numede+1)
            print*,'   node ',ndptrs(1),(nodlst(i,ndptrs(1)),i=1,3)
            print*,'   node ',ndptrs(numede+1),
     :                          (nodlst(i,ndptrs(numede+1)),i=1,3)
            ierr = -500
            return
         END IF
      ELSE IF( GTNDFL(BIGLST,NODLST,NOD1,2) .AND.
     :         GTNDFL(BIGLST,NODLST,NOD2,2) ) THEN
c         print*,'+++ Warning: attempted to shrink geometry edge +++'
c         print*,iedg,nod1,nod2
         GOTO 500
      ELSE IF( GTNDFL(BIGLST,NODLST,NOD1,2) .AND.
     :          GTNDFL(BIGLST,NODLST,NOD2,4) ) THEN
c         print*,'+++ Warning: attempted to shrink geom-split nodes +++'
         GOTO 500
      ELSE IF( GTNDFL(BIGLST,NODLST,NOD1,4) .AND.
     :          GTNDFL(BIGLST,NODLST,NOD2,2) ) THEN
c         print*,'+++ Warning: attempted to shrink split-geom nodes +++'
         GOTO 500
      ELSE IF( GTNDFL(BIGLST,NODLST,NOD1,4) .AND.
     :          GTNDFL(BIGLST,NODLST,NOD2,4) ) THEN
c         print*,'+++ Warning: attempted to shrink split geom edge +++'
         GOTO 500
      ELSE IF( GTNDFL(BIGLST,NODLST,NOD1,2) .OR.
     :         GTNDFL(BIGLST,NODLST,NOD1,4) ) THEN
         nod1 = 0
c         print*,'--- Info: shrinking surface edge to node 1 ---'
      ELSE IF( GTNDFL(BIGLST,NODLST,NOD2,2) .OR.
     :         GTNDFL(BIGLST,NODLST,NOD2,4) ) THEN
         nod2 = 0
c         print*,'--- Info: shrinking surface edge to node 2 ---'
      END IF
C
C - find a list of all elements, edges and nodes around the two nodes
C - that will shrink to one - or the one node that will move to the other
C
 50   CONTINUE
C
C - need to make sure we remove any test nodes that were added...
C
      IF( NOD1 .EQ. 0 ) THEN
C
c        if( debug ) print*,'   Collapsing to a node ',inod,nodo1
        IF( DIDADD ) CALL RMVND( BIGLST, NODLST, INOD )
        INOD = NODO1
        DIDADD = .TRUE.
        RATIO  = 0.0
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** EDGTST: GOT ERROR FROM RMVND'
            RETURN
         END IF
C
      ELSE IF( NOD2 .EQ. 0 ) THEN
C
c        if( debug ) print*,'   Collapsing to a node ',inod,nodo2
        IF( DIDADD ) CALL RMVND( BIGLST, NODLST, INOD )
        INOD = NODO2
        DIDADD = .TRUE.
        RATIO  = 1.0
C
        IF( IERR .NE. 0 ) THEN
           PRINT*,'*** EDGTST: GOT ERROR FROM RMVND'
           RETURN
        END IF
C
      END IF
C
      IF( NUMED1+NUMED2 .GT. 0 ) GOTO 987
C
      NUMED1 = 0
      NUMEL1 = 0
C
C      NUMED2 = 0
C      NUMEL2 = 0
C
      IF( NOD1 .GT. 0 ) THEN
C
         CALL NDCNCT( BIGLST, NODLST, NODO1, IELM,
     :                ELPTS2, EDPTS2, NDPTS2, NUMED1, NUMEL1 )
C
         IF( IERR .NE. 0 ) THEN
            print*,'*** EDGTST: GOT ERROR FROM NDCNCT'
            print*,'   node,element: ',nodo1,ielm
            RETURN
         END IF
C
c         if( debug ) then
c            print*,'   Back from NDCNCT 1 ...',numed1,numel1
c            print*,'      els: ',(elpts2(i),i=1,numel1)
c            print*,'      eds: ',(edpts2(i),i=1,numed1)
c            print*,'      nds: ',(ndpts2(i),i=1,numed1)
c         end if
         DO I = 1, NUMEL1
            CALL SETFLG( BIGLST, NODLST, ELPTS2(I), 6, .TRUE. )
         END DO
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** EDGTST: GOT ERROR FROM SETFLG'
            RETURN
         END IF
C
         NUMEDG = 0
C
         DO I = 1, NUMED1
            IF( GETFLG(BIGLST,NODLST,EDPTS2(I),6) ) GOTO 90
            IF( EDPTS2(I) .NE. IEDG ) THEN
               NUMEDG = NUMEDG + 1
               NDPTS2(NUMEDG) = NDPTS2(I)
               EDPTS2(NUMEDG) = EDPTS2(I)
               CALL STNDFL( BIGLST, NODLST, NDPTS2(NUMEDG), 7, .TRUE. )
            END IF
         END DO
C
         NUMED1 = NUMEDG
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** EDGTST: GOT ERROR FROM STNDFL/GETFLG'
            RETURN
         END IF
C
      END IF
C
      NUMELS = NUMEL1
      NUMEDG = NUMED1
      NUMEL2 = 0
      NUMED2 = 0
C
      IF( NOD2 .GT. 0 ) THEN
C
         CALL NDCNCT( BIGLST, NODLST, NODO2, IELM, ELPTS2(NUMEL1+1),
     :                EDPTS2(NUMED1+1), NDPTS2(NUMED1+1),
     :                NUMED2, NUMEL2 )
C
         IF( IERR .NE. 0 ) THEN
            print*,'*** EDGTST: GOT ERROR FROM NDCNCT'
            print*,'   node,element: ',nodo1,ielm
            RETURN
         END IF
C
c         if( debug ) then
c            print*,'   Back from NDCNCT 2 ...',numed2,numel2
c            print*,'      els: ',(elpts2(i),i=1+numel1,numel1+numel2)
c            print*,'      eds: ',(edpts2(i),i=numed1+1,numed1+numed2)
c            print*,'      nds: ',(ndpts2(i),i=numed1+1,numed1+numed2)
c         end if
         DO I = NUMEL1+1, NUMEL1+NUMEL2
            IF( .NOT. GETFLG( BIGLST, NODLST, ELPTS2(I), 6 ) ) THEN
               NUMELS = NUMELS + 1
               ELPTS2(NUMELS) = ELPTS2(I)
c            ELSE
c               if( debug ) print*,'   Repeated element ',elpts2(i)
c               CALL SETFLG( BIGLST, NODLST, ELPTS2(I), 5, .TRUE. )
            END IF
         END DO
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** EDGTST: GOT ERROR FROM GETFLG'
            RETURN
         END IF
C
         NUMEL2 = NUMELS - NUMEL1
C
         DO I = NUMED1+1, NUMED1+NUMED2
            IF( GETFLG(BIGLST,NODLST,EDPTS2(I),6) ) GOTO 90
            IF( EDPTS2(I) .NE. IEDG ) THEN
               NUMEDG = NUMEDG + 1
               NDPTS2(NUMEDG) = NDPTS2(I)
               EDPTS2(NUMEDG) = EDPTS2(I)
            END IF
         END DO
C
         NUMED2 = NUMEDG - NUMED1
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** EDGTST: GOT ERROR FROM GETFLG'
            RETURN
         END IF
C
      END IF
C
 987  CONTINUE
C
      IF( .NOT. DIDADD ) THEN
C
         XX = NODLST(1,NODO2)*RATIO + NODLST(1,NODO1)*(1.0-RATIO)
         YY = NODLST(2,NODO2)*RATIO + NODLST(2,NODO1)*(1.0-RATIO)
         ZZ = NODLST(3,NODO2)*RATIO + NODLST(3,NODO1)*(1.0-RATIO)
C
         INEW = -1
         INELM = INT(NODLST(16,NODO1))
         CALL ADDND( BIGLST, NODLST,
     :               ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :               SZENLS, SZNELS, NNOD,   NELM,
     :               NODX,   NODY,   NODZ,   ORGMTX,
     :               XX, YY, ZZ, INOD, INEW, INELM, .TRUE. )
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** EDGTST: GOT ERROR FROM ADDND (2)'
            RETURN
         END IF
c
         if( abs(xx-badndx).lt.1e-5.and.abs(yy-badndy).lt.1e-5) then
            print*,'*** EDGTST: ABOUT TO ADD BAD NODE (2)'
            print*, xx, yy, zz, ratio
            print*,NODLST(1,NODO1),NODLST(2,NODO1),NODLST(3,NODO1)
            print*,NODLST(1,NODO2),NODLST(2,NODO2),NODLST(3,NODO2)
         end if
c
         didadd = .true.
c
         if( inelm .lt. 0 ) then
           print*,'---+++ EDGTST(3): Got node outside mesh +++---'
           inelm = -inelm
           if( isgeom ) then
             print*,'  A geometry edge: ',ratio
           else if( isinsd ) then
             print*,'  An internal edge: ',ratio
           else
             print*,'  A plane surface edge: ',ratio
           end if
           print*,'node 1: ',(nodlst(i,nodo1),i=1,3)
           print*,'node 2: ',(nodlst(i,nodo2),i=1,3)
         end if
C
      END IF
C
C
C - work out the current functional of all the connected edges and
C - elements, as well as the predicted functional when the changes have
C - been made
C
      F2 = 0.0
      F3 = 0.0
      F4 = 0.0
      F5 = 1E+32
      F6 = 0.0
      F7 = F1
      F8 = 0.0
      F9 = 1E+32
      FA = 0.0
      CNTLFT = 0
      CNT1 = 0
      CNT2 = 0
      CNT3 = 0
C
      IF( NOD1 .GT. 0 ) THEN
C
       IF( NOD2 .GT. 0 .AND. RATIO .EQ. 0.5 ) THEN
          F5 = 0.0
          F9 = 0.0
       END IF
C
       DO I = 1, NUMEL1
C
         IF( GETFLG(BIGLST,NODLST,ELPTS2(I),5) ) THEN
            F3 = 1E+31
            F5 = 1E+31
            GOTO 60
         END IF
C         CALL SETFLG( BIGLST, NODLST, ELPTS2(I), 6, .FALSE. )
         CALL ELNODS( BIGLST, NODLST, ELPTS2(I), NODS, NXT, LST,
     :               .FALSE. )
C
         IF( IERR .NE. 0 ) THEN
            print*,'*** EDGTST: GOT ERROR FROM ELNODS'
            print*,'   element: ',elpts2(i)
            RETURN
         END IF
C
         F0 = ELMFNC( BIGLST, NODLST, ELPTS2(I),
     :                NODS(1), NODS(2), NODS(3), NODS(4), VOL1 )
         F2 = F2 + F0
         F1 = MAX( F1, F0 )
C
C - include elements with edge in the count (but don't count second time later)
C
         CNT3 = CNT3 + 1
C
         IF( NODS(1) .EQ. NODO2 ) GOTO 60
         IF( NODS(2) .EQ. NODO2 ) GOTO 60
         IF( NODS(3) .EQ. NODO2 ) GOTO 60
         IF( NODS(4) .EQ. NODO2 ) GOTO 60
C
         IF( NODS(1) .EQ. NODO1 ) THEN
            J = 1
         ELSE IF( NODS(2) .EQ. NODO1 ) THEN
            J = 2
         ELSE IF( NODS(3) .EQ. NODO1 ) THEN
            J = 3
         ELSE IF( NODS(4) .EQ. NODO1 ) THEN
            J = 4
         ELSE
            STOP
         END IF
C
C         DO J = 1, 4
C           IF( NODS(J) .EQ. NODO2 ) THEN
C              GOTO 60
C           ELSE IF( NODS(J) .EQ. NODO1 ) THEN
C              NODS(J) = INOD
C           END IF
C         END DO
C
         IF( F3 .LT. 1E+30 ) THEN
C
            NODS(J) = INOD
C
            F0 = ELMFNC( BIGLST, NODLST, 0,
     :                   NODS(1), NODS(2), NODS(3), NODS(4), VOL2 )
            F4 = F4 + F0
            F3 = MAX( F3 , F0 )
            CNTLFT = CNTLFT + 1
C
            IF( IERR .NE. 0 ) THEN
               PRINT*,'*** EDGTST: GOT ERROR FROM ELMFNC'
               RETURN
            END IF
C
c         if( vol1*vol2 .le. 0.0 ) print*,'INVERTED ELEMENT ',elpts2(i)
c         if( debug ) print*,'   element ',elpts2(i),
c     :          vol1,vol2,fncorg,fncdif
c         IF( ABS(VOL2) .LT. ORGRDS ) GOTO 90
            IF( VOL2/VOL1 .LE. 1E-4 ) F3 = 1E+31
C
         END IF
C
         IF( F5 .LT. 1E+30 ) THEN
C
            NODS(J) = NODO2
C
            F0 = ELMFNC( BIGLST, NODLST, 0,
     :                   NODS(1), NODS(2), NODS(3), NODS(4), VOL2 )
            F6 = F6 + F0
            F5 = MAX( F5 , F0 )
            CNT1 = CNT1 + 1
C
            IF( IERR .NE. 0 ) THEN
               PRINT*,'*** EDGTST: GOT ERROR FROM ELMFNC'
               RETURN
            END IF
C
            IF( VOL2/VOL1 .LE. 1E-4 ) F5 = 1E+31
C
         END IF
C
 60      CONTINUE
C
       END DO
C
      END IF
C
      FB = F2
C
      IF( NOD2 .GT. 0 ) THEN
C
       DO I = NUMEL1+1, NUMEL1+NUMEL2
C
         IF( GETFLG(BIGLST,NODLST,ELPTS2(I),5) ) THEN
            F3 = 1E+31
            F9 = 1E+31
            GOTO 80
         END IF
C
         CALL ELNODS( BIGLST, NODLST, ELPTS2(I), NODS, NXT, LST,
     :               .FALSE. )
C
         IF( IERR .NE. 0 ) THEN
            print*,'EDGTST: GOT ERROR FROM ELNODS/GETFLG'
            print*,'   element: ',elpts2(i)
            RETURN
         END IF
C
         F0 = ELMFNC( BIGLST, NODLST, ELPTS2(I),
     :                NODS(1), NODS(2), NODS(3), NODS(4), VOL1 )
         F8 = F8 + F0
         F7 = MAX( F7, F0 )
C
         IF( NODS(1) .EQ. NODO1 ) GOTO 80
         IF( NODS(2) .EQ. NODO1 ) GOTO 80
         IF( NODS(3) .EQ. NODO1 ) GOTO 80
         IF( NODS(4) .EQ. NODO1 ) GOTO 80
C
C - the elements around the edge were not included when finding the old
C - functional sum, as well as when counting original elements
C
         FB = FB + F0
         CNT3 = CNT3 + 1
C
         IF( NODS(1) .EQ. NODO2 ) THEN
            J = 1
         ELSE IF( NODS(2) .EQ. NODO2 ) THEN
            J = 2
         ELSE IF( NODS(3) .EQ. NODO2 ) THEN
            J = 3
         ELSE IF( NODS(4) .EQ. NODO2 ) THEN
            J = 4
         ELSE
            STOP
         END IF
C         DO J = 1, 4
C           IF( NODS(J) .EQ. NODO1 ) THEN
C              GOTO 80
C           ELSE IF( NODS(J) .EQ. NODO2 ) THEN
C              NODS(J) = INOD
C              K = J
C           END IF
C         END DO
C
         IF( F3 .LT. 1E+30 ) THEN
C
            NODS(J) = INOD
C
            F0 = ELMFNC( BIGLST, NODLST, 0,
     :                   NODS(1), NODS(2), NODS(3), NODS(4), VOL2 )
            F4 = F4 + F0
            F3 = MAX( F3, F0 )
            CNTLFT = CNTLFT + 1
C
            IF( IERR .NE. 0 ) THEN
               PRINT*,'*** EDGTST: GOT ERROR FROM ELMFNC'
               RETURN
            END IF
C
c         if( vol1*vol2 .le. 0.0 ) print*,'INVERTED ELEMENT ',elpts2(i)
c         if( debug ) print*,'   element ',elpts2(i),
c     :          vol1,vol2,fncorg,fncdif
c         IF( ABS(VOL2) .LT. ORGRDS ) GOTO 90
            IF( VOL2/VOL1 .LE. 1E-4 ) F3 = 1E+31
C
          END IF
C
          IF( F9 .LT. 1E+30 ) THEN
C
            NODS(J) = NODO1
C
            F0 = ELMFNC( BIGLST, NODLST, 0,
     :                   NODS(1), NODS(2), NODS(3), NODS(4), VOL2 )
            FA = FA + F0
            F9 = MAX( F9, F0 )
            CNT2 = CNT2 + 1
C
            IF( IERR .NE. 0 ) THEN
               PRINT*,'*** EDGTST: GOT ERROR FROM ELMFNC'
               RETURN
            END IF
C
c         if( vol1*vol2 .le. 0.0 ) print*,'INVERTED ELEMENT ',elpts2(i)
c         if( debug ) print*,'   element ',elpts2(i),
c     :          vol1,vol2,fncorg,fncdif
c         IF( ABS(VOL2) .LT. ORGRDS ) GOTO 90
            IF( VOL2/VOL1 .LE. 1E-4 ) F9 = 1E+31
C
          ELSE IF( F5 .GE. 1E+30 .AND. F3 .GE. 1E+30 ) THEN
             GOTO 90
          END IF
C
 80       CONTINUE
C
       END DO
C
       IF( NOD1 .LE. 0 ) THEN
          FB = F8
          CNT3 = NUMEL2
       END IF
C
      END IF
C
      SHRFNC = F3 - MAX(F1,F7)
      SHRAVE = F4 / CNTLFT - FB / CNT3
C
      IF( NOD2 .GT. 0 .AND. NOD1 .GT. 0 ) THEN
C
         FA = FA/CNT2 - F8/NUMEL2
         F9 = F9 - F7
C
         IF( F9 .LT. FNCDIF ) THEN
            IWHAT = 2
            FNCDIF = F9
            FNCORG = F7
            BSTX = NODLST(1,NODO1)
            BSTY = NODLST(2,NODO1)
            BSTZ = NODLST(3,NODO1)
            BSTRAT = 0.0
            BSTINE = INT(NODLST(16,NODO1))
         END IF
C
         IF( FA .LT. AVEDIF .AND. F9 .LT. 0.0 ) THEN
            AWHAT = 2
            AVEDIF = FA
            AVEORG = F8/NUMEL2
            BSTXA = NODLST(1,NODO1)
            BSTYA = NODLST(2,NODO1)
            BSTZA = NODLST(3,NODO1)
            BSTRAA = 0.0
            BSTINA = INT(NODLST(16,NODO1))
         END IF
C
         F6 = F6/CNT2 - F2/NUMEL1
         F5 = F5 - F1
C
         IF( F5 .LT. FNCDIF ) THEN
            IWHAT = 2
            FNCDIF = F5
            FNCORG = F1
            BSTX = NODLST(1,NODO2)
            BSTY = NODLST(2,NODO2)
            BSTZ = NODLST(3,NODO2)
            BSTRAT = 1.0
            BSTINE = INT(NODLST(16,NODO2))
         END IF
C
         IF( F6 .LT. AVEDIF .AND. F5 .LT. 0.0 ) THEN
            AWHAT = 2
            AVEDIF = F6
            AVEORG = F2/NUMEL1
            BSTXA = NODLST(1,NODO2)
            BSTYA = NODLST(2,NODO2)
            BSTZA = NODLST(3,NODO2)
            BSTRAA = 1.0
            BSTINA = INT(NODLST(16,NODO2))
         END IF
C
      END IF
C
C - always shrink an edge that is far too small,
C - and prevent a large edge from being shrunk
C
c      IF( ORGLEN .LT. 0.2 ) THEN
c         print*,'Weighting collapse of edge ',iedg,orglen
c         if( shrfnc .gt. 0.0 .and. shrfnc .lt. f0*50 ) then
c            shrfnc = -f0/5
c            if( shrfnc .lt. fncdif ) then
c               iwhat = 2
c               fncdif = shrfnc
c               fncorg = f0
c            end if
c         else if( shrfnc*50 .lt. fncdif ) then
c            IWHAT = 2
c            FNCORG = F0
c            FNCDIF = SHRFNC*10
c            IF( FNCDIF .GT. 0.0 .AND. SHRFNC .LT. FNCORG*50 )
c     :                                        FNCDIF = -FNCORG/5
c         end if
c      ELSE IF( SHRFNC .LT. FNCDIF .AND. ORGLEN .LT. 2.0 ) THEN
C
      IF( SHRFNC .LT. FNCDIF ) THEN
         IWHAT  = 2
         FNCDIF = SHRFNC
         FNCORG = MAX(F1,F7)
         BSTX   = XX
         BSTY   = YY
         BSTZ   = ZZ
         BSTRAT = RATIO
         BSTINE = INELM
      END IF
C
      IF( SHRAVE .LT. AVEDIF .AND. SHRFNC .LT. 0.0 ) THEN
         AWHAT  = 2
         AVEDIF = SHRAVE
         AVEORG = MAX(F2,F8)
         BSTXA  = XX
         BSTYA  = YY
         BSTZA  = ZZ
         BSTRAA = RATIO
         BSTINA = INELM
      END IF
c
c      if( debug ) then
c        print*,'   shrink edge functionals: ',F1,F3
c        print*,'   bits: ',F4,F3,F2,F1,F0
c      end if
C
  90  DO I = 1, NUMED1
         CALL STNDFL( BIGLST, NODLST, NDPTS2(I), 7, .FALSE. )
      END DO
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT ERROR FROM STNDFL'
         RETURN
      END IF
C
      DO I = 1, NUMEL1
         CALL SETFLG( BIGLST, NODLST, ELPTS2(I), 6, .FALSE. )
      END DO
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** EDGTST: GOT ERROR FROM SETFLG'
         RETURN
      END IF
C
 500  CONTINUE
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMSHR = TIMSHR + TIMEND - TIMSTT
      TOTSHR = TOTSHR + 1
C
  501 CONTINUE
C
C - remove the node that was added in the middle of the edge
C - (if it wasn't removed when checking edge shrinking).
C
      IF( NOD1 .NE. 0 .AND. NOD2 .NE. 0 .and. didadd ) THEN
C
         CALL RMVND( BIGLST, NODLST, INOD )
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** EDGTST: GOT ERROR FROM RMVND'
            RETURN
         END IF
C
      END IF
c
c      debug = .false.
C
      RETURN
      END
C
