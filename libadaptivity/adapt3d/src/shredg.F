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
      SUBROUTINE SHREDG( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   IEDG, ELPTRS, EDPTRS, NDPTRS, NUMEDE,
     :                   ORGRAT, X, Y, Z, INELM, SURFAC, INTRNL )
C-----------------------------------------------------------------------
C
C - This subroutine shrinks (collapses) the edge supplied (IEDG), and
C - also performs the necessary removal and joining of affected elements.
C - The necessary information in the arrays ELPTRS, EDPTRS, NDPTRS, as
C - well as the number of surrounding elements (NUMEDE) and the surface
C - value (SURFAC, if there is one) can be filled in easily by calling
C - subroutine EDGELS for the edge IEDG before calling this subroutine.
C
C-----------------------------------------------------------------------
      use write_log
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
      INTEGER IEDG, NUMEDE, ELPTRS(3,NUMEDE), EDPTRS(3,NUMEDE+1),
     :        NDPTRS(NUMEDE+1), SURFAC, INTRNL, INELM
C
      REAL ORGRAT, X, Y, Z
C
      INTEGER EDO1, EDO2, ED1, ED2, NODS(4), EDGLST(6), CONELM(4),
     :        IEDG1, IEDG2, IEDG3, IEDG4, IEDG5, LSTEL1, LSTEL2, INOD,
     :        NODS1(4), NDO, INEW, IFLAG, NXT, LST, NXTEL1, NXTEL2, I,
     :        IPOS, JPOS, IELM, J, K, NUMEDG, NUMELM, edpts2(maxede),
     :        elpts2(maxede), ndpts2(maxede), NUMED2, NUMEL2, IREG,
     :        nodt(4), IFAC, IED1, IED2, IED3
C
      REAL RATIO, elmfnc, v1, v2, q, edgfnc, rad
C
      LOGICAL GETFLG, ISGEOM, ISINSD, GTNDFL, ADDIT, goterr
C
      CALL ADAPT_SECOND( TIMSTT )
C
      numedg = numede
      if( surfac .gt. 0 ) numedg = numede+1
C
      CALL EDGINF( BIGLST, NODLST, IEDG, NODS, INEW, IFLAG, NXT, LST )
c
      goterr = debug
c      goterr = .true.
c      if( iedg .eq. 36128 ) goterr = .true.
c      if( iedg .eq. 944 ) goterr = .true.
c      if( iedg .eq. 26682 ) goterr = .true.
c      if( iedg .eq. 14519 ) goterr = .true.
c      if( iedg .eq. 26692 ) goterr = .true.
c      if( iedg .eq. 14425 ) goterr = .true.
c      if( iedg .eq. 55596 ) goterr = .true.
c      do i = 1, numede
c         if( elptrs(1,i) .eq. 30935 ) goterr = .true.
c         if( elptrs(2,i) .eq. 30935 ) goterr = .true.
c         if( elptrs(3,i) .eq. 30935 ) goterr = .true.
c         if( elptrs(1,i) .eq. 37104 ) goterr = .true.
c         if( elptrs(2,i) .eq. 37104 ) goterr = .true.
c         if( elptrs(3,i) .eq. 37104 ) goterr = .true.
c      end do
c      if( edptrs(2,numedg) .eq. 952 ) goterr = .true.
c      if( edptrs(1,numedg) .eq. 952 ) goterr = .true.
c      if( elptrs(2,numedg) .eq. 34117 ) goterr = .true.
c      if( elptrs(3,numedg) .eq. 34117 ) goterr = .true.
c
c
      if( goterr ) then
        ewrite(3,*) ' '
c        x = edgfnc( BIGLST, NODLST, nods(1), nods(2), q )
        ewrite(3,*)
     :       'Into SHREDG... ',ierr
        ewrite(3,*) 
     :       '   IEDG,NUMEDE,SURFAC: ',IEDG,NUMEDE,SURFAC
        ewrite(3,*) 
     :       '         X,Y,Z: ',X,Y,Z
        ewrite(3,*) 
     :       '   elptrs: ',(elptrs(1,i),i=1,numede)
        ewrite(3,*) 
     :       '   edptrs: ',(edptrs(3,i),i=1,numede)
        ewrite(3,*) 
     :       '   ndptrs: ',(ndptrs(i)  ,i=1,numedg)
        ewrite(3,*) 
     :       '   up els: ',(elptrs(2,i),i=1,numede)
        ewrite(3,*) 
     :       '   lw els: ',(elptrs(3,i),i=1,numede)
        ewrite(3,*) 
     :       '   up eds: ',(edptrs(1,i),i=1,numedg)
        ewrite(3,*) 
     :       '   lw eds: ',(edptrs(2,i),i=1,numedg)
      end if
  1   format( a, 30i7 )
C
      IF( IERR .NE. 0 ) RETURN
C
      RATIO = ORGRAT
C
      NXTEL1 = 0
      NXTEL2 = 0
C
      IF( SURFAC .GT. 0  ) THEN
         EDO1 = EDPTRS(1,NUMEDE+1)
         EDO2 = EDPTRS(2,NUMEDE+1)
         NDO  = NDPTRS(NUMEDE+1)
         LSTEL1 = SURFAC
         LSTEL2 = SURFAC
      ELSE
         EDO1 = EDPTRS(1,NUMEDE)
         EDO2 = EDPTRS(2,NUMEDE)
         NDO  = NDPTRS(NUMEDE)
         LSTEL1 = 0
         LSTEL2 = 0
      END IF
c
      if( goterr ) print*,'   edg nodes: ',nods(1),nods(2)
C
      ISGEOM = GETFLG( BIGLST, NODLST, IEDG, 2 )
      if( ierr .ne. 0 ) then
         print*,'*** SHREDG: Got an error from GETFLG (1)'
         return
      end if
      ISINSD = GETFLG( BIGLST, NODLST, IEDG, 3 ) .AND. .NOT. ISGEOM
      if( ierr .ne. 0 ) then
         print*,'*** SHREDG: Got an error from GETFLG (2)'
         return
      end if
C
      if( goterr ) then
        if( isgeom ) then
          print*,'   edge is geometry'
        else if( isinsd ) then
          print*,'   edge is inside'
        else
          print*,'   edge is on surface'
        end if
      end if
C
      IF( SURFAC .GT. 0 .AND. ISINSD ) THEN
         PRINT*,'*** SHREDG: SURFACE EDGE INCONSISTENCY'
         IERR = -111
         RETURN
      END IF
C
      IF( ISINSD ) THEN
         IF( .NOT. GTNDFL(BIGLST,NODLST,NODS(1),3) .AND.
     :       .NOT. GTNDFL(BIGLST,NODLST,NODS(2),3) ) THEN
            print*,'*** SHREDG: attempted to shrink internal edge +++'
            print*,iedg,nods(1),nods(2)
            ierr = -500
            RETURN
         ELSE IF( .NOT. GTNDFL(BIGLST,NODLST,NODS(1),3) ) THEN
            RATIO = 0.0
c            print*,'--- Info: shrinking internal edge to node 1 ---'
         ELSE IF( .NOT. GTNDFL(BIGLST,NODLST,NODS(2),3) ) THEN
            RATIO = 1.0
c            print*,'--- Info: shrinking internal edge to node 2 ---'
         END IF
      ELSE IF( ISGEOM ) THEN
         IF( GTNDFL(BIGLST,NODLST,NODS(1),2) .AND.
     :       GTNDFL(BIGLST,NODLST,NODS(2),2) ) THEN
            print*,'*** SHREDG: attempted to shrink geometry edge +++'
            print*,iedg,nods(1),nods(2)
            ierr = -500
            RETURN
         ELSE IF( GTNDFL(BIGLST,NODLST,NODS(1),2) .AND.
     :            GTNDFL(BIGLST,NODLST,NODS(2),4) ) THEN
            RATIO = 0.0
c            print*,'--- Info: shrinking geometry edge to node 1 ---'
         ELSE IF( GTNDFL(BIGLST,NODLST,NODS(2),2) .AND.
     :            GTNDFL(BIGLST,NODLST,NODS(1),4) ) THEN
            RATIO = 1.0
c            print*,'--- Info: shrinking geometry edge to node 2 ---'
         ELSE IF( GTNDFL(BIGLST,NODLST,NODS(1),4) .AND.
     :            GTNDFL(BIGLST,NODLST,NODS(2),4) ) THEN
c            print*,'--- INFO: shrinking a pre-split geometry edge ---'
         ELSE
            print*,'*** SHREDG: INCONSISTENT NODE FLAGS ON GEOM EDGE'
            print*,iedg,nods(1),nods(2)
            ierr = -500
            return
         END IF
      ELSE IF( GTNDFL(BIGLST,NODLST,NODS(1),2) .AND.
     :         GTNDFL(BIGLST,NODLST,NODS(2),2) ) THEN
         print*,'*** SHREDG: attempted to shrink geometry edge +++'
         print*,iedg,nods(1),nods(2)
         ierr = -500
         RETURN
      ELSE IF( GTNDFL(BIGLST,NODLST,NODS(1),4) .AND.
     :         GTNDFL(BIGLST,NODLST,NODS(2),4) ) THEN
         print*,'*** SHREDG: attempted to shrink split-split nodes +++'
         print*,iedg,nods(1),nods(2)
         ierr = -500
         RETURN
      ELSE IF( GTNDFL(BIGLST,NODLST,NODS(1),2) .AND.
     :         GTNDFL(BIGLST,NODLST,NODS(2),4) ) THEN
         print*,'*** SHREDG: attempted to shrink geom-split nodes +++'
         print*,iedg,nods(1),nods(2)
         ierr = -500
         RETURN
      ELSE IF( GTNDFL(BIGLST,NODLST,NODS(1),4) .AND.
     :         GTNDFL(BIGLST,NODLST,NODS(2),2) ) THEN
         print*,'*** SHREDG: attempted to shrink split-geom nodes +++'
         print*,iedg,nods(1),nods(2)
         ierr = -500
         RETURN
      ELSE IF( GTNDFL(BIGLST,NODLST,NODS(1),2) .OR.
     :         GTNDFL(BIGLST,NODLST,NODS(1),4) ) THEN
         RATIO = 0.0
c         print*,'--- Info: shrinking surface edge to node 1 ---'
      ELSE IF( GTNDFL(BIGLST,NODLST,NODS(2),2) .OR.
     :         GTNDFL(BIGLST,NODLST,NODS(2),4) ) THEN
         RATIO = 1.0
c         print*,'--- Info: shrinking surface edge to node 2 ---'
      END IF
C
c      X = NODLST(1,NODS(1))*(1.0-RATIO) + NODLST(1,NODS(2))*RATIO
c      Y = NODLST(2,NODS(1))*(1.0-RATIO) + NODLST(2,NODS(2))*RATIO
c      Z = NODLST(3,NODS(1))*(1.0-RATIO) + NODLST(3,NODS(2))*RATIO
C
      IF( RATIO .EQ. 1.0 ) THEN
         INOD = NODS(2)
         IPOS = 2
         JPOS = 1
         ED1  = EDO2
         ED2  = EDO1
         X = NODLST(1,NODS(2))
         Y = NODLST(2,NODS(2))
         Z = NODLST(3,NODS(2))
      ELSE IF( RATIO .EQ. 0.0 ) THEN
         INOD = NODS(1)
         IPOS = 1
         JPOS = 2
         ED1  = EDO1
         ED2  = EDO2
         X = NODLST(1,NODS(1))
         Y = NODLST(2,NODS(1))
         Z = NODLST(3,NODS(1))
      ELSE
         INEW = -1
c         INELM = INT(NODLST(16,NODS(1)))
         CALL ADDND( BIGLST, NODLST,
     :               ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :               SZENLS, SZNELS, NNOD,   NELM,
     :               NODX,   NODY,   NODZ,   ORGMTX,
     :               X, Y, Z, INOD, INEW, INELM, .FALSE. )
C
         IF( IERR .NE. 0 ) RETURN
c
         if( inelm .lt. 0 ) then
            print*,'---+++ SHREDG: Got node outside element +++---'
            inelm = -inelm
            if( isgeom ) then
               print*,'  A geometry edge: ',ratio
            else if( isinsd ) then
               print*,'  An internal edge: ',ratio
            else
               print*,'  A plane surface edge: ',ratio
            end if
            print*,'node 1: ',(nodlst(i,nods(1)),i=1,3)
            print*,'node 2: ',(nodlst(i,nods(2)),i=1,3)
         end if
c
         IF( ISGEOM ) THEN
            CALL STNDFL( BIGLST, NODLST, INOD, 4, .TRUE. )
         ELSE IF( ISINSD ) THEN
            CALL STNDFL( BIGLST, NODLST, INOD, 3, .TRUE. )
         END IF
         IPOS = 2
         JPOS = 1
         ED1  = EDO2
         ED2  = EDO1
      END IF
C
      IF( GTNDFL(BIGLST,NODLST,NODS(JPOS),2) ) THEN
         print*,'*** SHREDG: attempted to remove geometry node +++'
         print*,iedg,nods(1),nods(2)
         ierr = -500
         RETURN
      ELSE IF( GTNDFL(BIGLST,NODLST,NODS(JPOS),4) .AND.
     :         .NOT. ISGEOM ) THEN
         print*,'*** SHREDG: attempted to remove geom split node +++'
         print*,iedg,nods(1),nods(2)
         ierr = -500
         RETURN
      END IF
c
      do i = 1, numede
         if( elptrs(jpos+1,i) .lt. 0 ) then
             print*,'*** INCONSISTENCY IN SHREDG: ',IEDG,isgeom,isinsd
             print*,' Usually mean an element is very nearly planar...'
             print*,'edge being shrunk: ',iedg,'  has nodes:'
             print*,' nod1: ',nods(1),gtndfl(BIGLST,NODLST,nods(1),2),
     :                       gtndfl(BIGLST,NODLST,nods(1),3),
     :                       gtndfl(BIGLST,NODLST,nods(1),4)
             print*,'  pos: ',(nodlst(j,nods(1)),j=1,3)
             print*,' nod2: ',nods(2),gtndfl(BIGLST,NODLST,nods(2),2),
     :                       gtndfl(BIGLST,NODLST,nods(2),3),
     :                       gtndfl(BIGLST,NODLST,nods(2),4)
             print*,'  pos: ',(nodlst(j,nods(2)),j=1,3)
             print*,'stopped on ',i,ipos,jpos,ratio
             PRINT*,'  elptrs: ',(ELPTRS(1,J),J=1,NUMEDE)
             PRINT*,'  up els: ',(ELPTRS(2,J),J=1,NUMEDE)
             PRINT*,'  lw els: ',(ELPTRS(3,J),J=1,NUMEDE)
             PRINT*,'  up eds: ',(EDPTRS(1,J),J=1,NUMEDG)
             PRINT*,'  lw eds: ',(EDPTRS(2,J),J=1,NUMEDG)
             PRINT*,'  op eds: ',(EDPTRS(3,J),J=1,NUMEDE)
             PRINT*,'  ndptrs: ',(NDPTRS(J),J=1,NUMEDG)
             print*,'  problem element: ',elptrs(1,i),
     :              ' has face on surf: ',elptrs(jpos+1,i)
             print*,'  conelm: ',(biglst(j,elptrs(1,i)),j=1,3),
     :                            biglst(4,elptrs(1,i))/2
             ipos = biglst(nxtbig,elptrs(1,i))
             jpos = biglst(nxtbig,ipos)
             print*,'  edglst: ',(biglst(j,ipos),j=1,4),
     :                           (biglst(j,jpos),j=1,2)
             IERR = -503
             RETURN
          END IF
       END DO
C
cc      IF( RATIO .EQ. 0.5 .AND. ISINSD ) THEN
cc         DO I = 1, NUMEDE
cc            IF( .NOT. GETFLG(EDPTRS(1,I),3) ) THEN
cc               PRINT*,'*** SHREDG: INCONSISTENCY IN UPPER EDGES'
cc               print*,edptrs(1,i)
cc               IERR = -501
cc               RETURN
cc            END IF
cc            IF( .NOT. GETFLG(EDPTRS(2,I),3) ) THEN
cc               PRINT*,'*** SHREDG: INCONSISTENCY IN LOWER EDGES'
cc               print*,'Edge: ',edptrs(2,i)
cc               print*,'Nodes: ',(biglst(j,edptrs(2,i)),j=1,2)
cc               print*,'Element: ',elptrs(3,i)
cc               call elminf( elptrs(3,i), conelm, edglst, iflag, inew,
cc     :                      nxt, lst )
cc               print*,'   conelm: ',conelm
cc               do j = 1, 6
cc                 print*,'   edge: ',edglst(j),biglst(1,edglst(j)),
cc     :    biglst(2,edglst(j)),getflg(edglst(j),3),getflg(edglst(j),2)
cc               end do
cc               call elnods( elptrs(3,i), nods, nxt, lst )
cc               print*,'   nodes: ',nods
cc               print*,'   nod1: ',(nodlst(j,nods(1)),j=1,3),
cc     :                            (gtndfl(nods(1),j),j=2,4)
cc               print*,'   nod2: ',(nodlst(j,nods(2)),j=1,3),
cc     :                            (gtndfl(nods(2),j),j=2,4)
cc               print*,'   nod3: ',(nodlst(j,nods(3)),j=1,3),
cc     :                            (gtndfl(nods(3),j),j=2,4)
cc               print*,'   nod4: ',(nodlst(j,nods(4)),j=1,3),
cc     :                            (gtndfl(nods(4),j),j=2,4)
cc               IERR = -502
cc               RETURN
cc            END IF
cc         END DO
cc      END IF
C
C - remove all elements around the shrinking edge
C
      DO I = 1, NUMEDE
c         if( goterr ) then
c            call elnods( BIGLST, NODLST, elptrs(1,i), nodt, nxt, lst )
c            print*,'   removing elm ',elptrs(1,i)
c     :         elmfnc( BIGLST, NODLST, elptrs(1,i),
c     :                 nodt(1),nodt(2),nodt(3),nodt(4),rad)
c         end if
         if( goterr ) print*,'   removing element ',elptrs(1,i)
         CALL RMVELE( BIGLST, NODLST, ELPTRS(1,I) )
c         if( goterr ) then
c            print*,'   conelm to 2478: ',(biglst(j,2478),j=1,4)
c         end if
      END DO
C
C - remove the edge itself
C
      if( goterr ) print*,'   removing edge: ',iedg,ratio
      CALL RMVEDG( BIGLST, NODLST, IEDG )
C
      K = 0
 100  IF( K .LT. NUMEDE ) THEN
         K = K + 1
         IF( ELPTRS(JPOS+1,K) .LE. 0 ) GOTO 100
      ELSE
         NUMEDG = 0
         NUMELM = 0
         if( goterr ) print*,'NDCNCT has not been called...'
         GOTO 110
      END IF
C
C - find all edges/elements/nodes connected to node that will be removed
C
      CALL NDCNCT( BIGLST, NODLST, NODS(JPOS), ELPTRS(JPOS+1,K),
     :             ELPTS2, EDPTS2, NDPTS2, NUMEDG, NUMELM )
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,*) '*** SHREDG: Got an error from NDCNCT (1)'
         RETURN
      END IF
c
      if( goterr ) then
        ewrite(3,*) 'Back from NDCNCT:',numedg
        ewrite(3,*) '   els: ', elpts2(1:numelm)
        ewrite(3,*) '   eds: ', edpts2(1:numedg)
        ewrite(3,*) '   nds: ', ndpts2(1:numedg)
      end if
C
c      inew = -1
c      call addnd( x, y, z, k, inew )
c
c      do i = 1, numelm
c         call elnods( BIGLST, NODLST, elpts2(i), nods1, nxt, lst )
c         q = elmfnc( BIGLST, NODLST, nods1(1), nods1(2), nods1(3), nods1(4), v1 )
c         if( nods1(1) .eq. nods(jpos) ) then
c            nods1(1) = k
c         else if( nods1(2) .eq. nods(jpos) ) then
c            nods1(2) = k
c         else if( nods1(3) .eq. nods(jpos) ) then
c            nods1(3) = k
c         else if( nods1(4) .eq. nods(jpos) ) then
c            nods1(4) = k
c         end if
c         q = elmfnc( BIGLST, NODLST, nods1(1), nods1(2), nods1(3), nods1(4), v2 )
c         if( goterr ) then
c            print*,'   element ',elpts2(i),v1,v2
c         end if
c         if( v2/v1 .le. 1e-5 ) then
c            print*,'   ELEMENT ',elpts2(i),v1,v2,ratio
c            ierr = -1
c            return
c         end if
c      end do
c
c      do i = 1, numede
c         if( elptrs(ipos+1,i) .gt. 0 ) then
c          call elnods( elptrs(ipos+1,i), nods1, nxt, lst )
c          q = elmfnc( nods1(1), nods1(2), nods1(3), nods1(4), v1 )
c          if( nods1(1) .eq. nods(jpos) ) then
c            nods1(1) = k
c          else if( nods1(2) .eq. nods(jpos) ) then
c            nods1(2) = k
c          else if( nods1(3) .eq. nods(jpos) ) then
c            nods1(3) = k
c          else if( nods1(4) .eq. nods(jpos) ) then
c            nods1(4) = k
c          end if
c          q = elmfnc( BIGLST, NODLST, nods1(1), nods1(2), nods1(3), nods1(4), v2 )
c          if( goterr ) then
c            print*,'   element ',elptrs(ipos+1,i),v1,v2
c          end if
c          if( v2/v1 .le. 1e-5 ) then
c            print*,'   ELEMENT ',elptrs(ipos+1,i),v1,v2,ratio
c            ierr = -1
c            return
c          end if
c         end if
c      end do
c
c      call rmvnd( BIGLST, NODLST, k )
c
C - replace node NODS(JPOS) with 'new' node INOD in all the connected edges
C
      DO I = 1, NUMEDG
         IEDG1 = EDPTS2(I)
         IF( BIGLST(1,IEDG1) .EQ. NODS(JPOS) ) THEN
            BIGLST(1,IEDG1) = INOD
         ELSE IF( BIGLST(2,IEDG1) .EQ. NODS(JPOS) ) THEN
            BIGLST(2,IEDG1) = INOD
         ELSE
            PRINT*,'*** SHREDG: EDGE DOES NOT HAVE NODE'
            PRINT*,'   edge,node: ',iedg1,nods(jpos)
            IERR = -211
            RETURN
         END IF
      END DO
C
C - now that nothing is left using this node we can remove it
C
 110  CONTINUE
      if( goterr ) print*,'   removing node ',nods(jpos),jpos,inod
      CALL RMVND( BIGLST, NODLST, NODS(JPOS) )
C
      NUMED2 = 0
      NUMEL2 = 0
C
      K = 0
 120  IF( K .LT. NUMEDE ) THEN
         K = K + 1
         IF( ELPTRS(IPOS+1,K) .LE. 0 ) GOTO 120
      ELSE
         if( goterr ) print*,'NDCNCT has not been called...'
         GOTO 130
      END IF
C
C - if the shrink is to either end of the removed edge, then we don't
C - need to bother with the next chunk...
C
      IF( RATIO .LT. 1.0 .AND. RATIO .GT. 0.0 ) THEN
C
C - find all edges/elements/nodes connected to NODS(IPOS)
C
         CALL NDCNCT( BIGLST, NODLST, NODS(IPOS), ELPTRS(IPOS+1,K),
     :                ELPTS2(NUMELM+1), EDPTS2(NUMEDG+1),
     :                NDPTS2(NUMEDG+1), NUMED2, NUMEL2 )
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** SHREDG: Got an error from NDCNCT (2)'
            RETURN
         END IF
c
         if( goterr ) then
           ewrite(3,*) 
     :           'Back from NDCNCT 2:',numed2,numel2
           ewrite(3,*) 
     :          '   els: ', elpts2(numelm+1:numelm+numel2)
           ewrite(3,*) 
     :          '   eds: ', edpts2(1+numedg:numedg+numed2)
           ewrite(3,*) 
     :          '   nds: ', ndpts2(1+numedg:numedg+numed2)
         end if
C
C - replace NODS(IPOS) with INOD in all connected edges
C
         DO I = 1+NUMEDG, NUMEDG+NUMED2
            IEDG1 = EDPTS2(I)
            IF( BIGLST(1,IEDG1) .EQ. NODS(IPOS) ) THEN
               BIGLST(1,IEDG1) = INOD
            ELSE IF( BIGLST(2,IEDG1) .EQ. NODS(IPOS) ) THEN
               BIGLST(2,IEDG1) = INOD
            ELSE
               PRINT*,'*** SHREDG: EDGE DOES NOT HAVE NODE'
               PRINT*,'   edge,node: ',iedg1,nods(ipos)
               IERR = -211
               RETURN
            END IF
         END DO
C
C - now we can get rid of the node
C
         CALL RMVND( BIGLST, NODLST, NODS(IPOS) )
C
      END IF
C
 130  IF( IERR .NE. 0 ) RETURN
C
C - 'refresh' all elements that were connected to NODS(IPOS)
C
      DO I = 1, NUMEL2
         K = ELPTS2(NUMELM+I)
         IF( K .GT. 0 ) THEN
           CALL ELMINF( BIGLST, NODLST, K, CONELM, EDGLST,
     :                  INEW, IFLAG, IREG, NXT, LST )
           IF( IERR .NE. 0 ) THEN
              PRINT*,'*** SHREDG: Got an error from ELMINF (1)'
              RETURN
           END IF
           if( goterr ) then
c              call elnods( BIGLST, NODLST, K, nodt, nxt, lst )
              print*,'   removing elm ',K
c     :          elmfnc( BIGLST, NODLST, k,
c     :                  nodt(1),nodt(2),nodt(3),nodt(4),rad)
           end if
           CALL RMVELE( BIGLST, NODLST, K )
           INEW = -1
           CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, K, INEW, IREG )
      if( ierr.eq.12345 ) then
         print*,'--- inconsistent oppedg from addele in shredg (1)'
         ierr = 0
      end if
           IF( IERR .NE. 0 ) THEN
              PRINT*,'*** SHREDG: Got an error from RMVELE/ADDELE (1)'
              RETURN
           END IF
           if( goterr ) then
c              call elnods( BIGLST, NODLST, K, nodt, nxt, lst )
              print*,'   added elm ',K
c     :          elmfnc( BIGLST, NODLST, k,
c     :                  nodt(1),nodt(2),nodt(3),nodt(4),rad)
c              biglst(7,K) = 0
c              nxt = biglst(nxtbig,K)
c              biglst(7,nxt) = 0
           end if
         END IF
      END DO
C
C      NODLST(1,INOD) = X
C      NODLST(2,INOD) = Y
C      NODLST(3,INOD) = Z
C
C - reconnect the 'upper' and 'lower' elements, and replace edges
C - that have vanished
C
      IF( SURFAC .GT. 0 ) THEN
         NUMEDG = NUMEDE + 1
      ELSE
         NUMEDG = NUMEDE
      END IF
C
      ED1 = EDPTRS(IPOS,NUMEDG)
      ED2 = EDPTRS(JPOS,NUMEDG)
C
c      print*,'low eds: ',(edptrs(ipos,i),i=1,numedg)
c      print*,'upp eds: ',(edptrs(jpos,i),i=1,numedg)
c      print*,'out eds: ',(edptrs(3,i),i=1,numede)
C
C - Swap out the removed edges from elements that had them
C
      DO I = 1, NUMELM
C
         IELM = ELPTS2(I)
C
c         IF( .NOT. GETFLG(BIGLST,NODLST,IELM,4) ) THEN
C
            NXT = BIGLST(NXTBIG,IELM)
            LST = BIGLST(NXTBIG,NXT)
C
            DO K = 1, NUMEDG
C
               DO J = 1, 4
                  IF( BIGLST(J,NXT) .EQ. EDPTRS(JPOS,K) ) THEN
c                     print*,ielm,conelm
c                     print*,edglst
c                     print*,edptrs(jpos,k),edptrs(ipos,k)
                     BIGLST(J,NXT) = EDPTRS(IPOS,K)
c                     if( getflg(BIGLST,NODLST,ielm,4) ) print*,'   Missed edge: ',
c     :                  edptrs(jpos,k),edglst(j),ielm,k,j
                  END IF
               END DO
C
               DO J = 1, 2
                  IF( BIGLST(J,LST) .EQ. EDPTRS(JPOS,K) ) THEN
c                     print*,ielm,conelm
c                     print*,edglst
c                     print*,edptrs(jpos,k),edptrs(ipos,k)
                     BIGLST(J,LST) = EDPTRS(IPOS,K)
c                     if( getflg(BIGLST,NODLST,ielm,4) ) print*,'   Missed edge: ',
c     :                  edptrs(jpos,k),edglst(j),ielm,k,j
                  END IF
               END DO
C
            END DO
C
c         else
c            print*,'   ignored element ',ielm,i
c         END IF
C
      END DO
C
C - replace the elements (just to move them to the end of the list)
C
      DO I = 1, NUMELM
C
         IELM = ELPTS2(I)
C
c         IF( .NOT. GETFLG(BIGLST,NODLST,IELM,4) ) THEN
C
            CALL ELMINF( BIGLST, NODLST, IELM, CONELM, EDGLST,
     :                   INEW, IFLAG, IREG, NXT, LST )
            IF( IERR .NE. 0 ) THEN
               PRINT*,'*** SHREDG: Got an error from ELMINF (2)'
               RETURN
            END IF
c            print*,'   checking element ',ielm,i
c            print*,'      conelm: ',conelm
c            print*,'      edglst: ',edglst
C
               if( goterr ) print*,'   removing element ',ielm
               CALL RMVELE( BIGLST, NODLST, IELM )
               INEW = -1
               CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST,
     :                      IELM, INEW, IREG )
      if( ierr.eq.12345 ) then
         print*,'--- inconsistent oppedg from addele in shredg (2)'
         ierr = 0
      end if
           IF( IERR .NE. 0 ) THEN
              PRINT*,'*** SHREDG: Got an error from RMVELE/ADDELE (2)'
              RETURN
           END IF
C
                ELPTS2(I) = IELM
c               if( goterr ) then
c                  call elnods( BIGLST, NODLST, ielm, nodt, nxt, lst )
c                  print*,'   added elm ',ielm,
c     :              elmfnc( BIGLST, NODLST, ielm,
c     :                      nodt(1),nodt(2),nodt(3),nodt(4),rad)
c                  biglst(7,IELM) = 0
c                  nxt = biglst(nxtbig,IELM)
c                  biglst(7,nxt) = 0
c               end if
               if( goterr ) print*,'   added element ',ielm
C
c         else
c            print*,'   ignored element ',ielm,i
c         END IF
C
      END DO
C
C
      DO I = 1, NUMEDE
C
         EDO1 = ED1
         EDO2 = ED2
C
         ED1 = EDPTRS(IPOS,I)
         ED2 = EDPTRS(JPOS,I)
C
         K = ELPTRS(IPOS+1,I)
C
C         IF( K .GT. 0 ) THEN
C
C           CALL ELMINF( BIGLST, NODLST, K, CONELM, EDGLST, INEW, IFLAG, NXT, LST )
C
C           CALL RMVELE( BIGLST, NODLST, K )
C
C           INEW = -1
C           CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, ELPTRS(IPOS+1,I), INEW )
C
C         END IF
C
         IELM = ELPTRS(JPOS+1,I)
C
ccc         CALL CHGEDS( BIGLST, NODLST, IELM, ED2, EDO2, ED1, EDO1 )
C
c         DO J = 1, 6
C
c            IF( EDGLST(J) .EQ. EDO2 ) THEN
c               EDGLST(J) = EDO1
c            ELSE IF( EDGLST(J) .EQ. ED2 ) THEN
c               EDGLST(J) = ED1
C            ELSE IF( BIGLST(1,EDGLST(J)) .EQ. NODS(JPOS) ) THEN
C               BIGLST(1,EDGLST(J)) = INOD
C            ELSE IF( BIGLST(2,EDGLST(J)) .EQ. NODS(JPOS) ) THEN
C               BIGLST(2,EDGLST(J)) = INOD
c            END IF
C
c         END DO
C
C         INEW = -1
C         CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM, INEW )
C
c         if( goterr ) then
c            print*,'   connecting elements: ',ielm,elptrs(ipos+1,i)
c            print*,'   conelm for ielm: ',(biglst(j,ielm),j=1,4)
c            print*,'   conelm for jelm: ',
c     :                        (biglst(j,elptrs(ipos+1,i)),j=1,4)
c         end if
c
         IFAC = 0
         IED1 = 0
c         print*,'--- connecting elms in SHREDG: ',ielm,k,
c     :                                            ed1,edo1,edptrs(1,i)
         CALL NWCNEL( BIGLST, NODLST, IELM, K, IFAC,
     :                ED1, EDO1, EDPTRS(3,I) )
C         CALL CNCTEL( BIGLST, NODLST, IELM, K )
         CALL SETFLG( BIGLST, NODLST, IELM, 3, .FALSE. )
         IF( K .GT. 0 ) CALL SETFLG( BIGLST, NODLST, K, 3, .FALSE. )
c
c         if( goterr ) then
c            print*,'   connected elements: ',ielm,elptrs(ipos+1,i)
c            print*,'   conelm for ielm: ',(biglst(j,ielm),j=1,4)
c            print*,'   conelm for jelm: ',
c     :                        (biglst(j,elptrs(ipos+1,i)),j=1,4)
c         end if
C
c         CALL SETFLG( BIGLST, NODLST, IELM, 4, .TRUE. )
C
         IF( IERR .NE. 0 ) THEN
            PRINT*,'*** SHREDG: Got an error from CNCTEL/SETFLG'
            print*,'  failed for i,ipos,jpos: ',i,ipos,jpos
            PRINT*,'  elptrs: ',(ELPTRS(1,J),J=1,NUMEDE)
            PRINT*,'  up els: ',(ELPTRS(2,J),J=1,NUMEDE)
            PRINT*,'  lw els: ',(ELPTRS(3,J),J=1,NUMEDE)
            PRINT*,'  up eds: ',(EDPTRS(1,J),J=1,NUMEDG)
            PRINT*,'  lw eds: ',(EDPTRS(2,J),J=1,NUMEDG)
            PRINT*,'  op eds: ',(EDPTRS(3,J),J=1,NUMEDE)
            PRINT*,'  ndptrs: ',(NDPTRS(J),J=1,NUMEDG)
            RETURN
         END IF
C
      END DO
C
C      ADDIT = .FALSE.
C
      DO I = 1, NUMEDE
         if( goterr )
     :         print*,'   removing edge ',edptrs(jpos,i)
         CALL RMVEDG( BIGLST, NODLST, EDPTRS(JPOS,I) )
c         CALL SETFLG( BIGLST, NODLST, ELPTRS(JPOS+1,I), 4, .FALSE. )
      END DO
C
      IF( SURFAC .GT. 0 ) CALL RMVEDG( BIGLST, NODLST,
     :                                 EDPTRS(JPOS,NUMEDE+1) )
C
c      DO I = 1, NUMEDG
C
c         IF( BIGLST(1,EDPTS2(I)) .EQ. NODS(JPOS) ) THEN
c            BIGLST(1,EDPTS2(I)) = INOD
c         ELSE IF( BIGLST(2,EDPTS2(I)) .EQ. NODS(JPOS) ) THEN
c            BIGLST(2,EDPTS2(I)) = INOD
c         END IF
C
c      END DO
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMCHG = TIMCHG + TIMEND - TIMSTT
      TOTCHG = TOTCHG + 1
c
      if( .not. cnsych ) return
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
                  PRINT*,'*** SHREDG: INCONSISTENT CONNECTED ELEMENTS'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,conelm(j)) .le. 0 ) then
               if( .not. addit ) then
                  print*,'*** SHREDG: INVALID CONNECTED ELEMENT'
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
         call elnods( BIGLST,NODLST, i, nods, nxt, lst, .false. )
         x = elmfnc( BIGLST,NODLST, 0,
     :               nods(1), nods(2), nods(3), nods(4), y )
         if( abs(y) .lt. 1e-8 ) then
            if( .not. addit ) then
               PRINT*,'*** SHREDG: SMALL ELEMENT'
               print*,'   element: ',i
               print*,'   conelm: ',conelm
               print*,'   edglst: ',edglst
               addit  = .true.
               goterr = .true.
            end if
            PRINT*,'   vol,fnc: ',y,x
         end if
         do j = 1, 6
            nods(1) = biglst(1,edglst(j))
            nods(2) = biglst(2,edglst(j))
            IF( edglst(J) .LE. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** SHREDG: INCONSISTENT EDGE LIST'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,edglst(j)) .ge. 0 ) then
               if( .not. addit ) then
                  PRINT*,'*** SHREDG: INVALID EDGE'
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
                    PRINT*,'*** SHREDG: INVALID GEOM-EDGE NODE 1 FLAGS'
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
                    PRINT*,'*** SHREDG: INVALID GEOM-EDGE NODE 2 FLAGS'
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
                    PRINT*,'*** SHREDG: INVALID SURF-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
               else if( gtndfl(BIGLST,NODLST,nods(2),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** SHREDG: INVALID SURF-EDGE NODE 2 FLAGS'
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
         print*,'*** SHREDG: INCONSISTENCY AT END OF BIGLST'
         print*,i,emtbig
         ierr = -10000
      end if
C
      RETURN
      END
C
