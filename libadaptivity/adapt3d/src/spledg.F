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
      SUBROUTINE SPLEDG( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   IEDG,   ELPTRS, EDPTRS, NDPTRS, NUMEDE,
     :                   X, Y, Z, INELM,  SURFAC, INTRNL,
     :                   NEWED1, NEWED2, NEWEL1, NEWEL2 )
C-----------------------------------------------------------------------
C
C - This subroutine splits the edge supplied (IEDG) and also performs
C - the necessary splitting of surrounding elements.
C - The necessary information in the arrays ELPTRS, EDPTRS, NDPTRS, as
C - well as the number of surrounding elements (NUMEDE) and the surface
C - value (SURFAC, if there is one) can be filled in easily by calling
C - subroutine EDGELS for the edge IEDG before calling this subroutine.
C
C - P.S. Gotta love the name!
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
     :        NDPTRS(NUMEDE+1), SURFAC, NEWED1, NEWED2, NEWEL1, NEWEL2,
     :        INTRNL, INELM
C
      REAL X, Y, Z
C
      INTEGER EDO1, EDO2, ED1, ED2, NODS(4), EDGLST(6), CONELM(4), J, K,
     :        IEDG1, IEDG2, IEDG3, IEDG4, IEDG5, LSTEL1, LSTEL2, INOD,
     :        NODS1(2), NDO, INEW, IFLAG, NXT, LST, NXTEL1, NXTEL2, I,
     :        ELREGS(NUMEDE), KK, nodt(4), SURFC1, SURFC2, WCHFAC
C
      REAL elmfnc, rad
C
      LOGICAL GETFLG, ISGEOM, ISINSD, GTNDFL, addit, goterr
C
      CALL ADAPT_SECOND( TIMSTT )
C
      goterr = DEBUG
c      if( intrnl .lt. 0 ) goterr = .true.
c      goterr = .true.
c      do i = 1, numede
c         if( edptrs(1,i) .eq. 952 ) goterr = .true.
c         if( edptrs(2,i) .eq. 952 ) goterr = .true.
c         if( edptrs(3,i) .eq. 952 ) goterr = .true.
c      end do
c      if( surfac .gt. 0 ) then
c        if( edptrs(2,numede+1) .eq. 952 ) goterr = .true.
c        if( edptrs(1,numede+1) .eq. 952 ) goterr = .true.
c      end if
c      do i = 1, numede
c         if( elptrs(1,i) .eq. 30935 ) goterr = .true.
c         if( elptrs(2,i) .eq. 30935 ) goterr = .true.
c         if( elptrs(3,i) .eq. 30935 ) goterr = .true.
c         if( elptrs(1,i) .eq. 37104 ) goterr = .true.
c         if( elptrs(2,i) .eq. 37104 ) goterr = .true.
c         if( elptrs(3,i) .eq. 37104 ) goterr = .true.
c      end do
c      if( iedg .eq. 7825 ) goterr = .true.
c      if( iedg .eq. 944 ) goterr = .true.
c      if( iedg .eq. 26682 ) goterr = .true.
c      if( iedg .eq. 14519 ) goterr = .true.
c      if( iedg .eq. 26692 ) goterr = .true.
c      if( iedg .eq. 30130 ) goterr = .true.
c      if( iedg .eq. 32454 ) goterr = .true.
c
      NEWED1 = 0
      NEWED2 = 0
C
      NXTEL1 = 0
      NXTEL2 = 0
C - just to get rid of compile warnings...
      NDO = 0
C
      if( ierr .eq. 0 ) then
c
      IF( SURFAC .GT. 0  ) THEN
C
         EDO1 = EDPTRS(1,NUMEDE+1)
         EDO2 = EDPTRS(2,NUMEDE+1)
         NDO  = NDPTRS(NUMEDE+1)
C
         surfc1 = wchfac( biglst, nodlst, elptrs(1,1),
     :                    edo1, edo2, iedg )
         surfc2 = wchfac( biglst, nodlst, elptrs(1,numede),
     :                    edptrs(1,numede), edptrs(2,numede), iedg )
C
         LSTEL1 = biglst(surfc1,elptrs(1,1))
         if( surfc1 .eq. 4 ) LSTEL1 = LSTEL1/2
         if(lstel1 .ge. 0) then
            print*,'*** SPLEDG: First face not a surface!',lstel1,surfc1
            lstel1 = -surfac
         endif
C
         LSTEL2 = biglst(surfc2,elptrs(1,numede))
         if( surfc2 .eq. 4 ) LSTEL2 = LSTEL2/2
         if(LSTEL2 .ge. 0) then
            print*,'*** SPLEDG: Last face not a surface!',lstel2,surfc2
            surfc2 = -surfac
         else
            surfc2 = lstel2
         endif
C
         LSTEL2 = LSTEL1
C
      ELSE
C
         EDO1 = EDPTRS(1,NUMEDE)
         EDO2 = EDPTRS(2,NUMEDE)
         NDO  = NDPTRS(NUMEDE)
         LSTEL1 = 0
         LSTEL2 = 0
C
      END IF
c
      end if
C
      if( goterr ) then
        PRINT*,' '
        PRINT*,'Into SPLEDG... ',ierr
        print*,'    IEDG,NUMEDE: ',IEDG,NUMEDE
        if(surfac.gt.0) then
           PRINT*,'   SURFAC,INTRNL: ',SURFAC,INTRNL,surfc1,surfc2
        else
           PRINT*,'   SURFAC,INTRNL: ',SURFAC,INTRNL
        end if
        print*,'   x,y,z: ',x,y,z
        print*,'   nodes: ',biglst(1,iedg),biglst(2,iedg)
        print*,'   upnd: ',(nodlst(i,biglst(1,iedg)),i=1,3)
        print*,'   lwnd: ',(nodlst(i,biglst(2,iedg)),i=1,3)
        ewrite(3,*)
     :       '   els: ',(elptrs(1,i),i=1,numede)
        ewrite(3,*)
     :       '   eds: ',(edptrs(3,i),i=1,numede)
        ewrite(3,*)
     :       '   nds: ',(ndptrs(i),i=1,numede+1)
        ewrite(3,*)
     :       '   upels: ',(elptrs(2,i),i=1,numede)
        ewrite(3,*)
     :       '   lwels: ',(elptrs(3,i),i=1,numede)
        ewrite(3,*)
     :       '   upeds: ',(edptrs(1,i),i=1,numede+1)
        ewrite(3,*)
     :       '   lweds: ',(edptrs(2,i),i=1,numede+1)
        do i = 1, numede+1
           print*,'   node',i,':',nodlst(1,ndptrs(i)),
     :            nodlst(2,ndptrs(i)),nodlst(3,ndptrs(i))
        end do
cc        pause
      end if
c
   1  format( a, 30i7 )
C
      IF( IERR .NE. 0 ) RETURN
C
      CALL EDGINF( BIGLST, NODLST, IEDG, NODS, INEW, IFLAG, NXT, LST )
C
      ISGEOM = GETFLG( BIGLST, NODLST, IEDG, 2 )
      ISINSD = GETFLG( BIGLST, NODLST, IEDG, 3 ) .AND. .NOT. ISGEOM
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
         PRINT*,'*** SPLEDG: SURFACE EDGE INCONSISTENCY'
         IERR = -111
         RETURN
      END IF
C
      DO I = 1, NUMEDE
         NXT = BIGLST(NXTBIG,ELPTRS(1,I))
         NXT = BIGLST(NXTBIG,NXT)
         ELREGS(I) = BIGLST(7,NXT)
         if( goterr ) then
c            call elnods( BIGLST, NODLST, elptrs(1,i), nodt, nxt, lst )
            print*,'   removing element ',elptrs(1,i)
c     :        elmfnc( BIGLST, NODLST, elptrs(1,i),
c     :                nodt(1),nodt(2),nodt(3),nodt(4),rad)
         end if
         CALL RMVELE( BIGLST, NODLST, ELPTRS(1,I) )
      END DO
C
      if( goterr ) print*,'   removing edge ',iedg
      CALL RMVEDG( BIGLST, NODLST, IEDG )
C
c      X = NODLST(1,NODS(1))*(1.0-RATIO) + NODLST(1,NODS(2))*RATIO
c      Y = NODLST(2,NODS(1))*(1.0-RATIO) + NODLST(2,NODS(2))*RATIO
c      Z = NODLST(3,NODS(1))*(1.0-RATIO) + NODLST(3,NODS(2))*RATIO
C
      if( goterr ) print*,'    regs: ',(elregs(i),i=1,numede)
      INEW = -1
C      INELM = INT(NODLST(16,NODS(1)))
      CALL ADDND( BIGLST, NODLST,
     :            ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :            SZENLS, SZNELS, NNOD,   NELM,
     :            NODX,   NODY,   NODZ,   ORGMTX,
     :            X, Y, Z, INOD, INEW, INELM, .FALSE. )
C
      IF( IERR .NE. 0 ) RETURN
C
      if( goterr ) print*,'   added node ',inod,' at ',x,y,z
c
      if( inelm .lt. 0 ) then
         print*,'---+++ SPLEDG: Got node outside element +++---'
         inelm = -inelm
         if( isgeom ) then
            print*,'  Geometry edge: ',x,y,z
         else if( isinsd ) then
            print*,'  Internal edge: ',x,y,z
         else
            print*,'  Plane surface edge: ',x,y,z
         end if
         print*,'node 1: ',(nodlst(i,nods(1)),i=1,3)
         print*,'node 2: ',(nodlst(i,nods(2)),i=1,3)
      end if
c
      NODS1(1) = NODS(1)
      NODS1(2) = INOD
C
      INEW = -1
      CALL ADDEDG( BIGLST, NODLST, NODS1, IEDG1, INEW )
      if( goterr ) print*,'   added edge ',iedg1,', nodes ',nods1
C
      NODS1(1) = NODS(2)
      NODS1(2) = INOD
C
      INEW = -1
      CALL ADDEDG( BIGLST, NODLST, NODS1, IEDG2, INEW )
      if( goterr ) print*,'   added edge ',iedg2,', nodes ',nods1
C
      IF( ISGEOM ) THEN
         CALL STNDFL( BIGLST, NODLST, INOD,  4, .TRUE. )
         CALL SETFLG( BIGLST, NODLST, IEDG1, 2, .TRUE. )
         CALL SETFLG( BIGLST, NODLST, IEDG2, 2, .TRUE. )
      ELSE IF( ISINSD ) THEN
         CALL STNDFL( BIGLST, NODLST, INOD,  3, .TRUE. )
         CALL SETFLG( BIGLST, NODLST, IEDG1, 3, .TRUE. )
         CALL SETFLG( BIGLST, NODLST, IEDG2, 3, .TRUE. )
      END IF
C
      NEWED1 = IEDG1
      NEWED2 = IEDG2
C
      NODS1(1) = INOD
      NODS1(2) = NDO
C
      INEW = -1
      CALL ADDEDG( BIGLST, NODLST, NODS1, IEDG5, INEW )
      if( goterr ) print*,'   added edge ',iedg5,', nodes ',nods1
C
      IF( ISINSD .OR. (INTRNL .LT. 0 .AND. SURFAC .LE. 0 ) ) THEN
         CALL SETFLG( BIGLST, NODLST, IEDG5, 3, .TRUE. )
C      ELSE IF( GTNDFL(BIGLST,NODLST,NDO,3) ) THEN
C         CALL SETFLG( BIGLST, NODLST, IEDG5, 3, .TRUE. )
      END IF
C
      IEDG4 = IEDG5
C
      IF( IERR .NE. 0 ) RETURN
C
      KK = ELREGS(NUMEDE)
C
      DO I = 1, NUMEDE-1
C
         K = ELREGS(I)
C
         IF( K .NE. KK ) CALL SETFLG( BIGLST, NODLST, IEDG4, 3,.FALSE. )
         if( k .ne. kk .and. isinsd ) then
            print*,'*** INCONSISTENCY IN SPLEDG ',iedg
            print*,(elregs(k),k=1,numede)
            print*,i,isgeom
            print*,gtndfl(BIGLST,NODLST,nods(1),2),
     :             gtndfl(BIGLST,NODLST,nods(1),3),
     :             gtndfl(BIGLST,NODLST,nods(1),4)
            print*,gtndfl(BIGLST,NODLST,nods(2),2),
     :             gtndfl(BIGLST,NODLST,nods(2),3),
     :             gtndfl(BIGLST,NODLST,nods(2),4)
            stop
         end if
         KK = K
C
         IEDG3 = IEDG4
C
         NODS1(1) = INOD
         NODS1(2) = NDPTRS(I)
C
         INEW = -1
         CALL ADDEDG( BIGLST, NODLST, NODS1, IEDG4, INEW )
         if( goterr ) print*,'   added edge ',iedg4,', nodes ',nods1
C
         CALL SETFLG( BIGLST, NODLST, IEDG4, 3, .TRUE. )
c         IF( ISINSD ) THEN
c            CALL SETFLG(BIGLST,NODLST,IEDG4,3,.TRUE.)
c         ELSE IF( SURFAC .EQ. 0 ) THEN
c            CALL SETFLG(BIGLST,NODLST,IEDG4,3,.TRUE.)
c         ELSE IF( GTNDFL(BIGLST,NODLST,NODS1(2),3) ) THEN
c            CALL SETFLG(BIGLST,NODLST,IEDG4,3,.TRUE.)
c         ELSE IF( I .GT. 1 ) THEN
c            CALL SETFLG(BIGLST,NODLST,IEDG4,3,.TRUE.)
c         END IF
C
         EDGLST(1) = EDO1
         EDGLST(2) = EDPTRS(1,I)
         EDGLST(3) = EDPTRS(3,I)
         EDGLST(4) = IEDG1
         EDGLST(5) = IEDG4
         EDGLST(6) = IEDG3
C
         CONELM(2) = LSTEL1
         CONELM(1) = ELPTRS(2,I)
         CONELM(3) = 0
         CONELM(4) = 0
C
         INEW = -1
         CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, LSTEL1, INEW, K )
      if( ierr.eq.12345 ) then
         print*,'^^^ inconsistent oppedg from addele in spledg (1) ^^^'
         ierr = 0
      else if(ierr .lt. 0) then
         print*,'*** SPLEDG: Got error from ADDELE (1): ',IERR
         return
      end if
         if( goterr ) then
c            call elnods( BIGLST,NODLST, lstel1, nodt, nxt, lst )
            print*,'   added elm ',lstel1
c     :        elmfnc( BIGLST, NODLST, lstel1,
c     :                nodt(1),nodt(2),nodt(3),nodt(4),rad)
c              biglst(7,lstel1) = 0
c              nxt = biglst(nxtbig,lstel1)
c              biglst(7,nxt) = 0
         end if
C
         EDGLST(1) = EDO2
         EDGLST(2) = EDPTRS(2,I)
         EDGLST(3) = EDPTRS(3,I)
         EDGLST(4) = IEDG2
         EDGLST(5) = IEDG4
         EDGLST(6) = IEDG3
C
         CONELM(2) = LSTEL2
         CONELM(1) = ELPTRS(3,I)
         CONELM(3) = LSTEL1
         CONELM(4) = 0
C
         INEW = -1
         CALL ADDELE( BIGLST,NODLST, CONELM, EDGLST, LSTEL2, INEW, K )
      if( ierr.eq.12345 ) then
         print*,'^^^ inconsistent oppedg from addele in spledg (2) ^^^'
         ierr = 0
      else if(ierr .lt. 0) then
         print*,'*** SPLEDG: Got error from ADDELE (2): ',IERR
         return
      end if
         if( goterr ) then
c            call elnods( BIGLST, NODLST, lstel2, nodt, nxt, lst )
            print*,'   added elm ',lstel2
c     :        elmfnc( BIGLST, NODLST, lstel2,
c     :                nodt(1),nodt(2),nodt(3),nodt(4),rad)
c              biglst(7,lstel2) = 0
c              nxt = biglst(nxtbig,lstel2)
c              biglst(7,nxt) = 0
         end if
c         if( goterr ) print*,'   added elm ',lstel2
C
         EDO1 = EDPTRS(1,I)
         EDO2 = EDPTRS(2,I)
C
         IF( I .EQ. 1 ) THEN
            NXTEL1 = LSTEL1
            NXTEL2 = LSTEL2
         END IF
C
      END DO
C
      I = NUMEDE
C
      K = ELREGS(I)
C
      IF( K .NE. KK ) CALL SETFLG( BIGLST, NODLST, IEDG4, 3, .FALSE. )
         if( k .ne. kk .and. isinsd ) then
            print*,'*** INCONSISTENCY IN SPLEDG ',iedg
            print*,(elregs(k),k=1,numede)
            print*,i,isgeom
            print*,gtndfl(BIGLST,NODLST,nods(1),2),
     :             gtndfl(BIGLST,NODLST,nods(1),3),
     :             gtndfl(BIGLST,NODLST,nods(1),4)
            print*,gtndfl(BIGLST,NODLST,nods(2),2),
     :             gtndfl(BIGLST,NODLST,nods(2),3),
     :             gtndfl(BIGLST,NODLST,nods(2),4)
            stop
         end if
      KK = K
C
      IEDG3 = IEDG4
C
      NODS1(1) = INOD
      NODS1(2) = NDPTRS(I)
C
      IF( SURFAC .GT. 0 ) THEN
         INEW = -1
         CALL ADDEDG( BIGLST, NODLST, NODS1, IEDG4, INEW )
         if( goterr ) print*,'   added edge ',iedg4,', nodes ',nods1
         IF( ISINSD ) THEN
            print*,'*** SPLEDG: INCONSISTENCY IN EDGE/SURFACE INFO'
            PRINT*,ISINSD,SURFAC
            IERR = -100
C            CALL SETFLG( BIGLST, NODLST, IEDG4, 3, .TRUE. )
            RETURN
         ELSE IF( GTNDFL(BIGLST,NODLST,NODS1(2),3) ) THEN
            print*,'*** SPLEDG: INCONSISTENCY IN NODE/SURFACE INFO'
            PRINT*,ISINSD,SURFAC
            IERR = -100
C            CALL SETFLG( BIGLST, NODLST, IEDG4, 3, .TRUE. )
            RETURN
         END IF
         NXTEL1 = SURFC2
         NXTEL2 = NXTEL1
      ELSE
         IEDG4  = IEDG5
c         NXTEL1 = ELPTRS(2,1)
c         NXTEL2 = ELPTRS(3,1)
      END IF
C
      EDGLST(1) = EDO1
      EDGLST(2) = EDPTRS(1,I)
      EDGLST(3) = EDPTRS(3,I)
      EDGLST(4) = IEDG1
      EDGLST(5) = IEDG4
      EDGLST(6) = IEDG3
C
      CONELM(4) = NXTEL1
      CONELM(1) = ELPTRS(2,I)
      CONELM(3) = 0
      CONELM(2) = LSTEL1
C
      INEW = -1
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, LSTEL1, INEW, K )
      if( ierr.eq.12345 ) then
         print*,'^^^ inconsistent oppedg from addele in spledg (3) ^^^'
         ierr = 0
      else if(ierr .lt. 0) then
         print*,'*** SPLEDG: Got error from ADDELE (3): ',IERR
         return
      end if
c         if( goterr ) then
c            call elnods( BIGLST, NODLST, lstel1, nodt, nxt, lst )
c            print*,'   added elm ',lstel1,
c     :        elmfnc( BIGLST, NODLST, lstel1,
c     :                nodt(1),nodt(2),nodt(3),nodt(4),rad)
c              biglst(7,lstel1) = 0
c              nxt = biglst(nxtbig,lstel1)
c              biglst(7,nxt) = 0
c         end if
      if( goterr ) print*,'   added elm ',lstel1
C
      EDGLST(1) = EDO2
      EDGLST(2) = EDPTRS(2,I)
      EDGLST(3) = EDPTRS(3,I)
      EDGLST(6) = IEDG3
      EDGLST(4) = IEDG2
      EDGLST(5) = IEDG4
C
      CONELM(4) = NXTEL2
      CONELM(1) = ELPTRS(3,I)
      CONELM(3) = LSTEL1
      CONELM(2) = LSTEL2
C
      INEW = -1
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, LSTEL2, INEW, K )
      if( ierr.eq.12345 ) then
         print*,'^^^ inconsistent oppedg from addele in spledg (4) ^^^'
         ierr = 0
      else if(ierr .lt. 0) then
         print*,'*** SPLEDG: Got error from ADDELE (4): ',IERR
         return
      end if
c         if( goterr ) then
c            call elnods( BIGLST, NODLST, lstel2, nodt, nxt, lst )
c            print*,'   added elm ',lstel2,
c     :        elmfnc( BIGLST, NODLST, lstel2,
c     :                nodt(1),nodt(2),nodt(3),nodt(4),rad)
c              biglst(7,lstel2) = 0
c              nxt = biglst(nxtbig,lstel2)
c              biglst(7,nxt) = 0
c         end if
      if( goterr ) print*,'   added elm ',lstel2
C
      NEWEL1 = LSTEL1
      NEWEL2 = LSTEL2
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMCHG = TIMCHG + TIMEND - TIMSTT
      TOTCHG = TOTCHG + 1
c
c      if( goterr ) pause
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
     :                INEW, IFLAG, K, NXT, LST )
         DO J = 1, 4
            IF( CONELM(J) .EQ. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** SPLEDG: INCONSISTENT CONNECTED ELEMENTS'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,conelm(j)) .le. 0 ) then
               if( .not. addit ) then
                  print*,'*** SPLEDG: INVALID CONNECTED ELEMENT'
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
         x = elmfnc( BIGLST, NODLST, 0,
     :               nods(1), nods(2), nods(3), nods(4), y )
         if( abs(y) .lt. 1e-8 ) then
            if( .not. addit ) then
               PRINT*,'*** SPLEDG: SMALL ELEMENT'
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
                  PRINT*,'*** SPLEDG: INCONSISTENT EDGE LIST'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,edglst(j)) .ge. 0 ) then
               if( .not. addit ) then
                  PRINT*,'*** SPLEDG: INVALID EDGE'
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
                    PRINT*,'*** SPLEDG: INVALID GEOM-EDGE NODE 1 FLAGS'
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
                    PRINT*,'*** SPLEDG: INVALID GEOM-EDGE NODE 2 FLAGS'
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
                    PRINT*,'*** SPLEDG: INVALID SURF-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
               else if( gtndfl(BIGLST,NODLST,nods(2),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** SPLEDG: INVALID SURF-EDGE NODE 2 FLAGS'
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
         print*,'*** SPLEDG: INCONSISTENCY AT END OF BIGLST'
         print*,i,emtbig
         ierr = -10000
      end if
C
      RETURN
      END
C
