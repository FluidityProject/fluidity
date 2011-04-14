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
      SUBROUTINE EDG2FC( BIGLST, NODLST, IEDG, ELPTRS, EDPTRS, NDPTRS,
     :                   IREG, IWHAT )
C-----------------------------------------------------------------------
C
C - This subroutine replaces the edge supplied (IEDG) with a face (for the
C - case IWHAT=3) or two faces and one edge (for the case IWHAT=+/-4), and
C - also performs the necessary removal and joining of affected elements.
C - If IWHAT=4 then the new edge joins nodes 1 and 3; if IWHAT=-4 then the
C - new edge joins nodes 2 and 4.
C - This subroutine assumes there are three/four elements around an edge,
C - as appropriate for the value of IWHAT supplied (and that there are
C - no surfaces - i.e. EDGELS returned SURFAC as zero), so EDG2FC should
C - only be called after verifying that NUMEDE=3/4 after calling EDGELS.
C - It also assumes that the edge to be replaced (IEDG) does intersect
C - (one of) the face(s) which will replace it. This should always be the
C - case if there are three elements surrounding the edge, and should be
C - the case if EDGTST has cleared these elements for changing.
C - The necessary information in the arrays ELPTRS, EDPTRS, NDPTRS can
C - be filled in easily by calling subroutine EDGELS for the edge IEDG
C - before calling this subroutine.
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
      INCLUDE 'blktim.i'
C
      INTEGER IEDG, ELPTRS(3,4), EDPTRS(3,4), NDPTRS(4), IWHAT, IREG
C
      INTEGER NODS(4), EDGLST(6), CONELM(4), NUMEDE, IEDNW,
     :        INEW, IFLAG, NXT, LST, IELM1, IELM2, IELM3, IELM4, i, j
c
      real dx, dy, dz, elmfnc, v1, v2, nx, ny, nz
C
      LOGICAL GETFLG, GTNDFL, addit, goterr
C
      CALL ADAPT_SECOND( TIMSTT )
C
      NUMEDE = ABS(IWHAT)
c
      goterr = debug
c      if( iwhat .ne. 3 ) goterr = .true.
c      if( iedg .eq. 7235 ) goterr = .true.
c      if( iedg .eq. 944 ) goterr = .true.
c      if( iedg .eq. 26682 ) goterr = .true.
c      if( iedg .eq. 14519 ) goterr = .true.
c      if( iedg .eq. 26692 ) goterr = .true.
c      if( iedg .eq. 14425 ) goterr = .true.
c      if( iedg .eq. 21987 ) goterr = .true.
C
      IF( NUMEDE .NE. 3 .AND. NUMEDE .NE. 4 ) RETURN
C
      CALL EDGINF( BIGLST, NODLST, IEDG, NODS, INEW, IFLAG, NXT, LST )
C
      if( goterr ) then
        PRINT*,' '
        PRINT*,'Into EDG2FC... ',ierr
        print*,'    IEDG: ',IEDG,nods(1),nods(2)
        print*,'    els: ',(elptrs(1,i),i=1,numede)
        print*,'    eds: ',(edptrs(3,i),i=1,numede)
        print*,'    nds: ',(ndptrs(i),i=1,numede)
        print*,'  upels: ',(elptrs(2,i),i=1,numede)
        print*,'  lwels: ',(elptrs(3,i),i=1,numede)
        print*,'  upeds: ',(edptrs(1,i),i=1,numede)
        print*,'  lweds: ',(edptrs(2,i),i=1,numede)
      end if
C
      IF( GETFLG( BIGLST, NODLST, IEDG, 2 ) ) THEN
         PRINT*,'+++ Warning: attempt to swap a geometry edge +++'
         RETURN
      ELSE IF( .NOT. GETFLG( BIGLST, NODLST, IEDG, 3 ) ) THEN
         PRINT*,'+++ Warning: attempt to swap a surface edge +++'
         RETURN
      END IF
c
c      if( goterr ) then
c        dx = elmfnc( BIGLST, NODLST, 0, ndptrs(1), ndptrs(2), ndptrs(3), nods(1), v1 )
c        dy = elmfnc( BIGLST, NODLST, 0, ndptrs(1), ndptrs(2), ndptrs(3), nods(2), v2 )
c        print*,'   in-spheres of new elements: ',v1,v2
c        print*,'   node1: ',(nodlst(i,ndptrs(1)),i=1,3)
c        print*,'   node2: ',(nodlst(i,ndptrs(2)),i=1,3)
c        print*,'   node3: ',(nodlst(i,ndptrs(3)),i=1,3)
c        print*,'   topnd: ',(nodlst(i,nods(1)),i=1,3)
c        print*,'   botnd: ',(nodlst(i,nods(2)),i=1,3)
c        call facnrm( BIGLST, NODLST, ndptrs(1), ndptrs(2), ndptrs(3), nx, ny, nz )
c        dx = nodlst(1,nods(1)) - nodlst(1,ndptrs(1))
c        dy = nodlst(2,nods(1)) - nodlst(2,ndptrs(1))
c        dz = nodlst(3,nods(1)) - nodlst(3,ndptrs(1))
c        v1 = dx*nx + dy*ny + dz*nz
c        dx = nodlst(1,nods(2)) - nodlst(1,ndptrs(1))
c        dy = nodlst(2,nods(2)) - nodlst(2,ndptrs(1))
c        dz = nodlst(3,nods(2)) - nodlst(3,ndptrs(1))
c        v2 = dx*nx + dy*ny + dz*nz
c        print*,'   direction cosines: ',v1,v2
c      end if
C
      CALL RMVELE( BIGLST, NODLST, ELPTRS(1,1) )
      if( goterr ) print*,'   removed element ',elptrs(1,1)
      CALL RMVELE( BIGLST, NODLST, ELPTRS(1,2) )
      if( goterr ) print*,'   removed element ',elptrs(1,2)
      CALL RMVELE( BIGLST, NODLST, ELPTRS(1,3) )
      if( goterr ) print*,'   removed element ',elptrs(1,3)
      IF( NUMEDE .EQ. 4 ) THEN
         CALL RMVELE( BIGLST, NODLST, ELPTRS(1,4) )
         if( goterr ) print*,'   removed element ',elptrs(1,4)
      END IF
C
      CALL RMVEDG( BIGLST, NODLST, IEDG )
      if( goterr ) print*,'   removed edge ',iedg
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IWHAT .EQ. 3 ) THEN
C
        EDGLST(1) = EDPTRS(3,1)
        EDGLST(2) = EDPTRS(3,2)
        EDGLST(3) = EDPTRS(3,3)
        EDGLST(4) = EDPTRS(1,1)
        EDGLST(5) = EDPTRS(1,2)
        EDGLST(6) = EDPTRS(1,3)
C
        CONELM(2) = ELPTRS(2,1)
        CONELM(4) = ELPTRS(2,2)
        CONELM(3) = ELPTRS(2,3)
        CONELM(1) = 0
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM1, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm1
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (1)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (1)'
           return
        end if
C
        EDGLST(1) = EDPTRS(3,1)
        EDGLST(2) = EDPTRS(3,2)
        EDGLST(3) = EDPTRS(3,3)
        EDGLST(4) = EDPTRS(2,1)
        EDGLST(5) = EDPTRS(2,2)
        EDGLST(6) = EDPTRS(2,3)
C
        CONELM(2) = ELPTRS(3,1)
        CONELM(4) = ELPTRS(3,2)
        CONELM(3) = ELPTRS(3,3)
        CONELM(1) = IELM1
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM2, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm2
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (2)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (2)'
           return
        end if
C
      ELSE IF( IWHAT .EQ. 4 ) THEN
C
        NODS(1) = NDPTRS(1)
        NODS(2) = NDPTRS(3)
C
        INEW = -1
        CALL ADDEDG( BIGLST, NODLST, NODS, IEDNW, INEW )
        CALL SETFLG( BIGLST, NODLST, IEDNW, 3, .TRUE. )
        if( goterr ) print*,'   added edge ',iednw,nods(1),nods(2)
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,2)
        EDGLST(3) = EDPTRS(3,3)
        EDGLST(4) = EDPTRS(1,1)
        EDGLST(5) = EDPTRS(1,2)
        EDGLST(6) = EDPTRS(1,3)
C
        CONELM(4) = ELPTRS(2,2)
        CONELM(3) = ELPTRS(2,3)
        CONELM(2) = 0
        CONELM(1) = 0
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM1, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm1
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (3)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (3)'
           return
        end if
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,2)
        EDGLST(3) = EDPTRS(3,3)
        EDGLST(4) = EDPTRS(2,1)
        EDGLST(5) = EDPTRS(2,2)
        EDGLST(6) = EDPTRS(2,3)
C
        CONELM(4) = ELPTRS(3,2)
        CONELM(3) = ELPTRS(3,3)
        CONELM(1) = IELM1
        CONELM(2) = 0
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM2, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm2
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (4)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (4)'
           return
        end if
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,4)
        EDGLST(3) = EDPTRS(3,1)
        EDGLST(4) = EDPTRS(1,3)
        EDGLST(5) = EDPTRS(1,4)
        EDGLST(6) = EDPTRS(1,1)
C
        CONELM(4) = ELPTRS(2,4)
        CONELM(3) = ELPTRS(2,1)
        CONELM(2) = IELM1
        CONELM(1) = 0
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM3, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm3
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (5)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (5)'
           return
        end if
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,4)
        EDGLST(3) = EDPTRS(3,1)
        EDGLST(4) = EDPTRS(2,3)
        EDGLST(5) = EDPTRS(2,4)
        EDGLST(6) = EDPTRS(2,1)
C
        CONELM(4) = ELPTRS(3,4)
        CONELM(3) = ELPTRS(3,1)
        CONELM(1) = IELM3
        CONELM(2) = IELM2
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM4, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm4
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (6)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (6)'
           return
        end if
C
      ELSE IF( IWHAT .EQ. -4 ) THEN
C
        NODS(1) = NDPTRS(2)
        NODS(2) = NDPTRS(4)
C
        INEW = -1
        CALL ADDEDG( BIGLST, NODLST, NODS, IEDNW, INEW )
        CALL SETFLG( BIGLST, NODLST, IEDNW, 3, .TRUE. )
        if( goterr ) print*,'   added edge ',iednw,nods(1),nods(2)
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (7)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (7)'
           return
        end if
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,1)
        EDGLST(3) = EDPTRS(3,2)
        EDGLST(4) = EDPTRS(1,4)
        EDGLST(5) = EDPTRS(1,1)
        EDGLST(6) = EDPTRS(1,2)
C
        CONELM(4) = ELPTRS(2,1)
        CONELM(3) = ELPTRS(2,2)
        CONELM(2) = 0
        CONELM(1) = 0
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM1, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm1
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (8)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (8)'
           return
        end if
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,1)
        EDGLST(3) = EDPTRS(3,2)
        EDGLST(4) = EDPTRS(2,4)
        EDGLST(5) = EDPTRS(2,1)
        EDGLST(6) = EDPTRS(2,2)
C
        CONELM(4) = ELPTRS(3,1)
        CONELM(3) = ELPTRS(3,2)
        CONELM(1) = IELM1
        CONELM(2) = 0
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM2, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm2
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (9)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (9)'
           return
        end if
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,3)
        EDGLST(3) = EDPTRS(3,4)
        EDGLST(4) = EDPTRS(1,2)
        EDGLST(5) = EDPTRS(1,3)
        EDGLST(6) = EDPTRS(1,4)
C
        CONELM(4) = ELPTRS(2,3)
        CONELM(3) = ELPTRS(2,4)
        CONELM(2) = IELM1
        CONELM(1) = 0
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM3, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm3
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (10)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (10)'
           return
        end if
C
        EDGLST(1) = IEDNW
        EDGLST(2) = EDPTRS(3,3)
        EDGLST(3) = EDPTRS(3,4)
        EDGLST(4) = EDPTRS(2,2)
        EDGLST(5) = EDPTRS(2,3)
        EDGLST(6) = EDPTRS(2,4)
C
        CONELM(4) = ELPTRS(3,3)
        CONELM(3) = ELPTRS(3,4)
        CONELM(1) = IELM3
        CONELM(2) = IELM2
C
        INEW = -1
        CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM4, INEW, IREG )
        if( goterr ) print*,'   added element ',ielm4
        if( ierr.eq.12345 ) then
           print*,'--- inconsistent oppedg from addele in edg2fc (11)'
           ierr = 0
        else if( ierr.ne.0 ) then
           print*,'*** EDG2FC: Got an error from ADDELE (11)'
           return
        end if
C
      END IF
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
         CALL ELMINF( BIGLST, NODLST, I, CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
         DO J = 1, 4
            IF( CONELM(J) .EQ. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** EDG2FC: INCONSISTENT CONNECTED ELEMENTS'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,conelm(j)) .le. 0 ) then
               if( .not. addit ) then
                  print*,'*** EDG2FC: INVALID CONNECTED ELEMENT'
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
         dx = elmfnc( BIGLST, NODLST, 0,
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
            nods(1) = biglst(1,edglst(j))
            nods(2) = biglst(2,edglst(j))
            IF( edglst(J) .LE. 0 ) THEN
               if( .not. addit ) then
                  PRINT*,'*** EDG2FC: INCONSISTENT EDGE LIST'
                  print*,'   element: ',i
                  print*,'   conelm: ',conelm
                  print*,'   edglst: ',edglst
                  addit  = .true.
                  goterr = .true.
               end if
            else if( biglst(4,edglst(j)) .ge. 0 ) then
               if( .not. addit ) then
                  PRINT*,'*** EDG2FC: INVALID EDGE'
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
                    PRINT*,'*** EDG2FC: INVALID GEOM-EDGE NODE 1 FLAGS'
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
                    PRINT*,'*** EDG2FC: INVALID GEOM-EDGE NODE 2 FLAGS'
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
                    PRINT*,'*** EDG2FC: INVALID SURF-EDGE NODE 1 FLAGS'
                    print*,'   element: ',i
                    print*,'   conelm: ',conelm
                    print*,'   edglst: ',edglst
                    addit  = .true.
                    goterr = .true.
                  end if
                  print*,'   invalid edge: ',edglst(j),nods(1),nods(2)
               else if( gtndfl(BIGLST,NODLST,nods(2),3) ) then
                  if( .not. addit ) then
                    PRINT*,'*** EDG2FC: INVALID SURF-EDGE NODE 2 FLAGS'
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
         print*,'*** EDG2FC: INCONSISTENCY AT END OF BIGLST'
         print*,i,emtbig
         ierr = -10000
      end if
C
      RETURN
      END
C
