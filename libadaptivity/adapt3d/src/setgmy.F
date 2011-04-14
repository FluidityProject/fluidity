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
      SUBROUTINE SETGMY( BIGLST, NODLST, SRFGMY, CLCGMY,
     :                   NSELM, SZSNLS, SNLBAS, SNLIST, SURFID )
C-----------------------------------------------------------------------
C
C - SRFGMY should be TRUE if the surface mesh is to be left intact
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INTEGER NSELM, SZSNLS, NGMYND, NFACND, NSPLND,
     :        NGMYED, NFACED
C
      INTEGER SNLBAS(NSELM+1), SNLIST(SZSNLS), SURFID(NSELM)
C
      LOGICAL SRFGMY, CLCGMY
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IPTR, CONELM(4), EDGLST(6), IELM, INEW, IFLAG, NXT, LST,
     :        I, ELPTRS(3,MAXEDE), EDPTRS(3*MAXEDE+3), NDPTRS(MAXEDE+1),
     :        NODS(4), SURFAC, NUMEDE, I1, I2, I3, I4, INTRNL, IREG,
     :        NUMEDG, NUMELM, IEDG, J, CNTEDG, CNTNOD, CNTELM,
     :        CNTGMY, CNTINT, CNTSPL, cnts(20), gmyed1, gmyed2,
     :        SURFC1, SURFC2, lstel1, lstel2, wchfac
C
      REAL NX1, NY1, NZ1, NX2, NY2, NZ2, DOT, CUTOFF, TMPVAL, cuts(20)
C
      LOGICAL GETFLG, GTNDFL, isgeom, issplt
C
      PARAMETER( CUTOFF = 0.999999 )
C
      cuts(1)  = 0.9
      cuts(2)  = 0.95
      cuts(3)  = 0.98
      cuts(4)  = 0.99
      cuts(5)  = 0.995
      cuts(6)  = 0.998
      cuts(7)  = 0.999
      cuts(8)  = 0.9995
      cuts(9)  = 0.9998
      cuts(10) = 0.9999
      cuts(11) = 0.99995
      cuts(12) = 0.99998
      cuts(13) = 0.99999
      cuts(14) = 0.999995
      cuts(15) = 0.999998
      cuts(16) = 0.999999
      cuts(17) = 0.9999995
      cuts(18) = 0.9999998
      cuts(19) = 0.9999999
      cuts(20) = 0.99999995
C
      DO I = 1, 20
         CNTS(I) = 0
      END DO
C
      ewrite(2,*) 
     :     'Setting the edge geometry information...'
C
C - Loop through edges, clearing bits four and eight of the flag
C
      CNTEDG = 0
      CNTELM = 0
      IPTR = STTBIG
C
  10  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
         NXT = BIGLST(NXTBIG,IPTR)
         CALL SETFLG( BIGLST, NODLST, IPTR, 4, .FALSE. )
         CALL SETFLG( BIGLST, NODLST, IPTR, 8, .FALSE. )
         CNTEDG = CNTEDG + 1
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
         CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
         CNTELM = CNTELM + 1
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 10
C
      ewrite(2,*) 
     :     'Edge flags cleared: ',cntedg
      ewrite(2,*) 
     :     'Elements counted: ',cntelm
C
C - Loop through nodes, setting bits six and nine of the flag
C - Note that bit six being set (TRUE) indicates the node has to be
C - checked, whereas bit nine being set indicates this node has not
C - been checked before.
C
      CNTNOD = 0
      IPTR = STTNOD
C
 15   CALL STNDFL(BIGLST,NODLST,IPTR,6,.TRUE.)
      CALL STNDFL(BIGLST,NODLST,IPTR,9,.TRUE.)
      NXT = INT(NODLST(NXTNOD,IPTR))
      CNTNOD = CNTNOD + 1
C
      IPTR = NXT
C
      IF( GTNDFL(BIGLST,NODLST,IPTR,1) ) GOTO 15
c
      ewrite(2,*) 
     :     'Node flags cleared: ',cntnod
C
C - Loop through elements, picking out edges
C
      CNTGMY = 0
      CNTINT = 0
      IPTR = STTBIG
C
  20  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
C
         NXT = BIGLST(NXTBIG,IPTR)
C
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
C
         CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
c         print*,'Checking element ',iptr
C
         DO I = 1, 6
C
            IF( .NOT. GETFLG(BIGLST,NODLST,EDGLST(I),4) ) THEN
C
C - for edges not already done, find the surrounding elements
C
               IEDG = EDGLST(I)
C
               CALL EDGELS( BIGLST, NODLST, IEDG, IPTR,
     :                      ELPTRS, EDPTRS, NDPTRS, NUMEDE,
     :                      SURFAC, INTRNL )
C
               IF( IERR .NE. 0 ) THEN
                  WRITE(0,*) '*** SETGMY: GOT ERROR FROM EDGELS'
                  RETURN
               END IF
c               print*,'Edge ',IEDG,numede,surfac,intrnl
               I1 = BIGLST(1,IEDG)
               I2 = BIGLST(2,IEDG)
C
C - and work out what sort of edge it is (internal/surface/geometry)
C
               CALL SETFLG( BIGLST, NODLST, IEDG, 2, .FALSE. )
               IF( SURFAC .GT. 0 ) THEN
C
                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .FALSE. )
C
C - this edge is on a surface, so find the surfaces at each end of
C - the list of elements returned by EDGELS
C
                  surfc1 = wchfac( biglst, nodlst, elptrs(1,1),
     :                    edptrs(numede*3+1), edptrs(numede*3+2), iedg )
c     :                    EDPTRS(1,NUMEDE+1), EDPTRS(2,NUMEDE+1), iedg )
                  if( ierr .ne. 0 ) then
                     WRITE(0,*) '***SETGMY: Got error from WCHFAC (1)'
                     return
                  end if
                  surfc2 = wchfac( biglst, nodlst, elptrs(1,numede),
     :                    edptrs(numede*3-2), edptrs(numede*3-1), iedg )
c     :                    edptrs(1,numede), edptrs(2,numede), iedg )
                  if( ierr .ne. 0 ) then
                     WRITE(0,*) '***SETGMY: Got error from WCHFAC (1)'
                     return
                  end if
C
                  LSTEL1 = biglst(surfc1,elptrs(1,1))
                  if( surfc1 .eq. 4 ) LSTEL1 = LSTEL1/2
                  if(lstel1 .ge. 0) then
                     WRITE(0,*) '*** SETGMY: First face not a surface!',
     :                                                  lstel1,surfc1
                     surfc1 = -surfac
                  else
                     surfc1 = lstel1
                  endif
C
                  LSTEL2 = biglst(surfc2,elptrs(1,numede))
                  if( surfc2 .eq. 4 ) LSTEL2 = LSTEL2/2
                  if(LSTEL2 .ge. 0) then
                     WRITE(0,*) '*** SETGMY: Last face not a surface!',
     :                                                 lstel2,surfc2
                     surfc2 = -surfac
                  else
                     surfc2 = lstel2
                  endif
C
                  IF( surfc1 .ne. surfc2 ) THEN
C - surface IDs are different, so must be geometry edge
                     CNTGMY = CNTGMY + 1
                     CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                     IF( SRFGMY ) THEN
                        IF( .NOT. GETFLG(BIGLST,NODLST,IEDG,6) ) 
     :                    CALL SETFLG( BIGLST, NODLST, IEDG, 6, .TRUE. )
                     END IF
                     CALL SETFLG( BIGLST, NODLST, IEDG, 8, .TRUE. )
                     BIGLST( 7, IEDG ) = 1
c                     print*,'Geometry edge: ',i1,i2,dot
c                  ELSE IF( SRFGMY ) THEN
c                     CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
c                     CNTGMY = CNTGMY + 1
                  ELSE IF( INTRNL .LT. 0 ) THEN
C - internal boundary, as well as surface, so must be geometry edge
                     CNTGMY = CNTGMY + 1
                     CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                     IF( SRFGMY ) THEN
                        IF( .NOT. GETFLG(BIGLST,NODLST,IEDG,6) ) 
     :                    CALL SETFLG( BIGLST, NODLST, IEDG, 6, .TRUE. )
                     END IF
                     CALL SETFLG( BIGLST, NODLST, IEDG, 8, .TRUE. )
                     BIGLST( 7, IEDG ) = 1
c                     print*,'Geometry edge: ',i1,i2,dot
c                  ELSE IF( SRFGMY ) THEN
c                     CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
c                     CNTGMY = CNTGMY + 1
                  ELSE
C - check angle between faces, to see if it should be a geometry edge
                     I1 = BIGLST(1,IEDG)
                     I2 = BIGLST(2,IEDG)
                     I3 = NDPTRS(NUMEDE)
                     I4 = NDPTRS(NUMEDE+1)
                     CALL FACNRM( BIGLST, NODLST, I1, I2, I3,
     :                            NX1, NY1, NZ1 )
                     CALL FACNRM( BIGLST, NODLST, I2, I1, I4,
     :                            NX2, NY2, NZ2 )
                     IF( IERR .NE. 0 ) THEN
                       WRITE(0,*)'*** SETGMY: GOT ERROR FROM FACNRM (1)'
                       RETURN
                     END IF
                     DOT = NX1*NX2 + NY1*NY2 + NZ1*NZ2
                     BIGLST( 7, IEDG ) = INT(ABS(DOT)*100000000)+1
                     call docnts( dot, cnts, cuts )
                     IF( SRFGMY ) THEN
                        IF( .NOT. GETFLG(BIGLST,NODLST,IEDG,6) ) 
     :                    CALL SETFLG( BIGLST, NODLST, IEDG, 6, .TRUE. )
                     END IF
                     IF( ABS(DOT) .LT. CUTOFF .AND. CLCGMY ) THEN
                        CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                        CNTGMY = CNTGMY + 1
c                        print*,'Geometry edge: ',i1,i2,dot
c                     else
c                        print*,'Surface edge: ',i1,i2,dot
                     END IF
                  END IF
               ELSE IF( INTRNL .LT. -1 ) THEN
C - internal edge, but got more than one internal surface
C - so check angle between them to see if it should be a geometry edge
                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .FALSE. )
C - first two nodes of internal faces are from the current edge
                  I1 = BIGLST(1,IEDG)
                  I2 = BIGLST(2,IEDG)
C - INTRNL and SURFAC point to the far nodes for the internal faces
                  I3 = NDPTRS(-INTRNL-1)
                  I4 = NDPTRS(-SURFAC)
                  CALL FACNRM( BIGLST, NODLST, I1, I2, I3,
     :                         NX1, NY1, NZ1 )
                  CALL FACNRM( BIGLST, NODLST, I2, I1, I4,
     :                         NX2, NY2, NZ2 )
                  IF( IERR .NE. 0 ) THEN
                     WRITE(0,*) '*** SETGMY: GOT ERROR FROM FACNRM (2)'
                     RETURN
                  END IF
                  DOT = NX1*NX2 + NY1*NY2 + NZ1*NZ2
                  BIGLST( 7, IEDG ) = INT(ABS(DOT)*100000000)+1
                  call docnts( dot, cnts, cuts )
                  IF( ABS(DOT) .LT. CUTOFF ) THEN
                     CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                     CNTGMY = CNTGMY + 1
c                     print*,'Int.geom. edge: ',i1,i2,dot
c                  else
c                     print*,'Int.srf. edge: ',i1,i2,dot
                  END IF
               ELSE IF( INTRNL .EQ. -1 ) THEN
C - more than two internal surfaces, so must be a geometry edge
                  CALL SETFLG( BIGLST, NODLST, IEDG, 2, .TRUE. )
                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .FALSE. )
                  CALL SETFLG( BIGLST, NODLST, IEDG, 8, .TRUE. )
                  BIGLST( 7, IEDG ) = 1
                  CNTGMY = CNTGMY + 1
c                  print*,'Geometry edge: ',i1,i2
               ELSE
C - not a geometry edge, just a regular internal edge
                  CALL SETFLG( BIGLST, NODLST, IEDG, 3, .TRUE. )
                  BIGLST( 7, IEDG ) = 0
                  CNTINT = CNTINT + 1
c                  print*,'Internal edge: ',i1,i2
               END IF
C
               CALL SETFLG( BIGLST, NODLST, IEDG, 4, .TRUE. )
C
            END IF
C
         END DO
C
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 20
C
      ewrite(2,*) 'Geometry edges: ',CNTGMY
      ewrite(2,*) 'Internal edges: ',CNTINT
C
      IPTR = STTBIG
C
      ewrite(2,*) 
     :     'Setting the node geometry information...'
C
      CNTGMY = 0
      CNTINT = 0
      CNTSPL = 0
C
C - Loop through elements, picking out nodes
C
  30  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
C
         NXT = BIGLST(NXTBIG,IPTR)
C
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
C
         CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, .false. )
c         print*,'Checking element ',iptr
C
         DO I = 1, 4
C
C - Note that bit six being set (TRUE) indicates the node has to be
C - checked, whereas bit nine being set indicates this node has not
C - been checked before.
C
            IF( GTNDFL(BIGLST,NODLST,NODS(I),6) ) THEN
C
C - for nodes not already done, find list of surrounding elements/edges/nodes
C
               CALL NDCNCT( BIGLST, NODLST, NODS(I), IPTR,
     :                      ELPTRS, EDPTRS, NDPTRS, NUMEDG, NUMELM )
C
               IF( IERR .NE. 0 ) RETURN
C
C - check to see if we have visited this node before, and decrement
C - appropriate counters if we have.
C
               IF( .NOT. GTNDFL(BIGLST,NODLST,NODS(I),9) ) THEN
                  ewrite(3,*) 
     :                 '   Check again for node ',nods(i),
     :                   gtndfl(biglst,nodlst,nods(i),2),
     :                   gtndfl(biglst,nodlst,nods(i),3),
     :                   gtndfl(biglst,nodlst,nods(i),4)
                  IF( GTNDFL(BIGLST,NODLST,NODS(I),2) )
     :                                 CNTGMY = CNTGMY - 1
                  IF( GTNDFL(BIGLST,NODLST,NODS(I),3) )
     :                                 CNTINT = CNTINT - 1
                  IF( GTNDFL(BIGLST,NODLST,NODS(I),4) )
     :                                 CNTSPL = CNTSPL - 1
               END IF
C
C - work out what sort of node it is from the types of edges.
C
               I2 = 0
               I3 = 0
               I4 = 0
C
               gmyed1 = 0
               gmyed2 = 0
c
               DO I1 = 1, NUMEDG
                  IF( GETFLG(BIGLST,NODLST,EDPTRS(I1),2) ) THEN
                     I2 = I2 + 1
                     if( gmyed1 .gt. 0 ) then
                        gmyed2 = i1
                     else
                        gmyed1 = i1
                     end if
                     IF( GETFLG(BIGLST,NODLST,EDPTRS(I1),3) ) THEN
                        WRITE(0,*)
     :                        '+++ SETGMY: GOT INTERNAL GEOMETRY EDGE?'
                        I3 = I3 + 1
                     END IF
                  ELSE IF( GETFLG(BIGLST,NODLST,EDPTRS(I1),3) ) THEN
                     I3 = I3 + 1
                  END IF
               END DO
C
C - set up the default node flags (and flag that we have visited it)
C
               CALL STNDFL( BIGLST, NODLST, NODS(I), 2, .FALSE. )
               CALL STNDFL( BIGLST, NODLST, NODS(I), 3, .FALSE. )
               CALL STNDFL( BIGLST, NODLST, NODS(I), 4, .FALSE. )
               CALL STNDFL( BIGLST, NODLST, NODS(I), 6, .FALSE. )
C
               IF( I2 .GT. 2 ) THEN
C
C - node has more than two geometry edges, so this node cannot be changed.
C
                  CALL STNDFL( BIGLST, NODLST, NODS(I), 2, .TRUE. )
c                  print*,'Geometry node: ',nods(i),(nodlst(j,nods(i)),
c     :                                                            j=1,3)
c                  PRINT*,'Geom,intrnl,tot edges: ',i2,i3,numedg
                  CNTGMY = CNTGMY + 1
               ELSE IF( I3 .EQ. NUMEDG ) THEN
C
C - node has all internal edges, so this must be an internal node.
C
                  CALL STNDFL( BIGLST, NODLST, NODS(I), 3, .TRUE. )
c                  print*,'Internal node: ',nods(i),(nodlst(j,nods(i)),
c     :                                                            j=1,3)
                  CNTINT = CNTINT + 1
               ELSE IF( I2 .EQ. 2 ) THEN
C
C - found two geometry edges, so this is a geometry edge splitter node.
C - (The edges should only lie on a line - this is not checked here.
C -  If they do not then a weird error has occurred!)
C
                  CALL gmylin( BIGLST, NODLST, NODS(I), IPTR,
     :                         ELPTRS, EDPTRS, NDPTRS, NUMEDG, NUMELM,
     :                         gmyed1, gmyed2, isgeom, issplt, nxt )
c
                  if( ierr .ne. 0 ) return
c
                  if( issplt ) then
                     CALL STNDFL( BIGLST, NODLST, NODS(I), 4, .TRUE. )
                     CNTSPL = CNTSPL + 1
                  else if( isgeom ) then
c                     PRINT*,'+++ SETGMY: Bad splitter -made geom node!',
c     :                      NODS(I)
                     CALL STNDFL( BIGLST, NODLST, NODS(I), 2, .TRUE. )
                     CNTGMY = CNTGMY + 1
                  else
c                     PRINT*,'+++ SETGMY: Bad splitter -made surf node!',
c     :                      NODS(I)
                  end if
c                  print*,'Splitter node: ',nods(i),(nodlst(j,nods(i)),
c     :                                                            j=1,3)
c                  PRINT*,'Geom,intrnl,tot edges: ',i2,i3,numedg
               ELSE IF( I2 .EQ. 1 ) THEN
C
C - found only one geometry edge: impossible!
C
                  call gmybad( biglst, nodlst, nods(i), iptr,
     :                         elptrs, edptrs, ndptrs, numedg, numelm,
     :                         gmyed1, isgeom, issplt, nxt )
c
                  if( ierr .ne. 0 ) return
c
                  if( issplt ) then
c                    PRINT*,'+++ SETGMY: Bad geometry -made split node!',
c     :                      NODS(I)
                     CALL STNDFL( BIGLST, NODLST, NODS(I), 4, .TRUE. )
                     CNTSPL = CNTSPL + 1
                  else if( isgeom ) then
c                     PRINT*,'+++ SETGMY: Bad geometry -made geom node!',
c     :                      NODS(I)
                     CALL STNDFL( BIGLST, NODLST, NODS(I), 2, .TRUE. )
                     CNTGMY = CNTGMY + 1
c                  else
c                     PRINT*,'+++ SETGMY: Bad geometry -made surf node!',
c     :                      NODS(I)
                  end if
c
               ELSE
C
C - this node lies on a surface (not all edges are internal).
C
               END IF
C
               CALL STNDFL( BIGLST, NODLST, NODS(I), 9, .FALSE. )
C
            END IF
C
         END DO
C
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 30
C
      ewrite(2,*) 'Geometry nodes: ',CNTGMY
      ewrite(2,*) 'Splitter nodes: ',CNTSPL
      ewrite(2,*) 'Internal nodes: ',CNTINT
c
      ewrite(2,*) 'Count of dot products:'
      do i = 1, 20
        ewrite(2,*) 
     :        'Greater than ',cuts(21-i),' :   ',cnts(21-i)
      end do
C
c      PRINT*,'Leaving SETGMY'
c      print*,'Geometry setup completed!'
C
      RETURN
      END
C
c
c
      subroutine docnts( dot, cnts, cuts )
c
      implicit none
c
      integer cnts(20)
      real dot, cuts(20)
c
      integer i
c
      i = 21
c
  777 i = i - 1
c
      if( abs(dot) .gt. cuts(i) ) then
         cnts(i) = cnts(i) + 1
         return
      else if( i .gt. 1 ) then
         goto 777
      end if
c
      return
      end
c
c
c
      subroutine gmybad( biglst, nodlst, inod, ielm,
     :                   elptrs, edptrs, ndptrs, numedg, numelm,
     :                   gmyed1, isgeom, issplt, nxtel )
c
      implicit none
c
      integer inod, ielm, gmyed1, numelm, numedg, nxtel,
     :        elptrs(numelm), edptrs(numedg), ndptrs(numedg)
c
      logical isgeom, issplt
c
      include 'blknew.i'
c
      include 'blkbig.i'
c
      include 'blkerr.i'
c
      integer i, j, k, l, a1, minj, conelm(4), edglst(6), ireg, iflag,
     :        inew, nxt, lst, hasedg
c
      logical getflg
c
      isgeom = .false.
      issplt = .false.
c
c - look for largest angle of surrounding edges (except geometry edge)
c
      minj = 200000000
      k    = 0
c
      do i = 1, numedg
c
         if( edptrs(i) .le. 0 ) then
            WRITE(0,*) '*** GMYBAD: FOUND BAD EDGE!',i
            WRITE(0,*) (edptrs(j),j=1,numedg)
            WRITE(0,*) (ndptrs(j),j=1,numedg)
            WRITE(0,*) gmyed1,numedg,numelm,inod,ielm
            ierr = -5678
            return
         end if
c
         if( i .ne. gmyed1 ) then
            j = biglst( 7, edptrs(i) )
            if( j .lt. minj .and. j .gt. 0 ) then
               minj = j
               k = i
            end if
         end if
c
      end do
c
      a1 = biglst( 7, edptrs(gmyed1) )
      if( a1 .lt. 99900000 .or.
     :    getflg(biglst,nodlst,edptrs(gmyed1),8) ) then
c         print*,'+++ GMYBAD: Force gmynd, edge: ',edptrs(gmyed1),a1
         if( k .eq. 0 ) then
            WRITE(0,*) '*** GMYBAD: COULD NOT FIND AN EDGE AS GEOMETRY'
            WRITE(0,*) 'edge: ',edptrs(gmyed1),a1
            WRITE(0,*) 'numedg,numelm: ',numedg,numelm
            WRITE(0,*) edptrs
            WRITE(0,*) ndptrs
            WRITE(0,*) elptrs
            do i = 1, numedg
               a1 = biglst(7,edptrs(i))
               WRITE(0,*) 'ed:',edptrs(i),a1,
     :                      getflg(biglst,nodlst,edptrs(i),2),
     :                      getflg(biglst,nodlst,edptrs(i),3),
     :                      getflg(biglst,nodlst,edptrs(i),8)
            end do
            ierr = -4020
            return
         end if
         minj = 1
      end if
c
      if( minj .lt. 99999995 .and. minj .gt. 0 ) then
         isgeom = .true.
         call setflg( biglst, nodlst, edptrs(k), 2, .TRUE.  )
         call setflg( biglst, nodlst, edptrs(k), 3, .FALSE. )
         call setflg( biglst, nodlst, edptrs(k), 8, .TRUE.  )
         call stndfl( biglst, nodlst, ndptrs(k), 6, .TRUE. )
         call stndfl( biglst, nodlst, inod,      6, .TRUE. )
         WRITE(0,*) '+++ GMYBAD: Turned into geometry node'
         WRITE(0,*) '  olded,newed:',edptrs(gmyed1),edptrs(k)
         WRITE(0,*) '  node,othnd: ',inod,ndptrs(k)
         j = 1
 350     if( hasedg(biglst,nodlst,elptrs(j),edptrs(k)).eq.0 ) then
            j = j + 1
            if( j .le. numelm ) goto 350
            WRITE(0,*) '*** GMYBAD: REACHED END OF ELEMENT LIST (1)'
            WRITE(0,*) 'els:',(elptrs(l),l=1,numelm)
            WRITE(0,*) 'eds:',(edptrs(l),l=1,numedg)
            WRITE(0,*) 'nds:',(ndptrs(l),l=1,numedg)
            WRITE(0,*) 'minj,k,gmyed1: ',minj,k,gmyed1
            WRITE(0,*) 'inod,ielm,numedg,numelm: ',
     :                                  inod,ielm,numedg,numelm
            WRITE(0,*) 'nodes of gmyed1: ',
     :                           (biglst(j,edptrs(gmyed1)),j=1,2)
            do j = 1, numelm
               call elminf(biglst,nodlst,elptrs(j),conelm,edglst,
     :                     inew,iflag,ireg,nxt,lst)
               WRITE(0,*) 'eds of',elptrs(j),' : ',edglst
               WRITE(0,*) 'cncted els: ',conelm
            end do
            ierr = -5555
            return
         end if
         call el2end( biglst, nodlst, elptrs(j), nxtel )
      else
         call setflg( biglst, nodlst, edptrs(gmyed1), 2, .FALSE. )
         call setflg( biglst, nodlst, edptrs(gmyed1), 8, .FALSE. )
         call stndfl( biglst, nodlst, ndptrs(gmyed1), 6, .TRUE. )
         j = 1
 450     if( hasedg(biglst,nodlst,elptrs(j),edptrs(gmyed1)).eq.0 ) then
            j = j + 1
            if( j .le. numelm ) goto 450
            WRITE(0,*) '*** GMYBAD: REACHED END OF ELEMENT LIST (2)'
            WRITE(0,*) (elptrs(l),l=1,numelm)
            WRITE(0,*) (edptrs(l),l=1,numedg)
            WRITE(0,*) (ndptrs(l),l=1,numedg)
            WRITE(0,*) gmyed1,numedg,numelm,inod,ielm,k,minj,j
            ierr = -5556
            return
         end if
         call el2end( biglst, nodlst, elptrs(j), nxtel )
         WRITE(0,*) '+++ GMYBAD: Turned into surface node'
         WRITE(0,*) '  edge,elm: ',edptrs(gmyed1),elptrs(j)
         WRITE(0,*) '  node,othnd: ',inod,ndptrs(gmyed1)
      end if
c
      return
      end
c
c
c
      subroutine gmylin( biglst, nodlst, inod, ielm,
     :                   elptrs, edptrs, ndptrs, numedg, numelm,
     :                   gmyed1, gmyed2, isgeom, issplt, nxtel )
c
      implicit none
c
      integer inod, ielm, gmyed1, gmyed2, numelm, numedg, nxtel,
     :        elptrs(numelm), edptrs(numedg), ndptrs(numedg)
c
      logical isgeom, issplt
c
      include 'blknew.i'
c
      include 'blkbig.i'
c
      include 'blkerr.i'
c
      integer i, j, k, l, a1, a2, inod1, inod2, minj, hasedg
c
      real dot, x1, y1, z1, x2, y2, z2
c
      logical getflg
c
      isgeom = .false.
      issplt = .false.
c
      inod1 = biglst( 1, edptrs(gmyed1) )
      if( inod1 .eq. inod ) inod1 = biglst( 2, edptrs(gmyed1) )
c
      inod2 = biglst( 1, edptrs(gmyed2) )
      if( inod2 .eq. inod ) inod2 = biglst( 2, edptrs(gmyed2) )
c
      x1 = nodlst(1,inod) - nodlst(1,inod1)
      y1 = nodlst(2,inod) - nodlst(2,inod1)
      z1 = nodlst(3,inod) - nodlst(3,inod1)
c
      dot = 1.0/sqrt(x1*x1 + y1*y1 + z1*z1)
      x1 = x1*dot
      y1 = y1*dot
      z1 = z1*dot
c
      x2 = nodlst(1,inod) - nodlst(1,inod2)
      y2 = nodlst(2,inod) - nodlst(2,inod2)
      z2 = nodlst(3,inod) - nodlst(3,inod2)
c
      dot = 1.0/sqrt(x2*x2 + y2*y2 + z2*z2)
      x2 = x2*dot
      y2 = y2*dot
      z2 = z2*dot
c
      dot = abs( x1*x2 + y1*y2 + z1*z2 )
c
      if( dot .gt. 0.99998 ) then
         issplt = .true.
         return
      else
         isgeom = .true.
         return
      end if
c
c - found edges are not in straight line so look for largest
c - angle of surrounding edges (other than the two geometry edges)
c
      minj = 200000000
      k    = 0
c
      do i = 1, numedg
c
         if( edptrs(i) .le. 0 ) then
            WRITE(0,*) '*** GMYLIN: FOUND BAD EDGE!',i
            WRITE(0,*) (edptrs(j),j=1,numedg)
            WRITE(0,*) (ndptrs(j),j=1,numedg)
            WRITE(0,*) gmyed1,gmyed2,numedg,numelm,inod,ielm
            ierr = -5679
            return
         end if
c
         if( i .ne. gmyed1 .and. i .ne. gmyed2 ) then
            j = biglst( 7, edptrs(i) )
            if( j .lt. minj .and. j .gt. 0 ) then
               minj = j
               k = i
            end if
         end if
c
      end do
c
      a1 = biglst( 7, edptrs(gmyed1) )
      a2 = biglst( 7, edptrs(gmyed2) )
      if( a1 .lt. 99900000 .or. a2 .lt. 99900000 .or.
     :    getflg(biglst,nodlst,edptrs(gmyed1),8) .or.
     :    getflg(biglst,nodlst,edptrs(gmyed2),8) ) then
         if( k .eq. 0 ) then
            WRITE(0,*) '*** GMYLIN: COULD NOT FIND AN EDGE AS GEOMETRY'
            WRITE(0,*) 'edges angle: ',dot
            WRITE(0,*) 'geom edge angles: ',a1,a2
            WRITE(0,*) 'edges: ',edptrs(gmyed1),edptrs(gmyed2)
            WRITE(0,*) 'numedg,numelm: ',numedg,numelm
            WRITE(0,*) edptrs
            WRITE(0,*) ndptrs
            WRITE(0,*) elptrs
            do i = 1, numedg
               a1 = biglst(7,edptrs(i))
               WRITE(0,*) 'ed:',edptrs(i),a1,
     :                      getflg(biglst,nodlst,edptrs(i),2),
     :                      getflg(biglst,nodlst,edptrs(i),3),
     :                      getflg(biglst,nodlst,edptrs(i),8)
            end do
            ierr = -4010
            return
         end if
         minj = 1
      end if
c
      if( minj .lt. 99999995 .and. minj .gt. 0 ) then
         isgeom = .true.
         call setflg( biglst, nodlst, edptrs(k), 2, .TRUE.  )
         call setflg( biglst, nodlst, edptrs(k), 3, .FALSE. )
         call setflg( biglst, nodlst, edptrs(k), 8, .TRUE.  )
         call stndfl( biglst, nodlst, ndptrs(k), 6, .TRUE. )
         j = 1
 350     if( hasedg(biglst,nodlst,elptrs(j),edptrs(k)).eq.0 ) then
            j = j + 1
            if( j .le. numelm ) goto 350
            WRITE(0,*) '*** GMYLIN: REACHED END OF ELEMENT LIST (1)'
            WRITE(0,*) (elptrs(l),l=1,numelm)
            WRITE(0,*) (edptrs(l),l=1,numedg)
            WRITE(0,*) (ndptrs(l),l=1,numedg)
            WRITE(0,*) gmyed1,numedg,numelm,inod,ielm,k,minj,j
            ierr = -5557
            return
         end if
         call el2end( biglst, nodlst, elptrs(j), nxtel )
         WRITE(0,*) '+++ GMYLIN: Turned into geometry node'
         WRITE(0,*) '  oldeds,newed:',edptrs(gmyed1),edptrs(gmyed2),
     :                            edptrs(k)
         WRITE(0,*) '  node,othnd,elm: ',inod,ndptrs(k),elptrs(j)
      else
         call setflg( biglst, nodlst, edptrs(gmyed1), 2, .FALSE. )
         call setflg( biglst, nodlst, edptrs(gmyed2), 2, .FALSE. )
         call stndfl( biglst, nodlst, ndptrs(gmyed1), 6, .TRUE. )
         call stndfl( biglst, nodlst, ndptrs(gmyed2), 6, .TRUE. )
         WRITE(0,*) '+++ GMYLIN: Turned into surface node'
         WRITE(0,*) '  edges:',edptrs(gmyed1),edptrs(gmyed2)
         WRITE(0,*) '  node,othnds: ',inod,ndptrs(gmyed1),ndptrs(gmyed2)
         j = 1
 450     if( hasedg(biglst,nodlst,elptrs(j),edptrs(gmyed1)).eq.0 ) then
            if( gmyed1 .ne. gmyed2 ) then
               if( hasedg(biglst,nodlst,elptrs(j),edptrs(gmyed2))
     :                                                     .gt.0 ) then
                  gmyed2 = gmyed1
                  call el2end( biglst, nodlst, elptrs(j), nxtel )
               end if
            end if
            j = j + 1
            if( j .le. numelm ) goto 450
            WRITE(0,*) '*** GMYLIN: REACHED END OF ELEMENT LIST (2)'
            WRITE(0,*) (elptrs(l),l=1,numelm)
            WRITE(0,*) (edptrs(l),l=1,numedg)
            WRITE(0,*) (ndptrs(l),l=1,numedg)
            WRITE(0,*) gmyed1,numedg,numelm,inod,ielm,k,minj,j
            ierr = -5558
            return
         end if
         call el2end( biglst, nodlst, elptrs(j), nxtel )
         if( gmyed1 .ne. gmyed2 ) then
            gmyed1 = gmyed2
            goto 450
         end if
      end if
c
      return
      end
c
c
c
      subroutine el2end( biglst, nodlst, ielm, nxtel )
c
      implicit none
c
      integer ielm, nxtel
c
      include 'blknew.i'
      include 'blkbig.i'
c
      integer lst, nxt, iptr2, iptr3
c
      lst   = biglst(lstbig,ielm)
      iptr2 = biglst(nxtbig,ielm)
      iptr3 = biglst(nxtbig,iptr2)
      nxt   = biglst(nxtbig,iptr3)
c
      if( nxtel .eq. ielm ) nxtel = nxt
      if( nxtel .eq. emtbig ) then
         nxtel = ielm
         if( nxt .eq. emtbig ) return
      end if
c
      if( ielm .eq. sttbig ) sttbig = nxt
c
      if( lst .gt. 0 ) biglst( nxtbig, lst ) = nxt
      if( nxt .gt. 0 ) biglst( lstbig, nxt ) = lst
c
      if( endbig .eq. iptr3 ) endbig = lst
c
      biglst( lstbig, ielm  ) = endbig
      biglst( nxtbig, iptr3 ) = emtbig
c
      if( emtbig .gt. 0 ) biglst( lstbig, emtbig ) = iptr3
      if( endbig .gt. 0 ) biglst( nxtbig, endbig ) = ielm
c
      endbig = iptr3
c
      return
      end
c
