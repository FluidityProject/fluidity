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
      SUBROUTINE MKADPT( BIGLST, NODLST,
     :                   CLCGMY, SZSNLS, SNLBAS, SNLIST, SURFID,
     :                   PRDNDS, NPRDND,
     :                   NODX, NODY, NODZ, ORGMTX, NNOD,
     :                   ELMREG, ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   NELM, NSELM, SZENLS, SZNELS, GEOM3D, SRFGMY,
     :                   PNOD, NHALO, SCATER, NGATH, GATHER )
C-----------------------------------------------------------------------
C
C - This forms the adaptive mesh data structures (currently kept in a common
C - block) given the fixed mesh information (in NODX/Y/Z, ENLIST, NELIST,
C - and EELIST etc.) and also the element material/region id (ELMREG).
C - The geometry of the adapted mesh will retain the geometry of these
C - element material/region boundaries. ORGMTX should contain the nine
C - entries for the metric of each node.
C - If the surface of the mesh should remain unchanged then set SRFGMY
C - to TRUE - otherwise the surface mesh will probably be changed.
C - The halo nodes are made untouchable, as are the halo elements, and
C - all edges of halo elements.
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INTEGER NELM, SZENLS, SZNELS, NNOD, NHALO, PNOD, NGATH,
     :        NSELM, SZSNLS, NPRDND
C
      INTEGER ELMSUR(NELM), NELBAS(NNOD+1), ELMREG(NELM),
     :        ENLIST(SZENLS), NELIST(SZNELS), ENLBAS(NELM+1),
     :        EELIST(NELM*4), SCATER(NHALO), SNLIST(SZSNLS),
     :        SNLBAS(NSELM+1), SURFID(NSELM), PRDNDS(NPRDND),
     :        GATHER(NGATH)
C
      REAL NODX(NNOD), NODY(NNOD), NODZ(NNOD), ORGMTX(9,NNOD)
C
      LOGICAL GEOM3D, SRFGMY, GTNDFL, CLCGMY
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      REAL XSHIFT, YSHIFT, ZSHIFT, XSCALE, YSCALE, ZSCALE
      COMMON / BLKSCL / XSHIFT, YSHIFT, ZSHIFT, XSCALE, YSCALE, ZSCALE
C
      LOGICAL BLDING
      COMMON / BLKBLD / BLDING
C
      INTEGER CONELM(4), EDGLST(6), TMPEDG(6), IBAS, IELM, NLN, ILN,
     :        IEDG, IPTR, INEW, I1, I2, I3, I4, NODS(4), TMPCON(4),
     :        NXT, LST, IFLAG, IREG, LOCEDG(2,6), IND1, IND2, JELM,
     :        I, J, K
c
      real getdtm, det, XMAX, YMAX, ZMAX
C
      LOGICAL FLGELM, tstdts
c
      tstdts = .true.
C
      BLDING = .TRUE.
C
      IS3DMS = GEOM3D
C
      NELEMS = 0
      NEDGES = 0
C
      STTBIG = 1
      ENDBIG = -1
C
      ewrite(2,*)
     :     'Setting up initial integer list pointers...',MAXBIG
C
      CALL INITLS( EMTBIG, BIGLST, SIZBIG, NXTBIG, LSTBIG, MAXBIG )
C
      ewrite(2,*) 
     :     'Setting up initial real list pointers...',MXNODS
C
      CALL INIRLS( EMTNOD, NODLST, SIZNOD, NXTNOD, LSTNOD, MXNODS )
C
      STTNOD = 1
      STFRND = 1
C
      ewrite(2,*) 
     :     'Creating node list...',debug
C
      XMAX = NODX(1)
      YMAX = NODY(1)
      ZMAX = NODZ(1)
      XSHIFT = XMAX
      YSHIFT = YMAX
      ZSHIFT = ZMAX
      DO ILN = 2, NNOD
         XSHIFT = MIN( XSHIFT, NODX(ILN) )
         YSHIFT = MIN( YSHIFT, NODY(ILN) )
         ZSHIFT = MIN( ZSHIFT, NODZ(ILN) )
         XMAX = MAX( XMAX, NODX(ILN) )
         YMAX = MAX( YMAX, NODY(ILN) )
         ZMAX = MAX( ZMAX, NODZ(ILN) )
      END DO
      XSCALE = (XMAX-XSHIFT)/500.0
      YSCALE = (YMAX-YSHIFT)/500.0
      ZSCALE = (ZMAX-ZSHIFT)/500.0
C
      DO ILN = 1, NNOD
C
C - node positions are shifted (to make corner at origin) and scaled
C - (to make domain dimensions 500 in all directions)
C
         NODLST( 1,ILN) = (NODX(ILN)-XSHIFT)/XSCALE
         NODLST( 2,ILN) = (NODY(ILN)-YSHIFT)/YSCALE
         NODLST( 3,ILN) = (NODZ(ILN)-ZSHIFT)/ZSCALE
         NODLST( 6,ILN) = 1.0
C
C - must scale metrics too...
C
         NODLST( 7,ILN) = ORGMTX(1,ILN)*XSCALE*XSCALE
         NODLST( 8,ILN) = ORGMTX(2,ILN)*XSCALE*YSCALE
         NODLST( 9,ILN) = ORGMTX(3,ILN)*XSCALE*ZSCALE
         NODLST(10,ILN) = ORGMTX(4,ILN)*YSCALE*XSCALE
         NODLST(11,ILN) = ORGMTX(5,ILN)*YSCALE*YSCALE
         NODLST(12,ILN) = ORGMTX(6,ILN)*YSCALE*ZSCALE
         NODLST(13,ILN) = ORGMTX(7,ILN)*ZSCALE*XSCALE
         NODLST(14,ILN) = ORGMTX(8,ILN)*ZSCALE*YSCALE
         NODLST(15,ILN) = ORGMTX(9,ILN)*ZSCALE*ZSCALE
C         DO J = 1, 9
C            NODLST(J+6,ILN) = ORGMTX(J,ILN)
C         END DO
         if( tstdts ) then
            det = getdtm( nodlst(7,iln) )
            if( nodlst(7,iln) .lt. 0.0 .or.
     :          nodlst(11,iln) .lt. 0.0 .or.
     :          nodlst(15,iln) .lt. 0.0 ) then
               WRITE(0,*) '+++ MKADPT: Warning, not +ve def: ',det,iln
               WRITE(0,*) (nodlst(j,iln),j=7,9)
               WRITE(0,*) (nodlst(j,iln),j=10,12)
               WRITE(0,*) (nodlst(j,iln),j=13,15)
            else if( det .lt. 1e-20 ) then
               WRITE(0,*) '+++ MKADPT: Warning, small det: ',det,iln
               WRITE(0,*) (nodlst(j,iln),j=7,9)
               WRITE(0,*) (nodlst(j,iln),j=10,12)
               WRITE(0,*) (nodlst(j,iln),j=13,15)
            end if
         end if
C
         NODLST(16,ILN) = FLOAT(NELIST(NELBAS(ILN)+1))
C
      END DO
C
      IF( NPRDND .GT. 1 ) THEN
C
         ewrite(2,*) 
     :        'Flagging periodic nodes... ',nprdnd
C
         DO I = 1, NPRDND
            J = PRDNDS(I)
            CALL STNDFL( BIGLST, NODLST, J, 7, .TRUE. )
c            print*,'  node: ',j,nodx(j),nody(j),nodz(j)
         END DO
C
      END IF
C
      IF( NNOD .GT. PNOD ) THEN
C
C - flag halo (receive) nodes (both flag bits 8 & 9)
C
         ewrite(2,*) 
     :        'Flagging halo nodes...'
C
         DO ILN = PNOD+1, NNOD
            CALL STNDFL( BIGLST, NODLST, ILN, 8, .TRUE. )
            CALL STNDFL( BIGLST, NODLST, ILN, 9, .TRUE. )
         END DO
C
C - unflag bit 9 for all nodes that are in receive list, so only nodes
C - left that can have flag bit 9 set (from above) are those which are
C - >PNOD, but which are not listed in the receive list.
C - Later on these nodes are removed from the adaptivity list of nodes,
C - and elements that have such nodes are not added to the list of
C - elements adaptivity will deal with (so these elements and nodes
C - will never be seen again after coming out of adaptivity...)
C - Really don't know the point of it all!
C
         DO ILN = 1, NHALO
            CALL STNDFL( BIGLST, NODLST, SCATER(ILN), 9, .FALSE. )
         END DO
C
C - flag send nodes (flag bit 10)
C
         DO ILN = 1, NGATH
            CALL STNDFL( BIGLST, NODLST, GATHER(ILN), 10, .TRUE. )
         END DO
C
      END IF
C
      TOPNOD = NNOD
C
      NFRTND = 0
      NUMNDS = NNOD
      ENDNOD = NNOD
      EMTNOD = NNOD+1
C
      IF( .NOT. GEOM3D ) THEN
C
         WRITE(0,*) 'Creating 2D element and edge list...'
         WRITE(0,*) 'Errr, well... eventually!'
         STOP
C
      ELSE
C
         ewrite(2,*)
     :        'Creating 3D element and edge list...'
C
         DO IELM = 1, NELM
            if( debug ) then
              ewrite(3, *)
     :           '>>> Adding element ',ielm
            end if
C
            IBAS = ENLBAS(IELM)
            NLN  = ENLBAS(IELM+1) - IBAS
C
            IF( NLN .NE. 4 ) THEN
               IERR = -12
               WRITE(0,*) '*** MKADPT: INCORRECT NODES IN 3D ELEMENTS'
               WRITE(0,*) 'ielm, nln: ',ielm,nln
               BLDING = .FALSE.
               RETURN
            END IF
C
            I1 = ENLIST(IBAS+1)
            I2 = ENLIST(IBAS+2)
            I3 = ENLIST(IBAS+3)
            I4 = ENLIST(IBAS+4)
C
            FLGELM = .FALSE.
C
C - for some reason, this ignores any element that includes these
C - 'weird' nodes (ones that are >PNOD, but not in receive list)
C
            IF( GTNDFL(BIGLST,NODLST,I1,9) ) THEN
               ELMREG(IELM) = 0
               GOTO 100
            ELSE IF( GTNDFL(BIGLST,NODLST,I2,9) ) THEN
               ELMREG(IELM) = 0
               GOTO 100
            ELSE IF( GTNDFL(BIGLST,NODLST,I3,9) ) THEN
               ELMREG(IELM) = 0
               GOTO 100
            ELSE IF( GTNDFL(BIGLST,NODLST,I4,9) ) THEN
               ELMREG(IELM) = 0
               GOTO 100
            END IF
C
C - 'fix' elements that have receive nodes, or have all send nodes (2nd level)
C
            IF( GTNDFL(BIGLST,NODLST,I1,8) ) THEN
               FLGELM = .TRUE.
            ELSE IF( GTNDFL(BIGLST,NODLST,I2,8) ) THEN
               FLGELM = .TRUE.
            ELSE IF( GTNDFL(BIGLST,NODLST,I3,8) ) THEN
               FLGELM = .TRUE.
            ELSE IF( GTNDFL(BIGLST,NODLST,I4,8) ) THEN
               FLGELM = .TRUE.
            ELSE IF( GTNDFL(BIGLST,NODLST,I1,10) ) THEN
               IF( GTNDFL(BIGLST,NODLST,I2,10) ) THEN
                  IF( GTNDFL(BIGLST,NODLST,I3,10) ) THEN
                     IF( GTNDFL(BIGLST,NODLST,I4,10) ) THEN
                        FLGELM = .TRUE.   ! also fix elems with all send nodes
                     END IF
                  END IF
               END IF
            END IF
C
            LOCEDG(1,1) = I1
            LOCEDG(2,1) = I2
C
            LOCEDG(1,2) = I2
            LOCEDG(2,2) = I3
C
            LOCEDG(1,3) = I3
            LOCEDG(2,3) = I1
C
            LOCEDG(1,4) = I2
            LOCEDG(2,4) = I4
C
            LOCEDG(1,5) = I3
            LOCEDG(2,5) = I4
C
            LOCEDG(1,6) = I1
            LOCEDG(2,6) = I4
C
            EDGLST(1) = 0
            EDGLST(2) = 0
            EDGLST(3) = 0
            EDGLST(4) = 0
            EDGLST(5) = 0
            EDGLST(6) = 0
C
            CONELM(1) = EELIST(IELM*4-3)
            CONELM(2) = EELIST(IELM*4-2)
            CONELM(3) = EELIST(IELM*4-1)
            CONELM(4) = EELIST(IELM*4)
            if( debug ) print*,'    nodes: ',i1,i2,i3,i4
            if( debug ) print*,'    connected elements: ',CONELM
C
            DO I = 1, 4
C
               IF( CONELM(I) .GT. IELM ) THEN
                  CONELM(I) = 0
               ELSE IF( CONELM(I) .GT. 0 ) THEN
                  CONELM(I) = ELMREG(CONELM(I))
               ELSE
                  CONELM(I) = CONELM(I) - 1
               END IF
C
            END DO
C
            DO I = 1, 6
C
               IND1 = LOCEDG(1,I)
               IND2 = LOCEDG(2,I)
C
               JELM = IELM
               if( debug ) print*,'   node 1: ',
     :                 nodx(ind1),nody(ind1),nodz(ind1)
               if( debug ) print*,'   node 2: ',
     :                 nodx(ind2),nody(ind2),nodz(ind2)
               CALL FNDEDG( IND1, IND2, ENLBAS, ENLIST, NELBAS, NELIST,
     :                      NNOD, NELM, SZENLS, SZNELS, JELM,
     :                      .FALSE., .TRUE. )
               if( debug ) print*,'   JELM,NODS: ',JELM,IND1,IND2
C
               IF( JELM .GT. 0 ) THEN
C
C - need to check if the element has actually been added to the list
C
   5              IF( ELMREG(JELM) .LE. 0 ) THEN
                     CALL FNDEDG( IND1, IND2, ENLBAS, ENLIST,
     :                      NELBAS, NELIST,
     :                      NNOD, NELM, SZENLS, SZNELS, JELM,
     :                      .TRUE., .TRUE. )
                     IF( JELM .GE. IELM ) THEN
                        JELM = 0
                     ELSE IF( JELM .GT. 0 ) THEN
                        GOTO 5
                     END IF
                  END IF
C
               END IF
C
               IF( JELM .GT. 0 ) THEN
C
                  JELM = ELMREG(JELM)
C
                  CALL ELMINF( BIGLST, NODLST, JELM, TMPCON, TMPEDG,
     :                         INEW, IFLAG, IREG, NXT, LST )
C
                  IF( IERR .NE. 0 ) THEN
                    WRITE(0,*)'*** MKADAPT: Got error from ELMINF ',ierr
                    RETURN
                  END IF
C
                  J = 0
  10              IF( J .LT. 6 ) THEN
                     J = J + 1
                     IF( IND1 .NE. BIGLST(1,TMPEDG(J)) .AND.
     :                   IND1 .NE. BIGLST(2,TMPEDG(J)) ) THEN
                        GOTO 10
                     ELSE IF( IND2 .NE. BIGLST(1,TMPEDG(J)) .AND.
     :                        IND2 .NE. BIGLST(2,TMPEDG(J)) ) THEN
                        GOTO 10
                     END IF
                     EDGLST(I) = TMPEDG(J)
c                     if( jelm .lt. 10 )
c     :                  print*,'   used edge (nds)',edglst(i),
c     :                             LOCEDG(1,I),LOCEDG(2,I)
                     if( debug )  print*,'   used edge (nds)',edglst(i),
     :                                       LOCEDG(1,I),LOCEDG(2,I)
                  ELSE
                     WRITE(0,*) '*** MKADPT: LOST COMMON EDGE'
                     WRITE(0,*) IELM,IND1,IND2,JELM
                     WRITE(0,*) TMPEDG
                     IERR = -1099
                     BLDING = .FALSE.
                     RETURN
                  END IF
C
               ELSE
C
                  INEW = -1
                  CALL ADDEDG( BIGLST, NODLST,
     :                         LOCEDG(1,I), EDGLST(I), INEW )
C
C - If both nodes are periodic then flag this edge so it cannot be changed
C
                  IF( GTNDFL( BIGLST, NODLST, IND1, 7 ) .AND.
     :                GTNDFL( BIGLST, NODLST, IND2, 7 ) ) THEN
                     CALL SETFLG( BIGLST, NODLST, EDGLST(I), 6, .TRUE. )
                     ewrite(3,*) 
     :                    '   Fixed edge ',edglst(i),ind1,ind2
                  END IF
C
                  IF( IERR .NE. 0 ) THEN
                     BLDING = .FALSE.
                     RETURN
                  END IF
C
c                  if( jelm .lt. 10 )
c     :                 print*,'   added edge (nds)',edglst(i),
c     :                            LOCEDG(1,I),LOCEDG(2,I)
                  if( debug ) print*,'   added edge (nds)',edglst(i),
     :                                   LOCEDG(1,I),LOCEDG(2,I)
C
               END IF
C
            END DO
C
            if( debug ) print*,'   about to add new element...',edglst
c            if( ielm.lt.10 )
c     :         print*,'   about to add new element...',edglst
            INEW = -1
            CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IPTR,
     :                   INEW, ELMREG(IELM) )
      if( ierr.eq.12345 ) then
         WRITE(0,*) '--- inconsistent oppedg from addele in mkadpt (1)'
         ierr = 0
      end if
C
            IF( IERR .NE. 0 ) THEN
               BLDING = .FALSE.
               RETURN
            END IF
            if( debug ) print*,'   added element ',iptr
c            if( jelm.lt.10 )
c     :         print*,'   added element ',iptr
C
            IF( FLGELM ) THEN
C
               CALL SETFLG( BIGLST, NODLST, IPTR, 5, .TRUE. )
               DO I = 1, 6
                  CALL SETFLG( BIGLST, NODLST, EDGLST(I), 6, .TRUE. )
               END DO
C
            END IF
C
            ELMREG(IELM)   = IPTR
C
  100       CONTINUE
C
         END DO
C
      END IF
C
      BLDING = .FALSE.
C
      TOPBIG = EMTBIG
C
      IF( NNOD .GT. PNOD ) THEN
C
C - unflag send nodes (just in case flag bit 10 is used elsewhere later)
C
         DO ILN = 1, NGATH
            CALL STNDFL( BIGLST, NODLST, GATHER(ILN), 10, .FALSE. )
         END DO
C
C - this appears to remove the 'weird' nodes (i.e. that are >PNOD,
C - but not in receive list) -not sure the point of it all...
C
         DO ILN = PNOD+1, NNOD
            IF( GTNDFL( BIGLST, NODLST, ILN, 9 ) ) THEN
               CALL RMVND( BIGLST, NODLST, ILN )
            END IF
         END DO
C
      END IF
C
      CALL SETGMY( BIGLST, NODLST, SRFGMY, CLCGMY,
     :             NSELM, SZSNLS, SNLBAS, SNLIST, SURFID )
C
      ewrite(2,*) ' '
      ewrite(2,*) 
     :     '++++++++++ Finished making mesh ++++++++++++'
      ewrite(2,*) 
     :     'STTBIG,EMTBIG,ENDBIG: ',STTBIG,EMTBIG,ENDBIG
      ewrite(2,*) 
     :     'NELEMS,NEDGES: ',NELEMS,NEDGES
      ewrite(2,*) ' '
      RETURN
      END
C
