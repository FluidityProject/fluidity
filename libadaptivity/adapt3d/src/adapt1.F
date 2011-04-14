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
      SUBROUTINE ADAPT1( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   SPLTED, CLPSED, ED2FC,  FC2ED,  SPLTEL, MVNODE,
     :                   FLGSET, FLGCHK, ALWAYS, DOTOP,  MINCHG, CHANGE)
C-----------------------------------------------------------------------
C
C - This performs one adaptive sweep of the mesh performing the desired
C - operations as follows:
C
C           SPLTED : split edges if true
C           CLPSED : collapse edges if true
C           EDSWAP : perform edge to face and edge to edge swapping if true
C           FC2ED  : perform face to edge swapping if true
C           SPLTEL : split elements (do not use this yet)
C           MVNODE : move nodes if true
C
C - The FLGSET and FLGCHK allow for selectivity in elements to be checked
C - for adaptivity. These integers can be from 14 to 29 inclusive. ADAPT1
C - will set flag FLGSET of an element before attempting to adapt. Hence,
C - when ADAPT1 returns, any elements that remain with flag FLGSET set have
C - not been changed. However, elements with flag FLGCHK set will not be
C - checked for adaptivity (although they may be changed due to a neighbouring
C - element being changed). So, if FLGSET and FLGCHK are alternated between
C - values from one call of ADAPT1 to another, it allows redundant checking
C - to be minimised.
C
C - ALWAYS is the amount of functional change to always perform as soon as
C - it is found (i.e. ignore other checking for that element). It should
C - be negative (eg. -1.0 or -2.0). DOTOP is the minimum element functional
C - to allow to be changed (eg. 0.15). MINCHG is the minimum relative change
C - allowed to a local functional (eg. 0.01). CHANGE is returned as the
C - sum of local functional changes - note: it does not represent the total
C - change in the global functional.
C
C - It is advised that ADAPT1 should be called first with just node movement
C - then with everything but node movement (SPLTEL is not yet functional and
C - is excluded from this advice). The flags (FLGSET and FLGCHK) should be
C - independent for each of the two sets of calls, eg. interchange 14 & 15
C - for the node movement call, and interchange 16 & 17 for the other.
C
C-----------------------------------------------------------------------
      use write_log
      use AdaptProgress
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkorg.i'
C
      INCLUDE 'blknxt.i'
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      REAL ALWAYS, CHANGE, DOTOP, MINCHG
C
      LOGICAL SPLTED, CLPSED, MVNODE, FC2ED, ED2FC, SPLTEL
C
      INTEGER IPTR, CONELM(4), EDGLST(6), INEW, IFLAG, NXT, LST,
     :        NDPTRS(MAXEDE+1), ELPTRS(3,MAXEDE), EDPTRS(3,MAXEDE+1),
     :        NDPTS2(2*MAXEDE), ELPTS2(MAXEDE*2), EDPTS2(2*MAXEDE),
     :        BSTPTR, BSTCEL, BSTWHT, NUMEDE, SURFAC, NEWEL1, NEWEL2,
     :        NEWED1, NEWED2, IWHAT, I, MYNXT, NXTID, NODS(4), LASTPT,
     :        FLGSET, FLGCHK, IREG, BSTREG, INTRNL, BSTINE, STDOUT,
     :        INELM, AWHAT, INELMA, cnd2vl, trace, iumnds, ednd1, ednd2,
     :        gothere, numel1, numel2, numed1, numed2, trcnd
C
      REAL BSTDIF, CURDIF, CURORG, BSTCUR, BSTORG, RATIO, BSTRAT,
     :     BSTX, BSTY, BSTZ, X, Y, Z, XA, YA, ZA,
     :     GF, SUMDIF, SUMORG, BSTSUM, RATIA,
     :     edgfnc, elmfnc, getdtm, det,
     :     a, b, c, d, e, f, h, s, t, u
C
      LOGICAL GETFLG, DUNONE, GTNDFL, CLPS2, SPLED2,
     :        shwchg, tstdts, goterr
C
      CHARACTER*7 STR1, STR2

c      DEBUG = .FALSE.
      shwchg = .false.
      tstdts = .false.
c      if( anisot .gt. 1 )
c      shwchg = .true.
      goterr = debug
c      goterr = .true.
C
      ADPTNG = .FALSE.
      mvcnt  = 0
      trace  = 0
      trcnd  = 0
c      trace  = 6116063
c      trcnd  = 2
      SPLED2 = SPLTED
C
      CALL INITIM
C
      CALL ADAPT_SECOND( TIMSTT )
C
      a = 0.0
      b = 0.0
      c = 1e+30
      d = 0.0
      e = 0.0
      f = 1e+30
      h = 0.0
      s = 0.0
      t = 0.0
      u = 0.0
C
      GF = 0.0
c
      if( goterr ) then
         ewrite(3,6) 'Into ADAPT1... ',ierr
         ewrite(3,7) '   always,dotop: ',always,dotop
      end if
C
      IF( IERR .NE. 0 ) RETURN
C
      EDGON = .TRUE.
C
      ELSLFT = 0
      TOPBIG = 0
C
C - make sure the appropriate flags are clear before sweeping the elements
C
      IPTR = STTBIG
C
 10   IF( BIGLST(4,IPTR) .LT. 0 ) THEN
c         ewrite(3,6)
c     :           'Clearing flags 4 & 8 for edge ',iptr
         CALL SETFLG( BIGLST, NODLST, IPTR, 4, .FALSE. )
         CALL SETFLG( BIGLST, NODLST, IPTR, 8, .FALSE. )
         inew  = biglst(1,iptr)
         iflag = biglst(2,iptr)
c         ewrite(3,8) 'Getting edge functional...'
         h = edgfnc(BIGLST,NODLST,inew,iflag,ratio)
c         ewrite(3,9)
c     :           'functional,length: ',h,ratio
c         a = max( a, h )
c         s = s + h
         b = max( abs(ratio), b )
         c = min( abs(ratio), c )
         TOPBIG = MAX(TOPBIG,IPTR)
         NXT = BIGLST(NXTBIG,IPTR)
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
c         ewrite(3,6)
c     :           'Setting flag 4 for element ',iptr
         edgon = .false.
         biglst(7,iptr) = 0
         nxt = biglst(nxtbig,iptr)
         biglst(7,nxt) = 0
         CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, .false. )
c         ewrite(3,6)
c     :           'Nodes: ',nods(1),nods(2),nods(3),nods(4)
         H = ELMFNC( BIGLST, NODLST, iptr,
     :               NODS(1), NODS(2), NODS(3), NODS(4), RATIO )
c         ewrite(3,9)
c     :           'functional,in-sphere: ',h,ratio
         GF = MAX( GF, H )
         u = u + h
         e = max( abs(ratio), e )
         f = min( abs(ratio), f )
         edgon = .true.
c         ewrite(3,'(A,I5,A)')
c     :          'Setting flag ',flgset,'...'
         IF( FLGSET .GT. 9 ) CALL SETFLG(BIGLST,NODLST,IPTR,
     :                                           FLGSET,.TRUE.)
c         ewrite(3,8) 'Clearing flags 7 and 3...'
         CALL SETFLG( BIGLST, NODLST, IPTR, 7, .FALSE. )
         CALL SETFLG( BIGLST, NODLST, IPTR, 3, .FALSE. )
         if( flgchk .gt. 9 ) then
            if( .not. getflg(BIGLST,NODLST,iptr,flgchk) ) then
               elslft = elslft + 1
            else
c               ewrite(3,8) 'Setting flag 7...'
               call setflg(BIGLST,NODLST,iptr,7,.true.)
            end if
         else
            ELSLFT = ELSLFT + 1
         end if
c         ewrite(3,8) 'finding max bounds...'
         TOPBIG = MAX( TOPBIG, IPTR )
         IPTR   = BIGLST( NXTBIG, IPTR )
         TOPBIG = MAX( TOPBIG, IPTR )
         IPTR   = BIGLST( NXTBIG, IPTR )
         TOPBIG = MAX( TOPBIG, IPTR )
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 10
C
c      IF( ELSLFT .NE. NELEMS ) THEN
c         PRINT*,'WARNING: NUMBER OF ELEMENTS DIFFERENT: ',
c     :           NELEMS,ELSLFT
c      END IF
C
      LSTLFT = (ELSLFT/100)*100
      NELADD = 0
      NELSUB = 0
      NELCHK = 0
c      lstlft = elslft
C
      ewrite(1,9) '>>> Min/max edges: ',c,b
      ewrite(1,9) '>>> Min/max in-spheres: ',f,e
c      print*,'>>> Edge / in-sphere functional sums: ',s,t
      ewrite(1,9)
     :       '>>> GLOBAL MESH FUNCTIONAL & element ave: ',GF,u/nelems
      gf = 1.0
c      pause
C
      EDGON  = .FALSE.
      CHANGE = 0.0
c
      SUMORG = 0.0
      CURORG = 0.0
      TOPNOD = 0
C
      IPTR = STTNOD
C
 15   CALL STNDFL( BIGLST, NODLST, IPTR, 6, .FALSE. )
      CALL STNDFL( BIGLST, NODLST, IPTR, 7, .FALSE. )
      CALL STNDFL( BIGLST, NODLST, IPTR, 9, .FALSE. )
      NXT = INT(NODLST(NXTNOD,IPTR))
C
      if( tstdts ) then
        det = getdtm( nodlst(7,iptr) )
        if( det .lt. 1e-20 .or. goterr ) then
          if( det .lt. 1e-20 ) then
             ewrite(0,2)
     :           '+++ ADAPT1: Got small determinant: ',iptr,det
          else
             ewrite(0,2) 'Node, determinant: ',iptr,det
          end if
          ewrite(0,'(5X,3(1PE13.5))')
     :      nodlst(7,iptr),nodlst(10,iptr),nodlst(13,iptr),
     :      nodlst(8,iptr),nodlst(11,iptr),nodlst(14,iptr),
     :      nodlst(9,iptr),nodlst(12,iptr),nodlst(15,iptr)
        end if
      end if
c
c      x = nodlst(1,iptr)
c      y = nodlst(2,iptr)
c      z = nodlst(3,iptr)
c
      TOPNOD = MAX( TOPNOD, IPTR )
      IPTR = NXT
C
      IF( GTNDFL(BIGLST,NODLST,IPTR,1) ) GOTO 15
c
c      print*,'   done clearing flags'
C
      ewrite(1,6) 
     :     'Top of BIGLST & NODLST:',TOPBIG,TOPNOD
C
      IF( IERR .NE. 0 ) RETURN
c
      if( .not. mvnode ) then
         ewrite(1,8) 
     :        '--- ADAPT1: Node movement OFF...'
      else if( minchg .lt. 0.0 ) then
         ewrite(1,8)
     :        '--- ADAPT1: Node movement in OVERSHOOT mode...'
      else
         ewrite(1,8) 
     :        '--- ADAPT1: Node movement in NORMAL mode...'
      end if
C
c      ELCTMN = (NELEMS/500)*500 - 500
c      ELCTMX = ELCTMN + 1000
C
c      EDCTMN = (NEDGES/500)*500 - 500
c      EDCTMX = EDCTMN + 1000
C
      ADPTNG = .TRUE.
C
C - start at the beginning...
C
      IPTR = STTBIG
C
      MYNXT = NXTID( IPTR )
c      print*,'   Got MYNXT: ',mynxt
C
 17   BSTPTR = 0
      BSTCEL = 0
      BSTWHT = 0
      BSTCUR = 0.0
      BSTDIF = 0.0
      BSTSUM = 0.0
      BSTORG = 0.0
      BSTINE = 0
C
      TELCHK = 0
      TEDCHK = 0
      TNDCHK = 0
C
      DUNONE = .FALSE.
      gothere = 0
C
C - ...and run through BIGLST picking out the elements and finding out
C - what we can do to each one (face for edge swap) or its edges
C - (edge to face swap/edge splitting/edge removal)
C
 20   LASTPT = IPTR
      if(mod(telchk+1, 100000).eq.0) then
         if(should_exit(NUMNDS)) then
            goto 900
         end if
      end if
C
      IF( SPLED2 ) THEN
         IF( NELEMS*3+NEDGES .GT. (MAXBIG*4/5) ) THEN
            ewrite(0,'(A,A)') 
     :           '+++ ADAPT1: WARNING! ',
     :          'BIGLST FILLING UP -Switched off edge splitting'
            SPLED2 = .FALSE.
         ELSE IF( NUMNDS .GT. (MXNODS*9/10) ) THEN
            ewrite(0,'(A,A)') 
     :           '+++ ADAPT1: WARNING! ',
     :          'NODLST FILLING UP -Switched off edge splitting'
            SPLED2 = .FALSE.
         END IF
      ELSE IF( SPLTED ) THEN
         IF( NELEMS*3+NEDGES .LT. (MAXBIG*15/20) .AND.
     :       NUMNDS .LT. (MXNODS*8/10) ) THEN
            ewrite(0,8) 
     :   '--- ADAPT1: BIGLST/NODLST now clear enough for edge split'
            SPLED2 = .TRUE.
         END IF
      END IF
C
      IF( BIGLST(4,IPTR) .LT. 0 ) THEN
C
         BSTCUR = 0.0
         BSTSUM = 0.0
         BSTWHT = 0
         NXT = BIGLST(NXTBIG,IPTR)
C
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
C
cc         if( getflg(BIGLST,NODLST,iptr,7) ) then
cc            print*,'================= ALREADY DONE THIS ELEMENT: ',IPTR
cc         end if
c
         if( trace .gt. 0 ) cnd2vl = biglst(trcnd,trace)
c
         BSTCUR = 0.0
         BSTSUM = 0.0
         BSTWHT = 0
C
         CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
C
         NEXTLS( MYNXT ) = NXT
c         print*,'   Set NEXTLS ',NXT,MYNXT
C
         IF( GETFLG(BIGLST,NODLST,IPTR,5) ) THEN
            if( goterr ) then
              ewrite(3,6)
     :               'Ignoring element (5)',iptr,nxt
            end if
            GOTO 600
         END IF
C
         IF( FLGCHK .GT. 9 ) THEN
            IF( GETFLG(BIGLST,NODLST,IPTR,FLGCHK) ) THEN
               if( goterr ) then
                 ewrite(3,6)
     :               'Ignoring element (chk)',iptr,nxt
               end if
               GOTO 600
            END IF
         END IF
c
         if( goterr ) then
           ewrite(3,6) 'On element ',iptr,nxt
         end if
c         print*,'   conelm: ',conelm
c         print*,'   edglst: ',edglst
C
         IF( IERR .NE. 0 ) THEN
            ewrite(3,8) 
     :           '*** ADAPT1: GOT ERROR FROM ELMINF/GETFLG'
           ewrite(3,6) 
     :           '   element: ',iptr
           RETURN
         END IF
C
ccc         NEXTLS( MYNXT ) = NXT
cccc         print*,'   Set NEXTLS ',NXT,MYNXT
C
         IF( FC2ED ) THEN
C
C - check all pairs with this element that have not already been checked
C
          DO I = 1, 4
C
           IF( CONELM(I) .GT. 0 ) THEN
             IF( .NOT. GETFLG(BIGLST,NODLST,CONELM(I),3) ) THEN
             IF( .NOT. GETFLG(BIGLST,NODLST,CONELM(I),5) ) THEN
c
               if( ierr.ne.0 ) then
                  ewrite(-1,8)
     :                   '*** ADAPT1: GOT ERROR FROM GETFLG 1'
                  ewrite(-1,6)  '   element: ',iptr
                  ewrite(-1,6) '   conelm: ',conelm
                  ewrite(-1,6) '   edglst: ',edglst
                  return
               end if
C
               if( goterr ) then
                 ewrite(3,6)
     :              '   Checking with element ',conelm(i)
               end if
               CALL ELSTST( BIGLST, NODLST, IPTR, IREG, CONELM, EDGLST,
     :                      CURDIF, CURORG, SUMDIF, SUMORG, I,
     :                      GF*DOTOP, ABS(MINCHG) )
               if( goterr ) then
                 ewrite(3,8) '   Done.'
               end if
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :          '*** CHANGED AFTER ELSTST: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
               IF( IERR .NE. 0 ) THEN
                  ewrite(-1,8)
     :                 '*** ADAPT1: GOT ERROR FROM ELSTST'
                  ewrite(-1,6) '   element: ',iptr
                  ewrite(-1,6) '   conelm: ',conelm
                  ewrite(-1,6) '   edglst: ',edglst
                  return
               END IF
C
               IF( CURORG .GT. GF*DOTOP ) THEN
                 IF( CURDIF .LT. -ABS(MINCHG) ) THEN
                   IF( CURDIF .LT. BSTCUR ) THEN
                     BSTPTR = IPTR
                     BSTCEL = CONELM(I)
                     BSTWHT = -1
                     BSTCUR = CURDIF
                     BSTORG = CURORG
                     DUNONE = .TRUE.
                     if( curdif .lt. ALWAYS ) goto 590
                     if( curorg+curdif .lt. gf*dotop ) goto 590
                   END IF
                 ELSE IF( BSTCUR .EQ. 0.0 .AND.
     :                    SUMDIF .LT. -ABS(MINCHG) ) THEN
                   IF( SUMDIF .LT. BSTSUM ) THEN
                     BSTPTR = IPTR
                     BSTCEL = CONELM(I)
                     BSTWHT = -1
                     BSTSUM = SUMDIF
                     BSTORG = SUMORG
                     DUNONE = .TRUE.
                     if( sumdif .lt. ALWAYS ) goto 590
                   END IF
c                 ELSE IF( CURDIF .LT. BSTDIF .AND. .NOT. DUNONE ) THEN
c                   BSTPTR = IPTR
c                   BSTCEL = CONELM(I)
c                   BSTWHT = -1
c                   BSTDIF = CURDIF
                 END IF
               END IF
C
             else if( goterr ) then
                ewrite(3,6)
     :                 '   Skipping pair check (5) ',conelm(i),i
             END IF
             else if( goterr ) then
                ewrite(3,6)
     :                 '   Skipping pair check (3) ',conelm(i),i
             END IF
           else if( goterr ) then
              ewrite(3,6)
     :                 '   Skipping check with surface ',conelm(i),i
           END IF
C
          END DO
C
         END IF
C
         IF( ED2FC .OR. SPLED2 .OR. CLPSED ) THEN
C
C - now check the edges of this element that have not been checked
C
          DO I = 1, 6
C
           IF( .NOT. GETFLG(BIGLST,NODLST,EDGLST(I),4) ) THEN
           IF( .NOT. GETFLG(BIGLST,NODLST,EDGLST(I),6) ) THEN
c
              if( ierr.ne.0 ) then
                ewrite(-1,8) 
     :                '*** ADAPT1: GOT ERROR FROM GETFLG 2'
                ewrite(-1,6) '   element: ',iptr
                ewrite(-1,6) '   conelm: ',conelm
                ewrite(-1,6) '   edglst: ',edglst
                return
              end if
c
              if(goterr) then
                ewrite(3,6)
     :             '   checking edge ', edglst(i),biglst(4,edglst(i)),
     :             biglst(1,edglst(i)),biglst(2,edglst(i)),i
              end if
C
             TEDCHK = TEDCHK + 1
C
             CALL EDGELS( BIGLST, NODLST, EDGLST(I), IPTR, ELPTRS,
     :                    EDPTRS, NDPTRS, NUMEDE, SURFAC, INTRNL )
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :           '*** CHANGED AFTER EDGELS: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
             if(goterr) then
               ewrite(3,'(A)') '   Done EDGELS.'
             end if
             IF( IERR .NE. 0 ) THEN
                ewrite(-1,8) 
     :               '*** ADAPT1: GOT ERROR FROM EDGELS'
                ewrite(-1,6)
     :                 '   elem,edge,srf: ',iptr,edglst(i),surfac
                ewrite(-1,6) '   conelm: ',conelm
                ewrite(-1,6) '   edglst: ',edglst
                RETURN
             END IF
C
             ratio  = 0.5
             numel1 = 0
             numel2 = 0
             numed1 = 0
             numed2 = 0
C
             iumnds = numnds
             CALL EDGTST( BIGLST, NODLST, ENLBAS, ENLIST,
     :                    NELBAS, NELIST, EELIST,
     :                    SZENLS, SZNELS, NNOD, NELM,
     :                    NODX, NODY, NODZ, ORGMTX,
     :                    SPLED2, CLPSED, ED2FC, MVNODE,
     :                    EDGLST(I), IPTR, ELPTRS, EDPTRS, NDPTRS,
     :                    NUMEDE, SURFAC, INTRNL,
     :                    ELPTS2, EDPTS2, NDPTS2,
     :                    NUMEL1, NUMEL2, NUMED1, NUMED2,
     :                    CURDIF, CURORG, IWHAT,
     :                    SUMDIF, SUMORG, AWHAT,
     :                    RATIO, X,  Y,  Z,  INELM,
     :                    RATIA, XA, YA, ZA, INELMA,
     :                    ALWAYS, GF*DOTOP, ABS(MINCHG) )
             if(goterr) then
               ewrite(3,9) 
     :            '   Done EDGTST:',curdif,sumdif
             end if
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(-1,6)
     :           '*** CHANGED AFTER EDGTST: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
             if(iumnds.ne.numnds) then
                ewrite(-1,8)
     :               '+++ ADAPT1: NUMNDS CHANGED AFTER EDGTST!!!'
                ewrite(-1,6)
     :               '    iumnds,numnds,elem: ',iumnds,numnds,iptr
                ewrite(-1,*)
     :               '    iwhat,awhat,rat: ',iwhat,awhat,ratio
                ewrite(-1,'(A,3(1PE13.5))')
     :               '    x,y,z: ',x,y,z
                ednd1 = biglst(1,edglst(i))
                ednd2 = biglst(2,edglst(i))
                ewrite(-1,'(A,3(1PE13.5),I8)') 
     :               '    nod1: ',
     :         (nodlst(iumnds,ednd1),iumnds=1,3), int(nodlst(6,ednd1))
                ewrite(-1,'(A,3(1PE13.5),I8)') 
     :               '    nod2: ',
     :           (nodlst(iumnds,ednd2),iumnds=1,3), int(nodlst(6,ednd2))
             endif
c
             IF( IERR .NE. 0 ) THEN
               ewrite(-1,8) 
     :               '*** ADAPT1: GOT ERROR FROM EDGTST 1'
               ewrite(-1,6) '   element,edge: ',iptr,edglst(i)
               ewrite(-1,6) '   conelm: ',conelm
               ewrite(-1,6) '   edglst: ',edglst
               RETURN
             END IF
c
C             if( edglst(i) .eq. 30926 ) then
C                print*,'IN ADAPT1...'
C                print*,'   curdif,curorg: ',curdif,curorg
C                print*,'   always,bstdif,dunone: ',always,bstdif,dunone
C             end if
C
             IF( CURORG .GT. GF*DOTOP ) THEN
               IF( CURDIF .LT. -ABS(MINCHG) ) THEN
                  IF( CURDIF .LT. BSTCUR ) THEN
                    BSTPTR = EDGLST(I)
                    BSTCEL = IPTR
                    BSTWHT = IWHAT
                    BSTCUR = CURDIF
                    BSTORG = CURORG
                    BSTRAT = RATIO
                    DUNONE = .TRUE.
                    BSTREG = INTRNL
                    BSTX   = X
                    BSTY   = Y
                    BSTZ   = Z
                    BSTINE = INELM
                    if( curdif .lt. ALWAYS ) goto 590
                    if( curorg+curdif .lt. gf*dotop ) goto 590
                  END IF
               ELSE IF( BSTCUR .EQ. 0.0 .AND.
     :                  SUMDIF .LT. -ABS(MINCHG) ) THEN
                  IF( SUMDIF .LT. BSTSUM ) THEN
                    BSTPTR = EDGLST(I)
                    BSTCEL = IPTR
                    BSTWHT = AWHAT
                    BSTSUM = SUMDIF
                    BSTORG = SUMORG
                    BSTRAT = RATIA
                    DUNONE = .TRUE.
                    BSTREG = INTRNL
                    BSTX   = XA
                    BSTY   = YA
                    BSTZ   = ZA
                    BSTINE = INELMA
                    if( sumdif .lt. ALWAYS ) goto 590
                  END IF
c               ELSE IF( CURDIF .LT. BSTDIF .AND. .NOT. DUNONE ) THEN
c                 BSTPTR = EDGLST(I)
c                 BSTCEL = IPTR
c                 BSTWHT = IWHAT
c                 BSTDIF = CURDIF
               END IF
             END IF
             goto 987
C
             CLPS2 = CLPSED
             IF( RATIO .NE. 0.5 ) CLPS2 = .FALSE.
C
             ratio = 0.1
C
             CALL EDGTST( BIGLST, NODLST, ENLBAS, ENLIST,
     :                    NELBAS, NELIST, EELIST,
     :                    SZENLS, SZNELS, NNOD, NELM,
     :                    NODX, NODY, NODZ, ORGMTX,
     :                    SPLED2, CLPS2, .FALSE., MVNODE,
     :                    EDGLST(I), IPTR, ELPTRS, EDPTRS, NDPTRS,
     :                    NUMEDE, SURFAC, INTRNL,
     :                    ELPTS2, EDPTS2, NDPTS2,
     :                    NUMEL1, NUMEL2, NUMED1, NUMED2,
     :                    CURDIF, CURORG, IWHAT,
     :                    SUMDIF, SUMORG, AWHAT,
     :                    RATIO, X,  Y,  Z,  INELM,
     :                    RATIA, XA, YA, ZA, INELMA,
     :                    ALWAYS, GF*DOTOP, ABS(MINCHG) )
             if( goterr ) then
               ewrite(3,9) 
     :            '   Done EDGTST:',curdif,sumdif
             end if
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :       '*** CHANGED AFTER 2nd EDGTST: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
c
             IF( IERR .NE. 0 ) THEN
               ewrite(-1,8) 
     :               '*** ADAPT1: GOT ERROR FROM EDGTST 2'
               ewrite(-1,6) 
     :              '   element,edge: ',iptr,edglst(i)
               ewrite(-1,6) 
     :              '   conelm: ',conelm
               ewrite(-1,6) 
     :              '   edglst: ',edglst
               RETURN
             END IF
c
C             if( edglst(i) .eq. 30926 ) then
C                print*,'IN ADAPT1...'
C                print*,'   curdif,curorg: ',curdif,curorg
C                print*,'   always,bstdif,dunone: ',always,bstdif,dunone
C             end if
C
             IF( CURORG .GT. GF*DOTOP ) THEN
               IF( CURDIF .LT. -ABS(MINCHG) ) THEN
                  IF( CURDIF .LT. BSTCUR ) THEN
                    BSTPTR = EDGLST(I)
                    BSTCEL = IPTR
                    BSTWHT = IWHAT
                    BSTCUR = CURDIF
                    BSTORG = CURORG
                    BSTRAT = RATIO
                    DUNONE = .TRUE.
                    BSTREG = INTRNL
                    BSTX   = X
                    BSTY   = Y
                    BSTZ   = Z
                    BSTINE = INELM
                    if( curdif .lt. ALWAYS ) goto 590
                    if( curorg+curdif .lt. gf*dotop ) goto 590
                  END IF
               ELSE IF( BSTCUR .EQ. 0.0 .AND.
     :                  SUMDIF .LT. -ABS(MINCHG) ) THEN
                  IF( SUMDIF .LT. BSTSUM ) THEN
                    BSTPTR = EDGLST(I)
                    BSTCEL = IPTR
                    BSTWHT = AWHAT
                    BSTSUM = SUMDIF
                    BSTORG = SUMORG
                    BSTRAT = RATIA
                    DUNONE = .TRUE.
                    BSTREG = INTRNL
                    BSTX   = XA
                    BSTY   = YA
                    BSTZ   = ZA
                    BSTINE = INELMA
                    if( sumdif .lt. ALWAYS ) goto 590
                  END IF
c               ELSE IF( CURDIF .LT. BSTDIF .AND. .NOT. DUNONE ) THEN
c                 BSTPTR = EDGLST(I)
c                 BSTCEL = IPTR
c                 BSTWHT = IWHAT
c                 BSTDIF = CURDIF
               END IF
             END IF
C
             ratio = 0.9
C
             CALL EDGTST( BIGLST, NODLST, ENLBAS, ENLIST,
     :                    NELBAS, NELIST, EELIST,
     :                    SZENLS, SZNELS, NNOD, NELM,
     :                    NODX, NODY, NODZ, ORGMTX,
     :                    SPLED2, CLPS2, .FALSE., MVNODE,
     :                    EDGLST(I), IPTR, ELPTRS, EDPTRS, NDPTRS,
     :                    NUMEDE, SURFAC, INTRNL,
     :                    ELPTS2, EDPTS2, NDPTS2,
     :                    NUMEL1, NUMEL2, NUMED1, NUMED2,
     :                    CURDIF, CURORG, IWHAT,
     :                    SUMDIF, SUMORG, AWHAT,
     :                    RATIO, X,  Y,  Z,  INELM,
     :                    RATIA, XA, YA, ZA, INELMA,
     :                    ALWAYS, GF*DOTOP, ABS(MINCHG) )
             if( goterr ) then
               ewrite(3,9) 
     :            '   Done EDGTST:',curdif,sumdif
             end if
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :       '*** CHANGED AFTER 3rd EDGTST: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
c
             IF( IERR .NE. 0 ) THEN
               ewrite(-1,8) 
     :               '*** ADAPT1: GOT ERROR FROM EDGTST 3'
               ewrite(-1,6) 
     :              '   element,edge: ',iptr,edglst(i)
               ewrite(-1,6) 
     :              '   conelm: ',conelm
               ewrite(-1,6) 
     :              '   edglst: ',edglst
               RETURN
             END IF
c
C             if( edglst(i) .eq. 30926 ) then
C                print*,'IN ADAPT1...'
C                print*,'   curdif,curorg: ',curdif,curorg
C                print*,'   always,bstdif,dunone: ',always,bstdif,dunone
C             end if
C
             IF( CURORG .GT. GF*DOTOP ) THEN
               IF( CURDIF .LT. -ABS(MINCHG) ) THEN
                  IF( CURDIF .LT. BSTCUR ) THEN
                    BSTPTR = EDGLST(I)
                    BSTCEL = IPTR
                    BSTWHT = IWHAT
                    BSTCUR = CURDIF
                    BSTORG = CURORG
                    BSTRAT = RATIO
                    DUNONE = .TRUE.
                    BSTREG = INTRNL
                    BSTX   = X
                    BSTY   = Y
                    BSTZ   = Z
                    BSTINE = INELM
                    if( curdif .lt. ALWAYS ) goto 590
                    if( curorg+curdif .lt. gf*dotop ) goto 590
                  END IF
               ELSE IF( BSTCUR .EQ. 0.0 .AND.
     :                  SUMDIF .LT. -ABS(MINCHG) ) THEN
                  IF( SUMDIF .LT. BSTSUM ) THEN
                    BSTPTR = EDGLST(I)
                    BSTCEL = IPTR
                    BSTWHT = AWHAT
                    BSTSUM = SUMDIF
                    BSTORG = SUMORG
                    BSTRAT = RATIA
                    DUNONE = .TRUE.
                    BSTREG = INTRNL
                    BSTX   = XA
                    BSTY   = YA
                    BSTZ   = ZA
                    BSTINE = INELMA
                    if( sumdif .lt. ALWAYS ) goto 590
                  END IF
c               ELSE IF( CURDIF .LT. BSTDIF .AND. .NOT. DUNONE ) THEN
c                 BSTPTR = EDGLST(I)
c                 BSTCEL = IPTR
c                 BSTWHT = IWHAT
c                 BSTDIF = CURDIF
               END IF
             END IF
C
C - flag that we have checked this edge
C
  987        IF( SPLTED .EQV. SPLED2 )
     :          CALL SETFLG( BIGLST, NODLST, EDGLST(I), 4, .TRUE. )
C
           else if( goterr ) then
             ewrite(3,6)
     :              '   Skipping check of edge (6) ',edglst(i),i
           END IF
           else if( goterr ) then
             ewrite(3,6)
     :              '   Skipping check of edge (4) ',edglst(i),i
           END IF
C
          END DO
C
         END IF
C
         IF( SPLTEL ) THEN
C
            write(0,8) 'SHOULD NOT USE SPLTEL AT THE MOMENT!!!'
            stop
c
          CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, .false. )
C
          IF( IERR .NE. 0 ) THEN
             ewrite(-1,8) 
     :            '*** ADAPT1: GOT ERROR FROM ELNODS (1)'
             ewrite(-1,6) 
     :            '   element: ',iptr
          END IF
C
          CALL ELMTST( BIGLST, NODLST,
     :                 ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                 SZENLS, SZNELS, NNOD,   NELM,
     :                 NODX,   NODY,   NODZ,   ORGMTX,
     :                 IPTR,   NODS,   CURDIF, CURORG, SUMDIF, SUMORG )
C
          IF( CURORG .GT. GF*DOTOP ) THEN
            IF( CURDIF .LT. -ABS(MINCHG) ) THEN
              IF( CURDIF .LT. BSTCUR ) THEN
                BSTPTR = IPTR
                BSTCEL = IPTR
                BSTWHT = -2
                BSTCUR = CURDIF
                BSTORG = CURORG
                DUNONE = .TRUE.
                if( curdif .lt. ALWAYS ) goto 590
                if( curorg+curdif .lt. gf*dotop ) goto 590
              END IF
           ELSE IF( BSTCUR .EQ. 0.0 .AND.
     :              SUMDIF .LT. -ABS(MINCHG) ) THEN
              IF( SUMDIF .LT. BSTSUM ) THEN
                BSTPTR = IPTR
                BSTCEL = IPTR
                BSTWHT = -2
                BSTSUM = SUMDIF
                BSTORG = SUMORG
                DUNONE = .TRUE.
                if( sumdif .lt. ALWAYS ) goto 590
              END IF
c           ELSE IF( CURDIF .LT. BSTDIF .AND. .NOT. DUNONE ) THEN
c             BSTPTR = IPTR
c             BSTCEL = IPTR
c             BSTWHT = -2
c             BSTDIF = CURDIF
            END IF
          END IF
C
         END IF
C
         IF( MVNODE .and. bstwht .eq. 0 ) THEN
C
          IF( .NOT. SPLTEL ) THEN
C
             CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST,.false. )
C
             IF( IERR .NE. 0 ) THEN
              ewrite(-1,8) 
     :               '*** ADAPT1: GOT ERROR FROM ELNODS (2)'
              ewrite(-1,6) 
     :             '   element: ',iptr
             END IF
C
          END IF
C
          DO I = 1, 4
C
            IF( .NOT. GTNDFL(BIGLST,NODLST,NODS(I),6) .AND.
     :          .NOT. GTNDFL(BIGLST,NODLST,NODS(I),2) .AND.
     :          .NOT. GTNDFL(BIGLST,NODLST,NODS(I),8) ) THEN
c            IF( .NOT. GTNDFL(BIGLST,NODLST,NODS(I),2) ) THEN
C
               TNDCHK = TNDCHK + 1
C
               IF( IERR .NE. 0 ) THEN
                  ewrite(-1,8)
     :                   '*** ADAPT1: GOT ERROR FROM GTNDFL (NODTST)'
                  ewrite(-1,6) '   element: ',iptr
               END IF
C
               if( goterr ) then
                  ewrite(3,2) 
     :                 'Checking move of node ',nods(i),
     :                   (nodlst(ierr,nods(i)),ierr=1,3)
                  ierr = 0
                  debug  = .true.
                  goterr = .false.
               end if
c
               CALL NODTST( BIGLST, NODLST,
     :                      ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                      SZENLS, SZNELS, NNOD, NELM,
     :                      NODX, NODY, NODZ, ORGMTX,
     :                      NODS(I), IPTR,
     :                      CURDIF, CURORG, X,  Y,  Z,  INELM,
     :                      SUMDIF, SUMORG, XA, YA, ZA, INELMA,
     :                      ALWAYS, GF*DOTOP, MINCHG )
c
               if( debug .and. .not.goterr ) then
                  debug = .false.
                  goterr = .true.
               end if
c
               if( goterr ) then
                 ewrite(3,9) 'Done: ',
     :                            curdif,sumdif,x,y,z,xa,ya,za
               end if
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :           '*** CHANGED AFTER NODTST: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
               IF( IERR .NE. 0 ) THEN
                 ewrite(-1,8) 
     :                 '*** ADAPT1: GOT ERROR FROM NODTST'
                 ewrite(-1,6) 
     :                '   element,node: ',iptr,nods(i)
                 ewrite(-1,6) 
     :                '   conelm: ',conelm
                 ewrite(-1,6) 
     :                '   edglst: ',edglst
                 ewrite(-1,6) 
     :                '    nodes: ',nods
                 RETURN
               END IF
C
               IF( CURORG .GT. GF*DOTOP ) THEN
                 IF( CURDIF .LT. -ABS(MINCHG) ) THEN
                   IF( CURDIF .LT. BSTCUR ) THEN
                     BSTPTR = NODS(I)
                     BSTCEL = IPTR
                     BSTWHT = 10
                     BSTCUR = CURDIF
                     BSTORG = CURORG
                     BSTX   = X
                     BSTY   = Y
                     BSTZ   = Z
                     BSTINE = INELM
                     DUNONE = .TRUE.
                     if( curdif .lt. always ) goto 590
                     if( curorg+curdif .lt. gf*dotop ) goto 590
                   END IF
                ELSE IF( BSTCUR .EQ. 0.0 .AND.
     :                   SUMDIF .LT. -ABS(MINCHG) ) THEN
                   IF( SUMDIF .LT. BSTSUM ) THEN
                     BSTPTR = NODS(I)
                     BSTCEL = IPTR
                     BSTWHT = 10
                     BSTSUM = SUMDIF
                     BSTORG = SUMORG
                     BSTX   = XA
                     BSTY   = YA
                     BSTZ   = ZA
                     BSTINE = INELMA
                     DUNONE = .TRUE.
                     if( sumdif .lt. ALWAYS ) goto 590
                   END IF
c                ELSE IF( CURDIF .LT. BSTDIF .AND. .NOT. DUNONE ) THEN
c                  BSTPTR = NODS(I)
c                  BSTCEL = IPTR
c                  BSTWHT = 10
c                  BSTDIF = CURDIF
                 END IF
               END IF
C
               CALL STNDFL( BIGLST, NODLST, NODS(I), 6, .TRUE. )
C
               IF( IERR .NE. 0 ) THEN
                  ewrite(-1,8)
     :                   '*** ADAPT1: GOT ERROR FROM STNDFL (NODTST)'
                  ewrite(-1,6) '   element: ',iptr
               END IF
C

            else if( goterr ) then
               ewrite(3,2) 
     :              'Not checking move of node',nods(i),
     :               (nodlst(ierr,nods(i)),ierr=1,3)
               ierr = 0
            END IF
C
          END DO
C
         else if( goterr ) then
            ewrite(3,8) 
     :           'Not moving nodes of this elem...'
         END IF
C
  590    CONTINUE
C
         ELSLFT = ELSLFT - 1
C
         NELCHK = NELCHK + 1
         TELCHK = TELCHK + 1
C
         IF( ABS(ELSLFT-LSTLFT) .EQ. 4000 .OR.
     :    (NELCHK .GE. (NELEMS/500)*1000 .and. nelems .gt. 10000) ) THEN
c - this check was really only here to debug problematic cases with a few
C - elements swapping back and forth right near the end of an adapt sweep
c     :      (elslft .le. 100 .and. nelems .gt. 1000) ) THEN
c            TOTADD = TOTADD + NELADD
c            TOTSUB = TOTSUB + NELSUB
c            TOTCHK = TOTCHK + NELCHK
            if( lstlft .ne. elslft ) gothere = 0
c
            ewrite(1,3) 'LEFT/TOTAL:',ELSLFT,NELEMS,
     :              '      ADD/SUB/CHK:',NELADD,NELSUB,NELCHK,
     :              '      NODES:',NUMNDS
c            if( elslft .eq. 3100 .and. nelems .eq. 12553 .and.
c     :          neladd .eq. 2858 .and. nelsub .eq. 2814 .and.
c     :          nelchk .eq. 3352 .and. numnds .eq. 2922 ) then
c               goterr = .true.
c            end if
            if( elslft.le.100 .and. mvnode .and. nelems.gt.1000 ) then
c              if( goterr ) then
                gothere = gothere + 1
                if( gothere .gt. 1000 ) then
                   ewrite(-1,8) 
     :                  '*** STUCK IN LOOP - GIVING UP...'
                   ewrite(-1,*) 
     :                  nelems,elslft,lstlft,gothere
                  IERR = -1999
                  RETURN
                end if
c              else
c                write(myout,8)
c     :                 '>>>>>>>>>> SWITCHING ON DEBUGGING!!! <<<<<<<<<<'
c                goterr  = .true.
c                gothere = 0
c              end if
            else
              gothere = 0
            end if
C
  3         format( a, 2i9, a, 3i8, a, i9 )
C
            LSTLFT = ELSLFT
            NELADD = 0
            NELSUB = 0
            NELCHK = 0
            mvcnt  = 0
C
         ELSE IF( (1.0*NELCHK-1200.0)/1000 .GT. 1.0*NELEMS ) THEN
C
            ewrite(0,8) 
     :           '+++ WARNING: STUCK IN A LOOP???'
            LSTLFT = ELSLFT
            ewrite(0,3) 
     :           'LEFT/TOTAL:',ELSLFT,NELEMS,
     :           '      ADD/SUB/CHK:',NELADD,NELSUB,NELCHK,
     :           '      NODES:',NUMNDS
            NELADD = 0
            NELSUB = 0
            NELCHK = 0
            goterr = .TRUE.
            stop
C
         ELSE IF( (1.0*NELCHK-600.0)/1000 .GT. NELEMS
     :                            .and. .not. goterr ) THEN
            ewrite(0,8) 
     :           '+++ Warning: loop is looking likely...'
            ewrite(0,3) 
     :           'LEFT/TOTAL:',ELSLFT,NELEMS,
     :           '      ADD/SUB/CHK:',NELADD,NELSUB,NELCHK,
     :           '      NODES:',NUMNDS
            goterr = .TRUE.
         ELSE IF( (1.0*NELCHK)/1000 .GT. 1.0*NELEMS
     :                            .and. .not. shwchg ) THEN
            ewrite(0,8) 
     :           '+++ Warning: may be stuck in a loop...'
            ewrite(0,3) 
     :           'LEFT/TOTAL:',ELSLFT,NELEMS,
     :           '      ADD/SUB/CHK:',NELADD,NELSUB,NELCHK,
     :           '      NODES:',NUMNDS
            shwchg = .TRUE.
         else if( mvcnt .gt. 1100 ) then
            ewrite(-1,8) 
     :           '*** ADAPT1: ABORTING - TOO MANY MOVES'
            IERR = -11111
            RETURN
         else if( mvcnt .gt. 1000 .and. .not. shwchg ) then
            ewrite(0,8)
     :             '---===*** Warning: STUCK IN MOVE LOOP??? ***===---'
            shwchg = .TRUE.
         END IF
C
  600    CONTINUE
C
C - flag that we have checked this element
C
         CALL SETFLG( BIGLST, NODLST, IPTR, 3, .TRUE. )
         CALL SETFLG( BIGLST, NODLST, IPTR, 7, .TRUE. )
C
C - in the case when edge-splitting is inhibited, because one of the lists
C - is getting full, we want to switch off the flag that will be used to
C - prevent the next sweep from checking this element, otherwise we are in
C - danger of ignoring elements that still need work
C
         IF( SPLTED .NEQV. SPLED2 .AND. FLGSET .GT. 9 )
     :        CALL SETFLG(BIGLST,NODLST,IPTR,FLGSET,.FALSE.)
C
         IF( IERR .NE. 0 ) THEN
            ewrite(-1,8)
     :             '*** ADAPT1: GOT ERROR FROM SETFLG (Done elm)'
            ewrite(-1,6) '   element: ',iptr
         END IF
C
c      ELSE IF( goterr ) THEN
C
c         print*,'Skipping over edge ',iptr,nxt
C
      END IF
C
      IF( BSTCUR .EQ. 0.0 .AND. BSTSUM .EQ. 0.0 ) THEN
C
         if( goterr ) then
           ewrite(3,6)
     :       'Nothing to do, skipping to next ',iptr,nxt
         end if
         IPTR = NXT
C
C - have we reached an empty block?
C
         IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 20
c
         if( ierr.ne.0 ) then
            ewrite(-1,8)
     :             '*** ADAPT1: GOT ERROR FROM GETFLG (reached empty?)'
            ewrite(-1,6) '   element: ',iptr
            ewrite(-1,6) '   conelm: ',conelm
            ewrite(-1,6) '   edglst: ',edglst
            return
         end if
C
c         IF( DUNONE ) THEN
            BSTWHT = 0
            BSTDIF = 0.0
c         END IF
C
      END IF
C
C - yes... now do the actual splitting/shrinking etc. that was best
C
      IF( BSTWHT .EQ. 10 ) THEN
C
         if( goterr .or. shwchg ) then
            ewrite(3,1)
     :           '--- Info: Moving node ---',bstptr,bstcel,
     :           bstcur,bstorg
            ewrite(3,9) 
     :           ' Orig:',(nodlst(ierr,bstptr),ierr=1,3)
            ierr = 0
         end if
C
         if( goterr .and. .not. debug ) then
            debug = .true.
            goterr = .false.
         end if
c
         CALL NODMOV( BIGLST, NODLST,
     :                ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                SZENLS, SZNELS, NNOD,   NELM,
     :                NODX,   NODY,   NODZ,   ORGMTX, MINCHG.LT.0.0,
     :                BSTPTR, BSTCEL, BSTX, BSTY, BSTZ, BSTINE )
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :           '*** CHANGED AFTER NODMOV: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
c         print*,'Done...'
         IF( IERR .NE. 0 ) THEN
            ewrite(-1,8) 
     :           '*** ADAPT1: GOT ERROR FROM NODMOV'
            ewrite(-1,6) 
     :           '   Node,Element: ',bstptr,bstcel
            RETURN
         END IF
c
         if( debug .and. .not. goterr ) then
            debug = .false.
            goterr = .true.
            ewrite(3,9) 
     :           ' Done:',(nodlst(ierr,bstptr),ierr=1,3)
            ierr = 0
         end if
C
         CALL STNDFL( BIGLST, NODLST, BSTPTR, 6, .TRUE. )
C
      ELSE IF( BSTWHT .EQ. -2 ) THEN
C
         if( goterr .or. shwchg ) then
c         print*,' '
            ewrite(3,2) 
     :           '--- Info: Splitting element ---',bstptr,
     :                  bstcur,bstorg
c         print*,'   Element: ',bstptr
c         print*,'   Reduction: ',bstdif,bstcur,bstorg
c         print*,' '
         end if
C
         CALL ELMDIV( BIGLST, NODLST,
     :                ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                SZENLS, SZNELS, NNOD,   NELM,
     :                NODX,   NODY,   NODZ,   ORGMTX,
     :                BSTPTR )
C
         IF( IERR .NE. 0 ) THEN
            ewrite(-1,8) 
     :           '*** ADAPT1: GOT ERROR FROM ELMDIV'
            ewrite(-1,6) 
     :           '   element: ',bstptr
            RETURN
         END IF
C
      ELSE IF( BSTWHT .EQ. -1 ) THEN
C
         if( goterr .or. shwchg ) then
c         print*,' '
            ewrite(3,1) 
     :           '--- Info: Swapping face for edge ---',
     :           bstptr,bstcel,bstcur,bstorg
c         print*,'   Elements: ',bstptr,bstcel
c         print*,'   Reduction: ',bstdif,bstcur,bstorg
c         print*,' '
         end if
C
         CALL FC2EDG( BIGLST, NODLST, BSTPTR, BSTCEL )
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :           '*** CHANGED AFTER FC2EDG: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
         IF( IERR .NE. 0 ) THEN
            ewrite(-1,8) 
     :           '*** ADAPT1: GOT ERROR FROM FC2EDG'
            ewrite(-1,6)
     :           '   elements: ',bstptr,bstcel
            RETURN
         END IF
C
      ELSE IF( BSTWHT .LE. 4 .AND. BSTWHT .NE. 0 ) THEN
C
         CALL EDGELS( BIGLST, NODLST, BSTPTR, BSTCEL, ELPTRS,
     :                EDPTRS, NDPTRS, NUMEDE, SURFAC, INTRNL )
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :          '*** CHANGED AFTER EDGELS2: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
          IF( IERR .NE. 0 ) THEN
             ewrite(-1,8) 
     :            '*** ADAPT1: GOT ERROR FROM EDGELS'
             ewrite(-1,6)
     :            '   edge,element: ',bstptr,bstcel
             RETURN
          END IF
C
         IF( BSTWHT .EQ. 1 ) THEN
C
         if( goterr .or. shwchg ) then
c         print*,' '
            ewrite(3,4)
     :           '--- Info: Splitting edge ---',bstptr,bstcel,
     :                 bstcur,bstorg,bstrat
c         print*,'   Edge: ',bstptr
c         print*,'   Reduction: ',bstdif,bstcur,bstorg
c         print*,' '
         end if
            CALL SPLEDG( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   BSTPTR, ELPTRS, EDPTRS, NDPTRS, NUMEDE,
     :                   BSTX, BSTY, BSTZ, BSTINE, SURFAC, INTRNL,
     :                   NEWED1, NEWED2, NEWEL1, NEWEL2 )
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :           '*** CHANGED AFTER SPLEDG: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
            IF( IERR .NE. 0 ) THEN
               ewrite(-1,8) 
     :              '*** ADAPT1: GOT ERROR FROM SPLEDG'
               ewrite(-1,6) 
     :              '   edge,element: ',bstptr,bstcel
               ewrite(-1,6) 
     :              '   numede,surfac: ',numede,surfac
               RETURN
            END IF
C
         ELSE IF( BSTWHT .EQ. 2 ) THEN
C
         if( goterr .or. shwchg ) then
c         print*,' '
            ewrite(3,1)
     :           '--- Info: Shrinking edge ---',
     :                bstptr,bstcel,bstcur,bstorg
c         print*,'   Edge: ',bstptr
c         print*,'   Reduction: ',bstdif,bstcur,bstorg
c         print*,' '
         end if
            CALL SHREDG( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   BSTPTR, ELPTRS, EDPTRS, NDPTRS, NUMEDE,
     :                   BSTRAT, BSTX, BSTY, BSTZ, BSTINE,
     :                   SURFAC, INTRNL )
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :           '*** CHANGED AFTER SHREDG: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
            IF( IERR .NE. 0 ) THEN
               ewrite(-1,8) 
     :              '*** ADAPT1: GOT ERROR FROM SHREDG'
               ewrite(-1,6) 
     :              '   edge,element: ',bstptr,bstcel
               ewrite(-1,6) 
     :              '   numede,surfac: ',numede,surfac
               RETURN
            END IF
C
         ELSE IF( BSTWHT .EQ. 3 .OR. ABS(BSTWHT) .EQ. 4 ) THEN
C
c         if( goterr .or. bstwht .ne. 3 .or. shwchg ) then
         if( goterr .or. shwchg ) then
c         print*,' '
            if( bstwht .eq. 3 ) then
               ewrite(3,4) 
     :              '--- Info: Swapping edge for face ---',
     :              bstptr,bstcel,bstcur,bstorg
            else
               ewrite(3,5) 
     :              '--- Info: Swapping edge for edge ---',
     :              bstwht,bstptr,bstcel,bstcur,bstorg
            end if
c         print*,'   Edge: ',bstptr
c         print*,'   Reduction: ',bstdif,bstcur,bstorg
c         print*,' '
         end if
            CALL EDG2FC( BIGLST, NODLST, BSTPTR, ELPTRS, EDPTRS, NDPTRS,
     :                   BSTREG, BSTWHT )
c
         if( trace .gt. 0 ) then
          if( biglst(trcnd,trace) .ne. cnd2vl ) then
            ewrite(0,6)
     :          '*** CHANGED AFTER EDG2FC: ',cnd2vl,biglst(trcnd,trace)
            cnd2vl = biglst(trcnd,trace)
cccc            pause
          end if
         end if
C
            IF( IERR .NE. 0 ) THEN
               ewrite(-1,8) 
     :              '*** ADAPT1: GOT ERROR FROM EDG2FC'
               ewrite(-1,6) 
     :              '   edge,element: ',bstptr,bstcel
               ewrite(-1,6) 
     :              '   numede,surfac: ',numede,surfac
               RETURN
            END IF
C
         END IF
C
      ELSE if( goterr ) then
C
         ewrite(3,8) 
     :        '--- Info: nothing done on this pass ---'
c         BSTDIF = 0.0
C
      END IF
c
  1   format( a, 2i8, 10(1pe13.5) )
  2   format( a,  i8, 10(1pe13.5) )
  4   format( a, 2i8, 10(1pe13.5) )
  5   format( a, i3, 2i8, 10(1pe13.5) )
  6   format( a, 10i8 )
  7   format( a, 10(1pe13.5) )
  8   format( a )
  9   format( a, 10(1pe16.8) )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( BSTCUR .NE. 0.0 .OR. BSTSUM .NE. 0.0 ) THEN
C
         if( goterr ) then
           ewrite(1,9)
     :          '--- Info: reduced by ',BSTCUR,BSTSUM
         end if
c
         IF( BSTCUR .NE. 0.0 ) THEN
            CHANGE = CHANGE + BSTCUR
         ELSE
            CHANGE = CHANGE + BSTSUM
         END IF
C
         IPTR = NEXTLS( MYNXT )
         if( iptr .ne. nxt .and. goterr ) then
           ewrite(3,6) 
     :           '   NXTELM CHANGED: ',NXT,IPTR
           ewrite(3,6)
     :            '   STTBIG,ENDBIG,EMTBIG:',STTBIG,ENDBIG,EMTBIG
         end if
C
C - have we reached an empty block?
C
c         if( totchg .gt. 100 ) goto 900
         IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 20
c
         if( ierr.ne.0 ) then
            ewrite(-1,8) 
     :           '*** ADAPT1: GOT ERROR FROM GETFLG 3'
            ewrite(-1,6) 
     :           '   element: ',iptr
            ewrite(-1,6) 
     :           '   conelm: ',conelm
            ewrite(-1,6) 
     :           '   edglst: ',edglst
            return
         end if
C
         BSTDIF = 0.0
C
      END IF
C
c      IF( BSTDIF .LT. 0.0 ) THEN
C
c         PRINT*,'--- Info: reduced by ',BSTDIF
c         CHANGE = CHANGE + BSTDIF
C
c         IPTR = NEXTLS( MYNXT )
c         if( iptr .ne. nxt ) then
c           print*,'   NXTELM CHANGED: ',NXT,IPTR
c           print*,'   STTBIG,ENDBIG,EMTBIG:',STTBIG,ENDBIG,EMTBIG
c         end if
C
C - have we reached an empty block?
C
c         IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 17
c
c         if( ierr.ne.0 ) then
c            print*,'*** ADAPT1: GOT ERROR FROM GETFLG 4'
c            print*,'   element: ',iptr
c            print*,'   conelm: ',conelm
c            print*,'   edglst: ',edglst
c            return
c         end if
C
c      END IF
C
 900  CALL RMVNXT( MYNXT )
C
c      TOTADD = TOTADD + NELADD
c      TOTSUB = TOTSUB + NELSUB
c      TOTCHK = TOTCHK + NELCHK
C
      TOPBIG = STTBIG
      IPTR = STTBIG
C
 950  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
         TOPBIG = MAX( TOPBIG, IPTR )
         IPTR = BIGLST( NXTBIG, IPTR )
         GOTO 950
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
         TOPBIG = MAX( TOPBIG, IPTR )
         IPTR = BIGLST( NXTBIG, IPTR )
         TOPBIG = MAX( TOPBIG, IPTR )
         IPTR = BIGLST( NXTBIG, IPTR )
         TOPBIG = MAX( TOPBIG, IPTR )
         IPTR = BIGLST( NXTBIG, IPTR )
         GOTO 950
      END IF
C
      TOPNOD = STTNOD
      IPTR = STTNOD
C
 960  IF( GTNDFL(BIGLST,NODLST,IPTR,1) ) THEN
         TOPNOD = MAX( TOPNOD, IPTR )
         IPTR = INT(NODLST(NXTNOD,IPTR))
         GOTO 960
      END IF
C
      ewrite(1,9) 
     :     '--- Info: sum of reductions:',CHANGE
      ewrite(1,6) 
     :     '--- Info: total elems checked:',TELCHK
      ewrite(1,6) 
     :     '--- Info: total edges checked:',TEDCHK
      ewrite(1,6) 
     :     '--- Info: total nodes checked:',TNDCHK
C
      ewrite(1,8) ' '
      ewrite(1,'(A,i8,A,i8,A,i8)') 
     :     '    Nodes:',NUMNDS,
     :     '       Elements:',NELEMS,
     :     '       Edges:',NEDGES
c      goterr = .true.
      if( goterr ) then
         ewrite(3,6) 
     :        'TOPBIG,TOPNOD: ',TOPBIG,TOPNOD
         ewrite(3,6) 
     :        'MAXBIG,MXNODS: ',MAXBIG,MXNODS
      end if
C
      CALL RATOPC( TOPBIG, MAXBIG, STR1 )
      CALL RATOPC( 3*NELEMS+NEDGES, TOPBIG, STR2 )
      ewrite(1,'(A,A,A,A)') 
     :     '   BIGLST usage: ',STR1,
     :     '   BIGLST efficiency: ',STR2
      CALL RATOPC( TOPNOD, MXNODS, STR1 )
      CALL RATOPC( NUMNDS, TOPNOD, STR2 )
      ewrite(1,'(A,A,A,A)') 
     :     '   NODLST usage: ',STR1,
     :     '   NODLST efficiency: ',STR2
C
      IF( SPLTED .AND. .NOT. SPLED2 ) THEN
         ewrite(0,8) 
     :        '+++ ADAPT1: WARNING! LISTS WERE FILLING UP!'
         ewrite(0,8)
     :  '   ** Edge splitting was switched off during adapt sweep **'
      END IF
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMTOT = TIMEND - TIMSTT
C
      CALL SHWTIM
C
      ADPTNG = .FALSE.
      
      RETURN
      END
C
C
C
      SUBROUTINE RATOPC( TOP, BOT, STR )
C
      IMPLICIT NONE
C
      INTEGER TOP, BOT
      CHARACTER*7 STR
C
      INTEGER INTGR, DCML
C
      INTGR = INT(10000.0*TOP/BOT+0.5)
      DCML  = MOD(INTGR,100)
      INTGR = INTGR/100
C
      WRITE(STR,'(I3,A,I2,A)') INTGR,'.',DCML,'%'
      IF( DCML .LT. 10 ) STR(5:5) = '0'
C
      RETURN
      END
c
c
c
      real function getdtm( a )
c
      real a(3,3)
c
      getdtm = a(1,1)*( a(2,2)*a(3,3) - a(2,3)*a(3,2) )
     :       + a(2,1)*( a(3,2)*a(1,3) - a(1,2)*a(3,3) )
     :       + a(3,1)*( a(1,2)*a(2,3) - a(2,2)*a(1,3) )
c
      return
      end
c
