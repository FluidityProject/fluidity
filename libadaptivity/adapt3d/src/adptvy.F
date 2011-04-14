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
      SUBROUTINE ADPTVY( INTARR, INTSIZ, RLARR,  RLSIZ,
     :                   GEOM3D, SRFGMY, USEQ,
     :                   NNOD,   NELM,   NSELM,  ABSOLUTEMXNODS,
     :                   SZENLS, ENLBAS, ENLIST, ELMREG,
     :                   CLCGMY, SZSNLS, SNLBAS, SNLIST, SURFID,
     :                   PRDNDS, NPRDND,
     :                   NODX,   NODY,   NODZ,
     :                   INTNND, INTNEL, INTSZL, INTENL, INTENB,
     :                   INTNDX, INTNDY, INTNDZ,
     :                   ORGMTX, OLDFLD, NFREE,  TOTFRE, NFIELD,
     :                   XPCTEL, NWNNOD, NWNELM, NWNSEL,
     :                   NWSZEN, NWSZSN, NWSZNN, NWNDLC, NWSROW,
     :                   NWENLB, NWENLS, NWSNLB, NWSNLS, NWSFID,
     :                   NWELRG, NWNODX, NWNODY, NWNODZ,
     :                   NEWMTX, NEWFLD,
     :                   BIGLST, NODLST,
     :                   DOTOP,  MINCHG, NSWEEP, MSHOPT, TWOSTG, TOGTHR,
     :                   GATHER, SCATER, NGATH,  NHALO,  PNOD,
     :                   ATOSEN, ATOREC, NPROC,debug_level,dbg,chcnsy)
C-----------------------------------------------------------------------
C
C - This subroutine takes the current mesh information (in ENLIST & ELMREG)
C - and uses the metrics (in ORGMTX) to adapt the mesh. The new mesh info
C - is returned in the free space (INTARR and RLARR), pointed to by various
C - pointers (NWENLB, NWENLS etc. - see below). Size info for the new mesh
C - is also returned (NWNNOD, NWNELM etc. - see below), including the size
C - for certain lists which may be created by a user's code (eg. node
C - connectivity lists, node to element list etc). The fields given in
C - OLDFLD are also interpolated onto the new mesh, and returned in the
C - real array (RLARR), pointed to by NEWFLD.
C - Note that the interpolated metrics are also returned in RLARR (pointed
C - to by NEWMTX) if required. This can be switched off by supplying
C - NEWMTX as a negative integer variable.
C
C Things still to be done:
C
C - Sort out the different formats used for mesh information (ENLIST/ENLBAS)
C   by Chris and Cass... we could have a switch, or OFFSET, sent down
C   which will flick between the two.
C
C-----------------------------------------------------------------------
      use write_log
      use AdaptProgress
      IMPLICIT NONE
C
      INCLUDE 'blkprm.i'
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkerr.i'
C
      REAL ANGLEFACT
      COMMON / BLKTMP /  ANGLEFACT
C
      INTEGER RLSIZ, INTSIZ, NNOD, NELM, NSELM, ABSOLUTEMXNODS, SZENLS,
     :        NWNNOD, NWNELM, NWNSEL, NWSZEN, NWSZSN, NWSZNN, SZSNLS,
     :        NWENLB, NWENLS, NWSNLB, NWSNLS, NWSFID, NWELRG, NWSROW,
     :        NWNDLC, NWNODX, NWNODY, NWNODZ, NEWFLD, NEWMTX,
     :        NGATH, NHALO, NPROC, PNOD, NSWEEP, NPRDND, debug_level
C
      INTEGER INTARR(INTSIZ), ENLBAS(NELM+1), ENLIST(SZENLS),
     :        SNLIST(SZSNLS), SNLBAS(NSELM+1), SURFID(NSELM),
     :        ELMREG(NELM), XPCTEL, PRDNDS(NPRDND),
     :        GATHER(NGATH), SCATER(NHALO), ATOSEN(NPROC), ATOREC(NPROC)
C
      REAL RLARR(RLSIZ), NODX(NNOD), NODY(NNOD), NODZ(NNOD),
     :     MINCHG, DOTOP
      integer, dimension(:), allocatable :: tmparr
C
      LOGICAL GEOM3D, SRFGMY, USEQ, TOGTHR, CLCGMY, MSHOPT(6), TWOSTG,
     :        dbg, chcnsy
C
C - these arrays and sizes are all based on the interpolation mesh
C - (which can be, and often is, identical to the mesh-to-adapt -see notes)
C
      INTEGER INTNND, INTNEL, INTSZL, TOTFRE, NFIELD
C
      INTEGER INTENB(INTNEL+1), INTENL(INTSZL), NFREE(NFIELD)
C
      REAL INTNDX(INTNND), INTNDY(INTNND), INTNDZ(INTNND),
     :     ORGMTX(9*INTNND), OLDFLD(INTNND*TOTFRE)
C
C-----------------------------------------------------------------------
C
C The following should be sent into the subroutine:
C
C - INTARR is the array of free integers
C - RLARR is the array of free reals
C - INTSIZ is the size of INTARR
C - RLSIZ is the size of RLARR
C - GEOM3D is a flag which is true if this is a 3D problem
C - SRFGMY is TRUE if the surface mesh should be kept intact during adapting
C - CLCGMY is TRUE if the geometry should be calculated, and ignore SNLIST
C
C The following give information about the initial mesh to be adapted:
C
C - NNOD is the number of (local+halo) nodes in the original mesh
C - NELM is the number of (local+halo) volume elements
C - NSELM is the number of (local+halo) surface elements
C - ABSOLUTEMXNODS is the absolute maximum allowed no. of nodes to return
C - SZENLS is the size of ENLIST ( = 4*NELM in 3D)
C - ENLBAS is the base pointer into the (local+halo) volume element-node list
C - ENLIST the (local+halo) volume element node list
C - ELMREG contains (local+halo) elementwise region information
C - SNLIST is the surface element node list
C - SZSNLS is the size of SNLIST
C - SNLBAS is the base pointer into SNLIST
C - SURFID is the ID list for surfaces (if to be used for geometry)
C - NODX is the (local+halo) X coordinates
C - NODY is the (local+halo) Y coordinates
C - NODZ is the (local+halo) Z coordinates
C
C The following are all for the interpolation mesh (i.e. the mesh which
C is used to interpolate the metric and fields during and after the
C adaptivity respectively). The interpolation mesh is very often the
C same as the initial mesh to be adapted, in which case simply send
C INTNND=NNOD, INTNEL=NELM, INTSZL=SZENLS, INTENL=ENLIST, INTENB=ENLBAS,
C INTNDX=NODX, INTNDY=NODY and INTNDZ=NODZ
C
C - INTNND is the no. of nodes in the interpolation mesh
C - INTNEL is the no. of elements in the interpolation mesh
C - INTSZL is the size of the element-node list for the interpolation mesh
C - INTENL is the element-node list for the interpolation mesh
C - INTENB is the base pointer for the interpolation mesh 
C - INTNDX is the X coordinates for the nodes of the interpolation mesh
C - INTNDY is the Y coordinates for the nodes of the interpolation mesh
C - INTNDZ is the Z coordinates for the nodes of the interpolation mesh
C - ORGMTX is the metrics returned by SPAERR (for the interpolation mesh)
C - OLDFLD is the field data (for the interpolation mesh)
C - NFREE contains a list of the number of freedoms for each field
C - TOTFRE is the sum of the freedoms over all the fields
C - NFIELD is the number of fields (so OLDFLD has length INTNND*TOTFRE)
C 
C Another couple of paramaters to be sent into the subroutine:
C
C - XPCTEL is the expected number of elements (guestimated by SPAERR)
C - current_debug_level is the debug level used for EWRITE statements
C
C The following contain the information needed for a parallel run:
C
C - GATHER is the gather array for local nodes info to be sent
C - SCATER is the scatter array for local node info to be received
C - NGATH is the size of the gather array, GATHER
C - NHALO is the number of '1st order' halo nodes (therefore size of SCATER)
C - PNOD is the number of strictly local nodes (but NNOD-PNOD <> NHALO)
C - ATOSEN is the base pointer into GATHER for each processor
C - ATOREC is the base pointer into SCATER for each processor
C
C The following sizes are returned by the subroutine:
C
C - NWNNOD is the number of nodes in the adapted mesh
C - NWNELM is the number of volume elements in the adapted mesh
C - NWNSEL is the number of surface elements in the adapted mesh
C - NWSZEN is the size of the new volume element node list
C - NWSZSN is the size of the new surface element node list
C - NWSZNN is the size of the new node connectivity list
C
C The following are returned as pointers into INTARR:
C
C - NWENLB is the base pointer into the new volume element node list pointer
C - NWENLS is the new volume element node list pointer
C - NWSNLB is the base pointer into the new surface element node list pointer
C - NWSNLS is the new surface element node list pointer
C - NWELRG is the new elementwise region information pointer
C
C The following are returned as pointers into RLARR:
C
C - NWNODX is the new X coordinates pointer
C - NWNODY is the new Y coordinates pointer
C - NWNODZ is the new Z coordinates pointer
C - NEWMTX is the interpolated metrics pointer
C - NEWFLD is the interpolated fields pointer
C
C The following allow some control over the level of adaptivity:
C
C - DOTOP is the minimum element functional to allow to be changed (min 0.1)
C - MINCHG is the minimum relative change allowed to a local functional (min 0.01)
C - NSWEEP is the number of adaptive sweeps to perform ( min 1 )
C - MSHOPT is an array of six logicals which switch on/off mesh operations
C          the order of entries is the same as those listed in adapt1.f:
C           1 : split edges if true
C           2 : collapse edges if true
C           3 : perform edge to face and edge to edge swapping if true
C           4 : perform face to edge swapping if true
C           5 : split elements (do not use this yet)
C           6 : move nodes if true
C - TWOSTG if true, two stages of adapting, with no refinement on first
C - TOGTHR if true, lumps node movement adaptivity in with connectivity changes
C
C
C ======= Extra comments - please take note! =======
C
C INTARR must be large enough to contain the following...
C
C - node-element list base pointer for old mesh: NNOD+1
C - node-element list for old mesh: SZENLS
C - element connectivity list for the old mesh: 4*NELM (in 3D)
C - the adapting edge/element info: >30*NWNELM, probably >35*NWNELM
C - the new elementwise region info: NWNELM
C - the new volume element-node list base pointer: NWNELM+1
C - the new volume element-node list: 4*NWNELM (in 3D)
C - the new surface element-node list base pointer: NWNSEL+1
C - the new surface element-node list: 3*NWNSEL (in 3D)
C
C RLARR must be large enough to contain the following...
C
C - the adapting node info: >17*NWNNOD, probably >20*NWNNOD
C - the new node coordinates: 3*NWNNOD
C - the interpolated field: TOTFRE*NWNNOD
C - the interpolated metrics (if required): 9*NWNNOD
C
C-----------------------------------------------------------------------
      INTEGER NELBAS, NELIST, SZNELS, EELIST, SIZOWN, SIZLOC, NSWPS,
     :        XPCTND, I, SIZMTX, BIGLST, NODLST, FLCHK1, FLCHK2, FLCHK3,
     :        FLCHK4, FLCHK5, FLCHK6, mxsfid
C
      REAL ALWAYS, CHANG1, CHANG2
      real xd, yd, zd, xdf, ydf, zdf
c
      INTEGER*8 INTSUM
      REAL*8 RLSUM
      LOGICAL NOCNCT, DUNSTG, EDSPLT, ELSPLT

      call set_debug_level(debug_level)

C     Initialise AdaptProgress module
      call initialise(NNOD, 0.5)

      ANGLEFACT = 0.0
C
      IF( MSHOPT(5) ) THEN
         MSHOPT(5) = .FALSE.
C         ewrite(3, '(a)')
C     :         '+++ ADPTVY: Element splitting unsupported - turned off'
      END IF
C
      DUNSTG = .FALSE.
      IF( .NOT. TWOSTG ) THEN
         DUNSTG = .TRUE.
      ELSE IF( .NOT. MSHOPT(1) .AND. .NOT. MSHOPT(5) ) THEN
C - must make sure there is something switched on for two-stage adapt
         MSHOPT(1) = .TRUE.
         ewrite(0, '(a)')
     :     '+++ ADPTVY: Turned on edge splitting for two stage adapt'
      END IF
C
C      PNOD = NNOD - NHALO
C
      debug  = dbg
      cnsych = chcnsy
      if( chcnsy ) then
         ewrite(1,*)
     :   '--- ADPTVY: Consistency checks don''t work - switching off...'
      end if
      cnsych = .false.
C
      do i = 1, nhalo
        if(scater(i).le.0) ierr=54545
      end do
      if(ierr.eq.54545) then
         ewrite(-1, '(a)')
     :        '*** ADPTVY: GOT BAD NODE IN SCATER'
      end if
c
      do i = 1, ngath
        if(gather(i).le.0) ierr=54546
      end do
      if(ierr.eq.54546) then
         ewrite(-1, '(a)') 
     :        '*** ADPTVY: GOT BAD NODE IN GATHER'
      end if
c
      if(nselm.gt.0) then
         mxsfid = surfid(1)
         do i = 2, nselm
           mxsfid = max(surfid(i),mxsfid)
         end do
      else
         mxsfid = 0
      end if
c
      if( debug .or. ierr.ne.0 ) then
         ewrite(3,'(a)')
     :        'In ADPTVY with following parameters:'
         ewrite(3,'(a,6L3)')
     :          '  GEOM3D,CLCGMY,SRFGMY,TOGTHR,TWOSTG,CHCNSY:',
     :             GEOM3D,CLCGMY,SRFGMY,TOGTHR,TWOSTG,CHCNSY
         ewrite(3,'(a,4i7)')
     :        '  NNOD,NELM,NSELM,PNOD:',
     :        NNOD,NELM,NSELM,PNOD
         ewrite(3,'(a,3i7)')
     :        '  NPROC,NGATH,NHALO:',
     :        NPROC,NGATH,NHALO
         ewrite(3,'(a,2(1pe13.5),i4)')
     :        '  DOTOP,MINCHG,NSWEEP:',
     :        DOTOP,MINCHG,NSWEEP
         ewrite(3,'(a,i6)')
     :        '  MXSFID:',mxsfid
         if( ierr .ne. 0 ) then
            nwnnod = -1
            goto 999
         end if
      else
         ewrite(3,'(a,2i6)') 
     :        '--- ADPTVY: NSELM,MXSFID=',
     :        nselm,mxsfid
      end if
C
      USEQLY = USEQ
C
C - initialise the powers of two used for setting bits in
C - GETFLG/SETFLG/GTNDFL/STNDFL etc.
C
      CALL INIPWR
C
c      intsum = 0
c      do i = 1, nelm+1
c         intsum = intsum + enlbas(i)
c      end do
c      print*,'Sum of ENLBAS: ',intsum
c
c      intsum = 0
c      do i = 1, szenls
c         intsum = intsum + enlist(i)
c      end do
c      print*,'Sum of ENLIST: ',intsum,szenls
c
c      rlsum = 0.0
c      do i = 1, nnod
c         rlsum = rlsum + abs(nodx(i))
c      end do
c      print*,'Sum of NODX: ',rlsum
c
c      rlsum = 0.0
c      do i = 1, nnod
c         rlsum = rlsum + abs(nody(i))
c      end do
c      print*,'Sum of NODY: ',rlsum
c
c      rlsum = 0.0
c      do i = 1, nnod
c         rlsum = rlsum + abs(nodz(i))
c      end do
c      print*,'Sum of NODZ: ',rlsum
C
C The basic format to follow is as follows:
C
      IERR = 0
C
C - the expected number of nodes is usually about a fifth the expected
C - number of elements
C
      XPCTND = XPCTEL / 5
C
      IS3DMS = GEOM3D
C
      NELBAS = 1
      NELIST = NELBAS + NNOD + 1
      SZNELS = SZENLS
      EELIST = NELIST + SZNELS
C
C - This forms NELIST, the nodal element list (i.e. elems round a node),
C - and its base pointer (NELBAS), given the element node list (ENLIST),
C - its base pointer (ENLBAS), its size (SZENLS) and the number of nodes
C - (NNOD) and elements (NELM). The size of NELIST (SZNELS) is also
C - returned but always turns out to be the same size as the element node
C - list (SZENLS).
C
      CALL MKNELS( NNOD, NELM, ENLBAS, ENLIST, SZENLS,
     :             INTARR(NELBAS), INTARR(NELIST), SZNELS, 1 )
C
      if( debug ) then
         ewrite(3,'(a)')
     :        'Made NELBAS: '
         ewrite(3,'(8i8)')
     :        (INTARR(NELBAS+i),i=0,min(15,nnod))
         ewrite(3,'(a)') 
     :        'Made NELIST: '
         ewrite(3,'(8i8)') 
     :        (INTARR(NELIST+i),i=0,intarr(nelbas+min(15,nnod)+1))
      end if
C
C - This forms the element connectivity list (EELIST) given the element
C - node list (ENLIST), its base pointer (ENLBAS) and its size (SZENLS),
C - and the info returned by the previous call to MKNELS containing the
C - nodal element list etc., plus info about whether this is a 3D problem
C - (GEOM3D - logical).
C - It also uses the surface IDs to set the correct negative values for
C - the element connectivity list.
C
      CALL MKESCL( NNOD, NELM, ENLBAS, ENLIST, SZENLS,
     :             INTARR(NELBAS), INTARR(NELIST), SZNELS,
     :             NSELM, SZSNLS, SNLBAS, SNLIST, SURFID,
     :             GEOM3D, INTARR(EELIST) )
C
      if( debug ) then
         ewrite(3,'(a)') 
     :        'Made ECLIST:'
         ewrite(3,'(8i7)') 
     :        (INTARR(EELIST+i),i=0,4*min(NELM,15)-1)
      end if
C
      ewrite(2,'(A)') 
     :     'Checking consistency of elements...'
      CALL CHKCNSTNCY( NNOD, NELM, ENLBAS, ENLIST, SZENLS,
     :                 GEOM3D, INTARR(EELIST), NODX, NODY, NODZ )
      ewrite(2,'(A)') 'Passed!'
C
      IF( IS3DMS ) THEN
         BIGLST = EELIST + 4*NELM
      ELSE
         BIGLST = EELIST + 3*NELM
      END IF
      MAXBIG = ( INTSIZ - BIGLST )/SIZBIG
C
      IF( MAXBIG .LE. NELM*3 ) THEN
         WRITE(0,'(a,i10)')
     :        '*** ADPTVY: NOT ENOUGH SPACE FOR ELEMENTS:',MAXBIG
         NWNNOD = -1
         IERR = -5001
         GOTO 999
      END IF
C
      ewrite(1,'(A,I9,A)')
     :     '--- ADPTVY: Space set aside for about',
     :     2*MAXBIG/9,' elements'
C
      IF( 2*MAXBIG/9 .LT. 1.2*XPCTEL ) then
         ewrite(0,'(a)')
     :      '+++ ADPTVY: WARNING: MAY NOT BE ENOUGH SPACE FOR ELEMENTS'
      end if
C
      NODLST = 1
      MXNODS = ( RLSIZ - NODLST - XPCTND*(TOTFRE+3) )/SIZNOD
C
      IF( MXNODS .LT. NNOD ) THEN
         WRITE(0,'(a,i8)')
     :         '*** ADPTVY: NOT ENOUGH SPACE FOR NODES:', MXNODS
         NWNNOD = -1
         IERR = -5002
         GOTO 999
      END IF
C
C - want to make sure it doesn't go beyond the absolute maximum allowed nodes
C - but we do build in a teeny bit of leeway (adapt1 complains once the no.
C - of nodes reaches 90% of MXNODS, as defined here, so the 10/9 factor below
C - actually allows adaptivity to potentially go right up to the absolute
C - maximum before complaining)
C
      IF( ABSOLUTEMXNODS .GT. 0 )
     :         MXNODS = min(MXNODS, ABSOLUTEMXNODS*10/9)
C
      ewrite(1,'(A,I10,A)') 
     :     '--- ADPTVY: Space set aside for',
     :     MXNODS,' nodes'
C
      IF( MXNODS .LT. 1.2*XPCTND ) then
         ewrite(0,'(a)')
     :        '+++ ADPTVY: WARNING: MAY NOT BE ENOUGH SPACE FOR NODES'
      end if
C
C - This forms the adaptive mesh data structures given the fixed mesh
C - information as above and also the element material/region id (ELMREG).
C - The geometry of the adapted mesh will retain the geometry of these
C - element material/region boundaries. ORGMTX should contain the nine
C - entries for the metric of each node.
C
      CALL MKADPT( INTARR(BIGLST), RLARR(NODLST),
     :             CLCGMY, SZSNLS, SNLBAS, SNLIST, SURFID,
     :             PRDNDS, NPRDND,
     :             NODX, NODY, NODZ, ORGMTX, NNOD, ELMREG,
     :             ENLBAS, ENLIST, INTARR(NELBAS), INTARR(NELIST),
     :             INTARR(EELIST), NELM, NSELM,
     :             SZENLS, SZNELS, GEOM3D, SRFGMY,
     :             PNOD, NHALO, SCATER, NGATH, GATHER )
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,'(a)') '*** ADPTVY: GOT ERROR FROM MKADPT'
         NWNNOD = -1
         GOTO 999
      END IF
C
      IF( MINCHG .LT. 0.01 ) MINCHG = 0.01
      IF( DOTOP  .LT. 0.1  ) DOTOP  = 0.1
      ALWAYS = -1.0
      NSWPS  = MAX( NSWEEP, 0 )
C
      IF( NSWPS .EQ. 0 ) THEN
         ewrite(1,'(a)') 
     :        '--- ADPTVY: Number of sweeps is zero'
         DUNSTG = .TRUE.
         GOTO 600
      ENDIF
C
C -switch off refinement if we are performing a two stage adapt
C
      IF( TWOSTG ) THEN
         EDSPLT = .FALSE.
         ELSPLT = .FALSE.
         DUNSTG = .FALSE.
         ewrite(1,'(a)')
     :      '--- ADPTVY: STARTING STAGE ONE WITH REFINEMENT OFF'
      ELSE
         EDSPLT = MSHOPT(1)
c - MSHOPT(5) (element splitting) is currently ignored!!
c         ELSPLT = MSHOPT(5)
         ELSPLT = .FALSE.
         DUNSTG = .TRUE.
      END IF
C
      NOCNCT = .NOT. ( EDSPLT .OR. MSHOPT(2) .OR. MSHOPT(3) .OR.
     :                 ELSPLT .OR. MSHOPT(4) )
C
      IF( NOCNCT .AND. .NOT. MSHOPT(6) ) THEN
         ewrite(1,'(a)')
     :       '--- ADPTVY: All mesh operations are switched off'
         GOTO 600
      END IF
C
  10  CONTINUE
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
C         CALL ADAPT1( INTARR(BIGLST), RLARR(NODLST),
C     :                ENLBAS, ENLIST, SZENLS,
C     :                INTARR(NELBAS), INTARR(NELIST), SZNELS,
C     :                INTARR(EELIST),
C     :                NODX,   NODY,   NODZ,   ORGMTX, NNOD,   NELM,
C     :                SPLTED, CLPSED, EDSWAP, FC2ED, SPLTEL, MVNODE,
C     :                FLGSET, FLGCHK, ALWAYS, DOTOP, MINCHG, CHANGE )
C
C - The FLGSET and FLGCHK allow for selectivity in elements to be checked
C - for adaptivity. These integers can be from 14 to 29 inclusive. ADAPT1
C - will set flag FLGSET of an element before attempting to adapt. Hence,
C - when ADAPT1 returns, any elements that remain with flag FLGSET set
C - have not been changed. However, elements with flag FLGCHK set will
C - not be checked for adaptivity (although they may be changed due to a
C - neighbouring element being changed). So, if FLGSET and FLGCHK are
C - alternated between values from one call of ADAPT1 to another, it
C - allows redundant checking to be minimised.
C
C - ALWAYS is the amount of functional change to always perform as soon
C - as it is found (i.e. ignore other checking for that element). It
C - should be negative (eg. -1.0 or -2.0). DOTOP is the minimum element
C - functional to allow to be changed (eg. 0.15). MINCHG is the minimum
C - relative change allowed to a local functional (eg. 0.01). CHANGE is
C - returned as the sum of local functional changes - note: it does not
C - represent the total change in the global functional.
C
C - It is advised that ADAPT1 should be called first with just node
C - movement then with everything but node movement (SPLTEL is not yet
C - functional and is excluded from this advice). The flags (FLGSET and
C - FLGCHK) should be independent for each of the two sets of calls,
C - eg. interchange 14 & 15 for the node movement call, and interchange
C - 16 & 17 for the other.
C
C
      FLCHK1 = 0
      FLCHK2 = 0
      FLCHK3 = 14
      FLCHK4 = 15
      FLCHK5 = 16
      FLCHK6 = 17
C
C - Do the correct number of adaptive sweeps
C
      DO I = 1, NSWPS/2+1
C
c         IF( ANGLEFACT .GE. 0.1 ) THEN
c            ANGLEFACT = 0.1
c         ELSE IF(I .GE. NSWPS/2 ) THEN
c            ANGLEFACT = 0.1
c            FLCHK1 = 0
c         ELSE IF(I .GT. NSWPS/4 ) THEN
c            ANGLEFACT = 0.08*(I-NSWPS/4)/(NSWPS/2-NSWPS/4)+0.02
c            FLCHK1 = 0
c         ELSE
c            ANGLEFACT = 0.0
c         END IF
C
         IF( .NOT. NOCNCT .OR. TOGTHR ) THEN
C
            IF( TOGTHR ) THEN
               IF( NOCNCT ) THEN
                  ewrite(1,'(A,I3,1PE13.5)')
     :              '--- Starting node move adapt of sweep',I*2-1,
     :              ANGLEFACT
               ELSE
                  ewrite(1,'(A,I3,1PE13.5)')
     :              '--- Starting connect & move adapt of sweep',I*2-1,
     :              ANGLEFACT
               END IF
            ELSE
               ewrite(1,'(A,I3,1PE13.5)')
     :              '--- Starting connectivity adapt of sweep',I*2-1,
     :              ANGLEFACT
            END IF
C
C - an adaptive sweep with no node movement (element division is also
C - inactive at the moment)
C
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
            CALL ADAPT1( INTARR(BIGLST), RLARR(NODLST),
     :                   ENLBAS, ENLIST,
     :                   INTARR(NELBAS), INTARR(NELIST), INTARR(EELIST),
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   EDSPLT, MSHOPT(2), MSHOPT(3), MSHOPT(4),
     :                   ELSPLT, TOGTHR .AND. MSHOPT(6),
     :                   FLCHK3, FLCHK1, ALWAYS, DOTOP, MINCHG, CHANG1 )
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
C
            IF( IERR .NE. 0 ) THEN
               WRITE(0,'(a)') '*** ADPTVY: ERROR FROM FIRST ADAPT1'
               NWNNOD = -1
               GOTO 999
            END IF
C
         END IF
C
c         goto 600
C
         IF( TOGTHR .OR. .NOT. MSHOPT(6) ) THEN
C
           CHANG2 = 0.0
C
         ELSE
C
            ewrite(1,'(A,I3,1PE13.5)')
     :           '--- Starting node move adapt of sweep',I*2-1,
     :              ANGLEFACT
C
C - an adaptive sweep with only node movement
C
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
            CALL ADAPT1( INTARR(BIGLST), RLARR(NODLST),
     :                   ENLBAS, ENLIST,
     :                   INTARR(NELBAS), INTARR(NELIST), INTARR(EELIST),
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.,
     :                   .TRUE.,
     :                   FLCHK5, FLCHK2, ALWAYS, DOTOP, MINCHG, CHANG2 )
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
C
            IF( IERR .NE. 0 ) THEN
              WRITE(0,'(a)') '*** ADPTVY: ERROR FROM SECOND ADAPT1'
              NWNNOD = -1
              GOTO 999
            END IF
C
         END IF
C
         FLCHK1 = 15
         FLCHK2 = 17
         FLCHK3 = 14
C
         IF( CHANG1 + CHANG2 .EQ. 0.0 ) THEN
c            IF( ANGLEFACT .GT. 0.0 ) GOTO 600
c            ANGLEFACT = 0.1
c            FLCHK3 = 0
             goto 600
         END IF
         IF( I*2 .GT. NSWPS ) GOTO 600
C
         IF( .NOT. NOCNCT .OR. TOGTHR ) THEN
C
            IF( TOGTHR ) THEN
               IF( NOCNCT ) THEN
                  ewrite(1,'(A,I3,1PE13.5)')
     :              '--- Starting node move adapt of sweep',I*2,
     :              ANGLEFACT
               ELSE
                  ewrite(1,'(A,I3,1PE13.5)')
     :              '--- Starting connect & move adapt of sweep',I*2,
     :              ANGLEFACT
               END IF
            ELSE
               ewrite(1,'(A,I3,1PE13.5)')
     :              '--- Starting connectivity adapt of sweep',I*2,
     :              ANGLEFACT
            END IF
C
C - an adaptive sweep with no node movement (element division is also
C - inactive at the moment)
C
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
            CALL ADAPT1( INTARR(BIGLST), RLARR(NODLST),
     :                   ENLBAS, ENLIST,
     :                   INTARR(NELBAS), INTARR(NELIST), INTARR(EELIST),
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   EDSPLT, MSHOPT(2), MSHOPT(3), MSHOPT(4),
     :                   ELSPLT, TOGTHR .AND. MSHOPT(6),
     :                   FLCHK4, FLCHK3, ALWAYS, DOTOP, MINCHG, CHANG1 )
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
C
            IF( IERR .NE. 0 ) THEN
               WRITE(0,'(a)') '*** ADPTVY: ERROR FROM THIRD ADAPT1'
               NWNNOD = -1
               GOTO 999
            END IF
C
         END IF
C
         IF( TOGTHR .OR. .NOT. MSHOPT(6) ) THEN
C
            CHANG2 = 0.0
C
         ELSE
C
            ewrite(1,'(A,I3,1PE13.5)')
     :            '--- Starting node move adapt of sweep',I*2,
     :              ANGLEFACT
C
C - an adaptive sweep with only node movement
C
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
            CALL ADAPT1( INTARR(BIGLST), RLARR(NODLST),
     :                   ENLBAS, ENLIST,
     :                   INTARR(NELBAS), INTARR(NELIST), INTARR(EELIST),
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.,
     :                   .TRUE.,
     :                   FLCHK6, FLCHK5, ALWAYS, DOTOP, MINCHG, CHANG2 )
            IF( I .LT. NSWPS/2 ) MINCHG = -MINCHG
C
            IF( IERR .NE. 0 ) THEN
              WRITE(0,'(a)') '*** ADPTVY: ERROR FROM FOURTH ADAPT1'
              NWNNOD = -1
              GOTO 999
            END IF
C
         END IF
C
         IF( CHANG1 + CHANG2 .EQ. 0.0 ) THEN
c            IF( ANGLEFACT .GT. 0.0 ) GOTO 600
c            ANGLEFACT = 0.1
c            FLCHK1 = 0
            goto 600
         END IF
         IF( I*2 .GE. NSWPS ) GOTO 600
C
      END DO
C
 600  CONTINUE
C
      IF( TWOSTG .AND. .NOT. DUNSTG ) THEN
         EDSPLT = MSHOPT(1)
c - MSHOPT(5) (element splitting) is currently ignored!!
c         ELSPLT = MSHOPT(5)
         ELSPLT = .FALSE.
         ewrite(1,'(a)')
     :       '--- ADPTVY: STARTING STAGE TWO WITH REFINEMENT ON'
         DUNSTG = .TRUE.
         NOCNCT = .NOT. ( EDSPLT .OR. MSHOPT(2) .OR. MSHOPT(3) .OR.
     :                    ELSPLT .OR. MSHOPT(4) )
         GOTO 10
      END IF
C
c      CALL CHKGMY( INTARR(BIGLST), RLARR(NODLST) )
C
C - This returns the size of the new adapted mesh; the new number of
C - nodes (NWNNOD), the new number of elements (NWNELM) and the new size
C - of the element node list (NWSZLS - which is always 4*NWNELM for 3D)
C - and the size of the node-to-node (node connectivity) list (NWSZNN).
C - It also returns the size of the new surface element node list (NWSZSN).
C
      CALL GETSIZ( INTARR(BIGLST), RLARR(NODLST),
     :             NWNNOD, NWNELM, NWSZEN, NWSZNN, NWNSEL, NWSZSN )
C
      ewrite(1,'(a)') 
     :     '--- ADPTVY: FINISHED ADAPTING MESH'
      ewrite(1,'(a,2i8)') 
     :     'TOPNOD,TOPBIG:',TOPNOD,TOPBIG
      ewrite(1,'(a,3i8)') 
     :     'NODS,ELEMS,SURFS:',NWNNOD,NWNELM,NWNSEL
C
      NWNODX = NODLST + (TOPNOD+1)*SIZNOD
      NWNODY = NWNODX + NWNNOD
      NWNODZ = NWNODY + NWNNOD
      NEWFLD = NWNODZ + NWNNOD
C
      IF( NEWFLD + NWNNOD*TOTFRE .GE. RLSIZ ) THEN
         WRITE(0,'(a)')
     :        '*** ADPTVY: REAL ARRAY TOO SMALL TO INTERPOLATE'
         NWNNOD = -1
         IERR = -5008
         GOTO 999
      END IF
C
C - This interpolates, from the old mesh into NEWFLD, the NFIELD fields in
C - OLDFLD with NFREE(I) components in the Ith field. TOTFRE is the sum
C - of NFREE(1) to NFREE(NFIELD). This allows a long vector with moments,
C - nodes and groups to be interpolated in one call by sending TOTFRE as
C - moments*groups, NFIELD as no. of groups and each NFREE(I) as the
C - number of moments.
C - OLDFLD has length TOTM*NNOD and NEWFLD should have length TOTM*NEWNND.
C - INTFLD can be called as many times as necessary for appropriate fields.
C
      ewrite(1,'(a)') 
     :     '--- ADPTVY: Interpolating fields...'
c      pause
      CALL INTFLD( INTARR(BIGLST), RLARR(NODLST),
     :             NODX, NODY, NODZ, OLDFLD,
     :             NFREE, NFIELD, TOTFRE,
     :             ENLBAS, ENLIST,
     :             INTARR(NELBAS), INTARR(NELIST), INTARR(EELIST),
     :             SZENLS, SZNELS, NNOD, NELM, NHALO, RLARR(NEWFLD) )
c
      if( debug ) then
         call chkint( intarr(biglst), rlarr(nodlst),
     :                nodx, nody, nodz, oldfld, nnod,
     :                enlbas, enlist, szenls,
     :                nelm, rlarr(newfld) )
      end if
C
C - At this point the old mesh info is no longer needed.
C
      IF( NEWMTX .LT. 0 ) THEN
         SIZMTX = 1
      ELSE
         SIZMTX = NWNNOD*9
      END IF
      NEWMTX = NEWFLD + NWNNOD*TOTFRE
C
      IF( NEWMTX + SIZMTX .GE. RLSIZ ) THEN
         WRITE(0,'(a)')
     :      '*** ADPTVY: REAL ARRAY TOO SMALL FOR NEW METRICS'
         NWNNOD = -1
         IERR = -5009
         GOTO 999
      END IF
C
      IF( NWNDLC .GE. 0 ) THEN
         SIZLOC = NWNNOD
      ELSE
         SIZLOC = 1
      ENDIF
C
      IF( NWSROW .GE. 0 ) THEN
         SIZOWN = NWNSEL
      ELSE
         SIZOWN = 1
      END IF
C
      NWNDLC = BIGLST + (TOPBIG+1)*SIZBIG
      NWSROW = NWNDLC + SIZLOC
      NWELRG = NWSROW + SIZOWN
      NWENLB = NWELRG + NWNELM
      NWENLS = NWENLB + NWNELM + 1
      NWSNLB = NWENLS + NWSZEN
      NWSNLS = NWSNLB + NWNSEL + 1
      NWSFID = NWSNLS + NWSZSN
C
C - This forms the new element node list (NWENLS) and its base pointer
C - (NWENLB) for the adapted mesh, as well as returning the element
C - material/region id in NWELRG for each element. The nodal positions
C - for the adapted mesh is returned in NODX, NODY and NODZ, and the
C - interpolated error metrics are returned in NEWMTX if SIZMTX is sent
C - as 9*NWNNOD (otherwise SIZMTX should be sent as 1).
C
      ewrite(1,'(a)') 
     :     '--- ADPTVY: Forming new fixed mesh data...'
C
      IF( NWSFID + NWNSEL .GE. INTSIZ ) THEN
c         WRITE(0,'(a)')
c     :        '*** ADPTVY: INTEGER ARRAY TOO SMALL FOR NEW MESH'
c         NWNNOD = -1
c         IERR = -5010
c         GOTO 999
c      END IF
C
        ewrite(0, '(a)') 
     :  '*** ADPTVY: INTEGER ARRAY TOO SMALL FOR NEW MESH, REALLOCATING'
        allocate( tmparr(nwsfid+nwnsel-nwndlc) )
C
        NWNDLC = 1
        NWSROW = NWNDLC + SIZLOC
        NWELRG = NWSROW + SIZOWN
        NWENLB = NWELRG + NWNELM
        NWENLS = NWENLB + NWNELM + 1
        NWSNLB = NWENLS + NWSZEN
        NWSNLS = NWSNLB + NWNSEL + 1
        NWSFID = NWSNLS + NWSZSN
c
      CALL MKFIXD( INTARR(BIGLST), RLARR(NODLST), .FALSE.,
     :             RLARR(NWNODX), RLARR(NWNODY), RLARR(NWNODZ),
     :             RLARR(NEWMTX), SIZMTX, PRDNDS, NPRDND,
     :             TMPARR(NWNDLC), SIZLOC, TMPARR(NWSROW), SIZOWN,
     :             TMPARR(NWELRG), TMPARR(NWENLB), TMPARR(NWENLS),
     :             TMPARR(NWSNLB), TMPARR(NWSNLS), TMPARR(NWSFID),
     :             NWNELM, NWNNOD, NWNSEL, NWSZEN, NWSZSN, GEOM3D,
     :             GATHER, SCATER, NGATH,  NHALO,
     :             ATOSEN, ATOREC, NPROC )
C

        NWNDLC = BIGLST
        NWSROW = NWNDLC + SIZLOC
        NWELRG = NWSROW + SIZOWN
        NWENLB = NWELRG + NWNELM
        NWENLS = NWENLB + NWNELM + 1
        NWSNLB = NWENLS + NWSZEN
        NWSNLS = NWSNLB + NWNSEL + 1
        NWSFID = NWSNLS + NWSZSN

c
c - this overwrites the adaptive mesh info in INTARR, but hopefully
c - nothing will use this anyway once we leave adptvy
c
        intarr(biglst:nwsfid+nwnsel-1) = tmparr(:)
C
        deallocate( tmparr )
c
      ELSE
C
      CALL MKFIXD( INTARR(BIGLST), RLARR(NODLST), .FALSE.,
     :             RLARR(NWNODX), RLARR(NWNODY), RLARR(NWNODZ),
     :             RLARR(NEWMTX), SIZMTX, PRDNDS, NPRDND,
     :             INTARR(NWNDLC), SIZLOC, INTARR(NWSROW), SIZOWN,
     :             INTARR(NWELRG), INTARR(NWENLB), INTARR(NWENLS),
     :             INTARR(NWSNLB), INTARR(NWSNLS), INTARR(NWSFID),
     :             NWNELM, NWNNOD, NWNSEL, NWSZEN, NWSZSN, GEOM3D,
     :             GATHER, SCATER, NGATH,  NHALO,
     :             ATOSEN, ATOREC, NPROC )
C
      END IF
C
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,'(a)') '*** ADPTVY: GOT ERROR FROM MKFIXD'
         NWNNOD = -1
         GOTO 999
      END IF
c
      if( debug ) then
c
         intsum = 0
         do i = 1, nwnelm+1
            intsum = intsum + intarr(nwenlb+i-1)
         end do
         ewrite(3,'(a,i8)') 
     :        'Sum of ENLBAS: ',intsum
c
         intsum = 0
         do i = 1, nwszen
            intsum = intsum + intarr(nwenls+i-1)
         end do
         ewrite(3,'(a,2i8)') 
     :        'Sum of ENLIST: ',intsum,nwszen
c
      end if

C     Finalize AdaptProgress module
      call finalize(NWNNOD)

      ewrite(1,'(a)') 
     :     '--- ADPTVY: Finished successfully'
c
  999 IF( IERR .NE. 0 ) THEN
         write(0,'(a,i7)') '*** ADPTVY: FINISHED WITH ERROR',IERR
      END IF
C
      RETURN
      END
C
C
C     
      LOGICAL FUNCTION ADAPT_FALSE()
      ADAPT_FALSE = .FALSE.
      RETURN
      END
C
C
      LOGICAL FUNCTION ADAPT_TRUE()
      ADAPT_TRUE = .TRUE.
      RETURN
      END
C
C
C
      INTEGER FUNCTION MESH_IS_3D( NEWVAL )
C
      IMPLICIT NONE
C
      INTEGER NEWVAL
C
      INCLUDE 'blknew.i'
C
      IF( IS3DMS ) THEN
        MESH_IS_3D = 1
      ELSE
        MESH_IS_3D = 0
      END IF
C
      IF( NEWVAL .GT. 0 ) THEN
        IS3DMS = .TRUE.
      ELSE IF( NEWVAL .EQ. 0 ) THEN
        IS3DMS = .FALSE.
      END IF
C
      RETURN
      END
C
C
      INTEGER FUNCTION DO_CONSISTENCY_CHECK( NEWVAL )
C
      IMPLICIT NONE
C
      INTEGER NEWVAL
C
      INCLUDE 'blkerr.i'
C
      IF( CNSYCH ) THEN
        DO_CONSISTENCY_CHECK = 1
      ELSE
        DO_CONSISTENCY_CHECK = 0
      END IF
C
      IF( NEWVAL .GT. 0 ) THEN
        CNSYCH = .TRUE.
      ELSE IF( NEWVAL .EQ. 0 ) THEN
        CNSYCH = .FALSE.
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE SET_ADAPT_DEBUG( NEWVAL )
      use write_log
      IMPLICIT NONE
C
      INTEGER NEWVAL
C
      INCLUDE 'blkerr.i'
C
      IF( NEWVAL .GT. 0 ) THEN
        DEBUG = .TRUE.
        ewrite(1,*) '+++ Debugging switched on'
      ELSE IF( NEWVAL .EQ. 0 ) THEN
        DEBUG = .FALSE.
        ewrite(1,*) '+++ Debugging switched off'
      END IF
C
      RETURN
      END
C
C
      INTEGER FUNCTION GET_ADAPT_ERR()
C
      IMPLICIT NONE
C
      INCLUDE 'blkerr.i'
C
      GET_ADAPT_ERR = IERR
C
      RETURN
      END
C
      SUBROUTINE SET_BIGLST_SIZE( NUMINTS )
      use write_log
      IMPLICIT NONE
C
      INTEGER NUMINTS
C
      INCLUDE 'blknew.i'
      INCLUDE 'blkprm.i'
      INCLUDE 'blkerr.i'
C
      MAXBIG = NUMINTS/SIZBIG
C
      ewrite(1,*) 
     :     '--- Set aside space for about ',MAXBIG*2/9,
     :     ' elements'
      RETURN
      END
C
C
      SUBROUTINE SET_NODLST_SIZE( NUMRLS )
      use write_log
      IMPLICIT NONE
C
      INTEGER NUMRLS
C
      INCLUDE 'blknew.i'
      INCLUDE 'blkprm.i'
      INCLUDE 'blkerr.i'
C
      MXNODS = NUMRLS/SIZNOD
      ewrite(1,*) 
     :     '--- Set aside space for ',MXNODS,' nodes'
C
      RETURN
      END
C
