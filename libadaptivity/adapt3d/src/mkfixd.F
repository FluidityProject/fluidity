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
      SUBROUTINE MKFIXD( BIGLST, NODLST, NDNMBR,
     :                   NODX, NODY, NODZ, ORGMTX, SIZMTX,
     :                   PRDNDS, NPRDND,
     :                   NODLOC, SIZLOC, SUROWN, SIZOWN,
     :                   ELMREG, ENLBAS, ENLIST,
     :                   SNLBAS, SNLIST, SURFID,
     :                   NELM, NNOD, NSELM, SZENLS, SZSNLS, GEOM3D,
     :                   GATHER, SCATER, NGATH, NHALO,
     :                   ATOSEN, ATOREC, NPROC )
C-----------------------------------------------------------------------
C
C - This subroutine forms the element node list (ENLIST) and its base
C - pointer (ENLBAS) for the new adapted mesh, as well as returning the
C - element material/region id in ELMREG for each element. The nodal
C - positions for the adapted mesh is returned in NODX, NODY and NODZ,
C - and the interpolated error metrics are returned in ORGMTX if SIZMTX
C - is sent as 9*NNOD (otherwise it should be sent as 1).
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INTEGER NELM, SZENLS, NNOD, SZSNLS, NSELM, SIZOWN, SIZLOC, SIZMTX,
     :        NHALO, NPROC, NGATH, NPRDND
C
      INTEGER ELMREG(NELM), ENLIST(SZENLS), ENLBAS(NELM+1), 
     :        SNLBAS(NSELM+1), SNLIST(SZSNLS), SUROWN(SIZOWN), 
     :        NODLOC(SIZLOC), GATHER(NGATH), SCATER(NHALO),
     :        ATOSEN(NPROC), ATOREC(NPROC), PRDNDS(NPRDND),
     :        SURFID(NSELM)
C
      REAL NODX(NNOD), NODY(NNOD), NODZ(NNOD), ORGMTX(SIZMTX)
C
      LOGICAL GEOM3D, NDNMBR
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
      INTEGER IPTR, CONELM(4), EDGLST(6), INOD, IELM, IFLAG, IREG,
     :        INEW, NXT, LST, NODS(4), IPOS, I, J, K, SPOS, SELM,
     :        NODCON(4), NXTC, SUM, CNTGMY, CNTINT, CNTSPL, JNOD, KNOD,
     :        minelv, maxelv, minelr, maxela, i1, i2
C
      REAL X(4), Y(4), Z(4), TETVOL, avol, arad, asp, m(9), areas(4),
     :     minvol, minrad, minrvl, minvrd, maxasp, maxavl, maxard,
     :     minras, minvas, L(6), maxvol, maxvrd, maxvas
      integer mxsfid, nodcnt, nusfid
C      
      integer, dimension(:), allocatable:: srfcnt
      real, dimension(:), allocatable:: srfavx, srfavy, srfavz
      integer, parameter:: MAX_SFID_STATS=1000
C
      REAL*8 VOLSUM, ASPAVE
C
      LOGICAL GTNDFL, FACS(4), getflg
C
      IF( IERR .NE. 0 ) RETURN
C
c      print*,'NNOD,NELM,NSELM,GEOM3D: ',NNOD,NELM,NSELM,IS3DMS
      IF( GEOM3D .NEQV. IS3DMS ) THEN
        WRITE(0,*) '*** MKFIXD: Inconsistent GEOM3D!'
        IF( GEOM3D ) THEN
          WRITE(0,*) "   Didn't do 3D adapt before, but GEOM3D now true"
        ELSE
          WRITE(0,*) '   Did 3D adapt before, but GEOM3D now false'
        END IF
        IERR = -7835
        RETURN
      END IF
c      print*,'SIZMTX: ',SIZMTX
C
      CNTGMY = 0
      CNTSPL = 0
      CNTINT = 0
      IPTR = STTNOD
      INOD = 0
      JNOD = NNOD + 1
      NODCNT = 0
C
  10  IF( GTNDFL(BIGLST,NODLST,IPTR,1) ) THEN
C
         if( gtndfl(biglst,nodlst,iptr,2) ) then
            cntgmy = cntgmy + 1
         else if( gtndfl(biglst,nodlst,iptr,4) ) then
            cntspl = cntspl + 1
         else if( gtndfl(biglst,nodlst,iptr,3) ) then
            cntint = cntint + 1
         end if
c
         nodcnt = nodcnt + 1
         if(nodcnt.gt.nnod) then
            WRITE(0,*) '*** MKFIXD: Too many nodes!'
            WRITE(0,*) nodcnt,nnod,iptr,topnod
            ierr = -9797
            return
         end if
C
         IF( NDNMBR ) THEN
            IF( GTNDFL(BIGLST,NODLST,IPTR,8) ) THEN
               JNOD = JNOD - 1
               KNOD = JNOD
            ELSE
               INOD = INOD + 1
               KNOD = INOD
            END IF
         ELSE
            KNOD = INT(NODLST(17,IPTR))
            if( knod .le. 0 ) then
               WRITE(0,*) '*** MKFIXD: Got node again!'
               WRITE(0,*) knod,iptr,inod,jnod,nnod
               ierr = -1987
               return
            else if( knod .gt. nnod ) then
               WRITE(0,*)'*** MKFIXD: Got out-of-range KNOD: ',knod,nnod
               ierr = -1987
               return
            end if
         END IF
C
C - rescale and shift the domain back to what it was
C
         NODX(KNOD) = NODLST(1,IPTR)*XSCALE + XSHIFT
         NODY(KNOD) = NODLST(2,IPTR)*YSCALE + YSHIFT
         NODZ(KNOD) = NODLST(3,IPTR)*ZSCALE + ZSHIFT
C
         IF( SIZMTX .EQ. 9*NNOD ) THEN
            I = KNOD*9-9
            ORGMTX(I+1) = NODLST( 7,IPTR)/XSCALE/XSCALE
            ORGMTX(I+2) = NODLST( 8,IPTR)/XSCALE/YSCALE
            ORGMTX(I+3) = NODLST( 9,IPTR)/XSCALE/ZSCALE
            ORGMTX(I+4) = NODLST(10,IPTR)/YSCALE/XSCALE
            ORGMTX(I+5) = NODLST(11,IPTR)/YSCALE/YSCALE
            ORGMTX(I+6) = NODLST(12,IPTR)/YSCALE/ZSCALE
            ORGMTX(I+7) = NODLST(13,IPTR)/ZSCALE/XSCALE
            ORGMTX(I+8) = NODLST(14,IPTR)/ZSCALE/YSCALE
            ORGMTX(I+9) = NODLST(15,IPTR)/ZSCALE/ZSCALE
C            DO I = 1, 9
C               ORGMTX(I+KNOD*9-9) = NODLST(6+I,IPTR)
C            END DO
         END IF
C
         NODLST(17,IPTR) = -FLOAT(KNOD)
C
         IF( SIZLOC .EQ. NNOD ) NODLOC(KNOD) = INT(NODLST(16,IPTR))
c         print*,'Node ',nodx(Knod),nody(Knod),nodz(Knod),iptr
         if(nodlst(nxtnod,iptr).eq.iptr) then
            WRITE(0,*) '*** MKFIXD: Node points to itself!'
            WRITE(0,*) iptr,inod,jnod,knod,nnod
            ierr = -1876
            return
         end if
         IPTR = NODLST(NXTNOD,IPTR)
         if(iptr.le.topnod .and. iptr.gt.0) goto 10
C
         if(nodcnt.lt.nnod) then
            WRITE(0,*) '*** MKFIXD: Reached end of node list too early'
            WRITE(0,*) nodcnt,nnod,iptr,topnod
            ierr = -9696
            return
         end if
C
      END IF
C
      IF( IERR .NE. 0 ) THEN
         WRITE(0,*) '*** MKFIXD: GOT ERROR FROM GTNDFL'
         RETURN
      END IF
c
      iptr = sttnod
      inod = 0
  11  if( gtndfl(biglst,nodlst,iptr,1) ) then
         inod = inod + 1
         knod = int(nodlst(17,iptr))
         if(knod.ge.0) then
            WRITE(0,*) '*** MKFIXD: Found +ve entry!'
            WRITE(0,*) nodlst(17,iptr),iptr,inod,nnod
            ierr = -17856
            return
         end if
         if(inod.gt.nnod) then
            WRITE(0,*) '*** MKFIXD: Beyond end of nodes!'
            WRITE(0,*) nodlst(17,iptr),iptr,inod,nnod
            ierr = -16847
            return
         end if
         knod = -knod
         if(knod.gt.nnod) then
            WRITE(0,*) '*** MKFIXD: Node number too big!'
            WRITE(0,*) nodlst(17,iptr),iptr,inod,nnod
            ierr = -13968
            return
         end if
         nodlst(17,iptr) = float(knod)
         iptr = nodlst(nxtnod,iptr)
         if(iptr.le.topnod .and. iptr.gt.0) goto 11
      end if
C
      IF( NPRDND .GT. 1 ) THEN
C
         DO I = 1, NPRDND
            j = prdnds(i)
            PRDNDS(I) = INT( NODLST(17,PRDNDS(I)) )
            ewrite(2,*) 
     :           '  periodic renumber: ',j,prdnds(i)
         END DO
C
      END IF
c
      ewrite(2,*) 
     :     'Finished new node data'
      ewrite(2,*) 
     :     'Geom,split,int: ',cntgmy,cntspl,cntint
C
      CNTINT = 0
      CNTGMY = 0
      IPTR = STTBIG
      IELM = 0
      IPOS = 0
      SELM = 0
      SPOS = 0
      ENLBAS(IELM+1) = IPOS
      SNLBAS(SELM+1) = SPOS
      mxsfid = 0
      allocate( srfcnt(0:MAX_SFID_STATS) )
      allocate( srfavx(0:MAX_SFID_STATS) )
      allocate( srfavy(0:MAX_SFID_STATS) )
      allocate( srfavz(0:MAX_SFID_STATS) )
      srfavx = 0.0
      srfavy = 0.0
      srfavz = 0.0
      srfcnt = 0
c      ewrite(3,*) 'cleared srfavs'
C
  20  IF( BIGLST(4,IPTR) .GT. 0 ) THEN
C
         IELM = IELM + 1
C
         IF( IELM .GT. NELM ) THEN
            WRITE(0,*) '*** MKFIXD: ERROR COUNTING ELEMENTS'
            WRITE(0,*) NELM,IELM,IPTR
            IERR = -8989
            deallocate(srfcnt, srfavx, srfavy, srfavz)
            RETURN
         END IF
C
         CALL ELMINF( BIGLST, NODLST, IPTR, CONELM, EDGLST,
     :                INEW, IFLAG, IREG, NXT, LST )
C
         ELMREG(IELM) = IREG
C
         CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, .TRUE. )
C
         i1 = 1
         i2 = 2
         IF( IS3DMS ) THEN
            X(1) = NODLST(1,NODS(1))
            X(2) = NODLST(1,NODS(2))
            X(3) = NODLST(1,NODS(3))
            X(4) = NODLST(1,NODS(4))
            Y(1) = NODLST(2,NODS(1))
            Y(2) = NODLST(2,NODS(2))
            Y(3) = NODLST(2,NODS(3))
            Y(4) = NODLST(2,NODS(4))
            Z(1) = NODLST(3,NODS(1))
            Z(2) = NODLST(3,NODS(2))
            Z(3) = NODLST(3,NODS(3))
            Z(4) = NODLST(3,NODS(4))
            IF( TETVOL( X, Y, Z ) .LT. 0.0 ) THEN
c               print*,'*** MKFIXD: Got inside-out tet...'
c               STOP
c               print*,ielm,iptr
c               print*,nods
c               print*,x
c               print*,y
c               print*,z
               i1 = 2
               i2 = 1
            END IF
         END IF
C
         ENLIST(IPOS+1) = NODLST(17,NODS(i1))
         ENLIST(IPOS+2) = NODLST(17,NODS(i2))
         ENLIST(IPOS+3) = NODLST(17,NODS(3))
C
         IPOS = IPOS + 3
         IF( IS3DMS ) THEN
            ENLIST(IPOS+1) = NODLST(17,NODS(4))
            IPOS = IPOS + 1
         END IF
         if( ielm .lt. 11 .and. debug ) then
            print*,'old elm ',iptr,(nods(i),i=1,4)
            print*,'elm ',ielm,(enlist(ipos+i-4),i=1,4)
         end if
C
         IF( IPOS .GT. SZENLS+1 ) THEN
            WRITE(0,*) '*** MKFIXD: BEYOND END OF ENLIST'
            WRITE(0,*) SZENLS,IPOS,IELM,NELM,IPTR
            IERR = -9898
            deallocate(srfcnt, srfavx, srfavy, srfavz)
            RETURN
         END IF
C
         ENLBAS(IELM+1) = IPOS
C
         IF( IS3DMS ) THEN
C
c            IF( CONELM(1) .LE. 0 .OR. CONELM(2) .LE. 0 .OR.
c     :          CONELM(3) .LE. 0 .OR. CONELM(4) .EQ. 0 ) THEN
C
               IF( CONELM(1) .LE. 0 ) THEN
                  SELM = SELM + 1
                  SNLIST(SPOS+i1) = NODLST(17,NODS(2))
                  SNLIST(SPOS+i2) = NODLST(17,NODS(1))
                  SNLIST(SPOS+3)  = NODLST(17,NODS(3))
c                  if( selm .lt. 11 ) then
c                     print*,'srf ',selm,(snlist(spos+i),i=1,3)
c                  end if
                  SPOS = SPOS + 3
                  SNLBAS(SELM+1) = SPOS
                  SURFID(SELM) = -CONELM(1)-1
                  IF( SIZOWN .EQ. NSELM ) SUROWN(SELM) = IELM
                  if( ielm.gt.nelm .or. ielm.le.0 )
     :                WRITE(0,*) 'MKFIXD: IELM too big!',ielm,nelm,selm
                  mxsfid = max(mxsfid,surfid(selm))
                  if( surfid(selm) .lt. 0 ) then
                     WRITE(0,*) 'BAD SURFID RANGE: ',surfid(selm),selm
                  else if( surfid(selm) .le. MAX_SFID_STATS ) then
c                     if(surfid(selm).gt.0) elmreg(ielm) = 3
                     srfcnt(surfid(selm))=srfcnt(surfid(selm))+1
                     srfavx(surfid(selm))=srfavx(surfid(selm))
     :                                   +x(1)+x(2)+x(3)
                     srfavy(surfid(selm))=srfavy(surfid(selm))
     :                                   +y(1)+y(2)+y(3)
                     srfavz(surfid(selm))=srfavz(surfid(selm))
     :                                   +z(1)+z(2)+z(3)
                  end if
               END IF
C
               IF( CONELM(2) .LE. 0 ) THEN
                  SELM = SELM + 1
                  SNLIST(SPOS+i1) = NODLST(17,NODS(1))
                  SNLIST(SPOS+i2) = NODLST(17,NODS(2))
                  SNLIST(SPOS+3)  = NODLST(17,NODS(4))
c                  if( selm .lt. 11 ) then
c                     print*,'srf ',selm,(snlist(spos+i),i=1,3)
c                  end if
                  SPOS = SPOS + 3
                  SNLBAS(SELM+1) = SPOS
                  SURFID(SELM) = -CONELM(2)-1
                  IF( SIZOWN .EQ. NSELM ) SUROWN(SELM) = IELM
                  if( ielm.gt.nelm .or. ielm.le.0 )
     :                WRITE(0,*) 'MKFIXD: IELM too big!',ielm,nelm,selm
                  mxsfid = max(mxsfid,surfid(selm))
                  if( surfid(selm) .lt. 0 ) then
                     WRITE(0,*) 'BAD SURFID RANGE: ',surfid(selm),selm
                  else if( surfid(selm) .le. MAX_SFID_STATS ) then
c                     if(surfid(selm).gt.0) elmreg(ielm) = 3
                     srfcnt(surfid(selm))=srfcnt(surfid(selm))+1
                     srfavx(surfid(selm))=srfavx(surfid(selm))
     :                                   +x(1)+x(2)+x(3)
                     srfavy(surfid(selm))=srfavy(surfid(selm))
     :                                   +y(1)+y(2)+y(3)
                     srfavz(surfid(selm))=srfavz(surfid(selm))
     :                                   +z(1)+z(2)+z(3)
                  end if
               END IF
C
               IF( CONELM(3) .LE. 0 ) THEN
                  SELM = SELM + 1
                  SNLIST(SPOS+i1) = NODLST(17,NODS(3))
                  SNLIST(SPOS+i2) = NODLST(17,NODS(1))
                  SNLIST(SPOS+3)  = NODLST(17,NODS(4))
c                  if( selm .lt. 11 ) then
c                     print*,'srf ',selm,(snlist(spos+i),i=1,3)
c                  end if
                  SPOS = SPOS + 3
                  SNLBAS(SELM+1) = SPOS
                  SURFID(SELM) = -CONELM(3)-1
                  IF( SIZOWN .EQ. NSELM ) SUROWN(SELM) = IELM
                  if( ielm.gt.nelm .or. ielm.le.0 )
     :                WRITE(0,*) 'MKFIXD: IELM too big!',ielm,nelm,selm
                  mxsfid = max(mxsfid,surfid(selm))
                  if( surfid(selm) .lt. 0 ) then
                     WRITE(0,*) 'BAD SURFID RANGE: ',surfid(selm),selm
                  else if( surfid(selm) .le. MAX_SFID_STATS ) then
c                     if(surfid(selm).gt.0) elmreg(ielm) = 3
                     srfcnt(surfid(selm))=srfcnt(surfid(selm))+1
                     srfavx(surfid(selm))=srfavx(surfid(selm))
     :                                   +x(1)+x(2)+x(3)
                     srfavy(surfid(selm))=srfavy(surfid(selm))
     :                                   +y(1)+y(2)+y(3)
                     srfavz(surfid(selm))=srfavz(surfid(selm))
     :                                   +z(1)+z(2)+z(3)
                  end if
               END IF
C
               IF( CONELM(4) .LE. 0 ) THEN
                  SELM = SELM + 1
c                  if(facs(4).eq.4) then
                     WRITE(0,*) '+++ MKFIXD: Eh??? Surface on face 4!'
                     WRITE(0,*) ielm,selm
                     WRITE(0,*) conelm
                     WRITE(0,*) nods
c                  end if
                  SNLIST(SPOS+i1) = NODLST(17,NODS(2))
                  SNLIST(SPOS+i2) = NODLST(17,NODS(3))
                  SNLIST(SPOS+3)  = NODLST(17,NODS(4))
c                  if( selm .lt. 11 ) then
c                     print*,'srf ',selm,(snlist(spos+i),i=1,3)
c                  end if
                  SPOS = SPOS + 3
                  SNLBAS(SELM+1) = SPOS
                  SURFID(SELM) = -CONELM(4)-1
                  IF( SIZOWN .EQ. NSELM ) SUROWN(SELM) = IELM
                  if( ielm.gt.nelm .or. ielm.le.0 )
     :                WRITE(0,*) 'MKFIXD: IELM too big!',ielm,nelm,selm
                  mxsfid = max(mxsfid,surfid(selm))
                  if( surfid(selm) .lt. 0 ) then
                     WRITE(0,*) 'BAD SURFID RANGE: ',surfid(selm),selm
                  else if( surfid(selm) .le. MAX_SFID_STATS ) then
c                     if(surfid(selm).gt.0) elmreg(ielm) = 3
                     srfcnt(surfid(selm))=srfcnt(surfid(selm))+1
                     srfavx(surfid(selm))=srfavx(surfid(selm))
     :                                   +x(1)+x(2)+x(3)
                     srfavy(surfid(selm))=srfavy(surfid(selm))
     :                                   +y(1)+y(2)+y(3)
                     srfavz(surfid(selm))=srfavz(surfid(selm))
     :                                   +z(1)+z(2)+z(3)
                  end if
               END IF
C
               IF( SELM .GT. NSELM ) THEN
                 WRITE(0,*)'*** MKFIXD: ERROR COUNTING SURFACE ELEMENTS'
                 WRITE(0,*) NSELM,SELM
                 STOP
               END IF
C
               IF( SPOS .GT. SZSNLS+1 ) THEN
                  WRITE(0,*) '*** MKFIXD: BEYOND END OF SNLIST'
                  WRITE(0,*) SZSNLS,SPOS
                  STOP
               END IF
C
c            END IF
C
C         ELSE
C
C
         END IF
C
         IPTR = NXT
         IF( IPTR .LE. TOPBIG .AND. IPTR .GT. 0 ) GOTO 20
C
         IF( IELM .NE. NELM ) THEN
           WRITE(0,*)'*** MKFIXD: Reached end of element list too early'
           WRITE(0,*) ielm,nelm,iptr,topbig
           ierr = -9595
           deallocate(srfcnt, srfavx, srfavy, srfavz)
           return
         end if
C
      ELSE IF( BIGLST(4,IPTR) .LT. 0 ) THEN
C
C         CALL EDGINF( BIGLST, NODLST, IPTR, NODS,
C     :                INEW, IFLAG, NXT, LST )
C
         if( getflg(biglst,nodlst,iptr,2) ) then
            cntgmy = cntgmy + 1
         else if( getflg(biglst,nodlst,iptr,3) ) then
            cntint = cntint + 1
         end if
C
         IPTR = BIGLST( NXTBIG, IPTR )
         IF( IPTR .LE. TOPBIG .AND. IPTR .GT. 0 ) GOTO 20
C
         IF( IELM .NE. NELM ) THEN
           WRITE(0,*)'*** MKFIXD: Reached end of element list too early'
           WRITE(0,*) ielm,nelm,iptr,topbig
           ierr = -9595
           deallocate(srfcnt, srfavx, srfavy, srfavz)
           return
         end if
C
      END IF
c
      ewrite(2,*) 
     :     'Finished new element data'
      ewrite(2,*) 
     :     'Geom ed, int ed: ',cntgmy,cntint
      nusfid = count(srfcnt>0)
      if(nusfid.lt.200) then
         ewrite(2,*)
     :       'Surface averages  ( max id :',mxsfid,' )'
         do i = 0, MAX_SFID_STATS
            if( srfcnt(i) .gt. 0 ) then
               ewrite(2,*) 
     :              i,srfcnt(i),srfavx(i)/srfcnt(i)/3,
     :              srfavy(i)/srfcnt(i)/3,srfavz(i)/srfcnt(i)/3
            end if
         end do
      else
         ewrite(2,*) 
     :        'Check in MKFIXD:  NSELM, MXSFID, NUSFID =',
     :                                       nselm,mxsfid,nusfid
      end if
C
      IF( NPROC .GT. 1 ) THEN
C
         ewrite(2,*) 
     :        'Working out new gather array...'
         DO I = 1, NGATH
            INOD = GATHER(I)
            INOD = INT(NODLST(17,INOD))
            if(inod.gt.nnod) then
               WRITE(0,*) '  old:',gather(i),'  new:',inod,
     :                                   ' ***TOO BIG!!!',nnod
            else if(inod.lt.1) then
               WRITE(0,*) '  old:',gather(i),'  new:',inod,
     :                                   ' ***TOO SMALL!!!'
            else if( debug ) then
              WRITE(0,*) '  old:',gather(i),'  new:',
     :                                     inod,nodlst(6,inod)
            end if
            GATHER(I) = INOD
         END DO
C
         ewrite(2,*) 
     :        'Working out new scatter array...'
         DO I = 1, NHALO
            INOD = SCATER(I)
            INOD = INT(NODLST(17,INOD))
            if(inod.gt.nnod) then
               WRITE(0,*) '  old:',scater(i),'  new:',inod,
     :                                   ' ***TOO BIG!!!',nnod
            else if(inod.lt.1) then
               WRITE(0,*) '  old:',scater(i),'  new:',inod,
     :                                   ' ***TOO SMALL!!!'
            else if( debug ) then
               ewrite(2,*) 
     :              '  old:',scater(i),'  new:',
     :              inod,nodlst(6,inod)
            end if
            SCATER(I) = INOD
         END DO
C
      END IF
      ewrite(2,*) 
     :     'Checking local node ordering...'
C
      do i = 2, 8
         m(i) = 0.0
      end do
      m(1) = 1.0
      m(5) = 1.0
      m(9) = 1.0
c
      MAXASP = 0.0
      ASPAVE = 0.0
      VOLSUM = 0.0
      MAXVOL = 0.0
      MINVOL = 1E+20
      MINRAD = 1E+20
C - to avoid compiler warnings...
      MINELR = 0
      MAXELA = 0
      MAXELV = 0
      MINELV = 0
C
      DO I = 1, NELM
         IPOS = 0
c         print*,'element ',i,enlbas(i)+1
         DO J = ENLBAS(I)+1, ENLBAS(I+1)
            IPOS = IPOS + 1
            if(enlist(j).gt.nnod.or.enlist(j).le.0)
     :         WRITE(0,*) 'MKFIXD: bad node: ',enlist(j),i,j,ipos
            X(IPOS) = NODX(ENLIST(J))
            Y(IPOS) = NODY(ENLIST(J))
            Z(IPOS) = NODZ(ENLIST(J))
c            print*,'   node ',ipos,enlist(j),x(ipos),y(ipos),z(ipos)
         END DO
         CALL MTETIN( X, Y, Z, M, AVOL, AREAS, L, ARAD, ASP )
c         print*,'   vol,rad,asp ',avol,arad,asp
         ASPAVE = ASPAVE + ASP
         VOLSUM = VOLSUM + AVOL
         IF( AVOL .LE. 0.0 ) THEN
            WRITE(0,*) '+++ FOUND NEGATIVE VOLUME!',avol
            WRITE(0,*) I,(ENLIST(J),J=ENLBAS(I)+1,ENLBAS(I+1))
            J = ENLBAS(I) + 1
            IPOS = ENLIST(J)
            ENLIST(J)   = ENLIST(J+1)
            ENLIST(J+1) = IPOS
         END IF
         IF( ASP .GT. MAXASP ) THEN
            MAXASP = ASP
            MAXELA = I
            MAXAVL = AVOL
            MAXARD = ARAD
         END IF
         IF( AVOL .LT. MINVOL ) THEN
            MINVOL = AVOL
            MINELV = I
            MINVRD = ARAD
            MINVAS = ASP
         END IF
         IF( AVOL .GT. MAXVOL ) THEN
            MAXVOL = AVOL
            MAXELV = I
            MAXVRD = ARAD
            MAXVAS = ASP
         END IF
         IF( ARAD .LT. MINRAD ) THEN
            MINRAD = ARAD
            MINELR = I
            MINRVL = AVOL
            MINRAS = ASP
         END IF
      END DO
      ewrite(2,*) 
     :     'VOLSUM,ASPAVE: ',VOLSUM,ASPAVE/NELM
      ewrite(2,*) 
     :     'Maximum asp,vol,rad: ',maxasp,maxavl,maxard
      if( debug ) then
        ipos = enlbas(maxela)
        inod = enlist(ipos+1)
        print*,'  node 1 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+2)
        print*,'  node 2 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+3)
        print*,'  node 3 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+4)
        print*,'  node 4 : ',nodx(inod),nody(inod),nodz(inod)
      end if
      ewrite(2,*) 
     :     'Maximum vol,rad,asp: ',maxvol,maxvrd,maxvas
      if( debug ) then
        ipos = enlbas(maxelv)
        inod = enlist(ipos+1)
        print*,'  node 1 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+2)
        print*,'  node 2 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+3)
        print*,'  node 3 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+4)
        print*,'  node 4 : ',nodx(inod),nody(inod),nodz(inod)
      end if
      ewrite(2,*) 
     :     'Minimum vol,rad,asp: ',minvol,minvrd,minvas
      if( debug ) then
        ipos = enlbas(minelv)
        inod = enlist(ipos+1)
        print*,'  node 1 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+2)
        print*,'  node 2 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+3)
        print*,'  node 3 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+4)
        print*,'  node 4 : ',nodx(inod),nody(inod),nodz(inod)
      end if
      ewrite(2,*) 
     :     'Minimum rad,vol,asp: ',minrad,minrvl,minras
      if( debug ) then
        ipos = enlbas(minelr)
        inod = enlist(ipos+1)
        print*,'  node 1 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+2)
        print*,'  node 2 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+3)
        print*,'  node 3 : ',nodx(inod),nody(inod),nodz(inod)
        inod = enlist(ipos+4)
        print*,'  node 4 : ',nodx(inod),nody(inod),nodz(inod)
      end if
c
c                do i = 1, nelm
c                   iptr = 0
c                   do j = enlbas(i)+1, enlbas(i+1)
c                      iptr = iptr + 1
c                      xxx(iptr) = nodx(enlist(j))
c                      yyy(iptr) = nody(enlist(j))
c                      zzz(iptr) = nodz(enlist(j))
c                   end do
c                   aa = tetvol( xxx, yyy, zzz )
c                   if( aa .lt. 0.0 ) then
c                      j = enlbas(i)+1
c                      iptr = enlist(j)
c                      enlist(j) = enlist(j+1)
c                      enlist(j+1) = iptr
c                   end if
c                end do
c
c      print*,'ENLBAS: ',(enlbas(i),i=1,nelm+1)
c      print*,'ENLIST: ',(enlist(i),i=1,szenls)
c
      ewrite(2,*) 
     :     'Finished checking local node ordering'
C
  6   format( a, i7, a )
  7   format( a, 3(1pe13.5) )
  8   format( 2i7, 3(1pe13.5) )
  9   format( a, 2i7 )
C
      deallocate(srfcnt, srfavx, srfavy, srfavz)
      RETURN
      END
C
