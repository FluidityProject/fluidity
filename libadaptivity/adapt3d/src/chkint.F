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
      SUBROUTINE CHKINT( BIGLST, NODLST,
     :                   NODX, NODY, NODZ, OLDFLD, NNOD,
     :                   ENLBAS, ENLIST, SZENLS,
     :                   NELM, NEWFLD )
C-----------------------------------------------------------------------
C
C - This subroutine checks the interpolation of the first three fields
C - against the new x,y,z node co-ordinates - because (in fluidtp) the
C - first three fields should interpolate these coordinates and should do
C - so to within round-off errors.
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INTEGER NNOD, SZENLS, SZNELS, NELM,
     :        ENLBAS(NELM+1), ENLIST(SZENLS)
C
      REAL NODX(NNOD), NODY(NNOD), NODZ(NNOD), OLDFLD(NNOD,3)
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      REAL NEWFLD(NUMNDS,3)
C
      INCLUDE 'blkerr.i'
C
      integer i, j, k, l, iptr, inod
c
      real diff, maxdif, xx(5), yy(5), zz(5), b(5), tetvol
c
      real*8 sum1, sum2
c
      logical gtndfl
c
      sum1 = 0.0
      sum2 = 0.0
c
      maxdif = 0.0
      iptr = sttnod
c
 10   if( gtndfl(biglst,nodlst,iptr,1) ) then
         inod = int(nodlst(17,iptr))
         sum1 = sum1 + abs(newfld(inod,1))
         sum2 = sum2 + abs(nodlst(1,iptr))
         diff = newfld(inod,1) - nodlst(1,iptr)
         if( abs(diff) .gt. abs(maxdif) ) then
            maxdif = diff
            k = iptr
            j = inod
         end if
         iptr = int(nodlst(nxtnod,iptr))
         goto 10
      end if
c
      ewrite(3,*) 
     :     'Sum of real,intrp X: ',sum2,sum1
      ewrite(3,*) 
     :     'Largest X difference,iptr,inod: ',maxdif,k,j
      ewrite(3,*) 
     :     'Node co-ords: ',(nodlst(i,k),i=1,3)
      ewrite(3,*) 
     :     'Interpolated: ',(newfld(j,i),i=1,3)
      j = int(nodlst(16,k))
      l = enlbas(j)
      do i = 1, 4
         xx(i) = nodx(enlist(l+i))
         yy(i) = nody(enlist(l+i))
         zz(i) = nodz(enlist(l+i))
      end do
      ewrite(3,*) 
     :     'Inside old element,vol: ',j,tetvol(xx,yy,zz)
      call vals3d( xx, yy, zz, 5,
     :             nodlst(1,k), nodlst(2,k), nodlst(3,k), b )
      ewrite(3,*) 
     :     'Bases: ',(b(i),i=1,4)
      ewrite(3,*) 
     :     'Sum of bases: ',b(1)+b(2)+b(3)+b(4)
c
      sum1 = 0.0
      sum2 = 0.0
c
      maxdif = 0.0
      iptr = sttnod
c
 20   if( gtndfl(biglst,nodlst,iptr,1) ) then
         inod = int(nodlst(17,iptr))
         sum1 = sum1 + abs(newfld(inod,2))
         sum2 = sum2 + abs(nodlst(2,iptr))
         diff = newfld(inod,2) - nodlst(2,iptr)
         if( abs(diff) .gt. abs(maxdif) ) then
            maxdif = diff
            k = iptr
            j = inod
         end if
         iptr = int(nodlst(nxtnod,iptr))
         goto 20
      end if
c
      ewrite(3,*) 
     :     'Sum of real,intrp Y: ',sum2,sum1
      ewrite(3,*) 
     :     'Largest Y difference,iptr,inod: ',maxdif,k,j
      ewrite(3,*) 
     :     'Node co-ords: ',(nodlst(i,k),i=1,3)
      ewrite(3,*) 
     :     'Interpolated: ',(newfld(j,i),i=1,3)
      j = int(nodlst(16,k))
      l = enlbas(j)
      do i = 1, 4
         xx(i) = nodx(enlist(l+i))
         yy(i) = nody(enlist(l+i))
         zz(i) = nodz(enlist(l+i))
      end do
      call vals3d( xx, yy, zz, 5,
     :             nodlst(1,k), nodlst(2,k), nodlst(3,k), b )
      ewrite(3,*) 
     :     'Inside old element,vol: ',j,tetvol(xx,yy,zz)
      ewrite(3,*) 
     :     'Bases: ',(b(i),i=1,4)
      ewrite(3,*) 
     :     'Sum of bases: ',b(1)+b(2)+b(3)+b(4)
c
      sum1 = 0.0
      sum2 = 0.0
c
      maxdif = 0.0
      iptr = sttnod
c
 30   if( gtndfl(biglst,nodlst,iptr,1) ) then
         inod = int(nodlst(17,iptr))
         sum1 = sum1 + abs(newfld(inod,3))
         sum2 = sum2 + abs(nodlst(3,iptr))
         diff = newfld(inod,3) - nodlst(3,iptr)
         if( abs(diff) .gt. abs(maxdif) ) then
            maxdif = diff
            k = iptr
            j = inod
         end if
         iptr = int(nodlst(nxtnod,iptr))
         goto 30
      end if
c
      ewrite(3,*) 
     :     'Sum of real,intrp Z: ',sum2,sum1
      ewrite(3,*) 
     :     'Largest Z difference,iptr,inod: ',maxdif,k,j
      ewrite(3,*) 
     :     'Node co-ords: ',(nodlst(i,k),i=1,3)
      ewrite(3,*) 
     :     'Interpolated: ',(newfld(j,i),i=1,3)
      j = int(nodlst(16,k))
      l = enlbas(j)
      do i = 1, 4
         xx(i) = nodx(enlist(l+i))
         yy(i) = nody(enlist(l+i))
         zz(i) = nodz(enlist(l+i))
      end do
      call vals3d( xx, yy, zz, 5,
     :             nodlst(1,k), nodlst(2,k), nodlst(3,k), b )
      ewrite(3,*) 
     :     'Inside old element,vol: ',j,tetvol(xx,yy,zz)
      ewrite(3,*) 
     :     'Bases: ',(b(i),i=1,4)
      ewrite(3,*) 
     :     'Sum of bases: ',b(1)+b(2)+b(3)+b(4)
C
      RETURN
      END
C
