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
      subroutine assval( x, y, z, nnod, enlb, enls, nelm, szenls,
     :                   ecls, nelb, nels, sznels, tmpfld,
     :                   elmreg, sigmat, nmat, field, bndlyr )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
c
      integer nnod, nelm, szenls, enlb(nelm+1), enls(szenls),
     :        ecls(4*nelm), elmreg(nelm), nmat, sznels, nelb(nnod+1),
     :        nels(sznels)
c
      real x(nnod), y(nnod), z(nnod), field(nnod), sigmat(nmat),
     :     tmpfld(nelm)
c
      logical bndlyr
c
      integer ielm, nln, ibas, iln, i1, i2, i3, i4, inod, conelm, i,
     :        conbas, connln, j1, j2, j3, its, added, jbas, ipos,
     :        setit, chked
c
      real m(9), vol, areas(4), lngths(6), rad, q, xx(4), yy(4), zz(4),
     :     qq(4), d
c
      logical missed, setval, debug
c
      debug = .false.
c
      m(1) = 1.0
      do i = 2, 8
         m(i) = 0.0
      end do
      m(9) = 1.0
      m(4) = 1.0
c
      do inod = 1, nnod
         field(inod) = 0.0
      end do
c
      do ielm = 1, nelm
c
         ibas = enlb(ielm)
         nln  = enlb(ielm+1) - ibas
c
         i1 = enls(ibas+1)
         i2 = enls(ibas+2)
         i3 = enls(ibas+3)
         i4 = enls(ibas+4)
c
         j1 = 1
         j2 = 2
         j3 = 0
c
         do i = 1, 4
c
            if( j3 .eq. 0 ) then
               j3 = 3
            else if( j3 .eq. 3 ) then
               j3 = 4
            else if( j2 .eq. 2 ) then
               j2 = 3
            else if( j1 .eq. 1 ) then
               j1 = 2
            end if
c
            conelm = ecls(4*ielm-4+i)
c
            if( conelm .gt. 0 ) then
c
               conbas = enlb(conelm)
               connln = enlb(conelm+1) - conbas
c
               if( elmreg(ielm) .ne. elmreg(conelm) ) then
                  field(enls(ibas+j1)) = 1.0
                  field(enls(ibas+j2)) = 1.0
                  field(enls(ibas+j3)) = 1.0
               end if
c
            else if( bndlyr ) then
c
               field(enls(ibas+j1)) = 1.0
               field(enls(ibas+j2)) = 1.0
               field(enls(ibas+j3)) = 1.0
c
            end if
c
         end do
c
      end do
c
      its = 0
c      return
c
 10   its = its + 1
      added = 0
      setit = 0
      chked = 0
c
      do 1000 inod = 1, nnod
c
c         if( inod .eq. 1 ) debug = .true.
         if( field(inod) .gt. 0.0 ) goto 1000
c         field(inod) = 1.0
c         goto 1000
c
         chked = chked + 1
c
         setval = .false.
c
         xx(4) = x(inod)
         yy(4) = y(inod)
         zz(4) = z(inod)
         if( debug ) print*,'   node,pos: ',inod,xx(4),yy(4),zz(4)
c
         missed = .false.
c
         nln = nelb(inod+1) - nelb(inod)
c
         do ipos = nelb(inod)+1, nelb(inod+1)
c
            ielm = nels(ipos)
c
            jbas = enlb(ielm)
c            nln  = enlb(ielm+1) - jbas
c
            i1 = enls(jbas+1)
            i2 = enls(jbas+2)
            i3 = enls(jbas+3)
            i4 = enls(jbas+4)
c
            j1 = 1
            j2 = 0
c
            if( i1 .ne. inod ) then
               xx(j1) = x(i1) - xx(4)
               yy(j1) = y(i1) - yy(4)
               zz(j1) = z(i1) - zz(4)
               qq(j1) = abs(field(i1))
               if( qq(j1) .ne. 0.0 ) j2 = j2 + 1
               j1 = j1 + 1
            end if
c
            if( i2 .ne. inod ) then
               xx(j1) = x(i2) - xx(4)
               yy(j1) = y(i2) - yy(4)
               zz(j1) = z(i2) - zz(4)
               qq(j1) = abs(field(i2))
               if( qq(j1) .ne. 0.0 ) j2 = j2 + 1
               j1 = j1 + 1
            end if
c
            if( i3 .ne. inod ) then
               xx(j1) = x(i3) - xx(4)
               yy(j1) = y(i3) - yy(4)
               zz(j1) = z(i3) - zz(4)
               qq(j1) = abs(field(i3))
               if( qq(j1) .ne. 0.0 ) j2 = j2 + 1
               j1 = j1 + 1
            end if
c
            if( i4 .ne. inod ) then
               xx(j1) = x(i4) - xx(4)
               yy(j1) = y(i4) - yy(4)
               zz(j1) = z(i4) - zz(4)
               qq(j1) = abs(field(i4))
               if( qq(j1) .ne. 0.0 ) j2 = j2 + 1
               j1 = j1 + 1
            end if
c
            if( j2 .gt. 0 ) then
c
              if( debug ) print*,'      elm: ',ielm,sigmat(elmreg(ielm))
              do i = 1, 3
                d = xx(i)*xx(i) + yy(i)*yy(i) + zz(i)*zz(i)
                d = sqrt(d)*sigmat(elmreg(ielm))
                if( debug ) then
                   print*,'        nod,val,d: ',i,j2,qq(i),d
c                   print*,'            pos: ',xx(i),yy(i),zz(i)
                end if
c                print*,'   d is ',d,inod
                if( d .lt. 50.0 )
     :                field(inod) = field(inod) - qq(i)*exp(-d)/j2
              end do
c
              if( j2 .ne. 3 ) missed = .true.
              setval = .true.
c
            end if
c
         end do
c
         if( setval ) then
c
            field(inod) = field(inod)/nln
            if( field(inod) .eq. 0.0 ) field(inod) = -1e-20
c
            if( .not. missed ) then
               field(inod) = abs(field(inod))
               setit = setit + 1
            end if
c
            added = added + 1
            if( debug ) print*,'     final val: ',field(inod)
c
         end if
c
 1000 continue
c
      print*,'It: ',its,'  added: ',added,'  set: ',
     :       setit,'  nonz: ',chked
c
      if( its .lt. 1000 .and. added .gt. 0 ) goto 10
c
      return
      end
c
