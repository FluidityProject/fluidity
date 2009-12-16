!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"
subroutine fprint_backtrace
end subroutine fprint_backtrace

program evmerge
! - This program merges two or more event input files and meshes and is used for FETCH

! ******************************
! * DOES NOT USE IMPLICIT NONE *
! ******************************

! This needs a spring clean, some variables are NOT declared and there are lots of continue statements 
! and GOTO :(

! BUG
! - doesnt seem to get the merged matxs file correct
! - quick fix is to just copy .1.ev matxs to .12.ev as .1.ev and .2.ev are
!   probably the same anyway (no merge)

  integer aux, infil, inmsh, inmat
  integer outfil, outmsh, outmat
  INTEGER DATFMT
  !
  parameter ( maxnod = 500000, maxelm = 500000, &
       &            maxcon = 500000, maxfil = 1000,&
       & infil  = 10, inmsh = 11, inmat  = 12,&
       & outfil = 20, outmsh = 21, outmat = 22, &
       & aux    = 31, eps = 1.0e-05 )
  !
  real x(maxnod), y(maxnod), z(maxnod), mixden(maxelm)
  !
  real plaval(4,10000)
  !
  integer nodes(10,maxelm), newnod(maxnod), newroo(maxelm),&
       &        elmtyp(maxelm), elmreg(maxelm), elmmat(maxelm), &
       &        elmsrc(maxelm), elmroo(maxelm), nln(maxelm),&
       &        connod(maxcon), contyp(maxcon),&
       &        oldelm(maxelm), newcon(maxnod), oldnod(maxnod),&
       &        mixnum(maxelm), mixcom(maxelm), nodend(maxfil)
  !
  integer case, iprob, negom, nadj, soltyp, restrt,&
       &              flxopt, srcopt, tmpopt, emat, tmp,&
       &              nselm, nvelm, nlist, npntr, nelst
  !
  integer nsctr, lsctr, ksctr, angtyp
  !
  logical raytrc, mixing, upscat, fiss, totmerg, samemat, readmix,&
       &        newmat, newmix, binmsh, binmat, dfiss, matelm, coupled
  !
  character margin*4, filnam*80, line*80, rest(1000)*80, trgtnam*80
  !
  character path*80, str1*20, str2*80, str3*20, str4*20, marg*10
  !
  character str5*10, str6*10, str7*10, card*80
  !
  character comlin(10000)*80, matlin(1000000)*80, srclin(10000)*80,&
       &          bmlin(10000)*80, vellin(10000)*80, alblin(10000)*80,&
       &          spclin(10000)*80, bgrlin(10000)*80, bealin(10000)*80,&
       &          detlin(10000)*80, plalin(10000)*80, rsplin(10000)*80,&
       &          rstlin(10000)*80, dellin(10000)*80, sorlin(10000)*80
  !
  character*11 mshfmt, matfmt
  !
  character*1  quote
  !
  data      quote / '''' /
  !
  write( *, * )
  write( *, * )'**************************************************'
  write( *, * )'*** Program evmerge - version 3.1a - Aug 2001  ***'
  write( *, * )'**************************************************'
  write( *, * )
  !
  call getenv( "EVENT_DATA", path )
  !
  kp = index( path, ' ' ) - 1
  !
  ! *************************************************************************************
  !
  ! - read target file name
  !
  ! *************************************************************************************
  !
  print*
  !
  write( *, '( a, $ )' )&
       &   'Target file name                              => '
  !
  read( *, '( a )' ) trgtnam
  !
  kt = index( trgtnam, ' ' ) - 1
  !
  print*
  !
  write( *, '( a, $ )' )&
       &   'Total material merge (0) or same material (1) => '
  !
  read( *, * ) iopt
  !
  samemat = iopt .eq. 1
  !
  totmerg = iopt .ne. 1
  !
  write( *, '( a, $ )' )&
       &   'Region material (0) or element material (1)   => '
  !
  read( *, * ) iopt2
  !
  matelm = iopt2 .eq. 1
  !
  open( outfil, file = path(1:kp)//trgtnam(1:kt), &
       &      status = 'unknown', form = 'formatted' )
  !
  open( aux, status = 'scratch', form = 'unformatted' )
  !
  kcr    = 0
  kmix   = 0
  kcom   = 0
  kelm   = 0
  knod   = 0
  kcon   = 0
  kroot  = 0
  kreg   = 0
  kmat   = 0
  ksor   = 0
  ifil   = 0
  emat   = 0
  kiso   = 0
  ksurfr = 0
  kplan  = 0
  ksrc   = 0
  kdet   = 0
  kbeam  = 0
  !
  klcom  = 0
  klmat  = 0
  klsrc  = 0
  klsor  = 0
  klalb  = 0
  klbgr  = 0
  klbea  = 0
  kldet  = 0
  klpla  = 0
  klrsp  = 0
  klrst  = 0
  kldel  = 0
  klspc  = 0

  KBEAMS = 0
  KDEN = 0
  KSRCS = 0
  KRESP = 0
  KDETEC = 0
  KALB = 0
  !
  maxmat = 0
  !
  irest  = 0
  !
  raytrc = .false.
  !
  mixing = .false.
  !
  fiss   = .false.
  dfiss  = .false.
  !
  newmat = .true.
  newmix = .true.
  !
  binmat = .false.
  binmsh = .false.
  !
  xmin   = +1.0e+30
  xmax   = -1.0e+30
  !
  ymin   = +1.0e+30
  ymax   = -1.0e+30
  !
  zmin   = +1.0e+30
  zmax   = -1.0e+30
  !
  ! *************************************************************************************
  !
  ! - beginning of mesh loop
  !
  ! *************************************************************************************
  !
1 print*
  !
  write( *, '( a, $ )' )&
       &   'File name (type 0 to finish)                  => '
  !
  read( *, '( a )' ) filnam
  !
  if = index( filnam, ' ' ) - 1
  !
  if( filnam(1:1) .eq. '0' .or. filnam(1:1) .eq. '0' ) go to 50
  !
  open( infil, file = path(1:kp)//filnam(1:if), &
       &     status = 'unknown', form = 'formatted' )
  !
  ifil = ifil + 1
  print*,'ifil:',ifil,infil,outfil,path(1:kp)//filnam(1:if)
  !
  read( infil, '(a)' ) line
  if( ifil .eq. 1 ) write( outfil, '(a)' ) line
  !
  read( infil, '(a)' ) line
  if( ifil .eq. 1 ) write( outfil, '(a)' ) line
  !
  read( infil, '(a, 8i5 )' ) marg, iprob, ngeom, nadj, soltyp,&
       &                           restrt, flxopt, srcopt, tmpopt
  !
  case = mod( iprob, 10 )
  !
  myprob = iprob - case
  !
  if( ifil .eq. 1 ) write( outfil, '(a, 8i5)' ) marg, iprob, ngeom,&
       &                     nadj, soltyp, restrt,&
       &                     flxopt, srcopt, tmpopt
  !
  read( infil, '(a)' ) line
  !
  if( ifil .eq. 1 ) write( outfil, '(a)' ) line
  !
  if( case .eq. 3 ) then
     !
     read( infil, '(a)' ) line
     if( ifil .eq. 1 ) write( outfil, '(a)' ) line
     !
  end if
  !
  read( infil, '( a, 4i5 )' ) marg, nsctr, lsctr, ksctr, angtyp
  !
  if( ifil .eq. 1 ) then
     !
     nsctr1 = nsctr + 1
     lsctr1 = lsctr + 1
     ksctr1 = ksctr + 1
     !
     write( outfil, '( a, 4i5 )' ) marg, nsctr, lsctr, ksctr,&
          &                                 angtyp
     !
  end if
  !
  if( nsctr1 .ne. nsctr+1 .or. ksctr1 .ne. ksctr+1 .or.&
       &    ksctr1 .ne. ksctr+1 ) then
     !
     print*
     print*, 'Angular parameters disagree !!!!'
     print*
     !
  end if
  ! 
  read( infil, '( a, 14i5 )' ) marg, ngrps, nbgrps, lump,&
       &                             nupsct, ndgrps, nngrps, nggrps
  !
  upscat = nupsct .ne. 0
  !
  if( ifil .eq. 1 ) then
     !
     kgrps  = ngrps
     kbgrps = nbgrps
     klump  = lump
     kupsct = nupsct
     kngrps = nngrps
     kggrps = nggrps
     !
     write( outfil, '( a, 14i5 )' ) marg, ngrps, nbgrps, lump,&
          &                              nupsct, ndgrps, nngrps, nggrps
     !
  end if
  !
  if( ngrps .ne. kgrps .or. nbgrps .ne. kbgrps .or.&
       &    klump .ne. lump  .or. kupsct .ne. nupsct .or.&
       &    kngrps .ne. nngrps .or. kggrps .ne. nggrps ) then
     !
     print*
     print*, 'Energy parameters disagree !!!!'
     print*
     !
     print*,  kgrps, ngrps
     print*,  kbgrps, nbgrps
     print*,  klump, lump
     print*,  kupsct, nupsct
     print*,  kngrps, nngrps
     print*,  kggrps, nggrps
     !
  end if
  !
  ! *************************************************************************************
  !
  ! - read mesh file name and open it
  !
  ! *************************************************************************************
  !
  read( infil, * ) str1, str2, mshfmt
  !
  print*, str1, str2(1:25), ' ', mshfmt
  !
  if( mshfmt(1:1) .eq. 'f' ) then
     !
     open( inmsh, file = path(1:kp)//'mesh/'//str2, &
          &         status = 'unknown', form = 'formatted' )
     !
  else
     !
     binmsh = .true.
     open( inmsh, file = path(1:kp)//'mesh/'//str2, &
          &         status = 'unknown', form = 'unformatted' )
     !
  end if
  !
  ! *************************************************************************************
  !
  ! - read and open material file name
  !
  ! *************************************************************************************
  !
  read( infil, * ) str1, str2, str3, matfmt
  !
  print*, str1, str2(1:25), ' ', matfmt
  !
  if( matfmt(1:1) .eq. 'f' ) then
     !
     open( inmat, file = path(1:kp)//'matxs/'//str2, &
          &        status = 'unknown', form = 'formatted' )
     !
  else
     !
     binmat = .true.
     open( inmat, file = path(1:kp)//'matxs/'//str2, &
          &        status = 'unknown', form = 'unformatted' )
     !
  end if
  !
  ! *************************************************************************************
  !
  ! - read material card
  !
  ! *************************************************************************************
  !
  read( infil, * ) marg, mcr, mtps, nmix,&
       &                 iht, ihs, ihm, datfmt, nden, nbeams,&
       &                 ndetec, nsrcs, nresp
  !
  print*,'====>',datfmt
  nds = ihm - ihs
  nds = min( nds, ngrps - 1 )
  jmix = nmix
  !
  if( samemat .and. ifil .ne. 1 ) then
     !
     mcr    = 0
     mtps   = 0
     nmix   = 0
     nden   = 0
     nbeams = 0
     nsrcs  = 0
     nresp  = 0
     !
     ndgrps = 0
     !
     newmix = .false.
     newmat = .false.
     !
  end if
  !
  if( matelm ) then
     !
     newmix = .true.
     !
     nmix = jmix
     !
  end if
  !
  ! *************************************************************************************
  !
  ! - read comments
  !
  ! *************************************************************************************
  !
  read( infil, '(a, i5)' ) marg, ncom
  !
  klcom = klcom + 1
  !
  comlin(klcom) = 'File '//filnam(1:if)
  !
  do icom = 1, ncom
     !
     klcom = klcom + 1
     !
     read( infil, '( a )' ) comlin(klcom)
     !
  end do
  !
  kcom = kcom + ncom + 1
  !
  print*, 'ncom   :', ncom, ' kcom   :', kcom
  !
  ! *************************************************************************************
  !
  ! - read mixing instructions
  !
  ! *************************************************************************************
  !
  if( nmix .gt. 0 .and. newmix ) then
     !
     mixing = .true.
     !
     maxnum = 0
     !
     do k = kmix+1, kmix+nmix
        !
        read( infil, '( a, 2i5, e15.7 )' ) marg, mixnum(k), &
             &                                         mixcom(k), mixden(k)
        !
        maxnum    = max( maxnum, mixnum(k) )
        !
        if( matelm ) then
           !
           mixnum(k) = mixnum(k) + mixnum(kmix)
           !
        else
           !
           mixnum(k) = mixnum(k) + kmat
           !
        end if
        !
     end do
     !
     print*, 'maxnum :', maxnum
     !
     niso = 0
     !
     do k = kmix+1, kmix+nmix
        !
        if( mixcom(k) .gt. 0 ) then
           !
           mixcom(k) = mixcom(k) - maxnum + kiso
           !
           if( samemat ) mixcom(k) = mixcom(k) - kiso
           !
        end if
        !
        niso = max( niso, mixcom(k) )
        !
     end do
     !
     kiso = niso
     !
  end if
  !
  ! - increment kmix
  !
  kmix = kmix + nmix
  !
  print*, 'nmix   :', nmix, ' kmix   :', kmix
  !
  newmix = totmerg
  !
  ! *************************************************************************************
  !
  ! - read mesh file - here!
  !
  ! *************************************************************************************
  !
  if( mshfmt(1:1) .eq. 'f' ) then
     !
     read( inmsh, '( 10x, 13i5, l5 )' ) nreg, nmat, nsor,nelm, nnod,&
          &                                  ncon, nalb, nopt, nplan,nprocs,&
          &                                  nprocno, nsurfr,tmp,coupled
     !
     if( coupled ) &
          !
          &      read( inmsh, '( 10x, 5i10 )' ) nselm, nvelm,&
          &                                     nlist, npntr, nelst
     !
  else
     !
     read( inmsh ) nreg,nmat, nsor,nelm, nnod,&
          &                 ncon, nalb, nopt, nplan, &
          &                 nprocs,nprocno,&
          &                 nsurfr,tmp,coupled
     !
     if( coupled )&
          !
          &     read( inmsh ) nselm, nvelm, nlist, npntr, nelst
     !
  end if
  !
  print*, 'nelm   : ', nelm, ' nnod   : ', nnod, ' ncon   : ', ncon
  !
  print*, 'nreg   : ', nreg, ' nmat   : ', nmat, ' nsor   : ', nsor
  !
  raytrc = nsurfr .ne. 0 .or. raytrc
  !
  if( raytrc ) write( *, * ) 'Ray-tracing is on ...'
  !
  ksurfr = ksurfr + nsurfr
  !
  print*, 'nsurfr :', nsurfr, ' ksurfr :', ksurfr
  !
  if( samemat .and. ifil .ne. 1 ) then
     !
     nsor = 0
     nalb = 0
     !
  end if
  !
  ! *************************************************************************************
  !
  ! - elements
  !
  ! *************************************************************************************
  !
  maxreg = 0
  maxroo = 0
  !
  imat   = 0
  isor   = 0
  !
  jmat   = emat
  !
  do  ielm = 1, nelm! Was loop 10
     !
     kelm = kelm + 1
     !
     if( mshfmt(1:1) .eq. 'f' ) then
        !
        read( inmsh, '( 5x, 15i5 )' ) iroot, elmtyp(kelm),&
             &                              elmreg(kelm), &
             &                              elmmat(kelm), elmsrc(kelm), &
             &                              ( nodes(i,kelm), i = 1, 10 )
        !
        do i = 1, 10
           !
           if( nodes(i,kelm) .ne. 0 ) then
              !
              nodes(i,kelm) = nodes(i,kelm) + knod
              nln(kelm) = i
              !
           end if
           !
        end do
        !
     else
        !
        read( inmsh ) iroot, elmtyp(kelm),&
             &                              elmreg(kelm), &
             &                     elmmat(kelm), elmsrc(kelm), nln(kelm),&
             &                      ( nodes(i,kelm), i = 1, nln(kelm) )
        !
        do i = 1, nln(kelm)
           !
           nodes(i,kelm) = nodes(i,kelm) + knod
           !
        end do
        !
     end if
     !
     elmroo(kelm) = iroot + kroot
     !
     ! - take this out so that averaging regions are not merged
     !         elmreg(kelm) = elmreg(kelm) + kreg
     !
     if( mod( elmtyp(kelm), 2 ) .eq. 1 ) then
        !
        emat = max( elmmat(kelm), emat )
        imat = max( elmmat(kelm), imat )
        isor = max( elmsrc(kelm), isor )
        !
        if( totmerg ) elmmat(kelm) = elmmat(kelm) + kmat
        !
        if( elmsrc(kelm) .ne. 0 )&
             &          elmsrc(kelm) = elmsrc(kelm) + ksor
        !
     end if
     !
     maxreg = max( elmreg(kelm), maxreg )
     maxroo = max( iroot, maxroo )
     !
  end do ! Was loop 10
  !
  print*, 'imat   :', imat, ' emat   :', emat
  !
  !      ksor  = ksor  + isor
  if( totmerg ) then
     kmat = kmat + imat
  else if( kmat .eq. 0 ) then
     kmat = kmat + nmat
  end if
  !
  !      if( totmerg .or. kmat .eq. 0 ) kmat  = kmat  + nmat
  !
  print*, 'nmat   :', nmat, ' kmat   :', kmat
  !
  kroot = kroot + maxroo
  kreg  = kreg  + maxreg
  !
  ! *************************************************************************************
  !
  ! - nodes
  !
  ! *************************************************************************************
  !
  do inod = 1, nnod
     !
     if( mshfmt(1:1) .eq. 'f' ) then
        !
        read( inmsh, '( 10x, 3e15.7 )' ) x(inod+knod), y(inod+knod), &
             &                                       z(inod+knod)
        !
     else
        !
        read( inmsh ) x(inod+knod), y(inod+knod), &
             &                                  z(inod+knod)
        !
     end if
     !
     xmin = min( xmin, x(inod+knod) )
     xmax = max( xmax, x(inod+knod) )
     !
     ymin = min( ymin, y(inod+knod) )
     ymax = max( ymax, y(inod+knod) )
     !
     zmin = min( zmin, z(inod+knod) )
     zmax = max( zmax, z(inod+knod) )
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - constraints
  !
  ! *************************************************************************************
  !
  do j = 1, ncon
     !
     kcon = kcon + 1
     !
     if( mshfmt(1:1) .eq. 'f' ) then
        !
        read( inmsh, '( 10x, 3i5 )' ) connod(kcon), contyp(kcon)
        !            print*,'cons ',kcon,connod(kcon),contyp(kcon),knod
        !
     else
        !
        read( inmsh ) connod(kcon), contyp(kcon)
        !            if( connod(kcon) .eq. 41638 ) then
        !              print*,'cons ',kcon,connod(kcon),contyp(kcon),knod
        !            end if
        !
     end if
     !
     connod(kcon) = connod(kcon) + knod
     if( contyp(kcon) .lt. 0 ) contyp(kcon) = contyp(kcon) - knod
     !
  end do
  !
  knod = knod + nnod
  !
  nodend(ifil) = knod
  !
  ! *************************************************************************************
  !
  ! - read materials
  !
  ! *************************************************************************************
  !
  if( totmerg .or. klmat .lt. 1 ) then
     !
     totfis = 0.0
     !
     ! - formatted input
     !
     do i = 1, mcr
        !
        init = 1
        last = ngrps
        !
        if( nds .eq. 1 ) last = 1
        !
        do igrp = 1, ngrps
           !
           do k = 1, iht, 4
              !
              klmat = klmat + 1
              !
              read( inmat, '( a )', end = 4 ) matlin(klmat)
              !
              if( k .eq. 1 ) then
                 read( matlin(klmat), '(10X,3E15.7)' ) a, b, c
                 totfis = totfis + c
              end if
              !
           end do
           !
           do l = 0, lsctr
              !
              do k = init, last, 4
                 !
                 k3 = min( k+3, last )
                 !
                 klmat = klmat + 1
                 !
                 read( inmat, '( a )', end = 4 ) matlin(klmat)
                 !
              end do
              !
           end do
           !
           init = last + 1
           last = last + ngrps
           !
           if( .not. upscat ) last = last - igrp
           if(   nds .eq. 0 ) last = init
           !
        end do
        !
     end do
     !
     ! - binary format
     !
     do i = 1, mtps
        !
        init = 1
        last = ngrps
        !
        if( nds .eq. 1 ) last = 1
        !
        do igrp = 1, ngrps
           !
           do k = 1, iht, 4
              !
              klmat = klmat + 1
              !
              read( inmat, end = 4 ) matlin(klmat)
              !
           end do
           !
           do l = 0, nsctr
              !
              do k = init, last, 4
                 !
                 k3 = min( k+3, last )
                 !
                 klmat = klmat + 1
                 !
                 read( inmat, end = 4 ) matlin(klmat)
                 !
              end do
              !
           end do
           !
           init = last + 1
           last = last + ngrps
           !
           if( .not. upscat ) last = last - igrp
           if(   nds .eq. 0 ) last = init
           !
        end do
        !
     end do
     !
4    continue
     !
     kcr  = kcr  + mcr
     !
     ktps = ktsp + mtps
     !
     print*, 'mcr    :', mcr, ' kcr    :', kcr
     !
     print*, 'mtps   :', mcr, ' ktps   :', kcr
     !
     print*, 'read', klmat, ' material cards'
     !
  end if
  !
  ! ******************************************************************************
  !
  ! -  read velocities
  !
  ! ******************************************************************************
  !
  print*,'case, newmat :',case,newmat

  if( case .eq. 3 .and. newmat ) then
     !
     if( matfmt(1:1) .eq. 'f' ) then
        !
        do ig = 1, ngrps, 4
           !
           if( ifil .eq. 1 ) then
              !
              klvel = klvel + 1
              !
              read( inmat, '( a )' ) vellin(klvel)
              !
           else
              !
              read( inmat, '( a )' )
              !
           end if
           !
        end do
        !
     else
        !
        !
     end if
     !
  end if
  !
  print*, 'read', klvel, ' velocity cards'
  !
  ! *************************************************************************************
  !
  ! - read sources
  !
  ! *************************************************************************************
  !
  do isor = 1, nsor
     !
     do igrp = 1, ngrps, 4
        !
        klsor = klsor + 1
        !
        if( matfmt(1:1) .eq. 'f' ) then
           !
           read( inmat, '( a )', end = 5 ) sorlin(klsor)
           !
        else
           !
           read( inmat, end = 5 ) sorlin(klsor)
           !
        end if
        !
     end do
     !
  end do
  !
5 continue
  !
  ksor = ksor + nsor
  !
  print*, 'nsor   :', nsor, ' ksor   :', ksor
  !
  print*, 'read', klsor, ' source cards'
  !
  ! *************************************************************************************
  !
  ! - read albedos
  !
  ! *************************************************************************************
  !
  do ialb = 1, nalb
     !
     do igrp = 1, ngrps, 4
        !
        klalb = klalb + 1
        !
        if( matfmt(1:1) .eq. 'f' ) then
           !
           read( inmat, '( a )', end = 6 ) alblin(klalb)
           !
        else
           !
           read( inmat, end = 5 ) alblin(klalb)
           !
        end if
        !
     end do
     !
  end do
  !
6 continue
  !
  kalb = kalb + nalb
  !
  print*, 'nalb   :', nalb, ' kalb   :', kalb
  !
  print*, 'read', klalb, ' albedo cards'
  !
  ! ******************************************************************************
  !
  ! -  read fission spectrum
  !
  ! ******************************************************************************
  !
  if( ( totfis.gt.0.0 .or. case .eq. 1 ) .and. newmat ) then
     !
     if( matfmt(1:1) .eq. 'f' ) then
        !
        do ig = 1, ngrps, 4
           !
           if( .not. fiss ) then
              !
              klspc = klspc + 1
              !
              read( inmat, '( a )', end = 7 ) spclin(klspc)
              !
           else
              !
              read( inmat, '( a )' )
              !
           end if
           !
        end do
        !
     else
        !
        !
     end if
     !
     if( .not. fiss ) then
        !
        print*, 'read fission spectrum'
        !
     else
        !
        print*, '*ignored* fission spectrum'
        !
     end if
     !
     fiss = .true.
     !
     print*, 'read', klspc, ' fission spectrum cards'
     !
  end if
  !
  ! ******************************************************************************
  !
  ! -  read delayed neutron parameters
  !
  ! ******************************************************************************
  !
  if( ndgrps .gt. 0 ) then
     !
     if( matfmt(1:1) .eq. 'f' ) then
        !
        do i = 1, ndgrps, 4
           !
           do ig = 1, ngrps+2
              !
              if( .not. dfiss ) then
                 !
                 kldel = kldel + 1
                 read( inmat, '( a )', end = 7 ) dellin(kldel)
                 !
              else
                 !
                 read( inmat, '( a )' )
                 !
              end if
              !
           end do
           !
        end do
        !
     else
        !
        !
     end if
     !
     if( .not. dfiss ) then
        !
        print*, 'read delayed groups: ',ndgrps
        !
     else
        !
        print*, '*ignored* delayed groups: ',ndgrps
        !
     end if
     !
     dfiss = .true.
     !
     print*, 'read', kldel, ' delayed parameters cards'
     !
  end if
  !
7 continue
  !
  ! ******************************************************************************
  !
  ! - read broad group structure
  !
  ! ******************************************************************************
  !
  do i = 1, nbgrps, 4
     !
     if( ifil .eq. 1 ) then
        !
        klbgr = klbgr + 1
        !
        read( infil, '( a )', end = 9 ) bgrlin(klbgr)
        !
     else
        !
        read( infil, '( a )', end = 9 )
        !
     end if
     !
  end do
  !
  print*, 'read', klbgr, ' broad groups cards'
  !
9 continue
  !
  ! *************************************************************************************
  !
  ! - read general source
  !
  ! *************************************************************************************
  !
  do isrc = 1, nsrcs
     !
     if( matfmt(1:1) .eq. 'f' ) then
        !
        read( inmat, '( a, 2i5, 3( i5, 5x, a ), i5, e15.7 )' )&
             &            marg, ityp1, ityp2, isor1, str5, isor2, str6,&
             &            isor3, str7, itv, time
        !
        klsrc = klsrc + 1
        !
        if( isor1 .gt. 0 ) isor1 = isor1 + ksor - nsor
        if( isor2 .gt. 0 ) isor2 = isor2 + ksor - nsor
        if( isor3 .gt. 0 ) isor3 = isor3 + ksor - nsor
        !
        write( srclin(klsrc),&
             &          '( a, 2i5, 3( i5, 5x, a ), i5, e15.7 )' ) marg,&
             &           ityp1, ityp2, isor1, str5, isor2, str6, isor3, str7,&
             &           itv, time
        !
        do i = 1, 4
           !
           klsrc = klsrc + 1
           !
           read( inmat, '( a )' ) srclin(klsrc)
           !
        end do
        !
     else
        !
        read( inmat )&
             &           ityp1, ityp2, isor, str5, isor2, str6, isor3, str7,&
             &           itv, time
        !
        klsrc = klsrc + 1
        !
        if( isor1 .gt. 0 ) isor1 = isor1 + ksor - nsor
        if( isor2 .gt. 0 ) isor2 = isor2 + ksor - nsor
        if( isor3 .gt. 0 ) isor3 = isor3 + ksor - nsor
        !
        write( srclin(klsrc),&
             &          '( 2i5, 3( i5, 5x, a ), i5, e15.7 )' )&
             &           ityp1, ityp2, isor1, str5, isor2, str6, isor3, str7,&
             &           itv, time
        !
        do i = 1, 4
           !
           klsrc = klsrc + 1
           !
           read( inmat ) srclin(klsrc)
           !
        end do
        !
     end if
     !
  end do
  !
  ksrcs = ksrcs + nsrcs
  !
  print*, 'nsrcs  :', nsrcs, ' ksrcs  :', ksrcs
  !
  print*, 'read', klsrc, ' general source cards'
  !
  ! *************************************************************************************
  !
  ! - read beams
  !
  ! *************************************************************************************
  !
  do ibeam = 1, nbeams
     !
     if( mshfmt(1:1) .eq. 'f' ) then
        !
        read( inmat,'( a, 2i5, a, e15.7, i5, e15.7, 2i5 )' )&
             &           ityp, isor, str4, radius, itv, time, isurf
        !
        klbea = klbea + 1
        !
        write( bealin(klbea),&
             &            '( a, 2i5, a, e15.7, i5, e15.7, 2i5 )' )&
             &           ityp, isor+ksor, str4, radius, itv, time, isurf
        !
     else
        !
        read( inmat ) ityp, isor, str5, radius, itv, time, isurf
        !
        klbea = klbea + 1
        !
        write( bealin(klbea), '( 2i5, a, e15.7, i5, e15.7, 2i5 )' )&
             &           ityp, isor+ksor, str5, radius, itv, time, isurf
        !
     end if
     !
     klbea = klbea + 1
     !
     read( inmat, '( a )' ) bealin(klbea)
     !
     klbea = klbea + 1
     !
     read( inmat, '( a )' ) bealin(klbea)
     !
     isurf = mod( isurf, 1000 )
     !
     if( mshfmt(1:1) .eq. 'f' ) then
        !
        do i = 1, isurf, 14
           !
           klbea = klbea + 1
           !
           read( inmat, '( a )' ) bealin(klbea)
           !
        end do
        !
     else
        !
        read( inmat ) bealin(klbea)
        !
     end if
     !
  end do
  !
  kbeams = kbeams + nbeams
  !
  print*, 'nbeams :', nbeams, ' kbeams :', kbeams
  !
  print*, 'read', klbea, ' beam cards'
  !
  ! *************************************************************************************
  !
  ! - read detectors
  !
  ! *************************************************************************************
  !
  print*, 'ndetec :', ndetec
  !
  do idet = 1, ndetec
     !
     do i = 1, 3
        !
        kldet = kldet + 1
        !
        read( infil, '( a )' ) detlin(kldet)
        !
     end do
     !
  end do
  !
  kdetec = kdetec + ndetec
  !
  print*, 'ndetec :', ndetec, ' kdetec :', kdetec
  !
  print*, 'read', kldet, ' detector cards'
  !
  ! *************************************************************************************
  !
  ! - read planes
  !
  ! *************************************************************************************
  !
  do iplan = 1, nplan
     !
     if( mshfmt(1:1) .eq. 'f' ) then
        !
        read( inmsh, '( a, 2i5 )', end = 15 ) marg, ibc, ipts
        !
        klpla = klpla + 1
        !
        write( plalin(klpla), '( a, 2i5 )' ) ibc, ipts
        !
     else
        !
        read( inmsh, end = 15 ) ibc, ipts
        !
        klpla = klpla + 1
        !
        plalin(klpla) = 'u'
        !
        plaval(1,klpla) = ibc
        !
        plaval(2,klpla) = ipts
        !
        plaval(3,klpla) =  0
        !
        plaval(4,klpla) = -1
        !
     end if
     !
     do i = 1, ipts+1
        !
        klpla = klpla + 1
        !
        if( mshfmt(1:1) .eq. 'f' ) then
           !
           read( inmsh, '( a )', end = 15 ) plalin(klpla)
           !
        else
           !
           plalin(klpla) = 'u'
           !
           read( inmsh, end = 15 ) plaval(1,klpla),&
                &                                 plaval(2,klpla),&
                &                                 plaval(3,klpla)
           !
           plaval(4,klpla) = 0
           !
        end if
        !
     end do
     !
  end do
  !
15 kplan = kplan + nplan
  !
  print*, 'nplan  :', nplan, ' kplan  :', kplan
  !
  print*, 'read', klpla, ' plane cards'
  !
  ! *************************************************************************************
  !
  ! - read detector responses
  !
  ! *************************************************************************************
  !
  do iresp = 1, nresp
     !
     do igrp = 1, ngrps, 4
        !
        klrsp = klrsp + 1
        !
        if( matfmt(1:1) .eq. 'f' ) then
           !
           read( inmat, ' (a )', end = 5 ) rsplin(klrsp)
           !
        else
           !
           read( inmat, end = 5 ) rsplin(klrsp)
           !
        end if
        !
     end do
     !
  end do
  !
  kresp = kresp + nresp
  !
  print*, 'nresp  :', nresp, ' kresp  :', kresp
  !
  print*, 'read', klrsp, ' response cards'
  !
  ! *************************************************************************************
  !
  ! - read rest of input file
  !
  ! *************************************************************************************
  !
  do i = 1, 2000
     !
     read( infil, '( a )', end = 45 ) card
     !
     if( ifil .eq. 1 ) then
        !
        klrst = klrst + 1
        !
        rstlin(klrst) = card
        !
     end if
     !
  end do
  !
45 continue
  !
  print*, 'read', klrst, ' remaining cards'
  !
  go to 1
  !
  ! *************************************************************************************
  !
  ! - check for duplicate nodes and delete superfluous elements
  !
  ! *************************************************************************************
  !
50 nfil = ifil
  !
  jnod = 1
  !
  xran = abs( xmax - xmin )
  yran = abs( ymax - ymin )
  zran = abs( zmax - zmin )
  !
  range = max( xran, yran, zran )
  !
  if( yran .lt. eps ) yran = 1.0
  if( zran .lt. eps ) zran = 1.0
  !
  print*,'Setting up nodes from first file...'
  do inod = 1, nodend(1)
     newnod(inod) = inod
     oldnod(inod) = inod
     newcon(inod) = 0
  end do
  !
  jnod = nodend(1)
  !
  if( nfil .gt. 1 ) then
     !
     do ifil = 2, nfil
        !
        print*,'Merging common nodes of file ',ifil
        !
        jnodo = jnod
        !
        do inod = nodend(ifil-1)+1, nodend(ifil)
           !
           do mnod = 1, jnodo
              !
              difx = abs(x(mnod) - x(inod))
              dify = abs(y(mnod) - y(inod))
              difz = abs(z(mnod) - z(inod))
              !
              diff = sqrt( difx*difx + dify*dify + difz*difz )
              !
              if( difx .lt. eps*xran .and.&
                   &                dify .lt. eps*yran .and. &
                   &                difz .lt. eps*zran       ) then
                 !
                 print*, 'inod ', inod, ' is equal to node ', mnod
                 !
                 newnod(inod) = -abs(newnod(mnod))
                 newnod(mnod) = -abs(newnod(mnod))
                 !
                 go to 60
                 !
              end if
              !
           end do
           !
           jnod = jnod + 1
           !
           newnod(inod) = jnod
           oldnod(jnod) = inod
           newcon(jnod) = 0
           !
60         continue
           !
        end do
        !
     end do
     !
  end if
  !
  !      do 60 inod = 2, knod
  !c
  !c         if( inod .gt. 1 ) then
  !c
  !         do 55 mnod = 1, inod - 1
  !c
  !            difx = abs(x(mnod) - x(inod))
  !            dify = abs(y(mnod) - y(inod))
  !            difz = abs(z(mnod) - z(inod))
  !c
  !            diff = sqrt( difx*difx + dify*dify + dif*difz )
  !c
  !            if( difx .lt. eps*xran .and.
  !     :          dify .lt. eps*yran .and. 
  !     :          difz .lt. eps*zran       ) then
  !c
  !               print*, 'inod ', inod, ' is equal to node ', mnod
  !c
  !               newnod(inod) = -abs(newnod(mnod))
  !               newnod(mnod) = -abs(newnod(mnod))
  !c
  !               go to 60
  !c
  !            end if
  !c
  !   55    continue
  !*
  !c         end if
  !c
  !         jnod = jnod + 1
  !c
  !         newnod(inod) = jnod
  !         oldnod(jnod) = inod
  !         newcon(jnod) = 0
  !c
  !   60 continue
  !
  print*,'Finished merging nodes...'
  print*,'  Original total nodes   : ',knod
  print*,'  Nodes left after merge : ',jnod
  !      print*, 'knod :', knod, ' jnod :', jnod
  !
  jelm = 0
  !
  ! - surface elements first
  !
  print*,'Merging surface elements...'
  !
  do  ielm = 1, kelm! Was loop 80
     
    if( mod(elmtyp(ielm),2) .eq. 1 ) cycle
     
    do  i = 1, nln(ielm)! Was loop 65
        if( newnod(nodes(i,ielm)) .gt. 0 ) go to 70
    end do ! Was loop 65
     
    ! **** WARNING!!! ..............................
    ! - OK, we've found this surface has all nodes 'deleted', but that does
    ! - not necessarily mean this element should be ignored
          
    cycle
     
70  jelm = jelm + 1
     
    oldelm(jelm) = ielm
     
  end do ! Was loop 80
  !
  ! - volume elements
  !
  print*,'Merging volume elements...'
  !
  if( matelm ) kmat = 0
  !
  do  ielm = 1, kelm! Was loop 180
     
    if( mod(elmtyp(ielm),2) .eq. 0 ) cycle
     
    do i = 1, nln(ielm)
        if( newnod(nodes(i,ielm)) .gt. 0 ) go to 170
    end do
    
    ! **** WARNING!!! ..............................
    ! - OK, we've found this element has all nodes 'deleted', but that does
    ! - not necessarily mean this element should be ignored
     
    cycle
     
170 jelm = jelm + 1
     
    if( matelm ) then
        
        kmat         = kmat + 1
        
        elmmat(ielm) = kmat
        
    end if
     
    oldelm(jelm) = ielm
     
  end do ! Was loop 180
  !
  print*,'Finished merging elements...'
  print*,'  Original total elements   : ',kelm
  print*,'  Elements left after merge : ',jelm
  !
  print*,'Sorting out constraints...'
  !
  do icon = 1, kcon
     !
     inod = newnod(connod(icon))
     !
     if( inod .ne. 0 ) then
        !
        inodo = oldnod(abs(inod))
        !
        if( inodo .ne. connod(icon) ) then
           if( inod .lt. 0 ) goto 190
           print*,'*** INCONSISTENT OLD/NEW: ', connod(icon),&
                &                                              inodo, inod
        end if
        !
        inod = abs(inod)
        !
        if( newcon(inod) .ne. 0 ) then
           if( contyp(icon) .eq. 0 ) then
              print*,'*** CONFUSED?& !?!? ',icon,connod(icon),&
              &                                        contyp(icon),newcon(inod)
              print*,'  node co-ords: ',x(inodo),y(inodo),z(inodo)
           else if( contyp(icon) .lt. 0 ) then
              print*,'*** REPEAT CNSTR: ',icon,connod(icon),inod,&
                   &                                        contyp(icon),newcon(inod)
              print*,'  node co-ords: ',x(inodo),y(inodo),z(inodo)
              if( newcon(inod) .lt. 0 ) then
                 print*,'*** BOTH PERIODIC!?!?! -using newest'
                 newcon(inod) = -abs(newnod(-contyp(icon)))
              else
                 print*,'Ignoring reflective, using new periodic'
                 lstcon = newcon(inod)
                 newcon(inod) = -abs(newnod(-contyp(icon)))
                 inodo = -newcon(inod)
                 if( newcon(inodo) .eq. 0 ) then
                    print*,'Moving reflective to other node...'
                    newcon(inodo) = lstcon
                 else
                    print*,'*** OTHER NODE ALREADY CONSTRAINED:'
                    print*,inodo,newcon(inodo),lstcon
                 end if
              end if
           else if( contyp(icon) .gt. 0 ) then
              print*,'*** REPEAT CNSTR: ',icon,connod(icon),inod,&
                   &                                        contyp(icon),newcon(inod)
              print*,'  node co-ords: ',x(inodo),y(inodo),z(inodo)
              if( newcon(inod) .lt. 0 ) then
                 print*,'Keeping periodic condition (first one)'
                 inodo = -newcon(inod)
                 if( newcon(inodo) .eq. 0 ) then
                    print*,'Putting reflective on other node...'
                    newcon(inodo) = contyp(icon)
                 else
                    print*,'*** OTHER NODE ALREADY CONSTRAINED:'
                    print*,inodo,newcon(inodo),contyp(icon)
                 end if
              else
                 print*,'*** EXTRA REFLECTIONS! -using newest'
                 newcon(inod) = contyp(icon)
              end if
           end if
        else if( contyp(icon) .lt. 0 ) then
           newcon(inod) = -abs(newnod(-contyp(icon)))
        else if( contyp(icon) .gt. 0 ) then
           newcon(inod) = contyp(icon)
        else
           print*,'*** BAD CNSTR: ',icon,connod(icon),contyp(icon)
           print*,'  node co-ords: ',x(inodo),y(inodo),z(inodo)
        end if
        !
        !            print*,'old/new/typ: ',icon,connod(icon),inod,contyp(icon)
        !
     else
        !
        inodo = connod(icon)
        print*,'constraint removed: ',icon,connod(icon),contyp(icon)
        print*,'  node co-ords: ',x(inodo),y(inodo),z(inodo)
        !
     end if
     !
190 end do
  !
  jcon = 0
  !
  do inod = 1, jnod
     !
     if( newcon(inod) .ne. 0 ) jcon = jcon + 1
     !
  end do
  !      print*,'jcon: ',jcon
  !
  ! *************************************************************************************
  !
  ! - Make sure use unformatted for >=10000 nodes
  !
  ! *************************************************************************************
  !
  if( jnod .gt. 9999 .and. .not. binmsh ) then
     binmsh = .true.
     mshfmt = 'unformatted'
     print*,'Switching to unformatted mesh file (-lots of nodes)'
  end if
  !
  ! *************************************************************************************
  !
  ! - output mesh and material cards for target file
  !
  ! *************************************************************************************
  !
  if( .not. binmsh )then
     !
     write( outfil, '( 3(a,2x) )' ) '''MESHFILE''',&
          &                                     ''''//trgtnam(1:kt)//'''',&
          &                                     quote//mshfmt//quote
     !
  else
     !
     write( outfil, '( 3(a,2x) )' ) '''MESHFILE''',&
          &                             ''''//trgtnam(1:kt)//'.bin'//'''',&
          &                                     quote//mshfmt//quote
     !
  end if
  !
  if( .not. binmat )then
     !
     write( outfil, '( 4(a,2x) )' ) '''MATFILE''',&
          &                                 ''''//trgtnam(1:kt)//'''',&
          &                            '''event''', quote//matfmt//quote
     !
  else
     !
     write( outfil, '( 4(a,2x) )' ) '''MATFILE''',&
          &                            ''''//trgtnam(1:kt)//'.bin'//'''',&
          &                            '''event''', quote//matfmt//quote
     !
  end if
  !
  ! *************************************************************************************
  !
  ! - output merged mesh
  !
  ! *************************************************************************************
  !
  print*, 'nsurfr :', nsurfr, ' raytrace : ', raytrc
  !
  if( matelm .and. samemat ) then
     !
     kiso = kcr
     !
  end if
  !
  if( .not. binmsh ) then
     !
     open( outmsh, file = path(1:kp)//'mesh/'//trgtnam(1:kt), &
          &        status = 'unknown', form = 'formatted' )
     !
     PRINT*,'---00000',kreg, kmat+kiso,&
          &                                     ksor, jelm, jnod, jcon,&
          &                                     kalb, nopt, kplan, nprocs,&
          &                                     nprocno, ksurfr,0,.false.
     !        write( outmsh, '(a, 5x, 14i5  )' ) 'SIZES', kreg, kmat+kiso,
     write( outmsh, '(a,5x, 13i5, 1X,l1 )') 'SIZES', kreg, kmat+kiso,&
          &                                     ksor, jelm, jnod, jcon,&
          &                                     kalb, nopt, kplan, nprocs,&
          &                                     nprocno, ksurfr,0,.false.
     !     :                                     nprocno, ksurfr
!!!!!!!!!!!!!!!! changed this Jeff + Chris
  else
     !
     open( outmsh, file = path(1:kp)//'mesh/'//trgtnam(1:kt)//'.bin', &
          &        status = 'unknown', form = 'unformatted' )
     !
     write( outmsh ) kreg, kmat+kiso, ksor, jelm, jnod, jcon,&
          &                  kalb, nopt, kplan, nprocs, nprocno, ksurfr
     !
  end if
  !
  print*,'Writing ',jelm,' merged elements...'
  !
  do  ielm = 1, jelm! Was loop 100
     !
     kelm = oldelm(ielm)
     !
     imat = elmmat(kelm)
     !
     if( .not. binmsh ) then
        !
        write( outmsh, '(a, 15i5 )' ) 'elem ',ielm,&
             &                                    elmtyp(kelm), elmreg(kelm),&
             &                                    elmmat(kelm), elmsrc(kelm),&
             &                  (abs(newnod(nodes(j,kelm))), j = 1, nln(kelm) )
        !
     else
        !
        write( outmsh ) ielm, elmtyp(kelm), elmreg(kelm),&
             &                      elmmat(kelm), elmsrc(kelm), nln(kelm),&
             &               (abs(newnod(nodes(j,kelm))), j = 1, nln(kelm) )
        !
     end if
     !
  end do ! Was loop 100
  !
  ! *************************************************************************************
  !
  ! - output merged nodes
  !
  ! *************************************************************************************
  !
  print*,'Writing ',jnod,' merged nodes...'
  !
  do inod = 1, jnod
     !
     knod = oldnod(inod)
     !
     if( .not. binmsh ) then
        !
        write( outmsh, '(a, i5, 3(1pe15.7), 4i5 )' ) 'node ', inod,&
             &                   x(knod), y(knod), z(knod), inod, 0, 1, 0
        !
     else
        !
        write( outmsh ) x(knod), y(knod), z(knod), inod, 0, 1, 0
        !
     end if
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output merged constraints
  !
  ! *************************************************************************************
  !
  icon = 0
  !
  print*,'Writing ',jcon,' constraints (in theory)...'
  !      pause
  !
  do inod = 1, jnod
     !
     if( newcon(inod) .ne. 0 ) then
        !
        icon = icon + 1
        !
        if( .not. binmsh ) then
           !
           write( outmsh, '( a, 3i5 )' ) 'cons ', icon, &
                &                                       inod, newcon(inod)
           !
        else
           !
           write( outmsh ) inod, newcon(inod)
           !
        end if
        !
        !            write( *, '( a, 3i8 )' ) 'CONS', icon, 
        !     :                               inod, newcon(inod)
        !
        !            if( newcon(inod) .eq. 10 ) then
        !               print*,'Corner node: ',icon,inod,newcon(inod)
        !               inodo = oldnod(inod)
        !               print*,'  co-ords: ',x(inodo),y(inodo),z(inodo)
        !            end if
        !
        if( newcon(inod) .lt. 0 ) then
           if( newcon(inod) .lt. -inod ) then
              print*,'*** PERIODIC CONSTRAINT WRONG WAY ROUND!'
              write( *, '( a, 3i8 )' ) 'CONS', icon, &
                   &                                     inod, newcon(inod)
           end if
        end if
        !
     end if
     !
  end do
  !
  if( jcon .ne. icon ) then
     print*,'Actually wrote ',icon,' constraints'
     print*,'  (Was supposed to be ',jcon,' )'
     print*,'  *** MESH FILE SIZES PARAM WILL BE INCONSISTENT ***'
  else
     print*,'Found correct no. of constraints: ',icon
  end if
  !
  ! *************************************************************************************
  !
  ! - open output material file
  !
  ! *************************************************************************************
  !
  if( .not. binmat ) then
     !
     open( outmat, file = path(1:kp)//'matxs/'//trgtnam(1:kt), &
          &         status = 'unknown', form = 'formatted' )
     !
  else
     !
     open( outmat,&
          &         file = path(1:kp)//'matxs/'//trgtnam(1:kt)//'.bin', &
          &         status = 'unknown', form = 'unformatted' )
     !
  end if
  !
  ! *************************************************************************************
  !
  ! - write material card
  !
  ! *************************************************************************************
  !
  if( kmix .lt. 10000 ) then
     write( outfil, '( a, 3x, 12i5 )' ) '''MATXS''', kcr, mtps, kmix,&
          &                           iht, ihs, ihm, datfmt, kden, kbeams,&
          &                           kdetec, ksrcs, kresp
  else if( kmix .lt. 1000000 ) then
     print*,'#####:',kcr, mtps, kmix,&
          &                           iht, ihs, ihm, datfmt, kden, kbeams,&
          &                           kdetec, ksrcs, kresp
     write( outfil, '( a, 3x, 12i7 )' ) '''MATXS''', kcr, mtps, kmix,&
          &                           iht, ihs, ihm, datfmt, kden, kbeams,&
          &                           kdetec, ksrcs, kresp
  else
     write( outfil, * ) '''MATXS''', kcr, mtps, kmix,&
          &                           iht, ihs, ihm, datfmt, kden, kbeams,&
          &                           kdetec, ksrcs, kresp
  end if
  !
  ! *************************************************************************************
  !
  ! - output comment cards
  !
  ! *************************************************************************************
  !
  write( outfil, '( a, i5 )' ) '''NCOM    ''', kcom
  !
  do ilcom = 1, klcom
     !
     write(  outfil, '( a )' ) comlin(ilcom)
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output mixing instructions
  !
  ! *************************************************************************************
  !
   print*, 'Writing ', kmix,' mixing instructions...'

   maxnum = 0

   do k = 1, kmix

      maxnum = max( maxnum, mixnum(k) )
         
      if( mixcom(k) .gt. 0 ) mixcom(k) = mixcom(k) + kmat

! such that formatting for KMIX .GE. 100,000 is correct - and changed format slightly
! to match getdat.f in EVENT
      if(kmix .ge. 100000) then
          write( outfil, '( a,i11,2i10,1pe15.7 )' )'MIX',k,mixnum(k),&
            &   mixcom(k),&
            &   mixden(k)   
      else
          write( outfil, '( a,i7,2i5,1pe15.7 )' )'MIX',k,mixnum(k),&
            &   mixcom(k),&
            &   mixden(k)         
      end if  
         

   end do

   if( maxnum .ne. kmat ) then

      print*, '******'

      print*, 'maxnum :', maxnum, ' kmat :', kmat
 
      print*, '******'

   end if

  !
  ! *************************************************************************************
  !
  ! - output merged materials
  !
  ! *************************************************************************************
  !
  print*, 'Writing ', klmat,' merged material lines...'
  !
  do ilmat = 1, klmat
     !
     card = matlin(ilmat)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     write( outmat, '( a )' ) card(1:i)
     !
  end do
  !
  ! ******************************************************************************
  !
  ! -  output velocities
  !
  ! ******************************************************************************
  !
  print*, 'Writing ', klvel,' velocity lines...'
  !
  do ilvel = 1, klvel
     !
     card = vellin(ilvel)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     if( .not. binmat ) then
        !
        write( outmat, '( a )' ) card(1:i)
        !
     else
        !
        write( outmat ) vellin(ilvel)
        !
     end if
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output merged sources
  !
  ! *************************************************************************************
  !
  print*, 'Writing ', klsor,' merged source lines...'
  !
  do ilsor = 1, klsor
     !
     card = sorlin(ilsor)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     if( .not. binmat ) then
        !
        write( outmat, '( a )' ) card(1:i)
        !
     else
        !
        write( outmat ) card(1:i)
        !
     end if
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output merged albedos
  !
  ! *************************************************************************************
  !
  print*, 'Writing ', klalb,' merged albedo lines...'
  !
  do ilalb = 1, klalb
     !
     card = alblin(ilalb)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     if( .not. binmat ) then
        !
        read( outmat, ' (a )' ) card(1:i)
        !
     else
        !
        read( outmat ) card(1:i)
        !
     end if
     !
  end do
  !
  ! ******************************************************************************
  !
  ! -  output fission spectrum
  !
  ! ******************************************************************************
  !
  print*, 'Writing ', klspc,' fission spectrum lines...'
  !
  do ilspc = 1, klspc
     !
     card = spclin(ilspc)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     if( .not. binmat ) then
        !
        write( outmat, '( a )' ) card(1:i)
        !
     else
        !
        write( outmat ) card(1:i)
        !
     end if
     !
  end do
  !
  ! ******************************************************************************
  !
  ! -  output delayed neutron parameters
  !
  ! ******************************************************************************
  !
  print*, 'Writing ', kldel,' delayed fission lines...'
  !
  do ildel = 1, kldel
     !
     card = dellin(ildel)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     write( outmat, '( a )' ) card(1:i)
     !
  end do
  !
  ! ******************************************************************************
  !
  ! - output broad group structure
  !
  ! ******************************************************************************
  !
  do ilbgr = 1, klbgr
     !
     card = bgrlin(ilbgr)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     write( outfil, '( a )' ) card(1:i)
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output merged general source
  !
  ! *************************************************************************************
  !
  print*, 'Writing ', klsrc,' merged general source lines...'
  !
  do ilsrc = 1, klsrc
     !
     card = srclin(ilsrc)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     if( .not. binmat ) then
        !
        write( outmat, '( a )' ) card(1:i)
        !
     else
        !
        write( outmat ) card(1:i)
        !
     end if
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output beams
  !
  ! *************************************************************************************
  !
  do ilbea = 1, klbea
     !
     card =  bealin(ilbea)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     if( .not. binmsh ) then
        !
        write( outmat, '( a )' ) card(1:i)
        !
     else
        !
        write( outmat ) card(1:i)
        !
     end if
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output detectors
  !
  ! *************************************************************************************
  !
  do ildet = 1, kldet
     !
     card = detlin(ildet)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     write( outfil, '( a )' ) card(1:i)
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output planes
  !
  ! *************************************************************************************
  !
  print*, 'Writing ', klpla,' merged plane lines...'
  !
  do ilpla = 1, klpla
     !
     if( plalin(klpla) .ne. 'u' ) then
        !
        write( outfil, '( a )' ) plalin(klpla)
        !
     else
        !
        if( plaval(4,ilpla) .eq. -1 ) then
           !
           ibc  = plaval(1,ilpla)
           !
           ipts = plaval(2,ilpla)
           !
           write( outmsh ) ibc, ipts
           !
        else
           !
           write( outmsh ) plaval(1,ilpla), plaval(2,ilpla),&
                &                         plaval(3,ilpla)
           !
        end if
        !
     end if
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output detector responses
  !
  ! *************************************************************************************
  !
  do ilrsp = 1, klrsp
     !
     card = rsplin(ilrsp)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     if( .not. binmat ) then
        !
        write( outmat, ' (a )' ) card(1:i)
        !
     else
        !
        write( outmat ) card(1:i)
        !
     end if
     !
  end do
  !
  ! *************************************************************************************
  !
  ! - output remainder of input file
  !
  ! *************************************************************************************
  !
  do ilrst = 1, klrst
     !
     card = rstlin(ilrst)
     !
     do i = 80, 1, -1
        !
        if( card(i:i) .ne. ' ' ) exit
        !
     end do
     !
     write( outfil, '( a )' ) card(1:i)
     !
  end do
  !
131 continue
  !
  close( outfil )
  close( outmsh )
  close( outmat )
  !
  write( *, * )
  write( *, * )'***************************************************'
  write( *, * )'***         Program evmerge finished ok         ***'
  write( *, * )'***************************************************'
  write( *, * )
  stop
end program evmerge
