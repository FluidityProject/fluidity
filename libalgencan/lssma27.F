#include "lma27ad.f"

C     lssinfo:
C
C     0: Success.
C     1: Matrix not positive definite.
C     2: Rank deficient matrix.
C     6: Insufficient space to store the linear system.
C     7: Insufficient double precision working space.
C     8: Insufficient integer working space.

C     ******************************************************************
C     ******************************************************************

      logical function lss(lsssub)

C     SCALAR ARGUMENTS
      character * 4 lsssub

      lsssub = 'MA27'
      lss    = .true.

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssini(sclsys,acneig,usefac)

      implicit none

C     SCALAR ARGUMENTS
      logical acneig,sclsys,usefac

#include "dim.par"
#include "ma27dat.com"
#include "machconst.com"

      lacneig = acneig
      lusefac = usefac
      lsclsys = sclsys

      call ma27id(icntl,cntl)

C     Suppress monitoring, warning and error messages
      icntl(1) = 0
      icntl(2) = 0

C     Set the appropriate options in the case that the factorization
C     should stop when a nonpositive eigenvalue is found
      if ( .not. lacneig ) then
          cntl(1) =   - 1.0d0
          cntl(3) = macheps23
      end if
      
C     Initialize scaling module
      if ( lsclsys ) then
          call sclini()
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssana(nsys,hnnz,hlin,hcol,hval,hdiag,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,lssinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz)

#include "dim.par"
#include "outtyp.com"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer sclinfo

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(file10_unit,9000) nsysmax,nsys
          end if
          
          lssinfo = 6
          return
      end if

      if ( lsclsys ) then
          call sclana(nsys,hnnz,hlin,hcol,hval,hdiag,sclinfo)
          
          if ( sclinfo .ne. 0 ) then
            ! SCALING FAILURE
            
              lsclsys = .false.
            
              if ( iprintctl(3) ) then
                  write(* ,9010)
                  write(file10_unit,9010)
              end if
          end if

      end if

      call ma27ad(nsys,hnnz,hlin,hcol,iw,wintmax,ikeep,iw1,nsteps,0,
     +icntl,cntl,info,ops)

      if ( info(1) .eq. 0 ) then
        ! SUCCESS

          lssinfo = 0
          return
      end if

C     UNHANDLED ERROR

      if ( iprintctl(3) ) then
          write(* ,9020) info(1)
          write(file10_unit,9020) info(1)
      end if

      stop

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSANA-MA27 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9010 format(/,1X,'LSSANA-MA27 WARNING: Scaling of the system failed.',
     +       /,1X,'Scaling is being turned off.')
 9020 format(/,1X,'LSSANA-MA27 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssfac(nsys,hnnz,hlin,hcol,hval,hdiag,d,pind,pval,
     +nneigv,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,pind,nneigv,lssinfo
      double precision pval

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz),d(nsys)

#include "dim.par"
#include "outtyp.com"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer i,idiag,kb,kd,kw,nblk,ncols,nrows,block,sclinfo
      double precision d1

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(file10_unit,9000) nsysmax,nsys
          end if
          
          lssinfo = 6
          return
      end if

      do i = 1,hnnz
          fact(i) = hval(i)
      end do
      
      do i = 1,nsys
          idiag = hdiag(i)
          fact(idiag) = fact(idiag) + d(i)
      end do
      
      if ( lsclsys ) then
          
          call sclsys(nsys,hnnz,hlin,hcol,fact,s,sclinfo)
          
          if ( sclinfo .ne. 0 ) then
            ! SCALING FAILURE
            
              lsclsys = .false.
              
              if ( iprintctl(3) ) then
                  write(* ,9010)
                  write(file10_unit,9010)
              end if
          end if
          
      end if
      
      call ma27bd(nsys,hnnz,hlin,hcol,fact,fnnzmax,iw,wintmax,ikeep,
     +nsteps,maxfrt,iw1,icntl,cntl,info)

      if ( info(1) .eq. 0 ) then

          d1 = hval(hdiag(1)) + d(1)

          if ( d1 .gt. 0.0d0 ) then
            ! SUCCESS

              lssinfo = 0
          else
            ! MATRIX IS NEGATIVE DEFINITE

              pind    = 1
              pval    = abs( d1 )
              
              lssinfo = 1
          end if

      else if ( info(1) .eq. 2 .or. info(1) .eq. -6 ) then
        ! MATRIX NOT POSITIVE DEFINITE

          pind    = info(2)
          pval    = abs( fact(info(16)) )
          lssinfo = 1
          
      else if ( info(1) .eq. 3 .or. info(1) .eq. -5 ) then
        ! RANK DEFICIENT MATRIX

          pind    = info(2)
          pval    = abs( fact(info(16)) )
          lssinfo = 2
          
      else if ( info(1) .eq. -4 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          if ( iprintctl(3) ) then
              write(* ,9020) fnnzmax,info(2)
              write(file10_unit,9020) fnnzmax,info(2)
          end if
          
          lssinfo = 7

      else if ( info(1) .eq. -3 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          if ( iprintctl(3) ) then
              write(* ,9030) wintmax,info(2)
              write(file10_unit,9030) wintmax,info(2)
          end if
          
          lssinfo = 8

      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9040) info(1)
              write(file10_unit,9040) info(1)
          end if
          
          stop

      end if

C     NUMBER OF NEGATIVE EIGENVALUES

      nneigv = info(15)
      
      if ( lusefac .and. lssinfo .eq. 0 ) then

C         USED WHEN A PARTIAL SYSTEM IS SOLVED

          latop  = info(9)

C         Define matrix D^{-1} (stored in sdiag)

          kb = 1
          kd = 0
          kw = 2
          
          do i = 1,nsys
              sdiag(i) = 0.0d0
          end do
          
          nblk = abs( iw(1) )
          
          do block = 1,nblk
              
              ncols = iw(kw)
              if ( ncols .gt. 0 ) then
                  kw    = kw + 1
                  nrows = iw(kw)
              else
                  ncols = -ncols
                  nrows = 1
              end if
              
              do i = 1,nrows
                  kd        = kd + 1
                  sdiag(kd) = fact(kb)
                  
                  kb = kb + ncols + 1 - i
              end do
              
              kw = kw + ncols + 1
          end do
          
          if ( lsclsys ) then
              do i = 1,nsys
                  sdiag(ikeep(i)) = sdiag(ikeep(i)) * s(i)**2
              end do
          end if

      end if

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSFAC-MA27 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9010 format(/,1X,'LSSFAC-MA27 WARNING: Scaling of the system failed. ',
     +       /,1X,'Scaling is being turned off.')
 9020 format(/,1X,'LSSFAC-MA27 WARNING: Insufficient double precision ',
     +            'working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSFAC-MA27 WARNING: Insufficient integer working ',
     +            'space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9040 format(/,1X,'LSSFAC-MA27 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssol(nsys,sol)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      double precision sol(nsys)

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer i

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if
      
      call ma27cd(nsys,fact,fnnzmax,iw,wintmax,w,maxfrt,sol,iw1,nsteps,
     +icntl,info)

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssoltr(job,nsys,sol)

      implicit none

C     SCALAR ARGUMENTS
      character * 1 job
      integer nsys

C     ARRAY ARGUMENTS
      double precision sol(nsys)

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALAR
      integer i,nblk

      if ( job .eq. 'T' .or. job .eq. 't' ) then
          
          if ( lsclsys ) then
              do i = 1,nsys
                  sol(i) = sol(i) * s(i)
              end do
          end if

          nblk = abs( iw(1) )
          
          call ma27qd(nsys,fact,fnnzmax,iw(2),wintmax-1,w,maxfrt,sol,
     +    iw1,nblk,latop,icntl)

          if ( lsclsys ) then
              do i = 1,nsys
                  sol(i) = sol(i) / s(i)
              end do
          end if

          do i = 1,nsys
              w(i) = sol(i)
          end do

          do i = 1,nsys
              sol(ikeep(i)) = w(i)
          end do

          do i = 1,nsys
              sol(i) = sol(i) * sqrt( sdiag(i) )
          end do

      else
          
          do i = 1,nsys
              sol(i) = sol(i) / sqrt ( sdiag(i) )
          end do

          do i = 1,nsys
              w(i) = sol(i)
          end do

          do i = 1,nsys
              sol(i) = w(ikeep(i))
          end do

          if ( lsclsys ) then
              do i = 1,nsys
                  sol(i) = sol(i) * s(i)
              end do
          end if

          nblk = abs( iw(1) )
          
          call ma27rd(nsys,fact,fnnzmax,iw(2),wintmax-1,w,maxfrt,sol,
     +    iw1,nblk,latop,icntl)

          if ( lsclsys ) then
              do i = 1,nsys
                  sol(i) = sol(i) * s(i)
              end do
          end if

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsspermvec(nsys,v)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      integer v(nsys)

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      integer wi(nsysmax)

      do i = 1,nsys
          wi(ikeep(i)) = v(i)
      end do
      
      do i = 1,nsys
          v(i) = wi(i)
      end do
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssunpermvec(nsys,v)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      integer v(nsys)

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      integer wi(nsysmax)

      do i = 1,nsys
          wi(i) = v(ikeep(i))
      end do
      
      do i = 1,nsys
          v(i) = wi(i)
      end do
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsspermind(hnnz,hlin,hcol)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer col,i,lin

      do i = 1,hnnz
          col = hcol(i)
          lin = hlin(i)
          
          hcol(i) = ikeep(col)
          hlin(i) = ikeep(lin)
      end do
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssunpermind(nsys,hnnz,hlin,hcol)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer col,i,lin

      do i =1,nsys
          invp(i)=0
      end do
      
      do i = 1,nsys
          invp(ikeep(i)) = i
      end do
      
      do i = 1,hnnz
          col = hcol(i)
          lin = hlin(i)
          
          hcol(i) = invp(col)
          hlin(i) = invp(lin)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      double precision function lssgetd(j)

      implicit none

C     SCALAR ARGUMENTS
      integer j

#include "dim.par"
#include "ma27dat.com"

      lssgetd = sqrt( 1.0d0 / sdiag(j) )

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssetrow(nsys)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer i,ii,j,jj,j1,j2,nblk,apos,iblk,ilvl,ipiv,ipos,irow,ist,
     +        liell,npiv

      integer    ifrlvl
      parameter (ifrlvl = 5)
      
C     Array of null elements

      do i = 1,nsys
          pos(i) = 0
      end do

C     Compute inverse permutation

      do i = 1,nsys
          invp(ikeep(i)) = i
      end do

C     Set sta and next

      nblk = abs( iw(1) )
      
      apos   = 1
      ipos   = 1
      j2     = 0
      iblk   = 0
      npiv   = 0

      do j = 1,nsys
          sta(j) = 0
      end do
      
      do j = 1,nsys
          curr(j) = 0
      end do

      do 140 irow = 1,nsys
          if ( npiv .gt. 0 ) go to 90
          
          iblk = iblk + 1
          
          if ( iblk .gt. nblk ) go to 150
          
          ipos  = j2 + 1
          liell = -iw(ipos+1)
          npiv  = 1
          
          if ( liell .le. 0 ) then
              liell = -liell
              ipos  = ipos + 1
              npiv  = iw(ipos+1)
          end if
          
          j1   = ipos + 1
          j2   = ipos + liell
          ilvl = min( npiv, 10 )
          
          if ( liell .lt. icntl(ifrlvl+ilvl) ) go to 90

          do jj = j1,j2
              
              j  = abs( iw(jj+1)+0 )
              
              if ( j .ne. 0 ) then
                  
                  if ( sta(j) .eq. 0 ) then
                      sta(j)   = jj
                      curr(j)  = jj
                      next(jj) = 0
                  else
                      next(curr(j)) = jj
                      next(jj)      = 0
                      curr(j)       = jj
                  end if
                  
              end if
          end do
          
          do ipiv = 1,npiv
              apos = apos + 1
              ist  = ipiv + 1
              apos = apos + liell - ist + 1
          end do
          
          npiv = 0
          go to 140
          
 90       npiv = npiv - 1
          apos = apos + 1
          j1   = j1 + 1
          
          ii   = iw(j1)
          
          if ( ii .ne. 0 ) then
              
              if ( sta(ii) .eq. 0 ) then
                  sta(ii)    = j1-1
                  curr(ii)   = j1-1
                  next(j1-1) = 0
              else
                  next(curr(ii)) = j1-1
                  next(j1-1)     = 0
                  curr(ii)       = j1-1
              end if
              
          end if
          
          apos = apos + j2 - j1 + 1
 140  continue
      
      latop = apos - 1
 150  continue
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssgetrow(nsys,idx,rownnz,rowind,rowval)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,idx,rownnz

C     ARRAY ARGUMENTS
      integer rowind(nsys)
      double precision rowval(nsys)

#include "dim.par"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer i,ii,j,jj,j1,j2,k,apos,iblk,ifr,ilvl,ipiv,ipos,icurr,ist,
     +        liell,npiv,iiw1,prev,strt,nblk,ind
      double precision sclf

      integer    ifrlvl
      parameter (ifrlvl = 5)

      rownnz = 0
      apos   = 1
      ipos   = 1
      j2     = 0
      iblk   = 0
      npiv   = 0
      nblk   = abs( iw(1) )
      ind    = invp(idx)
      icurr  = sta(ind)

      if ( lsclsys ) then
          sclf = s(ind)
      else
          sclf = 1.0d0
      end if

 140  if ( icurr .ne. 0 ) then
          
          if ( npiv .gt. 0 ) go to 90
          
          iblk = iblk + 1
          
          if ( iblk .gt. nblk ) go to 150
          
          ipos      = j2 + 1
          iw1(iblk) = ipos
          liell     = -iw(ipos+1)
          npiv      = 1
          
          if ( liell .le. 0 ) then
              liell = -liell
              ipos  = ipos + 1
              npiv  = iw(ipos+1)
          end if
          
          j1   = ipos + 1
          j2   = ipos + liell
          ilvl = min( npiv, 10 )
          
          if ( liell .lt. icntl(ifrlvl+ilvl) ) go to 90
          
          ifr  = 0
          iiw1 = 0
          prev = icurr
          strt = j1
          
 100      if ( ( icurr .ne. 0 ) .and. ( icurr .ge. j1 ) .and.
     +         ( icurr .le. j2 ) ) then
              
              do jj = strt,icurr-1
                  
                  j   = abs( iw(jj+1)+0 )
                  ifr = ifr + 1
                  
                  if ( pos(j) .eq. 0 ) then
                      w(ifr) = 0.0d0
                  else
                      w(ifr) = rowval(pos(j))
                  end if
              end do
              
              j    = abs( iw(icurr+1)+0 )
              ifr  = ifr + 1
              iiw1 = ifr
              
              if ( pos(j) .eq. 0 ) then
                  w(ifr) = 0.0d0
              else
                  w(ifr) = rowval(pos(j))
              end if
              
              strt  = icurr + 1
              prev  = icurr
              icurr = next(icurr)
              go to 100
          end if
          
          if ( iiw1 .eq. 0 ) then
              go to 10
          end if
          
          do jj = prev+1,j2
              
              j   = abs( iw(jj+1)+0 )
              ifr = ifr + 1
              
              if ( ( j .eq. ind ) .and. ( j .ne. 0 ) ) then
                  iiw1 = ifr
              end if
              
              if ( pos(j) .eq. 0 ) then
                  w(ifr) = 0.0d0
              else
                  w(ifr) = rowval(pos(j))
              end if
          end do
          
          do ipiv = 1,npiv
              apos = apos + 1
              ist  = ipiv + 1
              
              if ( ipiv .eq. iiw1 ) then
                  k        = apos
                  w(ist-1) = w(ist-1) + sclf
                  
                  do j = ist,liell
                      w(j) = w(j) - sclf * fact(k)
                      k    = k + 1
                  end do
              end if
              
              apos = apos + liell - ist + 1
          end do
          
          ifr = 0
          do jj = j1,j2
              
              j   = abs( iw(jj+1)+0 )
              ifr = ifr + 1
              
              if ( pos(j) .eq. 0 ) then
                  rownnz         = rownnz + 1
                  rowind(rownnz) = j
                  rowval(rownnz) = w(ifr)
                  pos(j)         = rownnz
              else
                  rowval(pos(j)) = w(ifr)
              end if
          end do
          
 10       npiv = 0
          go to 140
          
 90       npiv = npiv - 1
          apos = apos + 1
          j1   = j1 + 1
          ii   = iw(j1)
          
          if ( ( icurr .ne. 0 ) .and. ( icurr .eq. j1-1 ) ) then
              k = apos
              
              if ( pos(ii) .eq. 0 ) then
                  rownnz          = rownnz + 1
                  rowind(rownnz)  = ii
                  rowval(rownnz)  = sclf
                  pos(ii)         = rownnz
              else
                  rowval(pos(ii)) = rowval(pos(ii)) + sclf
              end if
              
              do j = j1,j2
                  
                  ii = abs( iw(j+1)+0 )
                  
                  if ( pos(ii) .eq. 0 ) then
                      rownnz          = rownnz + 1
                      rowind(rownnz)  = ii
                      rowval(rownnz)  = - sclf * fact(k)
                      pos(ii)         = rownnz
                  else
                      rowval(pos(ii)) = rowval(pos(ii)) - sclf * fact(k)
                  end if
                  k = k + 1
              end do
              
              icurr = next(icurr)
          end if
          
          apos = apos + j2 - j1 + 1
          
          go to 140
      end if
      
 150  continue
      
      if ( lsclsys ) then
          do i = 1,rownnz
              rowval(i) = rowval(i) / s(rowind(i))
          end do
      end if

      do j = 1,rownnz
          pos(rowind(j)) = 0
      end do
      
      do i = 1,rownnz
          rowval(i) = rowval(i) / sqrt( sdiag(idx) )
      end do
      
      do i = 1,rownnz
          rowind(i) = ikeep(rowind(i))
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssafsol(nsys,hnnz,hlin,hcol,hval,hdiag,d,sol,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,lssinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz),d(nsys),sol(nsys)

#include "dim.par"
#include "outtyp.com"
#include "ma27dat.com"

C     LOCAL SCALARS
      integer i,nstepstmp,maxfrttmp
      double precision opstmp
      
C     LOCAL ARRAYS
      integer iwtmp(wintmax),iw1tmp(2*nsysmax),ikeeptmp(3*nsysmax)

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM
          
          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(file10_unit,9000) nsysmax,nsys
          end if
          
          lssinfo = 6
          return
      end if

      call ma27ad(nsys,hnnz,hlin,hcol,iwtmp,wintmax,ikeeptmp,iw1tmp,
     +nstepstmp,0,icntl,cntl,info,opstmp)

      if ( info(1) .eq. 0 .or. info(1) .eq. 1 ) then
        ! SUCCESS

          lssinfo = 0

      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) info(1)
              write(file10_unit,9030) info(1)
          end if
          
          stop
          
      end if
      
      do i = 1,hnnz
          fact(i) = hval(i)
      end do
      
      do i = 1,nsys
          fact(hdiag(i)) = fact(hdiag(i)) + d(i)
      end do

      if ( lsclsys ) then
          do i = 1,hnnz
              fact(i) = fact(i) * s(hlin(i)) * s(hcol(i))
          end do
      end if
      
      call ma27bd(nsys,hnnz,hlin,hcol,fact,fnnzmax,iwtmp,wintmax,
     +ikeeptmp,nstepstmp,maxfrttmp,iw1tmp,icntl,cntl,info)

      if ( info(1) .eq. 0 .or. info(1) .eq. 1 ) then
        ! SUCCESS

          lssinfo = 0

      else if ( info(1) .eq. 2 .or. info(1) .eq. -6 ) then
        ! MATRIX NOT POSITIVE DEFINITE
          
          lssinfo = 1

      else if ( info(1) .eq. 3 .or. info(1) .eq. -5 ) then
        ! RANK DEFICIENT MATRIX

          lssinfo = 2

      else if ( info(1) .eq. -4 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          lssinfo = 7

          if ( iprintctl(3) ) then
              write(* ,9010) fnnzmax,info(2)
              write(file10_unit,9010) fnnzmax,info(2)
          end if
          
          return
          
      else if ( info(1) .eq. -3 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          lssinfo = 8
          
          if ( iprintctl(3) ) then
              write(* ,9020) wintmax,info(2)
              write(file10_unit,9020) wintmax,info(2)
          end if
          
          return
          
      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) info(1)
              write(file10_unit,9030) info(1)
          end if
          
          stop
          
      end if

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if

      call ma27cd(nsys,fact,fnnzmax,iwtmp,wintmax,w,maxfrttmp,sol,
     +iw1tmp,nstepstmp,icntl,info)

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSAFSOL-MA27 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9010 format(/,1X,'LSSAFSOL-MA27 WARNING: Insufficient double ',
     +            'precision working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9020 format(/,1X,'LSSAFSOL-MA27 WARNING: Insufficient integer ',
     +            'working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSAFSOL-MA27 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end

