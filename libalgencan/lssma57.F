#include "lma57ad.f"

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

      lsssub = 'MA57'
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
#include "ma57dat.com"
#include "machconst.com"

      lacneig = acneig
      lusefac = usefac
      lsclsys = sclsys

      call ma57id(cntl,icntl)

C     Suppress monitoring, warning and error messages
      icntl(5) = 0

C     Chooses AMD pivot ordering using MC47
C     icntl(6) = 0

      if ( lsclsys ) then
          icntl(15) = 1
      else
          icntl(15) = 0
      end if
      
      if ( .not. lacneig ) then
          icntl(7)  = 2
          icntl(8)  = 1
          
          cntl(2)   = macheps23
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
#include "ma57dat.com"

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(file10_unit,9000) nsysmax,nsys
          end if
          
          lssinfo = 6
          return
      end if

      call ma57ad(nsys,hnnz,hlin,hcol,wintmax,keep,iwork,icntl,info,
     +rinfo)

      if ( info(1) .eq. 0 ) then
        ! SUCCESS

          lssinfo = 0
          return
      end if

C     UNHANDLED ERROR

      if ( iprintctl(3) ) then
          write(* ,9030) info(1)
          write(file10_unit,9030) info(1)
      end if
      
      stop

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSANA-MA57 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSANA-MA57 ERROR: Unhandled error ',I16,'.',
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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer i,idiag
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

      do i = 1,nsys
          idiag       = hdiag(i)
          w(i)        = hval(idiag)
          hval(idiag) = hval(idiag) + d(i)
      end do

      call ma57bd(nsys,hnnz,hval,fact,fnnzmax,ifact,wintmax,wintmax,
     +keep,iwork,icntl,cntl,info,rinfo)

      if ( lsclsys ) then
          do i = 1,nsys
              s(i) = fact(fnnzmax-nsys-1+i)
          end do    
      end if

      do i = 1,nsys
          idiag       = hdiag(i)
          hval(idiag) = w(i)
      end do
      
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

      else if ( info(1) .eq. 5 .or. info(1) .eq. -6 ) then
        ! MATRIX NOT POSITIVE DEFINITE

          pind    = info(2)
          pval    = abs( rinfo(20) )
          
          lssinfo = 1

      else if ( info(1) .eq. 4 .or. info(1) .eq. -5 ) then
        ! RANK DEFICIENT MATRIX

          pind    = info(2)
          pval    = abs( rinfo(20) )

          lssinfo = 2

      else if ( info(1) .eq. -3 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          if ( iprintctl(3) ) then
              write(* ,9010) fnnzmax,info(17)
              write(file10_unit,9010) fnnzmax,info(17)
          end if
          
          lssinfo = 7
          
      else if ( info(1) .eq. -4 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          if ( iprintctl(3) ) then
              write(* ,9020) wintmax,info(18)
              write(file10_unit,9020) wintmax,info(18)
          end if

          lssinfo = 8

      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) info(1)
              write(file10_unit,9030) info(1)
          end if

          stop

      end if

C     NUMBER OF NEGATIVE EIGENVALUES
      nneigv = info(24)

      if ( lusefac .and. lssinfo .eq. 0 ) then

C        Define matrix D^{-1} (stored in sdiag)
          
          do i = 1,nsys
              w(i) = 0.0d0
          end do
          
          do i = 1,nsys
              
              w(i) = 1.0d0

              call ma57cd(3,nsys,fact,fnnzmax,ifact,wintmax,1,w,nsys,
     +        work,nsysmax,iwork,icntl,info)

              sdiag(i) = w(i)
              w(i)     = 0.0d0
          end do
          
      end if

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSFAC-MA57 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9010 format(/,1X,'LSSFAC-MA57 WARNING: Insufficient double precision ',
     +            'working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9020 format(/,1X,'LSSFAC-MA57 WARNING: Insufficient integer working ',
     +            'space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSFAC-MA57 ERROR: Unhandled error ',I16,'.',
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
#include "ma57dat.com"

      call ma57cd(1,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +nsysmax,iwork,icntl,info)

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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer i

      if ( job .eq. 'T' .or. job .eq. 't' ) then

          call ma57cd(2,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +    nsysmax,iwork,icntl,info)

          do i = 1,nsys
              sol(i) = sol(i) * sqrt( sdiag(i) )
          end do
          
          do i = 1,nsys
              work(i) = sol(i)
          end do

          do i = 1,nsys
              sol(keep(i)) = work(i)
          end do

      else

          do i = 1,nsys
              work(i) = sol(i)
          end do

          do i = 1,nsys
              sol(keep(i)) = work(i)
          end do

          do i = 1,nsys
              sol(i) = sol(i) * sqrt( sdiag(i) )
          end do
          
          call ma57cd(4,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +    nsysmax,iwork,icntl,info)

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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer i

      do i = 1,nsys
          iwork(keep(i)) = v(i)
      end do
      
      do i = 1,nsys
          v(i) = iwork(i)
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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer i

      do i = 1,nsys
          iwork(i) = v(keep(i))
      end do
      
      do i = 1,nsys
          v(i) = iwork(i)
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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer col,i,lin

      do i = 1,hnnz
          col = hcol(i)
          lin = hlin(i)
          
          hcol(i) = keep(col)
          hlin(i) = keep(lin)
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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer col,i,lin

      do i = 1,nsys
          invp(i)=0
      end do
      
      do i = 1,nsys
          invp(keep(i)) = i
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
#include "ma57dat.com"

      lssgetd = sqrt( 1.0d0 / sdiag(invp(j)) )

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssetrow(nsys)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

#include "dim.par"
#include "ma57dat.com"

C     LOCAL SCALARS
      integer apos,i,nblk,iblk,ipiv,ix,iwpos,j1,j2,k,ncols,nrows,totlen,
     +        nsys2

C     Count the number of elements of fact and ifact used to form each
C     row of L^T

      totlen  = 0
      nsys2   = 2*nsys
      
      do i = 1,nsys
          iwork(nsys+i) = 0
      end do

      do i = 1,nsys
          iwork(nsys2+i) = 0
      end do

      do i = 1,nsys
          iwork(i) = 0
      end do

      apos  = 1
      iwpos = 4
      nblk  = ifact(3)
      
      do iblk = 1,nblk
          ncols = ifact(iwpos)
          nrows = ifact(iwpos+1)
          iwpos = iwpos + 2
          j1    = iwpos
          j2    = iwpos + nrows - 1
          
          do ipiv = 1,nrows
              apos = apos + 1
              ix   = abs( ifact(j1) )
              k    = apos
              
              if ( j1+1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k    = max( k, k+j2-j1 )
              apos = k
              j1   = j1 + 1
          end do
          
          j2 = iwpos + ncols - 1
          
          do ipiv = 1,nrows-1,2
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k  = apos+ncols-nrows
              ix = abs( ifact(iwpos+ipiv) )
              
              if ( j1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end do
          
          if ( mod(nrows,2) .eq. 1 ) then
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end if
          
          iwpos = iwpos + ncols
          
      end do
      
      iwork(nsys+1) = 1
      
      do i = 2,nsys
          iwork(nsys+i) = iwork(nsys+i-1) + iwork(nsys2+i-1)
      end do
      
      do i = 1,nsys
          iwork(i) = iwork(nsys+i)
      end do

C     Store the indices of elements of fact and ifact used to form each
C     row of L^T

      apos  = 1
      iwpos = 4
      nblk  = ifact(3)
      
      do iblk = 1,nblk
          ncols = ifact(iwpos)
          nrows = ifact(iwpos+1)
          iwpos = iwpos + 2
          j1    = iwpos
          j2    = iwpos + nrows - 1
          
          do ipiv = 1,nrows
              apos = apos + 1
              ix   = abs( ifact(j1) )
              k    = apos
              
              if ( j1+1 .le. j2 ) then
                  posfac(iwork(ix))   = j1+1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k   = max( k, k+j2-j1 )
              apos = k
              j1   = j1 + 1
          end do
          
          j2 = iwpos + ncols - 1
          
          do ipiv = 1,nrows-1,2
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  posfac(iwork(ix))   = j1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = -k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k  = apos+ncols-nrows
              ix = abs( ifact(iwpos+ipiv) )
              
              if ( j1 .le. j2 ) then
                  posfac(iwork(ix))   = j1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = -k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end do
          
          if ( mod(nrows,2) .eq. 1 ) then
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  posfac(iwork(ix))   = j1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = -k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end if
          
          iwpos = iwpos + ncols
          
      end do
      
C     Array of null elements

      do i = 1,nsys
          iwork(i) = 0
      end do

C     Compute inverse permutation

      do i = 1,nsys
          invp(keep(i)) = i
      end do

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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer i,j,j1,j2,k,ind,ista,iend,ix
      double precision sclf

      ind  = invp(idx)
      ista = iwork(nsys+ind)
      iend = ista + iwork(2*nsys+ind) - 1

      if ( lsclsys ) then
          sclf = s(ind)
      else
          sclf = 1.0d0
      end if

      rownnz         = 1
      rowind(rownnz) = ind
      rowval(rownnz) = sclf
      iwork(ind)     = rownnz

      do i = ista,iend,3
          j1 = posfac(i)
          j2 = posfac(i+1)
          k  = posfac(i+2)
          
          if ( k .gt. 0 ) then
              
              do j = j1,j2
                  ix = abs( ifact(j) )
                  
                  if ( iwork(ix) .eq. 0 ) then
                      if ( fact(k) .ne. 0.0d0 ) then
                          rownnz         = rownnz + 1
                          rowind(rownnz) = ix
                          rowval(rownnz) = sclf * fact(k)
                          iwork(ix)      = rownnz
                      end if
                  else
                      rowval(iwork(ix)) = 
     +                rowval(iwork(ix)) + sclf * fact(k)
                  end if
                  
                  k = k + 1
              end do
          else
              k = abs( k )
              
              do j = j1,j2
                  ix = abs( ifact(j) )
                  
                  if ( iwork(ix) .eq. 0 ) then
                      if ( fact(k) .ne. 0.0d0 ) then
                          rownnz         = rownnz + 1
                          rowind(rownnz) = ix
                          rowval(rownnz) = - sclf * fact(k)
                          iwork(ix)      = rownnz
                      end if
                  else
                      rowval(iwork(ix)) = 
     +                rowval(iwork(ix)) - sclf * fact(k)
                  end if
                  
                  k = k + 1
              end do
          end if
      end do

      if ( lsclsys ) then
          do i = 1,rownnz
              rowval(i) = rowval(i) / s(rowind(i))
          end do
      end if

C     Set used positions of iwork back to 0

      do i = 1,rownnz
          iwork(rowind(i)) = 0
      end do

      do i = 1,rownnz
          rowval(i) = rowval(i) / sqrt( sdiag(ind) )
      end do

      do i = 1,rownnz
          rowind(i) = keep(rowind(i))
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
#include "ma57dat.com"

C     LOCAL SCALARS
      integer i,tmp

C     LOCAL ARRAYS
      integer keeptmp(wintmax),infotmp(40)
      double precision rinfotmp(20)

      tmp       = icntl(15)
      icntl(15) = 0

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          lssinfo = 6

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(file10_unit,9000) nsysmax,nsys
          end if
          
          go to 500
          
      end if

      call ma57ad(nsys,hnnz,hlin,hcol,wintmax,keeptmp,iwork,icntl,
     +infotmp,rinfotmp)

      if ( ( infotmp(1) .eq. 0 ) .or. ( infotmp(1) .eq. 1 ) ) then
        ! SUCCESS

          lssinfo = 0

      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) infotmp(1)
              write(file10_unit,9030) infotmp(1)
          end if
          
          stop
          
      end if
      
      do i = 1,nsys
          hval(hdiag(i)) = hval(hdiag(i)) + d(i)
      end do

      if ( lsclsys ) then
          do i = 1,hnnz
              hval(i) = hval(i) * s(hlin(i)) * s(hcol(i))
          end do
      end if

      call ma57bd(nsys,hnnz,hval,fact,fnnzmax,ifact,wintmax,wintmax,
     +keeptmp,iwork,icntl,cntl,infotmp,rinfotmp)

      if ( lsclsys ) then
          do i = 1,hnnz
              hval(i) = hval(i) * s(hlin(i)) * s(hcol(i))
          end do
      end if

      do i = 1,nsys
          hval(hdiag(i)) = hval(hdiag(i)) - d(i)
      end do

      if ( ( infotmp(1) .eq. 0 ) .or. ( infotmp(1) .eq. 1 ) ) then
        ! SUCCESS

          lssinfo = 0

      else if ( info(1) .eq. 5 .or. info(1) .eq. -6 ) then
        ! MATRIX NOT POSITIVE DEFINITE

          lssinfo = 1

      else if ( infotmp(1) .eq. 4 .or. infotmp(1) .eq. -5 ) then
        ! RANK DEFICIENT MATRIX

          lssinfo = 2
         
      else if ( infotmp(1) .eq. -3 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          lssinfo = 7

          if ( iprintctl(3) ) then
              write(* ,9010) fnnzmax,infotmp(17)
              write(file10_unit,9010) fnnzmax,infotmp(17)
          end if

          go to 500
          
      else if ( infotmp(1) .eq. -4 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          lssinfo = 8

          if ( iprintctl(3) ) then
              write(* ,9020) wintmax,infotmp(18)
              write(file10_unit,9020) wintmax,infotmp(18)
          end if
          
          go to 500
          
      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) infotmp(1)
              write(file10_unit,9030) infotmp(1)
          end if
          
          stop
          
      end if

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if
      
      call ma57cd(1,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +nsysmax,iwork,icntl,infotmp)

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if

 500  continue

      icntl(15) = tmp

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSAFSOL-MA57 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9010 format(/,1X,'LSSAFSOL-MA57 WARNING: Insufficient double ',
     +            'precision working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9020 format(/,1X,'LSSAFSOL-MA57 WARNING: Insufficient integer ',
     +            'working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSAFSOL-MA57 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end
