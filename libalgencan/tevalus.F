C     ADD SLACKS

C     ******************************************************************
C     ******************************************************************

      subroutine tinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"
#include "slacks.com"
#include "outtyp.com"
#include "algparam.com"

C     LOCAL SCALARS
      integer j
      double precision dum

C     LOCAL ARRAYS
      double precision c(mmax)

      call uinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

      if ( slacks ) then

          nws = n

          call usetp(nws,x)

          if ( fccoded ) then
              call uevalfc(nws,x,dum,m,c,inform)
              if ( inform .lt. 0 ) return

          else
              do j = 1,m
                  if ( .not. equatn(j) ) then
                      call uevalc(nws,x,j,c(j),inform)
                      if ( inform .lt. 0 ) return
                  end if
              end do
          end if

          do j = 1,m
              if ( equatn(j) ) then
                  slaind(j) = - 1
              else
                  equatn(j) = .true.

                  n         = n + 1
                  slaind(j) = n

                  l(n)      = - 1.0d+20
                  u(n)      =   0.0d0
                  x(n)      = max( l(n), min( c(j), u(n) ) )
              end if
          end do

          if ( n .eq. nws ) slacks = .false.

          if ( iprintctl(2) ) then
              write(* ,100) n - nws
              write(file10_unit,100) n - nws
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Number of added slack variables   : ',I7)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"
#include "slacks.com"

C     LOCAL SCALARS
      integer j

      if  ( slacks ) then
          n = nws

          do j = 1,m
              if ( slaind(j) .ne. - 1 ) then
                  equatn(j) = .false.
              end if
          end do

          slacks = .false.
      end if

      call uendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "slacks.com"

      if ( .not. slacks ) then
          call uevalf(n,x,f,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalf(nws,x,f,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

#include "dim.par"
#include "slacks.com"

C     LOCAL SCALARS
      integer i

      if ( .not. slacks ) then
          call uevalg(n,x,g,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalg(nws,x,g,inform)
          if ( inform .lt. 0 ) return

          do i = nws + 1,n
              g(i) = 0.0d0
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "slacks.com"

      if ( .not. slacks ) then
          call uevalh(n,x,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalh(nws,x,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "slacks.com"

C     LOCAL SCALARS
      integer sind

      if ( .not. slacks ) then
          call uevalc(n,x,ind,c,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalc(nws,x,ind,c,inform)
          if ( inform .lt. 0 ) return

          sind = slaind(ind)
          if ( sind .ne. - 1 ) c = c - x(sind)
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

#include "dim.par"
#include "slacks.com"

C     LOCAL SCALARS
      integer sind

      if ( .not. slacks ) then
          call uevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevaljac(nws,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          sind = slaind(ind)
          if ( sind .ne. - 1 ) then
              jcnnz = jcnnz +  1
              jcvar(jcnnz) = sind
              jcval(jcnnz) = - 1.0d0
          end if
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "slacks.com"

      if ( .not. slacks ) then
          call uevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalhc(nws,x,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

#include "dim.par"
#include "slacks.com"

      if ( .not. slacks ) then
          call uevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalhl(nws,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

#include "dim.par"
#include "slacks.com"

C     LOCAL SCALARS
      integer i

      if ( .not. slacks ) then
          call uevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalhlp(nws,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

          do i = nws + 1,n
              hp(i) = 0.0d0
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

#include "dim.par"
#include "slacks.com"

C     LOCAL SCALARS
      integer j,sind

      if ( .not. slacks ) then
          call uevalfc(n,x,f,m,c,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalfc(nws,x,f,m,c,inform)
          if ( inform .lt. 0 ) return

          do j = 1,m
              sind = slaind(j)
              if ( sind .ne. - 1 ) c(j) = c(j) - x(sind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

#include "dim.par"
#include "slacks.com"

C     LOCAL SCALARS
      integer i,j,sind

      if ( .not. slacks ) then
          call uevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalgjac(nws,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = nws + 1,n
              g(i) = 0.0d0
          end do

          do j = 1,m
              sind = slaind(j)
              if ( sind .ne. - 1 ) then
                  jcnnz = jcnnz +  1
                  jcfun(jcnnz) = j
                  jcvar(jcnnz) = sind
                  jcval(jcnnz) = - 1.0d0
              end if
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tsetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "slacks.com"

      if ( .not. slacks ) then
          call usetp(n,x)

      else
          call usetp(nws,x)
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tunsetp()

      implicit none

      call uunsetp()

      end
