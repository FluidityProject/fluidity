C     ADD REMOVED FIXED VARIABLES

C     ******************************************************************
C     ******************************************************************

      subroutine uinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"
#include "fixvar.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vinip

      call vinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

C     Eliminate fixed variables (l=u) and save their values on y

      if ( rmfixv ) then

          yind(0) = n

          n = 0
          do i = 1,yind(0)
              if ( l(i) .lt. u(i) ) then
                  n = n + 1
                  yind(n) = i
                  ycor(i) = n
              else
                  y(i) = l(i)
                  ycor(i) = 0
              end if
          end do

          do i = 1,n
              x(i) = x(yind(i))
              l(i) = l(yind(i))
              u(i) = u(yind(i))
          end do

          if ( n .eq. yind(0) ) rmfixv = .false.

          if ( iprintctl(2) ) then
              write(* ,100) yind(0) - n
              write(file10_unit,100) yind(0) - n
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Number of removed fixed variables : ',I7)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(*),lambda(m),u(*),x(*)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer i,ind

C     EXTERNAL SUBROUTINES
      external vendp

C     Restore original x, l, u and n

      if ( rmfixv ) then
          do i = yind(0),1,-1
              ind = ycor(i)
              if ( ind .ne. 0 ) then
                  l(i) = l(ind)
                  u(i) = u(ind)
                  x(i) = x(ind)
              else
                  l(i) = y(i)
                  u(i) = y(i)
                  x(i) = y(i)
              end if
          end do

          n = yind(0)

          rmfixv = .false.
      end if

      call vendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalf

      if ( .not. rmfixv ) then
          call vevalf(n,x,f,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalf: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalf(yind(0),y,f,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

#include "dim.par"
#include "fixvar.com"

C     ARRAY ARGUMENTS
      double precision g(*),x(*)

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalg

      if ( .not. rmfixv ) then
          call vevalg(n,x,g,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalg: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalg(yind(0),y,g,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              g(i) = g(yind(i))
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer col,i,j,lin

C     EXTERNAL SUBROUTINES
      external vevalh

      if ( .not. rmfixv ) then
          call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalh: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalh(yind(0),y,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,hnnz
              lin = ycor(hlin(i))
              col = ycor(hcol(i))
              if ( lin .ne. 0 .and. col .ne. 0 ) then
                  j = j + 1
                  hlin(j) = lin
                  hcol(j) = col
                  hval(j) = hval(i)
              end if
          end do

          hnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalc

      if ( .not. rmfixv ) then
          call vevalc(n,x,ind,c,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalc: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalc(yind(0),y,ind,c,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

#include "dim.par"
#include "fixvar.com"

C     ARRAY ARGUMENTS
      integer jcvar(*)
      double precision x(*),jcval(*)

C     LOCAL SCALARS
      integer i,j,var

C     EXTERNAL SUBROUTINES
      external vevaljac

      if ( .not. rmfixv ) then
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevaljac: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevaljac(yind(0),y,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,jcnnz
              var = ycor(jcvar(i))
              if ( var .ne. 0 ) then
                  j = j + 1
                  jcvar(j) = var
                  jcval(j) = jcval(i)
              end if
          end do

          jcnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer col,i,j,lin

C     EXTERNAL SUBROUTINES
      external vevalhc

      if ( .not. rmfixv ) then
          call vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalhc: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalhc(yind(0),y,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,hnnz
              lin = ycor(hlin(i))
              col = ycor(hcol(i))
              if ( lin .ne. 0 .and. col .ne. 0 ) then
                  j = j + 1
                  hlin(j) = lin
                  hcol(j) = col
                  hval(j) = hval(i)
              end if
          end do

          hnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer col,i,j,lin

C     EXTERNAL SUBROUTINES
      external vevalhl

      if ( .not. rmfixv ) then
          call vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalhl: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalhl(yind(0),y,m,lambda,sf,sc,hlin,hcol,hval,hnnz,
     +    inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,hnnz
              lin = ycor(hlin(i))
              col = ycor(hcol(i))
              if ( lin .ne. 0 .and. col .ne. 0 ) then
                  j = j + 1
                  hlin(j) = lin
                  hcol(j) = col
                  hval(j) = hval(i)
              end if
          end do

          hnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(*),lambda(m),p(n),sc(m),x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision w(nmax)

C     EXTERNAL SUBROUTINES
      external vevalhlp

      if ( .not. rmfixv ) then
          call vevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalhlp: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          do i = 1,yind(0)
              w(i) = 0.0d0
          end do

          do i = 1,n
              w(yind(i)) = p(i)
          end do

          call vevalhlp(yind(0),y,m,lambda,sf,sc,w,hp,gothl,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              hp(i) = hp(yind(i))
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalfc

      if ( .not. rmfixv ) then
          call vevalfc(n,x,f,m,c,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevaljac: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalfc(yind(0),y,f,m,c,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer i,j,var

C     EXTERNAL SUBROUTINES
      external vevalgjac

      if ( .not. rmfixv ) then
          call vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevaljac: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalgjac(yind(0),y,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              g(i) = g(yind(i))
          end do

          j = 0
          do i = 1,jcnnz
              var = ycor(jcvar(i))
              if ( var .ne. 0 ) then
                  j = j + 1
                  jcfun(j) = jcfun(i)
                  jcvar(j) = var
                  jcval(j) = jcval(i)
              end if
          end do

          jcnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine usetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "fixvar.com"

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vsetp

      if ( .not. rmfixv ) then
          call vsetp(n,x)
          return
      end if

      yset = .true.

      do i = 1,n
          y(yind(i)) = x(i)
      end do

      call vsetp(yind(0),y)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uunsetp()

      implicit none

#include "dim.par"
#include "fixvar.com"

C     EXTERNAL SUBROUTINES
      external vunsetp

      yset = .false.

      call vunsetp()

      end
