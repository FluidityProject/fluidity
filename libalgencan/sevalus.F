C     SCALE OBJECTIVE FUNCTION AND CONSTRAINTS

C     ******************************************************************
C     ******************************************************************

      subroutine sinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"
#include "outtyp.com"
#include "scaling.com"
#include "algparam.com"

C     LOCAL SCALARS
      integer i,fun,j,jcnnz,nbds,neq,file50_unit
      double precision scmax
      logical connected

C     LOCAL ARRAYS
      integer jcfun(jcnnzmax),jcvar(jcnnzmax)
      double precision g(nmax),jcval(jcnnzmax)

      neq = 0
      do j = 1,m
          if ( equatn(j) ) neq = neq + 1
      end do

      nbds = 0
      do i = 1,n
          if ( l(i) .gt. - 1.0d+20 ) nbds = nbds + 1
          if ( u(i) .lt.   1.0d+20 ) nbds = nbds + 1
      end do

      if ( iprintctl(2) ) then
          write(* ,100) n,neq,m-neq,nbds
          write(file10_unit,100) n,neq,m-neq,nbds
      end if

      call tinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

C     Write classification line of final model

      if ( iprintctl(6) ) then
          file50_unit_loop: do file50_unit=10, 99
    
            inquire(unit=file50_unit, opened=connected)
    
            if (.not.connected) exit file50_unit_loop
    
          end do file50_unit_loop
          
          open(file50_unit,file='class-tabline.out')
          write(file50_unit,400) n,neq,m-neq,nbds
          close(file50_unit)
      end if

C     Scaling

      usf = 1.0d0
      do j = 1,m
          usc(j) = 1.0d0
      end do

      if ( scale ) then

          if ( m .eq. 0 ) then
              sf = 1.0d0
              if ( iprintctl(2) ) then
                  write(* ,200) 1.0d0 / sf
                  write(file10_unit,200) 1.0d0 / sf
              end if

              return
          end if

          call tsetp(n,x)

          if ( gjaccoded ) then

              call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
              if ( inform .lt. 0 ) return

C             Scale constraints

              do j = 1,m
                  sc(j) = 1.0d0
              end do

              do i = 1,jcnnz
                  fun = jcfun(i)

                  sc(fun) = max( sc(fun), abs( jcval(i) ) )
              end do

          else

              call tevalg(n,x,g,inform)
              if ( inform .lt. 0 ) return

C             Scale constraints

              do j = 1,m
                  sc(j) = 1.0d0

                  call tevaljac(n,x,j,jcvar,jcval,jcnnz,inform)
                  if ( inform .lt. 0 ) return

                  do i = 1,jcnnz
                      sc(j) = max( sc(j), abs( jcval(i) ) )
                  end do
              end do

          end if

C         Scale objective function

          sf = 1.0d0
          do i = 1,n
              sf = max( sf, abs( g(i) ) )
          end do

C         Report scaling factors

          scmax = 0.0d0
          do j = 1,m
              scmax = max( scmax, sc(j) )
          end do

          if ( iprintctl(2) ) then
              write(* ,300) 1.0d0 / sf,1.0d0 / scmax
              write(file10_unit,300) 1.0d0 / sf,1.0d0 / scmax
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Number of variables               : ',I7,
     +       /,1X,'Number of equality constraints    : ',I7,
     +       /,1X,'Number of inequality constraints  : ',I7,
     +       /,1X,'Number of bound constraints       : ',I7)

 200  format(/,1X,'Objective function scale factor   : ',1P,D7.1,
     +       /,1X,'The scaling feature was mainly developed for ',
     +            'constrained problems. For',/,1X,'unconstrained and ',
     +            'bound-constrained problem, please, set the ',
     +            'optimality',/,1X,'tolerance (related to the ',
     +            'sup-norm of the projected gradient of the',/,1X,
     +            'objective function) with a convenient value.')

 300  format(/,1X,'Objective function scale factor   : ',1P,D7.1,
     +       /,1X,'Smallest constraints scale factor : ',1P,D7.1)

 400  format(  1X,I6,1X,I6,1X,I6,1X,I6)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"
#include "scaling.com"

C     LOCAL SCALARS
      integer i

      if ( scale ) then
          do i = 1,m
              lambda(i) = lambda(i) * sf / sc(i)
          end do

          scale = .false.
      end if

      call tendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "scaling.com"

      call tevalf(n,x,f,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) f = f / sf

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

#include "dim.par"
#include "scaling.com"

C     LOCAL SCALARS
      integer i

      call tevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,n
              g(i) = g(i) / sf
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "scaling.com"

C     LOCAL SCALARS
      integer i

      call tevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,hnnz
              hval(i) = hval(i) / sf
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "scaling.com"

      call tevalc(n,x,ind,c,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) c = c / sc(ind)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

#include "dim.par"
#include "scaling.com"

C     LOCAL SCALARS
      integer i

      call tevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,jcnnz
              jcval(i) = jcval(i) / sc(ind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "scaling.com"

C     LOCAL SCALARS
      integer i

      call tevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,hnnz
              hval(i) = hval(i) / sc(ind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhl(n,x,m,lambda,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),x(n)

#include "dim.par"
#include "scaling.com"

      if ( scale ) then
          call tevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call tevalhl(n,x,m,lambda,usf,usc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhlp(n,x,m,lambda,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),x(n)

#include "dim.par"
#include "scaling.com"

      if ( scale ) then
          call tevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

      else
          call tevalhlp(n,x,m,lambda,usf,usc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

#include "dim.par"
#include "scaling.com"

C     LOCAL SCALARS
      integer j

      call tevalfc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          f = f / sf

          do j = 1,m
              c(j) = c(j) / sc(j)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

#include "dim.par"
#include "scaling.com"

C     LOCAL SCALARS
      integer i

      call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,n
              g(i) = g(i) / sf
          end do

          do i = 1,jcnnz
              jcval(i) = jcval(i) / sc(jcfun(i))
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ssetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

      call tsetp(n,x)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sunsetp()

      implicit none

      call tunsetp()

      end
