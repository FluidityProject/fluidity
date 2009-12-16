C     ******************************************************************
C     ******************************************************************

      subroutine algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,
     +equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,iprint,m,n,ncomp
      double precision cnorm,epsfeas,epsopt,f,nlpsupn,snorm

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"
#include "machconst.com"
#include "algconst.par"
#include "counters.com"
#include "outtyp.com"
#include "algparam.com"
#include "scaling.com"
#include "slacks.com"
#include "fixvar.com"

C     LOCAL SCALARS
      integer alinfo,geninfo,i,iter,j,maxit,nwcalls,nwtotit,outiter,
     +        solinfo,totiter,msqcalls,msqtotit,file20_unit
      double precision cnormb,cnormu,cnormub,fb,fu,fub,nlpsupnb
      real time
      logical :: connected
 
C     LOCAL ARRAYS
      double precision nl(nmax),c(mmax),rho(mmax)
      real(kind=4) :: dum(2)

C     DATA STATEMENTS
      data dum/0.0,0.0/

C     EXTERNAL FUNCTIONS AND SUBROUTINES
      external fparam,checkd,auglag,gencan

C     ==================================================================
C     Start timing
C     ==================================================================

      time = dtime(dum)

C     ==================================================================
C     Initialization
C     ==================================================================

C     Set machine-dependent constants

      bignum    = 1.0d+99
      macheps   = 1.0d-16
      macheps12 = sqrt( macheps )
      macheps13 = macheps ** ( 1.0d0 / 3.0d0 )
      macheps23 = macheps ** ( 2.0d0 / 3.0d0 )

C     Set global counters

      fcnt    = 0
      efcnt   = 0
      efccnt  = 0
      egcnt   = 0
      egjccnt = 0
      ehcnt   = 0
      ehlcnt  = 0
      ehlpcnt = 0

      do j = 1,m
          eccnt(j)  = 0
          ejccnt(j) = 0
          ehccnt(j) = 0
      end do

C     ==================================================================
C     Set default values for algoritmic parameters
C     ==================================================================

C     Set user-provided subroutines indicators

      fcoded    = coded(1)
      gcoded    = coded(2)
      hcoded    = coded(3)
      ccoded    = coded(4)
      jaccoded  = coded(5)
      hccoded   = coded(6)
      fccoded   = coded(7)
      gjaccoded = coded(8)
      hlcoded   = coded(9)
      hlpcoded  = coded(10)

      innercall = .false.
      useustp   = .false.

C     Set indicator of whether the true Hessian of the Lagrangian can be
C     computed or not

      truehl = .false.
      if ( hlcoded .or. ( hcoded .and. ( hccoded .or. m .eq. 0 ) )) then
          truehl = .true.
      end if

C     Hessian-vector product strategy: HAPPRO, INCQUO or TRUEHL (TRUEHL
C     is the default option. If the proper subroutines were not coded by
C     the user, then HAPPRO is used instead.)

      if ( truehl .or. hlpcoded ) then
          hptype = 'TRUEHL'
      else
          hptype = 'HAPPRO'
      end if

C     Ignore objective function (to only find a feasible point by 
C     minimizing 1/2 of the squared infeasibility)
      ignoref = .false.

C     Inner-to-the-face minimization algorithm (CG is the default option)
      avoidds = .true.

C     Skip acceleration step
      skipacc = .true.

C     Scaling of linear systems
      sclsys = .false.

C     Slacks for inequality constraints
      slacks = .false.

C     Remove fixed variables (with identical lower and upper bounds)
      rmfixv = .true.

C     Scale objective function and constraints
      if ( m .gt. 0 ) then
          scale = .true.
      else
          scale = .false.
      end if

C     Main output control (silent-mode?)

      iprintctl(1) = .true. ! Banner
      iprintctl(2) = .true. ! Parameters and problem processing
      iprintctl(3) = .true. ! Warnings and errors messages
      iprintctl(4) = .true. ! Screen-mirror file algencan.out
      iprintctl(5) = .true. ! Solution file solution.txt
      iprintctl(6) = .true. ! Statistic files with table lines

      file10_unit_loop: do file10_unit=10, 99

        inquire(unit=file10_unit, opened=connected)

        if (.not.connected) exit file10_unit_loop

      end do file10_unit_loop
      
      open(file10_unit,err=100,file='.silent',status='old')
      close(file10_unit)

      do i = 1,6
          iprintctl(i) = .false.
      end do

      iprint = 0

 100  continue

      if ( iprintctl(4) ) then
          open(unit=file10_unit,file='algencan.out',status='replace')
      else
          open(unit=file10_unit,                    status='scratch')
      end if

C     ==================================================================
C     Set solver arguments using the specification file
C     ==================================================================

      call fparam(epsfeas,epsopt,iprint,ncomp)

C     Outer and inner iterations output detail

      iprintout = iprint / 10
      iprintinn = mod( iprint, 10 )

C     Error tracker

      inform = 0

C     ==================================================================
C     Initialize problem data structures
C     ==================================================================

      call sinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

      nprint = min( n, ncomp )
      mprint = min( m, ncomp )

C     ==================================================================
C     Call the solver
C     ==================================================================

C     ALGENCAN for PNL problems

      if ( .not. ignoref .and. m .gt. 0 ) then
          call auglag(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,f,c,
     +    cnorm,snorm,nl,nlpsupn,fu,cnormu,fub,cnormub,fb,cnormb,
     +    nlpsupnb,outiter,totiter,nwcalls,nwtotit,msqcalls,msqtotit,
     +    alinfo,inform)

          solinfo = alinfo

C     GENCAN for box-constrained problems

      else
          maxit = 1000

C         Used in feasibility problems (ignoref=true). With lambda=0 and
C         rho=1, to minimize 1/2 of the squared infeasibility coincides
C         with minimizing the augmented Lagrangian.
          do j = 1,m
              lambda(j) = 0.0d0
              rho(j)    = 1.0d0
          end do

          call gencan(n,x,l,u,m,lambda,equatn,linear,rho,epsfeas,epsopt,
     +    maxit,iter,f,nl,nlpsupn,cnorm,cnormu,geninfo,inform)

          solinfo  = geninfo

          outiter  = 0
          totiter  = iter
          nwcalls  = 0
          nwtotit  = 0
          msqcalls = 0
          msqtotit = 0

          if ( scale ) then
              fu = f * sf
          else
              fu = f
          end if

          fb       = f
          fub      = fu
          cnormb   = cnorm
          cnormub  = cnormu
          nlpsupnb = nlpsupn
      end if

      if ( inform .lt. 0 ) return

C     Close output file

      close(file10_unit)

C     ==================================================================
C     End problem data structures
C     ==================================================================

      call sendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     Stop timing
C     ==================================================================

      time = dtime(dum)
      time = dum(1)

      if ( iprintout .eq. 1 ) then
          write(*,9000) time
      end if

C     ==================================================================
C     Write statistics
C     ==================================================================

      if ( iprintctl(6) ) then
          file20_unit_loop: do file20_unit=10, 99
    
            inquire(unit=file20_unit, opened=connected)
    
            if (.not.connected) exit file20_unit_loop
    
          end do file20_unit_loop
          
          open(file20_unit,file='algencan-tabline.out')
          write(file20_unit,9040) fu,cnormu,f,cnorm,nlpsupn,fub,cnormub,fb,
     +                   cnormb,nlpsupnb,inform,solinfo,n,m,outiter,
     +                   totiter,fcnt,nwcalls,nwtotit,msqcalls,msqtotit,
     +                   time
          close(file20_unit)
      end if

C     ==================================================================
C     NON-EXECUTABLE STATEMENTS
C     ==================================================================

 9000 format(/,1X,'Total CPU time in seconds: ',F8.2)
 9040 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,I3,1X,I1,1X,I6,1X,I6,1X,I2,1X,I7,1X,I7,1X,I2,1X,I7,
     +       1X,I7,1X,I7,0P,F8.2)


      end
