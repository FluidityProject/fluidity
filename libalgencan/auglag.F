C     ******************************************************************
C     ******************************************************************

      subroutine auglag(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,
     +f,c,cnorm,snorm,nl,nlpsupn,fu,cnormu,fub,cnormub,fb,cnormb,
     +nlpsupnb,outiter,totiter,nwcalls,nwtotit,msqcalls,msqtotit,alinfo,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      integer alinfo,inform,m,msqcalls,msqtotit,n,nwcalls,nwtotit,
     +        outiter,totiter
      double precision cnorm,cnormb,cnormu,cnormub,epsfeas,epsopt,f,fb,
     +        fu,fub,nlpsupn,nlpsupnb,snorm

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision c(m),l(n),lambda(m),nl(n),rho(m),u(n),x(n)

C     Solves the nonlinear programming problem
C
C     min f(x)
C
C     subject to
C
C             c_j(x)  = 0, j in E,
C             c_j(x) <= 0, j in I,
C             l <= x <= u,
C
C     where E is the set of indices of the equality constraints, I is
C     the set of indices of the inequality constraints, and there are
C     n variables and m constraints, using the method of multipliers
C     described in
C
C     R. Andreani, E. G. Birgin, J. M. Mart�nez and M. L. Schuverdt,
C     "On Augmented Lagrangian methods with general lower-level
C     constraints", SIAM Journal on Optimization 18, pp. 1286-1309, 2007.

C     alinfo:
C
C     0: Feasibility, optimality and complementarity satisfied
C     1: Maximum number of algencan iterations reached
C     2: Infeasible problem?

#include "dim.par"
#include "machconst.com"
#include "algconst.par"
#include "algparam.com"
#include "outtyp.com"
#include "scaling.com"

C     LOCAL SCALARS
      character innstp
      logical scaletmp,connected
      integer geninfo,i,iter,maxit,msqiter,nwinfo,nwiter,outitnp,
     +        rhoincr,xind,file20_unit
      double precision al,cnorma,epsfeas12,epsfeas14,epsopt12,epsopt14,
     +        epsopk,fa,nlpi,nlpsupna,rsupn,snormb,snormprev,xsupn

C     LOCAL ARRAYS
      double precision lambar(mmax),lambdaa(mmax),lambdab(mmax),
     +        sigma(mmax),sigpre(mmax),xa(nmax),xb(nmax)

C     DATA BLOCKS
      character genstp(0:9)
      data genstp(0) /'C'/
      data genstp(1) /'M'/
      data genstp(2) /' '/
      data genstp(3) /'F'/
      data genstp(4) /'G'/
      data genstp(5) /'P'/
      data genstp(6) /'U'/
      data genstp(7) /'S'/
      data genstp(8) /' '/
      data genstp(9) /' '/

C     ==================================================================
C     Print initial information
C     ==================================================================

      if ( iprintout .ge. 1 ) then
          write(* ,1000) n,m
          write(file10_unit,1000) n,m

          if ( iprintout .ge. 4 .and. nprint .ne. 0 ) then
              write(* ,1010) nprint,(l(i),i=1,nprint)
              write(* ,1020) nprint,(u(i),i=1,nprint)

              write(file10_unit,1010) nprint,(l(i),i=1,nprint)
              write(file10_unit,1020) nprint,(u(i),i=1,nprint)
          end if

      end if

C     ==================================================================
C     Initialization
C     ==================================================================

C     Counters

      outiter = 0
      outitnp = 0

      totiter = 0

      nwcalls = 0
      nwtotit = 0

      msqcalls = 0
      msqtotit = 0

C     Constants

      epsopt12  = sqrt( epsopt    )
      epsopt14  = sqrt( epsopt12  )
      epsfeas12 = sqrt( epsfeas   )
      epsfeas14 = sqrt( epsfeas12 )

C     Project initial point

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Compute current point norm

      xsupn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind  = i
          end if
      end do

C     Compute objective function and constraints

      call sevalobjc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

C     Set penalty parameters and compute maximum rho

      if ( rhoauto ) then
          call comprhoini(f,m,c,equatn,rho)
      end if

      rsupn = 0.0d0
      do i = 1,m
          rsupn = max( rsupn, rho(i) )
      end do

C     Compute safeguarded Lagrange multipliers

      do i = 1,m
          lambar(i) = max( lammin, min( lambda(i), lammax ) )
      end do

C     Compute complementarity and feasibility violations

      snormprev = bignum

      cnorm  = 0.0d0
      snorm  = 0.0d0
      do i = 1,m
          sigpre(i) = bignum

          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( c(i) ) )
              sigma(i) = c(i)

          else
              cnorm = max( cnorm, c(i) )
              sigma(i) = max( c(i), - lambar(i) / rho(i) )
          end if

          snorm = max( snorm, abs( sigma(i) ) )
      end do

C     Compute unscaled objective function and constraints norm

      if ( scale ) then
          fu = f * sf

          cnormu = 0.0d0
          do i = 1,m
              if ( equatn(i) ) then
                  cnormu = max( cnormu, abs( sc(i) * c(i) ) )
              else
                  cnormu = max( cnormu, sc(i) * c(i) )
              end if
          end do

      else
          fu     = f
          cnormu = cnorm
      end if

C     Compute continuous projected Lagrangian gradient norm

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

      nlpsupn = 0.0d0
      do i = 1,n
          nlpi = x(i) - nl(i)
          if ( l(i) .le. nlpi .and. nlpi .le. u(i) ) then
              nlpi = - nl(i)
          else
              nlpi = max( l(i), min( nlpi, u(i) ) ) - x(i)
          end if
          nlpsupn = max( nlpsupn, abs( nlpi ) )
      end do

C     Save best solution

      fub      = fu
      cnormub  = cnormu

      fb       = f
      cnormb   = cnorm
      snormb   = snorm
      nlpsupnb = nlpsupn

      do i = 1,n
          xb(i) = x(i)
      end do

      do i = 1,m
          lambdab(i) = lambda(i)
      end do

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Print information of this iteration
C     ==================================================================


      if ( iprintout .eq. 1 ) then

          if ( outiter .gt. 0 ) then
              innstp = genstp(geninfo)
          else
              innstp = ' '
          end if

          if ( .not. scale ) then
              if ( mod(outiter,10) .eq. 0 ) then
                  write(* ,1030)
                  write(file10_unit,1030)
              end if
              write(* ,1040) outiter,f,cnorm,snorm,nlpsupn,xsupn,rsupn,
     +                       totiter,innstp,nwcalls,nwtotit
              write(file10_unit,1040) outiter,f,cnorm,snorm,nlpsupn,xsupn,rsupn,
     +                       totiter,innstp,nwcalls,nwtotit
          else
              if ( mod(outiter,10) .eq. 0 ) then
                  write(* ,1050)
                  write(file10_unit,1050)
              end if
              write(* ,1060) outiter,fu,cnormu,f,cnorm,snorm,nlpsupn,
     +                       rsupn,totiter,innstp,nwcalls,nwtotit
              write(file10_unit,1060) outiter,fu,cnormu,f,cnorm,snorm,nlpsupn,
     +                       rsupn,totiter,innstp,nwcalls,nwtotit
          end if

      else if ( iprintout .ge. 2 ) then
          write(* ,1070) outiter
          write(file10_unit,1070) outiter

          if ( iprintout .ge. 3 ) then
              write(* ,1080) totiter,nwcalls,nwtotit
              write(file10_unit,1080) totiter,nwcalls,nwtotit
          end if

          if ( .not. scale ) then
              write(* ,1090) f,cnorm
              write(file10_unit,1090) f,cnorm
          else
              write(* ,1100) f,fu,cnorm,cnormu
              write(file10_unit,1100) f,fu,cnorm,cnormu
          end if

          write(* ,1110) snorm,nlpsupn,xind,xsupn,rsupn
          write(file10_unit,1110) snorm,nlpsupn,xind,xsupn,rsupn

          if ( iprintout .ge. 4 ) then
              if ( nprint .ne. 0 ) then
                  write(* ,1120) nprint,(x(i),i=1,nprint)
                  write(file10_unit,1120) nprint,(x(i),i=1,nprint)
              end if

              if ( mprint .ne. 0 ) then
                  write(* ,1130) mprint,(lambda(i),i=1,mprint)
                  write(* ,1140) mprint,(rho(i),i=1,mprint)

                  write(file10_unit,1130) mprint,(lambda(i),i=1,mprint)
                  write(file10_unit,1140) mprint,(rho(i),i=1,mprint)
              end if
          end if

      end if

C     Save intermediate data for crash report

      if ( iprintctl(6) ) then
          file20_unit_loop: do file20_unit=10, 99
    
            inquire(unit=file20_unit, opened=connected)
    
            if (.not.connected) exit file20_unit_loop
    
          end do file20_unit_loop
          
          open(file20_unit,file='algencan-tabline.out')
          write(file20_unit,1400) fu,cnormu,f,cnorm,nlpsupn,fub,cnormub,fb,
     +                   cnormb,nlpsupnb,inform,9,n,m,outiter,totiter,0,
     +                   nwcalls,nwtotit,msqcalls,msqtotit,999.9d0
          close(file20_unit)
      end if

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     Test feasibility, optimality and complementarity

      if ( max( snorm, cnormu ) .le. epsfeas .and.
     +     nlpsupn .le. epsopt ) then
          alinfo = 0

          if ( iprintout .ge. 1 ) then
              write(*, 1300)
              write(file10_unit,1300)
          end if

          return
      end if

C     Test whether the number of iterations is exhausted

      if ( outiter .ge. maxoutit ) then
          alinfo = 1

          if ( iprintout .ge. 1 ) then
              write(*, 1310)
              write(file10_unit,1310)
          end if

          return
      end if

C     Test whether the problem seems to be infeasible

      if ( snorm .gt. epsfeas .and. snorm .ge. snormprev ) then
          outitnp = outitnp + 1
      else
          outitnp = 0
      end if

      if ( outitnp .ge. maxoutitnp .or. rsupn .gt. rhomax ) then

          alinfo = 2

          if ( iprintout .ge. 1 ) then
              write(*, 1320)
              write(file10_unit,1320)
          end if

          return
      end if

C     ==================================================================
C     Near the solution, try to solve the KKT system by Newton's method
C     ==================================================================

      if ( .not. skipacc .and. truehl .and. ( nwcalls .gt. 0 .or.
     +     ( snorm .le. epsfeas12 .and. nlpsupn .le. epsopt12 ) .or.
     +     ( snorm .le. epsfeas14 .and. nlpsupn .le. epsopt14 .and.
     +       outiter .gt. 0 .and. geninfo .ne. 0 ) ) ) then

          nwcalls = nwcalls + 1

          do i = 1,n
              xa(i) = x(i)
          end do

          do i = 1,m
              lambdaa(i) = lambar(i)
          end do

          scaletmp = scale

          if ( scale ) then
              do i = 1,m
                  lambdaa(i) = lambdaa(i) * sf / sc(i)
              end do

              scale = .false.
          end if

          call newtonkkt(n,xa,l,u,m,lambdaa,equatn,linear,epsfeas,
     +    epsopt,fa,cnorma,nlpsupna,nwiter,msqiter,nwinfo,inform)

          nwtotit = nwtotit  + nwiter

          if ( nwinfo .ge. 6 ) then
              skipacc = .true.
          end if

          if ( msqiter .gt. 0 ) then
              msqcalls = msqcalls + 1
              msqtotit = msqtotit + msqiter
          end if

          scale = scaletmp

          if ( inform .lt. 0 ) return

C         Save best solution
 
          if ( ( cnormub .gt. epsfeas .and. cnorma .le. cnormub ) .or. 
     +         ( cnormub .le. epsfeas .and. cnorma .le. epsfeas .and. 
     +           fa .le. fub ) ) then

              fub      = fa
              cnormub  = cnorma

              fb       = fa
              cnormb   = cnorma
              snormb   = 0.0d0
              nlpsupnb = nlpsupna

              do i = 1,n
                  xb(i) = xa(i)
              end do

              do i = 1,m
                  lambdab(i) = lambdaa(i)
              end do

          end if

C         Test feasibility, optimality and complementarity

          if ( cnorma .le. epsfeas .and. nlpsupna .le. epsopt ) then

              fu      = fa
              cnormu  = cnorma

              f       = fa
              cnorm   = cnorma
              snorm   = 0.0d0
              nlpsupn = nlpsupna

              do i = 1,n
                  x(i) = xa(i)
              end do

              do i = 1,m
                  lambda(i) = lambdaa(i)
              end do

              alinfo = 0

              if ( iprintout .ge. 1 ) then
                  write(*, 1300)
                  write(file10_unit,1300)
              end if

              return

          end if

      end if

C     ==================================================================
C     Iteration
C     ==================================================================

      outiter = outiter + 1

C     ==================================================================
C     Solve the augmented Lagrangian subproblem
C     ==================================================================

C     Set optimality requeriment for the subproblem

      if ( outiter .gt. 1 .and. 
     +     snorm .le. epsfeas12 .and. nlpsupn .le. epsopt12 ) then
          epsopk = max( epsopt, 0.1d0 * epsopk )
      else
          epsopk = sqrt( epsopt )
      end if

      if ( outiter .eq. 1 ) then
          maxit =   10
      else
          maxit = 1000
      end if

C     Call the inner-solver

      call gencan(n,x,l,u,m,lambar,equatn,linear,rho,epsfeas,epsopk,
     +maxit,iter,al,nl,nlpsupn,cnorm,cnormu,geninfo,inform)

      totiter = totiter + iter

      if ( inform .lt. 0 ) return

C     ==================================================================
C     Prepare for the next iteration
C     ==================================================================

C     Compute current point norm

      xsupn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind  = i
          end if
      end do

C     Compute objective function and constraints

      call sevalobjc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

C     Compute complementarity and feasibility violations

      snormprev = snorm

      cnorm  = 0.0d0
      snorm  = 0.0d0
      do i = 1,m
          sigpre(i) = sigma(i)

          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( c(i) ) )
              sigma(i) = c(i)

          else
              cnorm = max( cnorm, c(i) )
              sigma(i) = max( c(i), - lambar(i) / rho(i) )
          end if

          snorm = max( snorm, abs( sigma(i) ) )
      end do

C     Compute unscaled objective function and constraints norm

      if ( scale ) then
          fu = f * sf

          cnormu = 0.0d0
          do i = 1,m
              if ( equatn(i) ) then
                  cnormu = max( cnormu, abs( sc(i) * c(i) ) )
              else
                  cnormu = max( cnormu, sc(i) * c(i) )
              end if
          end do

      else
          fu     = f
          cnormu = cnorm
      end if

C     Update Lagrange multipliers approximation

      do i = 1,m
          call evaldpdy(c(i),rho(i),lambar(i),equatn(i),lambda(i))
          lambar(i) = max( lammin, min( lambda(i), lammax ) )
      end do

C     Compute continuous projected Lagrangian gradient norm

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

      nlpsupn = 0.0d0
      do i = 1,n
          nlpi = x(i) - nl(i)
          if ( l(i) .le. nlpi .and. nlpi .le. u(i) ) then
              nlpi = - nl(i)
          else
              nlpi = max( l(i), min( nlpi, u(i) ) ) - x(i)
          end if
          nlpsupn = max( nlpsupn, abs( nlpi ) )
      end do

C     Save best solution

      if ( ( cnormub .gt. epsfeas .and. cnormu .le. cnormub ) .or. 
     +     ( cnormub .le. epsfeas .and. cnormu .le. epsfeas .and. 
     +       f .lt. fb ) ) then

          fub      = fu
          cnormub  = cnormu

          fb       = f
          cnormb   = cnorm
          snormb   = snorm
          nlpsupnb = nlpsupn

          do i = 1,n
              xb(i) = x(i)
          end do

          do i = 1,m
              lambdab(i) = lambda(i)
          end do

      end if

C     Update penalty parameters

      if ( outiter .eq. 1 ) then

          call comprhoini(f,m,c,equatn,rho)

          if ( iprintout .ge. 5 ) then
              write(*, 1200)
              write(file10_unit,1200)
          end if

      else

C         Note that |sigma(i)| <= eps implies c(i) <= eps.

          if ( rhoiden ) then

              if ( max( snorm, cnormu ) .gt. epsfeas .and.
     +             snorm .gt. rhofrac * snormprev ) then

                  do i = 1,m
                      rho(i) = rhomult * rho(i)
                  end do

                  if ( iprintout .ge. 5 ) then
                      write(*, 1210) rhomult
                      write(file10_unit,1210) rhomult
                  end if

              else

                  if ( iprintout .ge. 5 ) then
                      write(*, 1230)
                      write(file10_unit,1230)
                  end if

              end if

          else

              rhoincr = 0
              do i = 1,m
                  if ( max(abs(sigma(i)), sc(i)*c(i)) .gt. epsfeas .and.
     +                 abs(sigma(i)) .gt. rhofrac*abs(sigpre(i)) ) then
                      rho(i)  = rhomult * rho(i)
                      rhoincr = rhoincr + 1
                  end if
              end do

              if ( iprintout .ge. 5 ) then

                  if ( rhoincr .ne. 0 ) then
                      write(*, 1220) rhoincr,m,rhomult
                      write(file10_unit,1220) rhoincr,m,rhomult
                  else
                      write(*, 1230)
                      write(file10_unit,1230)
                  end if

              end if

          end if

      end if

C     Compute maximum rho

      rsupn = 0.0d0
      do i = 1,m
          rsupn = max( rsupn, rho(i) )
      end do

C     ==================================================================
C     Iterate
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 1000 format(/,1X,'Entry to ALGENCAN.',
     +       /,1X,'Number of variables  : ',I7,
     +       /,1X,'Number of constraints: ',I7)
 1010 format(/,1X,'Lower bounds  (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1020 format(/,1X,'Upper bounds  (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))

 1030 format(/,1X,'It',15X,'f',4X,'cnrm',4X,'snrm',2X,'nlpnrm',4X,
     +         'xnrm',1X,'penalty',2X,'innit',1X,'nwcalls',2X,'nwit')
 1040 format(I3,1P,D16.8,5(1X,1P,D7.1),1X,I5,A1,2(1X,I6))

 1050 format(/,1X,'It',4X,'unscaled f & cnrm',6X,'scaled f & cnrm',
     +         3X,'snrm',1X,'nlpnrm',1X,'penalty',1X,'innit',2X,'nwkkt')
 1060 format(I3,2(1P,D14.6,1X,1P,D6.0),3(1X,1P,D6.0),1X,I5,A1,1X,I2,1X,
     +         I3)

 1070 format(/,1X,'ALGENCAN OUTER ITERATION                      = ',I7)

 1080 format(/,1X,'Up-to-now total number of iterations          = ',I7,
     +       /,1X,'Up-to-now acceleration trials                 = ',I7,
     +       /,1X,'Up-to-now total number of Newton iterations   = ',I7)

 1090 format(/,1X,'Functional value                              = ',
     +             1P,D24.16,
     +       /,1X,'Sup-norm of constraints                       = ',
     +             17X,1P,D7.1)

 1100 format(/,1X,'Functional value          (scaled = ',1P,D8.1,
     +            ') = ',1P,D24.16,
     +       /,1X,'Sup-norm of constraints   (scaled = ',1P,D8.1,
     +            ') = ',17X,1P,D7.1)

 1110 format(  1X,'Sup-norm of complementarity-feasibility       = ',
     +             17X,1P,D7.1,
     +       /,1X,'Sup-norm of the Lagrangian projected gradient = ',
     +             17X,1P,D7.1,
     +       /,1X,'Sup-norm of x (attained at x_{', I7,'})       = ',
     +             17X,1P,D7.1,
     +       /,1X,'Largest penalty parameter                     = ',
     +             17X,1P,D7.1)

 1120 format(/,1X,'Current point (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1130 format(/,1X,'Updated Lagrange multipliers (first ',I7,
     +         1X,'components):',
     +       /,6(1X,1P,D11.4))
 1140 format(/,1X,'Updated penalty parameters (first ',I7,
     +         1X,'components):',
     +       /,6(1X,1P,D11.4))

 1200 format(/,1X,'Penalty parameters are re-initiated after the ',
     +            'resolution of the first',/,1X,'subproblem.')
 1210 format(/,1X,'The desired infeasibility was not achieved. ',
     +            'The penalty',/,1X,'parameter will be ',
     +            'increased multiplying by rhofrac = ',1PD11.4,'.')
 1220 format(/,1X,'The desired infeasibility was not achieved in ',I7,
     +         1X,'over ',I7,/,1X,'constraints.',/,1X,' Penalty ',
     +            'parameters will be increased multiplying by ',
     +            'rhofrac = ',1PD11.4,'.')
 1230 format(/,1X,'Desired feasibility improvement was achieved, ',
     +            'so, penalty parameters will not',/,1X,'be ',
     +            'modified.')

 1300 format(/,1X,'Flag of ALGENCAN: Solution was found.')
 1310 format(/,1X,'Flag of ALGENCAN: Maximum of iterations reached.')
 1320 format(/,1X,'Flag of ALGENCAN: The problem seems to be ',
     +            'infeasible.')

 1400 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,I3,1X,I1,1X,I6,1X,I6,1X,I2,1X,I7,1X,I7,1X,I2,1X,I7,
     +       1X,I7,1X,I7,0P,F8.2)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine comprhoini(f,m,c,equatn,rho)

C     SCALAR ARGUMENTS
      integer m
      double precision f

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision c(m),rho(m)

C     Consider the Augmented Lagrangian function
C
C         al = f(x) + \sum_i P(lambda_i,c_i(x),rho),
C
C     where
C
C         P(lambda,y,rho) = y ( lambda + 0.5 rho y ),
C
C     If c_i(x) is an equality constraint or lambda + rho y > 0, and
C
C         P(lambda,y,rho) = - 0.5 lambda^2 / rho,
C
C     otherwise.
C
C     Assuming that lambda_i = 0 for all i, it is clear that
C
C         P(lambda_i,c_i(x),rho) = 0.5 rho c_i(x)^2
C
C     and that the value of  rho that balances f(x) and
C
C     \sum_i P(lambda_i,c_i(x),rho) is given by
C
C     rho = f(x) / ( 0.5 \sum_i c_i(x)^2 ).

#include "machconst.com"

C     LOCAL SCALARS
      integer i
      double precision rhoini,sumc

      sumc = 0.0d0
      do i = 1,m
          if ( equatn(i) .or. c(i) .gt. 0.0d0 ) then
              sumc = sumc + 0.5d0 * c(i) ** 2
          end if
      end do

      rhoini = 10.0d0 * max( 1.0d0, abs( f ) ) / max( 1.0d0, sumc )

      rhoini = max( macheps12, min( rhoini, 1.0d0 / macheps12 ) )

      do i = 1,m
          rho(i) = rhoini
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision cnorm,cnormu
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision x(n)

#include "dim.par"
#include "graddat.com"
#include "scaling.com"

C     LOCAL SCALARS
      integer i

C     Subroutine sevalnal is always called by GENCAN with the same point
C     before calling this subroutine. Then, constraints are computed and
C     saved.

C     Compute infeasibility

      cnorm = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( c(i) ) )
          else
              cnorm = max( cnorm, c(i) )
          end if
      end do

      if ( scale ) then
          cnormu = 0.0d0
          do i = 1,m
              if ( equatn(i) ) then
                  cnormu = max( cnormu, abs( sc(i) * c(i) ) )
              else
                  cnormu = max( cnormu, sc(i) * c(i) )
              end if
          end do

      else
          cnormu = cnorm
      end if

      end
