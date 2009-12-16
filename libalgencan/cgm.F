C     ******************************************************************
C     ******************************************************************

      subroutine cgm(nind,x,m,lambda,rho,equatn,linear,l,u,g,delta,
     +eps,maxit,precond,d,rbdnnz,rbdind,rbdtype,iter,cginfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      character * 6 precond
      integer cginfo,inform,iter,m,maxit,nind,rbdnnz
      double precision delta,eps

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer rbdind(nind)
      character rbdtype(nind)
      double precision d(nind),g(nind),l(nind),lambda(m),rho(m),u(nind),
     +        x(*)

C     This subroutine implements the Conjugate Gradients method for 
C     minimizing the quadratic approximation q(d) of L(x,lambda,rho) 
C     at x
C
C     q(d) = 1/2 d^T H d + g^T d,
C
C     where H is an approximation of the Hessian matrix of the 
C     Augmented Lagrangian and g is its gradient vector,
C
C     subject to || d || <= delta and l <= x + d <= u.
C
C     In the constraint ''|| d || <= delta'', the norm will be the
C     Euclidian-norm if the input parameter trtype is equal to 0, and
C     it will be the sup-norm if trtype is equal to 1.
C
C     The method returns an approximation d of the solution such that
C 
C     (a) ||H d + g||_2 <= eps * ||g||_2, 
C
C     (b) ||d|| = delta or x + d is in the boundary of the box, or
C
C     (c) ( p such that p^t H p = 0 ) and ( d = - amax g if such p was 
C         found during the first CG iteration or the current point d 
C         of CG if such p was found in any other iteration ).
C
C     inform   integer
C              termination parameter:
C
C              0 = convergence with ||H d + g||_2 <= eps * ||g||_2;
C
C              1 = convergence to the boundary of ||d|| <= delta;
C
C              2 = convergence to the boundary of l <= x + d <= u;
C
C              3 = stopping with d = dk  such that <gk,dk> <= - theta 
C                  ||gk||_2 ||dk||_2 and <gk,d_{k+1}> > - theta 
C                  ||gk||_2 ||d_{k+1}||_2;
C
C              4 = not enough progress of the quadratic model during
C                  maxitnqmp iterations, i.e., during maxitnqmp 
C                  iterations | q - qprev | <= macheps * max( | q |, 1 )
C
C              6 = very similar consecutive iterates, for two 
C                  consecutive iterates x1 and x2 we have that
C
C                  | x2(i) - x1(i) | <= macheps * max ( | x1(i) |, 1 )
C
C                  for all i.
C
C              7 = stopping with p such that p^T H p = 0 and g^T p = 0;
C
C              8 = too many iterations;

#include "dim.par"
#include "machconst.com"
#include "algconst.par"
#include "algparam.com"
#include "outtyp.com"

C     LOCAL SCALARS
      character * 6 prectmp
      logical goth,gotp,negcur,restarted,samep
      integer i,itertmp,itnqmp
      double precision alpha,abox,amax,amaxbds,amaxbdsi,amaxdel,
     +        bestprog,bbeta,currprog,dnorm2,gnorm2,gtd,gtp,pnorm2,ptd,
     +        pthp,ptr,q,qprev,rnorm2,ztrprev,ztr,znorm2

C     LOCAL ARRAYS
      double precision p(nmax),hp(nmax),r(nmax),z(nmax)

C     ==================================================================
C     Initialization
C     ==================================================================

      restarted = .false.

 001  continue

      goth = .false.
      gotp = .false.

C     gnorm2 = norm2s(nind,g)
      gnorm2 = 0.0d0
      do i = 1,nind
          gnorm2 = gnorm2 + g(i) ** 2
      end do

      iter     =      0
      itnqmp   =      0
      qprev    = bignum
      bestprog =  0.0d0

      do i = 1,nind
          d(i) = 0.0d0
          r(i) =  g(i)
      end do

      q        =  0.0d0
      gtd      =  0.0d0
      dnorm2   =  0.0d0
      rnorm2   = gnorm2

      ztr      =  0.0d0

C     ==================================================================
C     Print initial information
C     ==================================================================

      if ( iprintinn .ge. 5 ) then
          write(*, 980) maxit,eps,delta,precond,hptype
          write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q

          write(file10_unit,980) maxit,eps,delta,precond,hptype
          write(file10_unit,984) iter,sqrt(rnorm2),sqrt(dnorm2),q
      end if

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     if ||r||_2 = ||H d + g||_2 <= eps * ||g||_2 then stop

      if ( iter .ne. 0 .and. 
     +   ( ( rnorm2 .le. eps ** 2 * gnorm2 .and. iter .ge. 4 ) .or. 
     +     ( rnorm2 .le. macheps ) ) ) then 

          cginfo = 0

          if ( iprintinn .ge. 5 ) then
              write(*, 990)
              write(file10_unit,990)
          end if
  
          go to 500

      end if

C     if the maximum number of iterations was achieved then stop

      if ( iter .ge. max(4, maxit) ) then

          cginfo = 8

          if ( iprintinn .ge. 5 ) then
              write(*, 998)
              write(file10_unit,998)
          end if
  
          go to 500

      end if

C     ==================================================================
C     Preconditioner
C     ==================================================================

      if ( precond .eq. 'NONE' ) then

          do i = 1,nind
              z(i) = r(i)
          end do

          ztrprev = ztr
          ztr     = rnorm2
          znorm2  = rnorm2

      else if ( precond .eq. 'QNCGNA' ) then

          call capplyhpre(nind,m,rho,equatn,gotp,r,z)

          ztrprev = ztr

          ztr = 0.0d0
          do i = 1,nind
              ztr = ztr + z(i) * r(i)
          end do

C         znorm2 = norm2s(nind,z)
          znorm2 = 0.0d0
          do i = 1,nind
              znorm2 = znorm2 + z(i) ** 2
          end do

      end if

C     ==================================================================
C     Compute direction
C     ==================================================================

      if ( iter .eq. 0 ) then

          do i = 1,nind
              p(i) = - z(i)
          end do

          ptr    = - ztr
          pnorm2 =   znorm2

      else

          bbeta = ztr / ztrprev

          do i = 1,nind
              p(i) = - z(i) + bbeta * p(i)
          end do

          if ( precond .eq. 'NONE' ) then

              pnorm2 = rnorm2 - 2.0d0 * bbeta * ( ptr + alpha * pthp ) 
     +               + bbeta ** 2 * pnorm2
              ptr = - rnorm2 + bbeta * ( ptr + alpha * pthp )

          else if ( precond .eq. 'QNCGNA' ) then

              ptr = 0.0d0
              pnorm2 = 0.0d0
              do i = 1,nind
                  ptr = ptr + p(i) * r(i)
                  pnorm2 = pnorm2 + p(i) ** 2
              end do

          end if

      end if

C     Force p to be a descent direction of q(d), i.e.,
C     <\nabla q(d), p> = <H d + g, p> = <r, p> \le 0.

      if ( ptr .gt. 0.0d0 ) then

          do i = 1,nind
              p(i) = - p(i)
          end do

          ptr = - ptr

      end if

C     ==================================================================
C     Compute p^T H p
C     ==================================================================

C     hp = H p

      call calchalp(nind,x,m,lambda,rho,equatn,linear,p,hp,goth,inform)
      if ( inform .lt. 0 ) return

C     Compute p^T hp

      pthp = 0.0d0
      do i = 1,nind
          pthp = pthp + p(i) * hp(i)
      end do 

C     ==================================================================
C     Compute maximum steps
C     ==================================================================

      amaxdel = bignum

      do i = 1,nind
          if ( p(i) .gt. 0.0d0 ) then
              amaxdel = min( amaxdel,  (   delta - d(i) ) / p(i) )
          else if ( p(i) .lt. 0.0d0 ) then
              amaxdel = min( amaxdel,  ( - delta - d(i) ) / p(i) )
          end if
      end do

      amaxbds = bignum

C     do i = 1,nind
C        if ( p(i) .gt. 0.0d0 ) then

C            amaxbdsi = ( u(i) - d(i) - x(i) ) / p(i)

C            if ( amaxbdsi .lt. amaxbds ) then
C                amaxbds    = amaxbdsi
C                rbdnnz     = 1
C                rbdind(1)  = i
C                rbdtype(1) = 'U'
C            else if ( amaxbdsi .eq. amaxbds ) then
C                rbdnnz          = rbdnnz + 1
C                rbdind(rbdnnz)  = i
C                rbdtype(rbdnnz) = 'U'
C            end if

C        else if ( p(i) .lt. 0.0d0 ) then

C            amaxbdsi = ( l(i) - d(i) - x(i) ) / p(i)

C            if ( amaxbdsi .lt. amaxbds ) then
C                amaxbds    = amaxbdsi
C                rbdnnz     = 1
C                rbdind(1)  = i
C                rbdtype(1) = 'L'
C            else if ( amaxbdsi .eq. amaxbds ) then
C                rbdnnz          = rbdnnz + 1
C                rbdind(rbdnnz)  = i
C                rbdtype(rbdnnz) = 'L'
C            end if

C        end if
C     end do

      amax = min( amaxdel, amaxbds )

C     ==================================================================
C     Compute the step
C     ==================================================================

      negcur = .false.

C     If p^T H p > 0 then take the conjugate gradients step

      if ( pthp .gt. 0.0d0 ) then

          alpha = min( amax, ztr / pthp )

C     Else, if we are at iteration zero then take the maximum 
C     positive step in the minus gradient direction

      else if ( iter .eq. 0 ) then

          alpha = amax

          negcur = .true.

C     Otherwise, stop at the current iterate

      else

          cginfo = 7

          if ( iprintinn .ge. 5 ) then
              write(*, 997)
              write(file10_unit,997)
          end if
  
          go to 500

      end if

C     ==================================================================
C     Test the angle condition
C     ==================================================================

      ptd = 0.0d0
      gtp = 0.0d0
      do i = 1,nind
          ptd = ptd + p(i) * d(i)
          gtp = gtp + g(i) * p(i)
      end do

C     These are gtd and dnorm2 for the new direction d which was not 
C     computed yet.
      gtd = gtd + alpha * gtp
      dnorm2 = dnorm2 + alpha ** 2 * pnorm2 + 2.0d0 * alpha * ptd

      if ( gtd .gt. 0.0d0 .or. 
     +     gtd ** 2 .lt. theta ** 2 * gnorm2 * dnorm2 ) then

          if ( precond .ne. 'NONE' .and. iter .eq. 0 ) then

              if ( iprintinn .ge. 5 ) then
                  write(*, 986)
                  write(file10_unit,986)
              end if

              restarted = .true.
              itertmp   = iter
              prectmp   = precond
              precond   = 'NONE'
              go to 001

          end if

          cginfo = 3

          if ( iprintinn .ge. 5 ) then
              write(*, 993)
              write(file10_unit,993)
          end if

          go to 500

      end if

C     ==================================================================
C     Compute the quadratic model functional value at the new point
C     ==================================================================

      qprev = q

      q = q + 0.5d0 * alpha ** 2 * pthp + alpha * ptr

C     ==================================================================
C     Compute new d
C     ==================================================================

      do i = 1,nind
          d(i) = d(i) + alpha * p(i)
      end do

C     ==================================================================
C     Compute the residual r = H d + g
C     ==================================================================

      do i = 1,nind
          r(i) = r(i) + alpha * hp(i)
      end do

C     rnorm2 = norm2s(nind,r)
      rnorm2 = 0.0d0
      do i = 1,nind
          rnorm2 = rnorm2 + r(i) ** 2
      end do

C     ==================================================================
C     Increment number of iterations
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     Print information of this iteration
C     ==================================================================

      if ( iprintinn .ge. 5 ) then
          write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q
          write(file10_unit,984) iter,sqrt(rnorm2),sqrt(dnorm2),q
      end if

C     ==================================================================
C     Test other stopping criteria
C     ==================================================================

C     Boundary of the "trust region"

      if ( alpha .eq. amaxdel ) then

          cginfo = 1

          if ( iprintinn .ge. 5 ) then
              if ( negcur ) then
                  write(*, 987)
                  write(file10_unit,987)
              end if

              write(*, 991)
              write(file10_unit,991)
          end if
  
          go to 500

      end if

C     Boundary of the box constraints

C     if ( alpha .eq. amaxbds ) then

C         cginfo = 2

C         if ( iprintinn .ge. 5 ) then
C             if ( negcur ) then
C                 write(*, 987)
C                 write(file10_unit,987)
C             end if

C             write(*, 992)
C             write(file10_unit,992)
C         end if
  
C         go to 500

C     end if

C     Small useful proportion

      abox = bignum

      do i = 1,nind
         if ( d(i) .gt. 0.0d0 ) then
             abox = min( abox, ( u(i) - x(i) ) / d(i) )
         else if ( d(i) .lt. 0.0d0 ) then
             abox = min( abox, ( l(i) - x(i) ) / d(i) )
         end if
      end do

      if ( abox .le. 0.1d0 ) then

          cginfo = 5

          if ( iprintinn .ge. 5 ) then
              write(* ,995)
              write(file10_unit,995)
          end if
  
          go to 500

      end if

C     Two consecutive iterates are too much close

      samep = .true.
      do i = 1,nind
         if ( abs( alpha * p(i) ) .gt. 
     +        macheps * max( 1.0d0, abs( d(i) ) ) ) then 
             samep = .false.
          end if
      end do

      if ( samep ) then

          cginfo = 6

          if ( iprintinn .ge. 5 ) then
              write(*, 996)
              write(file10_unit,996)
          end if
  
          go to 500

      end if

C     Many iterations without good progress of the quadratic model 

      currprog = qprev - q
      bestprog = max( currprog, bestprog )

      if ( currprog .le. epsnqmp * bestprog ) then

          itnqmp = itnqmp + 1

          if ( itnqmp .ge. maxcgitnp ) then
              cginfo = 4

              if ( iprintinn .ge. 5 ) then
                  write(*, 994)
                  write(file10_unit,994)
              end if

              go to 500
          endif

      else
          itnqmp = 0
      endif

C     ==================================================================
C     Iterate
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     ==================================================================
C     Return
C     ==================================================================

 500  continue

C     Print final information

      if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
          write(*, 985) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
          write(file10_unit,985) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
      end if

      if ( restarted ) then
          iter = iter + itertmp
          precond = prectmp
      end if

      return

C     Non-executable statements

 980  format(/,5X,'Conjugate Gradients (maxit = ',I7,',',1X,'eps = ',
     +            1P,D7.1,',',1X,'delta = ',1P,D7.1,')',
     +       /,5X,'(Preconditioner: ',A6,',',1X,
     +            'Hessian-vector product type: ',A6,')')
 984  format(  5X,'CG iter = ',I7,' rnorm = ',1P,D10.4,' dnorm = ',
     +            1P,D10.4,' q = ',1P,D11.4)
 985  format(/,5X,'Truncated Newton direction (first ',I7,
     +            ' components): ',/,1(5X,6(1X,1P,D11.4)))
 986  format(  5X,'The first CG-PREC iterate did not satisfy the ',
     +            'angle condition. CG will be restarted without ',
     +            'preconditioner)')
 987  format(  5X,'p such that p^T H p = 0 was found. ',
     +            'Maximum step was taken.')

 990  format(  5X,'Flag of CG: Convergence with small residual.')
 991  format(  5X,'Flag of CG: Convergence to the trust region ',
     +            'boundary.')
C992  format(  5X,'Flag of CG: Convergence to the box boundary.')
 993  format(  5X,'Flag of CG: The next CG iterate will not satisfy ',
     +            'the angle condition.')
 994  format(  5X,'Flag of CG: Not enough progress in the quadratic ',
     +            'model.')
 995  format(  5X,'Flag of CG: The maximum step to remain within the ',
     +            'box is smaller than 0.1)')
 996  format(  5X,'Flag of CG: Very near consecutive iterates.')
 997  format(  5X,'Flag of CG: p such that p^T H p = 0 was found.')
 998  format(  5X,'Flag of CG: Maximum number of CG iterations ',
     +            'reached.')

      end
