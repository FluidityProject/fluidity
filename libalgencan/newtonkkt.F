C     ******************************************************************
C     ******************************************************************

      subroutine newtonkkt(n,xo,lo,uo,m,lambdao,equatno,linearo,epsfeas,
     +epsopt,f,cnorm,nlnorm,iter,msqiter,accinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer accinfo,inform,iter,m,msqiter,n
      double precision cnorm,epsfeas,epsopt,f,nlnorm

C     ARRAY ARGUMENTS
      logical equatno(m),linearo(m)
      double precision lo(n),lambdao(m),uo(n),xo(n)

#include "dim.par"
#include "machconst.com"
#include "algparam.com"
#include "outtyp.com"
#include "graddat.com"
#include "hessdat.com"

C     accinfo:
C
C     0: KKT system solved.
C     1: Ignored constraints were violated.
C     2: After correcting the Lagrange multipliers signal, optimality
C        was lost.
C     3: Maximum number of iterations reached.
C     4: Newton seems to be diverging.
C     5: Singular Jacobian.
C     6: Insufficient space to store the KKT linear system.
C     7: Insufficient double precision working space for linear solver.
C     8: Insufficient integer working space for linear solver.

C     PARAMETERS
      integer mtotmax,ntotmax
      parameter ( mtotmax  = mmax + 2 * nmax )
      parameter ( ntotmax  = mmax + 3 * nmax )

C     LOCAL SCALARS
      integer cind,col,divit,fun,hlnnz,i,itmp,j,jcnnz,k,lssinfo,lin,
     +        maxit,minp,mnop,mtot,nbds,nineq,ninn,nneigv,nnon,nsys,
     +        ntot,pind,sind,var,vind
      double precision cnorm0,cnormprev,cnormref,epsact,epsadd,nlnorm0,
     +        nlnormprev,nlnormref,pval,sgnvio,val

C     LOCAL ARRAYS
      logical equatn(mtotmax),linear(mmax)
      character constt(2*nmax),status(nmax)
      integer consti(2*nmax),hdiag(ntotmax+mtotmax),inn(ntotmax),
     +        inp(mtotmax),jcfun(jcnnzmax),slaind(mtotmax)
      double precision adddiag(ntotmax+mtotmax),b(ntotmax+mtotmax),
     +        l(ntotmax),lambda(mtotmax),nl(ntotmax),r(mtotmax),
     +        t(nmax),u(ntotmax),x(ntotmax)

C     ==================================================================
C     PRESENTATION
C     ==================================================================

      if ( iprintout .ge. 1 ) then
          write(* ,8000)
          write(file10_unit,8000)
      end if

C     ==================================================================
C     INITIALIZE
C     ==================================================================

      iter    =  0
      divit   =  0
      maxit   = 10
      msqiter =  0

      epsact = sqrt( epsfeas )
      epsadd = macheps12

      call lssini(sclsys,.true.,.false.)

C     ==================================================================
C     SET INITIAL POINT
C     ==================================================================

      do i = 1,n
          l(i) = lo(i)
          u(i) = uo(i)
          x(i) = max( lo(i), min( xo(i), uo(i) ) )
          if ( lo(i) .eq. uo(i) ) then
              write(*,*) 'There is a fixed variable: ',i,lo(i),xo(i)
          end if
      end do

      do j = 1,m
          lambda(j) = lambdao(j)
          equatn(j) = equatno(j)
          linear(j) = linearo(j)
      end do

C     ==================================================================
C     COMPUTE CONSTRAINTS AND GRADIENT OF THE LAGRANGIAN
C     ==================================================================

      call ssetp(n,x)

C     Objective function and constraints

      call sevalobjc(n,x,f,m,r,inform)
      if ( inform .lt. 0 ) return

C     Gradient of the Lagrangian (and save Jacobian of the constraints)

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     STRUCTURES FOR NEW VARIABLES AND (BOUND) CONSTRAINTS
C     SET SLACKS VARIABLES VALUES
C     ==================================================================

C     Relate constraints to slacks

      nineq = 0
      do j = 1,m
          if ( .not. equatn(j) ) then
              nineq     = nineq + 1
              sind      = n + nineq
              slaind(j) = sind
              x(sind)   = sqrt( 2.0d0 * max( 0.0d0, - r(j) ) )
          else
              slaind(j) = 0
          end if
      end do

      nbds = 0
      do i = 1,n
          if ( l(i) .gt. - 1.0d+20 ) then
              nbds = nbds + 1
              constt(nbds) = 'L'
              consti(nbds) =  i
              cind         = m + nbds
              sind         = n + nineq + nbds
              slaind(cind) = sind
              r(cind)      = l(i) - x(i)
              x(sind)      = sqrt( 2.0d0 * max( 0.0d0, - r(cind) ) )
              r(cind)      = r(cind) + 0.5d0 * x(sind) ** 2
              equatn(cind) = .true.
          end if

          if ( u(i) .lt. 1.0d+20 .and. u(i) .ne. l(i) ) then
              nbds = nbds + 1
              constt(nbds) = 'U'
              consti(nbds) =  i
              cind         = m + nbds
              sind         = n + nineq + nbds
              slaind(cind) = sind
              r(cind)      = x(i) - u(i)
              x(sind)      = sqrt( 2.0d0 * max( 0.0d0, - r(cind) ) )
              r(cind)      = r(cind) + 0.5d0 * x(sind) ** 2
              equatn(cind) = .true.
          end if
      end do

      mtot = m + nbds
      ntot = n + nineq + nbds

C     ==================================================================
C     MAIN LOOP
C     ==================================================================

 100  continue

C     ==================================================================
C     SET ACTIVE CONSTRAINTS AND VARIABLES
C     ==================================================================

C     A fixed variable will be fixed forever.

      minp = 0
      mnop = 0

      ninn = 0
      nnon = 0

C     Set active variables

      do i = 1,n
           if ( x(i) .le. l(i) ) then
               status(i) = 'L'
               x(i) = l(i)

               nnon   = nnon + 1
               inn(i) = ntot + 1 - nnon

           else if ( x(i) .ge. u(i) ) then

               status(i) = 'U'
               x(i) = u(i)

               nnon   = nnon + 1
               inn(i) = ntot + 1 - nnon

           else
              status(i) = 'F'
              ninn   = ninn + 1
              inn(i) = ninn
           end if
      end do

C     Set active (regular) constraints and their slacks

      do j = 1,m
          if ( equatn(j) ) then
C             Active equality constraint
              minp   = minp + 1
              inp(j) = minp
          else
              sind = slaind(j)

              if ( r(j) .ge. - epsact ) then
C                 Active inequality constraint and its slack
                  minp      = minp + 1
                  inp(j)    = minp

                  ninn      = ninn + 1
                  inn(sind) = ninn
              else
C                 Deactivate inequality constraint and its slack
                  mnop      = mnop + 1
                  inp(j)    = mtot + 1 - mnop

                  nnon      = nnon + 1
                  inn(sind) = ntot + 1 - nnon
              end if
          end if
      end do

C     Set active bound constraints and their slacks

      do i = 1,nbds
          cind = m + i
          vind = consti(i)
          sind = slaind(cind)

          if ( status(vind) .eq. 'F' ) then
              minp      = minp + 1
              inp(cind) = minp

              ninn      = ninn + 1
              inn(sind) = ninn

          else
              mnop      = mnop + 1
              inp(cind) = mtot + 1 - mnop

              nnon      = nnon + 1
              inn(sind) = ntot + 1 - nnon

              if ( constt(i) .eq. status(vind) ) then
                  x(sind) = 0.0d0
              else
                  x(sind) = sqrt( 2.0d0 * ( u(vind) - l(vind) ) )
              end if
          end if
      end do

      nsys = ninn + minp

C     ==================================================================
C     SAVE NORMS TO CHECK IMPROVEMENT
C     ==================================================================

      if ( iter .eq. 0 ) then
          nlnormprev = bignum
          cnormprev  = bignum
      else
          nlnormprev = nlnorm
          cnormprev  = cnorm
      end if

C     ==================================================================
C     COMPUTE OBJECTIVE FUNCTION AND CONSTRAINTS
C     ==================================================================

      call ssetp(n,x)

C     Objective function and constraints

      call sevalobjc(n,x,f,m,r,inform)
      if ( inform .lt. 0 ) return

C     Add slacks effect

C     Regular constraints

      do j = 1,m
          sind = slaind(j)
          if ( sind .ne. 0 ) then
              r(j) = r(j) + 0.5d0 * x(sind) ** 2
          end if
      end do

C     Bound constraints

      do i = 1,nbds
          cind = m + i
          sind = slaind(cind)
          vind = consti(i)
          if ( constt(i) .eq. 'L' ) then
              r(cind) = l(vind) - x(vind)
          else
              r(cind) = x(vind) - u(vind)
          end if

          r(cind) = r(cind) + 0.5d0 * x(sind) ** 2
      end do

C     Constraints norm

      cnorm = 0.0d0
      do j = 1,mtot
          cnorm = max( cnorm, abs( r(j) ) )
      end do

      if ( iter .eq. 0 ) then
          cnorm0 = cnorm
      end if

C     ==================================================================
C     COMPUTE FIRST DERIVATIVES
C     ==================================================================

C     Gradient of the objective function and Jacobian of the constraints

      if ( gjaccoded ) then

          call sevalgjac(n,x,t,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

C         First derivative related to the slacks of the regular
C         constraints

          do j = 1,m
              sind = slaind(j)

              if ( sind .ne. 0 ) then
                  jcnnz = jcnnz + 1

                  jcfun(jcnnz) = j
                  jcvar(jcnnz) = sind
                  jcval(jcnnz) = x(sind)
              end if
          end do

          call coo2csr(m,jcnnz,jcfun,jcvar,jcval,jclen,jcsta)

      else

C         Gradient of the objective function

          call sevalg(n,x,t,inform)
          if ( inform .lt. 0 ) return

C         Jacobian of regular constraints with the slacks effect

          k = 0

          do j = 1,m
              jcsta(j) = k + 1

              call sevaljac(n,x,j,jcvar(k+1),jcval(k+1),jclen(j),inform)
              if ( inform .lt. 0 ) return

              k = k + jclen(j)

C             First derivative related to the slacks of the regular
C             constraints

              sind = slaind(j)

              if ( sind .ne. 0 ) then
                  jclen(j) = jclen(j) + 1

                  jcvar(k+1) = sind
                  jcval(k+1) = x(sind)

                  k = k + 1
              end if
          end do

          jcnnz = k

      end if

C     Bound constraints with the slacks effect

      do i = 1,nbds
          cind = m + i
          sind = slaind(cind)
          vind = consti(i)

          jcsta(cind) = jcnnz + 1
          jclen(cind) = 2

          jcvar(jcnnz+1) = vind
          if ( constt(i) .eq. 'L' ) then
              jcval(jcnnz+1) = - 1.0d0
          else
              jcval(jcnnz+1) =   1.0d0
          end if

          jcvar(jcnnz+2) = sind
          jcval(jcnnz+2) = x(sind)

          jcnnz = jcnnz + jclen(cind)
      end do

C     Gradient of the Lagrangian

      do i = 1,ntot
          nl(i) = 0.0d0
      end do

C     Gradient of the objective function

      do i = 1,n
          nl(i) = nl(i) + t(i)
      end do

C     Effect of the Jacobian of regular constraints

      do j = 1,m
          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
          end do
      end do

C     Set Lagrange multipliers related to bound constraints

      do i = 1,nbds
          cind = m + i
          vind = consti(i)
          if ( status(vind) .ne. 'F' ) then
              if ( constt(i) .eq. status(vind) ) then
                  if ( constt(i) .eq. 'L' ) then
                      lambda(cind) =   nl(vind)
                  else
                      lambda(cind) = - nl(vind)
                  end if
              else
                  lambda(cind) = 0.0d0
              end if
          else
              lambda(cind) = 0.0d0
          end if
      end do

C     Effect of Jacobian of bound constraints

      do j = m + 1,m + nbds
          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
          end do
      end do

C     Gradient of the Lagrangian norm

      nlnorm = 0.0d0
      do i = 1,ntot
          nlnorm = max( nlnorm, abs( nl(i) ) )
      end do

      if ( iter .eq. 0 ) then
          nlnorm0 = nlnorm
      end if

C     ==================================================================
C     WRITE INFORMATION OF THE CURRENT POINT
C     ==================================================================

      if ( iprintout .ge. 1 ) then
          write(* ,8010) iter,f,cnorm,nlnorm
          write(file10_unit,8010) iter,f,cnorm,nlnorm
      end if

C     ==================================================================
C     TEST STOPPING CRITERIA
C     ==================================================================

      if ( cnorm .le. epsfeas .and. nlnorm .le. epsopt ) then
        ! THE POINT SATISFIES FEASIBILITY AND OPTIMALITY

        ! (Ignored constraints and Lagrange multipliers signal
        !  constraint must be checked)

          go to 400
      end if

      if ( iter .ge. maxit ) then
        ! MAXIMUM NUMBER OF ITERATIONS EXCEEDED

          accinfo = 3

          if ( iprintout .ge. 1 ) then
              write(* ,9030)
              write(file10_unit,9030)
          end if

          go to 500
      end if

      cnormref  = max(cnorm0 ,max(cnormprev ,max(epsfeas, 1.0d+01)))
      nlnormref = max(nlnorm0,max(nlnormprev,max(epsopt,  1.0d+01)))

      if ( cnorm .gt. cnormref  .or. nlnorm .gt. nlnormref .or.
     +     cnorm .eq. cnormprev .or. nlnorm .eq. nlnormprev ) then
        ! IT SEEMS TO BE DIVERGING

          divit = divit + 1

          if ( divit .ge. 3 ) then

              accinfo = 4

              if ( iprintout .ge. 1 ) then
                  write(* ,9040)
                  write(file10_unit,9040)
              end if

              go to 500

          end if

      else
          divit = 0
      end if

C     ==================================================================
C     DO AN ITERATION
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     COMPUTE SECOND DERIVATIVES
C     ==================================================================

C     Hessian of the Lagrangian

      call sevalhl(n,x,m,lambda,hlin,hcol,hval,hlnnz,inform)
      if ( inform .lt. 0 ) return

C     Second derivatives related to slacks of the regular constraints

      do j = 1,m
          sind = slaind(j)

          if ( sind .ne. 0 ) then
              hlnnz = hlnnz + 1

              hlin(hlnnz) = sind
              hcol(hlnnz) = sind
              hval(hlnnz) = lambda(j)
          end if
      end do

C     Second derivatives related to slacks of the bound constraints

      do i = 1,nbds
          cind = m + i
          sind = slaind(cind)

          hlnnz = hlnnz + 1

          hlin(hlnnz) = sind
          hcol(hlnnz) = sind
          hval(hlnnz) = lambda(cind)
      end do

C     ==================================================================
C     ASSEMBLE THE JACOBIAN OF THE KKT SYSTEM
C     ==================================================================

      hnnz = 0

      do i = 1,nsys
          hdiag(i) = 0
      end do

C     Hessian of the Lagrangian

      do i = 1,hlnnz
          if ( hlin(i) .ge. hcol(i) ) then

              lin = inn(hlin(i))
              col = inn(hcol(i))
              val = hval(i)

              if ( lin .le. ninn .and. col .le. ninn ) then
                  if ( val .ne. 0.0d0 ) then
C                     A(lin,col) = A(lin,col) + val
                      hnnz = hnnz + 1
                      hlin(hnnz) = lin
                      hcol(hnnz) = col
                      hval(hnnz) = val
                      if ( lin .eq. col ) hdiag(lin) = hnnz
                  end if
              end if

          end if
      end do

C     Jacobian of the constraints

      do j = 1,mtot
          do i = jcsta(j),jcsta(j) + jclen(j) - 1

              fun = inp(j)
              var = inn(jcvar(i))
              val = jcval(i)

              if ( var .le. ninn .and. fun .le. minp ) then
                  if ( val .ne. 0.0d0 ) then
C                     A(fun+ninn,var) = A(fun+ninn,var) + val
                      hnnz = hnnz + 1
                      hlin(hnnz) = fun + ninn
                      hcol(hnnz) = var
                      hval(hnnz) = val
                  end if
              end if

          end do
      end do

      do i = 1,nsys
          if ( hdiag(i) .eq. 0 ) then
              hnnz = hnnz + 1
              hlin(hnnz) = i
              hcol(hnnz) = i
              hval(hnnz) = 0.0d0

              hdiag(i) = hnnz
          end if
      end do

C     ==================================================================
C     ANALYSE SPARSITY PATTERN
C     ==================================================================

      call lssana(nsys,hnnz,hlin,hcol,hval,hdiag,lssinfo)

      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          accinfo = 6
          go to 500

      end if

C     ==================================================================
C     SOLVE THE NEWTONIAN SYSTEM
C     ==================================================================

C     ==================================================================
C     COMPUTE REGULARIZATION
C     ==================================================================

 200  continue

      do i = 1,ninn
          adddiag(i) = epsadd
      end do

      do i = ninn + 1,nsys
          adddiag(i) = - macheps12
      end do

C     ==================================================================
C     FACTORIZE THE JACOBIAN OF THE NEWTONIAN SYSTEM
C     ==================================================================

      call lssfac(nsys,hnnz,hlin,hcol,hval,hdiag,adddiag,pind,pval,
     +nneigv,lssinfo)

      if ( lssinfo .eq. 0 .or. lssinfo .eq. 1 ) then

          if ( nneigv .ne. minp ) then
            ! WRONG INERTIA (SEE NOCEDAL AND WRIGHT)

C             Lemma 16.3 [pg. 447]: Suppose the the Jacobian of the
C             constraints has full rank and that the reduced Hessian
C             Z^T H Z is positive definite. Then the Jacobian of the
C             KKT system has ninn positive eigenvalues, minp negative
C             eigenvalues, and no zero eigenvalues.

C             Note that at this point we know that the matrix has no
C             zero eigenvalues. nneigv gives the number of negative
C             eigenvalues.

              epsadd = max( macheps12, epsadd * 10.0d0 )

              if ( iprintout .ge. 1 ) then
                  itmp = ninn+minp-nneigv
                  write(* ,8090) ninn,minp,itmp,nneigv,epsadd
                  write(file10_unit,8090) ninn,minp,itmp,nneigv,epsadd
              end if

              go to 200
          end if

      else if ( lssinfo .eq. 2 ) then
        ! SINGULAR JACOBIAN

          epsadd = max( macheps12, epsadd * 10.0d0 )

          if ( iprintout .ge. 1 ) then
              write(* ,8080) epsadd
              write(file10_unit,8080) epsadd
          end if

          if ( epsadd .le. 1.0d+20 ) then
              go to 200
          end if

          accinfo = 5

          if ( iprintout .ge. 1 ) then
              write(* ,9050)
              write(file10_unit,9050)
          end if

          go to 500

      else if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          accinfo = 6
          go to 500

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          accinfo = 7
          go to 500

      else ! if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          accinfo = 8
          go to 500

      end if

C     ==================================================================
C     SOLVE TRIANGULAR SYSTEMS
C     ==================================================================

C     SET RHS

      do i = 1,ntot
          if ( inn(i) .le. ninn ) then
              b(inn(i)) = - nl(i)
          end if
      end do

      do j = 1,mtot
          if ( inp(j) .le. minp ) then
              b(ninn+inp(j)) = - r(j)
          end if
      end do

C     SOLVE THE EQUATIONS

      call lsssol(nsys,b)

C     ==================================================================
C     UPDATE x AND lambda
C     ==================================================================

      do i = 1,ntot
          if ( inn(i) .le. ninn ) then
              x(i) = x(i) + b(inn(i))
          end if
      end do

      do j = 1,mtot
          if ( inp(j) .le. minp ) then
              lambda(j) = lambda(j) + b(ninn+inp(j))
          end if
      end do

C     ==================================================================
C     ITERATE
C     ==================================================================

      go to 100

C     ==================================================================
C     END OF MAIN LOOP
C     ==================================================================

C     ==================================================================
C     CHECK FEASIBILITY CONSIDERING ALL CONSTRAINTS
C     ==================================================================

 400  continue

      cnorm = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( r(i) ) )
          else
              cnorm = max( cnorm, r(i) )
          end if
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,8020) cnorm
          write(file10_unit,8020) cnorm
      end if

      if ( cnorm .gt. epsfeas ) then
          accinfo = 1

          if ( iprintout .ge. 1 ) then
              write(* ,9010)
              write(file10_unit,9010)
          end if

          go to 500
      end if

C     ==================================================================
C     CHECK LAGRANGE MULTIPLIERS SIGNAL (RELATED TO BOUND CONSTRAINTS)
C     (Lagrange mutipliers of inequality constraints must be tested too.)
C     ==================================================================

      sgnvio = 0.0d0
      do i = m + 1,mtot
          sgnvio = max( sgnvio, - lambda(i) )
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,8030) sgnvio
          write(file10_unit,8030) sgnvio
      end if

      if ( sgnvio .eq. 0.0d0 ) then
          accinfo = 0

          if ( iprintout .ge. 1 ) then
              write(* ,9000)
              write(file10_unit,9000)
          end if

          go to 500
      end if

C     TRY TO CORRECT MULTIPLIERS SIGNAL

      if ( iprintout .ge. 1 ) then
          write(* ,8040)
          write(file10_unit,8040)
      end if

      call minsq(n,t,m,nbds,lambda,constt,consti,status,msqiter,inform)
      if ( inform .lt. 0 ) return

C     Compute signal violation

      sgnvio = 0.0d0
      do i = m + 1,mtot
         sgnvio = max( sgnvio, - lambda(i) )
      end do

C     Compute optimality

      do i = 1,n
          nl(i) = t(i)
      end do

      do i = n + 1,ntot
          nl(i) = 0.0d0
      end do

      do j = 1,mtot
          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
          end do
      end do

      nlnorm = 0.0d0
      do i = 1,ntot
          nlnorm = max( nlnorm, abs( nl(i) ) )
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,8050) sgnvio,nlnorm
          write(file10_unit,8050) sgnvio,nlnorm
      end if

      if ( nlnorm .le. epsopt ) then
          accinfo = 0

          if ( iprintout .ge. 1 ) then
              write(* ,9000)
              write(file10_unit,9000)
          end if

          go to 500

      else
          accinfo = 2

          if ( iprintout .ge. 1 ) then
              write(* ,9020)
              write(file10_unit,9020)
          end if

          go to 500
      end if

C     ==================================================================
C     SET SOLUTION
C     ==================================================================

 500  continue

      do i = 1,n
          xo(i) = x(i)
      end do

      do j = 1,m
          lambdao(j) = lambda(j)
      end do

C     ==================================================================
C     NON-EXECUTABLE STATEMENTS
C     ==================================================================

 8000 format(/,' NEWTON-KKT scheme in action!')
 8010 format(/,' NEWTON-KKT Iteration',42X,' = ',5X,I6,
     +       /,' Objective function value',38X,' = ',1PD11.4,
     +       /,' Maximal violation of constraints',30X,' = ',1PD11.4,
     +       /,' Sup-norm of the gradient of the Lagrangian',20X,' = ',
     +           1PD11.4)
 8020 format(/,' Maximal violation of constraints (considering all',
     +         ' constraints) = ',1PD11.4)
 8030 format(/,' Maximal violation of Lagrange multipliers',
     +         ' non-negativity',6X,' = ',1PD11.4)
 8040 format(/,' GENCAN is being called to find the right',
     +         ' multipliers.')
 8050 format(/,' Maximal violation of Lagrange multipliers',
     +         ' non-negativity',6X,' = ',1PD11.4,
     +       /,' Sup-norm of the gradient of the Lagrangian',20X,' = ',
     +           1PD11.4)

 8080 format(/,' Singular Jacobian.',
     +         ' epsadd was increased to ',1PD11.4)
 8090 format(/,' Wrong Jacobian inertia. ',
     +       /,' Desired POS = ',I7,' NEG = ',I7,
     +       /,' Actual  POS = ',I7,' NEG = ',I7,
     +       /,' epsadd was increased to ',1PD11.4)

 9000 format(/,' Flag of NEWTON-KKT = KKT system solved!')
 9010 format(/,' Flag of NEWTON-KKT = Ignored constraints were',
     +         ' violated.')
 9020 format(/,' Flag of NEWTON-KKT = After correcting the Lagrange',
     +         ' multipliers signal,',/,' optimality was lost.')
 9030 format(/,' Flag of NEWTON-KKT = Maximum of iterations reached.')
 9040 format(/,' Flag of NEWTON-KKT = Newton seems to be diverging.')
 9050 format(/,' Flag of NEWTON-KKT = Singular Jacobian.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsq(n,t,m,nbds,lambda,constt,consti,status,iter,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,iter,m,n,nbds

C     ARRAY ARGUMENTS
      character constt(nbds),status(n)
      integer consti(nbds)
      double precision t(n),lambda(m+nbds)

#include "dim.par"
#include "algparam.com"
#include "graddat.com"
#include "hessdat.com"

C     The problem to be solved is
C
C         min 0.5 || b + A^T lambda ||_2^2
C
C         subject to
C
C         lambda_i >= 0, i = m+1,...,m+nabds.
C
C     Columns a_i of A^T are multiplied by
C
C         an_i = 1 / max( 1, ||a_i||_infty ).
C
C     So, if we define D = diag(an_1,...,an_{m+nabds}), the problem
C     can be rewrittem as
C
C         min 0.5 || b + A^T D D^{-1} lambda ||_2^2
C
C         subject to
C
C         an_i lambda_i >= 0, i = m+1,...,m+nabds.
C
C     Subtituting an_i lambda_i by lambda_i the problem to be solved
C     becomes
C
C         min 0.5 || b + A^T D lambda ||_2^2
C
C         subject to
C
C         lambda_i >= 0, i = m+1,...,m+nabds.

C     PARAMETERS
      integer mtotmax,ntotmax
      parameter ( mtotmax  = mmax + 2 * nmax )
      parameter ( ntotmax  = mmax + 3 * nmax )

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      logical avoiddstmp,dum2,dum3,useustptmp
      integer cind,geninfo,i,j,k,maxit,nabds,vind
      double precision dum1,dum4,dum5,dum6,eps,msqf,msqnlpsupn

C     LOCAL ARRAYS
      integer ind(mtotmax)
      double precision l(mtotmax),msqnl(mtotmax),scol(mtotmax),
     +        u(mtotmax),x(mtotmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

C     RHS

      do i = 1,n
          b(i) = t(i)
      end do

C     Matrix

      hnnz = 0

      do j = 1,m
          ind(j) = j

          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              if ( jcvar(i) .le. n ) then
                  hnnz = hnnz + 1
                  hcol(hnnz) = j
                  hlin(hnnz) = jcvar(i)
                  hval(hnnz) = jcval(i)
              end if
          end do

          x(j) = lambda(j)
          l(j) = - 1.0d+20
          u(j) =   1.0d+20
      end do

      nabds = 0
      do j = 1,nbds
          cind = m + j
          vind = consti(j)

          if ( constt(j) .eq. status(vind) ) then
              nabds = nabds + 1
              k = m + nabds

              ind(k) = cind

              do i = jcsta(cind),jcsta(cind) + jclen(cind) - 1
                  if ( jcvar(i) .le. n ) then
                      hnnz = hnnz + 1
                      hcol(hnnz) = k
                      hlin(hnnz) = jcvar(i)
                      hval(hnnz) = jcval(i)
                  end if
              end do

              x(k) = lambda(cind)
              l(k) = 0.0d0
              u(k) = 1.0d+20
          end if
      end do

C     Dimensions

      nrows = n
      ncols = m + nabds

C     Columns scaling

      do j = 1,ncols
          scol(j) = 1.0d0
      end do

      do i = 1,hnnz
          scol(hcol(i)) = max( scol(hcol(i)), abs( hval(i) ) )
      end do

      do j = 1,ncols
          scol(j) = 1.0d0 / scol(j)
      end do

      do i = 1,hnnz
          hval(i) = hval(i) * scol(hcol(i))
      end do

      do i = 1,ncols
          x(i) = x(i) / scol(i)
      end do

C     Call the solver

      innercall = .true.

      avoiddstmp = avoidds
      useustptmp = useustp

      avoidds = .true.
      useustp = .true.

      eps   = 1.0d-16
      maxit = 200

      call gencan(ncols,x,l,u,0,dum1,dum2,dum3,dum4,0.0d0,eps,
     +maxit,iter,msqf,msqnl,msqnlpsupn,dum5,dum6,geninfo,inform)

      innercall = .false.

      avoidds = avoiddstmp
      useustp = useustptmp

C     Copy unscaled solution to lambda

      do i = 1,ncols
          lambda(ind(i)) = x(i) * scol(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsqf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "hessdat.com"

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision p(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          p(i) = b(i)
      end do

      do i = 1,hnnz
          p(hlin(i)) = p(hlin(i)) + x(hcol(i)) * hval(i)
      end do

      f = 0.0d0
      do i = 1,nrows
          f = f + p(i) ** 2
      end do

      f = 0.5d0 * f

      f = 1.0d+08 * f

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsqg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

#include "dim.par"
#include "hessdat.com"

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision p(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          p(i) = b(i)
      end do

      do i = 1,hnnz
          p(hlin(i)) = p(hlin(i)) + x(hcol(i)) * hval(i)
      end do

      do i = 1,ncols
          g(i) = 0.0d0
      end do

      do i = 1,hnnz
          g(hcol(i)) = g(hcol(i)) + p(hlin(i)) * hval(i)
      end do

      do i = 1,ncols
          g(i) = 1.0d+08 * g(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsqhp(n,x,p,hp,goth,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer inform,n

C     ARRAY ARGUMENTS
      double precision hp(n),p(n),x(n)

#include "dim.par"
#include "hessdat.com"

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision tmp(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          tmp(i) = 0.0d0
      end do

      do i = 1,hnnz
          tmp(hlin(i)) = tmp(hlin(i)) + p(hcol(i)) * hval(i)
      end do

      do i = 1,ncols
          hp(i) = 0.0d0
      end do

      do i = 1,hnnz
          hp(hcol(i)) = hp(hcol(i)) + tmp(hlin(i)) * hval(i)
      end do

      do i = 1,ncols
          hp(i) = 1.0d+08 * hp(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      logical function minsqstop(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "machconst.com"
#include "hessdat.com"

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i
      double precision pnorm

C     LOCAL ARRAYS
      double precision p(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          p(i) = b(i)
      end do

      do i = 1,hnnz
          p(hlin(i)) = p(hlin(i)) + x(hcol(i)) * hval(i)
      end do

      pnorm = 0.0d0
      do i = 1,nrows
          pnorm = max( pnorm, abs( p(i) ) )
      end do

      minsqstop = .false.
      if ( pnorm .le. macheps12 ) then
          minsqstop = .true.
      end if

      end
