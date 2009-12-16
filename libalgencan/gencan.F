C     ******************************************************************
C     ******************************************************************

      subroutine gencan(n,x,l,u,m,lambda,equatn,linear,rho,epsfeas,
     +epsopt,maxit,iter,f,g,gpsupn,cnorm,cnormu,geninfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer geninfo,inform,iter,m,maxit,n
      double precision cnorm,cnormu,epsfeas,epsopt,f,gpsupn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision g(n),l(n),lambda(m),rho(m),u(n),x(n)

C     Solves the box-constrained minimization problem
C
C            Minimize f(x) subject to l <= x <= u
C     
C     using the method described in 
C
C     E. G. Birgin and J. M. Martinez, ''Large-scale active-set box-
C     constrained optimization method with spectral projected 
C     gradients'', Computational Optimization and Applications 23, pp. 
C     101-125, 2002.  

C     geninfo:
C
C     0: Small continuous-projected-gradient norm
C     1: Maximum number of gencan iterations reached
C     3: Lack of progress in the objective function value
C     4: Lack of progress in the continuous-projected-gradient norm
C     5: Lack of progress in the current point
C     6: Unbounded objective function?
C     7: Too small backtracking step. Wrong gradient?

#include "dim.par"
#include "machconst.com"
#include "algconst.par"
#include "algparam.com"
#include "counters.com"
#include "outtyp.com"
#include "rspace.com"
#include "itetyp.com"
#include "sydat.com"

C     LOCAL SCALARS
      character * 6 precond
      character * 2 innittype
      logical forceoi,memfail,vustop,connected
      integer cgcnt,cginfo,cgiter,cgmaxit,chcnt,cinit,i,innit,iscnt,
     +        itnfp,itngp,itnxp,lsinfo,mfit,nind,nindprev,nwcnt,lvfit,
     +        rbdnnz,trcnt,trinfo,xind,cisit,file20_unit
      double precision acgeps,amax,bcgeps,cgeps,cgdel,dsupn,fdiff,fplus,
     +        gieucn,gisupn,gpeucn,gpi,lamspg,mslamb,newtrdel,ssupn,tmp,
     +        trdel,tsmall,ysupn,xeucn,xsupn

C     LOCAL ARRAYS
      integer rbdind(nmax)
      character rbdtype(nmax)
      double precision d(nmax),gplus(nmax),xplus(nmax)

C     DATA BLOCKS
      character * 19 ittext(0:5)
      integer BDSREACHED,FSTORDPNT,SECORDPNT,SMALLSTEP,SMALLTRRAD,
     +        UNBOUNDED,UNDEFSTEP

      data ittext(0) /'(Initial point)   '/
      data ittext(1) /'(SPG step)        '/
      data ittext(2) /'(Truncated Newton)'/
      data ittext(3) /'(MS system)       '/
      data ittext(4) /'(Trust region)    '/
      data ittext(5) /'(Inner SPG)       '/

      data BDSREACHED /2/
      data FSTORDPNT  /6/
      data SECORDPNT  /7/
      data SMALLSTEP  /3/
      data SMALLTRRAD /4/
      data UNBOUNDED  /2/
      data UNDEFSTEP  /5/

C     EXTERNAL
      logical sstop

C     EXTERNAL SUBROUTINES
      external sevalal

C     ==================================================================
C     Initialization
C     ==================================================================

C     Set some initial values:

C     to record a memory failure in the direct solver
      memfail = .false.

C     for Conjugate Gradients
      precond = 'QNCGNA'

C     for the first inner-to-the-face minimization algorithm
      innittype = 'CG'

C     for testing lack of progress checking f, g and x
      itnfp = 0
      itngp = 0
      itnxp = 0

C     to force a leaving-face iteration when lack of progress is 
C     detected within a face
      forceoi = .false.

C     for calculating More-Sorensen's direction
      mslamb =  0.0d0

C     for counting number of iterations as well as inner-to-the-face
C     and leaving-face iterations
      iter  = 0
      cgcnt = 0
      nwcnt = 0
      trcnt = 0
      iscnt = 0
      chcnt = 0
      lvfit = 0
      innit = 0
      cinit = 0
      cisit = 0

C     just to print "Initial point" in the first ouput
      ittype = 0

C     Print problem information

      if ( iprintinn .ge. 1 ) then
          write(* ,1000) n
          write(file10_unit,1000) n
      end if

      if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
          write(* ,1010) nprint,(l(i),i=1,nprint)
          write(* ,1020) nprint,(u(i),i=1,nprint)
          write(* ,1030) nprint,(x(i),i=1,nprint)

          write(file10_unit,1010) nprint,(l(i),i=1,nprint)
          write(file10_unit,1020) nprint,(u(i),i=1,nprint)
          write(file10_unit,1030) nprint,(x(i),i=1,nprint)
      end if

C     Project initial guess. If the initial guess is infeasible, 
C     projection puts it into the box.

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Compute function and gradient at the initial point

      call sevalal(n,x,m,lambda,rho,equatn,linear,f,inform)
      if ( inform .lt. 0 ) return

      call sevalnal(n,x,m,lambda,rho,equatn,linear,g,inform)
      if ( inform .lt. 0 ) return

C     Compute x Euclidian and sup norms, continuous-project-gradient 
C     Euclidian and sup norms, internal gradient Euclidian norm. Set
C     nind as the number of free variables and save in array ind their 
C     identifiers. Also set variable sameface (same face) indicating 
C     that the "previous iterate" does not belong to the current face.

      sameface = .false.

      nind   = 0
      xsupn  = 0.0d0
      xeucn  = 0.0d0
      gpsupn = 0.0d0
      gpeucn = 0.0d0
      gisupn = 0.0d0
      gieucn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind = i
          end if
          xeucn = xeucn + x(i) ** 2
          gpi = x(i) - g(i)
          if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
              gpi = - g(i)
          else
              gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
          end if
          gpsupn = max( gpsupn, abs( gpi ) )
          gpeucn = gpeucn + gpi ** 2
          if ( x(i) - l(i) .gt. macheps23 * max( 1.0d0, l(i) ) .and. 
     +         u(i) - x(i) .gt. macheps23 * max( 1.0d0, u(i) ) ) then
              gisupn    = max( gisupn, abs( gpi ) )
              gieucn    = gieucn + gpi ** 2
              nind      = nind + 1
              ind(nind) = i
          end if
      end do
      xeucn  = sqrt( xeucn  )
      gpeucn = sqrt( gpeucn )
      gieucn = sqrt( gieucn )

C     Compute trut-region radius, in case betra is used as inner solver

      trdel    = max( trdelmin, trdelini * max( 1.0d0, xeucn ) )
      newtrdel = 0.0d0
  
C     To be used in feasibility problems, compute constraints norm

      call sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)
      if ( inform .lt. 0 ) return

C     Initial spectral steplength

C     Compute a small step and set the point at which the auxiliary
C     gradient will be computed

      if ( gpsupn .ne. 0.0d0 ) then 
          tsmall = macheps12 * max( 1.0d0, xsupn / gpsupn )
      else
          tsmall = 0.0d0
      end if

      do i = 1,n
          gpi = x(i) - g(i)
          if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
              gpi = - g(i)
          else
              gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
          end if
          s(i) = x(i) + tsmall * gpi
      end do

C     Compute the gradient at the auxiliary point

      call ievalnalu(n,s,m,lambda,rho,equatn,linear,.false.,y,inform)
      if ( inform .lt. 0 ) return

C     Compute s = x_{1/2} - x_0 and y = g_{1/2} - g_0

      sts = 0.0d0
      sty = 0.0d0
      ssupn = 0.0d0
      ysupn = 0.0d0
      yeucn = 0.0d0
      do i = 1,n
          s(i) = s(i) - x(i)
          y(i) = y(i) - g(i)
          sts  = sts + s(i) ** 2
          sty  = sty + s(i) * y(i)
          ssupn = max( ssupn, abs( s(i) ) )
          ysupn = max( ysupn, abs( y(i) ) )
          yeucn = yeucn + y(i) ** 2
      end do
      seucn = sqrt( sts )
      yeucn = sqrt( yeucn )

C     Compute a linear relation between gpsupn and cgeps, i.e.,
C     scalars a and b such that 
c
C         a * log10( ||g_P(x_ini)|| ) + b = log10(cgeps_ini) and
c
C         a * log10( ||g_P(x_fin)|| ) + b = log10(cgeps_fin),
c
C     where cgeps_ini and cgeps_fin are provided. Note that if 
C     cgeps_ini is equal to cgeps_fin then cgeps will be always 
C     equal to cgeps_ini and cgeps_fin.

      acgeps = log10( cgepsf / cgepsi ) / log10( cggpnf / gpsupn )
      bcgeps = log10( cgepsi ) - acgeps * log10( gpsupn )

C     Print initial information

      if ( iprintinn .eq. 1 ) then
          if ( mod(iter,10) .eq. 0 ) then
              write(* ,1040)
              write(file10_unit,1040)
          end if
          write(* ,1050) iter,f,gpsupn,xsupn
          write(file10_unit,1050) iter,f,gpsupn,xsupn

      else if ( iprintinn .ge. 2 ) then
          write(* ,1070) iter,ittext(ittype)
          write(file10_unit,1070) iter,ittext(ittype)

          if ( iprintinn .ge. 3 ) then
              write(* ,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
              write(file10_unit,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
          end if

          write(* ,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn
          write(file10_unit,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn

          if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
              write(* ,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(* ,1110) nprint,(x(i),i=1,nprint)
              write(* ,1130) nprint,(g(i),i=1,nprint)
              write(* ,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
              write(file10_unit,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(file10_unit,1110) nprint,(x(i),i=1,nprint)
              write(file10_unit,1130) nprint,(g(i),i=1,nprint)
              write(file10_unit,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
          end if
      end if

C     Save intermediate data for crash report

      if ( iprintctl(6) ) then
          file20_unit_loop: do file20_unit=10, 99
    
            inquire(unit=file20_unit, opened=connected)
    
            if (.not.connected) exit file20_unit_loop
    
          end do file20_unit_loop
          
          open(file20_unit,file='gencan-tabline.out')
          write(file20_unit,1400) f,0.0d0,f,0.0d0,gpsupn,f,0.0d0,f,0.0d0,gpsupn,
     +                   inform,9,n,m,0,iter,fcnt,0,0,0,0,999.9d0
          close(file20_unit)
      end if

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     Test user-provided stopping criterion

      if ( useustp ) then

          vustop = sstop(n,x,m,lambda,rho,equatn,linear,inform)
          if ( inform .lt. 0 ) return

          if ( vustop ) then
              geninfo = 8

              if ( iprintinn .ge. 1 ) then
                  write(*, 1280)
                  write(file10_unit,1280)
              end if

              return
          end if

      end if

C     Test whether the continuous-projected-gradient sup-norm
C     is small enough to declare convergence

      if ( .not. ignoref ) then
          if ( gpsupn .le. epsopt ) then
              geninfo = 0

              if ( iprintinn .ge. 1 ) then
                  write(*, 1200)
                  write(file10_unit,1200)
              end if

              return
          end if

      else
          if ( cnormu .le. epsfeas ) then
              geninfo = 0

              if ( iprintinn .ge. 1 ) then
                  write(*, 1200)
                  write(file10_unit,1200)
              end if

              return
          end if
      end if

C     Test whether the number of iterations is exhausted

      if ( iter .ge. maxit ) then
          geninfo = 1

          if ( iprintinn .ge. 1 ) then
              write(*, 1210)
              write(file10_unit,1210)
          end if

          return
      end if

C     ==================================================================
C     Test stopping criteria related to lack of progress
C     ==================================================================

      if ( iter .gt. 0 ) then

C     Test whether we have performed many iterations without 
C     moving from the current point, by checking the functional value

      if ( abs( fdiff ) .le. macheps * max( 1.0d0, abs( f ) ) ) then

          itnfp = itnfp + 1

          if ( itnfp .ge. maxinnitnp ) then

              if ( itnfp .gt. cinit ) then
                  geninfo = 3

                  if ( iprintinn .ge. 1 ) then
                      write(*, 1230)
                      write(file10_unit,1230)
                  end if

                  return

              else
                  forceoi = .true.
              end if
          end if

      else
          itnfp = 0
      end if

C     Test whether we have performed many iterations without 
C     moving from the current point, by checking their gradients

      if ( ysupn .le. macheps * max( 1.0d0, gpsupn ) .or.
     +     yeucn .le. macheps * max( 1.0d0, gpeucn ) ) then

          itngp = itngp + 1

          if ( itngp .ge. maxinnitnp ) then

              if ( itngp .gt. cinit ) then
                  geninfo = 4

                  if ( iprintinn .ge. 1 ) then
                      write(*, 1240)
                      write(file10_unit,1240)
                  end if

                  return

              else
                  forceoi = .true.
              end if
          end if

      else
          itngp = 0
      end if

C     Test whether we have performed many iterations without 
C     moving from the current point, by checking the step norm

      if ( ssupn .le. macheps * max( 1.0d0, xsupn ) .or.
     +     seucn .le. macheps * max( 1.0d0, xeucn ) ) then

          itnxp = itnxp + 1

          if ( itnxp .ge. maxinnitnp ) then

              if ( itnxp .gt. cinit ) then

                  geninfo = 5

                  if ( iprintinn .ge. 1 ) then
                      write(*, 1250)
                      write(file10_unit,1250)
                  end if

                  return

              else
                  forceoi = .true.
              end if
          end if

      else
          itnxp = 0
      end if

      end if

C     ==================================================================
C     Iteration
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     Compute new iterate
C     ==================================================================

C     We abandon the current face if the norm of the internal gradient
C     (here, internal components of the continuous projected gradient)
C     is smaller than eta times the norm of the continuous 
C     projected gradient. Using eta = 0.1 is a rather conservative 
C     strategy in the sense that internal iterations are preferred over 
C     SPG iterations. Replace eta = 0.1 by other tolerance in (0,1) if 
C     you find it convenient. 

      if ( gieucn .le. eta * gpeucn .or. 
     +     gisupn .le. eta * gpsupn .or. forceoi ) then

C         ==============================================================
C         Some constraints should be abandoned. Compute the new iterate 
C         using an SPG iteration
C         ==============================================================

          forceoi = .false.
          lvfit = lvfit + 1
          cinit = 0

C         Compute spectral steplength

          if ( sty .le. 0.0d0 ) then
              lamspg = max( 1.0d0, xsupn / gpsupn )
          else
              lamspg = sts / sty
          end if
          lamspg = min( lspgma, max( lspgmi, lamspg ) )

C         Perform safeguarded quadratic interpolation along the 
C         spectral continuous projected gradient

          call spgls(n,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg,
     +    xplus,fplus,tmp,d,sevalal,lsinfo,inform) 
          if ( inform .lt. 0 ) return

          call sevalnal(n,xplus,m,lambda,rho,equatn,linear,gplus,inform)
          if ( inform .lt. 0 ) return

C         Set iteration type

          ittype = 1

      else

C         ==============================================================
C         The new iterate will belong to the closure of the current face
C         ==============================================================

          innit = innit + 1
          cinit = cinit + 1

C         Shrink the point, its gradient and the bounds

          call shrink(nind,x)
          call shrink(nind,g)
          call shrink(nind,l)
          call shrink(nind,u)

C         Save values of fixed variables for further evaluations

          do i = 1,n - nind
              xcomplement(i) = x(nind + i)
          end do

          nt = n

C         Compute the new point using a trust-region approach

          if ( innittype .eq. 'TR' ) then

              call betra(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,
     +        trdel,newtrdel,mslamb,epsopt,xeucn,xsupn,gpsupn,d,rbdind,
     +        rbdtype,xplus,fplus,gplus,trcnt,iscnt,chcnt,memfail,
     +        trinfo,inform)

              if ( inform .lt. 0 ) return

              if ( memfail ) then
                  mfit = iter
              end if
          
          end if

C         Compute the descent direction solving the newtonian system

          if ( innittype .eq. 'DS' ) then

C             Solve the Newtonian system
C
C                 ( H + rho A^T A ) x = b
C
C             by solving the Martinez-Santos system
C
C                 H x + A^t y = b        
C                 A x - y/rho = 0
C
C             with a direct solver.

              call newtd(nind,x,l,u,g,m,rho,equatn,d,memfail,inform)
              if ( inform .lt. 0 ) return

              nwcnt = nwcnt + 1

              if ( memfail ) then
                  mfit = iter
              end if

C             Set iteration type

              ittype = 3

          end if

          if ( innittype .eq. 'CG' .or. 
     +       ( innittype .eq. 'DS' .and. memfail ) .or.
     +       ( innittype .eq. 'TR' .and. memfail ) .or.
     +       ( innittype .eq. 'TR' .and. trinfo .eq. UNDEFSTEP ) ) then

C             Compute "trust-region radius"

              if ( iter .eq. 1 ) then
                  cgdel = max( 100.0d0, 100.0d0 * xsupn ) 
              else
                  cgdel = max( delmin, 10.0d0 * ssupn )
              end if

C             Set conjugate gradient stopping criteria. 

              cgmaxit = min( 2 * nind, 10000 )
              cgeps   = 10.0d0 ** ( acgeps * log10( gpsupn ) + bcgeps )
              cgeps   = max( cgepsf, min( cgepsi, cgeps ) )

C             Call Conjugate Gradients to solve the Newtonian system

              call cgm(nind,x,m,lambda,rho,equatn,linear,l,u,g,cgdel,
     +        cgeps,cgmaxit,precond,d,rbdnnz,rbdind,rbdtype,cgiter,
     +        cginfo,inform)

              cgcnt = cgcnt + cgiter

              if ( inform .lt. 0 ) return

C             Set iteration type

              ittype = 2

          end if

          if ( ittype .le. 3 ) then

C             Compute maximum step

              if ( ittype .eq. 2 .and. cginfo .eq. BDSREACHED ) then
                  amax = 1.0d0
                  
              else
                  call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)
              end if

C             Perform line search

              call tnls(nind,x,l,u,m,lambda,rho,equatn,linear,f,g,amax,
     +        d,rbdnnz,rbdind,rbdtype,xplus,fplus,gplus,lsinfo,inform)
              if ( inform .lt. 0 ) return

          end if
              
C         Set iteration type based on the lenght of the maximum step

          dsupn = 0.0d0 
          do i = 1,nind
              dsupn = max( dsupn, abs( d(i) ) )
          end do
          
          if ( avoidds .or. 
     +         ( memfail .and. iter - mfit .le. 3 ) ) then
              innittype = 'CG'
              
          else
              if ( innittype .eq. 'DS' .and. 
     +             amax .le. macheps12 * dsupn ) then
                  innittype = 'CG'
              elseif ( cisit .ge. 3 ) then
                  innittype = 'CG'
              else
                  innittype = 'DS'
              end if
          end if
          
C         Expand the point, its gradient and the bounds. Also expand 
C         xplus and gplus.

          call expand(nind,x)
          call expand(nind,g)
          call expand(nind,l)
          call expand(nind,u)

          call expand(nind,xplus)
          call expand(nind,gplus)
      end if

C     Compute s = xplus - x and y = gplus - g

      fdiff = fplus - f

      sts = 0.0d0
      sty = 0.0d0
      ssupn = 0.0d0
      ysupn = 0.0d0
      yeucn = 0.0d0
      do i = 1,n
          s(i) = xplus(i) - x(i)
          y(i) = gplus(i) - g(i)
          sts  = sts + s(i) ** 2
          sty  = sty + s(i) * y(i)
          ssupn = max( ssupn, abs( s(i) ) )
          ysupn = max( ysupn, abs( y(i) ) )
          yeucn = yeucn + y(i) ** 2
      end do
      seucn = sqrt( sts )
      yeucn = sqrt( yeucn )

C     Set new point

      f = fplus

      do i = 1,n
          x(i) = xplus(i)
          g(i) = gplus(i)
      end do

C     ==================================================================
C     Prepare for the next iteration 
C     ==================================================================

C     This adjustment/projection is ''por lo que las putas pudiera''

c      do i = 1,n
c          if ( x(i) - macheps23 * max( 1.0d0, abs( x(i) ) ) .le. 
c     +         l(i) ) then
c              x(i) = l(i)
c          else if ( x(i) + macheps23 * max( 1.0d0, abs( x(i) ) ) .ge. 
c     +              u(i) ) then
c              x(i) = u(i)
c          end if
c      end do

C     Compute continuous-project-gradient Euclidian and Sup norms,
C     internal gradient Euclidian norm, and store in nind the number of
C     free variables and in array ind their identifiers. Verify whether 
C     x and xprev belong to the same face.

      sameface = .true.
      nindprev = nind

      nind   = 0
      xsupn  = 0.0d0
      xeucn  = 0.0d0
      gpsupn = 0.0d0
      gpeucn = 0.0d0
      gisupn = 0.0d0
      gieucn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind = i
          end if
          xeucn = xeucn + x(i) ** 2
          gpi = x(i) - g(i)
          if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
              gpi = - g(i)
          else
              gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
          end if
          gpsupn = max( gpsupn, abs( gpi ) )
          gpeucn = gpeucn + gpi ** 2
          if ( x(i) - l(i) .gt. macheps23 * max( 1.0d0, l(i) ) .and. 
     +         u(i) - x(i) .gt. macheps23 * max( 1.0d0, u(i) ) ) then
              gisupn    = max( gisupn, abs( gpi ) )
              gieucn    = gieucn + gpi ** 2
              nind      = nind + 1
              if ( nind .gt. nindprev .or. ind(nind) .ne. i ) then
                  sameface = .false.
              end if
              ind(nind) = i
          end if
      end do
      xeucn  = sqrt( xeucn  )
      gpeucn = sqrt( gpeucn )
      gieucn = sqrt( gieucn )

      if ( ittype .eq. 5 .and. sameface ) then
          cisit = cisit + 1
      else
          cisit = 0
      end if

C     To be used in feasibility problems, compute constraints norm

      call sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)
      if ( inform .lt. 0 ) return

C     Print information of this iteration

      if ( iprintinn .eq. 1 ) then
          if ( mod(iter,10) .eq. 0 ) then
              write(* ,1040)
              write(file10_unit,1040)
          end if
          write(* ,1060) iter,f,gpsupn,xsupn,ssupn,nind,lvfit,innit,fcnt
          write(file10_unit,1060) iter,f,gpsupn,xsupn,ssupn,nind,lvfit,innit,fcnt

      else if ( iprintinn .ge. 2 ) then
          write(* ,1070) iter,ittext(ittype)
          write(file10_unit,1070) iter,ittext(ittype)

          if ( iprintinn .ge. 3 ) then
              write(* ,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
              write(file10_unit,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
          end if

          write(* ,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn
          write(file10_unit,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn

          if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
              write(* ,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(* ,1110) nprint,(x(i),i=1,nprint)
              write(* ,1130) nprint,(g(i),i=1,nprint)
              write(* ,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
              write(file10_unit,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(file10_unit,1110) nprint,(x(i),i=1,nprint)
              write(file10_unit,1130) nprint,(g(i),i=1,nprint)
              write(file10_unit,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
          end if
      end if

C     Save intermediate data for crash report

      if ( iprintctl(6) ) then
          open(file20_unit,file='gencan-tabline.out')
          write(file20_unit,1400) f,0.0d0,f,0.0d0,gpsupn,f,0.0d0,f,0.0d0,gpsupn,
     +                   inform,9,n,m,0,iter,fcnt,0,0,0,0,999.9d0
          close(file20_unit)
      end if

C     ==================================================================
C     Test line-search and trust-region related stopping criteria
C     ==================================================================

      if ( ittype .ge. 4 ) then
          
C         Test whether the functional value is unbounded

          if ( trinfo .eq. UNBOUNDED ) then
              geninfo = 6
              
              if ( iprintinn .ge. 1 ) then
                  write(*, 1260)
                  write(file10_unit,1260)
              end if
              
              return
          end if

C         Test whether the line search did a very small step
          
          if ( trinfo .eq. SMALLSTEP ) then
              forceoi = .true.
          end if

C         Test whether the trust-region radius became too small
          
          if ( trinfo .eq. SMALLTRRAD ) then
              forceoi = .true.
          end if

C         Test whether the point is stationary in the face
          
          if ( trinfo .eq. FSTORDPNT ) then
              forceoi = .true.
          end if

          if ( trinfo .eq. SECORDPNT ) then
              forceoi = .true.
          end if

      else

C         Test whether the functional value is unbounded

          if ( lsinfo .eq. UNBOUNDED ) then
              geninfo = 6
              
              if ( iprintinn .ge. 1 ) then
                  write(*, 1260)
                  write(file10_unit,1260)
              end if
              
              return
          end if
          
C         Test whether the line search did a very small step
          
          if ( lsinfo .eq. SMALLSTEP ) then
              if ( ittype .eq. 1 ) then
                  geninfo = 7
                  
                  if ( iprintinn .ge. 1 ) then
                      write(*, 1270)
                      write(file10_unit,1270)
                  end if
                  
                  return
              else
                  forceoi = .true.
              end if
          end if

      end if

C     ==================================================================
C     Iterate 
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 1000 format(/,5X,'Entry to GENCAN.',
     +       /,5X,'Number of variables: ',I7)

 1010 format(/,5X,'Lower bounds  (first ',I7,' components):',
     +       /,(5X, 6(1X,1P,D11.4)))
 1020 format(/,5X,'Upper bounds  (first ',I7,' components):',
     +       /,(5X, 6(1X,1P,D11.4)))
 1030 format(/,5X,'Initial point (first ',I7,' components):',
     +       /,(5X, 6(1X,1P,D11.4)))

 1040 format(/,5X,'It',16X,'f',2X,'gpsupn',3X,'xsupn',3X,'ssupn',4X,
     +            'nind',3X,'lvfit',3X,'innit',4X,'fcnt')
 1050 format(     I7,1X,1P,D16.8,2(1X,1P,D7.1))
 1060 format(     I7,1X,1P,D16.8,3(1X,1P,D7.1),4(1X,I7))

 1070 format(/,5X,'GENCAN ITERATION                  ',12X,'= ',I7,1X,
     +         A19)
 1080 format(/,5X,'Leaving-face iterations           ',12X,'= ',I7,
     +       /,5X,'Inner-to-the-face iterations      ',12X,'= ',I7,
     +       /,5X,'Functional evaluations            ',12X,'= ',I7,
     +       /,5X,'Conjugate gradient iterations     ',12X,'= ',I7,
     +       /,5X,'Newtonian system factorizations   ',12X,'= ',I7,
     +       /,5X,'Trust-region iterations           ',12X,'= ',I7,
     +       /,5X,'Inner SPG iterations              ',12X,'= ',I7,
     +       /,5X,'Trust-region matrix factorizations',12X,'= ',I7)

 1090 format(/,5X,'Functional value                              = ',
     +             1P,D24.16,
     +       /,5X,'Sup-norm of the continuous projected gradient = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1,
     +       /,5X,'Sup-norm of the internal projection of gp     = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1,
     +       /,5X,'Sup-norm of x (attained at x_{', I7,'})       = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1,
     +       /,5X,'Sup-norm of x - x_{prev}                      = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1)
 1100 format(/,5X,'Current free variables (first ',I7,', total ',
     +            'number ',I7,'): ',
     +       /,(5X, 6(1X,I7)))
 1110 format(/,5X,'Current point (first ',I7, ' components): ',
     +       /,(5X, 6(1X,1P,D11.4)))
 1130 format(/,5X,'Current gradient (first ',I7,' components): ',
     +       /,(5X, 6(1X,1P,D11.4)))
 1140 format(/,5X,'Current continuous projected gradient (first ',I7, 
     +         5X,'components): ',
     +       /,(5X, 6(1X,1P,D11.4)))

 1200 format(/,5X,'Flag of GENCAN: Solution was found.',/)
 1210 format(/,5X,'Flag of GENCAN: Maximum of iterations reached.',/)
 1230 format(/,5X,'Flag of GENCAN: Lack of progress in the ',
     +            'functional value.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1240 format(/,5X,'Flag of GENCAN: Lack of progress in the gradient ',
     +            'norm.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1250 format(/,5X,'Flag of GENCAN: Lack of progress in the current ',
     +            'point.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1260 format(/,5X,'Flag of GENCAN: Objective function seems to be ',
     +            'unbounded.',/)
 1270 format(/,5X,'Flag of GENCAN: Too small step in the line search.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1280 format(/,5X,'Flag of GENCAN: User-provided stopping criterion ',
     +            'satisfied.',/)

 1400 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,I3,1X,I1,1X,I6,1X,I6,1X,I2,1X,I7,1X,I7,1X,I2,1X,I7,
     +       1X,I7,1X,I7,0P,F8.2)

      end


C     ******************************************************************
C     ******************************************************************

      subroutine compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)

      implicit none

C     SCALAR ARGUMENTS
      integer nind,rbdnnz
      double precision amax

C     ARRAY ARGUMENTS
      integer rbdind(nind)
      character rbdtype(nind)
      double precision d(nind),l(nind),u(nind),x(nind)

C     Compute maximum step amax > 0 such that l <= x + amax d <= u.

#include "dim.par"
#include "machconst.com"

C     LOCAL SCALARS
      integer i
      double precision amaxi

      rbdnnz = 0

      amax = bignum

      do i = 1,nind
         if ( d(i) .gt. 0.0d0 ) then

             amaxi = ( u(i) - x(i) ) / d(i)

             if ( amaxi .lt. amax ) then
                 amax       = amaxi
                 rbdnnz     = 1
                 rbdind(1)  = i
                 rbdtype(1) = 'U'
             else if ( amaxi .eq. amax ) then
                 rbdnnz          = rbdnnz + 1
                 rbdind(rbdnnz)  = i
                 rbdtype(rbdnnz) = 'U'
             end if

         else if ( d(i) .lt. 0.0d0 ) then

             amaxi = ( l(i) - x(i) ) / d(i)

             if ( amaxi .lt. amax ) then
                 amax      = amaxi
                 rbdnnz     = 1
                 rbdind(1)  = i
                 rbdtype(1) = 'L'
             else if ( amaxi .eq. amax ) then
                 rbdnnz          = rbdnnz + 1
                 rbdind(rbdnnz)  = i
                 rbdtype(rbdnnz) = 'L'
             end if

         end if
      end do

      end

