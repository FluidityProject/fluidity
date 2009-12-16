C     *****************************************************************
C     *****************************************************************

      subroutine betra(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,
     +trdelta,newdelta,mslamb,epsg,xeucn,xsupn,gpsupn,d,rbdind,rbdtype,
     +xtrial,ftrial,gtrial,triter,ispgiter,chcnt,memfail,betinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical memfail
      integer betinfo,chcnt,inform,ispgiter,m,n,nind,triter
      double precision epsg,f,ftrial,gpsupn,mslamb,newdelta,trdelta,
     +        xeucn,xsupn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      character rbdtype(nind)
      integer rbdind(nind)
      double precision d(nind),g(n),gtrial(n),l(n),lambda(m),rho(m),
     +        u(n),x(n),xtrial(n)

C     Solves the "unconstrained" inner problem that arises from the 
C     active-set method to solve box-constrained minimization problem
C
C            Minimize f(x) subject to l <= x <= u
C     
C     described in 
C
C     M. Andretta, E. G. Birgin e J. M. Mart�nez. "Practical active-set 
C     Euclidian trust-region method with spectral projected gradients 
C     for bound-constrained minimization". Optimization 54, pp. 
C     305-325, 2005.

C     betinfo:
C
C     0: Sufficient decrease of the function
C     1: x hits the boundary
C     2: Unbounded objective function?
C     3: Either the search direction or the step lenght is too small
C     4: Trust-region radius too small
C     5: Undefined direction
C     6: x is a first-order stationary point close to the boundary
C     7: x is a second-order stationary point

#include "dim.par"
#include "algconst.par"
#include "counters.com"
#include "hessdat.com"
#include "itetyp.com"
#include "machconst.com"
#include "outtyp.com"
#include "sydat.com"

C     LOCAL SCALARS
      logical frstit,pd,samep
      integer dinfo,extinfo,i,col,lin,index,lsinfo,rbdnnz
      double precision amax,ared,dbound,deucn,gsupn,lamspg,phi,pred,
     +        stpl,tmp

C     LOCAL ARRAYS
      integer hdiag(nmax)

C     EXTERNAL SUBROUTINES
      external calcal

C     Print presentation information

      if ( iprintinn .ge. 5 ) then
          write(*, 1000)
          write(file10_unit,1000)
      end if
      
C     Initialization

      if ( .not. sameface ) then
          trdelta = max( trdelta, newdelta )
      end if
      
      newdelta = 0.0d0
      
      memfail  = .false.

      frstit   = .true.

      ittype   = 4

C     ==================================================================
C     Compute distance to the boundary
C     ==================================================================
      
C     step 1: calculate the distance between x and the boundary.
C             dbound is the largest positive delta such that the ball 
C             centered at x with radius delta is still inside the box 
C             (set of constraints).

      dbound = bignum

      do i = 1,nind
          dbound = min( dbound, x(i) - l(i) )
      end do
      
      do i = 1,nind
          dbound = min( dbound, u(i) - x(i) )
      end do

C     Calculate infinite-norm of gradient.

      gsupn = 0.0d0
      do i = 1,nind
          gsupn = max( gsupn, abs( g(i) ) )
      end do

C     ==================================================================
C     Close to the boundary: perform inner SPG iteration
C     ==================================================================
      
C     step 2: close to the boundary: if the gradient is null then 
C             algorithm stops. Else, perform an inner SPG iteration.

      if ( dbound .lt. 2.0d0 * trdelmin ) then

          if ( gsupn .le. epsg ) then

              do i = 1,n
                  xtrial(i) = x(i)
                  gtrial(i) = g(i)
              end do
              
              ftrial  = f

              betinfo = 6
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9010) 
                  write(file10_unit,9010) 
              end if
              
              return
          end if

C         Inner SPG iteration

          if ( iprintinn .ge. 5 ) then
              write(*, 2070)
              write(file10_unit,2070)
          end if

          ittype = 5

C         Compute spectral steplength
          
          if ( sty .le. 0.0d0 ) then
              lamspg = max( 1.0d0, xsupn / gpsupn )
          else
              lamspg = sts / sty
          end if
          lamspg = min( lspgma, max( lspgmi, lamspg ) )
          
          ispgiter = ispgiter + 1

          call spgls(nind,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg,
     +    xtrial,ftrial,stpl,d,calcal,lsinfo,inform)
          
          if ( inform .lt. 0 ) return

          if ( lsinfo .eq. 0 .or. lsinfo .eq. 1 ) then
              betinfo = 0
          elseif ( lsinfo .eq. 2 ) then
              betinfo = 2
          elseif ( lsinfo .eq. 3 ) then
              betinfo = 3
          end if

          if ( betinfo .eq. 0 ) then

              if ( stpl .ge. 1.0d0 .and. extrp2 ) then

C                 Calculate amax, which is the biggest positive scalar 
C                 such that x + amax*d is still in the box bounded by l 
C                 and u.

                  call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)

                  call extrapolation(nind,x,l,u,m,lambda,rho,equatn,
     +            linear,g,xtrial,ftrial,gtrial,d,stpl,amax,rbdnnz,
     +            rbdind,rbdtype,fmin,beta,etaext,maxextrap,extinfo,
     +            inform)

                  if ( inform .lt. 0 ) return

                  if ( extinfo .eq. 2 ) then
                      betinfo = 2
                  else
                      betinfo = 0
                  end if

                  if ( betinfo .eq. 2 ) then
                      if ( iprintinn .ge. 5 ) then
                          write(*, 9050) 
                          write(file10_unit,9050) 
                      end if
                      
                      return
                  else
                      if ( stpl .ge. amax ) then
                          
                          betinfo = 1
                          
                          if ( iprintinn .ge. 5 ) then
                              write(*, 9030) 
                              write(file10_unit,9030) 
                          end if
                          
                      else
                          
                          betinfo = 0
                          
                          if ( iprintinn .ge. 5 ) then
                              write(*, 9000) 
                              write(file10_unit,9000) 
                          end if
                      end if

                      return
                  end if
              end if
          end if
          
C         Compute the gradient at the new iterate
          
          call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,gtrial,
     +    inform)
          if ( inform .lt. 0 ) return

          if ( betinfo .eq. 2 ) then

              if ( iprintinn .ge. 5 ) then
                  write(*, 9050) 
                  write(file10_unit,9050) 
              end if
              
              return
              
          elseif ( betinfo .eq. 3 ) then 

              if ( iprintinn .ge. 5 ) then
                  write(*, 9060) 
                  write(file10_unit,9060) 
              end if
              
              return
              
          else
              
C             If new point x is in the boundary, the inner algorithm 
C             stops and returns x as solution.
              
              rbdnnz = 0
      
              do i = 1,nind
                  if ( xtrial(i) .le. l(i) + 
     +                 macheps23 * max( abs( l(i) ), 1.0d0 ) ) then
                      xtrial(i) = l(i)
                      rbdnnz = rbdnnz + 1
                  else if ( xtrial(i) .ge. u(i) - macheps23 * 
     +                     max( abs( u(i) ), 1.0d0 ) ) then
                      xtrial(i) = u(i)
                      rbdnnz = rbdnnz + 1
                  end if
              end do
              
              if ( rbdnnz .gt. 0 ) then
                  
                  betinfo = 1
                  
                  if ( iprintinn .ge. 5 ) then
                      write(*, 9030) 
                      write(file10_unit,9030) 
                  end if

              else
                  
                  betinfo = 0
                  
                  if ( iprintinn .ge. 5 ) then
                      write(*, 9000) 
                      write(file10_unit,9000) 
                  end if

              end if

              return

          end if
          
      end if

C     ==================================================================
C     Far from the boundary: perform trust-region iteration
C     ==================================================================
      
      triter = triter + 1

C     step 3: far from the boundary, solve trust-region subproblem. 
C             Evaluate function Hessian at x.

      call calchal(nind,x,m,lambda,rho,equatn,linear,hlin,hcol,hval,
     +hnnz,inform)
      if ( inform .lt. 0 ) return
      
      do i = 1,nind
          hdiag(i) = 0
      end do

      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)

          if ( lin .eq. col ) then
              if ( hdiag(lin) .eq. 0 ) then
                  hdiag(lin) = i
              else
                  hval(hdiag(lin)) = hval(hdiag(lin)) + hval(i)
                  hval(i) = 0.0d0
              end if 
          end if
      end do

      do i = 1,nind
         if ( hdiag(i) .eq. 0 ) then
            hnnz       = hnnz + 1            
            hlin(hnnz) = i
            hcol(hnnz) = i
            hval(hnnz) = 0.0d0
            hdiag(i)   = hnnz
         end if
      end do

C     step 4: solve the trust-region subproblem using More-Sorensen's 
C             algorithm to minimize "exactly" quadratics subjected to 
C             balls. 

C     If trust-region radius is too small, the inner algorithm stops.

 100  continue

      if ( trdelta .lt. macheps * max( 1.0d0, xeucn ) ) then

          do i = 1,n
              xtrial(i) = x(i)
              gtrial(i) = g(i)
          end do
          
          ftrial  = f

          betinfo = 4
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9040) 
              write(file10_unit,9040) 
          end if
          
          return
      end if
      
      call moresor(nind,g,hnnz,hlin,hcol,hval,hdiag,trdelta,mssig,
     +0.0d0,mseps,msmaxit,mslamb,pd,d,chcnt,memfail,dinfo)

      if ( memfail ) then
          return
      end if

C     If maximum allowed number of MEQB iterations is achieved, another
C     direction d is calculated.

      if ( dinfo .eq. 5 ) then 
          
          if ( iprintinn .ge. 5 ) then
              write(*, 2000)
              write(file10_unit,2000)
          end if

          call dogleg(nind,g,hnnz,hlin,hcol,hval,pd,trdelta,d,dinfo)

      end if
      
C     If both internal gradient and Hessian matrix are null, subroutines 
C     MEQB and dogleg stop with dinfo = 0 and then the inner algorithm 
C     stops declaring "second-order stationary point".

      if ( dinfo .eq. 0 ) then

          do i = 1,n
              xtrial(i) = x(i)
              gtrial(i) = g(i)
          end do
          
          ftrial  = f

          betinfo = 7
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9020) 
              write(file10_unit,9020) 
          end if
          
          return
      end if
      
C     Print direction

      if ( iprintinn .ge. 5 ) then
          write(*, 1020) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
          write(file10_unit,1020) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
      end if

C     Evaluate the quadratic model of the objective function at d.

      call squad(nind,d,g,hnnz,hlin,hcol,hval,phi)

C     If the value of the quadratic model at d is 0 it means that x is a 
C     second-order stationary point. In this case, inner algorithm stops
C     declaring this.
      
      if ( ( abs( phi ) .le. phieps ) .and. ( gsupn .le. epsg ) ) then

          do i = 1,n
              xtrial(i) = x(i)
              gtrial(i) = g(i)
          end do
          
          ftrial  = f

          betinfo = 7
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9020) 
              write(file10_unit,9020) 
          end if
          
          return
      end if
      
C     Calculate predicted decrease of objective function
      
      pred = abs( phi )
      
C     Calculate d Euclidian-norm
      
      deucn = 0.0d0
      do i = 1,nind
          deucn = deucn + d(i)**2
      end do
      deucn = sqrt( deucn )
      
C     To avoid NaN and Inf directions

      if ( .not. ( deucn .le. bignum ) ) then
      
          trdelta = 2.5d-1 * trdelta

          if ( iprintinn .ge. 5 ) then
              write(*, 2060) 
              write(file10_unit,2060) 
          end if

          if ( trdelta .lt. macheps * max( 1.0d0, xeucn ) ) then

              do i = 1,n
                  xtrial(i) = x(i)
                  gtrial(i) = g(i)
              end do
              
              ftrial  = f

              betinfo = 5
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9070) 
                  write(file10_unit,9070) 
              end if
              
              return
          end if

          triter  = triter + 1
          frstit  = .false.

          go to 100

      end if

C     Calculate point xtrial = x + d.
      
      do i = 1,nind
          xtrial(i) = x(i) + d(i)
      end do

      stpl = 1.0d0

C     Verify if xtrial is inside de box. If not, interior is set to
C     false and amax is set to the biggest positive scalar such that
C     x + amax*d is inside the box.

      call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)
      
C     ==================================================================
C     Point on the boundary
C     ==================================================================      
      
C     If xtrial is not interior to the box, xtrial = x + d is replaced
C     by xtrial = x + amax*d. Now xtrial is definitely interior. Actually,
C     it is in the boundary. If the objective function decreases in
C     xtrial, the inner algorithm stops with xtrial as a solution.
C     Otherwise, a new trust-region radius trdelta is chosen (smaller than
C     dbound) and a new quadratic model minimizer is calculated (which
C     is necessarily interior because of the choice of trdelta).

      if ( amax .le. 1.0d0 ) then 
          
          if ( iprintinn .ge. 5 ) then
              write(*, 2010) 
              write(file10_unit,2010) 
          end if
          
          do i = 1,nind
              xtrial(i) = x(i) + amax * d(i)
          end do
          
          stpl = amax

C         Set x(i) to l(i) or u(i) for the indices i that got to the 
C         boundary (to prevent errors).
          
          do i = 1,rbdnnz
              index = rbdind(i)
              if ( rbdtype(i) .eq. 'L' ) then
                  xtrial(index) = l(index)
              elseif ( rbdtype(i) .eq. 'U' ) then
                  xtrial(index) = u(index)
              end if
          end do
          
          call calcal(nind,xtrial,m,lambda,rho,equatn,linear,ftrial,
     +    inform)
          if ( inform .lt. 0 ) return
          
C         Print functional value

          if ( iprintinn .ge. 5 ) then
              write(*, 1030) trdelta,ftrial,fcnt
              write(file10_unit,1030) trdelta,ftrial,fcnt
          end if

C         Test whether f is very small
      
          if ( ftrial .le. fmin ) then
              
              call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,
     +        gtrial,inform)
              if ( inform .lt. 0 ) return

              betinfo = 2
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9050) 
                  write(file10_unit,9050) 
              end if
              
              return
          end if      
          
C         If the new point x + d is too close to the previous point x, 
C         inner algorithm stops

          samep = .true.
          do i = 1,nind
              if ( abs( x(i) - xtrial(i) ) .gt. 
     +             macheps * max( abs( xtrial(i) ), 1.0d0 ) ) then
                  samep = .false.
              end if
          end do
          
          if ( samep .and. ftrial - f .le. macheps23 * abs(f) ) then
              
              betinfo = 3
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9060) 
                  write(file10_unit,9060) 
              end if
              
              return
          end if

C         Test if function value decreases at xtrial
        
          if ( ( ftrial .le. f ) .or. 
     +         ( ( deucn .le. macheps23 * xeucn ) .and.
     +         ( ftrial .le. f + macheps23 * abs( f ) ) ) ) then
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 2030) 
                  write(file10_unit,2030) 
              end if
              
              if ( extrp4 ) then
                  
                  call extrapolation(nind,x,l,u,m,lambda,rho,equatn,
     +            linear,g,xtrial,ftrial,gtrial,d,stpl,amax,rbdnnz,
     +            rbdind,rbdtype,fmin,beta,etaext,maxextrap,extinfo,
     +            inform)

                  if ( inform .lt. 0 ) return

                  if ( extinfo .eq. 2 ) then
                      betinfo = 2
                  else
                      betinfo = 0
                  end if

C                 Update the trust-region radius (which may or may not
C                 be used)

                  if ( frstit ) then
                      trdelta = 2.0d0 * max( trdelta, stpl * deucn )
                  end if
                  
                  if ( betinfo .eq. 2 ) then
                      if ( iprintinn .ge. 5 ) then
                          write(*, 9050) 
                          write(file10_unit,9050) 
                      end if
                      
                      return
                  end if
                  
                  betinfo = 1
                  
C                 Calculate actual reduction of objective function

                  ared = f - ftrial
                  
                  go to 200
                  
              else
                  
C                 Update the trust-region radius (which may or may not
C                 be used)
               
                  trdelta = 2.0d0 * trdelta
                  
C                 Calculate actual reduction of objective function

                  ared = f - ftrial
                  
C                 Compute gtrial.
                  
                  call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,
     +            gtrial,inform)
                  if ( inform .lt. 0 ) return
                  
                  betinfo = 1
                  
                  go to 200
              end if
              
          else
              tmp      = trdelmin + msrho * 
     +                   ( ( dbound / (1.0d0 + mssig ) ) - trdelmin )
              newdelta = trdelta
              trdelta  = max( trdelmin, tmp )

              triter   = triter + 1
              frstit   = .false.

              if ( iprintinn .ge. 5 ) then
                  write(*, 2040) 
                  write(file10_unit,2040) 
              end if

              go to 100
          end if
      end if

C     ==================================================================
C     Point interior to the box
C     ==================================================================
      
C     step 5: in this case xtrial is inside the box. Acceptance or
C             rejection of the trust-region subproblem solution.

      call calcal(nind,xtrial,m,lambda,rho,equatn,linear,ftrial,inform)
      if ( inform .lt. 0 ) return

C     Print functional value

      if ( iprintinn .ge. 5 ) then
          write(*, 1030) trdelta,ftrial,fcnt
          write(file10_unit,1030) trdelta,ftrial,fcnt
      end if

C     Test whether f is very small
      
      if ( ftrial .le. fmin ) then
          
          call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,gtrial,
     +    inform)
          if ( inform .lt. 0 ) return

          betinfo = 2
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9050) 
              write(file10_unit,9050) 
          end if
          
          return
      end if
      
C     If the new point x + d is too close to the previous point x, inner
C     algorithm stops

      samep = .true.
      do i = 1,nind
          if ( abs( x(i) - xtrial(i) ) .gt. 
     +         macheps * max( abs( xtrial(i) ), 1.0d0 ) ) then
              samep = .false.
          end if
      end do
      
      if ( samep .and. ftrial - f .le. macheps23 * abs(f) ) then
          
          betinfo = 3
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9060) 
              write(file10_unit,9060) 
          end if
          
          return
      end if
      
C     Calculate actual reduction of objective function
      
      ared = f - ftrial
      
C     If there is not sufficient decrease of the function, the 
C     trust-region radius is decreased and the new quadratic model 
C     minimizer will be calculated. 

      if ( iprintinn .ge. 5 ) then
          write(*, 1010) deucn,pred,ared
          write(file10_unit,1010) deucn,pred,ared
      end if

      if ( deucn .le. macheps23 * xeucn .and.
     +     ared  .le. macheps23 * abs( f ) ) then
         
          call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,gtrial,
     +    inform)
          if ( inform .lt. 0 ) return
          
          go to 200

      end if
      
      if ( ( pred .le. phieps .or. ared .ge. tralpha*pred ) .and.
     +     ( pred .gt. phieps .or. f .ge. ftrial ) ) then

C         If extrapolation at step 5 is not to be performed, point 
C         xtrial is accepted.

          if ( extrp5 ) then
          
              call extrapolation(nind,x,l,u,m,lambda,rho,equatn,linear,
     +        g,xtrial,ftrial,gtrial,d,stpl,amax,rbdnnz,rbdind,rbdtype,
     +        fmin,beta,etaext,maxextrap,extinfo,inform)
              
              if ( inform .lt. 0 ) return
              
              if ( extinfo .eq. 2 ) then
                  betinfo = 2
              else
                  betinfo = 0
              end if
              
              if ( frstit ) then
                  trdelta = max( trdelta, stpl * deucn )
              end if
              
              if ( betinfo .eq. 2 ) then
                  if ( iprintinn .ge. 5 ) then
                      write(*, 9050) 
                      write(file10_unit,9050) 
                  end if
                  
                  return
              end if
              
C             Update actual reduction of objective function
              
              ared = f - ftrial
              
          else
              
              call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,
     +        gtrial,inform)
              if ( inform .lt. 0 ) return
              
          end if

      else

          trdelta = 2.5d-1 * deucn
          triter  = triter + 1
          frstit  = .false.

          if ( iprintinn .ge. 5 ) then
              write(*, 2050) 
              write(file10_unit,2050) 
          end if
          
          go to 100

      end if
      
C     ==================================================================
C     Prepare for next call to this routine
C     ==================================================================
      
C     Update the trust-region radius (which may or may not be used). 
C     This update can only be done when the current iteration was a 
C     trust-region iteration (and not an inner SPG one).

 200  continue

      if ( ared .lt. 2.5d-1 * pred ) then
          trdelta = max( 2.5d-1 * deucn, trdelmin )
      else
          if ( ( ared .ge. 0.5d0 * pred ) .and.
     +         ( deucn .ge. trdelta - macheps23 * 
     +         max( trdelta, 1.0d0 ) ) ) then
              trdelta = max( 2.0d0 * trdelta, trdelmin )
          else
              trdelta = max( trdelta, trdelmin )
          end if
      end if
      
C     If new point x is in the boundary, the inner algorithm stops
C     and returns x as solution.

      if ( stpl .ge. amax ) then
          
          betinfo = 1
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9030) 
              write(file10_unit,9030) 
          end if
          
      else
          
          betinfo = 0
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9000) 
              write(file10_unit,9000) 
          end if
      end if
      
C     Non-executable statements
      
 1000 format(/,5X,'Trust-region iteration')
 1010 format(  5X,'deucn = ',1P,D7.1,' pred = ',1P,D11.4,
     +            ' ared = ',1P,D11.4)
 1020 format(/,5X,'Trust-region direction (first ',I7,
     +            ' components): ',/,1(5X,6(1X,1P,D11.4)))
 1030 format(/,5X,'delta = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)

 2000 format(/,5X,'Since the direction More-Sorensen calculation ',
     +       /,5X,'failed, dogleg direction will be computed.')
 2010 format(/,5X,'x+d is not interior to the box.')
 2030 format(  5X,'f(x+d) < f(x), x+d is accepted.')
 2040 format(  5X,'f(x+d) >= f(x), a new direction d will be computed.')
 2050 format(  5X,'x+d did not obtain suficcient functional reduction. '
     +       /,5X,'A new direction d will be computed.')
 2060 format(/,5X,'Direction is undefined. ',
     +            'A new direction d will be computed.')
 2070 format(/,5X,'Point close to the boundary. ',
     +            'An inner SPG iteration will be used.')

 9000 format(  5X,'Flag of TR: Sufficient decrease of function.')
 9010 format(  5X,'Flag of TR: ',
     +            'First-order stationary point close to boundary.')
 9020 format(/,5X,'Flag of TR: Second-order stationary point.')
 9030 format(  5X,'Flag of TR: Point on the boundary.')      
 9040 format(  5X,'Flag of TR: Trust-region radius too small.')
 9050 format(  5X,'Flag of TR: Unbounded objective function?')
 9060 format(  5X,'Flag of TR: Very similar consecutive points.')
 9070 format(  5X,'Flag of TR: Undefined direction.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine squad(nred,x,g,hnnz,hlin,hcol,hval,phi)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,nred
      double precision phi
      
C     ARRAY ARGUMENTS
      integer hcol(hnnz),hlin(hnnz)
      double precision g(nred),hval(hnnz),x(nred)

C     Evaluates the quadratic model phi(x) = 1/2 x^T H x + g^T x.

#include "dim.par"

C     LOCAL SCALARS
      integer col,i,lin
      
C     LOCAL ARRAYS
      double precision wd(nmax)

      do i = 1,nred
          wd(i) = 0.0d0
      end do
      
      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          
          wd(lin) = wd(lin) + hval(i) * x(col)
          if ( lin .ne. col ) then
              wd(col) = wd(col) + hval(i) * x(lin)
          end if
      end do
      
      phi = 0.0d0
      do i = 1,nred
          phi = phi + wd(i) * x(i)
      end do
      
      phi = phi * 0.5d0
      
      do i = 1,nred
          phi = phi + g(i) * x(i)
      end do
      
      end      
      
C     *****************************************************************
C     *****************************************************************

C     Algorithm that finds a unconstrained minimizer of objective 
C     function inside the box of constraints, hits the boundary 
C     (obtaining function decrease), or finds an interior point where 
C     the objective function has sufficient decrease (compared to its 
C     value at x). Extrapolation may be done.
C
C     When the current point x is "close to" the boundary, a Spectral 
C     Projected Gradient (SPG) iteration is used to calculate the new 
C     point. If this new point is at the boundary, the algorithm stops. 
C     Otherwise, a new iteration begins.
C
C     When x is "far from" the boundary, trust-region radius is 
C     determined and d is calculated using More-Sorensen algorithm to 
C     solve the trust-region subproblem (which is to find a minimizer a 
C     to a function quadratic model provided that the minimizer's 
C     Euclidian-norm is smaller than a given delta). The new point is
C     xtrial = x + d. 
C
C     If xtrial lies outside the box of constraints, it is truncated on 
C     the boundary. This new y on the boundary will be candidate to be a 
C     solution. If function value at new xtrial is smaller than function 
C     value at x, inner algorithm stops with xtrial. Otherwise, the 
C     trust-region radius is decreased so that the new solution d' to 
C     the trust-region subproblem makes x + d' be interior to the box.
C     More-Sorensen algorithm is used to calculate d' too.
C
C     If xtrial lies inside the box, sufficient decrease of objective 
C     function is tested. If it is true, xtrial is accepted as a solution
C     candidate. If xtrial in on the boundary, inner algorithm stops and 
C     if it is interior, a new iteration begins. If sufficient decrease is 
C     not obtained, trust-region radius is decreased and a new quadratic 
C     model minimizer is calculated (as in a classical trust-region 
C     algorithm for unconstrained minimization).
C
C     If the user wants, after calculating the candidate solution xtrial, 
C     extrapolation may be performed. For this, set extrpi to true, 
C     where i is the step of the inner algorithm that can call 
C     extrapolation procedure. If extrp4 is true, extrapolation will be
C     tried after xtrial hit the boundary. And if extrp5 is true,
C     extrapolation will be tried after xtrial is calculated by
C     trust-region algorithm, when xtrial is interior and provides
C     sufficient decrease of objective function.
C
C     If gradient at current point is null, inner algorithm stops 
C     declaring "first-order stationary point". If quadratic model 
C     minimum is 0, inner algorithm stops declaring "second-order 
C     stationary point".
C
C     M. Andretta, E. G. Birgin and J. M. Martinez, ''Practical active-set
C     Euclidian trust-region method with spectral projected gradients for
C     bound-constrained minimization'', Optimization 54, pp. 305-325, 2005.
C
C     On Entry
C      
C     n        integer
C              dimension of full space
C
C     nind     integer
C              dimension of reduced space
C
C     x        double precision x(n)
C              initial point, interior to the current face
C
C     l        double precision l(n)
C              lower bounds on x
C
C     u        double precision u(n)
C              upper bounds on x
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C     equatn   logical equatn(m)
C     linear   logical linear(m)
C              These five parameters are not used nor modified by 
C              BETRA and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if BETRA is 
C              being used to solve the bound-constrained subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters.
C              equatn is logical array that, for each constraint, 
C              indicates whether the constraint is an equality constraint
C              (.true.) or an inequality constraint (.false.). Finally,
C              linear is logical array that, for each constraint, 
C              indicates whether the constraint is a linear constraint
C              (.true.) or a nonlinear constraint (.false.)
C
C     f        double precision
C              objective function value at x
C
C     g        double precision g(n)
C              gradient at x
C
C     trdelta  double precision 
C              trust-region radius
C
C     newdelta double precision
C              trust-region radius is set to the maximum between newdelta 
C              and trdelta
C     
C     mslamb   double precision
C              value that More-Sorensen algorithm calculates to find 
C              the trust-region subproblem solution (MEQB)
C
C     epsg     double precision
C              allowed error for projected gradient norm
C
C     xeucn    double precision
C              x Euclidian norm
C
C     xsupn    double precision
C              x sup-norm
C
C     gpeucn   double precision
C              projected-gradient Euclidian norm
C
C     On Return
C      
C     trdelta  double precision
C              updated trut-region radius
C
C     newdelta double precision
C              when the trust-region radius trdelta is decreased so that 
C              the point x + d fit the current face, newdelta is set to 
C              the previous value of trdelta. Otherwise, it is set to 0
C
C     mslamb   double precision
C              updated value for next iteration (see entry parameter)
C
C     d        double precision d(nind)
C              direction computed such that xtrial = x + d
C
C     rbdind   integer rbdind(n)
C              indices of variables that reached their bounds
C     
C     rbdtype  character rbdtype(n)
C              if variable rbdind(i) reached its lower bound, 
C              rbdtype(i) = 'L'. If variable rbdind(i) reached its upper
C              bound, rbdtype(i) = 'U'. 
C
C     xtrial   double precision xtrial(n)
C              solution candidate, with inform as described bellow
C
C     ftrial   double precision
C              objective function value at xtrial
C
C     gtrial   double precision gtrial(n)
C              gradient at xtrial
C
C     triter   integer
C              number of trust-region iterations
C
C     ispgiter integer
C              number of inner SPG iterations
C
C     chcnt    integer
C              number of Cholesky decompositions
C
C     memfail  logical
C              true iff linear solver failed because of lack of memory
C
C     betinfo  integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C      
C              0 = Sufficient decrease of the function;
C
C              1 = x hits the boundary;
C
C              2 = Unbounded objective function?
C
C              3 = Either the search direction or the step lenght is too
C                  small;
C
C              4 = Trust-region radius too small;
C
C              5 = Undefined direction;
C
C              6 = x is a first-order stationary point close to the
C                  boundary;
C
C              7 = x is a second-order stationary point;
C
C     inform   0 = no error occurred;
C             <0 = error in some function, gradient or Hessian routine.
