C     *****************************************************************
C     *****************************************************************

      subroutine tnls(nind,x,l,u,m,lambda,rho,equatn,linear,f,g,amax,d,
     +rbdnnz,rbdind,rbdtype,xp,fp,gp,lsinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,lsinfo,m,nind,rbdnnz
      double precision amax,f,fp

C     ARRAY ARGUMENTS
      integer rbdind(nind)
      character rbdtype(nind)
      logical equatn(m),linear(m)
      double precision d(nind),g(*),gp(*),l(nind),lambda(m),rho(m),
     +        u(nind),x(*),xp(*)

C     This subroutine computes a line search in direction d.
C
C     lsinfo:
C
C     At the first trial:
C
C     5: x + amax d is at the boundary and f(x + amax d) is smaller than f

C     Extrapolation:
C
C     2: Unbounded objective function?
C     4: beta-condition holds. No extrapolation is done
C     6: Maximum number of extrapolations reached
C     7: Similar consecutive projected points
C     8: Not-well-defined objective function
C     9: Functional value increases

C     Backtracking:
C
C     0: Armijo satisfied
C     1: Small step with functional value similar to the current one
C     2: Unbounded objective function?
C     3: Too small backtracking step. Wrong gradient?

#include "dim.par"
#include "machconst.com"
#include "algconst.par"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical boundary
      integer i
      double precision alpha,dsupn,gtd,xsupn

C     EXTERNAL SUBROUTINES
      external calcal

C     ==================================================================
C     ==================================================================

C     Test Armijo condition and beta condition. Decide between accept 
C     the first trial, extrapolate or backtrack. 

C     ==================================================================
C     ==================================================================

C     ------------------------------------------------------------------
C     Compute directional derivative, dsupn and xsupn
C     ------------------------------------------------------------------

      gtd = 0.0d0
      dsupn = 0.0d0
      xsupn = 0.0d0
      do i = 1,nind
          gtd = gtd + g(i) * d(i)
          dsupn = max( dsupn, abs( d(i) ) )
          xsupn = max( xsupn, abs( x(i) ) )
      end do

      if ( iprintinn .ge. 6 ) then
          write(*, 100) xsupn,amax,dsupn
          write(file10_unit,100) xsupn,amax,dsupn
      end if

C     ------------------------------------------------------------------
C     Compute first trial (projected point)
C     ------------------------------------------------------------------

      alpha = 1.0d0

      boundary = .false.
      do i = 1,nind
          xp(i) = x(i) + d(i)

          if ( xp(i) .lt. l(i) .or. xp(i) .gt. u(i) ) then
              boundary = .true.
              xp(i) = max( l(i), min( xp(i), u(i) ) )
          end if
      end do

      if ( amax .eq. 1.0d0 ) then
          boundary = .true.

          do i = 1,rbdnnz
              if ( rbdtype(i) .eq. 'L' ) then
                  xp(rbdind(i)) = l(rbdind(i))
              else if ( rbdtype(i) .eq. 'U' ) then
                  xp(rbdind(i)) = u(rbdind(i))
              end if
          end do
      end if

      call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

      if ( .not. boundary ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 110) fp,fcnt
              write(file10_unit,110) fp,fcnt
          end if

      else
          if ( iprintinn .ge. 6 ) then
              write(*, 120) fp,fcnt
              write(file10_unit,120) fp,fcnt
          end if
      end if

C     ------------------------------------------------------------------
C     The first trial is an interior point.
C     ------------------------------------------------------------------

      if ( .not. boundary ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 140)
              write(file10_unit,140)
          end if

C         Armijo condition holds.

          if ( fp .le. f + alpha * gamma * gtd ) then

              if ( iprintinn .ge. 6 ) then
                  write(*, 150)
                  write(file10_unit,150)
              end if

              go to 1000

          end if

C         Armijo condition does not hold. We will do backtracking.  

          if ( iprintinn .ge. 6 ) then
              write(* ,180)
              write(file10_unit,180)
          end if

          go to 2000

      end if

C     ------------------------------------------------------------------
C     The first trial is at the boundary.
C     ------------------------------------------------------------------

      if ( iprintinn .ge. 6 ) then
          write(*, 190) 
          write(file10_unit,190) 
      end if

C     Function value is smaller than at the current point. We will 
C     extrapolate.

      if ( fp .lt. f ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 200) 
              write(file10_unit,200) 
          end if

          go to 1000

      end if

C     Discard the projected point and consider x + amax d

      if ( iprintinn .ge. 6 ) then
          write(*, 210) 
          write(file10_unit,210) 
      end if

      alpha = amax

      do i = 1,nind
          xp(i) = x(i) + alpha * d(i)
      end do

      do i = 1,rbdnnz
          if ( rbdtype(i) .eq. 'L' ) then
              xp(rbdind(i)) = l(rbdind(i))
          else if ( rbdtype(i) .eq. 'U' ) then
              xp(rbdind(i)) = u(rbdind(i))
          end if
      end do

      call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

      if ( iprintinn .ge. 6 ) then
          write(*, 130) alpha,fp,fcnt
          write(file10_unit,130) alpha,fp,fcnt
      end if

C     Function value is smaller than or equal to (or even just a little 
C     bit greater than) at the current point. Line search is over.

      if ( fp .lt. f + macheps23 * abs(f) ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 220) 
              write(file10_unit,220) 
          end if

          call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
          if ( inform .lt. 0 ) return

          lsinfo = 5

          if ( iprintinn .ge. 6 ) then
              write(*, 900)
              write(file10_unit,900)
          end if

          return

      end if

C     Function value is greater than at the current point. We will 
C     do backtracking.

      if ( iprintinn .ge. 6 ) then
          write(*, 230) 
          write(file10_unit,230) 
      end if

      go to 2000

C     ==================================================================
C     ==================================================================

C     Extrapolation

C     ==================================================================
C     ==================================================================

 1000 continue

      call extrapolation(nind,x,l,u,m,lambda,rho,equatn,linear,g,xp,fp,
     +gp,d,alpha,amax,rbdnnz,rbdind,rbdtype,fmin,beta,etaext,maxextrap,
     +lsinfo,inform)

      return

C     ==================================================================
C     ==================================================================

C     Backtracking

C     ==================================================================
C     ==================================================================

 2000 continue

      call backtracking(nind,x,m,lambda,rho,equatn,linear,f,xsupn,d,gtd,
     +dsupn,alpha,fp,xp,calcal,lsinfo,inform)
      if ( inform .lt. 0 ) return

      call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     End of backtracking
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 100  format(/,5X,'TN Line search (xsupn = ',1P,D7.1,', amax = ',
     +             1P,D7.1,', dsupn = ',1P,D7.1,')')
 110  format(  5X,'Unitary step    F = ',1P,D24.16,' FE = ',I7)
 120  format(  5X,'Projected point F = ',1P,D24.16,' FE = ',I7)
 130  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 140  format(  5X,'The first trial is an interior point.')
 150  format(  5X,'Armijo condition holds.')
 180  format(  5X,'Armijo condition does not hold. We will backtrack.')
 190  format(  5X,'The first trial is at the boundary.')
 200  format(  5X,'Function value at the boundary is smaller than at ',
     +            'the current point.',/,5X,'We will extrapolate.')
 210  format(  5X,'Discarding projected point. We will now consider x ',
     +            '+ amax d.')
 220  format(  5X,'Function value at the boundary is smaller than or ',
     +            'equal to than at the',/,5X,'current point. Line ',
     +            'search is over.')
 230  format(  5X,'Function value at the boundary is greater than at ',
     +            'the current point.',/,5X,'We will backtrack.')
 900  format(  5X,'Flag of TN Line search: First trial accepted.')

      end
