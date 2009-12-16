C     ******************************************************************
C     ******************************************************************

      subroutine spgls(n,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg,
     +xp,fp,alpha,d,evalaldim,lsinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,lsinfo,m,n
      double precision alpha,f,fp,lamspg

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision d(n),g(n),l(n),lambda(m),rho(m),u(n),x(n),xp(n)

C     SUBROUTINE ARGUMENTS
      external evalaldim

C     This subroutine computes a line search in the Spectral Projected 
C     Gradient direction.
C
C     lsinfo:
C
C     0: Armijo satisfied
C     1: Small step with functional value similar to the current one
C     2: Unbounded objective function?
C     3: Too small backtracking step. Wrong gradient?

#include "dim.par"
#include "outtyp.com"
#include "counters.com"

C     LOCAL SCALARS
      integer i
      double precision dsupn,gtd,xsupn

C     ------------------------------------------------------------------
C     Compute search direction, directional derivative, dsupn, xsupn and 
C     first trial
C     ------------------------------------------------------------------

      gtd   = 0.0d0
      dsupn = 0.0d0
      xsupn = 0.0d0
      do i = 1,n
          d(i)  = - lamspg * g(i)
          xp(i) = x(i) + d(i)
          if ( xp(i) .lt. l(i) .or. xp(i) .gt. u(i) ) then
              xp(i) = max( l(i), min( xp(i), u(i) ) )
              d(i)  = xp(i) - x(i)
          end if
          gtd   = gtd + g(i) * d(i)
          dsupn = max( dsupn, abs( d(i) ) )
          xsupn = max( xsupn, abs( x(i) ) )
      end do

      if ( iprintinn .ge. 6 ) then
          write(* ,100) xsupn,lamspg,dsupn
          write(file10_unit,100) xsupn,lamspg,dsupn
      end if
      
      call evalaldim(n,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

      alpha = 1.0d0
      
      if ( iprintinn .ge. 6 ) then
          write(*, 110) alpha,fp,fcnt
          write(file10_unit,110) alpha,fp,fcnt
      end if
      
C     ==================================================================
C     Backtracking
C     ==================================================================

      call backtracking(n,x,m,lambda,rho,equatn,linear,f,xsupn,d,gtd,
     +dsupn,alpha,fp,xp,evalaldim,lsinfo,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     End of backtracking
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 100  format(/,5X,'SPG Line search (xsupn = ',1P,D7.1,1X,'SPGstep= ',
     +             1P,D7.1,1X,'dsupn = ',1P,D7.1,')')
 110  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)

      end
