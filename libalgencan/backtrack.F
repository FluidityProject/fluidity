      subroutine backtracking(dim,x,m,lambda,rho,equatn,linear,f,xsupn,
     +d,gtd,dsupn,alpha,fp,xp,evalaldim,btinfo,inform)

C     SCALAR ARGUMENTS
      integer btinfo,dim,inform,m
      double precision alpha,dsupn,f,fp,gtd,xsupn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision d(dim),lambda(m),rho(m),x(dim),xp(*)

C     SUBROUTINE ARGUMENTS
      external evalaldim

C     Backtracking with quadratic interpolation.
C
C     btinfo:
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
      logical smallstep
      integer i,interp
      double precision atmp

      interp = 0

 2010 continue

C     Test Armijo condition

      if ( fp .le. f + alpha * gamma * gtd ) then

C         Finish backtracking with the current point
   
          btinfo = 0

          if ( iprintinn .ge. 6 ) then
              write(*, 900)
              write(file10_unit,900)
          end if

          return

      end if

C     Test if we obtained a functional value similar to the current one
C     associated to a very small step
      
      if ( alpha * dsupn .le. macheps23 * xsupn .and.
     +     fp - f        .le. macheps23 * abs(f) ) then

C         Finish backtracking with the current point
          
          btinfo = 1
          
          if ( iprintinn .ge. 6 ) then
              write(*, 910)
              write(file10_unit,910)
          end if
          
          return
          
      end if

C     Test f going to -inf

      if ( fp .le. fmin ) then

C         Finish backtracking with the current point

          btinfo = 2

          if ( iprintinn .ge. 6 ) then
              write(*, 920)
              write(file10_unit,920)
          end if

          return

      end if

C     Compute new step

      interp = interp + 1

      atmp = ( - gtd * alpha ** 2 ) / 
     +       ( 2.0d0 * ( fp - f - alpha * gtd ) )

      if ( atmp .ge. sigma1 * alpha .and. 
     +     atmp .le. sigma2 * alpha ) then
          alpha = atmp
      else
          alpha = alpha / etaint
      end if

C     Compute new trial point

      do i = 1,dim
          xp(i) = x(i) + alpha * d(i)
      end do

      call evalaldim(dim,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

C     Print information of this iteration

      if ( iprintinn .ge. 6 ) then
          write(*, 110) alpha,fp,fcnt
          write(file10_unit,110) alpha,fp,fcnt
      end if

C     Test whether the step is too small

      smallstep = .true.
      do i = 1,dim
          if ( alpha * abs( d(i) ) .gt. 
     +         macheps * max( 1.0d0, abs( x(i) ) ) ) then
              smallstep = .false.
          end if
      end do

      if ( interp .ge. mininterp .and. smallstep .and. 
     +     fp - f .le. macheps23 * abs(f) ) then

C         Finish backtracking with the current point

          btinfo = 3

          if ( iprintinn .ge. 6 ) then
              write(*, 930)
              write(file10_unit,930)
          end if
  
          return

      end if

      go to 2010

C     NON-EXECUTABLE STATEMENTS

 110  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 900  format(  5X,'Flag of backtracking: Armijo condition holds.')
 910  format(  5X,'Flag of backtracking: Small step with similar ',
     +            'functional value.')
 920  format(  5X,'Flag of backtracking: Unbounded objective function?')
 930  format(  5X,'Flag of backtracking: Too small backtracking step.')

      end
