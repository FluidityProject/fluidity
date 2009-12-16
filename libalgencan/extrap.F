C     *****************************************************************
C     *****************************************************************

      subroutine extrapolation(nind,x,l,u,m,lambda,rho,equatn,linear,g,
     +xp,fp,gp,d,alpha,amax,rbdnnz,rbdind,rbdtype,fmin,beta,etaext,
     +maxextrap,extinfo,inform)

C     SCALAR ARGUMENTS
      integer extinfo,inform,m,maxextrap,nind,rbdnnz
      double precision alpha,amax,beta,etaext,fmin,fp

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      character rbdtype(nind)
      integer rbdind(nind)
      double precision d(nind),g(*),gp(*),l(nind),lambda(m),rho(m),
     +        u(nind),x(*),xp(*)
      
C     Performs extrapolation.
C
C     extinfo:
C
C     0: Success
C     2: Unbounded objective function?
C     4: beta-condition holds. No extrapolation is done
C     6: Maximum number of extrapolations reached
C     7: Similar consecutive projected points
C     8: Not-well-defined objective function
C     9: Functional value increases

#include "dim.par"
#include "counters.com"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical projected,samep
      integer extrap,i
      double precision atmp,fbext,ftmp,gptd,gtd

C     LOCAL ARRAYS
      double precision xbext(nmax),xtmp(nmax)

      extinfo = 0

      extrap  = 0

C     Compute directional derivative

      gtd = 0.0d0
      do i = 1,nind
          gtd = gtd + g(i) * d(i)
      end do
      
      call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
      if ( inform .lt. 0 ) return
      
      gptd = 0.0d0
      do i = 1,nind
          gptd = gptd + gp(i) * d(i)
      end do

C     Beta condition holds. No extrapolation is done.

      if ( gptd .ge. beta * gtd ) then
          
          if ( iprintinn .ge. 6 ) then
              write(*, 110)
              write(file10_unit,110)
          end if
          
          extinfo = 4
          
          return
          
      end if

C     Beta condition does not holds. We will extrapolate.

      if ( iprintinn .ge. 6 ) then
          write(* ,120)
          write(file10_unit,120)
      end if
      
C     Save f and x before extrapolation to return in case of a
C     not-well-defined objective function at an extrapolated point

      fbext = fp

      do i = 1,nind
          xbext(i) = xp(i)
      end do

 1010 continue

C     Test f going to -inf

      if ( fp .le. fmin ) then

C         Finish the extrapolation with the current point

          if ( extrap .ne. 0 ) then
              call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
              if ( inform .lt. 0 ) return
          end if

          extinfo = 2

          if ( iprintinn .ge. 6 ) then
              write(*, 910)
              write(file10_unit,910)
          end if

          return

      end if

C     Test if the maximum number of extrapolations was exceeded

      if ( extrap .ge. maxextrap ) then

C         Finish the extrapolation with the current point

          if ( extrap .ne. 0 ) then
              call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
              if ( inform .lt. 0 ) return
          end if

          extinfo = 6

          if ( iprintinn .ge. 6 ) then
              write(*, 930)
              write(file10_unit,930)
          end if

          return

      end if

C     Chose new step 

      extrap = extrap + 1

      if ( alpha .lt. amax .and. etaext * alpha .gt. amax ) then
          atmp = amax
      else
          atmp = etaext * alpha
      end if

C     Compute new trial point

      do i = 1,nind
          xtmp(i) = x(i) + atmp * d(i)
      end do

      if ( atmp .eq. amax ) then
          do i = 1,rbdnnz
              if ( rbdtype(i) .eq. 'L' ) then
                  xtmp(rbdind(i)) = l(rbdind(i))
              else if ( rbdtype(i) .eq. 'U' ) then
                  xtmp(rbdind(i)) = u(rbdind(i))
              end if
          end do
      end if

C     Project

      if ( atmp .gt. amax ) then

          projected = .false.
          do i = 1,nind
              if ( xtmp(i) .lt. l(i) .or. xtmp(i) .gt. u(i) ) then
                  projected = .true.
                  xtmp(i) = max( l(i), min( xtmp(i), u(i) ) )
              end if
          end do

      end if

C     Test if this is not the same point as the previous one. This test 
C     is performed only when xtmp is in fact a projected point.

      if ( projected ) then

          samep = .true.
          do i = 1,nind
              if ( abs( xtmp(i) - xp(i) ) .gt. 
     +             macheps * max( 1.0d0, abs( xp(i) ) ) ) then
                  samep = .false.
              end if
          end do

          if ( samep ) then

C             Finish the extrapolation with the current point

              call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
              if ( inform .lt. 0 ) return

              extinfo = 7

              if ( iprintinn .ge. 6 ) then
                  write(*, 940)
                  write(file10_unit,940)
              end if

              return

          end if

      end if

C     Evaluate function

      call calcal(nind,xtmp,m,lambda,rho,equatn,linear,ftmp,inform)

C     If the objective function is not well defined in an extrapolated 
C     point, we discard all the extrapolated points and return to a 
C     safe region (to the last point before the extrapolation)

      if ( inform .lt. 0 ) then

          fp = fbext

          do i = 1,nind
              xp(i) = xbext(i)
          end do

          call csetp(nind,xp)

          call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
          if ( inform .lt. 0 ) return

          extinfo = 8

          if ( iprintinn .ge. 6 ) then
              write(*, 950)
              write(file10_unit,950)
          end if

          return

      end if

C     Print information of this iteration

      if ( iprintinn .ge. 6 ) then
          write(*, 100) atmp,ftmp,fcnt
          write(file10_unit,100) atmp,ftmp,fcnt
      end if

C     If the functional value decreases then set the current point and 
C     continue the extrapolation

      if ( ftmp .lt. fp ) then

          alpha = atmp

          fp = ftmp

          do i = 1,nind
              xp(i) = xtmp(i)
          end do

          go to 1010

      end if

C     If the functional value does not decrease then discard the last 
C     trial and finish the extrapolation with the previous point

      call csetp(nind,xp)

      call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
      if ( inform .lt. 0 ) return

      extinfo = 9

      if ( iprintinn .ge. 6 ) then
          write(*, 960)
          write(file10_unit,960)
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 110  format(  5X,'Beta condition also holds. ',
     +            'No extrapolation is done.')
 120  format(  5X,'Beta condition does not hold. We will extrapolate.')

 910  format(  5X,'Flag of Extrapolation: Unbounded objective ',
     +            'function?')
 930  format(  5X,'Flag of Extrapolation: Maximum of consecutive ',
     +            'extrapolations reached.')
 940  format(  5X,'Flag of Extrapolation: Very similar consecutive ',
     +            'projected points.')
 950  format(  5X,'Flag of Extrapolation: Not-well-defined objective ',
     +            'function in an extrapolated point.')
 960  format(  5X,'Flag of Extrapolation: Functional value increased ',
     +            'when extrapolating.')

      end
