C     *****************************************************************
C     *****************************************************************

      subroutine checkd(n,l,u,m,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision l(n),u(n)

C     This subrotutine checks the user supplied first and second
C     derivatives subroutines (evalg, evalh, evaljac and evalhc) for
C     computing the objective function gradient and Hessian and the
C     constraints gradients and Hessians, respectively.

#include "dim.par"
#include "algparam.com"
#include "outtyp.com"

C     LOCAL SCALARS
      character answer
      integer i,j
      double precision drand,seed,smalll,smallu

C     LOCAL ARRAYS
      double precision x(nmax)

C     EXTERNAL FUNCTIONS
      external drand

      call vunsetp()

C     SET A RANDOM POINT 

      seed = 123456.0d0
      do i = 1,n
          smalll = max( l(i), - 10.0d0 )
          smallu = min( u(i),   10.0d0 )
          if ( .not. smalll .lt. smallu ) then
              smalll = l(i)
              smallu = u(i)
          end if
          x(i) = smalll + ( smallu - smalll ) * drand(seed)
      end do
 
      write(* ,100)
      write(file10_unit,100)

      do i = 1,n
          write(* ,110) i,x(i)
          write(file10_unit,110) i,x(i)
      end do

C     CHECK OBJECTIVE FUNCTION GRADIENT

      if ( .not. gcoded ) then
          write(* ,160) 'evalg'
          write(file10_unit,160) 'evalg'

          go to 1000
      end if

      write(* ,120)
      write(file10_unit,120)

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkg(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK CONSTRAINTS GRADIENTS

 1000 continue

      if ( .not. jaccoded ) then
          write(* ,160) 'evaljac'
          write(file10_unit,160) 'evaljac'

          go to 1020
      end if

      j = 1

 1010 if ( j .le. m ) then

          write(* ,130) j
          write(file10_unit,130) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              go to 1020

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkjac(n,x,j,inform)
              if ( inform .lt. 0 ) return
          end if

          j = j + 1

          go to 1010

      end if

C     CHECK HESSIAN OF THE OBJECTIVE FUNCTION

 1020 continue

      if ( .not. hcoded ) then
          write(* ,160) 'evalh'
          write(file10_unit,160) 'evalh'

          go to 1030
      end if

      write(* ,140)
      write(file10_unit,140)

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkh(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK HESSIANS OF THE CONSTRAINTS

 1030  continue

      if ( .not. hccoded ) then
          write(* ,160) 'evalhc'
          write(file10_unit,160) 'evalhc'

          go to 1050
      end if

      j = 1

 1040 if ( j .le. m ) then

          write(* ,150) j
          write(file10_unit,150) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              return

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkhc(n,x,j,inform)
              if ( inform .lt. 0 ) return
          end if

          j = j + 1

          go to 1040

      end if

C     CHECK HESSIAN OF THE LAGRANGIAN

 1050 continue

      if ( .not. hlcoded ) then
          write(* ,160) 'evalhl'
          write(file10_unit,160) 'evalhl'

          go to 1060
      end if

      write(* ,*) 'Test of evalhl not implemented yet!'
      write(file10_unit,*) 'Test of evalhl not implemented yet!'

C     CHECK HESSIAN OF THE LAGRANGIAN TIMES A VECTOR

 1060 continue

      if ( .not. hlpcoded ) then
          write(* ,160) 'evalhlp'
          write(file10_unit,160) 'evalhlp'

          return
      end if

      write(* ,*) 'Test of evalhlp not implemented yet!'
      write(file10_unit,*) 'Test of evalhlp not implemented yet!'

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Derivatives will be tested at the random point: ')
 110  format(  1X,'x(',I6,') = ',1P,D15.8)
 120  format(/,1X,'Check the gradient of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 130  format(/,1X,'Check the gradient of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')
 140  format(/,1X,'Check the Hessian matrix of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 150  format(/,1X,'Check the Hessian of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')
 160  format(/,1X,'Skipping test of uncoded ',A7,' subroutine.')

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkg(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalg for 
C     computing the gradient of the objective function using central
C     finite differences with two different discretization steps.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer i
      double precision fminus,fplus,gdiff1,gdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      double precision g(nmax)

      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

      write(* ,100)
      write(file10_unit,100)

      maxerr = 0.0d0

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

          write(* ,110) i,g(i),gdiff1,gdiff2,tmp
          write(file10_unit,110) i,g(i),gdiff1,gdiff2,tmp

          maxerr = max( maxerr, tmp )

      end do

      write(* ,120) maxerr
      write(file10_unit,120) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of the objective function.',
     +       /,1X,'Index',13X,'evalg',2X,'Central diff (two different ',
     +            'steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkh(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalh for 
C     computing the Hessian of the objective function using central 
C     finite differences with two different discretization steps.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
          write(*, 100) nsmax,nsmax
          write(file10_unit,100) nsmax,nsmax

          return
      end if

C     Compute the gradient of the objective function at x

      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

C     Compute the Hessian of the objective function at x and save in a
C     dense matrix

      call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

C     Test column by column

      write(* ,200)
      write(file10_unit,200)

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vevalg(n,x,gplus1,inform)
          if ( inform .lt. 0 ) return

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vevalg(n,x,gplus2,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp

          write(* ,210) j
          write(file10_unit,210) j

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
              if ( elem   .ne. 0.0d0 .or. 
     +             hdiff1 .ne. 0.0d0 .or. 
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      write(* ,220)
                      write(file10_unit,220)
                  end if
                  write(* ,230) i,elem,hdiff1,hdiff2,tmp
                  write(file10_unit,230) i,elem,hdiff1,hdiff2,tmp
              end if
              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( nullcol ) then
              write(* ,240)
              write(file10_unit,240)
          else
              write(* ,250) maxcoe(j)
              write(file10_unit,250) maxcoe(j)
          end if

      end do

      write(* ,*)
      write(file10_unit,*)

      do j = 1,n
          write(* ,260) j,maxcoe(j)
          write(file10_unit,260) j,maxcoe(j)
      end do

      write(* ,270) maxerr
      write(file10_unit,270) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKH uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of the objective function column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',13X,'evalh',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkjac(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evaljac for 
C     computing the gradients of the constraints using central finite 
C     differences with two different discretization steps.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical nullcol
      integer i,jcnnz
      double precision cminus,cplus,jacdiff1,jacdiff2,maxerr,step1,
     +        step2,tmp

C     LOCAL ARRAYS
      integer jcvar(nmax)
      double precision g(nmax),jcval(nmax)

C     COMPUTE THE GRADIENT OF THE CONSTRAINT AND SAVE IT INTO A DENSE 
C     VECTOR

      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     COMPARE WITH CENTRAL FINITE DIFFERENCES

      write(* ,100) ind
      write(file10_unit,100) ind

      maxerr = 0.0d0

      nullcol = .true.

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          jacdiff1 = ( cplus - cminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          jacdiff2 = ( cplus - cminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - jacdiff1 ), abs( g(i) - jacdiff2 ) )

          if ( g(i)     .ne. 0.0d0 .or. 
     +         jacdiff1 .ne. 0.0d0 .or. 
     +         jacdiff2 .ne. 0.0d0 ) then
              if ( nullcol ) then
                  nullcol = .false.
                  write(* ,110)
                  write(file10_unit,110)
              end if
              write(* ,120) i,g(i),jacdiff1,jacdiff2,tmp
              write(file10_unit,120) i,g(i),jacdiff1,jacdiff2,tmp
          end if

          maxerr = max( maxerr, tmp )
      end do

      if ( nullcol ) then
          write(* ,130)
          write(file10_unit,130)
      else
          write(* ,140) maxerr
          write(file10_unit,140) maxerr
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of constraints ',I5,'.')
 110  format(/,1X,'Index',11X,'evaljac',2X,'Central diff (two ',
     +            'different steps)',4X,'Absolute error')
 120  format(  1X,I5,4(3X,1P,D15.8))
 130  format(  1X,'All the elements of this gradient are null.')
 140  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************
      subroutine checkhc(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalhc for 
C     computing the Hessians of the constraints using finite
C     differences.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz,jcnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2),jcvar(nsmax)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),jcval(nsmax),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
          write(*, 100) nsmax,nsmax
          write(file10_unit,100) nsmax,nsmax

          return
      end if

C     Compute the gradient of constraint ind at x and save it in a 
C     dense vector

      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     Compute the Hessian of constraint ind at x and save it in a 
C     dense matrix

      call vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

      write(* ,200) ind
      write(file10_unit,200) ind

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

C         Compute the gradient of constraint ind at xplus1 and 
C         save in a dense vector

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus1(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus1(jcvar(i)) = jcval(i)
          end do

C         Compute the gradient of constraint ind at xplus2 and 
C         save in a dense vector

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus2(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus2(jcvar(i)) = jcval(i)
          end do

          x(j) = tmp

          write(* ,210) j
          write(file10_unit,210) j

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )

              if ( elem   .ne. 0.0d0 .or.
     +             hdiff1 .ne. 0.0d0 .or. 
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      write(* ,220)
                      write(file10_unit,220)
                  end if
                  write(* ,230) i,elem,hdiff1,hdiff2,tmp
                  write(file10_unit,230) i,elem,hdiff1,hdiff2,tmp
              end if

              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( nullcol ) then
              write(* ,240)
              write(file10_unit,240)
          else
              write(* ,250) maxcoe(j)
              write(file10_unit,250) maxcoe(j)
          end if

      end do

      write(* ,*)
      write(file10_unit,*)

      do j = 1,n
          write(* ,260) j,maxcoe(j)
          write(file10_unit,260) j,maxcoe(j)
      end do

      write(* ,270) maxerr
      write(file10_unit,270) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKHC uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of constraint ',I5,' column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',12X,'evalhc',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end

