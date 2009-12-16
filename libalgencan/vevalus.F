C     SAFE CALL TO THE USER PROVIDED SUBROUTINES

C     ******************************************************************
C     ******************************************************************

      subroutine vinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(6),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "dim.par"

C     LOCAL SCALARS
      integer i

C     INTRINSEC FUNCTIONS
      intrinsic max,min

C     EXTERNAL SUBROUTINES
      external checkd

C     Avoid huge bounds

      do i = 1,n
          l(i) = max( l(i), - 1.0d+20 )
          u(i) = min( u(i),   1.0d+20 )
      end do

C     Project initial guess

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Check derivatives

      if ( checkder ) then
          call checkd(n,l,u,m,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n,file20_unit
      logical :: connected

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

#include "outtyp.com"

C     LOCAL SCALARS
      integer i

      if ( iprintctl(2) ) then

C         Save solution
          file20_unit_loop: do file20_unit=10, 99
    
            inquire(unit=file20_unit, opened=connected)
    
            if (.not.connected) exit file20_unit_loop
    
          end do file20_unit_loop

          open(file20_unit,file='solution.txt')

C         Point

          write(file20_unit,100)
          do i = 1,n
              write(file20_unit,200) i,x(i)
          end do

C         Lagrange multipliers

          if ( m .gt. 0 ) then
              write(file20_unit,300)
              do i = 1,m
                  write(file20_unit,200) i,lambda(i)
              end do
          end if

          close(file20_unit)

      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,'FINAL POINT:',//,2X,'INDEX',16X,'X(INDEX)')
 200  format(I7,1P,D24.16)
 300  format(/,'FINAL ESTIMATION OF THE LAGRANGE MULTIPLIERS: ',
     +       //,2X,'INDEX',11X,'LAMBDA(INDEX)')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalf,reperr

      call evalf(n,x,f,flag)

      efcnt = efcnt + 1
      fcnt  = fcnt  + 1

      if ( flag .ne. 0 ) then
          inform = - 90
          call reperr(inform)
          return
      end if

      if ( .not. IsANumber(f) ) then
          if ( iprintctl(3) ) then
              write(* ,100)
              write(* ,200) f
              write(file10_unit,100)
              write(file10_unit,200) f
          end if

          inform = - 90
          call reperr(inform)
          return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALF WARNING: The objective function value ',
     +            'computed by the user-supplied subroutine EVALF is ',
     +            '+Inf, -Inf or NaN.')

 200  format(/,1X,'Value: ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

#include "dim.par"
#include "algparam.com"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalg,ievalg,reperr

      if ( gcoded ) then
          call evalg(n,x,g,flag)

          egcnt = egcnt + 1

          if ( flag .ne. 0 ) then
              inform = - 92
              call reperr(inform)
              return
          end if

          do i = 1,n
              if ( .not. IsANumber(g(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,200) n,i,g(i)
                      write(file10_unit,100)
                      write(file10_unit,200) n,i,g(i)

                      inform = - 92
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ievalg(n,x,g,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALG WARNING: There is an element whose value is ',
     +            '+Inf, -Inf or NaN in the gradient of the objective ',
     +            'function computed by the user-supplied subroutine ',
     +            'EVALG.')

 200  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

#include "machconst.com"

C     LOCAL SCALARS
      integer j
      double precision fminus,fplus,step,tmp

C     INTRINSIC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external setp,vevalf

      do j = 1,n
          tmp  = x(j)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(j) = tmp + step
          call setp(n,x)
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp - step
          call setp(n,x)
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          g(j) = ( fplus - fminus ) / ( 2.0d0 * step )

          x(j) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalh,reperr

      call evalh(n,x,hlin,hcol,hval,hnnz,flag)

      ehcnt = ehcnt + 1

      if ( flag .ne. 0 ) then
          inform = - 94
          call reperr(inform)
          return
      end if

      do i = 1,hnnz
          if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +         hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. 
     +         hcol(i) .gt. hlin(i) ) then

              if ( iprintctl(3) ) then
                  write(* ,100)
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(file10_unit,100)
                  write(file10_unit,300) n,i,hlin(i),hcol(i),hval(i)
              end if

              hlin(i) = 1
              hcol(i) = 1
              hval(i) = 0.0d0
          end if

          if ( .not. IsANumber(hval(i)) ) then
              if ( iprintctl(3) ) then
                  write(* ,200) 
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(file10_unit,200) 
                  write(file10_unit,300) n,i,hlin(i),hcol(i),hval(i)

                  inform = - 94
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALH WARNING: There is an element out of range, ',
     +            'or in the upper triangle, of the Hessian of the ',
     +            'objetive function computed by the user-supplied ',
     +            'subroutine EVALH. It will be ignored.')

 200  format(/,1X,'VEVALH WARNING: There is an element whose value is ',
     +            '+Inf, -Inf or NaN in the Hessian of the objetive ',
     +            'function computed by the user-supplied subroutine ',
     +            'EVALH.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalc,reperr

      call evalc(n,x,ind,c,flag)

      eccnt(ind) = eccnt(ind) + 1

      if ( flag .ne. 0 ) then
          inform = - 91
          call reperr(inform)
          return
      end if

      if ( .not. IsANumber(c) ) then
          if ( iprintctl(3) ) then
              write(* ,100) ind
              write(* ,200) c
              write(file10_unit,100) ind
              write(file10_unit,200) c

              inform = - 91
              call reperr(inform)
              return
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALC WARNING: The value of constraint ',I16,' ',
     +            'computed by the user-supplied subroutine EVALC is ',
     +            '+Inf, -Inf or NaN.')

 200  format(/,1X,'Value: ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

#include "dim.par"
#include "algparam.com"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evaljac,ievaljac,reperr

      if ( jaccoded ) then
          call evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

          ejccnt(ind) = ejccnt(ind) + 1

          if ( flag .ne. 0 ) then
              inform = - 93
              call reperr(inform)
              return
          end if

          do i = 1,jcnnz
              if ( jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then

                  if ( iprintctl(3) ) then
                      write(* ,100) ind
                      write(* ,300) n,i,jcvar(i),jcval(i)
                      write(file10_unit,100) ind
                      write(file10_unit,300) n,i,jcvar(i),jcval(i)
                  end if

                  jcvar(i) = 1
                  jcval(i) = 0.0d0
              end if

              if ( .not. IsANumber(jcval(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,200) ind
                      write(* ,300) n,i,jcvar(i),jcval(i)
                      write(file10_unit,200) ind
                      write(file10_unit,300) n,i,jcvar(i),jcval(i)

                      inform = - 93
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ievaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALJAC WARNING: There is an element out of ',
     +            'range in the gradient of constraint ',I16,' ',
     +            'computed by the user-supplied subroutine EVALJAC. ',
     +            'It will be ignored.')

 200  format(/,1X,'VEVALJAC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the gradient of constraint ',
     +            I16,'computed by the user-supplied subroutine ',
     +            'EVALJAC.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Variable : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision jcval(n),x(n)

#include "dim.par"
#include "machconst.com"

C     LOCAL SCALARS
      integer j
      double precision cminus,cplus,step,tmp

C     INTRINSEC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external setp,vevalc

      jcnnz = 0

      do j = 1,n
          tmp  = x(j)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(j) = tmp + step
          call setp(n,x)
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp - step
          call setp(n,x)
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          jcvar(jcnnz + 1) = j
          jcval(jcnnz + 1) = ( cplus - cminus ) / ( 2.0d0 * step )

          if ( abs( jcval(jcnnz + 1) ) .gt. 0.0d0 ) then
              jcnnz = jcnnz + 1
          end if

          x(j) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

#include "dim.par"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhc,reperr

      call evalhc(n,x,ind,hlin,hcol,hval,hnnz,flag)

      ehccnt(ind) = ehccnt(ind) + 1

      if ( flag .ne. 0 ) then
          inform = - 95
          call reperr(inform)
          return
      end if

      do i = 1,hnnz
          if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +         hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. 
     +         hcol(i) .gt. hlin(i) ) then

              if ( iprintctl(3) ) then
                  write(* ,100) ind
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(file10_unit,100) ind
                  write(file10_unit,300) n,i,hlin(i),hcol(i),hval(i)
              end if

              hlin(i) = 1
              hcol(i) = 1
              hval(i) = 0.0d0
          end if

          if ( .not. IsANumber(hval(i)) ) then
              if ( iprintctl(3) ) then
                  write(* ,200) ind
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(file10_unit,200) ind
                  write(file10_unit,300) n,i,hlin(i),hcol(i),hval(i)

                  inform = - 95
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHC WARNING: There is an element out of range ',
     +            'or in the upper triangle of the Hessian of ',
     +            'constraint ',I16,' computed by the user-supplied ',
     +            'subroutine EVALHC. It will be ignored.')

 200  format(/,1X,'VEVALHC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the Hessian of constraint ',
     +            I16,' computed by the user-supplied subroutine ',
     +            'EVALHC.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

#include "dim.par"
#include "algparam.com"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhl,ievalhl,reperr

      if ( hlcoded ) then

          call evalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,flag)

          ehlcnt = ehlcnt + 1

          if ( flag .ne. 0 ) then
              inform = - 96
              call reperr(inform)
              return
          end if

          do i = 1,hnnz
              if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +             hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. 
     +             hcol(i) .gt. hlin(i) ) then

                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                      write(file10_unit,100)
                      write(file10_unit,300) n,i,hlin(i),hcol(i),hval(i)
                  end if

                  hlin(i) = 1
                  hcol(i) = 1
                  hval(i) = 0.0d0
              end if

              if ( .not. IsANumber(hval(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,200)
                      write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                      write(file10_unit,200)
                      write(file10_unit,300) n,i,hlin(i),hcol(i),hval(i)

                      inform = - 96
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else if ( hcoded .and. ( hccoded .or. m .eq. 0 ) ) then
         call ievalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
         if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHL WARNING: There is an element out of range, ',
     +            'or in the upper triangle, of the',/,1X,'Hessian of ',
     +            'the Lagrangian computed by the user-supplied ',
     +            'subroutine EVALHL. It',/,1X,'will be ignored.')

 200  format(/,1X,'VEVALHL WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the',/,1X,'Hessian of the ',
     +            'Lagrangian computed by the user-supplied ',
     +            'subroutine EVALHL.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

#include "dim.par"

C     LOCAL SCALARS
      integer col,con,hnnztmp,i,ind,itmp,j,lin,nextj,rnnz
      double precision val

C     LOCAL ARRAYS
      integer hcon(hnnzmax),pos(nmax),rind(nmax),stlin(nmax)
      double precision rval(nmax)

C     EXTERNAL SUBROUTINES
      external vevalh,vevalhc

C     ==================================================================
C     COMPUTE HESSIANS
C     ==================================================================

C     COMPUTE HESSIAN OF THE OBJECTIVE FUNCTION

      call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

C     For each element of the Hessian of the objective function,
C     set constraint index as zero
      do i = 1,hnnz
          hval(i) = hval(i) / sf
          hcon(i) = 0
      end do

      if ( m .eq. 0 ) return

C     COMPUTE HESSIANS OF THE CONSTRAINTS

      ind = 0

      do j = 1,m
C         Compute Hessian of constraint j
          call vevalhc(n,x,j,hlin(hnnz+ind+1),hcol(hnnz+ind+1),
     +    hval(hnnz+ind+1),hnnztmp,inform)
          if ( inform .lt. 0 ) return

C         For each element of the Hessian, set constraint as j
          do i = hnnz + ind + 1,hnnz + ind + hnnztmp
              hval(i) = hval(i) / sc(j)
              hcon(i) = j
          end do
         
          ind = ind + hnnztmp
      end do

      if ( ind .eq. 0 ) return

      hnnz = hnnz + ind

C     ==================================================================
C     SET ROW LINKED LISTS
C     ==================================================================

C     Initialize pointers to the first element of each row
      do i = 1,n
         stlin(i) = 0
      end do

C     Set row linked lists
      do i = 1,hnnz
         lin        = hlin(i)
         itmp       = stlin(lin)
         stlin(lin) = i
         hlin(i)    = itmp
      end do

C     ==================================================================
C     BUILD HESSIAN OF THE LAGRANGIAN ROW BY ROW
C     ==================================================================

C     Initialize array pos
      do i = 1,n
          pos(i) = 0
      end do

      do i = 1,n
C         Initialize the i-th row of the Hessian of the Lagrangian
          rnnz = 0

C         Process the i-th row of all the Hessians
          j = stlin(i)

 10       if ( j .ne. 0 ) then

C             Process element (i,hcol(j)) of the Hessian of constraint 
C             hcon(j) (Hessian of the objective function if hcon(j)=0)

              col = hcol(j)
              con = hcon(j)
              if ( con .eq. 0 ) then
                  val = hval(j)
              else
                  val = hval(j) * lambda(con)
              end if

              if ( pos(col) .ne. 0 ) then
                  rval(pos(col)) = rval(pos(col)) + val

              else
                  rnnz           = rnnz + 1
                  pos(col)       = rnnz
                  rind(pos(col)) = col
                  rval(pos(col)) = val
              end if

C             Get next element in the i-th row linked list
              j = hlin(j)
              go to 10
          end if

C         Clean array pos
          do j = 1,rnnz
              pos(rind(j)) = 0
          end do

C         Set i-th row of hl (over the i-th rows of the Hessians) 
C         and mark remaining elements to be deleted
          j = stlin(i)

 20       if ( j .ne. 0 ) then
              nextj = hlin(j)

              if ( rnnz .ne. 0 ) then
                  hlin(j) = i
                  hcol(j) = rind(rnnz)
                  hval(j) = rval(rnnz)
                  rnnz    = rnnz - 1
              else
                  hlin(j) = 0
              end if

              j = nextj
              go to 20
          end if

      end do

C     Eliminate remaining elements (marked with hlin(j)=0)
      j = 1

 30   if ( j .le. hnnz ) then
          if ( hlin(j) .eq. 0 ) then
              if ( j .ne. hnnz ) then
                  hlin(j) = hlin(hnnz)
                  hcol(j) = hcol(hnnz)
                  hval(j) = hval(hnnz)
              end if
              hnnz = hnnz - 1
          else
              j = j + 1
          end if

          go to 30
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

#include "dim.par"
#include "algparam.com"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhlp,ievalhlp,reperr

      if ( hlpcoded ) then
          call evalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,flag)

          ehlpcnt = ehlpcnt + 1

          if ( flag .ne. 0 ) then
              inform = - 97
              call reperr(inform)
              return
          end if

          do i = 1,n
              if ( .not. IsANumber(hp(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,200) n,i,hp(i)
                      write(file10_unit,100)
                      write(file10_unit,200) n,i,hp(i)

                      inform = - 97
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else if ( truehl ) then
          call ievalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHLP WARNING: There is an element in the ',
     +            'product of the Hessian of the Lagrangian by a ',
     +            'computed by the user-supplied subroutine EVALHLP ',
     +            'whose value is +Inf, -Inf or NaN.')

 200  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

#include "dim.par"
#include "hessdat.com"

C     LOCAL SCALARS
      integer col,i,lin
      double precision val

C     EXTERNAL SUBROUTINES
      external vevalhl

      if ( .not. gothl ) then
          gothl = .true.

          call vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      do i = 1,n
          hp(i) = 0.0d0
      end do

      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          val = hval(i)

          hp(lin) = hp(lin) + p(col) * val 

          if ( lin .ne. col ) then
              hp(col) = hp(col) + p(lin) * val 
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

#include "dim.par"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUCNTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalfc,reperr

      call evalfc(n,x,f,m,c,flag)

      efccnt = efccnt + 1
      fcnt   = fcnt   + 1

      if ( flag .ne. 0 ) then
          inform = - 98
          call reperr(inform)
          return
      end if

      if ( .not. IsANumber(f) ) then
          if ( iprintctl(3) ) then
              write(* ,100)
              write(* ,300) f
              write(file10_unit,100)
              write(file10_unit,300) f

              inform = - 98
              call reperr(inform)
              return
          end if
      end if

      do i = 1,m
          if ( .not. IsANumber(c(i)) ) then
              if ( iprintctl(3) ) then
                  write(* ,200)
                  write(* ,400) n,m,i,c(i)
                  write(file10_unit,200)
                  write(file10_unit,400) n,m,i,c(i)

                  inform = - 98
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALFC WARNING: The objective function value ',
     +            'computed by the user-supplied',/,1X,'subroutine ',
     +            'EVALFC is +Inf, -Inf or NaN.')

 200  format(/,1X,'VEVALFC WARNING: The value of a constraint ',
     +            'computed by the user-supplied',/,1X,'subroutine ',
     +            'EVALFC is +Inf, -Inf or NaN.')

 300  format(/,1X,'Value: ',1P,D24.16)

 400  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Number of constraints : ',I16,
     +       /,1X,'Constraint            : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

#include "dim.par"
#include "algparam.com"
#include "counters.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUCNTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalgjac,ievalgjac,reperr

      if ( gjaccoded ) then
          call evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

          egjccnt = egjccnt + 1

          if ( flag .ne. 0 ) then
              inform = - 99
              call reperr(inform)
              return
          end if

          do i = 1,jcnnz
              if ( jcfun(i) .lt. 1 .or. jcfun(i) .gt. m .or.
     +             jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then

                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,300) n,m,i,jcfun(i),jcvar(i),jcval(i)
                      write(file10_unit,100)
                      write(file10_unit,300) n,m,i,jcfun(i),jcvar(i),jcval(i)
                  end if

                  jcfun(i) = 1
                  jcvar(i) = 1
                  jcval(i) = 0.0d0
              end if

              if ( .not. IsANumber(jcval(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,200)
                      write(* ,300) n,m,i,jcfun(i),jcvar(i),jcval(i)
                      write(file10_unit,200)
                      write(file10_unit,300) n,m,i,jcfun(i),jcvar(i),jcval(i)

                      inform = - 99
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ievalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALGJAC WARNING: There is an element out of ',
     +            'range in the gradient of the objective function or ',
     +            'in the Jacobian of the constraints computed by the ',
     +            'user-supplied subroutine EVALGJAC. It will be ',
     +            'ignored.')

 200  format(/,1X,'VEVALGJAC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the',/,1X,'gradient of the ',
     +            'objective function or in the Jacobian of the ',
     +            'constraints',/,1X,'computed by the user-supplied ',
     +            'subroutine EVALGJAC.')

 300  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Number of constraints : ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Constraint            : ',I16,
     +       /,1X,'Variable              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

#include "dim.par"
#include "machconst.com"

C     LOCAL SCALARS
      integer i,j
      double precision fminus,fplus,step,tmp

C     LOCAL ARRAYS
      double precision cminus(mmax),cplus(mmax)

C     INTRINSIC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external setp,vevalfc

      jcnnz = 0

      do i = 1,n
          tmp  = x(i)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(i) = tmp + step
          call setp(n,x)
          call vevalfc(n,x,fplus,m,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step
          call setp(n,x)
          call vevalfc(n,x,fminus,m,cminus,inform)
          if ( inform .lt. 0 ) return

          do j = 1,m
              jcfun(jcnnz + 1) = j
              jcvar(jcnnz + 1) = i
              jcval(jcnnz + 1) = ( cplus(j) - cminus(j) ) /
     +                           ( 2.0d0 * step )

              if ( abs( jcval(jcnnz + 1) ) .gt. 0.0d0 ) then
                  jcnnz = jcnnz + 1
              end if
          end do

          g(i) = ( fplus - fminus ) / ( 2.0d0 * step )

          x(i) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine reperr(inform)

C     SCALAR ARGUMENTS
      integer inform

#include "outtyp.com"

      if ( iprintctl(3) ) then

          if ( inform .eq. -90 ) then
              write(* ,100) 'EVALF' 
              write(file10_unit,100) 'EVALF' 

          else if ( inform .eq. -91 ) then
              write(* ,100) 'EVALC' 
              write(file10_unit,100) 'EVALC' 

          else if ( inform .eq. -92 ) then
              write(* ,100) 'EVALG' 
              write(file10_unit,100) 'EVALG' 

          else if ( inform .eq. -93 ) then
              write(* ,100) 'EVALJAC' 
              write(file10_unit,100) 'EVALJAC' 

          else if ( inform .eq. -94 ) then
              write(* ,100) 'EVALH' 
              write(file10_unit,100) 'EVALH' 

          else if ( inform .eq. -95 ) then
              write(* ,100) 'EVALHC' 
              write(file10_unit,100) 'EVALHC' 

          else if ( inform .eq. -96 ) then
              write(* ,100) 'EVALHL' 
              write(file10_unit,100) 'EVALHL' 

          else if ( inform .eq. -97 ) then
              write(* ,100) 'EVALHLP' 
              write(file10_unit,100) 'EVALHLP' 

          else if ( inform .eq. -98 ) then
              write(* ,100) 'EVALFC' 
              write(file10_unit,100) 'EVALFC' 

          else if ( inform .eq. -99 ) then
              write(* ,100) 'EVALGJAC' 
              write(file10_unit,100) 'EVALGJAC' 
          end if

      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'*** There was an error in the user supplied ',
     +            'subroutine ',A10,' ***',/)

      end

C     ******************************************************************
C     ******************************************************************

      logical function IsANumber(x)

      implicit none

C     SCALAR ARGUMENTS
      double precision x

#include "machconst.com"

      IsANumber = .true.

      if ( .not. abs( x ) .le. bignum ) then
          IsANumber = .false.
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vsetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

#include "dim.par"
#include "graddat.com"

C     EXTERNAL SUBROUTINES
      external setp

      gotc = .false.

      call setp(n,x)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vunsetp()

      implicit none

#include "dim.par"
#include "graddat.com"

C     EXTERNAL SUBROUTINES
      external unsetp

      gotc = .false.

      call unsetp()

      end
