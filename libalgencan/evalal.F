C     *****************************************************************
C     *****************************************************************

      subroutine sevalobjc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision f
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

#include "algparam.com"

      call ssetp(n,x)

      if ( fccoded ) then
          call sevalobjcb(n,x,f,m,c,inform)

      else ! if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then
          call sevalobjca(n,x,f,m,c,inform)
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalobjcb(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision f
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     COMPUTE OBJECTIVE FUNTION AND CONSTRAINTS

      call sevalfc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalobjca(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision f
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     LOCAL SCALARS
      integer j

C     COMPUTE OBJECTIVE FUNTION

      call sevalf(n,x,f,inform)
      if ( inform .lt. 0 ) return

C     COMPUTE CONSTRAINTS

      do j = 1,m

C         COMPUTE THE j-TH CONSTRAINT
          call sevalc(n,x,j,c(j),inform)
          if ( inform .lt. 0 ) return

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalal(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

#include "dim.par"
#include "graddat.com"
#include "algparam.com"

      if ( innercall ) then
          call minsqf(n,x,al,inform)
          return
      end if

      call ssetp(n,x)

      if ( fccoded ) then
          call sevalalb(n,x,m,lambda,rho,equatn,linear,al,inform)
      else
          call sevalala(n,x,m,lambda,rho,equatn,linear,al,inform)
      end if

      gotc = .true.

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalalb(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer j
      double precision f,p

C     COMPUTE OBJECTIVE FUNTION AND CONSTRAINTS

      call sevalfc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

C     COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)

      al = f

      do j = 1,m

C         ADD P(c_j, rho_j, lambda_j)
          call evalp(c(j),rho(j),lambda(j),equatn(j),p)
          al = al + p

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalala(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer j
      double precision f,p

C     COMPUTE OBJECTIVE FUNTION

      call sevalf(n,x,f,inform)
      if ( inform .lt. 0 ) return

C     COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)

      al = f

      do j = 1,m

C         COMPUTE j-TH CONSTRAINT
          call sevalc(n,x,j,c(j),inform)
          if ( inform .lt. 0 ) return

C         ADD P(c_j, rho_j, lambda_j)
          call evalp(c(j),rho(j),lambda(j),equatn(j),p)
          al = al + p

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnl(n,x,m,lambda,equatn,linear,nl,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

#include "algparam.com"

      if ( gjaccoded ) then
          call sevalnlb(n,x,m,lambda,equatn,linear,nl,inform)
      else
          call sevalnla(n,x,m,lambda,equatn,linear,nl,inform)
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnlb(n,x,m,lambda,equatn,linear,nl,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer i,j,jcnnz

C     LOCAL ARRAYS
      integer jcfun(jcnnzmax)

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND JACOBIAN
C     OF CONSTRAINTS

      call sevalgjac(n,x,nl,m,jcfun,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = nl(i)
      end do

C     CONVERT JACOBIAN OF CONSTRAINTS FROM COORDINATE FORMAT TO
C     COMPRESSED SPARSE ROW FORMAT

      call coo2csr(m,jcnnz,jcfun,jcvar,jcval,jclen,jcsta)

C     COMPUTE \nabla L = \nabla f + \sum_j lambda_j * \nabla c_j

      constrc = .false.

      do j = 1,m

          if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then

C             ADD lambda_j * \nabla c_j

              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
              end do

              if ( .not. linear(j) ) then
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      g(jcvar(i)) = g(jcvar(i)) + lambda(j) * jcval(i)
                  end do
              end if

              constrc = .true.

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnla(n,x,m,lambda,equatn,linear,nl,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer i,ind,j

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION

      call sevalg(n,x,nl,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = nl(i)
      end do

C     COMPUTE \nabla L = \nabla f + \sum_j lambda_j * \nabla c_j

      ind = 0

      constrc = .false.

      do j = 1,m

          if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then

              jcsta(j) = ind + 1

C             COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
              call sevaljac(n,x,j,jcvar(ind+1),jcval(ind+1),jclen(j),
     +        inform)
              if ( inform .lt. 0 ) return

              ind = ind + jclen(j)

C             ADD lambda_j * \nabla c_j
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
              end do

              if ( .not. linear(j) ) then
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      g(jcvar(i)) = g(jcvar(i)) + lambda(j) * jcval(i)
                  end do
              end if

              constrc = .true.

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnal(n,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

#include "dim.par"
#include "graddat.com"
#include "algparam.com"

      if ( innercall ) then
          call minsqg(n,x,nal,inform)
          return
      end if

      if ( gjaccoded ) then
          call sevalnalb(n,x,m,lambda,rho,equatn,linear,nal,inform)
      else
          call sevalnala(n,x,m,lambda,rho,equatn,linear,nal,inform)
      end if

      gotc = .true.

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnalb(n,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer j
      double precision dum

C     COMPUTE CONSTRAINTS

      if ( .not. gotc ) then
          call sevalfc(n,x,dum,m,c,inform)
          if ( inform .lt. 0 ) return
      end if

      do j = 1,m
C         COMPUTE dP/dc
          call evaldpdy(c(j),rho(j),lambda(j),equatn(j),dpdc(j))
      end do

C     COMPUTE GRADIENT OF THE LAGRANGIAN WITH DPDC INSTEAD OF LAMBDA
      call sevalnlb(n,x,m,dpdc,equatn,linear,nal,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnala(n,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer j

      if ( .not. gotc ) then
          do j = 1,m
C             COMPUTE THE j-TH CONSTRAINT
              call sevalc(n,x,j,c(j),inform)
              if ( inform .lt. 0 ) return
          end do
      end if

      do j = 1,m
C         COMPUTE dP/dc
          call evaldpdy(c(j),rho(j),lambda(j),equatn(j),dpdc(j))
      end do

C     COMPUTE GRADIENT OF THE LAGRANGIAN WITH DPDC INSTEAD OF LAMBDA
      call sevalnla(n,x,m,dpdc,equatn,linear,nal,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalp(y,rho,lambda,equatn,p)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision lambda,p,rho,y

      if ( equatn ) then
          p  = y * ( lambda + 0.5d0 * rho * y )
      else
          if ( lambda + rho * y .ge. 0.0d0 ) then
              p = y * ( lambda + 0.5d0 * rho * y )
          else
              p = - 0.5d0 * lambda ** 2 / rho
          end if
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evaldpdy(y,rho,lambda,equatn,dpdy)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision y,rho,lambda,dpdy

      if ( equatn ) then
          dpdy = lambda + rho * y
      else
          dpdy = max( 0.0d0, lambda + rho * y )
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalnalu(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical ignlin
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

C     This subroutine computes the gradient of the Augmented Lagrangian 
C     function at a point xp, which is near to x, taking care of the
C     non-differentiability. The Augmented Lagrangian gradient must be
C     previously computed at x.

#include "algparam.com"

      if ( innercall ) then
          call minsqg(n,xp,nalp,inform)
          return
      end if

      call ssetp(n,xp)

      if ( gjaccoded ) then
          call ievalnalub(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +    inform)
      else
          call ievalnalua(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +    inform)
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalnalub(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical ignlin
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer i,j,jcpnnz
      double precision dpdcp,dum

C     LOCAL ARRAYS
      integer jcpfun(jcnnzmax),jcplen(mmax),jcpsta(mmax),
     +        jcpvar(jcnnzmax)
      double precision cp(mmax),jcpval(jcnnzmax)

C     COMPUTE CONSTRAINTS AT xp

      call sevalfc(n,xp,dum,m,cp,inform)
      if ( inform .lt. 0 ) return

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND JACOBIAN
C     OF CONSTRAINTS AT xp

      call sevalgjac(n,xp,nalp,m,jcpfun,jcpvar,jcpval,jcpnnz,inform)
      if ( inform .lt. 0 ) return

C     CONVERT JACOBIAN OF CONSTRAINTS FROM COODINATE FORMAT TO
C     COMPRESSED SPARSE ROW FORMAT

      call coo2csr(m,jcpnnz,jcpfun,jcpvar,jcpval,jcplen,jcpsta)

C     COMPUTE \nabla L = \nabla f + \sum_j dPdc * dcdx

      do j = 1,m

          if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and.
     +         .not. ( ignlin .and. linear(j) ) ) then

C             COMPUTE dP/dc

              dpdcp = lambda(j) + rho(j) * cp(j)

              if ( dpdcp .ne. 0.0d0 ) then

C                 ADD dPdc * dcdx

                  do i = jcpsta(j),jcpsta(j) + jcplen(j) - 1
                      nalp(jcpvar(i)) = 
     +                nalp(jcpvar(i)) + dpdcp * jcpval(i)
                  end do

              end if

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalnalua(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +inform)

      implicit none

C     This subroutine computes the gradient of the Augmented Lagrangian 
C     function at a point xp, which is near to x, taking care of the
C     non-differentiability. The Augmented Lagrangian gradient must be
C     previously computed at x.

C     SCALAR ARGUMENTS
      logical ignlin
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer i,j,jcpnnz
      double precision cp,dpdcp

C     LOCAL ARRAYS
      integer jcpvar(nmax)
      double precision jcpval(nmax)

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AT xp

      call sevalg(n,xp,nalp,inform)
      if ( inform .lt. 0 ) return

C     COMPUTE \nabla L = \nabla f + \sum_j dPdc * dcdx

      do j = 1,m

          if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and.
     +         .not. ( ignlin .and. linear(j) ) ) then

C             COMPUTE THE i-TH CONSTRAINT
              call sevalc(n,xp,j,cp,inform)
              if ( inform .lt. 0 ) return

C             COMPUTE dP/dc
              dpdcp = lambda(j) + rho(j) * cp

              if ( dpdcp .ne. 0.0d0 ) then

C                 COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
                  call sevaljac(n,xp,j,jcpvar,jcpval,jcpnnz,inform)
                  if ( inform .lt. 0 ) return

C                 ADD dPdc * dcdx
                  do i = 1,jcpnnz
                      nalp(jcpvar(i)) = 
     +                nalp(jcpvar(i)) + dpdcp * jcpval(i)
                  end do

              end if

          end if

      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhal(n,x,m,lambda,rho,equatn,linear,hallin,halcol,
     +halval,halnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer halnnz,inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer hallin(*),halcol(*)
      double precision lambda(m),halval(*),rho(m),x(n)

C     This subroutine computes the Hessian of the augmented Lagrangian.

#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer i,lin,j,k,l,var

C     LOCAL ARRAYS
      integer stlin(nmax)
      double precision r(nmax)

      call sevalhl(n,x,m,dpdc,hallin,halcol,halval,halnnz,inform)
      if ( inform .lt. 0 ) return

      if ( m .eq. 0 ) return

C     PUT MATRIX INTO A ROW-LINKED LIST

      do i = 1,n
         r(i) = 0.0d0
      end do

      do i = 1,n
          stlin(i) = 0
      end do
      
      do i = 1,halnnz
          lin = hallin(i)
          k   = stlin(lin)
          stlin(lin) = i
          hallin(i)  = k
      end do

C     ADD \sum_j \rho_j * \nabla c_j \nabla c_j^t

      do j = 1,m
         
          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             ADD \rho_j * \nabla c_j \nabla c_j^t

              do k = jcsta(j),jcsta(j) + jclen(j) - 1
               
                  var = jcvar(k)

C                 PUT ROW jcvar(k) INTO A DENSE VECTOR

                  lin = stlin(var)
 10               if ( lin .ne. 0 ) then
                      r(halcol(lin)) = r(halcol(lin)) + halval(lin)
                      lin = hallin(lin)
                      go to 10
                  end if
               
C                 ADD VALUE

                  do l = jcsta(j),jcsta(j) + jclen(j) - 1
                      if ( jcvar(l) .le. var ) then
                          r(jcvar(l)) = 
     +                    r(jcvar(l)) + rho(j) * jcval(l) * jcval(k)
                      end if
                  end do
               
C                 UPDATE VALUES IN HALVAL
               
                  lin = stlin(var)
 20               if ( lin .ne. 0 ) then
                      halval(lin) = r(halcol(lin))
                      r(halcol(lin)) = 0.0d0
                      lin = hallin(lin)
                      go to 20
                  end if
               
C                 INSERT NEW ELEMENTS IN HAL REPRESENTATION

                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      l = jcvar(i)
                      if ( r(l) .ne. 0.0d0 ) then
                          halnnz = halnnz + 1
                          halval(halnnz) = r(l)
                          halcol(halnnz) = l
                          hallin(halnnz) = stlin(var)
                          stlin(var) = halnnz
                          r(l) = 0.0d0
                      end if
                  end do
 
              end do

          end if

      end do
  
C     PUT MATRIX BACK INTO COORDINATE SQUEME

      do i = 1,n

          lin = stlin(i)
 30       if ( lin .ne. 0 ) then
              k = hallin(lin)
              hallin(lin) = i
              lin = k
              go to 30
          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,gothl,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(n),lambda(m),p(n),rho(m),x(n)

#include "dim.par"
#include "graddat.com"
#include "algparam.com"

C     LOCAL SCALARS
      integer i,j
      double precision atp

      if ( innercall ) then
          call minsqhp(n,x,p,hp,gothl,inform)
          return
      end if

C     --------------------------------------------------------------
C     Hessian approximation
C     --------------------------------------------------------------

      if ( hptype .eq. 'HAPPRO' .and. constrc ) then

          call applyhapp(n,m,rho,equatn,gothl,p,hp)

C     --------------------------------------------------------------
C     Incremental quotients
C     --------------------------------------------------------------

      else if ( hptype .eq. 'INCQUO' .or. hptype .eq. 'HAPPRO' ) then

          call ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)
          if ( inform .lt. 0 ) return

C     --------------------------------------------------------------
C     True Hessian
C     --------------------------------------------------------------

      else if ( hptype .eq. 'TRUEHL' ) then

C         Compute Hessian of Lagrangian times p using dpdc 
C         instead of lambda

          call sevalhlp(n,x,m,dpdc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

C         Add rho A^T A

          do j = 1,m

              if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

                  atp = 0.0d0
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      atp = atp + jcval(i) * p(jcvar(i))
                  end do

                  atp = atp * rho(j)

                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
                  end do

              end if

          end do

      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(n),lambda(m),p(n),rho(m),x(n)

C     Computes an approximation of the product of the Hessian of the 
C     Augmented Lagrangian times a vector using incremental quotients.

#include "machconst.com"
#include "dim.par"
#include "graddat.com"

C     LOCAL SCALARS
      integer i,j
      double precision atp,psupn,step,xsupn

C     LOCAL  ARRAYS
      double precision xp(nmax)

C     --------------------------------------------------------------
C     Set auxiliary point
C     --------------------------------------------------------------

      xsupn = 0.0d0
      psupn = 0.0d0
      do i = 1,n
          xsupn = max( xsupn, abs( x(i) ) )
          psupn = max( psupn, abs( p(i) ) )
      end do

      step = macheps12 * max( xsupn / psupn, 1.0d0 )

      do i = 1,n
          xp(i) = x(i) + step * p(i)
      end do

C     --------------------------------------------------------------
C     Compute gradient of the augmented Lagrangian at xp considering 
C     the same constraints considered at x and ignoring linear 
C     constraints
C     --------------------------------------------------------------

      call ievalnalu(n,xp,m,lambda,rho,equatn,linear,.true.,hp,inform)
      if ( inform .lt. 0 ) return

C     --------------------------------------------------------------
C     Compute gradients difference
C     --------------------------------------------------------------

      do i = 1,n
          hp(i) = ( hp(i) - g(i) ) / step
      end do

C     --------------------------------------------------------------
C     Add contribution of linear constraints
C     --------------------------------------------------------------

      do j = 1,m

          if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and. 
     +         linear(j) ) then

C             Compute inner product <a,p>
              atp = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  atp = atp + jcval(i) * p(jcvar(i))
              end do

              atp = atp * rho(j)

C             Add rho * atp * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
              end do

          end if

      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine coo2csr(m,nnz,alin,acol,aval,alen,asta)

      implicit none

C     SCALAR ARGUMENTS
      integer m,nnz

C     ARRAY ARGUMENTS
      integer acol(nnz),alen(m),alin(nnz),asta(m)
      double precision aval(nnz)

C     This subroutines converts a matrix from coordinate format to
C     compressed sparse row format.

C     LOCAL SCALARS
      integer i,j,col,coltmp,lin,lintmp
      double precision val,valtmp

      do i = 1,m
          alen(i) = 0
      end do

      do i = 1,nnz
          lin = alin(i)
          alen(lin) = alen(lin) + 1
      end do

      asta(1) = 1
      do i = 2,m
          asta(i) = asta(i-1) + alen(i-1)
      end do

      do i = 1,nnz

          val = aval(i)
          col = acol(i)
          lin = alin(i)

          alin(i) = - 1

 10       if ( lin .ge. 0 ) then

              j = asta(lin)
              asta(lin) = j + 1

              valtmp = aval(j)
              coltmp = acol(j)
              lintmp = alin(j)

              aval(j) = val
              acol(j) = col
              alin(j) = - 1

              val = valtmp
              col = coltmp
              lin = lintmp

              go to 10

          end if

      end do

      do i = 1,m
          asta(i) = asta(i) - alen(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      logical function sstop(n,x,m,lambda,rho,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

#include "algparam.com"

C     EXTERNAL FUNCTIONS
      logical minsqstop

      if ( innercall ) then
          sstop = minsqstop(n,x,inform)
          return
      end if

      end
