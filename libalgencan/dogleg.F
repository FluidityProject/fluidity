C     ******************************************************************
C     ******************************************************************

      subroutine dogleg(n,g,hnnz,hlin,hcol,hval,pd,delta,p,dlinfo)
      
      implicit none

C     SCALAR ARGUMENTS
      logical pd
      integer dlinfo,hnnz,n
      double precision delta
      
C     ARRAY ARGUMENTS
      integer hcol(hnnz),hlin(hnnz)
      double precision g(n),hval(hnnz),p(n)

C     Compute approximate minimizer of the quadratic model using Dogleg
C     method when the Hessian is positive definite and Cauchy point
C     otherwise.

C     dlinfo:
C
C     0: successfull exit. Both H and g are null;
C     1: successfull exit.

#include "dim.par"
#include "machconst.com"
#include "outtyp.com"

C     LOCAL SCALARS
      integer col,i,lin
      double precision a,b,c,coef,d,delta2,geucn,geucn2,gthg,pbeucn2,
     +        pueucn,pueucn2,putpb,r

C     LOCAL ARRAYS
      double precision pu(nmax)
C     Presentation

      if ( iprintinn .ge. 5 ) then
          write(* ,1000)
          write(file10_unit,1000)
      end if

C     Initialization

      dlinfo = 1
      delta2 = delta**2

C     If H is not positive definite, compute Cauchy point

      if ( .not. pd ) then 
          go to 100
      end if
      
      pbeucn2 = 0.0d0
      do i = 1,n
          pbeucn2 = pbeucn2 + p(i)**2
      end do
      
C     If Newton step is inside the trust region, this step is taken
      
      if ( pbeucn2 .le. delta2 ) then
          go to 500
      end if
      
C     If Newton step is outside the trust region, compute the
C     unconstrained minimizer pu of the quadratic function
      
C     Compute g^T H g
      
      do i = 1,n
          pu(i) = 0.0d0
      end do
      
      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          
          pu(lin) = pu(lin) + hval(i) * g(col)
          if ( lin .ne. col ) then
              pu(col) = pu(col) + hval(i) * g(lin)
          end if
      end do
      
      gthg = 0.0d0
      do i = 1,n
          gthg = gthg + g(i) * pu(i)
      end do
      
C     Compute g^T g.
      
      geucn2 = 0.0d0
      do i = 1,n
          geucn2 = geucn2 + g(i)**2
      end do
      
C     Compute pu = coef * g.
      
      coef = - geucn2 / gthg 
      do i = 1,n
          pu(i) = coef * g(i)
      end do
      
C     If uncontrained minimizer is outside the trust region, it is
C     truncated at the border
      
      pueucn2 = 0.0d0
      do i = 1,n
          pueucn2 = pueucn2 + pu(i)**2
      end do
      pueucn = sqrt( pueucn2 )
      
      if ( pueucn2 .ge. delta2 ) then
          r = delta / pueucn
          do i = 1,n
              p(i) = r * pu(i)
          end do
          go to 500
      end if
      
C     Compute step length in directions pu and pb. Direction p is a
C     linear combination of pu and pb. To compute the step length in
C     each direction, we have to solve (in r):
C     \| pu + (r - 1) (pb - pu) \|^2 = delta^2.
      
      putpb = 0.0d0
      do i = 1,n
          putpb = putpb + pu(i) * p(i)
      end do
      
      a = pbeucn2 + pueucn2 - 2.0d0 * putpb
      b = - 2.0d0 * pbeucn2 - 4.0d0 * pueucn2 + 6.0d0 * putpb
      c = pbeucn2 + 4.0d0 * pueucn2 - 4.0d0 * putpb - (delta**2)
      d = b**2 - 4.0d0 * a * c
      
      if ( ( 2.0d0 * abs( a ) .lt. macheps23 ) .or. 
     +     ( d .lt. 0.0d0 ) ) then
          go to 200
      end if
      
      r = (- b + sqrt( d )) / (2.0d0 * a)
      
      if ( ( r .lt. 1.0d0 ) .or. ( r .gt. 2.0d0 ) ) then
          r = (- b - sqrt( d )) / (2.0d0 * a)
      end if
      r = r - 1.0d0
      
      do i = 1,n
          p(i) = r * p(i) + (1.0d0 - r) * pu(i)
      end do
      
      go to 500
      
C     Compute Cauchy point, when H is not positive definite

 100  continue
      
C     Compute g^T H g

      do i = 1,n
          pu(i) = 0.0d0
      end do
      
      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          
          pu(lin) = pu(lin) + hval(i) * g(col)
          if ( lin .ne. col ) then
              pu(col) = pu(col) + hval(i) * g(lin)
          end if
      end do
      
      gthg = 0.0d0
      do i = 1,n
          gthg = gthg + g(i) * pu(i)
      end do
      
C     Compute g^T g
      
      geucn2 = 0.0d0
      do i = 1,n
          geucn2 = geucn2 + g(i)**2
      end do
      geucn = sqrt( geucn2 )
      
C     Compute step length coef
      
 200  if ( ( abs( gthg ) .le. macheps23 ) .and. 
     +     ( geucn .le. macheps23 ) ) then

          dlinfo = 0
          
          if ( iprintinn .ge. 5 ) then
              write(* ,1020)
              write(file10_unit,1020)
          end if

          return
      end if
      
      if ( ( gthg .le. 0.0d0 ) .or. 
     +     ( geucn2*geucn .lt. delta*gthg ) ) then
          coef = - delta / geucn
      else
          coef = - geucn2 / gthg 
      end if

C     Compute p = coef * g
      
      do i = 1,n
          p(i) = coef * g(i)
      end do

C     Termination

 500  continue

      if ( iprintinn .ge. 5 ) then
          write(* ,1010)
          write(file10_unit,1010)
      end if

      if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
          write(*, 1030) min0(n,nprint),(p(i),i=1,min0(n,nprint))
          write(file10_unit,1030) min0(n,nprint),(p(i),i=1,min0(n,nprint))
      end if

 1000 format(/,5X,'Computation of Dogleg direction.')
 1010 format(  5X,'Dogleg computed successfully.')
 1020 format(  5X,'Null direction was computed.')
 1030 format(/,5X,'Dogleg direction (first ',I7,' components): ',
     +       /,1(5X,6(1X,1P,D11.4)))

      end
      
C     ******************************************************************
C     ******************************************************************

C     Compute approximate minimizer of the quadratic model using Dogleg
C     method when the Hessian is positive definite and Cauchy point
C     otherwise.

C     On Entry
C
C     n        integer
C              dimension
C
C     g        double precision g(n)
C              vector used to define the quadratic function
C
C     hnnz     integer
C              number of nonzero elements of H
C
C     hlin     integer hlin(hnnz)
C              row indices of nonzero elements of H
C
C     hcol     integer hcol(hnnz)
C              column indices of nonzero elements of H
C
C     hval     double precision hval(hnnz)
C              nonzero elements of H, that is,
C              H(hlin(i),hcol(i)) = hval(i). Since H is symmetric, just
C              one element H(i,j) or H(j,i) must appear in its sparse
C              representation. If more than one pair corresponding to
C              the same position of H appears in the sparse
C              representation, the multiple entries will be summed.
C              If any hlin(i) or hcol(i) is out of range, the entry will
C              be ignored
C     
C     pd       logical
C              indicates if the last Cholesky decomposition of moresor 
C              was successfull. That is, if the last matrix used by 
C              moresor was positive definite
C     
C     delta    double precision
C              trust-region radius
C     
C     On Return
C     
C     p        double precision p(n)      
C              solution to problem
C              minimize     psi(w)
C              subjected to ||w|| <= delta
C      
C     dlinfo   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C
C              0 = successfull exit. Both H and g are null;
C
C              1 = successfull exit.
