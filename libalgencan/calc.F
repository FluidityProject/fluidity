C     ******************************************************************
C     ******************************************************************

      subroutine calcal(nind,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,nind
      double precision al

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(*)

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Compute augmented Lagrangian

      call sevalal(nt,x,m,lambda,rho,equatn,linear,al,inform)
      if ( inform .lt. 0 ) return

C     Shrink x to the reduced space

      call shrink(nind,x)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calcnal(nind,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,nind

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(*),rho(m),x(*)

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Compute the gradient of the augmented Lagrangian

      call sevalnal(nt,x,m,lambda,rho,equatn,linear,nal,inform)
      if ( inform .lt. 0 ) return

C     Shrink x and nal to the reduced space

      call shrink(nind,x)
      call shrink(nind,nal)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calchal(nind,x,m,lambda,rho,equatn,linear,hlin,hcol,
     +hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,nind

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer hcol(*),hlin(*)
      double precision lambda(m),hval(*),rho(m),x(*)

C     This subroutine computes the Hessian of the augmented Lagrangian
C     in the reduced space.

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer col,i,k,lin

C     LOCAL ARRAYS
      integer wi(nmax)

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Compute the Hessian of the augmented Lagrangian

      call sevalhal(nt,x,m,lambda,rho,equatn,linear,hlin,hcol,hval,hnnz,
     +inform)
      if ( inform .lt. 0 ) return

C     Shrink x to the reduced space

      call shrink(nind,x)

C     Shrink representation of H
      
      do i = 1,nt
         wi(i) = 0
      end do

      do i = 1,nind
         wi(ind(i)) = i
      end do
      
      k = 0

      do i = 1,hnnz
          lin = wi(hlin(i))
          col = wi(hcol(i))

          if ( lin .ne. 0 .and. col .ne. 0 ) then
              k = k + 1
              hlin(k) = lin
              hcol(k) = col
              hval(k) = hval(i)
          end if
      end do

      hnnz = k

      end


C     ******************************************************************
C     ******************************************************************

      subroutine calchalp(nind,x,m,lambda,rho,equatn,linear,p,hp,gothl,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,nind

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(*),lambda(m),p(*),rho(m),x(*)

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Complete p with zeroes

      do i = 1,nt - nind
          p(nind+i) = 0.0d0
      end do

C     Expand x and p to the full space

      call expand(nind,x)
      call expand(nind,p)

C     Compute the Hessian-vector product

      call sevalhalp(nt,x,m,lambda,rho,equatn,linear,p,hp,gothl,inform)
      if ( inform .lt. 0 ) return

C     Shrink x, p and hp to the reduced space

      call shrink(nind,x)
      call shrink(nind,p)
      call shrink(nind,hp)
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine capplyhpre(nind,m,rho,equatn,gotp,r,z)

      implicit none

C     SCALAR ARGUMENTS
      logical gotp
      integer m,nind

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision r(*),rho(m),z(*)

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer i

C     Complete r with zeroes

      do i = nind + 1,nt
          r(i) = 0.0d0
      end do

C     Expand r to the full space

      call expand(nind,r)

C     Solve P z = r

      call applyhpre(nt,m,rho,equatn,gotp,r,z)

C     Shrink r and z to the reduced space

      call shrink(nind,r)
      call shrink(nind,z)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine csetp(nind,x)

      implicit none

C     SCALAR ARGUMENTS
      integer nind

C     ARRAY ARGUMENTS
      double precision x(*)

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Set point

      call ssetp(nt,x)

C     Shrink x

      call shrink(nind,x)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine shrink(nind,v)

      implicit none

C     This subroutine shrinks vector v from the full dimension space 
C     (dimension n) to the reduced space (dimension nind).

C     SCALAR ARGUMENTS
      integer nind

C     ARRAY ARGUMENTS
      double precision v(*)

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = 1,nind
           indi = ind(i)
           if ( i .ne. indi ) then
               tmp     = v(indi)
               v(indi) = v(i)
               v(i)    = tmp
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine expand(nind,v)

      implicit none

C     This subroutine expands vector v from the reduced space 
C     (dimension nind) to the full space (dimension n).

C     SCALAR ARGUMENTS
      integer nind

C     ARRAY ARGUMENTS
      double precision v(*)

#include "dim.par"
#include "rspace.com"

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = nind,1,- 1
          indi = ind(i)
          if ( i .ne. indi ) then
              tmp     = v(indi)
              v(indi) = v(i)
              v(i)    = tmp
          end if
      end do
     
      end
 
