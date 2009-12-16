C     ******************************************************************
C     ******************************************************************

      program algencanma

      implicit none

#include "dim.par"

C     LOCAL SCALARS
      logical checkder
      integer inform,iprint,m,n,ncomp
      double precision cnorm,epsfeas,epsopt,f,nlpsupn,snorm

C     LOCAL ARRAYS
      logical coded(10),equatn(mmax),linear(mmax)
      double precision l(nmax),lambda(mmax),u(nmax),x(nmax)

C     EXTERNAL SUBROUTINES
      external algencan,param

C     SET UP PROBLEM DATA

      call param(epsfeas,epsopt,iprint,ncomp)

      call inip(n,x,l,u,m,lambda,equatn,linear,coded,checkder)

      call algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
     +linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)

      call endp(n,x,l,u,m,lambda,equatn,linear)

      stop

      end

C     ******************************************************************
C     ******************************************************************

      subroutine param(epsfeas,epsopt,iprint,ncomp)

C     SCALAR ARGUMENTS
      integer iprint,ncomp
      double precision epsfeas,epsopt

      epsfeas  = 1.0d-08
      epsopt   = 1.0d-08

      iprint   = 10
      ncomp    = 6

      end

