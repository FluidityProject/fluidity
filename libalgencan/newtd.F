C     ******************************************************************
C     ******************************************************************

      subroutine newtd(nind,x,l,u,g,m,rho,equatn,d,memfail,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical memfail
      integer inform,nind,m

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision d(nind),g(*),l(*),rho(m),u(*),x(*)

C     This subroutine solves the Newtonian system
C
C             ( H + rho A^T A ) x = b
C
C     by solving the Martinez-Santos system
C
C             H x + A^t y = b
C             A x - y/rho = 0

#include "dim.par"
#include "machconst.com"
#include "algparam.com"
#include "outtyp.com"
#include "itetyp.com"

C     PARAMETERS
      integer dimmax
      parameter ( dimmax = nmax + mmax )

C     COMMON ARRAYS
      double precision mindiag(nmax)

C     LOCAL SCALARS
      integer dim,hnnz,i,iter,lssinfo,nneigv,pind
      double precision adsupn,diff,pval

C     LOCAL ARRAYS
      integer hdiag(dimmax),hlin(hnnzmax),hcol(hnnzmax)
      double precision adddiag(dimmax),hval(hnnzmax),sol(dimmax)

C     COMMON BLOCKS
      common /diadat/ mindiag
      save   /diadat/

C     ------------------------------------------------------------------
C     Presentation
C     ------------------------------------------------------------------

      if ( iprintinn .ge. 5 ) then
          write(* ,1000)
          write(file10_unit,1000)
      end if

C     ------------------------------------------------------------------
C     Initialization
C     ------------------------------------------------------------------

      iter    = 0
      memfail = .false.

      call lssini(sclsys,.true.,.false.)

C     ------------------------------------------------------------------
C     Compute ML matrix
C     ------------------------------------------------------------------

      call mlsyst(nind,x,g,m,rho,equatn,hlin,hcol,hval,hnnz,hdiag,sol,
     +dim,inform)

      if ( inform .lt. 0 ) return

C     ------------------------------------------------------------------
C     Analyse sparsity pattern
C     ------------------------------------------------------------------

      call lssana(dim,hnnz,hlin,hcol,hval,hdiag,lssinfo)

      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      end if

C     ------------------------------------------------------------------
C     Main loop
C     ------------------------------------------------------------------

 100  continue

      iter = iter + 1

C     ------------------------------------------------------------------
C     Compute regularization
C     ------------------------------------------------------------------

      if ( iter .eq. 1 ) then

          if ( sameface .and. ittype .eq. 3 ) then

              do i = 1,nind
                  mindiag(i) = 0.1d0 * mindiag(i)
              end do

          else

              do i = 1,nind
                  if ( g(i) .eq. 0.0d0 ) then
                      mindiag(i) = 0.0d0
                  else
                      if ( g(i) .gt. 0.0d0 ) then
                          diff = x(i) - l(i)
                      else if ( g(i) .lt. 0.0d0 ) then
                          diff = u(i) - x(i)
                      end if
                      mindiag(i) = abs( g(i) / diff )
                  end if
              end do

          end if

      else

          do i = 1,nind
              if ( mindiag(i) .eq. 0.0d0 ) then
                  mindiag(i) = macheps23
              else
                  mindiag(i) = 10.0d0 * mindiag(i)
              end if
          end do

      end if

      do i = 1,nind
          adddiag(i) = max( macheps23, mindiag(i) - hval(hdiag(i)) )
      end do

      do i = nind + 1,dim
          adddiag(i) = 0.0d0
      end do

      adsupn = 0.0d0
      do i = 1,dim
          adsupn = max( adsupn, adddiag(i) )
      end do

      if ( iprintinn .ge. 5 ) then
          write(* ,1010) adsupn
          write(file10_unit,1010) adsupn
      end if

C     ------------------------------------------------------------------
C     Factorize matrix
C     ------------------------------------------------------------------

      call lssfac(dim,hnnz,hlin,hcol,hval,hdiag,adddiag,pind,pval,
     +nneigv,lssinfo)

      if ( lssinfo .eq. 0 .or. lssinfo .eq. 1 ) then

          if ( nneigv .ne. dim - nind ) then
C           ! WRONG INERTIA (SEE NOCEDAL AND WRIGHT)

C             Lemma 16.3 [pg. 447]: Assume that the Jacobian of the
C             constraints has full rank and that the reduced Hessian
C             Z^T H Z is positive definite. Then the Jacobian of the
C             KKT system has n positive eigenvalues, m negative
C             eigenvalues, and no zero eigenvalues.

C             Note that at this point we know that the matrix has no
C             zero eigenvalues. nneigv gives the number of negative
C             eigenvalues.

              if ( iprintinn .ge. 5 ) then
                  write(* ,1020) nneigv,dim - nind
                  write(file10_unit,1020) nneigv,dim - nind
              end if

              go to 100

          else

              if ( iprintinn .ge. 5 ) then
                  write(* ,1030)
                  write(file10_unit,1030)
              end if

          end if

      else if ( lssinfo .eq. 2 ) then
        ! SINGULAR JACOBIAN

          if ( iprintinn .ge. 5 ) then
              write(* ,1040)
              write(file10_unit,1040)
          end if

          go to 100

      else if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else ! if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if

C     ------------------------------------------------------------------
C     Solve
C     ------------------------------------------------------------------

      call lsssol(dim,sol)

      do i = 1,nind
          d(i) = sol(i)
      end do

      if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
          write(*, 1050) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
          write(file10_unit,1050) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
      end if

C     NON-EXECUTABLE STATEMENTS

 1000 format(/,5X,'Sparse factorization of the ML system.')
 1010 format(  5X,'Maximum value added to the diagonal: ',1P,D24.16)
 1020 format(  5X,'ML-matrix with wrong inertia.',
     +       /,5X,'Actual number of negative eigenvalues  = ',I16,'.',
     +       /,5X,'Desired number of negative eigenvalues = ',I16,'.')
 1030 format(  5X,'Direct solver finished successfully.')
 1040 format(  5X,'ML-matrix numerically singular.')
 1050 format(/,5X,'Newton direction (first ',I7,' components): ',
     +       /,1(5X,6(1X,1P,D11.4)))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine mlsyst(nind,x,nal,m,rho,equatn,ulin,ucol,uval,unnz,
     +udiag,b,dim,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer dim,inform,m,nind,unnz

C     ARRAY ARGUMENTS
      logical equatn(m)
      integer ucol(*),udiag(*),ulin(*)
      double precision b(*),nal(nind),rho(m),uval(*),x(*)

#include "dim.par"
#include "graddat.com"
#include "rspace.com"

C     LOCAL SCALARS
      integer col,i,j,k,lin,var

C     LOCAL ARRAYS
      integer wi(nmax)

C     This subrotuine is called from the reduced space.

C     MATRIX

C     Compute Hessian of the Lagrangian

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

      call expand(nind,x)

      call sevalhl(nt,x,m,dpdc,ulin,ucol,uval,unnz,inform)
      if ( inform .lt. 0 ) return

      call shrink(nind,x)

C     Preparation for shrink (wi indicates, for each free variable x_i,
C     its rank within the set of free variables. wi(i)=0 if x_i is not
C     a free variable)

      do i = 1,nt
         wi(i) = 0
      end do

      do i = 1,nind
         wi(ind(i)) = i
      end do

C     Shrink Hessian of the Lagrangian and set diagonal-elements indices

      k = 0

      do i = 1,nind
          udiag(i) = 0
      end do

      do i = 1,unnz
          lin = wi(ulin(i))
          col = wi(ucol(i))

          if ( lin .ne. 0 .and. col .ne. 0 ) then
              k = k + 1
              ulin(k) = lin
              ucol(k) = col
              uval(k) = uval(i)

              if ( lin .eq. col ) then
                  udiag(lin) = k
              end if
          end if
      end do

      do i = 1,nind
          if ( udiag(i) .eq. 0 ) then
              k = k + 1
              ulin(k) = i
              ucol(k) = i
              uval(k) = 0.0d0

              udiag(i) = k
          end if
      end do

C     Shrink Jacobian and add diagonal matrix - 1.0 / rho

      dim = nind

      do j = 1,m
          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

              dim = dim + 1

              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  var = wi(jcvar(i))

                  if ( var .ne. 0 ) then
                      k = k + 1
                      ulin(k) = dim
                      ucol(k) = var
                      uval(k) = jcval(i)
                  end if
              end do

              k = k + 1
              ulin(k) = dim
              ucol(k) = dim
              uval(k) = - 1.0d0 / rho(j)

              udiag(dim) = k
          end if
      end do

      unnz = k

C     RHS

      do i = 1,nind
          b(i) = - nal(i)
      end do

      do i = nind + 1,dim
          b(i) = 0.0d0
      end do

      end
