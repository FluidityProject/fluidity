C     lssinfo:
C
C     0: Success.
C     1: Matrix not positive definite.
C     2: Rank deficient matrix.
C     6: Insufficient space to store the linear system.
C     7: Insufficient double precision working space.
C     8: Insufficient integer working space.

C     ******************************************************************
C     ******************************************************************

      logical function lss(lsssub)

C     SCALAR ARGUMENTS
      character * 4 lsssub

      lss = .false.
      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssini(sclsys,acneig,usefac)

      implicit none

C     SCALAR ARGUMENTS
      logical acneig,sclsys,usefac

C     All parameters of this routine are entry parameters:
C
C     sclsys: set it equal to TRUE to scale the linear system.

C     acneig: set it equal to TRUE to stop the factorization if the
C             matrix of coefficients is not positive definite.

C     usefac: sset it equal to TRUE if the Cholesky factor is 
C             explicitely needed.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssana(nsys,hnnz,hlin,hcol,hval,hdiag,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,lssinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssfac(nsys,hnnz,hlin,hcol,hval,hdiag,d,pind,pval,
     +nneigv,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,pind,nneigv,lssinfo
      double precision pval

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz),d(nsys)

C     This routine computes the factorization of a matrix H + D. H is
C     given in the coordinate format by the arrays hlin, hcol and hval, 
C     all with size hnnz. hdiag(i) stores the position in hval that 
C     holds the i-th diagonal elements of H (if such element is null 
C     and is not present in the coordinate representation, hdiag(i) 
C     should be set to 0). The diagonal matrix D is representend by an 
C     array of nsys elements.
C
C     On return, if the user chose to stop the factorization whenever 
C     the matrix H + D is not positive definite, pind will hold the 
C     pivot index where the factorization failed and pval will hold the
C     pivot value. nneigv will hold the number of negative eigenvalues
C     of H + D.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssol(nsys,sol)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      double precision sol(nsys)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssoltr(job,nsys,sol)

      implicit none

C     SCALAR ARGUMENTS
      character * 1 job
      integer nsys

C     ARRAY ARGUMENTS
      double precision sol(nsys)

C     Given the factorization R^T R = P(H + D)P^T, this routine will 
C     solve: P^T R^T x = b, if job = 'T' or 't' and R P x = b otherwise.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsspermvec(nsys,v)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      integer v(nsys)

C     Permute integer array v using the same permutation used in the 
C     factorization.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssunpermvec(nsys,v)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      integer v(nsys)

C     Permute integer array v using the inverse permutation used in the 
C     factorization.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsspermind(hnnz,hlin,hcol)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)

C     Permute the contents of vectors hlin and hcol using the same 
C     permutation used in the factorization.

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssunpermind(nsys,hnnz,hlin,hcol)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)

C     Permute the contents of vectors hlin and hcol using the inverse
C     permutation used in the factorization.

      end

C     ******************************************************************
C     ******************************************************************

      double precision function lssgetd(j)

      implicit none

C     SCALAR ARGUMENTS
      integer j

C     Returns the j-th element of the diagonal of factor R (where R is 
C     such that R^T R = P(H + D)P^T).

      lssgetd = 0.0d0

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssetrow(nsys)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssgetrow(nsys,idx,rownnz,rowind,rowval)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,idx,rownnz

C     ARRAY ARGUMENTS
      integer rowind(nsys)
      double precision rowval(nsys)

C     Stores in a sparse vector, represented by arrays rowind and 
C     rowval, the idx-th row of factor R (where R is such that R^T R = 
C     P(H + D)P^T).

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssafsol(nsys,hnnz,hlin,hcol,hval,hdiag,d,sol,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,lssinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz),d(nsys),sol(nsys)

C     Solves a linear system (H' + D')x = b. 

      end
