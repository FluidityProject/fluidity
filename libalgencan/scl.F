C     sclinfo:
C
C     0: Success.
C     1: Invalid row or column index.
C     2: Duplicate entry.
C     6: Insufficient space to store the input matrix.
C     7: Insufficient double precision working space.

C     ******************************************************************
C     ******************************************************************

      logical function scl(sclsub)

C     SCALAR ARGUMENTS
      character * 4 sclsub

      scl = .false.
      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclini()

      implicit none

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclana(nsys,hnnz,hlin,hcol,hval,hdiag,sclinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,sclinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclsys(nsys,hnnz,hlin,hcol,hval,s,sclinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,sclinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)
      double precision hval(hnnz),s(nsys)

      end
