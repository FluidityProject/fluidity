C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/
