C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/
