C     COMMON SCALARS
      double precision hlspg,hstds

C     COMMON ARRAYS
      double precision hds(nmax)

C     COMMON BLOCKS
      common /happdata/ hlspg,hds,hstds
      save   /happdata/
