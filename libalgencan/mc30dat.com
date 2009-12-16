C     COMMON ARRAYS
      double precision w(4*nsysmax)

C     COMMON BLOCKS
      common /scldat/ w
      save   /scldat/
