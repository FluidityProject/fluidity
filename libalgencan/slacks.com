C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/
