C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/
