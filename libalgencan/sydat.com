C     COMMON SCALARS
      double precision seucn,sts,sty,yeucn

C     COMMON ARRAYS
      double precision s(nmax),y(nmax)

C     COMMON BLOCKS
      common /sydata/ s,y,seucn,yeucn,sts,sty
      save   /sydata/
