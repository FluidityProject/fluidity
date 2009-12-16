C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/
