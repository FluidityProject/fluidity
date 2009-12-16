C     COMMON SCALARS
      double precision plspg,psmdyty

C     COMMON ARRAYS
      double precision pdiag(nmax),psmdy(nmax)

C     COMMON BLOCKS
      common /hpredata/ pdiag,psmdy,plspg,psmdyty
      save   /hpredata/
