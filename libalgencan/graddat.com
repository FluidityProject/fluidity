C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

