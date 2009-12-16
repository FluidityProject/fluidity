C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/
