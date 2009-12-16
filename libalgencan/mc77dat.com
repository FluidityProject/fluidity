C     COMMON ARRAYS
      integer icntl(10),iw(nsysmax),info(10)
      double precision cntl(10),dw(hnnzmax),rinfo(10)

C     COMMON BLOCKS
      common /scldat/ cntl,dw,rinfo,icntl,iw,info
      save   /scldat/
