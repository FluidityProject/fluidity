C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
