C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
