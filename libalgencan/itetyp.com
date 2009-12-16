C     COMMON SCALARS
      logical sameface
      integer ittype

C     COMMON BLOCKS
      common /itedat/ sameface,ittype
      save   /itedat/
