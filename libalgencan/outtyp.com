C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint,file10_unit

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint
      common /outdat/ mprint,file10_unit
      save   /outdat/

