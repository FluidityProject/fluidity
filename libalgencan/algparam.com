C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
