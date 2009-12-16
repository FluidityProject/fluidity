!*---------------------------------------------------------------------*
!*These lines of code contain the parameter statements for the global 
!*definitions of many important parameters. 
!*---------------------------------------------------------------------*  
      real tol,rk,rkd,t0
      integer imax,jmax,mx,mxpoi,mxele,mxbou,nd
      parameter ( imax=10, jmax=10)
      parameter ( mx=100, mxpoi=mx, mxele=200, mxbou=mx/5, nd=4 )
      parameter ( tol=1.0e-6, rk=0.1, rkd = 0.0 )
      parameter ( T0=0.3024 )
      integer ntimemax
      parameter (ntimemax = 10)
      integer istate
      integer subistate
      integer nvar
      parameter (nvar=4)
      parameter ( istate = mxpoi*nvar)
      parameter ( subistate = mxpoi)
      integer nrsnapshots 
      parameter (nrsnapshots = 3 )
      integer SnapNdT
      parameter (SnapNdT=5 )
