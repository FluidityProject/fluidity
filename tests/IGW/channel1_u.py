def val(X,t):
   from math import exp
   g = [0.0,0.0]
   g[0] = 0.0
   g[1] = exp(-X[0]/0.1)*exp(-(X[1]-5.0)**2)
   return g