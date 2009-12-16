def val(X,t):
   from math import exp
   X0 = 5.0 - t
   sig = 1.0
   g = [0.0,0.0]
   g[0] = exp(-X[1]/0.1)*exp(-(X[0]-X0)**2/(sig*sig))
   return g
