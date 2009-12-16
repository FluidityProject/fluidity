def val(X,t):
   from math import exp
   r0 = 0.5
   g = [0.0,0.0]
   r = ((X[0])**2 + (X[1])**2)**0.5
   if(r>0.0001):
      g[0] = -exp((r-r0)/0.1)*X[0]/r*X[1]/r
      g[1] = exp((r-r0)/0.1)*X[0]/r*X[0]/r
   return g