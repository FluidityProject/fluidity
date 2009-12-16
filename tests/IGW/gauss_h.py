def val(X,t):
   from math import exp, pi, cos
   x = X[0]
   y = X[1]
   rsq = (X[0]-0.5)**2 + (X[1]-0.5)**2
   r0sq = 0.3**2
   g = exp(-rsq/r0sq)
   if x>0.999 or x<0.001 or y>0.999 or y<0.001:
      g = 0.0
   
   return g
