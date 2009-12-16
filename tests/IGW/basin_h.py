def val(X,t):
   from math import exp
   r0 = 0.5
   r = ((X[0])**2 + (X[1]**2))**0.5
   if(r>0.0001):
      return exp((r-r0)/0.1)*X[0]/r
   else:
      return 0.0
