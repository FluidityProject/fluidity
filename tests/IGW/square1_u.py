def val(X,t):
   from math import exp
   print "Hello, world"
   x = X[0]
   y = X[1]
   g = [None, None]
   g[0] =  0.555555555556*(y - 0.5)*exp(2.77777777778*(-(y - 0.5)**2 - (x - 0.5)**2))
   g[1] = -0.555555555556*(x - 0.5)*exp(2.77777777778*(-(y - 0.5)**2 - (x - 0.5)**2))
   return g
