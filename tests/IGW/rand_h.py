def val(X,t):
   from random import random

   x = X[0]
   y = X[1]
   
   g = random()
   if x>0.999 or x<0.001 or y>0.999 or y<0.001:
      g = 0.0
   
   return g
