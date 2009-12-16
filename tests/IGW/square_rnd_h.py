def val(X,t):
   from random import random
   if X[0]<0.0001 or X[0]>0.9999 or X[1]<0.0001 or X[1]>0.9999:
      return 0.0
   else:
      return random()
