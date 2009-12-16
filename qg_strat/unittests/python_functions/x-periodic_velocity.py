def val(X,t):
   from math import cos, sin, pi
   L=1.0e+6
   k=(2.0*pi)/L
   l=(2.0*pi)/L
   return [-l*cos(k*X[0])*cos(l*X[1]), -k*sin(k*X[0])*sin(l*X[1]), 0.0]
