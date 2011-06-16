from math import sin, cos, pi
import numpy

theta = 0.5
h = 1.0
g = 1.0
d0 = 1.0
r = 6371220.0

def eta_src(X, t):
  x = X[0]
  y = X[1]
  return 0

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([g/r, -g*x*sin(y)/r])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([0, 0])

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  return sin(x)
