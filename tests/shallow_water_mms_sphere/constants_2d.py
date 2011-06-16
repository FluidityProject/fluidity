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
  return d0*cos(1.0/1000*pi*t + x)/r

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([1.0/1000*pi*cos(1.0/1000*pi*t + x)*cos(y), -g*sin(y)/r])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([cos(y)*sin(x+2*pi*t/2000), 0])

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  return cos(y)
