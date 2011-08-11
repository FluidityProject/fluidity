from math import sin, cos, pi
import numpy

theta = 0.5
h = 1.0
g = 1.0
d0 = 1.0

def eta_src(X, t):
  x = X[0]
  y = X[1]
  return   2*pi*d0*cos(2*(t + y)*pi) + 2*pi*d0*cos(2*(t + x)*pi) - 2*pi*cos(2*(t + y)*pi) - 2*pi*cos(2*(t + x)*pi)

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([-2*pi*g*cos(2*(t + x)*pi) + 2*pi*cos(2*(t + x)*pi),-2*pi*g*cos(2*(t + y)*pi) + 2*pi*cos(2*(t + y)*pi),0])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([sin(2*(t + x)*pi), sin(2*(t + y)*pi), 0.0])

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  return -sin(2*(t + y)*pi) - sin(2*(t + x)*pi)
