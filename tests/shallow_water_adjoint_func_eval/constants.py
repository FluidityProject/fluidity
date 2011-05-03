from math import sin, cos, pi
import numpy

d0 = 1.0
theta = 0.5
h = 1.0
functional = 0.5
dfunctional = 1

def eta_src(x, t):
  return 0.0

def u_src(x, t):
  return numpy.array([0.0,0,0])

def u_exact(x, t):
  return numpy.array([sin(2*pi*(x+t)), 0.0, 0.0])

def eta_exact(x, t):
  return h-sin(2*pi*(x+t))
