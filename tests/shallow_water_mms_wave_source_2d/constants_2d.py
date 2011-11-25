from math import sin, cos, pi
import numpy

d0 = 0.1
theta = 0.5
h = 1.0
g = 0.43
functional = 2
dfunctional = 4

def eta_src(X, t):
  x = X[0]
  y = X[1]
  return -2*pi*d0*sin(2*pi*y + t) + 2*pi*d0*cos(2*pi*x + t) - cos(2*pi*x + t)

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([-2*pi*g*cos(2*pi*x + t) + cos(2*pi*x + t), -sin(2*pi*y + t), 0])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([sin(2*pi*x + t), cos(2*pi*y + t), 0.0])

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  return -sin(2*pi*x + t)
