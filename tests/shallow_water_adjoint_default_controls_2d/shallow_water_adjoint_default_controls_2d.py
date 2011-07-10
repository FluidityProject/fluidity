from math import sin, cos, pi
import numpy

d0 = 0.5
theta = 0.5
h = 1.0
g = 0.1
functional = 2
dfunctional = 4

def eta_src(X, t):
  x = X[0]
  y = X[1]
  return 2*(t + 1)*pi*d0*h*cos(2*pi*x) - 2*pi*d0*sin(2*pi*y) + h*cos(2*pi*x)

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([ -2*(t + 1)*pi*g*h*sin(2*pi*x) + h*sin(2*pi*x),0,0])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([(t + 1)*h*sin(2*pi*x), cos(2*pi*x) + cos(2*pi*y), 0])

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  return (t + 1)*h*cos(2*pi*x)
