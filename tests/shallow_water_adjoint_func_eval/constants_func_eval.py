from math import sin, cos, pi
import numpy

d0 = 5.5
theta = 0.5
h = 1.0
g = 9.81
functional = 0.5
dfunctional = 1


def eta_src(X, t):
  x = X[0]
  return 2*pi*d0*cos(2*pi*x)

def u_src(X, t):
  x = X[0]
  return numpy.array([2*pi*g*cos(2*pi*x)+1,0,0])

def u_exact(X, t):
  x = X[0]
  return numpy.array([sin(2*pi*x)+t, 0.0, 0.0])

def eta_exact(X, t):
  x = X[0]
  return sin(2*pi*x)

def functional_vector_eta(X,t):
  x = X[0]
  return X[0]*(X[0]-1)+t

def functional_vector_u(X,t):
  x = X[0]
  return [sin(2*pi*X[0])+t, 0.0, 0.0]
