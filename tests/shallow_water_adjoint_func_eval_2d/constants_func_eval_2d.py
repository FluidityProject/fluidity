from math import sin, cos, pi
import numpy

d0 = 1.0
theta = 0.5
h = 1.0
g = 9.81
functional = 2
dfunctional = 4

def eta_src(X, t):
  return 2*pi*d0*cos(2*pi*X[0] + t) - cos(2*pi*X[0] + t)

def u_src(X, t):
  return numpy.array([-g*2*pi*cos(2*pi*X[0] + t) + cos(2*pi*X[0] + t),0,0])

def u_exact(X, t):
  return numpy.array([sin(2*pi*X[0] + t), 0.0, 0.0])

def eta_exact(X, t):
  return -sin(2*pi*X[0] + t)

def functional_vector_eta(X,t):
  return X[0]*(X[0]-1)+t

def functional_vector_u(X,t):
  return [sin(2*pi*X[0])+t, 0.0, 0.0]
