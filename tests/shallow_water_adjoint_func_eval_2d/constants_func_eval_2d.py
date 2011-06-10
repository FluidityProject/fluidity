from math import sin, cos, pi
import numpy

d0 = 0.9
theta = 0.5
g = 0.8
functional = 2
dfunctional = 4

def eta_src(X, t):
  x = X[0]
  y = X[1]
  return 2*pi*d0*g*cos(2*(t + y)*pi) + 2*pi*d0*g*cos(2*(t + x)*pi) - 2*pi*cos(2*(t + y)*pi) - 2*pi*cos(2*(t + x)*pi)

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([0,0,0])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([g*sin(2*(t + x)*pi),  g*sin(2*(t + y)*pi), 0.0])

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  return -sin(2*(t + y)*pi) - sin(2*(t + x)*pi)

def functional_vector_eta(X,t):
  x = X[0]
  y = X[1]
  return sin(2*pi*(x+t))  

def functional_vector_u(X,t):
  x = X[0]
  y = X[1]
  if t == 0:
    return [sin(2*pi*x), sin(2*pi*x), 0.0]
  else:  
    # Note: The velocity source term currently has to be zero, because the cartesian <-> local projection methods are not implemented in python.
    #       Therefore, the dot product of the functional_vector_u will always be zero, so we just choose functional_vector_u = 0 here.
    return [0.0, 0.0, 0.0]
