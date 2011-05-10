from math import sin, cos, pi
import numpy

d0 = 1.0
theta = 0.5
h = 1.0
g = 1.0
functional = 2
dfunctional = 4

def eta_src(X, t):
  x = X[0]
  y = X[1]
  return sin(2*pi*X[0]) 

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([0,0,0])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([sin(2*pi*x) + t, 0.0, 0.0])

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  return 1.0

def functional_vector_eta(X,t):
  x = X[0]
  y = X[1]
  #return sin(2*pi*x)  
  return 1.0 

def functional_vector_u(X,t):
  x = X[0]
  y = X[1]
  return [0.0, 0.0, 0.0]
