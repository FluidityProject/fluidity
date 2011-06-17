from math import sin, cos, pi
import numpy

theta = 0.5
h = 1.0
g = 1.0
d0 = 1.0
r = 6371220.0

K = 0.242334 # equilibrium amplitude of tidal constituent
sigma = 1.40519e-04 # frequency of tidal constituent

def eta_src(X, t):
  x = X[0]
  y = X[1]
  return -K*sigma*sin(sigma*t + 2*x)*cos(y)**2 - 2*K*d0*cos(sigma*t + 2*x)*cos(y)/r

def u_src(X, t):
  x = X[0]
  y = X[1]
  return numpy.array([-K*sigma*cos(sigma*t + 2*x)*cos(y)**2 - 2*K*g*sin(sigma*t + 2*x)*cos(y)/r, -2*K*g*sin(y)*cos(sigma*t + 2*x)*cos(y)/r])

def u_exact(X, t):
  x = X[0]
  y = X[1]
  # We prescribe the equilibrium tide response for the M2 constituent
  return [-cos(y)**2*K*sin(sigma*t+2*x), 0] # only considring semodiurnal constituents 

def eta_exact(X, t):
  x = X[0]
  y = X[1]
  # We prescribe the equilibrium tide response for the M2 constituent
  return cos(y)**2*K*cos(sigma*t+2*x) # only considring semodiurnal constituents 
