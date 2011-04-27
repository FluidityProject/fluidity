from math import sin, cos, pi
import numpy
def eta_source(x, t):
  return -2*(t - 10)*pi*t*cos(2*pi*x) - (t - 5)*cos(2*pi*x) - t*cos(2*pi*x)

def u_source(x, t):
  return numpy.array([2*(t - 5)*pi*t*sin(2*pi*x) - (t - 10)*sin(2*pi*x) - t*sin(2*pi*x),0,0])

def u_init(x):
  return numpy.array([0.0, 0.0, 0.0])
def eta_init(x):
  return 0.0
