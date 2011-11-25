from math import sin, cos, pi
import numpy

theta = 0.5
h = 1.0
g = 1.0
d0 = 1.0

def eta_src(x, t):
  return  2*pi*d0*h*cos(2*(t + x)*pi) - 2*pi*h*cos(2*(t + x)*pi) 

def u_src(x, t):
  return numpy.array([-2*pi*g*h*cos(2*(t + x)*pi) + 2*pi*h*cos(2*(t + x)*pi),0,0])

def u_exact(x, t):
  return numpy.array([h*sin(2*(t + x)*pi), 0.0, 0.0])

def eta_exact(x, t):
  return -h*sin(2*(t + x)*pi)
