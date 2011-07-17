from math import sin, cos, pi
import numpy

d0 = 0.5
theta = 0.5
h = 1.0
g = 0.1
functional = 2
dfunctional = 4

def eta_src(x, t):
  return 2*(t + 1)*d0*pi*h*cos(2*pi*x) + h*cos(2*pi*x)

def u_src(x, t):
  return numpy.array([-2*(t + 1)*pi*g*h*sin(2*pi*x) + h*sin(2*pi*x),0,0])

def u_exact(x, t):
  return numpy.array([(t+1)*h*sin(2*pi*x), 0.0, 0.0])

def eta_exact(x, t):
  return (t+1)*h*cos(2*pi*x)
