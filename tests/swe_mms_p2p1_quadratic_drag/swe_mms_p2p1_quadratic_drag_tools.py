from math import sin, cos, tanh, pi, sqrt

def u(x):
    return cos(x[1])*sin(x[0])

def v(x):
    return -cos(x[0])*sin(x[1])

def h(x):
    return sin(x[0])*sin(x[1])

def forcing_u(x):
    return cos(x[0])*cos(x[1])**2*sin(x[0]) + cos(x[0])*sin(x[0])*sin(x[1])**2 + 1.20*cos(x[1])*sin(x[0]) + 9.80*cos(x[0])*sin(x[1]) + 0.00250*sqrt(cos(x[1])**2*sin(x[0])**2 + cos(x[0])**2*sin(x[1])**2)*cos(x[1])*sin(x[0])/(sin(x[0])*sin(x[1]) + 20.0)

def forcing_v(x):
    return cos(x[0])**2*cos(x[1])*sin(x[1]) + cos(x[1])*sin(x[0])**2*sin(x[1]) + 9.80*cos(x[1])*sin(x[0]) - 1.20*cos(x[0])*sin(x[1]) - 0.00250*sqrt(cos(x[1])**2*sin(x[0])**2 + cos(x[0])**2*sin(x[1])**2)*cos(x[0])*sin(x[1])/(sin(x[0])*sin(x[1]) + 20.0)

def velocity(x):
   return [u(x), v(x)]

def forcing_velocity(x):
   return [forcing_u(x), forcing_v(x)]

