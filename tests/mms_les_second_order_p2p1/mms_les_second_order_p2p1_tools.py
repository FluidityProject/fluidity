from math import sin, cos, tanh, pi, sqrt

def u(X):
    return 0.600*sin(X[1]) + cos(X[0]) + 3.00

def v(X):
    return X[1]*sin(X[0])

def p(X):
    return sin(X[0]*X[1]/pi) + sin(X[0]) + cos(X[1]) - 1.00

def rho(X):
    return 2.50

def forcing_u(X):
    return 1.50*X[1]*sin(X[0])*cos(X[1]) - 0.160*(-0.300*sin(X[1]) + 0.500*cos(X[0]))*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) - (0.0400*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 0.700)*(-0.600*sin(X[1]) + cos(X[0])) + (0.0800*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 1.40)*cos(X[0]) - (1.50*sin(X[1]) + 2.50*cos(X[0]) + 7.50)*sin(X[0]) - 0.160*((0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*X[1]*sin(X[0]) - 2*sin(X[0])*cos(X[0]))*sin(X[0])/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + X[1]*cos(X[0]*X[1]/pi)/pi + cos(X[0])

def forcing_v(X):
    return (0.0400*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 0.700)*X[1]*sin(X[0]) + (1.50*sin(X[1]) + 2.50*cos(X[0]) + 7.50)*X[1]*cos(X[0]) + 2.50*X[1]*sin(X[0])**2 - 0.320*(-0.300*sin(X[1]) + 0.500*cos(X[0]))*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*sin(X[0])/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 0.0800*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*((0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*X[1]*sin(X[0]) - 2*sin(X[0])*cos(X[0]))/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + X[0]*cos(X[0]*X[1]/pi)/pi - sin(X[1])

def velocity(X):
   return [u(X), v(X)]

def forcing_velocity(X):
   return [forcing_u(X), forcing_v(X)]

