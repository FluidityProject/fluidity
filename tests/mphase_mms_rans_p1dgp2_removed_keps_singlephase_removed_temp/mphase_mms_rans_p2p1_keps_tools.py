from math import sin, cos, tanh, pi

def u(X):
    return cos(X[0]) + 0.600*sin(X[1]) + 2.50

def v(X):
    return X[1]*sin(X[0])

def p(X):
    return cos(X[1]) + sin(X[0]) + sin(X[0]*X[1]/pi) - 1.00

def rho(X):
    return 10.0

def forcing_u(X):
    return 6.00*X[1]*cos(X[1])*sin(X[0]) - (10.0*cos(X[0]) + 6.00*sin(X[1]) + 25.0)*sin(X[0]) + X[1]*cos(X[0]*X[1]/pi)/pi + 1.70*cos(X[0]) + 0.420*sin(X[1]) - 7.07106781000000

def forcing_v(X):
    return X[1]*(10.0*cos(X[0]) + 6.00*sin(X[1]) + 25.0)*cos(X[0]) + 10.0*X[1]*sin(X[0])**2 + 0.700*X[1]*sin(X[0]) + X[0]*cos(X[0]*X[1]/pi)/pi - sin(X[1]) - 7.07106781000000

def velocity(X):
   return [u(X), v(X)]

def forcing_velocity(X):
   return [forcing_u(X), forcing_v(X)]
