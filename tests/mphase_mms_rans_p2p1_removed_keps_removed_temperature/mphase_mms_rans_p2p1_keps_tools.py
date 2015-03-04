from math import sin, cos, tanh, pi

def u1(X):
    return cos(X[0]) + 0.600*sin(X[1]) + 2.50

def v1(X):
    return X[1]*sin(X[0])

def u2(X):
    return cos(X[0]) + 0.600*sin(X[1]) + 2.50

def v2(X):
    return X[1]*sin(X[0])

def p(X):
    return cos(X[1]) + sin(X[0]) + sin(X[0]*X[1]/pi) - 1.00

def rho1(X):
    return 10.0

def rho2(X):
    return 10.0

def forcing_u1(X):
    return 4.20*X[1]*cos(X[1])*sin(X[0]) - (7.00*cos(X[0]) + 4.20*sin(X[1]) + 17.5)*sin(X[0]) + 0.700*X[1]*cos(X[0]*X[1]/pi)/pi + 1.19*cos(X[0]) + 0.294*sin(X[1]) - 4.94974746700000

def forcing_v1(X):
    return X[1]*(7.00*cos(X[0]) + 4.20*sin(X[1]) + 17.5)*cos(X[0]) + 7.00*X[1]*sin(X[0])**2 + 0.490*X[1]*sin(X[0]) + 0.700*X[0]*cos(X[0]*X[1]/pi)/pi - 0.700*sin(X[1]) - 4.94974746700000

def forcing_u2(X):
    return 1.80*X[1]*cos(X[1])*sin(X[0]) - (3.00*cos(X[0]) + 1.80*sin(X[1]) + 7.50)*sin(X[0]) + 0.300*X[1]*cos(X[0]*X[1]/pi)/pi + 0.510*cos(X[0]) + 0.126*sin(X[1]) - 2.12132034300000

def forcing_v2(X):
    return X[1]*(3.00*cos(X[0]) + 1.80*sin(X[1]) + 7.50)*cos(X[0]) + 3.00*X[1]*sin(X[0])**2 + 0.210*X[1]*sin(X[0]) + 0.300*X[0]*cos(X[0]*X[1]/pi)/pi - 0.300*sin(X[1]) - 2.12132034300000

def velocity1(X):
   return [u1(X), v1(X)]

def forcing_velocity1(X):
   return [forcing_u1(X), forcing_v1(X)]

def velocity2(X):
   return [u2(X), v2(X)]

def forcing_velocity2(X):
   return [forcing_u2(X), forcing_v2(X)]
