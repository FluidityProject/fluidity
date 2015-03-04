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

def temperature1(X):
    return -1.30*cos(2.10*X[1]) - 1.80*sin(1.70*X[0]) + 3.70*sin(1.30*X[0]*X[1]/pi) + 5.20

def temperature2(X):
    return -1.30*cos(2.10*X[1]) - 1.80*sin(1.70*X[0]) + 3.70*sin(1.30*X[0]*X[1]/pi) + 5.20

def rho1(X):
    return -13.0*cos(2.10*X[1]) - 18.0*sin(1.70*X[0]) + 37.0*sin(1.30*X[0]*X[1]/pi) + 62.0

def rho2(X):
    return -13.0*cos(2.10*X[1]) - 18.0*sin(1.70*X[0]) + 37.0*sin(1.30*X[0]*X[1]/pi) + 62.0

def forcing_u1(X):
    return 0.600*X[1]*(-9.10*cos(2.10*X[1]) - 12.6*sin(1.70*X[0]) + 25.9*sin(1.30*X[0]*X[1]/pi) + 43.4)*cos(X[1])*sin(X[0]) - (cos(X[0]) + 0.600*sin(X[1]) + 2.50)*(-9.10*cos(2.10*X[1]) - 12.6*sin(1.70*X[0]) + 25.9*sin(1.30*X[0]*X[1]/pi) + 43.4)*sin(X[0]) + 0.700*X[1]*cos(X[0]*X[1]/pi)/pi + 1.19*cos(X[0]) + 6.43467170710000*cos(2.10*X[1]) + 8.90954544060000*sin(1.70*X[0]) - 18.3140656279000*sin(1.30*X[0]*X[1]/pi) + 0.294*sin(X[1]) - 30.6884342954000

def forcing_v1(X):
    return X[1]*(cos(X[0]) + 0.600*sin(X[1]) + 2.50)*(-9.10*cos(2.10*X[1]) - 12.6*sin(1.70*X[0]) + 25.9*sin(1.30*X[0]*X[1]/pi) + 43.4)*cos(X[0]) + X[1]*(-9.10*cos(2.10*X[1]) - 12.6*sin(1.70*X[0]) + 25.9*sin(1.30*X[0]*X[1]/pi) + 43.4)*sin(X[0])**2 + 0.490*X[1]*sin(X[0]) + 0.700*X[0]*cos(X[0]*X[1]/pi)/pi + 6.43467170710000*cos(2.10*X[1]) + 8.90954544060000*sin(1.70*X[0]) - 18.3140656279000*sin(1.30*X[0]*X[1]/pi) - 0.700*sin(X[1]) - 30.6884342954000

def forcing_u2(X):
    return 0.600*X[1]*(-3.90*cos(2.10*X[1]) - 5.40*sin(1.70*X[0]) + 11.1*sin(1.30*X[0]*X[1]/pi) + 18.6)*cos(X[1])*sin(X[0]) - (cos(X[0]) + 0.600*sin(X[1]) + 2.50)*(-3.90*cos(2.10*X[1]) - 5.40*sin(1.70*X[0]) + 11.1*sin(1.30*X[0]*X[1]/pi) + 18.6)*sin(X[0]) + 0.300*X[1]*cos(X[0]*X[1]/pi)/pi + 0.510*cos(X[0]) + 2.75771644590000*cos(2.10*X[1]) + 3.81837661740000*sin(1.70*X[0]) - 7.84888526910000*sin(1.30*X[0]*X[1]/pi) + 0.126*sin(X[1]) - 13.1521861266000

def forcing_v2(X):
    return X[1]*(cos(X[0]) + 0.600*sin(X[1]) + 2.50)*(-3.90*cos(2.10*X[1]) - 5.40*sin(1.70*X[0]) + 11.1*sin(1.30*X[0]*X[1]/pi) + 18.6)*cos(X[0]) + X[1]*(-3.90*cos(2.10*X[1]) - 5.40*sin(1.70*X[0]) + 11.1*sin(1.30*X[0]*X[1]/pi) + 18.6)*sin(X[0])**2 + 0.210*X[1]*sin(X[0]) + 0.300*X[0]*cos(X[0]*X[1]/pi)/pi + 2.75771644590000*cos(2.10*X[1]) + 3.81837661740000*sin(1.70*X[0]) - 7.84888526910000*sin(1.30*X[0]*X[1]/pi) - 0.300*sin(X[1]) - 13.1521861266000

def forcing_temperature1(X):
    return X[1]*(4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*sin(X[0]) + (4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0]))*(cos(X[0]) + 0.600*sin(X[1]) + 2.50) + 6.25300000000000*X[0]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 + 6.25300000000000*X[1]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 - 5.73300000000000*cos(2.10*X[1]) - 5.20200000000000*sin(1.70*X[0])

def forcing_temperature2(X):
    return X[1]*(4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*sin(X[0]) + (4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0]))*(cos(X[0]) + 0.600*sin(X[1]) + 2.50) + 6.25300000000000*X[0]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 + 6.25300000000000*X[1]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 - 5.73300000000000*cos(2.10*X[1]) - 5.20200000000000*sin(1.70*X[0])

def velocity1(X):
   return [u1(X), v1(X)]

def forcing_velocity1(X):
   return [forcing_u1(X), forcing_v1(X)]

def velocity2(X):
   return [u2(X), v2(X)]

def forcing_velocity2(X):
   return [forcing_u2(X), forcing_v2(X)]
