from math import sin, cos, tanh, pi

def u(X):
    return cos(X[0]) + 0.600*sin(X[1]) + 2.50

def v(X):
    return X[1]*sin(X[0])

def p(X):
    return cos(X[1]) + sin(X[0]) + sin(X[0]*X[1]/pi) - 1.00

def temperature(X):
    return -1.30*cos(2.10*X[1]) - 1.80*sin(1.70*X[0]) + 3.70*sin(1.30*X[0]*X[1]/pi) + 5.20

def rho(X):
    return -13.0*cos(2.10*X[1]) - 18.0*sin(1.70*X[0]) + 37.0*sin(1.30*X[0]*X[1]/pi) + 62.0

def forcing_u(X):
    return 0.600*X[1]*(-13.0*cos(2.10*X[1]) - 18.0*sin(1.70*X[0]) + 37.0*sin(1.30*X[0]*X[1]/pi) + 62.0)*cos(X[1])*sin(X[0]) - (cos(X[0]) + 0.600*sin(X[1]) + 2.50)*(-13.0*cos(2.10*X[1]) - 18.0*sin(1.70*X[0]) + 37.0*sin(1.30*X[0]*X[1]/pi) + 62.0)*sin(X[0]) + X[1]*cos(X[0]*X[1]/pi)/pi + 1.70*cos(X[0]) + 9.19238815300000*cos(2.10*X[1]) + 12.7279220580000*sin(1.70*X[0]) - 26.1629508970000*sin(1.30*X[0]*X[1]/pi) + 0.420*sin(X[1]) - 43.8406204220000

def forcing_v(X):
    return X[1]*(cos(X[0]) + 0.600*sin(X[1]) + 2.50)*(-13.0*cos(2.10*X[1]) - 18.0*sin(1.70*X[0]) + 37.0*sin(1.30*X[0]*X[1]/pi) + 62.0)*cos(X[0]) + X[1]*(-13.0*cos(2.10*X[1]) - 18.0*sin(1.70*X[0]) + 37.0*sin(1.30*X[0]*X[1]/pi) + 62.0)*sin(X[0])**2 + 0.700*X[1]*sin(X[0]) + X[0]*cos(X[0]*X[1]/pi)/pi + 9.19238815300000*cos(2.10*X[1]) + 12.7279220580000*sin(1.70*X[0]) - 26.1629508970000*sin(1.30*X[0]*X[1]/pi) - sin(X[1]) - 43.8406204220000

def forcing_temperature(X):
    return X[1]*(4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*sin(X[0]) + (4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0]))*(cos(X[0]) + 0.600*sin(X[1]) + 2.50) + 6.25300000000000*X[0]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 + 6.25300000000000*X[1]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 - 5.73300000000000*cos(2.10*X[1]) - 5.20200000000000*sin(1.70*X[0])

def velocity(X):
   return [u(X), v(X)]

def forcing_velocity(X):
   return [forcing_u(X), forcing_v(X)]
