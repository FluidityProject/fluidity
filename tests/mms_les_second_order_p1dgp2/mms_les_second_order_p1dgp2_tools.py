from math import sin, cos, tanh, pi, sqrt

def nu_T(X):
    return 4*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2)

def u(X):
    return 0.600*sin(X[1]) + cos(X[0]) + 3.00

def v(X):
    return X[1]*sin(X[0])

def forcing_u(X):
    return -16*(-0.300*sin(X[1]) + 0.500*cos(X[0]))*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) - (4*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 1.00)*(-0.600*sin(X[1]) + cos(X[0])) + 2*(4*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 1.00)*cos(X[0]) - 16*((0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*X[1]*sin(X[0]) - 2*sin(X[0])*cos(X[0]))*sin(X[0])/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 0.600*sin(X[1]) + cos(X[0]) + 3.00

def forcing_v(X):
    return (4*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 1.00)*X[1]*sin(X[0]) - 32*(-0.300*sin(X[1]) + 0.500*cos(X[0]))*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*sin(X[0])/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + X[1]*sin(X[0]) + 8*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*((0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*X[1]*sin(X[0]) - 2*sin(X[0])*cos(X[0]))/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2)

def U(X):
   return [u(X), v(X)]

def forcing_U(X):
   return [forcing_u(X), forcing_v(X)]

