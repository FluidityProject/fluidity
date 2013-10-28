from math import sin, cos, tanh, pi, e, sqrt

def u(X):
    return cos(X[0]) + cos(X[1])

def v(X):
    return X[1]*sin(X[0])

def nu(X):
    return -1.20*sin(0.300*X[0]*X[1]) + 1.20*sin(1.70*X[0]) + 1.40*cos(1.10*X[1]) + 4.00

def forcing_u(X):
    return -(cos(X[0]) - cos(X[1]))*(-1.20*sin(0.300*X[0]*X[1]) + 1.20*sin(1.70*X[0]) + 1.40*cos(1.10*X[1]) + 4.00) - (X[1]*cos(X[0]) - sin(X[1]))*(-0.360*X[0]*cos(0.300*X[0]*X[1]) - 1.54*sin(1.10*X[1])) + (-0.720*X[1]*cos(0.300*X[0]*X[1]) + 4.08*cos(1.70*X[0]))*sin(X[0]) + (-2.40*sin(0.300*X[0]*X[1]) + 2.40*sin(1.70*X[0]) + 2.80*cos(1.10*X[1]) + 8.00)*cos(X[0]) + cos(X[0]) + cos(X[1])

def forcing_v(X):
    return (-1.20*sin(0.300*X[0]*X[1]) + 1.20*sin(1.70*X[0]) + 1.40*cos(1.10*X[1]) + 4.00)*X[1]*sin(X[0]) - (X[1]*cos(X[0]) - sin(X[1]))*(-0.360*X[1]*cos(0.300*X[0]*X[1]) + 2.04*cos(1.70*X[0])) - (-0.720*X[0]*cos(0.300*X[0]*X[1]) - 3.08*sin(1.10*X[1]))*sin(X[0]) + X[1]*sin(X[0])

def U(X):
    return [u(X), v(X)]

def forcing_U(X):
    return [forcing_u(X), forcing_v(X)]
