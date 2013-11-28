from math import sin, cos, tanh, pi, sqrt

def nu_T(X):
    return 0.800*cos(0.3183098861837907*X[0]*X[1]) + 0.400*sin(X[0]) + 0.300*sin(X[1]) + 2.00

def Ri(X):
    return 1

def c(X):
    return -X[1]

def forcing_c(X):
    return -0.415838293991294*X[0]*sin(0.3183098861837907*X[0]*X[1]) + 0.489897948556636*cos(X[1])

def u(X):
    return X[1]

def v(X):
    return 0

def U(X):
   return [u(X), v(X)]

