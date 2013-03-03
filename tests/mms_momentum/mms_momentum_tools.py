from math import sin, cos, tanh, pi

def u(X):
    return sin(X[0]**2 + X[1]**2) + 0.500

def v(X):
    return 0.100*cos(X[0]**2 + X[1]**2) + 0.0500

def p(X):
    return 0.200*cos(X[0] + X[1]) + 0.300

def rho(X):
    return 1.00

def ie(X):
    return 0.500*cos(X[0] + X[1]) + 0.750

def forcing_u(X):
    return 2*(0.100*cos(X[0]**2 + X[1]**2) + 0.0500)*X[1]*cos(X[0]**2 + X[1]**2) + 4*(sin(X[0]**2 + X[1]**2) + 0.500)*X[0]*cos(X[0]**2 + X[1]**2) - 0.200*(sin(X[0]**2 + X[1]**2) + 0.500)*X[1]*sin(X[0]**2 + X[1]**2) + 2.80*X[0]**2*sin(X[0]**2 + X[1]**2) + 2.80*X[1]**2*sin(X[0]**2 + X[1]**2) - 0.200*sin(X[0] + X[1]) - 2.80*cos(X[0]**2 + X[1]**2)

def forcing_v(X):
    return 2*(0.100*cos(X[0]**2 + X[1]**2) + 0.0500)*X[0]*cos(X[0]**2 + X[1]**2) - 0.400*(0.100*cos(X[0]**2 + X[1]**2) + 0.0500)*X[1]*sin(X[0]**2 + X[1]**2) - 0.200*(sin(X[0]**2 + X[1]**2) + 0.500)*X[0]*sin(X[0]**2 + X[1]**2) + 0.280*X[0]**2*cos(X[0]**2 + X[1]**2) + 0.280*X[1]**2*cos(X[0]**2 + X[1]**2) - 0.200*sin(X[0] + X[1]) + 0.280*sin(X[0]**2 + X[1]**2)

def forcing_rho(X):
    return 2*X[0]*cos(X[0]**2 + X[1]**2) - 0.200*X[1]*sin(X[0]**2 + X[1]**2)

def forcing_ie(X):
    return 2*(0.200*cos(X[0] + X[1]) + 0.300)*X[0]*cos(X[0]**2 + X[1]**2) - 0.200*(0.200*cos(X[0] + X[1]) + 0.300)*X[1]*sin(X[0]**2 + X[1]**2) - 0.500*(0.100*cos(X[0]**2 + X[1]**2) + 0.0500)*sin(X[0] + X[1]) - 0.500*(sin(X[0]**2 + X[1]**2) + 0.500)*sin(X[0] + X[1])

def velocity(X):
   return [u(X), v(X)]

def forcing_velocity(X):
   return [forcing_u(X), forcing_v(X)]

