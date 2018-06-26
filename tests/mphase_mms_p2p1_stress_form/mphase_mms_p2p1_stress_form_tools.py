from math import sin, cos, tanh, pi

def u1(X):
    return -X[0]*cos(X[1]) + 0.250*cos(X[0])*cos(X[1])

def v1(X):
    return sin(X[1])

def u2(X):
    return cos(X[1])*sin(X[0])

def v2(X):
    return -cos(X[0])*sin(X[1]) + sin(X[0])*sin(X[1])

def vfrac1(X):
    return 0.800

def vfrac2(X):
    return 0.200

def p(X):
    return 1.00*cos(1.00*X[1]) + 1.00*sin(1.00*X[0]) + 1.00*sin(1.00*X[0]*X[1]/pi) - 1.16128283994290

def rho1(X):
    return 2.50

def rho2(X):
    return 0.500

def forcing_u1(X):
    return (-2.00*X[0]*cos(X[1]) + 0.500*cos(X[0])*cos(X[1]))*(-0.250*cos(X[1])*sin(X[0]) - cos(X[1])) - 0.560*X[0]*cos(X[1]) + 0.326666666666667*cos(X[0])*cos(X[1]) + 2.00*(X[0]*sin(X[1]) - 0.250*cos(X[0])*sin(X[1]))*sin(X[1]) + 0.800*X[1]*cos(1.00*X[0]*X[1]/pi)/pi + 0.800*cos(1.00*X[0])

def forcing_v1(X):
    return 2.00*cos(X[1])*sin(X[1]) - 0.0466666666666667*sin(X[0])*sin(X[1]) + 0.800*X[0]*cos(1.00*X[0]*X[1]/pi)/pi - 0.800*sin(1.00*X[1]) + 0.560*sin(X[1])

def forcing_u2(X):
    return 0.100*cos(X[0])*cos(X[1])**2*sin(X[0]) - (-0.100*cos(X[0])*sin(X[1]) + 0.100*sin(X[0])*sin(X[1]))*sin(X[0])*sin(X[1]) - 0.0466666666666666*cos(X[0])*cos(X[1]) + 0.280*cos(X[1])*sin(X[0]) + 0.200*X[1]*cos(1.00*X[0]*X[1]/pi)/pi + 0.200*cos(1.00*X[0])

def forcing_v2(X):
    return 0.100*(cos(X[0])*sin(X[1]) + sin(X[0])*sin(X[1]))*cos(X[1])*sin(X[0]) - (cos(X[0])*cos(X[1]) - cos(X[1])*sin(X[0]))*(-0.100*cos(X[0])*sin(X[1]) + 0.100*sin(X[0])*sin(X[1])) - 0.280*cos(X[0])*sin(X[1]) + 0.326666666666667*sin(X[0])*sin(X[1]) + 0.200*X[0]*cos(1.00*X[0]*X[1]/pi)/pi - 0.200*sin(1.00*X[1])

def velocity1(X):
   return [u1(X), v1(X)]

def velocity2(X):
   return [u2(X), v2(X)]

def forcing_velocity1(X):
   return [forcing_u1(X), forcing_v1(X)]

def forcing_velocity2(X):
   return [forcing_u2(X), forcing_v2(X)]

