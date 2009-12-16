def val(X,t):
    from math import sin, cos, pi
    eps0 = 0.001
    return cos(pi*X[0]*eps0)*cos(pi*X[1]*eps0*2.0)*pi*pi*eps0*eps0*5.0 + \
                cos(2*pi*X[0]*eps0)*cos(3*pi*X[1]*eps0)*cos(pi*X[2]) * \
                pi*pi*(13*eps0*eps0 + 1)
