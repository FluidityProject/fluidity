def val(X,t):
    from math import sin, cos, pi, sqrt
    
    g = cos(pi*X[0])*cos(2*pi*X[1]) + \
        cos(2*pi*X[0])*cos(3*pi*X[1])*cos(pi*X[2])
    
    R = sqrt( (X[0]-0.5)**2 + (X[1]-0.5)**2 + (X[2]-0.5)**2)

    return g
