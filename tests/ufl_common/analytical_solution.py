from numpy import cos, pi

def helmholtz(X,t):
    n = 8
    return -cos(X[0]*pi*n)*cos(X[1]*pi*n)

def helmholtz_initial(X,t):
    n = 8
    lmbda = 1
    coeff = -(lmbda+2*(n**2)*pi**2)
    sins = cos(X[0]*pi*n)*cos(X[1]*pi*n)
    return coeff*sins
