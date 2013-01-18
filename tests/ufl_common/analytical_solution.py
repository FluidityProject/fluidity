from numpy import sin, cos, exp, matrix, pi, sqrt

def advection_diffusion(x, t):
    dx= (matrix(x)-matrix((0.45 + t,0.5)))
    r=sqrt(dx*dx.T)
    D = 0.1 # Diffusivity
    A = 0.1 # Normalisation
    return A*(exp((-r**2)/(4*D*t))/(4*pi*D*t))

def helmholtz(X,t):
    n = 8
    return -cos(X[0]*pi*n)*cos(X[1]*pi*n)

def helmholtz_initial(X,t):
    n = 8
    lmbda = 1
    return -(lmbda+2*(n**2)*pi**2) * cos(X[0]*pi*n)*cos(X[1]*pi*n)

def helmholtz_sin(X,t):
    n = 8
    return -sin(X[0]*pi*n)*sin(X[1]*pi*n)

def helmholtz_initial_sin(X,t):
    n = 8
    lmbda = 1
    return -(lmbda+2*(n**2)*pi**2) * sin(X[0]*pi*n)*sin(X[1]*pi*n)
