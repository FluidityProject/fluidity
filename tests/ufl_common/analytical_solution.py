from numpy import cos, exp, matrix, pi, sqrt

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
    coeff = -(lmbda+2*(n**2)*pi**2)
    sins = cos(X[0]*pi*n)*cos(X[1]*pi*n)
    return coeff*sins
