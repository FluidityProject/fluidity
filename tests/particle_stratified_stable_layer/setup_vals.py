from math import log, sqrt, pi

# List of constant values for simulations:
kappa = 1.
Ra    = 3.0e05

lda = 2.
u_0 = (lda**(7./3.)/ ( (1. + lda**4)**(2./3.))) * ((Ra/(2.*sqrt(pi)))**(2./3.))
Q   = 2. * sqrt(lda / (pi * u_0))

