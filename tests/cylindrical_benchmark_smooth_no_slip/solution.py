import numpy
from math import sqrt, atan2, cos, sin

nu = 1.
g = 1.
R1, R2 = 1.22, 2.22
alpha = R1/R2
n = 2
k = 2

if n<=1:
    raise NotImplemented()

# velocity solution: coefficients for n, -n, n+2, -n+2 and k+3 powers of r
A = 0.5*((alpha**(k + n + 3) + alpha**(2*n))*(k + n + 1)*(n + 1) - (alpha**(k + n + 1) + alpha**(2*n + 2))*(k + n + 3)*n - (alpha**(k + 3*n + 3) + 1)*(k - n + 1))*R2**(-n + 3)*g*n/(((alpha**(n + 1) - alpha**(n - 1))**2*n**2 - (alpha**(2*n) - 1)**2)*((k + 3)**2 - n**2)*((k + 1)**2 - n**2)*nu)
B = -0.5*((alpha**(k + 3*n + 3) + alpha**(2*n))*(k - n + 1)*(n - 1) - (alpha**(k + 3*n + 1) + alpha**(2*n + 2))*(k - n + 3)*n + (alpha**(k + n + 3) + alpha**(4*n))*(k + n + 1))*R2**(n + 3)*g*n/(((alpha**(n + 1) - alpha**(n - 1))**2*n**2 - (alpha**(2*n) - 1)**2)*((k + 3)**2 - n**2)*((k + 1)**2 - n**2)*nu)
C = 0.5*((alpha**(k + n + 1) + alpha**(2*n))*(k + n + 3)*(n - 1) - (alpha**(k + n + 3) + alpha**(2*n - 2))*(k + n + 1)*n + (alpha**(k + 3*n + 1) + 1)*(k - n + 3))*R2**(-n + 1)*g*n/(((alpha**(n + 1) - alpha**(n - 1))**2*n**2 - (alpha**(2*n) - 1)**2)*((k + 3)**2 - n**2)*((k + 1)**2 - n**2)*nu)
D = -0.5*((alpha**(k + 3*n + 1) + alpha**(2*n))*(k - n + 3)*(n + 1) - (alpha**(k + 3*n + 3) + alpha**(2*n - 2))*(k - n + 1)*n - (alpha**(k + n + 1) + alpha**(4*n))*(k + n + 3))*R2**(n + 1)*g*n/(((alpha**(n + 1) - alpha**(n - 1))**2*n**2 - (alpha**(2*n) - 1)**2)*((k + 3)**2 - n**2)*((k + 1)**2 - n**2)*nu)
E = R2**(-k)*g*n/(((k + 3)**2 - n**2)*((k + 1)**2 - n**2)*nu)

# pressure: coefficients for n, -n, and k+1
G = -4*nu*C*(n+1)
H = -4*nu*D*(n-1)
F = -g*(k + 1)*R2**(-k)/((k+1)**2-n**2)


def u_r(r, phi):
    dpsi_dphi = n*cos(n*phi)*(A*r**n+B*r**(-n)+C*r**(n+2)+D*r**(-n+2)+E*r**(k+3))
    return -dpsi_dphi/r

numpy.testing.assert_almost_equal(u_r(R1, 0.), 0.)
numpy.testing.assert_almost_equal(u_r(R2, 0.), 0.)

def u_phi(r, phi):
    dpsi_dr = sin(n*phi)*(A*n*r**(n-1) + B*-n*r**(-n-1) + C*(n+2)*r**(n+1) + D*(-n+2)*r**(-n+1)+E*(k+3)*r**(k+2))
    return dpsi_dr

numpy.testing.assert_almost_equal(u_phi(R1, 0.), 0.)
numpy.testing.assert_almost_equal(u_phi(R2, 0.), 0.)

def p(r, phi):
    return (G*r**n + H*r**(-n) + F*r**(k+1))*cos(n*phi)

def delta_rho(r, phi):
    return r**k * cos(n*phi) / R2**k

def get_cartesian_solution(X):
    r = sqrt(X[0]**2+X[1]**2)
    phi = atan2(X[1], X[0])
    ur = u_r(r,phi)
    ut = u_phi(r,phi)
    return [ur*X[0]/r - ut*X[1]/r, ur*X[1]/r + ut*X[0]/r]

def get_cartesian_pressure_solution(X):
    r = sqrt(X[0]**2+X[1]**2)
    phi = atan2(X[1], X[0])
    return p(r, phi)

def delta_rho_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2)
    phi = atan2(X[1], X[0])
    return delta_rho(r, phi)
