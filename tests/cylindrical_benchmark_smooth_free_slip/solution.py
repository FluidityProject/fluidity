from __future__ import division
import numpy
from math import sqrt, atan2, cos, sin

nu = 1.
g = 1.
# outer and inner radius R_+ and R_-
Rm, Rp = 1.22, 2.22
alpha = Rm/Rp
n = 2
k = 2

if n<=1:
    raise NotImplemented()

# velocity solution: coefficients for n, -n, n+2, -n+2 and k+3 powers of r
A = -0.25*(alpha**2 - alpha**(k + n + 3))*Rp**(-n + 3)*g/((alpha + alpha**n)*(alpha**n - alpha)*(k + n + 1)*(k - n + 3)*nu)
B = 0.25*Rp**(n + 3)*(alpha**(k + n + 3) - alpha**(2*n + 2))*g/((alpha**(n + 1) + 1)*(alpha**(n + 1) - 1)*(k + n + 3)*(k - n + 1)*nu)
C = -0.25*Rp**(-n + 1)*(alpha**(k + n + 3) - 1)*g/((alpha**(n + 1) + 1)*(alpha**(n + 1) - 1)*(k + n + 3)*(k - n + 1)*nu)
D = -0.25*Rp**(n + 1)*(alpha**(k + n + 3) - alpha**(2*n))*g/((alpha + alpha**n)*(alpha**n - alpha)*(k + n + 1)*(k - n + 3)*nu)
E = g*n/(((k + 3)**2 - n**2)*((k + 1)**2 - n**2)*Rp**k*nu)

# pressure: coefficients for n, -n, and k+1
G = -4*nu*C*(n+1)
H = -4*nu*D*(n-1)
F = -g*(k + 1)*Rp**(-k)/((k+1)**2-n**2)


def u_r(r, phi):
    dpsi_dphi = n*cos(n*phi)*(A*r**n+B*r**(-n)+C*r**(n+2)+D*r**(-n+2)+E*r**(k+3))
    return -dpsi_dphi/r

# some sanity checks:
# no-normal flow
numpy.testing.assert_almost_equal(u_r(Rm, 0.), 0.)
numpy.testing.assert_almost_equal(u_r(Rp, 0.), 0.)

def u_phi(r, phi):
    dpsi_dr = sin(n*phi)*(A*n*r**(n-1) + B*-n*r**(-n-1) + C*(n+2)*r**(n+1) + D*(-n+2)*r**(-n+1)+E*(k+3)*r**(k+2))
    return dpsi_dr



def p(r, phi):
    return (G*r**n + H*r**(-n) + F*r**(k+1))*cos(n*phi)

def normal_stress(r, phi):
    dpsi_dphi = n*cos(n*phi)*(A*r**n + B*r**(-n) + C*r**(n+2) + D*r**(-n+2) + E*r**(k+3))
    dpsi_drdphi = n*cos(n*phi)*(A*n*r**(n-1) + B*-n*r**(-n-1) + C*(n+2)*r**(n+1) + D*(-n+2)*r**(-n+1) + E*(k+3)*r**(k+2))
    tau_rr = 2*nu*(dpsi_dphi/r**2 - dpsi_drdphi/r)
    return tau_rr - p(r,phi)

def normal_stress_cartesian(X):
  r = sqrt(X[0]**2+X[1]**2)
  phi = atan2(X[1], X[0])
  return normal_stress(r, phi)

def delta_rho(r, phi):
    return r**k * cos(n*phi) / Rp**k

def velocity_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2)
    phi = atan2(X[1], X[0])
    ur = u_r(r,phi)
    ut = u_phi(r,phi)
    return [ur*X[0]/r - ut*X[1]/r, ur*X[1]/r + ut*X[0]/r]

def pressure_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2)
    phi = atan2(X[1], X[0])
    return p(r, phi)

def delta_rho_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2)
    phi = atan2(X[1], X[0])
    return delta_rho(r, phi)
