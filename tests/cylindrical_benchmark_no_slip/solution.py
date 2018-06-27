import numpy
from math import sqrt, atan2, cos, sin, pi

nu = 1.
g = 1.
Rp, Rm = 2.22, 1.22
rp = Rm + 0.5
n = 2

if n<=1:
    raise NotImplemented()

alpha_pm = numpy.array([Rp/rp, Rm/rp])
alpha_mp = numpy.array([Rm/rp, Rp/rp])
alpha, beta = alpha_pm
pm = numpy.array([1, -1])
mp = -pm

# velocity solution: coefficients for n, -n, n+2, and -n+2 power of r
A_pm = 0.125*(((alpha**2 - beta**2)*n + (n + 1)*pm + alpha**(-2*n) - beta**(-2*n))*(n - 1) + (alpha**2*beta**(-2*n) - alpha**(-2*n)*beta**2)*n - (n**2*(beta/alpha)**(2*mp) - (beta/alpha)**(2*n*pm))*pm)*g*rp**(-n + 2)/((n**2*(alpha/beta - beta/alpha)**2 - ((beta/alpha)**(-n) - (beta/alpha)**n)**2)*(n - 1)*nu)
B_pm = 0.125*(((alpha**2 - beta**2)*n + (n - 1)*pm - alpha**(2*n) + beta**(2*n))*(n + 1) - (alpha**2*beta**(2*n) - alpha**(2*n)*beta**2)*n - (n**2*(beta/alpha)**(2*mp) - (beta/alpha)**(2*mp*n))*pm)*g*rp**(n + 2)/((n**2*(alpha/beta - beta/alpha)**2 - ((beta/alpha)**(-n) - (beta/alpha)**n)**2)*(n + 1)*nu)
C_pm = -0.125*((n**2*(beta/alpha)**(2*pm) - (beta/alpha)**(2*n*pm))*mp - (mp*(n - 1) - n*(1/alpha**2 - 1/beta**2) + alpha**(-2*n) - beta**(-2*n))*(n + 1) + n*(alpha**(-2*n)/beta**2 - beta**(-2*n)/alpha**2))*g*rp**(-n)/((n**2*(alpha/beta - beta/alpha)**2 - ((beta/alpha)**(-n) - (beta/alpha)**n)**2)*(n + 1)*nu)
D_pm = -0.125*((n**2*(beta/alpha)**(2*pm) - (beta/alpha)**(2*mp*n))*mp - (mp*(n + 1) - n*(1/alpha**2 - 1/beta**2) - alpha**(2*n) + beta**(2*n))*(n - 1) - n*(alpha**(2*n)/beta**2 - beta**(2*n)/alpha**2))*g*rp**n/((n**2*(alpha/beta - beta/alpha)**2 - ((beta/alpha)**(-n) - (beta/alpha)**n)**2)*(n - 1)*nu)


# pressure solution: coefficients for n and -n
G_pm = -4*nu*C_pm*(n+1)
H_pm = -4*nu*D_pm*(n-1)


def u_r(r, phi):
    dpsi_dphi = n*cos(n*phi)*(A_pm*r**n+B_pm*r**(-n)+C_pm*r**(n+2)+D_pm*r**(-n+2))
    return -dpsi_dphi/r

numpy.testing.assert_almost_equal(u_r(Rm, 0.)[1], 0.)
numpy.testing.assert_almost_equal(u_r(Rp, 0.)[0], 0.)
numpy.testing.assert_almost_equal(u_r(rp, 0.)[0], u_r(rp, 0.)[1])

def u_phi(r, phi):
    dpsi_dr = sin(n*phi)*(A_pm*n*r**(n-1) + B_pm*-n*r**(-n-1) + C_pm*(n+2)*r**(n+1) + D_pm*(-n+2)*r**(-n+1))
    return dpsi_dr

numpy.testing.assert_almost_equal(u_phi(rp, pi/(2*n))[0], u_phi(rp, pi/(2*n))[1])

numpy.testing.assert_almost_equal(u_phi(Rp, 0.), 0.)
numpy.testing.assert_almost_equal(u_phi(Rm, 0.), 0.)

def p(r, phi):
    return (G_pm*r**n + H_pm*r**(-n))*cos(n*phi)

def get_cartesian_solution(X, i):
  # i==0: upper mantle, i==1: lower mantle
  r = sqrt(X[0]**2+X[1]**2)
  phi = atan2(X[1], X[0])
  ur = u_r(r, phi)[i]
  ut = u_phi(r, phi)[i]
  return [ur*X[0]/r - ut*X[1]/r, ur*X[1]/r + ut*X[0]/r]

def get_cartesian_pressure_solution(X, i):
  # i==0: upper mantle, i==1: lower mantle
  r = sqrt(X[0]**2+X[1]**2)
  phi = atan2(X[1], X[0])
  return p(r, phi)[i]
