
from fluidity_tools import stat_parser as stat
from math import log, pi, cos, sin, cosh, sinh, exp
from numpy import sqrt, interp
from scipy.integrate import quad

wavelengthfactor = 1.0
dx = 1./20. # nond

eta0 = 1000. # m
xi0 = 1000. # m
rho0 = 4500. # kg/m**3
rhou = 2*rho0 # kg/m**3
mu0 = 1.e21  # Pas
g = 10.      # m/s**2
D = 3.e6     # m

#eta0 = 1./3000. # m
#xi0 = 1./3000. # m
#rho0 = 1.0 # kg/m**3
#rhou = 2.0 # kg/m**3
#mu0 = 1.0  # Pas
#g = 1.0      # m/s**2
#D = 1.0     # m

def wavenumber():
  wavelength = wavelengthfactor*D
  return 2.0*pi/wavelength

def nond_wavenumber():
  wavelength = wavelengthfactor
  return 2.0*pi/wavelength

def t0_eta():
  delta_rho = (rhou-rho0)*g
  rhog = rho0*g
  mu = mu0
  k = wavenumber()
  return -1.0*(((delta_rho*k**2*mu + k**2*mu*rhog)*D*sinh(D*k)**2 - (delta_rho*k**2*mu + k**2*mu*rhog)*D*cosh(D*k)**2 \
          - (delta_rho*k*mu + k*mu*rhog)*sinh(D*k)*cosh(D*k) - sqrt((delta_rho**2*k**2 + 2*delta_rho*k**2*rhog \
          + k**2*rhog**2)*D**2*cosh(D*k)**4 - 2*(delta_rho**2*k + 2*delta_rho*k*rhog + k*rhog**2)*D*sinh(D*k)**3*cosh(D*k) \
          + 2*(delta_rho**2*k + 2*delta_rho*k*rhog + k*rhog**2)*D*sinh(D*k)*cosh(D*k)**3 \
          + ((delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2)*D**2 - 4*delta_rho*rhog)*sinh(D*k)**4 \
          - (2*(delta_rho**2*k**2 + k**2*rhog**2)*D**2 - delta_rho**2 - 2*delta_rho*rhog \
          - rhog**2)*sinh(D*k)**2*cosh(D*k)**2)*k*mu)/(delta_rho*rhog*sinh(D*k)**2))

def t0_xi():
  delta_rho = (rhou-rho0)*g
  rhog = rho0*g
  mu = mu0
  k = wavenumber()
  return -1.0*(((delta_rho*k**2*mu + k**2*mu*rhog)*D*sinh(D*k)**2 - (delta_rho*k**2*mu + k**2*mu*rhog)*D*cosh(D*k)**2 \
               - (delta_rho*k*mu + k*mu*rhog)*sinh(D*k)*cosh(D*k) + sqrt((delta_rho**2*k**2 + 2*delta_rho*k**2*rhog \
               + k**2*rhog**2)*D**2*cosh(D*k)**4 - 2*(delta_rho**2*k + 2*delta_rho*k*rhog \
               + k*rhog**2)*D*sinh(D*k)**3*cosh(D*k) + 2*(delta_rho**2*k + 2*delta_rho*k*rhog \
               + k*rhog**2)*D*sinh(D*k)*cosh(D*k)**3 + ((delta_rho**2*k**2 - 2*delta_rho*k**2*rhog \
               + k**2*rhog**2)*D**2 - 4*delta_rho*rhog)*sinh(D*k)**4 - (2*(delta_rho**2*k**2 + k**2*rhog**2)*D**2 \
               - delta_rho**2 - 2*delta_rho*rhog - rhog**2)*sinh(D*k)**2*cosh(D*k)**2)*k*mu)/(delta_rho*rhog*sinh(D*k)**2))

def t0():
  return min(t0_eta(), t0_xi())

def nond_factor():
  return rho0*g*D*t0()/mu0

def nond_eta0():
  return eta0/D

def nond_xi0():
  return xi0/D

def nond_F(x,t):
  k = nond_wavenumber()
  F0 = nond_eta0()
  G0 = nond_xi0()
  delta_rho = (rhou-rho0)*nond_factor()/rho0
  rhog = nond_factor() # use this as a proxy for the nondimensional factorisation
  return ((0.5*(exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k + k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog - sqrt((delta_rho**2\
- 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k + 2*delta_rho*k*rhog\
  + k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog + delta_rho**2 + 2*delta_rho*rhog +\
rhog**2)*sinh(k)**2)*k))/((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho -\
rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 +\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)))\
- exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k + k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog + sqrt((delta_rho**2 -\
  2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k + 2*delta_rho*k*rhog\
+ k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog + delta_rho**2 + 2*delta_rho*rhog +\
rhog**2)*sinh(k)**2)*k))/((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho -\
rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 +\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))))*G0 + (((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho -\
rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 +\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k -\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)))\
+ 1)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k + k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog + sqrt((delta_rho**2 -\
2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k + 2*delta_rho*k*rhog +\
k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog + delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**2)*k)) -\
(k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k\
+ k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog - sqrt((delta_rho**2 - 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 +\
2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k + 2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog +\
delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 -\
sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4\
+ 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 -\
2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) +\
4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2\
- 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) -\
  k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2\
+ (delta_rho*k - k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 +\
delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 -\
2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 -\
4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) +\
2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k -\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))))*F0))*cos(k*x)

def nond_G(x,t):
  k = nond_wavenumber()
  F0 = nond_eta0()
  G0 = nond_xi0()
  delta_rho = (rhou-rho0)*nond_factor()/rho0
  rhog = nond_factor() # use this as a proxy for the nondimensional factorisation
  return ((((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k + k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog -\
sqrt((delta_rho**2 - 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k +\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog + delta_rho**2 + 2*delta_rho*rhog +\
rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho -\
rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 +\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k -\
k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)))\
- (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
  rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k + k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog +\
sqrt((delta_rho**2 - 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k +\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog + delta_rho**2 + 2*delta_rho*rhog +\
rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho -\
rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 +\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k -\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))))*G0 + 2*(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho -\
rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 +\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k -\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)))\
+ 1)*(k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k + k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog +\
sqrt((delta_rho**2 - 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k +\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog + delta_rho**2 + 2*delta_rho*rhog +\
rhog**2)*sinh(k)**2)*k))/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 +\
sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4\
+ 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 -\
2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) +\
4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2\
- 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) -\
  k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)**2*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k +\
k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 + k**2*rhog - sqrt((delta_rho**2 - 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 +\
2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k + 2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) + (4*delta_rho*k**2*rhog +\
delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k - k*rhog)*cosh(k)**2 -\
sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4\
+ 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 -\
2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) +\
4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2\
- 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) -\
  k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2\
+ (delta_rho*k - k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 +\
delta_rho**2*k**2*cosh(k)**4 - 2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 -\
2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 -\
4*delta_rho*k*rhog*sinh(k)**3*cosh(k) + 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) +\
2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k -\
k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))*((delta_rho - rhog)*sinh(k)*cosh(k) - (delta_rho*k - k*rhog)*sinh(k)**2 + (delta_rho*k -\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 -\
2*delta_rho*k**2*rhog*sinh(k)**4 + 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 - 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
+ 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
delta_rho**2*sinh(k)**2*cosh(k)**2 - 4*delta_rho*rhog*sinh(k)**4 + 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))))*F0))*cos(k*x)


def nond_amp(t):
  return max(nond_F(0.0,t), nond_G(0.0,t))

def numerical_amp(statn, t):
  return interp([t], statn["ElapsedTime"]["value"], statn["Fluid"]["FreeSurface"]["max"])[0]

def nond_error_amp(statn, t):
  return abs(nond_amp(t)-numerical_amp(statn, t))


