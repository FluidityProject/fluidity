#!/usr/bin/env python3

import scipy.integrate
import sys

(x, y) = [float(z) for z in sys.argv[1:3]]
t_start = 0.0
t_end = float(sys.argv[3])

def val(X, t):
  from math import exp, sin, cos
  ##L=1.e6; F=0.1; rho=1000.0; H=200.0; gamma=1.e-6;
  ##A = F*L/(pi*gamma*rho*H) 
  A = 159154.94309189534
  ##alpha = beta/gamma
  ##1/alpha = 100km=1e5; alpha = 1.e-5
  ##beta=1e-11
  ##zplus = -alpha/2 + sqrt((alpha**2)/4 + (pi/L)**2)
  zplus = 9.0504906000698326e-07
  ##zminus = -alpha/2 - sqrt((alpha**2)/4 + (pi/L)**2)
  zminus = -1.0905049060006982e-05  
  ##p=(1-exp(L*zminus))/(exp(L*zplus)-exp(L*zminus))  
  p =  0.40451761482981791
  ## q = 1-p
  q = 0.59548238517018204
  ## pi/L = 3.1415926535897933e-06
  ##psi = A*sin(3.1415926535897933e-06*y)*(p*exp(x*zplus) + q*exp(x*zminus) - 1.0)
  ## u = d psi/d y
  ## A*pi/L = 0.5
  u = 0.5*cos(3.1415926535897933e-06*X[1])*(p*exp(X[0]*zplus) + q*exp(X[0]*zminus) - 1.0)
  ## v = -d psi/d x
  v = -A*sin(3.1415926535897933e-06*X[1])*(p*zplus*exp(X[0]*zplus) + q*zminus*exp(X[0]*zminus))
  return [u,v]

def func(X, t):
  [u, v] = val(X, t)
  return [-u, -v]

print("func([x, y], t_start): ", func([x, y], t_start))

pos = scipy.integrate.odeint(func, [x, y], [t_start, t_end])
print("pos[0, :]: ", pos[0, :])
print("pos[-1, :]: ", pos[-1, :])
