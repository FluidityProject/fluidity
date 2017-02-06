import numpy
from math import sqrt, atan2, cos, sin

Ra = 1.
R = numpy.array([2.22, 1.22])
Rr = R[::-1]
rp = R[1]+0.5
n = 2

if n>1:
    E = -Ra*rp**(-n)/(8*(n+1)) * \
            (Rr**(2*n+2)-rp**(2*n+2))/(R[0]**(2*n+2)-R[1]**(2*n+2))
    F = -Ra*R**(2*n)*rp**(-n)/(8*(n-1)) * \
            (Rr**2*rp**(2*n)-Rr**(2*n)*rp**2)/(R[0]**2*R[1]**(2*n)-R[0]**(2*n)*R[1]**2)
    C = -F/R**(2*n-2)
    D = -E*R**(2*n+2)
else:
    raise NotImplemented()


def u_r(r, theta):
  dpsi_dtheta = n*cos(n*theta)*(C*r**n+D*r**(-n)+E*r**(n+2)+F*r**(-n+2))
  return -dpsi_dtheta/r

def u_theta(r, theta):
  dpsi_dr = sin(n*theta)*(C*n*r**(n-1) + D*-n*r**(-n-1) + E*(n+2)*r**(n+1) + F*(-n+2)*r**(-n+1))
  return dpsi_dr

def get_cartesian_solution(X, i):
  r = sqrt(X[0]**2+X[1]**2)
  theta = atan2(X[1], X[0])
  ur = u_r(r,theta)[i]
  ut = u_theta(r,theta)[i]
  return [ur*X[0]/r - ut*X[1]/r, ur*X[1]/r + ut*X[0]/r]
