from __future__ import division
import scipy.special
import numpy
from math import sqrt, atan2, cos, sin, tan, acos, pi

nu = 1.
g = 1.
# outer and inner radius R_+ and R_-
Rp, Rm = 2.22, 1.22
alpha = Rm/Rp
l = 2
k = 3
m = l

# in the below:
# phi is the longitude (0<=phi<=2*pi)
# theta is the colatitude (0<=theta<=pi)
# NOTE: scipy uses the "mathematical" definition of phi and theta
# which is the reverse of the "physical" convention we use

# coefficients for Pl(r) = A*r**l + B*r**(-l-1) + C*r**(l+2) + D*r**(-l+1) + E*r**(k+3)
Gamma = ((alpha**(l + 1) + alpha**(l - 3))*(2*l + 1)**2 - 2*alpha**(l - 1)*(2*l + 3)*(2*l - 1) - 4*alpha**(3*l) - 4*alpha**(-l - 2))*(k + l + 4)*(k + l + 2)*(k - l + 3)*(k - l + 1)
A = ((alpha**(k + 2) + alpha**(l - 1))*(k + l + 2)*(2*l + 3) - (alpha**k + alpha**(l + 1))*(k + l + 4)*(2*l + 1) - 2*(alpha**(k + 2*l + 3) + alpha**(-l - 2))*(k - l + 1))*Rp**(-l + 3)*g/(Gamma*nu)
B = ((alpha**(k + 2*l + 1) + alpha**(l + 1))*(k - l + 3)*(2*l + 1) - (alpha**(k + 2*l + 3) + alpha**(l - 1))*(k - l + 1)*(2*l - 1) - 2*(alpha**(k + 2) + alpha**(3*l))*(k + l + 2))*Rp**(l + 4)*g/(Gamma*nu)
C = -((alpha**(k + 2) + alpha**(l - 3))*(k + l + 2)*(2*l + 1) - (alpha**k + alpha**(l - 1))*(k + l + 4)*(2*l - 1) - 2*(alpha**(k + 2*l + 1) + alpha**(-l - 2))*(k - l + 3))*Rp**(-l + 1)*g/(Gamma*nu)
D = -((alpha**(k + 2*l + 1) + alpha**(l - 1))*(k - l + 3)*(2*l + 3) - (alpha**(k + 2*l + 3) + alpha**(l - 3))*(k - l + 1)*(2*l + 1) - 2*(alpha**k + alpha**(3*l))*(k + l + 4))*Rp**(l + 2)*g/(Gamma*nu)
E = g/(Rp**k*(k + l + 4)*(k + l + 2)*(k - l + 3)*(k - l + 1)*nu)

# coefficients for pressure solution: p = G r^l + H r^{-l-1} + K r^{k+1}
G = -2*nu*(l+1)*(2*l+3)*C
H = -2*nu*l*(2*l-1)*D
K = -g*(k+2)/((k+1)*(k+2)-l*(l+1))/Rp**k


def Pl(r):
    return A*r**l + B*r**(-l-1) + C*r**(l+2) + D*r**(-l+1) + E*r**(k+3)


def dPldr(r):
    return l*A*r**(l-1) + (-l-1)*B*r**(-l-2) + (l+2)*C*r**(l+1) + (-l+1)*D*r**-l + (k+3)*E*r**(k+2)


# some sanity checks:
# no-normal flow
assert abs(Pl(Rm)) < 1e-12
assert abs(Pl(Rp)) < 1e-12


def Y(m, l, theta, phi):
    # everywhere we take the real part of Y, corresponding to the cos(m phi) part of the solution
    return scipy.special.sph_harm(m, l, phi, theta).real


def dYdphi(m, l, theta, phi):
    # except in theta derivatives (of odd order)
    return -m * scipy.special.sph_harm(m, l, phi, theta).imag


def dYdtheta(m, l, theta, phi):
    # this is from http://functions.wolfram.com/Polynomials/SphericalHarmonicY/20/01/01/
    # which can be derived from (x^2-1)d/dx P^ml = sqrt(1-x^2) P^(m+1)l -mx P^ml
    dydt = m/tan(theta) * Y(m, l, theta, phi)
    if m < l:
        # for m==l, Y^(m+1)_l=0
        # note we fiddle with phi to obtain the desired exp(im phi)
        # despite raising m to m+1
        dydt += sqrt((l-m)*(l+m+1)) * Y(m+1, l, theta, phi*m/(m+1))
    return dydt


def P(r, theta, phi):
    return Pl(r) * Y(m, l, theta, phi)


def pressure(r, theta, phi):
    return (G*r**l + H*r**(-l-1) + K*r**(k+1))*Y(m, l, theta, phi)


def u_theta(r, theta, phi):
    # u_theta = -1/r d/dtheta d/dr (r P)
    #         = -1/r dP/dtheta - d/dtheta d/dr P
    #         = -(1/r Pl + dPl/dr) dY/dtheta
    return -(Pl(r)/r + dPldr(r)) * dYdtheta(m, l, theta, phi)


def u_phi(r, theta, phi):
    # u_phi = -1/(r sin(theta)) d/dphi d/dr (r P)
    #       = -1/(r sin(theta)) * (dP/dphi + r d/dphi d/dr P)
    #       = -1/sin(theta) * (Pl/r + dPl/dr) * dY/dphi
    return -(Pl(r)/r + dPldr(r)) / sin(theta) * dYdphi(m, l, theta, phi)


def u_r(r, theta, phi):
    # u_r = 1/r \Lambda^2 P = -1/r l(l+1) P
    return -l*(l+1)*Pl(r)*Y(m, l, theta, phi)/r


def tau_rr(r, theta, phi):
    # tau_rr = 2 nu d/dr (1/r \Lambda^2 P)
    #        = -2 nu l(l+1) [1/r dP/dr - 1/r^2 P]
    return -2*nu*l*(l+1)*(dPldr(r) - Pl(r)/r)*Y(m, l, theta, phi)/r


def delta_rho(r, theta, phi):
    return r**k * Y(m, l, theta, phi) / Rp**k


def velocity_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    if theta < 1e-7*pi or theta > pi*(1.-1e-7):
        # workaround pole problem by averaging in 4 points near the pole
        dx = 1e-6*r
        return tuple(numpy.mean([velocity_cartesian((x, y, X[2])) for x, y in [[dx, dx], [-dx, dx], [dx, -dx], [-dx, -dx]]], axis=0))
    phi = atan2(X[1], X[0])
    ur = u_r(r, theta, phi)
    uth = u_theta(r, theta, phi)
    uph = u_phi(r, theta, phi)
    costh = cos(theta)
    req = sqrt(X[0]**2+X[1]**2)
    return (X[0]/r*ur + X[0]/req*costh*uth - X[1]/req*uph,
            X[1]/r*ur + X[1]/req*costh*uth + X[0]/req*uph,
            X[2]/r*ur - sin(theta)*uth)


def pressure_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return pressure(r, theta, phi)


def delta_rho_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return delta_rho(r, theta, phi)


def normal_stress_cartesian(X, i):
    # i==0: "Upper Mantle" (above anomaly)
    # i==1: below anomaly
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return tau_rr(r, theta, phi)[i]
