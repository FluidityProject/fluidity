import scipy.special
import numpy
from math import sqrt, atan2, cos, sin, tan, acos, pi

nu = 1.
g = 1.
# outer and inner radius R_+ and R_-
Rp, Rm = 2.22, 1.22
# radial height of anomaly: r'
rp = Rm + 0.5
l = 2
m = l

# in the below:
# phi is the longitude (0<=phi<=2*pi)
# theta is the colatitude (0<=theta<=pi)
# NOTE: scipy uses the "mathematical" definition of phi and theta
# which is the reverse of the "physical" convention we use

alpha_pm = numpy.array([Rp/rp, Rm/rp])
alpha_mp = numpy.array([Rm/rp, Rp/rp])
alpha, beta = alpha_pm
R_pm = alpha_pm*rp
pm = numpy.array([1, -1])
mp = -pm

# coefficients for Pl(r) = A_pm*r**l + B_pm*r**(-l-1) + C_pm*r**(l+2) + D_pm*r**(-l+1)
# since the alpha's are length-2 arrays (each entry corresponding to one halve of the domain)
# so will the coefficients and thus the solution
A_pm = -0.5*(alpha_mp**(2*l - 1) - 1)*g*pm*rp**(-l + 2)/((alpha_mp**(2*l - 1) - alpha_pm**(2*l - 1))*(2*l + 1)*(2*l - 1)*nu)
B_pm = -0.5*(alpha_mp**(-2*l - 3) - 1)*g*pm*rp**(l + 3)/((alpha_mp**(-2*l - 3) - alpha_pm**(-2*l - 3))*(2*l + 3)*(2*l + 1)*nu)
C_pm = 0.5*(alpha_mp**(2*l + 3) - 1)*g*pm*rp**(-l)/((alpha_mp**(2*l + 3) - alpha_pm**(2*l + 3))*(2*l + 3)*(2*l + 1)*nu)
D_pm = 0.5*(alpha_mp**(-2*l + 1) - 1)*g*pm*rp**(l + 1)/((alpha_mp**(-2*l + 1) - alpha_pm**(-2*l + 1))*(2*l + 1)*(2*l - 1)*nu)


# coefficients for pressure solution: p = G_pm r^l + H_pm r^{-l-1}
G_pm = -2*nu*(l+1)*(2*l+3)*C_pm
H_pm = -2*nu*l*(2*l-1)*D_pm

def Pl(r):
    return A_pm*r**l + B_pm*r**(-l-1) + C_pm*r**(l+2) + D_pm*r**(-l+1)

def dPldr(r):
    return l*A_pm*r**(l-1) + (-l-1)*B_pm*r**(-l-2) + (l+2)*C_pm*r**(l+1) + (-l+1)*D_pm*r**-l

# some sanity checks:
# no-normal flow
assert all(numpy.abs(Pl(R_pm))<1e-12)
# continuity of velocity at r=rp:
assert abs(Pl(rp)[0]-Pl(rp)[1])<1e-12
assert abs(dPldr(rp)[0]-dPldr(rp)[1])<1e-12


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
    if m<l:
        # for m==l, Y^(m+1)_l=0
        # note we fiddle with phi to obtain the desired exp(im phi)
        # despite raising m to m+1
        dydt += sqrt((l-m)*(l+m+1)) * Y(m+1, l, theta, phi*m/(m+1))
    return dydt

def P(r, theta, phi):
    return Pl(r) * Y(m, l, theta, phi)

def pressure(r, theta,phi):
    return (G_pm*r**l + H_pm*r**(-l-1))*Y(m, l, theta, phi)

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

def velocity_cartesian(X, i):
    # i==0: "Upper Mantle" (above anomaly)
    # i==1: below anomaly
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    if theta<1e-7*pi or theta>pi*(1.-1e-7):
        # workaround pole problem by averaging in 4 points near the pole
        dx = 1e-6*r
        return tuple(numpy.mean([velocity_cartesian((x,y,X[2])) for x,y in [[dx,dx],[-dx,dx],[dx,-dx],[-dx,-dx]]],axis=0))
    phi = atan2(X[1], X[0])
    ur = u_r(r,theta, phi)[i]
    uth = u_theta(r,theta, phi)[i]
    uph = u_phi(r, theta, phi)[i]
    costh = cos(theta)
    req = sqrt(X[0]**2+X[1]**2)
    return ( X[0]/r*ur + X[0]/req*costh*uth - X[1]/req*uph,
             X[1]/r*ur + X[1]/req*costh*uth + X[0]/req*uph,
             X[2]/r*ur - sin(theta)*uth )

def pressure_cartesian(X, i):
    # i==0: "Upper Mantle" (above anomaly)
    # i==1: below anomaly
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return pressure(r, theta, phi)[i]

