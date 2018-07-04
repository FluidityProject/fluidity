import scipy.special
import numpy
from math import sqrt, atan2, cos, sin, tan, acos, pi

nu = 1.
g = 1.
# outer and inner radius (these are a1 and a2 in [1])
R1 = 1.22
R2 = 2.22
alpha = R1/R2
l = 2
k = 3
m = l

# in the below:
# phi is the longitude (0<=phi<=2*pi)
# theta is the colatitude (0<=theta<=pi)
# NOTE: scipy uses the "mathematical" definition of phi and theta
# which is the reverse of the "physical" convention we use

# coefficients for Pl(r) = A*r**l + B*r**(-l-1) + C*r**(l+2) + D*r**(-l+1) + E*r**(k+3)
A = 0.5*R2**(-l + 3)*(alpha**(k + 3) - alpha**(-l + 1))*g/((alpha**l - alpha**(-l + 1))*(k + l + 2)*(k - l + 3)*(2*l + 1)*nu)
B = -0.5*R2**(l + 4)*(alpha**(k + 4) - alpha**(l + 3))*g/((alpha**(-l) - alpha**(l + 3))*(k + l + 4)*(k - l + 1)*(2*l + 1)*nu)
C = 0.5*R2**(-l + 1)*(alpha**(k + 4) - alpha**(-l))*g/((alpha**(-l) - alpha**(l + 3))*(k + l + 4)*(k - l + 1)*(2*l + 1)*nu)
D = -0.5*R2**(l + 2)*(alpha**(k + 3) - alpha**l)*g/((alpha**l - alpha**(-l + 1))*(k + l + 2)*(k - l + 3)*(2*l + 1)*nu)
E = R2**(-k)*g/((k + l + 4)*(k + l + 2)*(k - l + 3)*(k - l + 1)*nu)

# coefficients for pressure solution: p = G r^l + H r^{-l-1} + K r^{k+1}
G = -2*nu*(l+1)*(2*l+3)*C
H = -2*nu*l*(2*l-1)*D
K = -g*(k+2)/((k+1)*(k+2)-l*(l+1))/R2**k

def Pl(r):
    return A*r**l + B*r**(-l-1) + C*r**(l+2) + D*r**(-l+1) + E*r**(k+3)

def dPldr(r):
    return l*A*r**(l-1) + (-l-1)*B*r**(-l-2) + (l+2)*C*r**(l+1) + (-l+1)*D*r**-l + (k+3)*E*r**(k+2)

# some sanity checks:
# eqn (106):
assert abs(Pl(R1))<1e-12
assert abs(Pl(R2))<1e-12


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
    return (G*r**l + H*r**(-l-1) + K*r**(k+1))*Y(m, l, theta, phi)

def u_theta(r, theta, phi):
    # starting from eqn (44) in [1]:
    # u_theta = -1/r d/dtheta d/dr (r P)
    #         = -1/r dP/dtheta - d/dtheta d/dr P
    #         = -(1/r Pl + dPl/dr) dY/dtheta
    return -(Pl(r)/r + dPldr(r)) * dYdtheta(m, l, theta, phi)

def u_phi(r, theta, phi):
    # starting from eqn (44) in [1]:
    # u_phi = -1/(r sin(theta)) d/dphi d/dr (r P)
    #       = -1/(r sin(theta)) * (dP/dphi + r d/dphi d/dr P)
    #       = -1/sin(theta) * (Pl/r + dPl/dr) * dY/dphi
    return -(Pl(r)/r + dPldr(r)) / sin(theta) * dYdphi(m, l, theta, phi)

def u_r(r, theta, phi):
    # starting from eqn(44) in [1]:
    # u_r = 1/r B^2 P = -1/r l(l+1) P
    # (see above eqn (102) in [1])
    return -l*(l+1)*Pl(r)*Y(m, l, theta, phi)/r

def delta_rho(r, theta, phi):
    return r**k * Y(m, l, theta, phi) / R2**k

def get_cartesian_solution(X):
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    if theta<1e-7*pi or theta>pi*(1.-1e-7):
        # workaround pole problem by averaging in 4 points near the pole
        dx = 1e-6*r
        return tuple(numpy.mean([get_cartesian_solution((x,y,X[2])) for x,y in [[dx,dx],[-dx,dx],[dx,-dx],[-dx,-dx]]],axis=0))
    phi = atan2(X[1], X[0])
    ur = u_r(r,theta, phi)
    uth = u_theta(r,theta, phi)
    uph = u_phi(r, theta, phi)
    costh = cos(theta)
    req = sqrt(X[0]**2+X[1]**2)
    return ( X[0]/r*ur + X[0]/req*costh*uth - X[1]/req*uph,
             X[1]/r*ur + X[1]/req*costh*uth + X[0]/req*uph,
             X[2]/r*ur - sin(theta)*uth )

def get_pressure_solution(X):
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return pressure(r, theta, phi)

def get_delta_rho(X):
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return delta_rho(r, theta, phi)
