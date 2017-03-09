import scipy.special
import numpy
from math import sqrt, atan2, cos, sin, tan, acos, pi

Ra = 1.
# outer and inner radius (these are a1 and a2 in [1])
R = numpy.array([2.22, 1.22])
# r', the radius of the anomaly
rp = R[1]+0.5
l = 2
m = l

# in the below:
# phi is the longitude (0<=phi<=2*pi)
# theta is the colatitude (0<=theta<=pi)
# NOTE: scipy uses the "mathematical" definition of phi and theta
# which is the reverse of the "physical" convention we use

# solutions based on free slip solutions in [1] N.M. Ribe, Analytical Approaches to Mantle Dynamics, in Mantle Dynamics, Vol. 7
# no slip solution here is arrived at by solving (106-109) in sympy, but with (107) replaced with:
#   dPl/dr(R_1, r') = dPl/dr(R_2, r') = 0
# (i.e. zero *first* derivatives)

# this horrendous mess, comes out of sympy:
R1, R2 = R
A = numpy.array([
      R1*Ra*rp**(-l + 1)*(-4*R1**(2*l + 1)*R2*l*rp**(2*l + 3) + 2*R1**(2*l + 1)*R2*rp**(2*l + 3) + 4*R1**(2*l + 1)*R2**(2*l + 2)*l**2*rp**2 + 4*R1**(2*l + 1)*R2**(2*l + 2)*l*rp**2 - 3*R1**(2*l + 1)*R2**(2*l + 2)*rp**2 - 4*R1**(2*l + 1)*R2**(2*l + 4)*l**2 + R1**(2*l + 1)*R2**(2*l + 4) + 4*R1**(2*l + 3)*R2*l*rp**(2*l + 1) + 2*R1**(2*l + 3)*R2*rp**(2*l + 1) - 4*R1**(2*l + 3)*R2**(2*l)*l**2*rp**2 - 4*R1**(2*l + 3)*R2**(2*l)*l*rp**2 - R1**(2*l + 3)*R2**(2*l)*rp**2 + 4*R1**(2*l + 3)*R2**(2*l + 2)*l**2 - R1**(2*l + 3)*R2**(2*l + 2) + 4*R2**(2*l + 2)*l*rp**(2*l + 3) - 2*R2**(2*l + 2)*rp**(2*l + 3) - 4*R2**(2*l + 4)*l*rp**(2*l + 1) - 2*R2**(2*l + 4)*rp**(2*l + 1) + 4*R2**(4*l + 3)*rp**2)/(2*(-16*R1*R2**(4*l + 3)*l**2 + 4*R1*R2**(4*l + 3) + 16*R1**(2*l)*R2**(2*l + 4)*l**4 + 16*R1**(2*l)*R2**(2*l + 4)*l**3 - 4*R1**(2*l)*R2**(2*l + 4)*l - R1**(2*l)*R2**(2*l + 4) - 32*R1**(2*l + 2)*R2**(2*l + 2)*l**4 - 32*R1**(2*l + 2)*R2**(2*l + 2)*l**3 + 32*R1**(2*l + 2)*R2**(2*l + 2)*l**2 + 8*R1**(2*l + 2)*R2**(2*l + 2)*l - 6*R1**(2*l + 2)*R2**(2*l + 2) + 16*R1**(2*l + 4)*R2**(2*l)*l**4 + 16*R1**(2*l + 4)*R2**(2*l)*l**3 - 4*R1**(2*l + 4)*R2**(2*l)*l - R1**(2*l + 4)*R2**(2*l) - 16*R1**(4*l + 3)*R2*l**2 + 4*R1**(4*l + 3)*R2)),
      R2*Ra*rp**(-l + 1)*(4*R1*R2**(2*l + 1)*l*rp**(2*l + 3) - 2*R1*R2**(2*l + 1)*rp**(2*l + 3) - 4*R1*R2**(2*l + 3)*l*rp**(2*l + 1) - 2*R1*R2**(2*l + 3)*rp**(2*l + 1) + 4*R1**(2*l)*R2**(2*l + 3)*l**2*rp**2 + 4*R1**(2*l)*R2**(2*l + 3)*l*rp**2 + R1**(2*l)*R2**(2*l + 3)*rp**2 - 4*R1**(2*l + 2)*R2**(2*l + 1)*l**2*rp**2 - 4*R1**(2*l + 2)*R2**(2*l + 1)*l*rp**2 + 3*R1**(2*l + 2)*R2**(2*l + 1)*rp**2 - 4*R1**(2*l + 2)*R2**(2*l + 3)*l**2 + R1**(2*l + 2)*R2**(2*l + 3) - 4*R1**(2*l + 2)*l*rp**(2*l + 3) + 2*R1**(2*l + 2)*rp**(2*l + 3) + 4*R1**(2*l + 4)*R2**(2*l + 1)*l**2 - R1**(2*l + 4)*R2**(2*l + 1) + 4*R1**(2*l + 4)*l*rp**(2*l + 1) + 2*R1**(2*l + 4)*rp**(2*l + 1) - 4*R1**(4*l + 3)*rp**2)/(2*(-16*R1*R2**(4*l + 3)*l**2 + 4*R1*R2**(4*l + 3) + 16*R1**(2*l)*R2**(2*l + 4)*l**4 + 16*R1**(2*l)*R2**(2*l + 4)*l**3 - 4*R1**(2*l)*R2**(2*l + 4)*l - R1**(2*l)*R2**(2*l + 4) - 32*R1**(2*l + 2)*R2**(2*l + 2)*l**4 - 32*R1**(2*l + 2)*R2**(2*l + 2)*l**3 + 32*R1**(2*l + 2)*R2**(2*l + 2)*l**2 + 8*R1**(2*l + 2)*R2**(2*l + 2)*l - 6*R1**(2*l + 2)*R2**(2*l + 2) + 16*R1**(2*l + 4)*R2**(2*l)*l**4 + 16*R1**(2*l + 4)*R2**(2*l)*l**3 - 4*R1**(2*l + 4)*R2**(2*l)*l - R1**(2*l + 4)*R2**(2*l) - 16*R1**(4*l + 3)*R2*l**2 + 4*R1**(4*l + 3)*R2)),
    ])
B = numpy.array([
      R1**(2*l + 2)*Ra*rp**(-l + 1)*(-4*R1**2*R2**(2*l)*l**2*rp**(2*l + 3) - 4*R1**2*R2**(2*l)*l*rp**(2*l + 3) - R1**2*R2**(2*l)*rp**(2*l + 3) + 4*R1**2*R2**(2*l + 2)*l**2*rp**(2*l + 1) + 8*R1**2*R2**(2*l + 2)*l*rp**(2*l + 1) + 3*R1**2*R2**(2*l + 2)*rp**(2*l + 1) - 4*R1**2*R2**(4*l + 3)*l - 2*R1**2*R2**(4*l + 3) + 4*R1**(2*l + 1)*R2*rp**(2*l + 3) - 4*R1**(2*l + 1)*R2**(2*l + 2)*l*rp**2 - 6*R1**(2*l + 1)*R2**(2*l + 2)*rp**2 + 4*R1**(2*l + 1)*R2**(2*l + 4)*l + 2*R1**(2*l + 1)*R2**(2*l + 4) + 4*R2**(2*l + 2)*l**2*rp**(2*l + 3) + 4*R2**(2*l + 2)*l*rp**(2*l + 3) - 3*R2**(2*l + 2)*rp**(2*l + 3) - 4*R2**(2*l + 4)*l**2*rp**(2*l + 1) - 8*R2**(2*l + 4)*l*rp**(2*l + 1) - 3*R2**(2*l + 4)*rp**(2*l + 1) + 4*R2**(4*l + 3)*l*rp**2 + 6*R2**(4*l + 3)*rp**2)/(2*(-16*R1*R2**(4*l + 3)*l**2 - 32*R1*R2**(4*l + 3)*l - 12*R1*R2**(4*l + 3) + 16*R1**(2*l)*R2**(2*l + 4)*l**4 + 48*R1**(2*l)*R2**(2*l + 4)*l**3 + 48*R1**(2*l)*R2**(2*l + 4)*l**2 + 20*R1**(2*l)*R2**(2*l + 4)*l + 3*R1**(2*l)*R2**(2*l + 4) - 32*R1**(2*l + 2)*R2**(2*l + 2)*l**4 - 96*R1**(2*l + 2)*R2**(2*l + 2)*l**3 - 64*R1**(2*l + 2)*R2**(2*l + 2)*l**2 + 24*R1**(2*l + 2)*R2**(2*l + 2)*l + 18*R1**(2*l + 2)*R2**(2*l + 2) + 16*R1**(2*l + 4)*R2**(2*l)*l**4 + 48*R1**(2*l + 4)*R2**(2*l)*l**3 + 48*R1**(2*l + 4)*R2**(2*l)*l**2 + 20*R1**(2*l + 4)*R2**(2*l)*l + 3*R1**(2*l + 4)*R2**(2*l) - 16*R1**(4*l + 3)*R2*l**2 - 32*R1**(4*l + 3)*R2*l - 12*R1**(4*l + 3)*R2)),
      R2**(2*l + 2)*Ra*rp**(-l + 1)*(-4*R1*R2**(2*l + 1)*rp**(2*l + 3) + 4*R1**(2*l)*R2**2*l**2*rp**(2*l + 3) + 4*R1**(2*l)*R2**2*l*rp**(2*l + 3) + R1**(2*l)*R2**2*rp**(2*l + 3) - 4*R1**(2*l + 2)*R2**2*l**2*rp**(2*l + 1) - 8*R1**(2*l + 2)*R2**2*l*rp**(2*l + 1) - 3*R1**(2*l + 2)*R2**2*rp**(2*l + 1) + 4*R1**(2*l + 2)*R2**(2*l + 1)*l*rp**2 + 6*R1**(2*l + 2)*R2**(2*l + 1)*rp**2 - 4*R1**(2*l + 2)*l**2*rp**(2*l + 3) - 4*R1**(2*l + 2)*l*rp**(2*l + 3) + 3*R1**(2*l + 2)*rp**(2*l + 3) - 4*R1**(2*l + 4)*R2**(2*l + 1)*l - 2*R1**(2*l + 4)*R2**(2*l + 1) + 4*R1**(2*l + 4)*l**2*rp**(2*l + 1) + 8*R1**(2*l + 4)*l*rp**(2*l + 1) + 3*R1**(2*l + 4)*rp**(2*l + 1) + 4*R1**(4*l + 3)*R2**2*l + 2*R1**(4*l + 3)*R2**2 - 4*R1**(4*l + 3)*l*rp**2 - 6*R1**(4*l + 3)*rp**2)/(2*(-16*R1*R2**(4*l + 3)*l**2 - 32*R1*R2**(4*l + 3)*l - 12*R1*R2**(4*l + 3) + 16*R1**(2*l)*R2**(2*l + 4)*l**4 + 48*R1**(2*l)*R2**(2*l + 4)*l**3 + 48*R1**(2*l)*R2**(2*l + 4)*l**2 + 20*R1**(2*l)*R2**(2*l + 4)*l + 3*R1**(2*l)*R2**(2*l + 4) - 32*R1**(2*l + 2)*R2**(2*l + 2)*l**4 - 96*R1**(2*l + 2)*R2**(2*l + 2)*l**3 - 64*R1**(2*l + 2)*R2**(2*l + 2)*l**2 + 24*R1**(2*l + 2)*R2**(2*l + 2)*l + 18*R1**(2*l + 2)*R2**(2*l + 2) + 16*R1**(2*l + 4)*R2**(2*l)*l**4 + 48*R1**(2*l + 4)*R2**(2*l)*l**3 + 48*R1**(2*l + 4)*R2**(2*l)*l**2 + 20*R1**(2*l + 4)*R2**(2*l)*l + 3*R1**(2*l + 4)*R2**(2*l) - 16*R1**(4*l + 3)*R2*l**2 - 32*R1**(4*l + 3)*R2*l - 12*R1**(4*l + 3)*R2)),
    ])
C = numpy.array([
      Ra*rp**(-l + 1)*(-4*R1*R2**(2*l)*l*rp**(2*l + 3) - 2*R1*R2**(2*l)*rp**(2*l + 3) + 4*R1*R2**(2*l + 2)*l*rp**(2*l + 1) + 6*R1*R2**(2*l + 2)*rp**(2*l + 1) - 4*R1*R2**(4*l + 3) + 4*R1**(2*l)*R2*l*rp**(2*l + 3) + 2*R1**(2*l)*R2*rp**(2*l + 3) - 4*R1**(2*l)*R2**(2*l + 2)*l**2*rp**2 - 8*R1**(2*l)*R2**(2*l + 2)*l*rp**2 - 3*R1**(2*l)*R2**(2*l + 2)*rp**2 + 4*R1**(2*l)*R2**(2*l + 4)*l**2 + 4*R1**(2*l)*R2**(2*l + 4)*l + R1**(2*l)*R2**(2*l + 4) - 4*R1**(2*l + 2)*R2*l*rp**(2*l + 1) - 6*R1**(2*l + 2)*R2*rp**(2*l + 1) + 4*R1**(2*l + 2)*R2**(2*l)*l**2*rp**2 + 8*R1**(2*l + 2)*R2**(2*l)*l*rp**2 + 3*R1**(2*l + 2)*R2**(2*l)*rp**2 - 4*R1**(2*l + 2)*R2**(2*l + 2)*l**2 - 4*R1**(2*l + 2)*R2**(2*l + 2)*l + 3*R1**(2*l + 2)*R2**(2*l + 2))/(-32*R1*R2**(4*l + 3)*l**2 - 64*R1*R2**(4*l + 3)*l - 24*R1*R2**(4*l + 3) + 32*R1**(2*l)*R2**(2*l + 4)*l**4 + 96*R1**(2*l)*R2**(2*l + 4)*l**3 + 96*R1**(2*l)*R2**(2*l + 4)*l**2 + 40*R1**(2*l)*R2**(2*l + 4)*l + 6*R1**(2*l)*R2**(2*l + 4) - 64*R1**(2*l + 2)*R2**(2*l + 2)*l**4 - 192*R1**(2*l + 2)*R2**(2*l + 2)*l**3 - 128*R1**(2*l + 2)*R2**(2*l + 2)*l**2 + 48*R1**(2*l + 2)*R2**(2*l + 2)*l + 36*R1**(2*l + 2)*R2**(2*l + 2) + 32*R1**(2*l + 4)*R2**(2*l)*l**4 + 96*R1**(2*l + 4)*R2**(2*l)*l**3 + 96*R1**(2*l + 4)*R2**(2*l)*l**2 + 40*R1**(2*l + 4)*R2**(2*l)*l + 6*R1**(2*l + 4)*R2**(2*l) - 32*R1**(4*l + 3)*R2*l**2 - 64*R1**(4*l + 3)*R2*l - 24*R1**(4*l + 3)*R2),
     -Ra*rp**(-l + 1)*(4*R1*R2**(2*l)*l*rp**(2*l + 3) + 2*R1*R2**(2*l)*rp**(2*l + 3) - 4*R1*R2**(2*l + 2)*l*rp**(2*l + 1) - 6*R1*R2**(2*l + 2)*rp**(2*l + 1) - 4*R1**(2*l)*R2*l*rp**(2*l + 3) - 2*R1**(2*l)*R2*rp**(2*l + 3) + 4*R1**(2*l)*R2**(2*l + 2)*l**2*rp**2 + 8*R1**(2*l)*R2**(2*l + 2)*l*rp**2 + 3*R1**(2*l)*R2**(2*l + 2)*rp**2 + 4*R1**(2*l + 2)*R2*l*rp**(2*l + 1) + 6*R1**(2*l + 2)*R2*rp**(2*l + 1) - 4*R1**(2*l + 2)*R2**(2*l)*l**2*rp**2 - 8*R1**(2*l + 2)*R2**(2*l)*l*rp**2 - 3*R1**(2*l + 2)*R2**(2*l)*rp**2 - 4*R1**(2*l + 2)*R2**(2*l + 2)*l**2 - 4*R1**(2*l + 2)*R2**(2*l + 2)*l + 3*R1**(2*l + 2)*R2**(2*l + 2) + 4*R1**(2*l + 4)*R2**(2*l)*l**2 + 4*R1**(2*l + 4)*R2**(2*l)*l + R1**(2*l + 4)*R2**(2*l) - 4*R1**(4*l + 3)*R2)/(-32*R1*R2**(4*l + 3)*l**2 - 64*R1*R2**(4*l + 3)*l - 24*R1*R2**(4*l + 3) + 32*R1**(2*l)*R2**(2*l + 4)*l**4 + 96*R1**(2*l)*R2**(2*l + 4)*l**3 + 96*R1**(2*l)*R2**(2*l + 4)*l**2 + 40*R1**(2*l)*R2**(2*l + 4)*l + 6*R1**(2*l)*R2**(2*l + 4) - 64*R1**(2*l + 2)*R2**(2*l + 2)*l**4 - 192*R1**(2*l + 2)*R2**(2*l + 2)*l**3 - 128*R1**(2*l + 2)*R2**(2*l + 2)*l**2 + 48*R1**(2*l + 2)*R2**(2*l + 2)*l + 36*R1**(2*l + 2)*R2**(2*l + 2) + 32*R1**(2*l + 4)*R2**(2*l)*l**4 + 96*R1**(2*l + 4)*R2**(2*l)*l**3 + 96*R1**(2*l + 4)*R2**(2*l)*l**2 + 40*R1**(2*l + 4)*R2**(2*l)*l + 6*R1**(2*l + 4)*R2**(2*l) - 32*R1**(4*l + 3)*R2*l**2 - 64*R1**(4*l + 3)*R2*l - 24*R1**(4*l + 3)*R2),
    ])
D = numpy.array([
     -R1**(2*l)*Ra*rp**(-l + 1)*(4*R1**2*R2**(2*l)*l**2*rp**(2*l + 3) - R1**2*R2**(2*l)*rp**(2*l + 3) - 4*R1**2*R2**(2*l + 2)*l**2*rp**(2*l + 1) - 4*R1**2*R2**(2*l + 2)*l*rp**(2*l + 1) + 3*R1**2*R2**(2*l + 2)*rp**(2*l + 1) + 4*R1**2*R2**(4*l + 3)*l - 2*R1**2*R2**(4*l + 3) - 4*R1**(2*l + 3)*R2*rp**(2*l + 1) + 4*R1**(2*l + 3)*R2**(2*l)*l*rp**2 + 2*R1**(2*l + 3)*R2**(2*l)*rp**2 - 4*R1**(2*l + 3)*R2**(2*l + 2)*l + 2*R1**(2*l + 3)*R2**(2*l + 2) - 4*R2**(2*l + 2)*l**2*rp**(2*l + 3) + R2**(2*l + 2)*rp**(2*l + 3) + 4*R2**(2*l + 4)*l**2*rp**(2*l + 1) + 4*R2**(2*l + 4)*l*rp**(2*l + 1) + R2**(2*l + 4)*rp**(2*l + 1) - 4*R2**(4*l + 3)*l*rp**2 - 2*R2**(4*l + 3)*rp**2)/(32*R1*R2**(4*l + 3)*l**2 - 8*R1*R2**(4*l + 3) - 32*R1**(2*l)*R2**(2*l + 4)*l**4 - 32*R1**(2*l)*R2**(2*l + 4)*l**3 + 8*R1**(2*l)*R2**(2*l + 4)*l + 2*R1**(2*l)*R2**(2*l + 4) + 64*R1**(2*l + 2)*R2**(2*l + 2)*l**4 + 64*R1**(2*l + 2)*R2**(2*l + 2)*l**3 - 64*R1**(2*l + 2)*R2**(2*l + 2)*l**2 - 16*R1**(2*l + 2)*R2**(2*l + 2)*l + 12*R1**(2*l + 2)*R2**(2*l + 2) - 32*R1**(2*l + 4)*R2**(2*l)*l**4 - 32*R1**(2*l + 4)*R2**(2*l)*l**3 + 8*R1**(2*l + 4)*R2**(2*l)*l + 2*R1**(2*l + 4)*R2**(2*l) + 32*R1**(4*l + 3)*R2*l**2 - 8*R1**(4*l + 3)*R2),
      R2**(2*l)*Ra*rp**(-l + 1)*(-4*R1*R2**(2*l + 3)*rp**(2*l + 1) + 4*R1**(2*l)*R2**2*l**2*rp**(2*l + 3) - R1**(2*l)*R2**2*rp**(2*l + 3) + 4*R1**(2*l)*R2**(2*l + 3)*l*rp**2 + 2*R1**(2*l)*R2**(2*l + 3)*rp**2 - 4*R1**(2*l + 2)*R2**2*l**2*rp**(2*l + 1) - 4*R1**(2*l + 2)*R2**2*l*rp**(2*l + 1) + 3*R1**(2*l + 2)*R2**2*rp**(2*l + 1) - 4*R1**(2*l + 2)*R2**(2*l + 3)*l + 2*R1**(2*l + 2)*R2**(2*l + 3) - 4*R1**(2*l + 2)*l**2*rp**(2*l + 3) + R1**(2*l + 2)*rp**(2*l + 3) + 4*R1**(2*l + 4)*l**2*rp**(2*l + 1) + 4*R1**(2*l + 4)*l*rp**(2*l + 1) + R1**(2*l + 4)*rp**(2*l + 1) + 4*R1**(4*l + 3)*R2**2*l - 2*R1**(4*l + 3)*R2**2 - 4*R1**(4*l + 3)*l*rp**2 - 2*R1**(4*l + 3)*rp**2)/(32*R1*R2**(4*l + 3)*l**2 - 8*R1*R2**(4*l + 3) - 32*R1**(2*l)*R2**(2*l + 4)*l**4 - 32*R1**(2*l)*R2**(2*l + 4)*l**3 + 8*R1**(2*l)*R2**(2*l + 4)*l + 2*R1**(2*l)*R2**(2*l + 4) + 64*R1**(2*l + 2)*R2**(2*l + 2)*l**4 + 64*R1**(2*l + 2)*R2**(2*l + 2)*l**3 - 64*R1**(2*l + 2)*R2**(2*l + 2)*l**2 - 16*R1**(2*l + 2)*R2**(2*l + 2)*l + 12*R1**(2*l + 2)*R2**(2*l + 2) - 32*R1**(2*l + 4)*R2**(2*l)*l**4 - 32*R1**(2*l + 4)*R2**(2*l)*l**3 + 8*R1**(2*l + 4)*R2**(2*l)*l + 2*R1**(2*l + 4)*R2**(2*l) + 32*R1**(4*l + 3)*R2*l**2 - 8*R1**(4*l + 3)*R2)
    ])

def Pl(r):
    # eqn (105) from [1]
    return A*r**l + B*r**(-l-1) + C*r**(l+2) + D*r**(-l+1)

def dPldr(r):
    return l*A*r**(l-1) + (-l-1)*B*r**(-l-2) + (l+2)*C*r**(l+1) + (-l+1)*D*r**-l

# some sanity checks:
# eqn (108) from [1]:
assert abs(Pl(rp)[0]-Pl(rp)[1])<1e-12
assert abs(dPldr(rp)[0]-dPldr(rp)[1])<1e-12
# eqn (106):
assert all(numpy.abs(Pl(R))<1e-12)
# no slip also requires:
assert all(numpy.abs(dPldr(R))<1e-12)

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

def get_cartesian_solution(X, i):
    # i==0: "Upper Mantle" (above anomaly)
    # i==1: below anomaly
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    if theta<1e-6*pi or theta>pi*(1.-1e-6):
        return 0., 0., 0.
    phi = atan2(X[1], X[0])
    ur = u_r(r,theta, phi)[i]
    uth = u_theta(r,theta, phi)[i]
    uph = u_phi(r, theta, phi)[i]
    costh = cos(theta)
    req = sqrt(X[0]**2+X[1]**2)
    return ( X[0]/r*ur + X[0]/req*costh*uth - X[1]/req*uph,
             X[1]/r*ur + X[1]/req*costh*uth + X[0]/req*uph,
             X[2]/r*ur - sin(theta)*uth )

def get_spherical_solution(X, i):
    # i==0: "Upper Mantle" (above anomaly)
    # i==1: below anomaly
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    if theta<1e-6*pi or theta>pi*(1.-1e-6):
        return 0., 0., 0.
    phi = atan2(X[1], X[0])
    ur = u_r(r,theta, phi)[i]
    uth = u_theta(r,theta, phi)[i]
    uph = u_phi(r, theta, phi)[i]
    return ur, uth, uph

def Y_cartesian(X):
    r = sqrt(X[0]**2+X[1]**2+X[2]**2)
    theta = acos(X[2]/r)
    phi = atan2(X[1], X[0])
    return Y(m, l, theta, phi)
