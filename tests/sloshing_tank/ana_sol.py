import math
from scipy.special import erf
from numpy import poly1d
from numpy import pi, sin, linspace
from numpy import exp, cos

# scipy erfc does not support complex numbers but erf does
def erfc(z):
        return 1.0-erf(z)  


def analytic_solution_simple(t, x, a_0, L, g, eta): # time, x_value, initial maximum perturbation (at the very left), wavelength, gravity, viscosity
        k= 2*pi/L # wave number 2pi/wavelength 
        gamma=2*eta*k**2
        a=a_0*cos(pi*x)        
        return a*exp(-gamma*t)



# Formulas were compared with the same formula hacked into matlab
# based on Motion of two superposed viscous fluids, prosperetti
# assumes that fluid in the tank has same viscosity as fluid above
def analytic_solution(t, x, a_0, L, g, eta): # time, x_value, initial maximum perturbation (at the very left), wavelength, gravity, viscosity
  debug=False
  eta=eta/2.0 
  k= 2.0*pi/L # wave number 2pi/wavelength 
  omega_0sq = g*k  # inviscid natural frequency        
  a=a_0*cos(pi*x)        

  p1 = poly1d([1.0,0.0,2*k**2*eta,4*k**3*eta**(1.5),eta**2*k**4+omega_0sq],r=0)

  p1roots = p1.r
  z1=p1roots[0]
  z2=p1roots[1]
  z3=p1roots[2]
  z4=p1roots[3]
  Z1=(z2-z1)*(z3-z1)*(z4-z1)
  Z2=(z1-z2)*(z3-z2)*(z4-z2)
  Z3=(z1-z3)*(z2-z3)*(z4-z3)
  Z4=(z1-z4)*(z2-z4)*(z3-z4)

  if debug:
                print('Calculate analytic solution:'        )
                print('t=', t)
                print('a=', a)
                print('k=', k)
                print('g=', g)
                print('eta=', eta)
                print('omega_0sq=', omega_0sq)
        
  t0=4*eta**2*k**4/(8*eta**2*k**4+omega_0sq)*a*erfc((eta*k**2*t)**0.5)
  t1=z1/Z1*omega_0sq*a/(z1**2-eta*k**2)*exp((z1**2-eta*k**2)*t)*erfc(z1*t**0.5)
  t2=z2/Z2*omega_0sq*a/(z2**2-eta*k**2)*exp((z2**2-eta*k**2)*t)*erfc(z2*t**0.5)
  t3=z3/Z3*omega_0sq*a/(z3**2-eta*k**2)*exp((z3**2-eta*k**2)*t)*erfc(z3*t**0.5)
  t4=z4/Z4*omega_0sq*a/(z4**2-eta*k**2)*exp((z4**2-eta*k**2)*t)*erfc(z4*t**0.5)

  a=t0+t1+t2+t3+t4
        
  if debug:
                print('a(t)=', a.real        )
  if (a.imag>0.000001):
                print('Warning: Imaginary part of a(t) is not zero!')
  return a.real


