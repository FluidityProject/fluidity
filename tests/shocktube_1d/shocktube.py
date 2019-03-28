# bit of python to compute analytical solution for the shock tube problem
# solution taken from "Mathematical and Computational Methods for Compressible Flow"
# by Miloslav Feistauer, Jiri Felcman and Ivan Straskraba, 2003
# eqn numbers below refer to this books
from math import sqrt
from scipy.optimize import newton

# ratio of specific heats
gamma=1.4

# left state:
ul=0.
rhol=1.0
iel=2.5

# right state:
ur=0.
rhor=0.2
ier=2.5

# this eos is assumed in the eqns below (so can't change this):
def p_eos(ie, rho):
  return rho*ie*(gamma-1.0)
  
def rho_eos(ie, p):
  return p/ie/(gamma-1.0)


# derived quantities, left:
pl=p_eos(iel,rhol)
El=rhol*iel
al=sqrt(gamma*pl/rhol)
# derived quantities, right:
pr=p_eos(ier,rhor)
Er=rhor*ier
ar=sqrt(gamma*pr/rhor)


def F1l(p):
  """Function F1l defines difference between us and ul:
     us=ul+F1l(p)    eqn. (3.1.159)"""
     
  # eqn. (3.1.160)
  if p<=pl:
    return 2.*al/(gamma-1.)*(1.-(p/pl)**((gamma-1.)/(2.*gamma)))
  else:
    return -(p-pl)*sqrt((2./(gamma+1.)/rhol)/(p+(gamma-1.)/(gamma+1.)*pl))
    
def F1lprime(p):
  # derivative of the above
  if p<=pl:
    return -al/gamma*(p/pl)**((gamma-1.)/(2.*gamma)-1.)
  else:
    return   -1.*sqrt((2./(gamma+1.)/rhol)/(p+(gamma-1.)/(gamma+1.)*pl)) + \
      -(p-pl)/2./sqrt((2./(gamma+1.)/rhol)/(p+(gamma-1.)/(gamma+1.)*pl))

def F3r(p):
  """Function F3r defines difference between us and ur:
     us=ur+F3r(p)    eqn. (3.1.161)"""
     
  # eqn. (3.1.162)
  if p<=pr:
    return -2.*ar/(gamma-1.)*(1.-(p/pr)**((gamma-1.)/(2.*gamma)))
  else:
    return (p-pr)*sqrt((2./(gamma+1.)/rhor)/(p+(gamma-1.)/(gamma+1.)*pr))
    
def F3rprime(p):
  # derivative of the above
  if p<=pr:
    return ar/gamma*(p/pr)**((gamma-1.)/(2.*gamma)-1.)
  else:
    return      sqrt((2./(gamma+1.)/rhor)/(p+(gamma-1.)/(gamma+1.)*pr)) + \
      (p-pr)/2./sqrt((2./(gamma+1.)/rhor)/(p+(gamma-1.)/(gamma+1.)*pr))

def F(p):
  # eqn (3.1.165)
  return F3r(p)-F1l(p)+ur-ul
  
def Fprime(p):
  return F3rprime(p)-F1lprime(p)
  
# inital guess:
p=(pl+pr)/2.

# solve eqn (3.1.164): F(p*)=0, to find p* the pressure between the u-a and u+a waves
# (is constant over contact discontinuity at the u wave)
ps=newton(F, p, fprime=Fprime)
print("should be zero: F(p*) = ", F(ps))
print("p* =", ps)
us=ul+F1l(ps)
print("u* =", us)
print("should be equal to:", ur+F3r(ps))

# the star region is divide by the contact discontinuity which gives a discontinuity in
# density rhosl/rhosr and internal energy iesl/iesr

# left:
if ps<pl: # rarefaction
  rhosl=rhol*(ps/pl)**(1/gamma) # eqn (3.1.172)
else: # shock
  rhosl=rhol*( (gamma-1.)/(gamma+1.)*pl/ps +1. )/(pl/ps+(gamma-1.)/(gamma+1.)) # eqn (3.1.176)
  s1=ul-al*sqrt((gamma+1.)/2./gamma*ps/pl + (gamma-1.)/2./gamma) # shock speed, eqn (3.1.177)
iesl=ps/rhosl/(gamma-1.)
asl=sqrt(gamma*ps/rhosl)

# right:
if ps<pr: # rarefaction
  rhosr=rhor*(ps/pr)**(1/gamma) # eqn (3.1.178)
else: # shock
  rhosr=rhor*( ps/pr + (gamma-1.)/(gamma+1.) )/( (gamma-1.)/(gamma+1.)*ps/pr +1. ) # eqn (3.1.182)
  s3=ur+ar*sqrt((gamma+1.)/2./gamma*ps/pr + (gamma-1.)/2./gamma) # shock speed, eqn (3.1.183)
iesr=ps/rhosr/(gamma-1.)
asr=sqrt(gamma*pr/rhosr)

def solution(x,t):
  if x/t<us:    
    # before the contact discontinuity:
    
    if ps<pl:
      # u-a is a rarefaction wave
      if x/t<ul-al: # left
        return (pl, ul, rhol)
      elif x/t<us-asl: # within the rarefaction wave
        p=pl*(2./(gamma+1.) + (gamma-1.)/(gamma+1.)/al*(ul-x/t))**(2*gamma/(gamma-1.)) # eqn (3.1.97)
        u=2./(gamma+1.)*(al + (gamma-1.)/2.*ul + x/t) # eqn (3.1.97)
        rho=rhol*(2./(gamma+1.) + (gamma-1.)/(gamma+1.)/al*(ul-x/t))**(2./(gamma-1.)) # eqn (3.1.97)
        return (p, u, rho)
      else: # after the rarefaction before the contact disc.
        return (ps,us,rhosl)
    else:
      # u-a is a shock wave
      if x/t<s1: # left
        return (pl, ul, rhol)
      else: # between u-a shock and contact disc.
        return (psl, us, rhosl)
        
  else:
    # after the contact discontinuity:
        
    if ps<pr:
      # u+a is a rarefaction wave
      if x/t>ur+ar: # right
        return (pr, ur, rhor)
      elif x/t>usr+asr: # within the rarefaction wave
        p=pr*(2./(gamma+1.) + (gamma-1.)/(gamma+1.)/ar*(ur-x/t))**(2*gamma/(gamma-1.)) # eqn (3.1.98)
        u=2./(gamma+1.)*(-ar + (gamma-1.)/2.*ur + x/t) # eqn (3.1.98)
        rho=rhor*(2./(gamma+1.) - (gamma-1.)/(gamma+1.)/ar*(ur-x/t))**(2./(gamma-1.)) # eqn (3.1.98)
        return (p, u, rho)
      else: # after the contact disc. before the rarefaction
        return (ps,us, rhosr)
    else:
      # u+a is a shock wave
      if x/t>s3: # right
        return (pr, ur, rhor)
      else: # between contact disc. and u-a shock
        return (ps, us, rhosr)

