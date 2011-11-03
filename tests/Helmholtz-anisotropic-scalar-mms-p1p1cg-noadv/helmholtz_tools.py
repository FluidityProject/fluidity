from math import pi, exp

# prognostic variable

def gaussian(XX):
   '''Gaussian bump'''
   x = XX[0]; y = XX[1]
   r=(x**2+y**2)**.5
   u=exp(-pi*r**2)
   return u

def velocity(XX):
   x = XX[0]; y = XX[1]
   return 0.0
