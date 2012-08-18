from math import pi, exp

# prognostic variable

def gaussian(XX):
   '''Gaussian bump'''
   x = XX[0]; y = XX[1]
   r=(x**2+y**2)**.5
   u=exp(-pi*r**2)
   return u

def gaussian_advected(XX):
   '''Gaussian bump advected for 1 sec'''
   x = XX[0]; y = XX[1]
   u = velocity(XX)[0]
   v = velocity(XX)[1]
   t = 1.0
   r=((x-u*t)**2+(y-v*t)**2)**.5
   u=exp(-pi*r**2)
   return u

def gaussian_1x1(XX):
   '''Gaussian bump'''
   x = XX[0]-0.5; y = XX[1]-0.5
   r=(x**2+y**2)**.5
   u=exp(-pi*3.*r**2)
   return u

def velocity(XX):
   x = XX[0]; y = XX[1]
   u = 0.0
   v = 0.0
   return [u, v]
