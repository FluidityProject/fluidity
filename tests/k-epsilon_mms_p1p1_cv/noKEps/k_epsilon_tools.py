from math import sin, cos, tanh

def velocity(XX):
   '''Velocity'''
   x = XX[0]; y = XX[1]
   u= 0.1*sin(x)*cos(y)
   v= -0.1*cos(x)*sin(y)
   return [u, v]

def pressure(XX):
   x = XX[0]; y = XX[1]
   return 1.0-cos(x)*cos(y)

def rho_p(XX):
   '''Density perturbation'''
   x = XX[0]; y = XX[1]
   return 0.1*cos(x)*cos(y)
   
def forcing_mom(XX):
   '''Forcing function: momentum'''
   x = XX[0]; y = XX[1]

   Su = (0.0100*sin(x)*sin(y)**2*cos(x) + 0.0100*sin(x)*cos(x)*cos(y)**2 + 1.20*sin(x)*cos(y))
   Sv = (0.0100*sin(x)**2*sin(y)*cos(y) + 0.0100*sin(y)*cos(x)**2*cos(y) + 0.800*sin(y)*cos(x) + 0.100*cos(x)*cos(y))
   return (Su, Sv)

def forcing_density(XX):
   '''Forcing function: temperature field'''
   x = XX[0]; y = XX[1]

   Srho_p = (-0.0100*sin(x)**2*cos(y)**2 + 0.0100*sin(y)**2*cos(x)**2)
   return Srho_p 
