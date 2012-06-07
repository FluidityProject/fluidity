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

def tke(XX):
   '''Turbulent kinetic energy'''
   x = XX[0]; y = XX[1]
   return y+1.0

def eps(XX):
   '''Turbulent dissipation'''
   x = XX[0]; y = XX[1]
   return y+1.0
   
def forcing_mom(XX):
   '''Forcing function: momentum'''
   x = XX[0]; y = XX[1]

   Su = (0.0100*sin(x)*sin(y)**2*cos(x) + 0.0100*sin(x)*cos(x)*cos(y)**2 + 0.200*(y + 2.00)*sin(x)*cos(y) + 0.100*sin(x)*sin(y) + sin(x)*cos(y))
   Sv = (0.0100*sin(x)**2*sin(y)*cos(y) + 0.0100*sin(y)*cos(x)**2*cos(y) - 0.200*(y + 2.00)*sin(y)*cos(x) + sin(y)*cos(x) + 0.200*cos(x)*cos(y))
   return (Su, Sv)

def forcing_density(XX):
   '''Forcing function: temperature field'''
   x = XX[0]; y = XX[1]

   Srho_p = (-0.0100*sin(x)**2*cos(y)**2 + 0.0100*sin(y)**2*cos(x)**2 + 0.200*(y + 1.00)*cos(x)*cos(y) + 0.100*sin(y)*cos(x))
   return Srho_p 

def forcing_k(XX):
   '''Forcing function: k'''
   x = XX[0]; y = XX[1]

   Sk = (-0.0400*(y + 1.00)*cos(x)**2*cos(y)**2 + 0.100*(y + 1.00)*sin(y)*cos(x) - 0.100*sin(y)*cos(x) + y)
   return Sk

def forcing_eps(XX):
   '''Forcing function: epsilon'''
   x = XX[0]; y = XX[1]

   if sin(x)*cos(y) == 0:
      Ce3 = 1.0
   else:
      Ce3 = tanh(abs(-sin(y)*cos(x))/abs((sin(x)*cos(y))))
   Se = (-0.0400*(y + 1.00)*cos(x)**2*cos(y)**2 + 0.100*(y + 1.00)*sin(y)*cos(x)*Ce3 - 0.100*sin(y)*cos(x) + y)
   return Se

def Pke(XX):
   x = XX[0]; y = XX[1]

   Pke = 0.0400*(y + 1.00)*cos(x)**2*cos(y)**2
   return Pke

def Gke(XX):
   x = XX[0]; y = XX[1]

   Gke = -0.100*(y + 1.00)*sin(y)*cos(x)
   return Gke

def Ake(XX):
   x = XX[0]; y = XX[1]

   Ake = 1.00
   return Ake

def Peps(XX):
   x = XX[0]; y = XX[1]

   Peps = 0.0400*(y + 1.00)*cos(x)**2*cos(y)**2
   return Peps

def Geps(XX):
   x = XX[0]; y = XX[1]

   if sin(x)*cos(y) == 0:
      Ce3 = 1.0
   else:
      Ce3 = tanh(abs(-sin(y)*cos(x))/abs((sin(x)*cos(y))))

   Geps = Ce3*(-0.100*(y + 1.00)*sin(y)*cos(x))
   return Geps

def Aeps(XX):
   x = XX[0]; y = XX[1]

   Aeps = 1.00
   return Aeps
