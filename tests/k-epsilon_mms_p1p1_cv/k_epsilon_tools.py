from math import sin, cos, tanh

def velocity(XX):
   '''Velocity'''
   x = XX[0]; y = XX[1]
   u= 0.1*sin(x)*cos(y)
   v= -0.1*cos(x)*sin(y)
   return [u, v]

def pressure(XX):
   x = XX[0]; y = XX[1]
   return cos(x)*cos(y)-1.0

def rho_p(XX):
   '''Density perturbation'''
   x = XX[0]; y = XX[1]
   return 0.1*cos(x)*cos(y)
   
def forcing_mom(XX):
   '''Forcing function: momentum'''
   x = XX[0]; y = XX[1]

   Su = (0.01*sin(x)*sin(y)**2*cos(x) 
         + 0.01*sin(x)*cos(x)*cos(y)**2 
         - 0.8*sin(x)*cos(y))
   Sv = (0.01*sin(x)**2*sin(y)*cos(y) 
         + 0.01*sin(y)*cos(x)**2*cos(y) 
         - 1.2*sin(y)*cos(x) 
         + 0.1*cos(x)*cos(y))
   return (Su, Sv)

def forcing_density(XX):
   '''Forcing function: temperature field'''
   x = XX[0]; y = XX[1]

   Srho_p = -0.01*sin(x)**2*cos(y)**2 + 0.01*sin(y)**2*cos(x)**2
   return Srho_p 

def forcing_k(XX):
   '''Forcing function: k'''
   x = XX[0]; y = XX[1]

#    S = 2.*(u_x**2 + v_y**2 + u_y*v_x) + u_y**2 + v_x**2

#    Pk = nut(XX)*S
#    Gk = nut(XX)/sigmaT*(rho_p_x*g[0] + rho_p_y*g[1])
#    Ak = e
#    Sk = (k_t 
#          + (u*k_x+v*k_y) 
#          - 1./sigmak*(nut(XX)*(k_xx+k_yy) + (nut_x*k_x+nut_y*k_y)) 
#          - Pk - Gk + Ak)
#    return Sk

def forcing_eps(XX):
   '''Forcing function: epsilon'''
   x = XX[0]; y = XX[1]

#    S = 2.*(u_x**2 + v_y**2 + u_y*v_x) + u_y**2 + v_x**2
#    Pe = e/k*nut(XX)*S
#    Ge = e/k*nut(XX)/sigmaT*(rho_p_x*g[0] + rho_p_y*g[1])
#    Ae = e**2/k
#    ce3 = tanh(u_y/u_x)

#    Se = (eps_t 
#          + (u*eps_x+v*eps_y) 
#          - 1./sigmae*(nut(XX)*(eps_xx+eps_yy) + (nut_x*eps_x+nut_y*eps_y))
#          - ce1*Pe - ce1*ce3*Ge + ce2*Ae)
#    return Se
