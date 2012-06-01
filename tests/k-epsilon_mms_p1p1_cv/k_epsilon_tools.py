from math import sin, cos, tanh

g = [0,-1.0]
nu = 1.0
sigmaT = 1.0
sigmak = 1.0
sigmae = 1.3
ce1 = 1.44 
ce2 = 1.92

# simulated variables

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
   
def tke(XX):
   '''Turbulent kinetic energy'''
   x = XX[0]; y = XX[1]
   return y+1.0
   
def eps(XX):
   '''Turbulent dissipation'''
   x = XX[0]; y = XX[1]
   return y+1.0
   
def nut(XX):
   x = XX[0]; y = XX[1]
   return 0.09*(y+1.0)
   
# gradients

def grad_u(XX):
   '''grad u'''
   x = XX[0]; y = XX[1]
   u_x = 0.1*cos(x)*cos(y)
   u_y = -0.1*sin(x)*sin(y)
   return [u_x, u_y]

def grad_v(XX):
   '''grad v'''
   x = XX[0]; y = XX[1]
   v_x = 0.1*sin(x)*sin(y)
   v_y = -0.1*cos(x)*cos(y)
   return [v_x, v_y]

def grad_p(XX):
   '''grad pressure'''
   x = XX[0]; y = XX[1]
   p_x= -1.0*sin(x)*cos(y)
   p_y= -1.0*cos(x)*sin(y)
   return [p_x, p_y]

def grad_rho_p(XX):
   '''grad density perturbation'''
   x = XX[0]; y = XX[1]
   p_x= -0.1*sin(x)*cos(y)
   p_y= -0.1*cos(x)*sin(y)
   return [p_x, p_y]
   
def grad_tke(XX):
   '''grad turbulent kinetic energy'''
   x = XX[0]; y = XX[1]
   k_x=0.0
   k_y=1.0
   return [k_x, k_y]
   
def grad_eps(XX):
   '''grad turbulent dissipation'''
   x = XX[0]; y = XX[1]
   eps_x=0.0
   eps_y=1.0
   return  [eps_x, eps_y]

def grad_nut(XX):
   x = XX[0]; y = XX[1]
   nut_x = 0.0
   nut_y = 0.09
   return  [nut_x, nut_y]

def grad2_u(XX):
   x = XX[0]; y = XX[1]
   u_xx = -1.0*velocity(XX)[0]
   u_yy = -1.0*velocity(XX)[0]
   u_xy = velocity(XX)[1]
   return [u_xx, u_yy, u_xy]

def grad2_v(XX):
   x = XX[0]; y = XX[1]
   v_xx = -1.0*velocity(XX)[1]
   v_yy = -1.0*velocity(XX)[1]
   v_xy = velocity(XX)[0]
   return [v_xx, v_yy, v_xy]

def grad2_rho_p(XX):
   x = XX[0]; y = XX[1]
   rho_p_xx = -rho_p(XX)
   rho_p_yy = -rho_p(XX)
   rho_p_xy = 0.1*sin(x)*sin(y)
   return  [rho_p_xx, rho_p_yy, rho_p_xy]

def grad2_tke(XX):
   x = XX[0]; y = XX[1]
   k_xx = 0.0
   k_yy = 0.0
   return [k_xx, k_yy]

def grad2_eps(XX):
   x = XX[0]; y = XX[1]
   eps_xx = 0.0
   eps_yy = 0.0
   return [eps_xx, eps_yy]

def forcing_mom(XX):
   '''Forcing function: momentum: correct expansion?'''
   x = XX[0]; y = XX[1]

   u = velocity(XX)[0]
   v = velocity(XX)[1]
   p = pressure(XX)

   p_x = grad_p(XX)[0]
   p_y = grad_p(XX)[1]
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   u_t = 0.0

   nut_x = grad_nut(XX)[0]
   nut_y = grad_nut(XX)[1]

   u_xx = grad2_u(XX)[0]
   u_yy = grad2_u(XX)[1]
   u_xy = grad2_u(XX)[2]
   v_xx = grad2_v(XX)[0]
   v_yy = grad2_v(XX)[1]
   v_xy = grad2_v(XX)[2]

   #print "source terms x: ", adv*rho*(u*u_x+v*u_y), p_x, -1.*(nu+nut(XX))*(u_xx+u_yy), -2.*nut_x*u_x, -1.*nut_y*(u_y+v_x)
   #print "source terms y: ", adv*rho*(u*v_x+v*v_y), p_y, -1.*(nu+nut(XX))*(v_xx+v_yy), -2.*nut_y*v_y, -1.*nut_x*(u_y+v_x)
   # Expansion from Eca paper:
   #Su = mass*rho*u_t + adv*rho*(u*u_x+v*u_y) + p_x - (nu+nut(XX))*(u_xx+u_yy) - 2.*nut_x*u_x - nut_y*(u_y+v_x)
   #Sv = mass*rho*u_t + adv*rho*(u*v_x+v*v_y) + p_y - (nu+nut(XX))*(v_xx+v_yy) - 2.*nut_y*v_y - nut_x*(u_y+v_x)

   # Correct expansion of viscous terms. Eca paper was wrong.
   Su = (u_t 
         + (u*u_x+v*u_y) 
         + p_x 
         - (nu+nut(XX))*(v_xx+v_yy) + nut_x*u_x+nut_y*u_y 
         - (rho_p(XX))*g[0])
   Sv = (u_t 
         + (u*v_x+v*v_y) 
         + p_y 
         - (nu+nut(XX))*(v_xx+v_yy) + nut_x*v_x+nut_y*v_y 
         - (rho_p(XX))*g[1])
   return (Su, Sv)

def forcing_density(XX):
   '''Forcing function: temperature field'''
   x = XX[0]; y = XX[1]

   u = velocity(XX)[0]
   v = velocity(XX)[1]

   p_x = grad_p(XX)[0]
   p_y = grad_p(XX)[1]
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   rho_p_x = grad_rho_p(XX)[0]
   rho_p_y = grad_rho_p(XX)[1]
   rho_p_xx = grad2_rho_p(XX)[0]
   rho_p_yy = grad2_rho_p(XX)[1]
   rho_p_t = 0.0

   nut_x = grad_nut(XX)[0]
   nut_y = grad_nut(XX)[1]

   ST = (rho_p_t 
         + (u*rho_p_x+v*rho_p_y) 
         - 1./sigmaT*(nut(XX)*(rho_p_xx+rho_p_yy) + (nut_x*rho_p_x+nut_y*rho_p_y)))
   return ST

def forcing_k(XX):
   '''Forcing function: k'''
   x = XX[0]; y = XX[1]

   u = velocity(XX)[0]
   v = velocity(XX)[1]
   k = tke(XX)
   e = eps(XX)

   k_x = grad_tke(XX)[0]
   k_y = grad_tke(XX)[1]
   k_xx = grad2_tke(XX)[0]
   k_yy = grad2_tke(XX)[1]
   k_t = 0.0

   nut_x = grad_nut(XX)[0]
   nut_y = grad_nut(XX)[1]
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   rho_p_x = grad_rho_p(XX)[0]
   rho_p_y = grad_rho_p(XX)[1]

   S = 2.*(u_x**2 + v_y**2 + u_y*v_x) + u_y**2 + v_x**2

   Pk = nut(XX)*S
   Gk = nut(XX)/sigmaT*(rho_p_x*g[0] + rho_p_y*g[1])
   Ak = e
   Sk = (k_t 
         + (u*k_x+v*k_y) 
         - 1./sigmak*(nut(XX)*(k_xx+k_yy) + (nut_x*k_x+nut_y*k_y)) 
         - Pk - Gk + Ak)
   return Sk

def forcing_eps(XX):
   '''Forcing function: epsilon'''
   x = XX[0]; y = XX[1]

   u = velocity(XX)[0]
   v = velocity(XX)[1]

   # Need to prevent k=0:
   kmin = 1.e-6
   k = max(tke(XX), kmin)
   e = eps(XX)

   eps_x = grad_eps(XX)[0]
   eps_y = grad_eps(XX)[1]
   eps_xx = grad2_eps(XX)[0]
   eps_yy = grad2_eps(XX)[1]
   eps_t = 0.0

   nut_x = grad_nut(XX)[0]
   nut_y = grad_nut(XX)[1]
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   rho_p_x = grad_rho_p(XX)[0]
   rho_p_y = grad_rho_p(XX)[1]

   S = 2.*(u_x**2 + v_y**2 + u_y*v_x) + u_y**2 + v_x**2
   Pe = e/k*nut(XX)*S
   Ge = e/k*nut(XX)/sigmaT*(rho_p_x*g[0] + rho_p_y*g[1])
   Ae = e**2/k
   ce3 = tanh(u_y/u_x)

   Se = (eps_t 
         + (u*eps_x+v*eps_y) 
         - 1./sigmae*(nut(XX)*(eps_xx+eps_yy) + (nut_x*eps_x+nut_y*eps_y))
         - ce1*Pe - ce1*ce3*Ge + ce2*Ae)
   return Se

def forcing_mom_noturb(XX):
   '''Forcing function: momentum'''
   x = XX[0]; y = XX[1]

   u = velocity(XX)[0]
   v = velocity(XX)[1]
   p = pressure(XX)

   p_x = grad_p(XX)[0]
   p_y = grad_p(XX)[1]
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   u_xx = grad2_u(XX)[0]
   u_yy = grad2_u(XX)[1]
   u_xy = grad2_u(XX)[2]
   v_xx = grad2_v(XX)[0]
   v_yy = grad2_v(XX)[1]
   v_xy = grad2_v(XX)[2]
   u_t = 0.0

   Su = u_t + (u*u_x+v*u_y) + p_x - nu*(u_xx+u_yy) - rho_p(XX)*g[0]
   Sv = u_t + (u*v_x+v*v_y) + p_y - nu*(v_xx+v_yy) - rho_p(XX)*g[1]
   return (Su, Sv)
