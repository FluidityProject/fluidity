from math import sin, cos, pi

# simulated variables

def velocity(XX):
   '''Velocity'''
   x = XX[0]; y = XX[1]
   u=2.*pi*(sin(pi*x))**3*sin(pi*y)*cos(pi*y)
   v=-3.*pi*(sin(pi*x))**2*cos(pi*x)*(sin(pi*y))**2
   return [u, v]

def pressure(XX):
   x = XX[0]; y = XX[1]
   return cos(pi*x)*cos(pi*y)

# gradients

def grad_u(XX):
   '''grad u'''
   x = XX[0]; y = XX[1]
   u_x = 6.*pi**2*(sin(pi*x))**2*cos(pi*x)*sin(pi*y)*cos(pi*y)
   u_y = -2.*pi**3*(sin(pi*x))**3*sin(pi*y)*cos(pi*y)
   return [u_x, u_y]

def grad_v(XX):
   '''grad v'''
   x = XX[0]; y = XX[1]
   v_x = 6.*pi**3*(sin(pi*x))**2*cos(pi*x)*(sin(pi*y))**2
   v_y = -6.*pi**2*(sin(pi*x))**2*cos(pi*x)*sin(pi*y)*cos(pi*y)
   return [v_x, v_y]

def grad_p(XX):
   '''grad pressure'''
   x = XX[0]; y = XX[1]
   p_x=-pi*sin(pi*x)+1
   p_y=-pi*cos(pi*x)+1
   return [p_x, p_y]

def grad2_u(XX):
   x = XX[0]; y = XX[1]
   u_xx = -12.*pi**4*(sin(pi*x))**2*cos(pi*x)*sin(pi*y)*cos(pi*y)
   u_yy = 2.*pi**5*(sin(pi*x))**3*sin(pi*y)*cos(pi*y)
   return [u_xx, u_yy]

def grad2_v(XX):
   x = XX[0]; y = XX[1]
   v_xx = -12.*pi**5*(sin(pi*x))**2*cos(pi*x)*(sin(pi*y))**2
   v_yy = 6.*pi**4*(sin(pi*x))**2*cos(pi*x)*sin(pi*y)*cos(pi*y)
   return [v_xx, v_yy]

def forcing_mom(XX):
   '''Forcing function: momentum'''
   x = XX[0]; y = XX[1]
   adv = 1.0; beta = 0.0; mass = 0.0; rho = 1.0
   nu = 1.0

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
   v_xx = grad2_v(XX)[0]
   v_yy = grad2_v(XX)[1]
   u_t = 0.0

   Su = mass*rho*u_t + adv*rho*(u*u_x+v*u_y) + p_x - nu*(u_xx+u_yy)
   Sv = mass*rho*u_t + adv*rho*(u*v_x+v*v_y) + p_y - nu*(v_xx+v_yy)
   return (Su, Sv)
