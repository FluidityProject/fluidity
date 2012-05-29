from math import sin, cos

# simulated variables

def velocity(XX):
   '''Velocity'''
   x = XX[0]; y = XX[1]
   u=sin(x)*cos(y)
   v=-1.0*cos(x)*sin(y)
   return [u, v]

def pressure(XX):
   x = XX[0]; y = XX[1]
   return cos(x)*cos(y)-1.0

# gradients

def grad_u(XX):
   '''grad u'''
   x = XX[0]; y = XX[1]
   u_x = cos(x)*cos(y)
   u_y = -1.0*sin(x)*sin(y)
   return [u_x, u_y]

def grad_v(XX):
   '''grad v'''
   x = XX[0]; y = XX[1]
   v_x = sin(x)*sin(y)
   v_y = -1.0*cos(x)*cos(y)
   return [v_x, v_y]

def grad_p(XX):
   '''grad pressure'''
   x = XX[0]; y = XX[1]
   p_x=-1.0*sin(x)*cos(y)
   p_y=-1.0*cos(x)*sin(y)
   return [p_x, p_y]

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

def forcing_mom(XX):
   '''Forcing function: momentum'''
   x = XX[0]; y = XX[1]
   adv = 1.0; beta = 0.0; mass = 0.0; rho = 1.0
   nu = 0.7

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

   Su = mass*rho*u_t + adv*rho*(u*u_x+v*u_y) + p_x - nu*(u_xx+u_yy)
   Sv = mass*rho*u_t + adv*rho*(u*v_x+v*v_y) + p_y - nu*(v_xx+v_yy)
   return (Su, Sv)

