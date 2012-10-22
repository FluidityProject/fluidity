y = var('y')
from math import pi

def function(phi_0, phi_x, phi_y, phi_xy, 
             f_sin_x, f_cos_x, f_sin_y, f_cos_y, f_sin_xy, f_cos_xy, 
             alpha_x, alpha_y, alpha_xy):
    
    f_0 = phi_0 
    f_x = phi_x*(f_sin_x*sin(alpha_x*x) + f_cos_x*cos(alpha_x*x)) 
    f_y = phi_y*(f_sin_y*sin(alpha_y*y) + f_cos_y*cos(alpha_y*y)) 
    f_xy = phi_xy*(f_sin_xy*sin(alpha_xy*x*y/pi) + f_cos_xy*cos(alpha_xy*x*y/pi)) 
    f = f_0 + f_x + f_y + f_xy
    return f

u = function(3.0, 1.0, 0.6, 0.0, 
             0.0, 1.0, 1.0, 0.0, 1.0, 0.0,
             1.0, 1.0, 1.0)
v = function(2.0, 0.25, 1.5, 0.0, 
             0.0, 1.0, 5.0, 0.0, 2.0, 0.0,
             1.0, 0.0, 1.0)
             
u = 1.0*(sin(x**2+y**2)+0.5)
v = 0.1*(cos(x**2+y**2)+0.5)

rho = 0.5*(sin(x*x + y*y) + 1.5)
ie = 4.0*x**2 + 1.0
gamma = 1.4

p = (gamma-1.0)*ie*rho

nu = 0.7

# Smagorinsky model
Cs = 0.1
n = 8
Delta2 = 0.0771062843836 # Delta squared - the area of each element in MMS_A.msh.

S_xx = diff(u,x)
S_xy = 0.5*(diff(u,y) + diff(v,x))
S_yy = diff(v,y)
S_yx = 0.5*(diff(u,y) + diff(v,x))
S_norm = sqrt(2*(S_xx**2 + S_xy**2 + S_yx**2 + S_yy**2))

nu_T = 2*rho*((Cs)**2)*S_norm*Delta2

#print "DIVERGENCE = ", diff(u,x) + diff(v,y)

#nu = nu+nu_T

#tau_xx = 2*nu*diff(u,x) - (2.0/3.0)*nu*(diff(u,x) + diff(v,y))
#tau_xy = nu*(diff(u,y) + diff(v,x))
#tau_yy = 2*nu*diff(v,y) - (2.0/3.0)*nu*(diff(u,x) + diff(v,y)) 
#tau_yx = nu*(diff(u,y) + diff(v,x))

tau_xx = nu*diff(u,x)            
tau_yy = nu*diff(v,y)
tau_xy = nu*diff(u,y)
tau_yx = nu*diff(v,x)  

Su = rho*u*diff(u,x) + rho*v*diff(u,y) - diff(tau_xx, x) - diff(tau_xy, y) + diff(p,x)  
Sv = rho*u*diff(v,x) + rho*v*diff(v,y) - diff(tau_yx, x) - diff(tau_yy, y) + diff(p,y)

Sie = rho*u*diff(ie,x) + rho*v*diff(ie,y) + p*(diff(u,x) + diff(v,y))

Srho = diff(u*rho,x) + diff(v*rho,y)
  
print 'from math import sin, cos, tanh, pi, sqrt'
print ''
print 'def u(X):'
print '    return', str(u).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v(X):'
print '    return', str(v).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def p(X):'
print '    return', str(p).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print '' 
print 'def rho(X):'
print '    return', str(rho).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def ie(X):'
print '    return', str(ie).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u(X):'
print '    return', str(Su).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v(X):'
print '    return', str(Sv).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_rho(X):'
print '    return', str(Srho).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_ie(X):'
print '    return', str(Sie).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity(X):'
print '   return [u(X), v(X)]'
print ''
print 'def forcing_velocity(X):'
print '   return [forcing_u(X), forcing_v(X)]'
print ''
