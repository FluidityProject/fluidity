y = var('y')

def function(phi_0, phi_x, phi_y, phi_xy, 
             f_sin_x, f_cos_x, f_sin_y, f_cos_y, f_sin_xy, f_cos_xy, 
             alpha_x, alpha_y, alpha_xy):
    
    f_0 = phi_0 
    f_x = phi_x*(f_sin_x*sin(alpha_x*x) + f_cos_x*cos(alpha_x*x)) 
    f_y = phi_y*(f_sin_y*sin(alpha_y*y) + f_cos_y*cos(alpha_y*y)) 
    f_xy = phi_xy*(f_sin_xy*sin(alpha_xy*x*y/pi) + f_cos_xy*cos(alpha_xy*x*y/pi)) 
    f = f_0 + f_x + f_y + f_xy
    return f
             
rho = 1.0 #0.5*(sin(x*x + y*y) + 1.5)

ie = 0.5*(cos(x + y) + 1.5)

u = 1.0*(sin(x**2+y**2)+0.5)
v = 0.1*(cos(x**2+y**2)+0.5)
 
gamma = 1.4
p = (gamma-1.0)*ie*rho

mu = 0.7

tau_xx = mu*diff(u,x)  #2*mu*diff(u,x) - (2.0/3.0)*mu*(diff(u,x) + diff(v,y))
tau_xy = mu*diff(u,y)  #mu*(diff(u,y) + diff(v,x))
tau_yy = mu*diff(v,y)  #2*mu*diff(v,y) - (2.0/3.0)*mu*(diff(u,x) + diff(v,y)) 
tau_yx = mu*diff(v,x)  #mu*(diff(u,y) + diff(v,x)) 

# Conservative form (?)
#rho*(u*u_x+v*u_y + beta*(u*u_x + u*v_y))
Su = diff(rho*u*u, x) + diff(rho*u*v, y) - diff(tau_xx, x) - diff(tau_xy, y) + diff(p,x)  
#Su = diff(rho*u*u, x) + diff(rho*u*v, y) + diff(p,x)  

#rho*(u*v_x+v*v_y + beta*(v*u_x + v*v_y))
Sv = diff(rho*v*u, x) + diff(rho*v*v, y) - diff(tau_yx, x) - diff(tau_yy, y) + diff(p,y)
#Sv = diff(rho*v*u, x) + diff(rho*v*v, y) + diff(p,y)

# Non-conservative form
#Su = rho*u*diff(u,x) + rho*v*diff(u,y) - diff(tau_xx, x) - diff(tau_xy, y) + diff(p,x)  
#Sv = rho*u*diff(v,x) + rho*v*diff(v,y) - diff(tau_yx, x) - diff(tau_yy, y) + diff(p,y)  

Srho = diff(u*rho,x) + diff(v*rho,y)

Sie = rho*u*diff(ie, x) + rho*v*diff(ie, y) + p*diff(u, x) + p*diff(v, y)

print 'from math import sin, cos, tanh, pi'
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
