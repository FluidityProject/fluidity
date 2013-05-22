y = var('y')

def function(phi_0, phi_x, phi_y, phi_xy, 
             f_sin_x, f_cos_x, f_sin_y, f_cos_y, f_sin_xy, f_cos_xy, 
             alpha_x, alpha_y, alpha_xy):
    
    f_0 = phi_0 
    f_x = phi_x*(f_sin_x*sin(alpha_x*x) + f_cos_x*cos(alpha_x*x)) 
    f_y = phi_y*(f_sin_y*sin(alpha_y*y) + f_cos_y*cos(alpha_y*y)) 
    f_xy = phi_xy*(f_sin_xy*sin(alpha_xy*x*y) + f_cos_xy*cos(alpha_xy*x*y)) 
    f = f_0 + f_x + f_y + f_xy
    return f

u = function(0.0, 1.0, 1.0, 0.0, 
             0.0, 1.0, 0.0, 1.0, 0.5, 0.0,
             1.0, 1.0, 2.0)
v = integral(-diff(u,x),y) # divergence free 
nu = function(4.0, 1.2, 1.4, -1.2, 
              1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
              1.7, 1.1, 0.3) 

tau_xx = 2*nu*diff(u,x)            
tau_xy = nu*(diff(u,y) + diff(v,x))
tau_yy = 2*nu*diff(v,y)            
tau_yx = nu*(diff(u,y) + diff(v,x))  

# Helmholtz problem  
Su = - (diff(tau_xx, x) + diff(tau_xy, y)) + u
Sv = - (diff(tau_yx, x) + diff(tau_yy, y)) + v
  
print 'from math import sin, cos, tanh, pi, e, sqrt'
print ''
print 'def u(X):'
print '    return', str(u).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v(X):'
print '    return', str(v).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def nu(X):'
print '    return', str(nu).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u(X):'
print '    return', str(Su).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v(X):'
print '    return', str(Sv).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def U(X):'
print '    return [u(X), v(X)]'
print ''
print 'def forcing_U(X):'
print '    return [forcing_u(X), forcing_v(X)]'
