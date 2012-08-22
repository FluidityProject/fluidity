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

u = function(0.5, 1.0, 0.3, 0.0, 
             0.0, 1.0, 0.5, 0.5, 1.0, 0.0,
             0.6, 0.75, 1.0)
p = function(-1.0, 1.0, 1.0, 1.0,
              1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
              1.0, 1.0, 1.0)
s1 = function(0.1, -0.01, -0.02, 0.07, 
              1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
              1.7, 2.1, 1.3)
s2 = function(0.1, 0.01, -0.01, -0.1, 
              1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
              1.4, 3.0, 0.6)
v = integral(-diff(u,x),y)  # divergence free  

nu = 1.0/(1.0 - (s1 + s2)/0.65)**1.625 # concentration dependent viscosity

s1_R = 0.33
s2_R = 0.66
rho = s1_R*s1 + s2_R*s2 

tau_xx = 2*nu*diff(u,x)            
tau_xy = nu*(diff(u,y) + diff(v,x))
tau_yy = 2*nu*diff(v,y)            
tau_yx = nu*(diff(u,y) + diff(v,x))  

s1_us = 0.33*(1.0 - (s1 + s2))**2.39 # hindered settling
s2_us = 0.66*(1.0 - (s1 + s2))**2.39 # hindered settling

Su = u*diff(u,x) + v*diff(u,y) - diff(tau_xx, x) - diff(tau_xy, y) - rho + diff(p,x)  
Sv = u*diff(v,x) + v*diff(v,y) - diff(tau_yx, x) - diff(tau_yy, y) - rho + diff(p,y)  

Ss1 = (u + s1_us)*diff(s1,x) + (v + s1_us)*diff(s1,y) - 1.0*(diff(s1, x, x) + diff(s1, y, y)) 
Ss2 = (u + s2_us)*diff(s2,x) + (v + s2_us)*diff(s2,y) - 1.0*(diff(s2, x, x) + diff(s2, y, y)) 
  
print 'from math import sin, cos, tanh, pi, e, sqrt'
print ''
print 'def u(X):'
print '    return', str(u).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v(X):'
print '    return', str(v).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def p(X):'
print '    return', str(p).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def s1(X):'
print '    return', str(s1).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def s2(X):'
print '    return', str(s2).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u(X):'
print '    return', str(Su).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v(X):'
print '    return', str(Sv).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_s1(X):'
print '    return', str(Ss1).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_s2(X):'
print '    return', str(Ss2).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def rho(X):'
print '    return', str(rho).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def nu(X):'
print '    return', str(nu).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity(X):'
print '    return [u(X), v(X)]'
print ''
print 'def forcing_velocity(X):'
print '    return [forcing_u(X), forcing_v(X)]'
