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

u = 1.0*(sin(x**2+y**2)+0.5)
v = 0.1*(cos(x**2+y**2)+0.5) 

rho = 0.5*(sin(x*x + y*y) + 1.5)

ie = 1.25*x*y + cos(x + y) + 20.0

gamma = 1.4
csq = 100.0
rho0 = 0.5
p = (gamma-1.0)*ie*rho
#p = csq*(rho - rho0)

nu = 0.7

g_x = 0.707106781
g_y = 0.707106781

tau_xx = nu*diff(u,x)            
tau_yy = nu*diff(v,y)
tau_xy = nu*diff(u,y)
tau_yx = nu*diff(v,x)  

Su = rho*u*diff(u,x) + rho*v*diff(u,y) - diff(tau_xx, x) - diff(tau_xy, y) - g_x*rho + diff(p,x)  
Sv = rho*u*diff(v,x) + rho*v*diff(v,y) - diff(tau_yy, y) - diff(tau_yx, x) - g_y*rho + diff(p,y)  

Srho = diff(u*rho,x) + diff(v*rho,y)

Sie = rho*u*diff(ie, x) + rho*v*diff(ie, y) + p*diff(u, x) + p*diff(v, y)
# - diff(0.5*diff(ie,x), x) - diff(0.5*diff(ie,y), y)

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
