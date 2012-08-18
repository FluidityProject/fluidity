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

u = function(2.5, 1.0, 0.6, 0.0, 
             0.0, 1.0, 1.0, 0.0, 1.0, 0.0,
             1.0, 1.0, 1.0)
p = function(-1.0, 1.0, 1.0, 1.0,
             1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
             1.0, 1.0, 1.0)
temperature = function(5.2, -1.8, -1.3, 3.7, 
               1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
               1.7, 2.1, 1.3)
ke = function(0.9, 0.9, 0.6, 0.4, 
             0.0, 1.0, 1.0, 0.0, 0.0, 1.0,
             0.6, 0.7, 0.8)
eps = function(8.2, -3.8, 4.3, 1.7, 
             1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
             0.7, 0.8, 0.6)
v = integral(-diff(u,x),y)  # divergence free

rho = 10.0*temperature + 10.0
nu_T = rho*(ke^2)/eps
nu = 0.7

tau_xx = 2*nu*diff(u,x)            
tau_xy = nu*(diff(u,y) + diff(v,x))
tau_yy = 2*nu*diff(v,y)            
tau_yx = nu*(diff(u,y) + diff(v,x))  

tau_xx_R = 2*nu_T*diff(u,x) - (2./3.)*ke*rho
tau_xy_R = nu_T*(diff(u,y) + diff(v,x))
tau_yy_R = 2*nu_T*diff(v,y) - (2./3.)*ke*rho
tau_yx_R = nu_T*(diff(u,y) + diff(v,x))

Su = rho*u*diff(u,x) + rho*v*diff(u,y) - diff(tau_xx, x) - diff(tau_xy, y) - diff(tau_xx_R, x) - diff(tau_xy_R, y) - rho + diff(p,x)  
Sv = rho*u*diff(v,x) + rho*v*diff(v,y) - diff(tau_yx, x) - diff(tau_yy, y) - diff(tau_yx_R, x) - diff(tau_yy_R, y) - rho + diff(p,y)  

Stemperature = u*diff(temperature,x) + v*diff(temperature,y) - (1.0 + nu_T)*(diff(temperature, x, x) + diff(temperature, y, y)) - diff(nu_T, x)*diff(temperature, x) -  diff(nu_T, y)*diff(temperature, y)

P = nu_T*(2*(diff(u,x)^2 + diff(v,y)^2 + diff(u,y)*diff(v,x)) + diff(u,y)^2 + diff(v,x)^2) - (2./3.)*ke*rho*(diff(u,x) + diff(v,y))

C3 = 1.0
g_x = 1.0
g_y = 1.0
B = nu_T*(g_x*diff(rho,x) + g_y*diff(rho,y))/rho

pr = 1
ab = 1
bo = 1

Ske = rho*u*diff(ke,x) + rho*v*diff(ke,y) - nu_T*(diff(ke, x, x) + diff(ke, y, y)) - diff(nu_T, x)*diff(ke, x) -  diff(nu_T, y)*diff(ke, y) - pr*P + ab*rho*eps - bo*B
Seps = rho*u*diff(eps,x) + rho*v*diff(eps,y) - nu_T*(diff(eps, x, x) + diff(eps, y, y)) - diff(nu_T, x)*diff(eps, x) -  diff(nu_T, y)*diff(eps, y) - pr*(eps/ke)*P + ab*rho*(eps^2/ke) - bo*C3*(eps/ke)*B
  
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
print 'def temperature(X):'
print '    return', str(temperature).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def rho(X):'
print '    return', str(rho).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print '' 
print 'def ke(X):'
print '    return', str(ke).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def eps(X):'
print '    return', str(eps).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u(X):'
print '    return', str(Su).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v(X):'
print '    return', str(Sv).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_temperature(X):'
print '    return', str(Stemperature).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_ke(X):'
print '    return', str(Ske).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_eps(X):'
print '    return', str(Seps).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P(X):'
print '    return', str(P).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P_eps(X):'
print '    return', str((eps/ke)*P).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def A(X):'
print '    return', str(rho*eps/ke).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B(X):'
print '    return', str(B).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B_eps(X):'
print '    return', str(C3*(eps/ke)*B).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def EV(X):'
print '    return', str(nu_T).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity(X):'
print '   return [u(X), v(X)]'
print ''
print 'def forcing_velocity(X):'
print '   return [forcing_u(X), forcing_v(X)]'
print ''
print 'def A_ke(X):'
print '   return A([X[0],X[1]])*ke([X[0],X[1]])'
print ''
print 'def A_eps(X):'
print '   return A([X[0],X[1]])*eps([X[0],X[1]])'
