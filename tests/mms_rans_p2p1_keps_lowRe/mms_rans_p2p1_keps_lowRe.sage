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

vel_scale = 1e1
u = vel_scale*function(0.5, 1.0, 0.3, 0.0, 
                       0.0, 1.0, 0.5, 0.5, 1.0, 0.0,
                       0.6, 0.75, 1.0)
p_prime = function(-1.0, 1.0, 1.0, 1.0,
                    1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                    1.0, 1.0, 1.0)
rho_scale = 1e2
rho = rho_scale*function(5.2, -1.8, -1.3, 3.7, 
                         1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                         1.7, 2.1, 1.3)
ke_scale = 1e0
eps_scale = 1e0
ke = ke_scale*function(1.0, 0.9, 0.6, 0.4, 
                       0.0, 1.0, 1.0, 0.0, 0.0, 1.0,
                       0.6, 0.7, 0.8)
eps = eps_scale*function(13.2, -3.8, 4.3, 0.0, 
                         1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
                         0.7, 0.3, 0.0)
v = integral(-diff(u,x),y)  # divergence free  
nu = 1.0

p = p_prime - 2./3.*ke

g_x = 0.707106781
g_y = 0.707106781

Y = 1.5
R_y = ke**0.5*Y/nu
Re_T = ke**2.0/(eps*nu)

f_mu = (1 - e**-(0.0165*R_y))**2.0*(1.0 + 20.5/Re_T)
f_1 = 1.0 + (0.05/f_mu)**3.0 
f_2 = 1.0 - e**(- Re_T**2.0)

nu_T = f_mu*ke^2/eps

tau_xx = 2*nu*diff(u,x)            
tau_xy = nu*(diff(u,y) + diff(v,x))
tau_yy = 2*nu*diff(v,y)            
tau_yx = nu*(diff(u,y) + diff(v,x))  

tau_xx_R = 2*nu_T*diff(u,x) - (2./3.)*ke
tau_xy_R = nu_T*(diff(u,y) + diff(v,x))
tau_yy_R = 2*nu_T*diff(v,y) - (2./3.)*ke
tau_yx_R = nu_T*(diff(u,y) + diff(v,x))

Su = u*diff(u,x) + v*diff(u,y) - diff(tau_xx, x) - diff(tau_xy, y) - diff(tau_xx_R, x) - diff(tau_xy_R, y) - g_x*rho + diff(p,x)  
Sv = u*diff(v,x) + v*diff(v,y) - diff(tau_yx, x) - diff(tau_yy, y) - diff(tau_yx_R, x) - diff(tau_yy_R, y) - g_y*rho + diff(p,y)  

Srho = u*diff(rho,x) + v*diff(rho,y) - (1.0 + nu_T)*(diff(rho, x, x) + diff(rho, y, y)) - diff(nu_T, x)*diff(rho, x) - diff(nu_T, y)*diff(rho, y)

P = nu_T*(2*(diff(u,x)^2 + diff(v,y)^2 + diff(u,y)*diff(v,x)) + diff(u,y)^2 + diff(v,x)^2) - (2./3.)*ke*(diff(u,x) + diff(v,y))

u_z = g_x*u + g_y*v
u_xy = ((u**2.0 + v**2.0) - u_z**2.0)**0.5
C3 = tanh(u_z/u_xy)
B = -nu_T*(g_x*diff(rho,x) + g_y*diff(rho,y))  
pr = 1.0
ab = 1.0
bo = 1.0

Ske = u*diff(ke,x) + v*diff(ke,y) - (nu + nu_T)*(diff(ke, x, x) + diff(ke, y, y)) - diff(nu_T, x)*diff(ke, x) -  diff(nu_T, y)*diff(ke, y) - pr*P + ab*eps - bo*B
Seps = u*diff(eps,x) + v*diff(eps,y) - (nu + nu_T)*(diff(eps, x, x) + diff(eps, y, y)) - diff(nu_T, x)*diff(eps, x) -  diff(nu_T, y)*diff(eps, y) - f_1*pr*(eps/ke)*P + f_2*ab*(eps^2/ke) - f_1*bo*C3*(eps/ke)*B
  
print 'from math import sin, cos, tanh, pi, e, sqrt'
print ''
print 'def u(X):'
print '    return', str(u).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v(X):'
print '    return', str(v).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def p(X):'
print '    return', str(p_prime).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def rho(X):'
print '    return', str(rho).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def ke(X):'
print '    return', str(ke).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def eps(X):'
print '    return', str(eps).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u(X):'
print '    return', str(Su).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v(X):'
print '    return', str(Sv).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_rho(X):'
print '    return', str(Srho).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_ke(X):'
print '    return', str(Ske).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_eps(X):'
print '    return', str(Seps).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P_ke(X):'
print '    return', str(P).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P_eps(X):'
print '    return', str(f_1*(eps/ke)*P).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def A_ke(X):'
print '    return', str(-eps).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def A_eps(X):'
print '    return', str(-eps^2.0/ke).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B_ke(X):'
print '    return', str(B).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B_eps(X):'
print '    return', str(f_1*C3*(eps/ke)*B).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def EV(X):'
print '    return', str(nu_T).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity(X):'
print '    return [u(X), v(X)]'
print ''
print 'def forcing_velocity(X):'
print '    return [forcing_u(X), forcing_v(X)]'
print ''
print 'def A_ke(X):'
print '    return', str(-eps).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def A_eps(X):'
print '    return', str(-f_2*eps**2/ke).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def d_eps(X):'
print '    return', str(diff(eps,y)).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def f_1(X):'
print '    return', str(f_1).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def f_2(X):'
print '    return', str(f_2).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def f_mu(X):'
print '    return', str(f_mu).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def Y(X):'
print '    return', str(Y).replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
