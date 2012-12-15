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
    
# NOTE: Phase 1 is the Fluid phase. Phase 2 is the Particle phase.

u1 = function(3.0, 1.0, 0.6, 0.0, 
             0.0, 1.0, 1.0, 0.0, 1.0, 0.0,
             1.0, 1.0, 1.0)
v1 = function(5.2, -1.8, -1.3, 3.7, 
               1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
               1.7, 2.1, 1.3)
u2 = function(0.9, 0.9, 0.6, 0.4, 
             0.0, 1.0, 1.0, 0.0, 0.0, 1.0,
             0.6, 0.7, 0.8)
v2 = function(-1.0, -2.0, 1.0, 1.1,
             1.0, 0.0, 0.0, 1.0, 0.25, 0.0,
             7.5, 1.0, 3.0)
             
vfrac2 = 0.5*(sin(x*x + y*y)) + 0.1 # The prognostic PhaseVolumeFraction should always be in the Particle phase (phase 2)
vfrac1 = 1.0 - vfrac2

rho1 = 0.5*(sin(x*x + y*y) + 1.5)
rho2 = 3.0

gamma = 1.4
rho0 = 0.5
csq = 100.0
p = csq*(rho1 - rho0)
            
mu1 = 0.7
mu2 = 0.9

g_x = 0.707106781
g_y = 0.707106781

tau_xx1 = 2*mu1*diff(u1,x) - (2.0/3.0)*mu1*(diff(u1,x) + diff(v1,y))
tau_xy1 = mu1*(diff(u1,y) + diff(v1,x))
tau_yy1 = 2*mu1*diff(v1,y) - (2.0/3.0)*mu1*(diff(u1,x) + diff(v1,y)) 
tau_yx1 = mu1*(diff(u1,y) + diff(v1,x))

tau_xx2 = 2*mu2*diff(u2,x) - (2.0/3.0)*mu2*(diff(u2,x) + diff(v2,y))
tau_xy2 = mu2*(diff(u2,y) + diff(v2,x))
tau_yy2 = 2*mu2*diff(v2,y) - (2.0/3.0)*mu2*(diff(u2,x) + diff(v2,y)) 
tau_yx2 = mu2*(diff(u2,y) + diff(v2,x))

Su1 = vfrac1*rho1*u1*diff(u1,x) + vfrac1*rho1*v1*diff(u1,y) - diff(vfrac1*tau_xx1, x) - diff(vfrac1*tau_xy1, y) - vfrac1*rho1*g_x + vfrac1*diff(p,x)  
Sv1 = vfrac1*rho1*u1*diff(v1,x) + vfrac1*rho1*v1*diff(v1,y) - diff(vfrac1*tau_yy1, y) - diff(vfrac1*tau_yx1, x) - vfrac1*rho1*g_y + vfrac1*diff(p,y)  

Su2 = vfrac2*rho2*u2*diff(u2,x) + vfrac2*rho2*v2*diff(u2,y) - diff(vfrac2*tau_xx2, x) - diff(vfrac2*tau_xy2, y) - vfrac2*rho2*g_x + vfrac2*diff(p,x)  
Sv2 = vfrac2*rho2*u2*diff(v2,x) + vfrac2*rho2*v2*diff(v2,y) - diff(vfrac2*tau_yy2, y) - diff(vfrac2*tau_yx2, x) - vfrac2*rho2*g_y + vfrac2*diff(p,y) 

Srho1 = diff(vfrac1*u1*rho1,x) + diff(vfrac1*v1*rho1,y) + rho1*diff(u2*vfrac2, x) + rho1*diff(v2*vfrac2, y)

Svfrac2 = diff(u2*vfrac2,x) + diff(v2*vfrac2,y)

print 'from math import sin, cos, tanh, pi, sqrt'
print ''
print 'def u1(X):'
print '    return', str(u1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v1(X):'
print '    return', str(v1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def p(X):'
print '    return', str(p).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def rho1(X):'
print '    return', str(rho1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''   
print 'def vfrac1(X):'
print '    return', str(vfrac1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def forcing_u1(X):'
print '    return', str(Su1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v1(X):'
print '    return', str(Sv1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_rho1(X):'
print '    return', str(Srho1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity1(X):'
print '   return [u1(X), v1(X)]'
print ''
print 'def forcing_velocity1(X):'
print '   return [forcing_u1(X), forcing_v1(X)]'
print ''
print 'def u2(X):'
print '    return', str(u2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v2(X):'
print '    return', str(v2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def rho2(X):'
print '    return', str(rho2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def vfrac2(X):'
print '    return', str(vfrac2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def forcing_u2(X):'
print '    return', str(Su2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v2(X):'
print '    return', str(Sv2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity2(X):'
print '   return [u2(X), v2(X)]'
print ''
print 'def forcing_velocity2(X):'
print '   return [forcing_u2(X), forcing_v2(X)]'
print ''
print 'def forcing_vfrac2(X):'
print '    return', str(Svfrac2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
