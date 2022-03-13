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

p = function(-1.0, 1.0, 1.0, 1.0,
             1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
             1.0, 1.0, 1.0)
# ensure the pressure solution has 0 mean:
pavg = N(integrate(integrate(p, x, 0, pi), y, 0, pi)/(pi^2))
p -= pavg

rho1 = 2.5
rho2 = 0.5
u1 = 0.25*cos(x)*cos(y) - x*cos(y)
v1 = sin(y)
u2 = sin(x)*cos(y)
v2 = sin(y)*sin(x) - cos(x)*sin(y)

vfrac1 = 0.8 #0.1 + y/20.0
vfrac2 = 1.0 - vfrac1 #1.0 - vfrac1

#print "DIVERGENCE = ", diff(vfrac1*u1,x) + diff(vfrac1*v1,y) + diff(vfrac2*u2,x) + diff(vfrac2*v2,y)

nu = 0.7

tau_xx1 = 2*nu*diff(u1,x) - (2.0/3.0)*nu*(diff(u1,x) + diff(v1,y))
tau_xy1 = nu*(diff(u1,y) + diff(v1,x))
tau_yy1 = 2*nu*diff(v1,y) - (2.0/3.0)*nu*(diff(u1,x) + diff(v1,y))
tau_yx1 = nu*(diff(u1,y) + diff(v1,x))

tau_xx2 = 2*nu*diff(u2,x) - (2.0/3.0)*nu*(diff(u2,x) + diff(v2,y))
tau_xy2 = nu*(diff(u2,y) + diff(v2,x))
tau_yy2 = 2*nu*diff(v2,y) - (2.0/3.0)*nu*(diff(u2,x) + diff(v2,y))
tau_yx2 = nu*(diff(u2,y) + diff(v2,x))

Su1 = vfrac1*rho1*u1*diff(u1,x) + vfrac1*rho1*v1*diff(u1,y) - diff(vfrac1*tau_xx1, x) - diff(vfrac1*tau_xy1, y) + vfrac1*diff(p,x)
Sv1 = vfrac1*rho1*u1*diff(v1,x) + vfrac1*rho1*v1*diff(v1,y) - diff(vfrac1*tau_yx1, x) - diff(vfrac1*tau_yy1, y) + vfrac1*diff(p,y)
Su2 = vfrac2*rho2*u2*diff(u2,x) + vfrac2*rho2*v2*diff(u2,y) - diff(vfrac2*tau_xx2, x) - diff(vfrac2*tau_xy2, y) + vfrac2*diff(p,x)
Sv2 = vfrac2*rho2*u2*diff(v2,x) + vfrac2*rho2*v2*diff(v2,y) - diff(vfrac2*tau_yx2, x) - diff(vfrac2*tau_yy2, y) + vfrac2*diff(p,y)

print 'from math import sin, cos, tanh, pi'
print ''
print 'def u1(X):'
print '    return', str(u1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v1(X):'
print '    return', str(v1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def u2(X):'
print '    return', str(u2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v2(X):'
print '    return', str(v2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def vfrac1(X):'
print '    return', str(vfrac1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def vfrac2(X):'
print '    return', str(vfrac2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def p(X):'
print '    return', str(p).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def rho1(X):'
print '    return', str(rho1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def rho2(X):'
print '    return', str(rho2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u1(X):'
print '    return', str(Su1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v1(X):'
print '    return', str(Sv1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u2(X):'
print '    return', str(Su2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v2(X):'
print '    return', str(Sv2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity1(X):'
print '   return [u1(X), v1(X)]'
print ''
print 'def velocity2(X):'
print '   return [u2(X), v2(X)]'
print ''
print 'def forcing_velocity1(X):'
print '   return [forcing_u1(X), forcing_v1(X)]'
print ''
print 'def forcing_velocity2(X):'
print '   return [forcing_u2(X), forcing_v2(X)]'
print ''
