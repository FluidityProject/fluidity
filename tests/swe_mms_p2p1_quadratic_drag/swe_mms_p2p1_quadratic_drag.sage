y = var('y')
from math import pi

h = sin(x)*sin(y)
h_mean = 20.0
u = cos(y)*sin(x)
v = integral(-diff(u,x),y)

print diff((h_mean+h)*u, x) + diff((h_mean+h)*v, y)

g = 9.8

nu = 0.6
tau_xx = 2*nu*diff(u,x) - (2.0/3.0)*nu*(diff(u,x) + diff(v,y))
tau_xy = nu*(diff(u,y) + diff(v,x))
tau_yy = 2*nu*diff(v,y) - (2.0/3.0)*nu*(diff(u,x) + diff(v,y))
tau_yx = nu*(diff(u,y) + diff(v,x))

C_D = 0.0025 # Drag coefficient
magnitude = sqrt(u*u + v*v)

Su = u*diff(u,x) + v*diff(u,y) + g*diff(h,x) - diff(tau_xx, x) - diff(tau_xy, y) + (C_D/(h_mean + h))*magnitude*u
Sv = u*diff(v,x) + v*diff(v,y) + g*diff(h,y) - diff(tau_yy, y) - diff(tau_yx, x) + (C_D/(h_mean + h))*magnitude*v

print 'from math import sin, cos, tanh, pi, sqrt'
print ''
print 'def u(X):'
print '    return', str(u).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''
print 'def v(X):'
print '    return', str(v).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''  
print 'def h(X):'
print '    return', str(h).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print '' 
print 'def forcing_u(X):'
print '    return', str(Su).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''
print 'def forcing_v(X):'
print '    return', str(Sv).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''
print 'def velocity(X):'
print '   return [u(X), v(X)]'
print ''
print 'def forcing_velocity(X):'
print '   return [forcing_u(X), forcing_v(X)]'
print ''
