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

# tracer
c = -y
# molecular diffusion
k = 1.0

# velocity
u = y
v = integral(-diff(u,x),y)  # divergence free

Ri = -diff(c,y)/diff(u, y)^2.0
Ri_c = 3.0

# eddy visc
nu_T = function(2.0, 0.4, 0.3, 0.8, 
                1.0, 0.0, 1.0, 0.0, 0.0, 1.0,
                1.0, 1.0, 1.0)
# prandtl
pr = 0.5

# calculate source term
k_T_x = k + nu_T/pr
k_T_y = k + ((1-Ri/Ri_c)^0.5)*nu_T/pr
D_cx = k_T_x*diff(c,x) 
D_cy = k_T_y*diff(c,y)
Sc = -(diff(D_cx, x) + diff(D_cy, y))
  
print 'from math import sin, cos, tanh, pi, sqrt'
print ''
print 'def nu_T(X):'
print '    return', str(nu_T).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def Ri(X):'
print '    return', str(Ri).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print '' 
print 'def c(X):'
print '    return', str(c).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print '' 
print 'def forcing_c(X):'
print '    return', str(Sc).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def u(X):'
print '    return', str(u).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v(X):'
print '    return', str(v).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def U(X):'
print '   return [u(X), v(X)]'
print ''
