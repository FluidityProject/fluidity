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

ie1 = 0.5*(cos(x + y) + 1.5)
ie2 = 2.5*(sin(x*x) + 1.5)
#ie2 = 3.0*y

#u1 = 0.25*(sin(2*x+(x**2)*(y**2))+0.5)
#v1 = 0.1*(cos(x**2+y**2)+0.5)

#u2 = 1.0*(sin(x**2+y**2)+0.25)
#v2 = 1.15*(cos(x*y)+0.75)

u1 = 1.0*(sin(x**2+y**2)+0.5)
v1 = 0.1*(cos(x**2+y**2)+0.5)
u2 = cos(x**2+y**2)+2.5*x
v2 = 0.5*x*y 

rho1 = 0.5*(sin(x*x + y*y) + 1.5)
rho2 = 3.0

ie1 = 1.25*x*y + cos(x + y)
ie2 = sin(x*y)


gamma = 1.4
rho0 = 0.5
csq = 100.0
#p = csq*(rho1 - rho0)
p = (gamma-1.0)*ie1*rho1
             
#vfrac1 = 0.5*cos(x)*cos(x) + 0.1
vfrac1 = 0.8
vfrac2 = 1.0 - vfrac1 
            
nu1 = 0.7
nu2 = 0.7

g_x = 0.707106781
g_y = 0.707106781

tau_xx1 = nu1*diff(u1,x)            
tau_yy1 = nu1*diff(v1,y)
tau_xy1 = nu1*diff(u1,y)
tau_yx1 = nu1*diff(v1,x)  

tau_xx2 = nu2*diff(u2,x)            
tau_yy2 = nu2*diff(v2,y)
tau_xy2 = nu2*diff(u2,y)
tau_yx2 = nu2*diff(v2,x)  

Su1 = vfrac1*rho1*u1*diff(u1,x) + vfrac1*rho1*v1*diff(u1,y) - diff(vfrac1*tau_xx1, x) - diff(vfrac1*tau_xy1, y) - vfrac1*g_x*rho1 + vfrac1*diff(p,x)  
Sv1 = vfrac1*rho1*u1*diff(v1,x) + vfrac1*rho1*v1*diff(v1,y) - diff(vfrac1*tau_yy1, y) - diff(vfrac1*tau_yx1, x) - vfrac1*g_y*rho1 + vfrac1*diff(p,y)  

Su2 = vfrac2*rho2*u2*diff(u2,x) + vfrac2*rho2*v2*diff(u2,y) - diff(vfrac2*tau_xx2, x) - diff(vfrac2*tau_xy2, y) - vfrac2*g_x*rho2 + vfrac2*diff(p,x)  
Sv2 = vfrac2*rho2*u2*diff(v2,x) + vfrac2*rho2*v2*diff(v2,y) - diff(vfrac2*tau_yy2, y) - diff(vfrac2*tau_yx2, x) - vfrac2*g_y*rho2 + vfrac2*diff(p,y) 

Srho1 = diff(vfrac1*u1*rho1,x) + diff(vfrac1*v1*rho1,y) + rho1*diff(u2*vfrac2, x) + rho1*diff(v2*vfrac2, y)

# Heat transfer term for the internal energy equations
k = 0.5
Cv1 = 500.0
Cv2 = 500.0
d = 0.1

Pr = Cv1*gamma*nu1/k
Re = rho1*d*(sqrt((u1 - u2)**2 + (v1 - v2)**2))/nu1
Nu = (7.0 - 10.0*vfrac1 + 5.0*vfrac1**2)*(1.0 + 0.7*(Re**0.2)*(Pr**(1.0/3.0))) + (1.33 - 2.4*vfrac1 + 1.2*vfrac1**2)*(Re**0.7)*(Pr**(1.0/3.0))
Q = (6.0*k*vfrac2*Nu/(d**2))*(ie2/Cv2 - ie1/Cv1)

heat_flux_x1 = (k/Cv1)*vfrac1*diff(ie1,x)
heat_flux_y1 = (k/Cv1)*vfrac1*diff(ie1,y)
heat_flux_x2 = (k/Cv2)*vfrac2*diff(ie2,x)
heat_flux_y2 = (k/Cv2)*vfrac2*diff(ie2,y)

Sie1 = vfrac1*rho1*u1*diff(ie1, x) + vfrac1*rho1*v1*diff(ie1, y) - diff(heat_flux_x1, x) - diff(heat_flux_y1, y) + vfrac1*p*diff(u1, x) + vfrac1*p*diff(v1, y) #- Q
Sie2 = vfrac2*rho2*u2*diff(ie2, x) + vfrac2*rho2*v2*diff(ie2, y) - diff(heat_flux_x2, x) - diff(heat_flux_y2, y) #+ Q

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
print 'def ie1(X):'
print '    return', str(ie1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
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
print 'def forcing_ie1(X):'
print '    return', str(Sie1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
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
print 'def ie2(X):'
print '    return', str(ie2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
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
print 'def forcing_ie2(X):'
print '    return', str(Sie2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity2(X):'
print '   return [u2(X), v2(X)]'
print ''
print 'def forcing_velocity2(X):'
print '   return [forcing_u2(X), forcing_v2(X)]'
print ''
