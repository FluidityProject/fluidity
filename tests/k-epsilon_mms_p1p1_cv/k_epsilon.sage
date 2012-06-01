y = var('y')
u = 0.1*sin(x)*cos(y)
v = -0.1*cos(x)*sin(y)
p = cos(x)*cos(y)-1.0
rho_p = 0.1*cos(x)*cos(y)
ke = y+1.0
eps = y+1.0

nuT = ke^2/eps

Su = u*diff(u,x) + v*diff(u,y) + diff(p,x) - (1.0+nuT)*(diff(u, x, x) + diff(u, y, y)) - diff(u, x)*diff(nuT, x) - diff(u, y)*diff(nuT, y)
Sy = u*diff(v,x) + v*diff(v,y) + diff(p,y) - (1.0+nuT)*(diff(v, x, x) + diff(v, y, y)) - diff(v, x)*diff(nuT, x) - diff(v, y)*diff(nuT, y) - rho_p*-1.0 

Srho_p = u*diff(rho_p,x) + v*diff(rho_p,y) - nuT*(diff(rho_p, x, x) + diff(rho_p, y, y)) - diff(nuT, x)*diff(rho_p, x) -  diff(nuT, y)*diff(rho_p, y)

S = 2*(diff(u,x)**2 + diff(v,y)**2 + diff(u,y)*diff(v,x)) + diff(u,y)**2 + diff(v,x)**2
Pke = nuT*S
Gke = nuT*diff(rho_p,y)
Ake = eps
Ske = u*diff(ke,x) + v*diff(ke,y) - nuT*(diff(ke, x, x) + diff(ke, y, y)) - diff(nuT, x)*diff(ke, x) -  diff(nuT, y)*diff(ke, y) - Pke - Gke + Ake

Peps = eps/ke*nuT*S
Geps = eps/ke*nuT*diff(rho_p,y)
Aeps = eps^2/ke
Ce3 = tanh(v/u)
Seps = u*diff(eps,x) + v*diff(eps,y) - nuT*(diff(eps, x, x) + diff(eps, y, y)) - diff(nuT, x)*diff(eps, x) -  diff(nuT, y)*diff(eps, y) - Peps - Ce3*Geps + Aeps

print str(Su).replace('e^', 'exp').replace('^', '**')
print str(Sy).replace('e^', 'exp').replace('^', '**')
print str(Srho_p).replace('e^', 'exp').replace('^', '**')
print str(Ske).replace('e^', 'exp').replace('^', '**')
print str(Seps).replace('e^', 'exp').replace('^', '**')

print str(Pke).replace('e^', 'exp').replace('^', '**')
print str(Gke).replace('e^', 'exp').replace('^', '**')
print str(Ake).replace('e^', 'exp').replace('^', '**')