y = var('y')
u = 0.1*sin(x)*cos(y)
v = -0.1*cos(x)*sin(y)
p = 1.0-cos(x)*cos(y)
rho_p = 0.1*cos(x)*cos(y)

Su = u*diff(u,x) + v*diff(u,y) - (diff(u, x, x) + diff(u, y, y)) + diff(p,x)
Sy = u*diff(v,x) + v*diff(v,y) - (diff(v, x, x) + diff(v, y, y)) + diff(p,y) + rho_p*-1.0 

Srho_p = u*diff(rho_p,x) + v*diff(rho_p,y)

print str(Su).replace('e^', 'exp').replace('^', '**').replace('000000000000', '')
print str(Sy).replace('e^', 'exp').replace('^', '**').replace('000000000000', '')
print str(Srho_p).replace('e^', 'exp').replace('^', '**').replace('000000000000', '')
