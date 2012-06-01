y = var('y')
u = 0.1*sin(x)*cos(y)
v = -0.1*cos(x)*sin(y)
p = cos(x)*cos(y)-1.0
rho_p = 0.1*cos(x)*cos(y)

Su = u*diff(u,x) + v*diff(u,y) + diff(p,x) - diff(u, x, x) - diff(u, y, y)
Sy = u*diff(v,x) + v*diff(v,y) + diff(p,y) - diff(v, x, x) - diff(v, y, y) - rho_p*-1.0 

Srho_p = u*diff(rho_p,x) + v*diff(rho_p,y)  

str(Su).replace('e^', 'exp').replace('^', '**')
str(Sy).replace('e^', 'exp').replace('^', '**')
str(Srho_p).replace('e^', 'exp').replace('^', '**')

print Su
print Sy
print Srho_p