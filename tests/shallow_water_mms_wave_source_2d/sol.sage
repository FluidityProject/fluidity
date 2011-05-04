t = var('t')
x = var('x')
y = var('y')
c = var('c')
k = var('k')
g = var('g') # gravity
H = var('H') # depth
a = 0
b = 1

eta = -c * sin(k*x+t)
u_x = c * sin(k*x+t)
u_y = 0.0

print "u_x: ", u_x
print "u_y: ", u_y
print "eta: ", eta
print "mom source x: ", g*diff(eta, x) + diff(u_x, t)
print "mom source y: ", g*diff(eta, y) + diff(u_y, t)
print "cont source: ", diff(eta, t) + diff(H*u_x, x) + diff(H*u_y, y)
