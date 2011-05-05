t = var('t')
x = var('x')
y = var('y')
k = var('k')
g = var('g') # gravity
H = var('H') # depth
a = 0
b = 1

eta = -sin(2*pi*x+t)
u_x = sin(2*pi*x+t)
u_y = cos(2*pi*y+t)

print "u_x: ", u_x
print "u_y: ", u_y
print "eta: ", eta
print "u_x source: ", g*diff(eta, x) + diff(u_x, t)
print "u_y source: ", g*diff(eta, y) + diff(u_y, t)
print "cont source: ", diff(eta, t) + diff(H*u_x, x) + diff(H*u_y, y)
