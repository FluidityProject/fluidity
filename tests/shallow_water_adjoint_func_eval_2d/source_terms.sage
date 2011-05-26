t = var('t')
x = var('x')
y = var('y')
g = var('g')
H = var('d0')

eta = -sin(2*pi*(x+t)) -sin(2*pi*(y+t))
u_x = g*sin(2*pi*(x+t))
u_y = g*sin(2*pi*(y+t))

u_x_src = g*diff(eta, x) + diff(u_x, t)
u_y_src = g*diff(eta, y) + diff(u_y, t)
eta_src = diff(eta, t) + diff(H*u_x, x) + diff(H*u_y, y)

print "u_x: ", u_x
print "u_y: ", u_y
print "eta: ", eta
print "u_x src: ", u_x_src
print "u_y src: ", u_y_src
print "eta_src: ", eta_src
