x = var('x')
y = var('y')
t = var('t')
h = var('h')
g = var('g')
H = var('d0')

eta = -h*sin(2*pi*(x+t)) - h*sin(2*pi*(y+t))
u_x = h*sin(2*pi*(x+t)) 
u_y = h*sin(2*pi*(y+t))

u_x_src = g*diff(eta, x) + diff(u_x, t)
u_y_src = g*diff(eta, y) + diff(u_y, t)
eta_src = diff(eta, t) + diff(H*u_x, x) + diff(H*u_y, y)

print "u_x: ", u_x.subs(h=1)
print "u_y: ", u_y.subs(h=1)
print "eta: ", eta.subs(h=1)

print "u_x src: ", u_x_src.subs(h=1)
print "u_y src: ", u_y_src.subs(h=1)
print "eta_src: ", eta_src.subs(h=1)
