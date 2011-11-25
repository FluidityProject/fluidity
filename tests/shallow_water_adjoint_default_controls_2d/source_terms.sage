t = var('t')
x = var('x')
y = var('y')
g = var('g')
h = var('h')
H = var('d0')

eta = (t+1)*h*cos(2*pi*x)
u_x = (t+1)*h*sin(2*pi*x)
u_y = cos(2*pi*x)+cos(2*pi*y)

u_x_src = g*diff(eta, x) + diff(u_x, t)
u_y_src = g*diff(eta, y) + diff(u_y, t)
eta_src = diff(eta, t) + diff(H*u_x, x) + diff(H*u_y, y)

print "u_x: ", u_x
print "u_y: ", u_y
print "eta: ", eta

print "u_x src: ", u_x_src
print "u_y src: ", u_y_src
print "eta_src: ", eta_src

J = integrate(integrate((eta.subs(t=1))**2, x, 0, 1), y, 0, 1)
print "J(t=1): ", J.subs(h=1)
print "diff(J(t=1), h).subs(h=1): ", diff(J, h).subs(h=1)
