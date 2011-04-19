t = var('t')
h = var('h')
c = var('c')

# --------------------------------------------------

u = h * sin(2*pi*x + c*t)
eta = h * cos(2*pi*x + c*t)

u_src = diff(u, t) + diff(eta, x)
eta_src = diff(eta, t) + diff(u, x)

print "u_src: ", u_src
print "eta_src: ", eta_src

J = integrate((eta.subs(t=1))**2, x, 0, 1)
print "J: ", J
print "dJ/dh: ", diff(J, h)
print "dJ/dc: ", diff(J, c)

# --------------------------------------------------

mu = sin(2*pi*x) + cos(2*pi*t)
lmbda = cos(2*pi*x) + sin(2*pi*t)

j0 = -diff(mu, t) - diff(lmbda, x)
j1 = -diff(lmbda, t) - diff(mu, x)

print "Velocity adjoint source term: ", j0
print "LayerThickness adjoint source term: ", j1
