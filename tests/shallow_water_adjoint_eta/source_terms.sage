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
print "--------------------------------------------------"

mu = 0
lmbda = t*(10-t)*(5-t)

j0 = diff(mu, t) - diff(lmbda, x)
j1 = diff(lmbda, t) - diff(mu, x)

print "Velocity adjoint source term: ", j0
print "LayerThickness adjoint source term: ", j1

# --------------------------------------------------
print "--------------------------------------------------"

mu = t*(10-t)*(4-t)
lmbda = 0

j0 = diff(mu, t) - diff(lmbda, x)
j1 = diff(lmbda, t) - diff(mu, x)

print "Velocity adjoint source term: ", j0
print "LayerThickness adjoint source term: ", j1
