t = var('t')
h = var('h')

# --------------------------------------------------

u = h*t*(10-t)*sin(2*pi*x)
eta = h*t*(5-t)*cos(2*pi*x)

u_src = diff(u, t) + diff(eta, x)
eta_src = diff(eta, t) + diff(u, x)

print "u_src: ", u_src
print "eta_src: ", eta_src

J = integrate((eta.subs(t=1))**2, x, 0, 1)
print "J: ", J.subs(h=1)
print "diff(J, h).subs(h=1): ", diff(J, h).subs(h=1)

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
