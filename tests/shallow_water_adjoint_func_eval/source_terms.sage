t = var('t')
g = var('g')
H = var('H')

# --------------------------------------------------

u = sin(2*pi*x) + t
eta = sin(2*pi*x)

u_src = diff(u, t) + g*diff(eta, x)
eta_src = diff(eta, t) + diff(H*u, x)

print "u_src: ", u_src
print "eta_src: ", eta_src
