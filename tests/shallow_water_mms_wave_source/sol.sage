t = var('t')
x = var('x')
c = var('c')
k = var('k')
a = 0
b = 1

eta = -c * sin(k*x+t)
u = c * sin(k*x+t)

print "u: ", u
print "eta: ", eta
print "mom source: ", diff(eta, x) + diff(u, t)
print "cont source: ", diff(eta, t) + diff(u, x)
