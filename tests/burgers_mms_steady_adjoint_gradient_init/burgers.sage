import numpy

u = sin(x) + cos(x)
v = 1.0

#f = u * diff(u, x) - v * diff(u, x, 2) # forward source
f = - v * diff(u, x, 2) # forward source
J = 0.5 * integrate(u**2, x, -10, 10)

print "Forward source term: ", str(f).replace('e^', 'exp').replace('^', '**')
print "True value of functional: ", float(J)
