import numpy

t = var('t')
T = 10.0

u = sin(x + t)
v = 1.0

f = diff(u, t) + u * diff(u, x) - v * diff(u, x, 2) # forward source
J = integrate(integrate(u, x, -10, 10), t, 0, T)

dt = 2.0
ntimesteps = int(numpy.ceil(10.0/dt))
partial_J = numpy.zeros(ntimesteps)
for i in range(ntimesteps):
  partial_J[i] = float(integrate(integrate(u, x, -10, 10), t, i*dt, (i+1)*dt))


print "Forward source term: ", str(f).replace('e^', 'exp').replace('^', '**')
print "Partial values of functional: ", partial_J
print "Sum of partial values: ", sum(partial_J)
print "True value of functional: ", float(J)

intu = integrate(u, x, -10, 10)
print "integrate(u, x, -10, 10): ", intu
for i in range(ntimesteps+1):
  print " .. (t=%f): %f" % (dt * i, intu.subs(t=dt * i))
