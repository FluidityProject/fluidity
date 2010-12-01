#!/usr/bin/env sage
t = var('t')
m = var('m')

a = -10
b = +10
v = 1.0

u = m * (x - a) * (x - b) # velocity, parameterised by a parameter m
T = 0.01 * (x - a) * (x - b) # target velocity is achieved when m = 1

# first, let's compute the source term for the forward model:
src = simplify(expand(diff(u, t) + u * diff(u, x) - v * diff(u, x, 2)))

print "Source term for forward model: ", src
print "\partial F/ \partial m: ", -1 * diff(src, m)

# now compute J(m)
# J(u) = 0.5 * \int_-10^+10 (u - t)**2 dx
J = (1/2) * integrate( (u - T)**2, x, a, b)

print "J as a functional of m: ", J
print "dJ/dm: ", diff(J, m)

print "J(m=0.005): ", J.subs(m=0.005)
print "dJ/dm|m=0.005: ", diff(J, m).subs(m=0.005)
