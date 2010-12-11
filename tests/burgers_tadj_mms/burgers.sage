#!/usr/bin/env sage
t = var('t')
m = var('m')

a = -10
b = +10
v = 1.0

dt = 1.0
theta = var('theta')

lambda1 = (x + 10) * x * (x - 10)
src = simplify(expand(lambda1 / dt - theta * v * diff(lambda1, x, 2)))

print "Source term at time t=1.0: ", src
lambda0 = simplify(expand(lambda1 / dt - (theta - 1.0) * v * diff(lambda1, x, 2)))
print "Adjoint solution at time t=0.0: ", lambda0
