t = var('t')

u = 1/sqrt(2*pi) * exp(-x**2 / 2)
v = 1.0


f = diff(u, t) + u * diff(u, x) - v * diff(u, [x, x])
s = str(f).replace('e^', 'exp').replace('^', '**')
print s

endpoints = (-10.0, 10.0)
du = diff(u, x)
bcs = [float(sign(y)*du(x=y)) for y in endpoints]
print bcs
