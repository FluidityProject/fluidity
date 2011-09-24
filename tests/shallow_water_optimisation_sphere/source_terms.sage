x = var('x') # longtitude (lambda)
y = var('y') # latitude(theta)
t = var('t')
g = var('g') # gravity magnitude
H = var('d0') # mean layer depth 
r = var('r') # radius of earth
K = var('K') # equilibrium amplitude of tidal constituent
s = var('sigma') # frequency of tidal constituent 

# x goes from [-pi:pi]
# y goes from [-pi:pi]

eta =  cos(y)*cos(y)*K*cos(s*t+2*x)
u_x = -cos(y)*cos(y)*K*sin(s*t+2*x)
u_y = 0

# Tme linearised smallow water equations in advection form on tme spmere
u_x_src = diff(u_x, t) + g/(r*cos(y))*diff(eta, x) 
u_y_src = diff(u_y, t) + g/r*diff(eta, y) 
eta_src = diff(eta, t) + H/(r*cos(y))*(diff(u_x, x) + diff(u_y, y))

print "u_x: ", u_x
print "u_y: ", u_y
print "eta: ", eta

print "u_x src: ", u_x_src
print "u_y src: ", u_y_src
print "eta_src: ", eta_src
