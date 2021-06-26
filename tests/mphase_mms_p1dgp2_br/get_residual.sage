# This Sage script obtains the residual when 
# the 'analytical' solutions for u and p are
# plugged in to the momentum equation.
# ctj10, March 2012

# Spatial variables x and y
x = var('x')
y = var('y')

# PhaseVolumeFraction. For single-phase, vfrac = 1.0. For multiphase, sum(vfrac_i) = 1.0.
vfrac = var('vfrac')
dvfrac = [var('dvfrac_x'), var('dvfrac_y')]
# Density.
rho = var('rho')
# Isotropic Viscosity.
mu = var('mu')

# Analytical solution
u = [sin(x)*cos(y), -sin(y)*cos(x)]
p = cos(x)*cos(y)

# For MMS tests, the analytical solution needs to be divergence-free.
# Let's check this first...
if(diff(u[0], x) + diff(u[1], y) != 0):
   print "\nWARNING: Analytical solution for this phase is NOT divergence-free!"
   print "If single-phase: choose a divergence-free solution."
   print "If multiphase: this may be ok; check that sum(vfrac_i * u_i) = 0.\n"
else:
   print "Analytical solution for this phase is divergence-free."

# Find the residual of the momentum equation:
# rho*u_t + rho*(u*u_x+v*u_y) + p_x - nu*u_xx - nu*u_yy

# Time independent problem, so the mass term is zero.
mass_x = 0
mass_y = 0

# Advection term in non-conservative form (i.e. beta = 0)
advection_x = rho*vfrac*(u[0]*diff(u[0], x) + u[1]*diff(u[0], y))
advection_y = rho*vfrac*(u[0]*diff(u[1], x) + u[1]*diff(u[1], y))

# Pressure gradient
pressure_gradient_x = vfrac*diff(p, x)
pressure_gradient_y = vfrac*diff(p, y)

# Stress term (assumes isotropic viscosity)
# -mu*diff(vfrac*diff(u)) = vfrac*diff(u, 2) + diff(u)*diff(vfrac)
stress_x = -mu*(vfrac*(diff(diff(u[0], x), x) + diff(diff(u[0], y), y))) -mu*(dvfrac[0]*diff(u[0], x) + dvfrac[1]*diff(u[0], y))
stress_y = -mu*(vfrac*(diff(diff(u[1], x), x) + diff(diff(u[1], y), y))) -mu*(dvfrac[0]*diff(u[1], x) + dvfrac[1]*diff(u[1], y))

momentum_x = mass_x + advection_x + pressure_gradient_x + stress_x
momentum_y = mass_y + advection_y + pressure_gradient_y + stress_y

print momentum_x
print momentum_y
