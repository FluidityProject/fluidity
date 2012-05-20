# This Sage script obtains the residual when 
# the 'analytical' solutions for u and p are
# plugged in to the momentum equation.
# ctj10, April 2012

# Spatial variables x and y
x = var('x')
y = var('y')

# PhaseVolumeFraction. For single-phase, vfrac = 1.0. For multiphase, sum(vfrac_i) = 1.0.
vfrac_c = var('vfrac_c')
vfrac_i = var('vfrac_i')
vfrac = var('vfrac')
dvfrac = [var('dvfrac_x'), var('dvfrac_y')]
# Isotropic Viscosity.
mu = var('mu')

# Analytical solutions
u_c = [1.0*(sin(x**2+y**2)+0.5), 0.1*(cos(x**2+y**2)+0.5)]
u_i = [cos(x**2+y**2)+2.5*x, 0.5*x*y]
u = u_i

# Note that this phase is compressible, so we have an analytical solution for density
rho_c = 0.5*(sin(x*x + y*y) + 1.5)
rho_i = 2.0 # Incompressible phase density
rho = rho_i

# Pressure depends on the density of the compressible phase, and is given by the compressible_gas equation of state
csq = var('csq') # Bulk sound speed squared, c^2
gamma = var('gamma') # Ratio of specific heats
rho0 = var('rho0') # Reference density
e = 1.25*x*y + cos(x + y) # Internal energy
p = csq*(rho_c - rho0) + (gamma - 1)*rho_c*e

# For compressible MMS tests, we also provide a source term to the Density field
# so that the continuity equation is satisfied.
print "\nSOURCE TERM FOR COMPRESSIBLE DENSITY FIELD (i.e. the compressible continuity equation's residual):"
print diff(u_c[0]*rho_c*vfrac_c, x) + diff(u_c[1]*rho_c*vfrac_c, y) + rho_c*diff(u_i[0]*vfrac_i, x) + rho_c*diff(u_i[1]*vfrac_i, y)
print "\n"

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

print "SOURCE TERM FOR COMPRESSIBLE MOMENTUM EQUATION (X COMPONENT):"
print momentum_x
print "\n"

print "SOURCE TERM FOR COMPRESSIBLE MOMENTUM EQUATION (Y COMPONENT):"
print momentum_y
print "\n"

print "SOURCE TERM FOR INTERNAL ENERGY EQUATION:"
print vfrac_c*rho_c*u_c[0]*diff(e, x) + vfrac_c*rho_c*u_c[1]*diff(e, y) + p*diff(vfrac_c*u_c[0], x) + p*diff(vfrac_c*u_c[1], y)
print "\n"
