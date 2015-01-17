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
u = [0.25*cos(x)*cos(y) - x*cos(y), sin(y)] # Velocity of the fluid phase
#u = [sin(x)*cos(y), sin(y)*sin(x) - cos(x)*sin(y)] # Velocity of particle phase 1
#u = [cos(x)*cos(y), sin(y)*sin(x) + sin(x)*sin(y)] # Velocity of particle phase 2
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

buoyancy_x = 0
buoyancy_y = 9.8*vfrac*rho

# Fluid-particle drag term
is_particle = True # Are we considering the particle phase, or the fluid phase?

vfrac_p = var('vfrac_p')
mu_f = var('mu_f')
d_p = var('d_p')
u_f = [0.25*cos(x)*cos(y) - x*cos(y), sin(y)]
#u_p = [sin(x)*cos(y), sin(y)*sin(x) - cos(x)*sin(y)] # Velocity of particle phase 1
u_p = [cos(x)*cos(y), sin(y)*sin(x) + sin(x)*sin(y)] # Velocity of particle phase 2

if(is_particle == True):
   fluid_particle_drag_x = -(18.0*mu_f*vfrac_p/(d_p**2))*(u_p[0] - u_f[0])
   fluid_particle_drag_y = -(18.0*mu_f*vfrac_p/(d_p**2))*(u_p[1] - u_f[1])
else:
   fluid_particle_drag_x = (18.0*mu_f*vfrac_p/(d_p**2))*(u_p[0] - u_f[0])
   fluid_particle_drag_y = (18.0*mu_f*vfrac_p/(d_p**2))*(u_p[1] - u_f[1])

# Particle-particle drag term
# An extension of the Syamlal (1985) correlation by Neri et al. (2003)
if(is_particle == True):

   # Particle phase parameters
   e = 0.6
   vfrac_i = 0.05
   vfrac_j = 0.15
   d_i = 1.0
   d_j = 0.1
   rho_i = 3.0
   rho_j = 2.0
   phi_i = 0.7
   phi_j = 0.4
   alpha = 0.9

   u_i = [sin(x)*cos(y), sin(y)*sin(x) - cos(x)*sin(y)] # Particle phase 1
   u_j = [cos(x)*cos(y), sin(y)*sin(x) + sin(x)*sin(y)] # Particle phase 2
   
   a = sqrt(d_j/d_i)
   X_i = vfrac_i/(vfrac_i + vfrac_j)
   
   if(X_i <= phi_i/(phi_i + (1.0-phi_i)*phi_j)):
      # "Maximum volume fraction of a random closely packed mixture" (Neri et al., 2003).
      vfrac_ij = ((phi_i - phi_j) + (1.0-a)*(1.0-phi_i)*phi_j)*((phi_i + (1.0-phi_j)*phi_i)/phi_i)*X_i + phi_j
   else:
      vfrac_ij = (1.0 - a)*(phi_i + (1.0 - phi_i)*phi_j)*(1 - X_i) + phi_i

   F = (3.0*vfrac_ij**(1.0/3.0) + (vfrac_i + vfrac_j)**(1.0/3.0)) / (2.0*(vfrac_ij**(1.0/3.0) - (vfrac_i + vfrac_j)**(1.0/3.0)))

   particle_particle_drag_x = F*alpha*(1.0+e)*rho_i*vfrac_i*rho_j*vfrac_j*( (d_i + d_j)**2 / (rho_i*d_i**3 + rho_j*d_j**3) )*sqrt((u_i[0] - u_j[0])**2 + (u_i[1] - u_j[1])**2) * (u_j[0] - u_i[0])
   particle_particle_drag_y = F*alpha*(1.0+e)*rho_i*vfrac_i*rho_j*vfrac_j*( (d_i + d_j)**2 / (rho_i*d_i**3 + rho_j*d_j**3) )*sqrt((u_i[0] - u_j[0])**2 + (u_i[1] - u_j[1])**2) * (u_j[1] - u_i[1])

else:
   particle_particle_drag_x = 0
   particle_particle_drag_y = 0

momentum_x = mass_x + advection_x + pressure_gradient_x + stress_x + buoyancy_x
momentum_y = mass_y + advection_y + pressure_gradient_y + stress_y + buoyancy_y

print momentum_x
print momentum_y

print "\n"

print fluid_particle_drag_x
print fluid_particle_drag_y

print "\n"

print particle_particle_drag_x
print particle_particle_drag_y

