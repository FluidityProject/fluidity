import pylab
import numpy
import shocktube
import vtktools
from math import sqrt
f = open('/home/christian/Documents/scripts/rcparams.py', 'r')
exec(f.read())

filename='mphase_dusty_gas_shock_tube_378.vtu' # This corresponds to approximately t = 3.77e-4 seconds (normalised time of tau = 4 in the paper by Miura & Glass (1982))
vt=vtktools.vtu(filename)

# Time and coordinates
t=vt.GetScalarField('Air::Time')[0]
xyz=vt.GetLocations()
x=xyz[:,0]

# Solution fields
p=vt.GetScalarField('Air::PressureProjected')
uvw_air=vt.GetVectorField('Air::Velocity')
u_air=uvw_air[:,0]
uvw_dust=vt.GetVectorField('Dust::Velocity')
u_dust=uvw_dust[:,0]
rho_air=vt.GetScalarField('Air::DensityProjected')
rho_dust=vt.GetScalarField('Dust::Density')
ie_air=vt.GetScalarField('Air::InternalEnergyProjected')
ie_dust=vt.GetScalarField('Dust::InternalEnergyProjected')
vfrac=vt.GetScalarField('Dust::PhaseVolumeFraction')



# Get the single phase results
filename='single_phase_frozen_flow_test_378.vtu' # This corresponds to approximately t = 3.77e-4 seconds (normalised time of tau = 4 in the paper by Miura & Glass (1982))
vt=vtktools.vtu(filename)

# Solution fields
single_p=vt.GetScalarField('PressureProjected')
single_uvw_air=vt.GetVectorField('Velocity')
single_u_air=single_uvw_air[:,0]
single_rho_air=vt.GetScalarField('DensityProjected')
single_ie_air=vt.GetScalarField('InternalEnergyProjected')


# Plot the numerical results.
# Note that we have divided the results by the reference velocity/pressure/density/internal energy to normalise them.
# x is normalised by the reference length.

# Pressure
pylab.figure(num=1, dpi=500, facecolor='w', edgecolor='k')
pylab.plot( x/0.0271002710027, p/1.013e5,'b-', label="Numerical (Fluidity, Multiphase)")
pylab.plot( x/0.0271002710027, single_p/1.013e5,'g-', label="Numerical (Fluidity, Single-phase)")

# Velocity of Air
pylab.figure(num=2, dpi=500, facecolor='w', edgecolor='k')
pylab.plot( x[:len(x)]/0.0271002710027, u_air/339.559734079,'b-', label="Numerical (Fluidity, Multiphase) - Air")

# Velocity of Particles
pylab.plot( x[:len(x)]/0.0271002710027, u_dust/339.559734079,'r-', label="Numerical (Fluidity, Multiphase) - Particles")

pylab.plot( x[:len(x)]/0.0271002710027, single_u_air/339.559734079,'g-', label="Numerical (Fluidity, Single-phase) - Air")

# Density of Air
pylab.figure(num=4, dpi=500, facecolor='w', edgecolor='k')
pylab.plot( x/0.0271002710027, rho_air/1.23,'b-', label="Numerical (Fluidity, Multiphase)")
pylab.plot( x/0.0271002710027, single_rho_air/1.23,'g-', label="Numerical (Fluidity, Single-phase)")

# Internal energy of Air
pylab.figure(num=5, dpi=500, facecolor='w', edgecolor='k')
pylab.plot( x/0.0271002710027, ie_air/205894.308943089,'b-', label="Numerical (Fluidity, Multiphase) - Air")

# Internal energy of Particles
pylab.plot( x/0.0271002710027, ie_dust/205894.308943089,'r-', label="Numerical (Fluidity, Multiphase) - Particles")

pylab.plot( x/0.0271002710027, single_ie_air/205894.308943089,'g-', label="Numerical (Fluidity, Single-phase) - Air")

# Multiply the PhaseVolumeFraction by the Dust Density to get the mass concentration, and then normalise.
pylab.figure(num=7, dpi=500, facecolor='w', edgecolor='k')
pylab.plot( x/0.0271002710027, vfrac*2500.0/1.23,'b-', label="Numerical (Fluidity, Multiphase)")
pylab.legend(loc=2)


# Now work out the frozen flow ('analytical') solutions
sol=numpy.array([shocktube.solution(xi,4.0) for xi in x/0.0271002710027])
gamma = 1.4
p=sol[:,0]
u=sol[:,1]/sqrt(gamma)
rho=sol[:,2]
ie=p/rho/(shocktube.gamma-1.0)/2.5
  
pylab.figure(1)
pylab.plot( x/0.0271002710027, p,'g--', label='Frozen flow')
pylab.legend(loc=1)
pylab.axis([-10, 10, 0, 12])
pylab.xlabel(r"Normalised position ($x$/$l$)")
pylab.ylabel(r"Normalised pressure ($p$/$p^{\mathrm{ref}}$)")
pylab.figtext(0.16, 0.16, "(a)")
pylab.savefig('mphase_dusty_gas_shock_tube_explicit_pressure.pdf')
  
pylab.figure(2)
pylab.plot( x/0.0271002710027, u,'g--', label='Frozen flow - Air')
pylab.legend(loc=1)
pylab.axis([-10, 10, -0.2, 1.4])
pylab.xlabel(r"Normalised position ($x$/$l$)")
pylab.ylabel(r"Normalised speed ($|\mathbf{u}_i|$/$|\mathbf{u}^{\mathrm{ref}}|$)")
pylab.figtext(0.16, 0.16, "(b)")
pylab.savefig('mphase_dusty_gas_shock_tube_explicit_velocity.pdf')

pylab.figure(4)
pylab.plot( x/0.0271002710027, rho,'g--', label='Frozen flow')
pylab.legend(loc=1)
pylab.axis([-10, 10, 0, 12])
pylab.xlabel(r"Normalised position ($x$/$l$)")
pylab.ylabel(r"Normalised density of air ($\rho_c$/$\rho^{\mathrm{ref}}$)")
pylab.figtext(0.16, 0.16, "(c)")
pylab.savefig('mphase_dusty_gas_shock_tube_explicit_density_air.pdf')

pylab.figure(5)
pylab.plot( x/0.0271002710027, ie,'g--', label='Frozen flow - Air')
pylab.legend(loc=1)
pylab.axis([-10, 10, 0.4, 2.0])
pylab.xlabel(r"Normalised position ($x$/$l$)")
pylab.ylabel(r"Normalised specific internal energy ($e_i$/$e^{\mathrm{ref}}$)")
pylab.figtext(0.16, 0.16, "(d)")
pylab.savefig('mphase_dusty_gas_shock_tube_explicit_ie.pdf')

pylab.figure(7)
pylab.legend(loc=1)
pylab.axis([-10, 10, 0, 3])
pylab.xlabel(r"Normalised position ($x$/$l$)")
pylab.ylabel(r"Normalised mass concentration ($\alpha_d\rho_d$/$\rho^{\mathrm{ref}}$)")
pylab.figtext(0.16, 0.16, "(e)")
pylab.savefig('mphase_dusty_gas_shock_tube_explicit_mass_concentration_particles.pdf')

