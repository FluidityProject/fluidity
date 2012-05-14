import pylab
import numpy
import shocktube
import vtktools
params = {'text.fontsize': 6,
         'legend.fontsize': 8,
         'xtick.labelsize': 8,
         'ytick.labelsize': 8,
         'lines.markersize': 6,
         'axes.titlesize': 'small'}
pylab.rcParams.update(params)

filename='mphase_dusty_gas_shock_tube_311.vtu' # This is at t = 3.11e-4 seconds (normalised time of tau = 4 in the paper by Miura & Glass (1982))
vt=vtktools.vtu(filename)

# Time and coordinates
t=vt.GetScalarField('Air::Time')[0]
xyz=vt.GetLocations()
x=xyz[:,0]

# Solution fields
p=vt.GetScalarField('Air::Pressure')
uvw_air=vt.GetVectorField('Air::Velocity')
u_air=uvw_air[:,0]
uvw_dust=vt.GetVectorField('Dust::Velocity')
u_dust=uvw_dust[:,0]
rho_air=vt.GetScalarField('Air::Density')
rho_dust=vt.GetScalarField('Dust::Density')
ie=vt.GetScalarField('Air::InternalEnergy')
vfrac=vt.GetScalarField('Dust::PhaseVolumeFraction')

# First normalise and plot the numerical results.
pylab.figure(1)
pylab.subplot(321)
pylab.plot( x/0.0271002710027, p/1.013e5,'b.', label="Numerical")

pylab.subplot(322)
print len(x)
print len(u_air)
pylab.plot( x[:len(x)]/0.0271002710027, u_air/339.559734079,'b.', label="Numerical (Air)")
pylab.plot( x[:len(x)]/0.0271002710027, u_dust/339.559734079,'r.', label="Numerical (Particles)")

pylab.subplot(323)
pylab.plot( x/0.0271002710027, rho_air/1.23,'b.', label="Numerical")

pylab.subplot(324)
pylab.plot( x/0.0271002710027, rho_dust,'b.', label="Numerical")
pylab.title('Density of Particles')
pylab.legend(loc=2)

pylab.subplot(325)
pylab.plot( x/0.0271002710027, ie/205894.308943,'b.', label="Numerical")

# Multiply the PhaseVolumeFraction by the Dust Density to get the mass concentration, and then normalise.
pylab.subplot(326)
pylab.plot( x/0.0271002710027, vfrac*2500.0/1.23,'b.', label="Numerical")
pylab.title('Mass Concentration from PhaseVolumeFraction of Particles')
pylab.legend(loc=2)
  
# Now work out the frozen flow ('analytical') solutions
sol=numpy.array([shocktube.solution(xi,4.0) for xi in x/0.0271002710027])
p=sol[:,0]
u=sol[:,1]
rho=sol[:,2]
ie=p/rho/(shocktube.gamma-1.0)
  
pylab.subplot(321)
pylab.plot( x/0.0271002710027, p*1.0e5/1.013,'-g', label='Frozen flow')
pylab.title('Normalised Pressure')
pylab.legend(loc=2)
  
pylab.subplot(322)
pylab.plot( x/0.0271002710027, u,'g-', label='Frozen flow (Air)')
pylab.title('Normalised Velocity')
pylab.legend(loc=2)

pylab.subplot(323)
pylab.plot( x/0.0271002710027, rho*1e5/1.23,'g-', label='Frozen flow')
pylab.title('Normalised Density of Air')
pylab.legend(loc=2)

pylab.subplot(325)
pylab.plot( x/0.0271002710027, ie*1e5/205894.308943,'g-', label='Frozen flow')
pylab.title('Normalised Internal Energy of Air')
pylab.legend(loc=2)

pylab.savefig('mphase_dusty_gas_shock_tube.png')

pylab.show()
