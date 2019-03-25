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

filename='single_phase_frozen_flow_test_378.vtu' # This corresponds to approximately t = 3.77e-4 seconds (normalised time of tau = 4 in the paper by Miura & Glass (1982))
vt=vtktools.vtu(filename)

# Time and coordinates
t=vt.GetScalarField('Time')[0]
xyz=vt.GetLocations()
x=xyz[:,0]

# Solution fields
p=vt.GetScalarField('Pressure')
uvw_air=vt.GetVectorField('Velocity')
u_air=uvw_air[:,0]
rho_air=vt.GetScalarField('Density')
ie_air=vt.GetScalarField('InternalEnergy')

# Plot the numerical results.
# Note that we have divided the results by the reference velocity/pressure/density/internal energy to normalise them.
# x is normalised by the reference length.
pylab.figure(1)
pylab.subplot(321)
pylab.plot( x/0.0271002710027, p/1.013e5,'b.', label="Numerical")

pylab.subplot(322)
print(len(x))
print(len(u_air))
pylab.plot( x[:len(x)]/0.0271002710027, u_air/286.980353992,'b.', label="Numerical (Air)")

pylab.subplot(323)
pylab.plot( x/0.0271002710027, rho_air/1.23,'b.', label="Numerical")

pylab.subplot(324)
pylab.plot( x/0.0271002710027, ie_air/82357.7235772,'b.', label="Numerical (Air)")

# Now work out the frozen flow ('analytical') solutions
sol=numpy.array([shocktube.solution(xi,4.0) for xi in x/0.0271002710027])
p=sol[:,0]
u=sol[:,1]
rho=sol[:,2]
ie=p/rho/(shocktube.gamma-1.0)
  
pylab.subplot(321)
pylab.plot( x/0.0271002710027, p,'-g', label='Frozen flow')
pylab.title('Normalised Pressure')
pylab.legend(loc=2)
  
pylab.subplot(322)
pylab.plot( x/0.0271002710027, u,'g-', label='Frozen flow (Air)')
pylab.title('Normalised Velocity')
pylab.legend(loc=2)

pylab.subplot(323)
pylab.plot( x/0.0271002710027, rho,'g-', label='Frozen flow')
pylab.title('Normalised Density of Air')
pylab.legend(loc=2)

pylab.subplot(324)
pylab.plot( x/0.0271002710027, ie,'g-', label='Frozen flow')
pylab.title('Normalised Internal Energy of Air')
pylab.legend(loc=2)

pylab.savefig('single_phase_frozen_flow_test.png')

pylab.show()
