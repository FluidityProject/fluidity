#!/usr/bin/env python3

# This sript plots up the RMS velocity:
import matplotlib.pyplot as pylab
import numpy
from fluidity_tools import stat_parser as stat
pylab.rc('font', size=32)

pylab.subplot(211)

statfile24="particle_rayleigh-taylor-mu10-24.stat"
statfile48="particle_rayleigh-taylor-mu10-48.stat"

# Data from paper:
cnd_data   = numpy.loadtxt("rms_data.txt")
particle_data = numpy.loadtxt("rms_data_4.txt")

#Scaling factor to account fo dimensions of model
scaling_factor = numpy.sqrt(1./0.9142)

# First plot 24x24 case:
pylab.plot(stat(statfile24)["ElapsedTime"]["value"][:],
           stat(statfile24)["Buoyant"]["Velocity%magnitude"]["l2norm"][:]*scaling_factor,'r', linestyle='-',lw=4.0, label="Particle-24")

# Next plot 48x48 case:
pylab.plot(stat(statfile48)["ElapsedTime"]["value"][:],
           stat(statfile48)["Buoyant"]["Velocity%magnitude"]["l2norm"][:]*scaling_factor,'b', linestyle='-',lw=4.0, label="Particle-48")

# Plot Vankekan benchmark value as line for comparison:
pylab.plot(cnd_data[:,0],cnd_data[:,1],'k-',lw=4.0, label="VK chain")
pylab.plot(particle_data[:,0],particle_data[:,1],'k--',lw=4.0, label="VK particle")

# Standard plot:python
pylab.ylabel(r"RMS Velocity (N-D)",size=38)
pylab.xlim(0,1000)
pylab.ylim(0.0,0.01)
pylab.ticklabel_format(axis='y', style='sci',scilimits=(-2,2),useOffset=False)
pylab.grid()
pylab.tick_params(axis='y',labelsize=32)
pylab.tick_params(axis='x',labelbottom=False)

pylab.subplot(212)
# Data from paper:
cnd_data   = numpy.loadtxt("entrainment_data.txt")
particle_data = numpy.loadtxt("entrainment_data_4.txt")

#Scaling factor to account fo dimensions of model
ld=0.9142
de=0.2
scaling_factor = ld*de

# First plot 24x24 case:
pylab.plot(stat(statfile24)["ElapsedTime"]["value"][:],
           stat(statfile24)["Buoyant"]["Entrainment"]["integral"][:]/scaling_factor,'r', linestyle='-',lw=4.0,label="Particle-24")

# Next plot 48x48 case:
pylab.plot(stat(statfile48)["ElapsedTime"]["value"][:],
           stat(statfile48)["Buoyant"]["Entrainment"]["integral"][:]/scaling_factor,'b', linestyle='-',lw=4.0,label="Particle-48")

# Plot Vankekan benchmark value as line for comparison:
pylab.plot(cnd_data[:,0],cnd_data[:,1],'k-',lw=4.0, label="VK chain")
pylab.plot(particle_data[:,0],particle_data[:,1],'k--',lw=4., label="VK particle")

# Standard plot:
pylab.xlabel(r"Time (N-D)",size=38)
pylab.ylabel(r"Entrainment (% of Material)",size=38)
pylab.xlim(0,1000)
pylab.ylim(0,1)
pylab.grid()
pylab.tick_params(axis='both',labelsize=32)
pylab.tick_params(axis='both',labelsize=32)
pylab.legend(loc='lower right',ncol=1, fontsize=36)

pylab.gcf().set_tight_layout(True)
pylab.gcf().set_size_inches(18, 20)
pylab.gcf().savefig('RT.png', bbox_inches='tight')
