#!/usr/bin/python3

# This script plots up the Nusselt number for both the 24x24 and 48x48 cases:
import pylab
from fluidity_tools import stat_parser as stat
import numpy

# Stafiles:
statfile12="particle_rayleigh-taylor-mu10-12.stat"
statfile24="particle_rayleigh-taylor-mu10-24.stat"

# Data from paper:
cnd_data   = numpy.loadtxt("entrainment_data.txt")
particle_data = numpy.loadtxt("entrainment_data_4.txt")

#Scaling factor to account fo dimensions of model
ld=0.9142
de=0.2
scaling_factor = ld*de

# First plot 24x24 case:
pylab.plot(stat(statfile12)["ElapsedTime"]["value"][:],
           stat(statfile12)["Buoyant"]["Entrainment"]["integral"][:]/scaling_factor,'r', linestyle='-',lw=4.0,label="Particle-12")

# Next plot 48x48 case:
pylab.plot(stat(statfile24)["ElapsedTime"]["value"][:],
           stat(statfile24)["Buoyant"]["Entrainment"]["integral"][:]/scaling_factor,'b', linestyle='-',lw=4.0,label="Particle-24")

# Plot benchmark value as line for comparison:
pylab.plot(cnd_data[:,0],cnd_data[:,1],'k-',lw=4.0, label="VK chain")
pylab.plot(particle_data[:,0],particle_data[:,1],'k--',lw=4., label="VK particle")

pylab.xlabel(r"Time")
pylab.ylabel(r"Entrainment")
pylab.xlim(0,1000)
pylab.ylim(0,1)
pylab.grid()
pylab.legend(loc='lower right',ncol=1)


pylab.savefig("Entrainment.png")
