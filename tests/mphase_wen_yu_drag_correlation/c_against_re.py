# Plots the drag coefficient against the particle Reynolds number.
from math import sqrt, pi
from numpy import *
import pylab
from fluidity_tools import stat_parser

C_wen_yu = zeros(200000)
C_stokes = zeros(200000)
particle_Re = arange(0.001, 1000, 0.005)
for i in range(0, len(particle_Re)):
   # Drag coefficients for the Wen & Yu and Stokes drag correlations respectively.
   C_wen_yu[i] = (24.0/particle_Re[i])*(1.0 + 0.15*particle_Re[i]**0.687)
   C_stokes[i] = (24.0/particle_Re[i])
   
s = stat_parser("./mphase_wen_yu_drag_correlation.stat")
numerical_particle_Re_wen_yu = s["Tephra"]["ParticleReynoldsNumber"]["max"][-1]
numerical_C_wen_yu = s["Tephra"]["DragCoefficient"]["max"][-1]

pylab.loglog(particle_Re, C_stokes, "-r", label="Stokes drag correlation")
pylab.loglog(particle_Re, C_wen_yu, "-g", label="Wen & Yu drag correlation")
pylab.loglog(numerical_particle_Re_wen_yu, numerical_C_wen_yu, "*b", label="Numerical result")
pylab.legend(loc=1)
pylab.xlabel("ParticleReynoldsNumber")
pylab.ylabel("DragCoefficient")

pylab.show()