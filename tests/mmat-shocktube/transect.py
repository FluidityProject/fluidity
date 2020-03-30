#!/usr/bin/env python3
import vtktools
from numpy import concatenate, arange, newaxis, ones, array
from pylab import plot, show

colxz=ones((39,1))*0.005
coly=arange(-0.95,0.99,0.05)[:,newaxis]
coordinates=concatenate((colxz,coly,colxz),1)
vtufile=vtktools.vtu('1material_shocktube_15.vtu')
newprobedpressure=vtktools.vtu.ProbeData(vtufile,coordinates,'Pressure')
newprobedvelocity=vtktools.vtu.ProbeData(vtufile,coordinates,'Velocity')
newprobeddensity=vtktools.vtu.ProbeData(vtufile,coordinates,'Density')
newprobedmatnrg=vtktools.vtu.ProbeData(vtufile,coordinates,'MaterialInternalEnergy')
newprobedmatdens=vtktools.vtu.ProbeData(vtufile,coordinates,'MaterialDensity')
newnormmatnrg=(newprobedmatnrg-min(newprobedmatnrg))/(max(newprobedmatnrg)-min(newprobedmatnrg))
plot(coordinates[:,1],newnormmatnrg,coordinates[:,1],newprobedmatdens,coordinates[:,1],newprobeddensity,coordinates[:,1],newprobedpressure,coordinates[:,1],newprobedvelocity[:,1])
show()
