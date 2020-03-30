#!/usr/bin/env python3
import vtktools
from numpy import concatenate, arange, newaxis, ones, array
from pylab import plot, show, figure, axis
import glob
import operator

vtus = glob.glob("1material_shocktube_*[0-9].vtu")

colxz=ones((39,1))*0.005
coly=arange(-0.95,0.99,0.05)[:,newaxis]
coordinates=concatenate((colxz,coly,colxz),1)

nums = []
for i in range(len(vtus)):
  nums.append([int(vtus[i].split(".vtu")[0].split("_")[-1]), i])
nums.sort(key=operator.itemgetter(0))

for i in range(len(vtus)):
  vtufile=vtktools.vtu(vtus[nums[i][-1]])
  newprobedpressure=vtktools.vtu.ProbeData(vtufile,coordinates,'Pressure')
  #newprobedvelocity=vtktools.vtu.ProbeData(vtufile,coordinates,'Velocity')
  newprobeddensity=vtktools.vtu.ProbeData(vtufile,coordinates,'Density')
  newprobedmatnrg=vtktools.vtu.ProbeData(vtufile,coordinates,'InternalEnergy')
  newnormmatnrg=(newprobedmatnrg-min(newprobedmatnrg))/(max(newprobedmatnrg)-min(newprobedmatnrg))
  #plot(coordinates[:,1],newnormmatnrg,coordinates[:,1],newprobeddensity,coordinates[:,1],newprobedpressure,coordinates[:,1],newprobedvelocity[:,1])
  figure(i)
  plot(coordinates[:,1],newnormmatnrg,coordinates[:,1],newprobeddensity,coordinates[:,1],newprobedpressure)
  axis([-1, 1, 0, 1])

show()
