#!/usr/bin/env python

from numpy import arange,concatenate
import os
import sys
sys.path.append('/data2/BP-CFD/BP/run/flu-150808/tools')
import vtktools
import math
from pylab import *
from matplotlib.ticker import MaxNLocator


filelist = sys.argv[1:]

print filelist

min_z = 0.0
max_z = -50.0
step_z = -0.25

pts = []
zz = []
for z in arange(min_z, max_z+step_z, step_z):
      pts.append((0.5, 0.5, z))
      zz.append(z)



for file in filelist:
   print file
   try:
     os.stat(file)
   except:
     print "No such file: %s" % file
     sys.exit(1)

   num = int(file.split(".vtu")[0].split('_')[-1])


   u=vtktools.vtu(file)
##   u.GetFieldNames()

   fig = figure()
   temp = u.ProbeData(vtktools.arr(pts), "Temperature")
   #vert_diff = u.ProbeData(vtktools.arr(pts), "VerticalDiffusivity")
   #vert_visc = u.ProbeData(vtktools.arr(pts), "VerticalViscosity")
   KineticEnergy = u.ProbeData(vtktools.arr(pts), "KineticEnergy")
   uvw = u.ProbeData(vtktools.arr(pts), "Velocity")


##   for i in range(0,len(zz)):
##      print i,zz[i],KineticEnergy[i]/2.0

   tke0 = 1.0e-5
   ggg = ((KineticEnergy/2.0)<tke0).nonzero()
   gggg = array(ggg[0])
   print gggg[0]-1,zz[gggg[0]-1],KineticEnergy[gggg[0]-1]/2.0
   print gggg[0],zz[gggg[0]],KineticEnergy[gggg[0]]/2.0
   print (zz[gggg[0]-1]-(zz[gggg[0]-1]-zz[gggg[0]])*(((KineticEnergy[gggg[0]]/2.0)-tke0)/((KineticEnergy[gggg[0]]/2.0)-(KineticEnergy[gggg[0]-1]/2.0))))

   ax = fig.add_subplot(221)
   ax.plot(temp*10*2.0e-4,zz,'b')
   ax.xaxis.set_major_locator(MaxNLocator(5))
   xlabel('Buoyancy')
   ylabel('Depth')

   ax = fig.add_subplot(222)
   ax.plot(uvw[:,0],zz,'b')
   ax.xaxis.set_major_locator(MaxNLocator(5))
   xlabel('U Velocity')
   ylabel('Depth')

   ax = fig.add_subplot(223)
   ax.plot(KineticEnergy/2.0,zz,'b')
   ax.xaxis.set_major_locator(MaxNLocator(5))
   xlabel('TKE: (q^2)/2')
   ylabel('Depth')

   ax = fig.add_subplot(224)
   #ax.plot(vert_diff,zz,'b--')
   #ax.plot(vert_visc,zz,'b')
   #ax.xaxis.set_major_locator(MaxNLocator(5))
   xlabel('Vertical Diffusivity (--) and Viscosity (-)')
   ylabel('Depth')

   show()
