#!/usr/bin/env python3

from numpy import arange,concatenate,array,argsort
import os
import sys
import vtktools
import math
from pylab import *
from matplotlib.ticker import MaxNLocator
import re 
from scipy.interpolate import UnivariateSpline

#### taken from http://www.codinghorror.com/blog/archives/001018.html  #######
def sort_nicely( l ): 
  """ Sort the given list in the way that humans expect. 
  """ 
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key ) 
##############################################################################


filelist = sys.argv[1:]

sort_nicely(filelist)

x0 = 166.6
y0 = 50.
tke0 = 1.0e-5

last_mld = 0
fig = figure() 
times = []
depths = []
for file in filelist:
   print(file)
   try:
     os.stat(file)
   except:
     print("No such file: %s" % file)
     sys.exit(1)

   num = int(file.split(".vtu")[0].split('_')[-1])

   u=vtktools.vtu(file)


   time = u.GetScalarField('Time')
   tt = time[0]
   kk = u.GetScalarField('GLSTurbulentKineticEnergy')
   pos = u.GetLocations()
   if (tt < 1800):
      continue

   xyzkk = []
   for i in range(0,len(kk)):
      if( abs(pos[i,0] - x0) < 0.1 and abs(pos[i,1] - y0) < 0.1):
         xyzkk.append((pos[i,0],pos[i,1],-pos[i,2],(kk[i])))

   xyzkkarr = vtktools.arr(xyzkk)
   III = argsort(xyzkkarr[:,2])
   xyzkkarrsort = xyzkkarr[III,:]
   ggg = ((xyzkkarrsort[:,3])>tke0).nonzero()
   LL = ggg[-1][-1]
   zz = -xyzkkarrsort[:,2]
   zza = zz[LL+1]
   kea = xyzkkarrsort[LL+1,3]
   zzb = zz[LL]
   keb = xyzkkarrsort[LL,3]

   # We're using quite low resolution, so a better estimate of the MLD is
   # somewhere between the 2 layers
   mld = (zzb+zza)/2
   if (last_mld == mld):
      continue

   times.append(tt/3600)
   depths.append(-1.0*mld)
   last_mld = mld

   

plot(times,depths,'b-')

times2 = arange(0, 10*60*60, 60)
Dm = 1.05*0.0098*(1.0/sqrt(0.01))*sqrt(times2);
plot(times2/3600.0,Dm,'k--')

show()
