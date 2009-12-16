#!/usr/bin/env python

from numpy import arange,concatenate,array,argsort
import os
import sys
import vtktools
import math
from pylab import *
from matplotlib.ticker import MaxNLocator
import re 


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



x0 = 500.0


tke0 = 1.0e-5

times = []
depths = []
for file in filelist:
   print file
   try:
     os.stat(file)
   except:
     print "No such file: %s" % file
     sys.exit(1)

   num = int(file.split(".vtu")[0].split('_')[-1])
   if (num<2):
     continue


   u=vtktools.vtu(file)


   time = u.GetScalarField('Time')
   tt = time[0]
   kk = u.GetScalarField('GLSTurbulentKineticEnergy')
   pos = u.GetLocations()
   
   xyzkk = []
   for i in range(0,len(kk)):
      if( abs(pos[i,0] - x0) < 0.1 ):
         xyzkk.append((pos[i,0],-pos[i,1],pos[i,2],(kk[i])))

   xyzkkarr = vtktools.arr(xyzkk)
   III = argsort(xyzkkarr[:,1])
   xyzkkarrsort = xyzkkarr[III,:]
   ggg = ((xyzkkarrsort[:,3])>tke0).nonzero()
   LL = ggg[-1][-1]
   zz = -xyzkkarrsort[:,1]
   zza = zz[LL+1]
   kea = xyzkkarrsort[LL+1,3]
   zzb = zz[LL]
   keb = xyzkkarrsort[LL,3]
   print tt/3600,(zza-(zza-zzb)*(((keb)-tke0)/((keb)-(kea))))
   times.append(tt/3600)
   depths.append(-1.0*(zza-(zza-zzb)*(((keb)-tke0)/((keb)-(kea)))))

   
fig = figure() 
plot(times,depths,'b-')

times2 = arange(0, 30*60*60, 60)
Dm = 1.05*1.0e-2*(1.0/sqrt(0.01))*sqrt(times2);
plot(times2/3600.0,Dm,'k--')




tke0 = 1.0e-4

times = []
depths = []
for file in filelist:
   print file
   try:
     os.stat(file)
   except:
     print "No such file: %s" % file
     sys.exit(1)

   num = int(file.split(".vtu")[0].split('_')[-1])
   if (num<2):
     continue


   u=vtktools.vtu(file)


   time = u.GetScalarField('Time')
   tt = time[0]
   kk = u.GetScalarField('GLSTurbulentKineticEnergy')
   pos = u.GetLocations()
   
   xyzkk = []
   for i in range(0,len(kk)):
      if( abs(pos[i,0] - x0) < 0.1 ):
         xyzkk.append((pos[i,0],-pos[i,1],pos[i,2],(kk[i])))

   xyzkkarr = vtktools.arr(xyzkk)
   III = argsort(xyzkkarr[:,1])
   xyzkkarrsort = xyzkkarr[III,:]
   ggg = ((xyzkkarrsort[:,3])>tke0).nonzero()
   LL = ggg[-1][-1]
   zz = -xyzkkarrsort[:,1]
   zza = zz[LL+1]
   kea = xyzkkarrsort[LL+1,3]
   zzb = zz[LL]
   keb = xyzkkarrsort[LL,3]
   print tt/3600,(zza-(zza-zzb)*(((keb)-tke0)/((keb)-(kea))))
   times.append(tt/3600)
   depths.append(-1.0*(zza-(zza-zzb)*(((keb)-tke0)/((keb)-(kea)))))

   
plot(times,depths,'r-')


show()
