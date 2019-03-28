from numpy import arange,concatenate,array,argsort,zeros
import os
import sys
import vtktools
import math
import re 
from math import sqrt
from scipy.interpolate import UnivariateSpline

#### taken from http://www.codinghorror.com/blog/archives/001018.html  #######
def sort_nicely( l ): 
  """ Sort the given list in the way that humans expect. 
  """ 
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key ) 
##############################################################################


# compute the mixed layer depth over time
def MLD(filelist):
  x0 = 0.
  tke0 = 1.0e-5
  last_mld = 0
  
  times = []
  depths = []
  Dm = []
  for file in filelist:
     try:
       os.stat(file)
     except:
       print("No such file: %s" % file)
       sys.exit(1)
     
     u=vtktools.vtu(file)
     time = u.GetScalarField('Time')
     tt = time[0]
     kk = u.GetScalarField('GLSTurbulentKineticEnergy')
     pos = u.GetLocations()
     # ignore first 4 hours of simulaiton
     if (tt < 14400):
       continue

     xyzkk = []
     for i in range(0,len(kk)):
       if( abs(pos[i,0] - x0) < 0.1 ):
         xyzkk.append((pos[i,0],-pos[i,1],pos[i,2],(kk[i])))

     xyzkkarr = vtktools.arr(xyzkk)
     III = argsort(xyzkkarr[:,1])
     xyzkkarrsort = xyzkkarr[III,:]
     # march down the column, grabbing the last value above tk0 and the first 
     # one less than tke0. Interpolate between to get the MLD
     kea = 1000
     keb = 0
     zza = 0
     zzb = 0
     for values in xyzkkarrsort:
        if (values[3] > tke0):
            kea = values[3]
            zza = -values[1]
        if (values[3] < tke0):
            keb = values[3]
            zzb = -values[1]
            break

     # the MLD is somewhere between these two values - let's estimate half way!
     mld = (zzb+zza)/2.
     if (last_mld == mld):
        continue

     times.append(tt/3600)
     depths.append(-1.0*mld)
     last_mld = mld
     Dm.append(1.05*0.00988211768*(1.0/sqrt(0.01))*sqrt(tt))


  return times, depths, Dm

