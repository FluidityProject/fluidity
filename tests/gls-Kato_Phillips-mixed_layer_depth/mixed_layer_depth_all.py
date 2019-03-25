#!/usr/bin/env python

from numpy import arange,concatenate,array,argsort
import os
import sys
import vtktools
import math
from pylab import *
from matplotlib.ticker import MaxNLocator
import re 
from scipy.interpolate import UnivariateSpline
import glob

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
     if (tt < 100):
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

     mld = zza
     if (last_mld == mld):
        continue

     times.append(tt/3600)
     depths.append(-1.0*mld)
     last_mld = mld

  return times, depths

path = sys.argv[1]

x0 = 0.
tke0 = 1.0e-5
files_to_look_through_ke = [
"Kato_Phillips-mld-k_e-CA",
"Kato_Phillips-mld-k_e-CB",
"Kato_Phillips-mld-k_e-GL",
"Kato_Phillips-mld-k_e-KC"
]
files_to_look_through_gen = [
"Kato_Phillips-mld-gen-CA",
"Kato_Phillips-mld-gen-CB",
"Kato_Phillips-mld-gen-GL",
"Kato_Phillips-mld-gen-KC"
]
files_to_look_through_kw = [
"Kato_Phillips-mld-k_w-CA",
"Kato_Phillips-mld-k_w-CB",
"Kato_Phillips-mld-k_w-GL",
"Kato_Phillips-mld-k_w-KC"
]
files_to_look_through_kkl = [
"Kato_Phillips-mld-k_kl-CA",
"Kato_Phillips-mld-k_kl-KC"
]
colours = ['r','g','b','#8000FF']

times2 = arange(0, 10, 0.1)
Dm = 1.05*1.0e-2*(1.0/sqrt(0.01))*sqrt((times2*60*60));

figke = figure(figsize=(9.172,4.5),dpi=90)
ax = figke.add_subplot(111)
i = 0
for simulation in files_to_look_through_ke:
    filelist = glob.glob(path+simulation+"*.vtu")
    sort_nicely(filelist)
    times, depths = MLD(filelist)
    ax.plot(times,depths,colours[i],label=simulation)
    i = i+1

ax.plot(times2,Dm,'k-',label='Analytical')
ax.set_ylim(ax.get_ylim()[::-1])
xlabel('Time (hours)')
ylabel('ML Depth (m)')
legend(loc=0)
savefig(path + '/ke.png', dpi=90,format='png')

figkw = figure(figsize=(9.172,4.5),dpi=90)
ax = figkw.add_subplot(111)
i = 0
for simulation in files_to_look_through_kw:
    filelist = glob.glob(path+simulation+"*.vtu")
    sort_nicely(filelist)
    times, depths = MLD(filelist)
    ax.plot(times,depths,colours[i],label=simulation)
    i = i+1

ax.plot(times2,Dm,'k-',label='Analytical')
ax.set_ylim(ax.get_ylim()[::-1])
xlabel('Time (hours)')
ylabel('ML Depth (m)')
legend(loc=0)
savefig(path + '/kw.png', dpi=90,format='png')


figgen = figure(figsize=(9.172,4.5),dpi=90)
ax = figgen.add_subplot(111)
i = 0
for simulation in files_to_look_through_gen:
    filelist = glob.glob(path+simulation+"*.vtu")
    sort_nicely(filelist)
    times, depths = MLD(filelist)
    ax.plot(times,depths,colours[i],label=simulation)
    i = i+1

ax.plot(times2,Dm,'k-',label='Analytical')
ax.set_ylim(ax.get_ylim()[::-1])
xlabel('Time (hours)')
ylabel('ML Depth (m)')
legend(loc=0)
savefig(path + '/gen.png', dpi=90,format='png')

figkkl = figure(figsize=(9.172,4.5),dpi=90)
ax = figkkl.add_subplot(111)
i = 0
for simulation in files_to_look_through_kkl:
    filelist = glob.glob(path+simulation+"*.vtu")
    sort_nicely(filelist)
    times, depths = MLD(filelist)
    ax.plot(times,depths,colours[i],label=simulation)
    i = i+1

ax.plot(times2,Dm,'k-',label='Analytical')
ax.set_ylim(ax.get_ylim()[::-1])
xlabel('Time (hours)')
ylabel('ML Depth (m)')
legend(loc=0)
savefig(path + '/kkl.png', dpi=90,format='png')


