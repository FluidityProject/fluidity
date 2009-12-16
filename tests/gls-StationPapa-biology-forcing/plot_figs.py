#!/usr/bin/env python

from fluidity_tools import stat_parser
from numpy import arange,concatenate,array,argsort,append,fromfile
from scipy.interpolate import interp1d
import os
import sys
import vtktools
import math
from pylab import *
from matplotlib.ticker import MaxNLocator
import re 
import string

#### taken from http://www.codinghorror.com/blog/archives/001018.html  #######
def sort_nicely( l ): 
  """ Sort the given list in the way that humans expect. 
  """ 
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key ) 
##############################################################################


filelist = sys.argv[1:]
files = []
for file in filelist:
   if (string.find(file, 'checkpoint') != -1):
     continue
   elif (string.find(file, '_0.vtu') != -1):
     continue
   else:
     files.append(file)

sort_nicely(files)     
NN = 120

x0 = -3358350+100
y0 = -2.35154e+06
tke0 = 1.0e-5
den0 = 0.03


times = zeros((size(files),NN))
depths = zeros((size(files),NN))
temps = zeros((size(files),NN))
speeds = zeros((size(files),NN))
density = zeros((size(files),NN))
salinity = zeros((size(files),NN))
phyto =  zeros((size(files),NN))
nutrients =  zeros((size(files),NN))
detritus =  zeros((size(files),NN))
zoo =  zeros((size(files),NN))
# temperature
mld_times = []
mld_depths = []
#density
mld_times2 = []
mld_depths2 = []
#TKE
mld_times3 = []
mld_depths3 = []


ii=0
for file in files:
   print file
   try:
     os.stat(file)
   except:
     print "No such file: %s" % file
     sys.exit(1)
   num = int(file.split(".vtu")[0].split('_')[-1])
   u=vtktools.vtu(file)
   time = u.GetScalarField('Time')
   tt = time[0]
   # skip is less than 15 minutes of simulated time
   if (tt < 900):
      continue
   temp = u.GetScalarField('Temperature')
   pos = u.GetLocations()
   vel = u.GetVectorField('Velocity')
   kk = u.GetScalarField('GLSTurbulentKineticEnergy')
   den = u.GetScalarField('Density')
   sal = u.GetScalarField('Salinity')
   phy = u.GetScalarField('Phytoplankton')
   z = u.GetScalarField('Zooplankton')
   nut = u.GetScalarField('Nutrient')
   det = u.GetScalarField('Detritus')
   xyzkkden = []
   for i in range(0,len(temp)):
      if( (abs(pos[i,0] - x0) < 0.1) & (abs(pos[i,1] - y0) < 0.1) ):
         xyzkkden.append((pos[i,0],pos[i,1],-pos[i,2] + 4.88594e+06,temp[i],vel[i,0],vel[i,1],vel[i,2],kk[i],1000*den[i],sal[i],phy[i],z[i],nut[i],det[i]))
   xyzkkarr = vtktools.arr(xyzkkden)
   III = argsort(xyzkkarr[:,2])
   xyzkkarrsort = xyzkkarr[III,:]
   sst = xyzkkarrsort[0,3]
   sden = xyzkkarrsort[0,8]
   #print sst
   #print sden
   ggg = ((xyzkkarrsort[:,3]) > (sst-0.1)).nonzero()
   denggg = ((xyzkkarrsort[:,8]) <= (sden+den0)).nonzero()
   gggkk = ((xyzkkarrsort[:,7]) > tke0).nonzero()
   
   # MLD based on temperature
   if( (size(ggg) > 0 ) ):
      LL = ggg[-1][-1]
      zz = xyzkkarrsort[:,2]
      if (LL+1 < size(zz)):
        zza = zz[LL+1]
        kea = xyzkkarrsort[LL+1,3]
      else:
        zza = zz[LL]
        kea = xyzkkarrsort[LL,3]
      zzb = zz[LL]
      keb = xyzkkarrsort[LL,3]
      mld_times.append(tt/(24*60*60))
      mld_depths.append((zza-(zza-zzb)*(((keb)-(sst-0.1))/((keb)-(kea)))))
      mld_depths[-1] = -1 * mld_depths[-1]
      if (mld_depths[-1] < -500.0):
          mld_depths[-1] = -500.0

   # MLD based on density       
   if( (size(denggg) > 0 ) ):
      LL = denggg[-1][-1]
      zz = xyzkkarrsort[:,2]
      if (LL+1 < size(zz)):
        zza = zz[LL+1]
        kea = xyzkkarrsort[LL+1,8]
      else:
        zza = zz[LL]
        kea = xyzkkarrsort[LL,8]
      zzb = zz[LL]
      keb = xyzkkarrsort[LL,8]
      mld_times2.append(tt/(24*60*60))
      mld_depths2.append((zza-(zzb-zza)*(((sden+den0)-kea)/(kea-keb))))
      mld_depths2[-1] = -1 * mld_depths2[-1]
      if (mld_depths2[-1] < -500.0):
          mld_depths2[-1] = -500.0

   # MLD based on TKE       
   if( (size(gggkk) > 0 ) ):
      LL = gggkk[-1][-1]
      zz = xyzkkarrsort[:,2]
      if (LL+1 < size(zz)):
        zza = zz[LL+1]
        kea = xyzkkarrsort[LL+1,7]
      else:
        zza = zz[LL]
        kea = xyzkkarrsort[LL,7]
      zzb = zz[LL]
      keb = xyzkkarrsort[LL,7]
      mld_times3.append(tt/(24*60*60))
      mld_depths3.append((zza-(zza-zzb)*((keb-tke0)/((keb)-(kea)))))
      mld_depths3[-1] = -1 * mld_depths3[-1]
      if (mld_depths3[-1] < -500.0):
          mld_depths3[-1] = -500.0

   for jj in arange(NN):
      times[ii,jj] = time[0]
      depths[ii,jj] = -xyzkkarrsort[jj,2]
      temps[ii,jj] = xyzkkarrsort[jj,3]
      density[ii,jj] = xyzkkarrsort[jj,8]/1000
      salinity[ii,jj] = xyzkkarrsort[jj,9]
      phyto[ii,jj] = xyzkkarrsort[jj,10]
      zoo[ii,jj] = xyzkkarrsort[jj,11]
      nutrients[ii,jj] = xyzkkarrsort[jj,12]
      detritus[ii,jj] = xyzkkarrsort[jj,13]
      speeds[ii,jj] = sqrt(xyzkkarrsort[jj,4]**2 + xyzkkarrsort[jj,5]**2)
   ii = ii+1

# grab data from stat file
file = file.split("_")[0]
file = file+".stat"
stat=stat_parser( file )
time = stat["ElapsedTime"]["value"]/(24*60*60)


###########
# PHYSICS #
###########

fig=figure()

# The code below plots a 2 column, 4 row plot of graphs
#     
#    temp             l2norm temp
#   
#    velocity         UML depth
#
#    Diff/Viscos      TKE
#
#    Density          Temp profile

ax = fig.add_subplot(4,2,1)
ax.plot(time,stat["Fluid"]["Temperature"]["max"],'b')
ax.plot(time,stat["Fluid"]["Temperature"]["min"],'r')
xlabel('Time (days)')
ylabel('Max temperature')

ax = fig.add_subplot(4,2,2)
ax.plot(time,stat["Fluid"]["Temperature"]["l2norm"],'b')
xlabel('Time (days)')
ylabel('l2norm temperature')

ax = fig.add_subplot(4,2,3)
ax.plot(time,stat["Fluid"]["Velocity%magnitude"]["max"],'k')
ax.plot(time,stat["Fluid"]["Velocity%1"]["max"],'b')
ax.plot(time,stat["Fluid"]["Velocity%2"]["max"],'r')
ax.plot(time,stat["Fluid"]["Velocity%1"]["min"],'b')
ax.plot(time,stat["Fluid"]["Velocity%2"]["min"],'r')
xlabel('Time (days)')
ylabel('Max speed, u, v')

ax = fig.add_subplot(4,2,4)
ax.plot(mld_times,mld_depths,'b')
ax.plot(mld_times2,mld_depths2,'r')
ax.plot(mld_times3,mld_depths3,'k')
xlabel('Time (days)')
ylabel('UML depth (T-b, D-r, TKE-k)')

ax = fig.add_subplot(4,2,5)
ax.plot(time,stat["Fluid"]["GLSVerticalDiffusivity"]["max"],'b')
ax.plot(time,stat["Fluid"]["GLSVerticalViscosity"]["max"],'r')
ax.set_ylim(0,1)
xlabel('Time (days)')
ylabel('Max KH (blue) and KM (red)')

ax = fig.add_subplot(4,2,6)
ax.plot(time,stat["Fluid"]["GLSTurbulentKineticEnergy"]["max"],'b')
xlabel('Time (days)')
ylabel('Max TKE')


ax = fig.add_subplot(4,2,7)
cs=ax.pcolor(times/86400.0,depths,density)
ax.plot(mld_times2,mld_depths2,'w')
pp=colorbar(cs)
ax.set_ylim(-150, 0)
pp.set_label("Density")
xlabel('Time (days)')
ylabel('Depth (m)')

ax = fig.add_subplot(4,2,8)
cs=ax.pcolor(times/86400.0,depths,temps)
ax.plot(mld_times,mld_depths,'w')
ax.set_ylim(-150, 0)
pp=colorbar(cs)
pp.set_label("Temperature")
xlabel('Time (days)')
ylabel('Depth (m)')

###########
# BIOLOGY #
###########

# PLots 3 rows by 2 columns
#
# Max Zoo/Plankton    Photo. Radiation
#
# Detritus            Nutrients
#
# Plankton            Zooplankton
#

fig_bio = figure()
ax = fig_bio.add_subplot(3,2,1)
ax.plot(time,stat["Fluid"]["Zooplankton"]["max"],'b')
ax.plot(time,stat["Fluid"]["Phytoplankton"]["max"],'r')
xlabel('Time (days)')
ylabel('Max Plankton (Z-b, P-r)')

ax = fig_bio.add_subplot(3,2,2)
ax.plot(time,stat["Fluid"]["PhotosyntheticRadiation"]["max"],'b')
xlabel('Time (days)')
ylabel('Max Photosynthetic Radiation (W/m^2)')

ax = fig_bio.add_subplot(3,2,3)
cs=ax.pcolor(times/86400.0,depths,detritus)
pp=colorbar(cs)
ax.set_ylim(-150, 0)
pp.set_label("Detritus")
xlabel('Time (days)')
ylabel('Depth (m)')

ax = fig_bio.add_subplot(3,2,4)
cs=ax.pcolor(times/86400.0,depths,nutrients)
pp=colorbar(cs)
ax.plot(mld_times3,mld_depths3,'w')
ax.set_ylim(-150, 0)
cs.set_clim(0, 0.015)
pp.set_label("Nutrients")
xlabel('Time (days)')
ylabel('Depth (m)')

ax = fig_bio.add_subplot(3,2,5)
cs=ax.pcolor(times/86400.0,depths,phyto)
ax.plot(mld_times3,mld_depths3,'w')
ax.set_ylim(-150, 0)
pp=colorbar(cs)
pp.set_label("Phytoplankton")
xlabel('Time (days)')
ylabel('Depth (m)')

ax = fig_bio.add_subplot(3,2,6)
cs=ax.pcolor(times/86400.0,depths,zoo)
ax.plot(mld_times3,mld_depths3,'w')
ax.set_ylim(-150, 0)
pp=colorbar(cs)
pp.set_label("Zooplankton")
xlabel('Time (days)')
ylabel('Depth (m)')

show()
