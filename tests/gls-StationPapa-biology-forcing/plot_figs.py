#!/usr/bin/env python

from fluidity_tools import stat_parser
import numpy
from numpy import arange,concatenate,array,argsort,append,fromfile
from scipy.interpolate import interp1d
import os
import sys
import vtktools
import math
import string
from pylab import *
from matplotlib.ticker import MaxNLocator
import re 
from datetime import datetime, timedelta
import numpy.ma as MA


###################################
#Usage:
# - edit the config below
# - run using: ./plot_figs ./*.vtu
#   where *.vtu points to a directory full of VTUs you want to plot
#   If the VTUs are in the same directory as the script, MAKE SURE you do ./*.vtu
#
# Figures will be saved as PNGs in the directory containing your VTUs
###################################


####################
#                  #
#      CONFIG      #
#                  #
####################
# change default font sizes
params = {
          'legend.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'font.size' : 14,
          'axes.labelsize' : 14,
          'font.size' : 14,
          'text.fontsize' : 14,
          'figure.subplot.bottom' : 0.3
}
rcParams.update(params)

# max number of times
TT = 3652

# Origin (and where data is taken from)
x0 = 0
y0 = 0
z0 = 0
# this is the maximum depth that will be plotted
z_end = 150 # metres from surface

# dz
# The scripts extracts data on a regular grid from the VTUs
# This specifies the spacing required. Too small a number
# will make this script run very, very slowly.
# Manual suggests 1m
res = 1

# location of GOTM output
gotm_output = 'ows_papa.nc'

# definitions of UML depth
tke0 = 1.0e-5
den0 = 0.03
temp0 = 0.1

# data file that contains measuered UML depths
uml_data = "mld_stationpapa.csv"

# start date/time
# YYY-MM-DD HH:MM:SS
start = "1970-01-01 00:00:00"
end = "1978-06-01 00:00:00"
# this is used to match real data to the simulation data

# clip data?
# This plots all subplots (of the same variable) to the same colour scale
clip_data = True
# This doesn't do anything in this script as there is on one dataset!


##############################################################################
##############################################################################
#   Main script
##############################################################################
##############################################################################

# this may need tweaking if you do < 1m resolution
FNN = int(z_end/res)

#########
## Some helpful functions
#########

#### taken from http://www.codinghorror.com/blog/archives/001018.html  #######
def sort_nicely( l ): 
  """ Sort the given list in the way that humans expect. 
  """ 
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key ) 
##############################################################################

def plot_2d_data_with_mld(fluidity_data,file_name,axis_label):

    fig = figure(figsize=(15,8),dpi=90)
    fig.subplots_adjust(hspace=0.6)
    ax = fig.add_subplot(111)
    adjusted = fluidity_data
    cs=ax.pcolormesh(real_dates,depths,adjusted)
    pp=colorbar(cs)
    ax.plot(date2num(real_dates_mld_kk),mld_depths_kk,'w', alpha=0.5)
    dateFmt = mpl.dates.DateFormatter('%m/%Y')
    ax.xaxis.set_major_formatter(dateFmt)
    monthsLoc = mpl.dates.MonthLocator(interval=3)
    ax.xaxis.set_major_locator(monthsLoc)
    ax.set_xlim(date2num(start_date),date2num(end_date))
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    ax.set_ylim(-z_end, 0)
    ax.set_xlim(start_date,end_date)
    pp.set_label(axis_label + ' (ICOM)')
    xlabel('Date (mm/yyyy)')
    ylabel('Depth (m)')
    savefig(working_dir+'/mld_' +file_name, dpi=90,format='png')

def plot_2d_data(fluidity_data,file_name,axis_label):

    fig = figure(figsize=(15,8),dpi=90)
    fig.subplots_adjust(hspace=0.6)
    ax = fig.add_subplot(111)
    adjusted = fluidity_data
    cs=ax.pcolormesh(real_dates,depths,adjusted)
    pp=colorbar(cs)
    dateFmt = mpl.dates.DateFormatter('%m/%Y')
    ax.xaxis.set_major_formatter(dateFmt)
    monthsLoc = mpl.dates.MonthLocator(interval=3)
    ax.xaxis.set_major_locator(monthsLoc)
    ax.set_xlim(date2num(start_date),date2num(end_date))
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    ax.set_ylim(-z_end, 0)
    ax.set_xlim(start_date,end_date)
    pp.set_label(axis_label + ' (ICOM)')
    xlabel('Date (mm/yyyy)')
    ylabel('Depth (m)')
    savefig(working_dir+'/' +file_name, dpi=90,format='png')

######################################################################

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

times = zeros((size(files),FNN))
depths = zeros((size(files),FNN))
temps = zeros((size(files),FNN))
speeds = zeros((size(files),FNN))
density = zeros((size(files),FNN))
salinity = zeros((size(files),FNN))
tke = zeros((size(files),FNN))
NN2 = zeros((size(files),FNN))
vert_diff = zeros((size(files),FNN))
vert_visc = zeros((size(files),FNN))
len_scale = zeros((size(files),FNN))
u_vel = zeros((size(files),FNN))
v_vel = zeros((size(files),FNN))
phytoplankton = zeros((size(files),FNN))
zooplankton = zeros((size(files),FNN))
detritus = zeros((size(files),FNN))
nutrients = zeros((size(files),FNN))

mld_times_kk = []
mld_depths_kk = []
mld_times_temp = []
mld_depths_temp = []
mld_times_den = []
mld_depths_den = []

# set up new coords, where we want data
# vertical profile from surface to bottom with 1m vertical res
coords = numpy.array([x0,y0,z0])
nPoints = 0
for z in arange(z0-1,z0-z_end-res,-res): 
    coords = numpy.vstack(((coords,[x0, y0, z])))
    nPoints = nPoints + 1

ii=0
working_dir=file.rsplit("/",1)[0]
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

   # grab the data e need from the VTUs and shove in an array
   temp = u.ProbeData(coords,'Temperature')
   vel = u.ProbeData(coords,'Velocity')
   kk = u.ProbeData(coords,'GLSTurbulentKineticEnergy')
   nn = u.ProbeData(coords,'GLSBuoyancyFrequency')
   den = u.ProbeData(coords,'Density')
   sal = u.ProbeData(coords,'Salinity')
   vdiff = u.ProbeData(coords,'GLSVerticalDiffusivity')
   vvisc = u.ProbeData(coords,'GLSVerticalViscosity')
   length = u.ProbeData(coords,'GLSLengthScale')
   pp = u.ProbeData(coords,'Phytoplankton')
   zp = u.ProbeData(coords,'Zooplankton')
   det = u.ProbeData(coords,'Detritus')
   nutr = u.ProbeData(coords,'Nutrient')
   xyz_data = []
   for i in range(0,len(temp)):
        xyz_data.append((coords[i,0],coords[i,1],-coords[i,2]+z0,temp[i],vel[i,0],vel[i,1],vel[i,2],kk[i],1000*den[i],sal[i],nn[i],vdiff[i],vvisc[i],length[i],pp[i],zp[i],det[i],nutr[i]))
   

   # sorted the array based on depth
   xyz_data = vtktools.arr(xyz_data)
   III = argsort(xyz_data[:,2])
   xyz_data_sorted = xyz_data[III,:]

   # Surface values
   sst = xyz_data_sorted[0,3]
   sden = xyz_data_sorted[0,8]

   # grab any values where the UML condition is met for temperature, density and TKE
   uml_temp = ((xyz_data_sorted[:,3]) > (sst-temp0)).nonzero()
   uml_den = ((xyz_data_sorted[:,8]) <= (sden+den0)).nonzero()
   uml_kk = ((xyz_data_sorted[:,7]) > tke0).nonzero()

   # work out ML depth based on temperature
   if( (size(uml_temp) > 0 ) ):
      LL = uml_temp[-1][-1]
      zz = xyz_data_sorted[:,2]
      if (LL+1 < size(zz)):
        zza = zz[LL+1]
        kea = xyz_data_sorted[LL+1,3]
      else:
        zza = zz[LL]
        kea = xyz_data_sorted[LL,3]
      zzb = zz[LL]
      keb = xyz_data_sorted[LL,3]
      mld_times_temp.append(tt/(24*60*60))
      mld_depths_temp.append((zza-(zza-zzb)*(((keb)-(sst-temp0))/((keb)-(kea)))))
      mld_depths_temp[-1] = mld_depths_temp[-1] * -1.
      if (mld_depths_temp[-1] < -500.0):
          mld_depths_temp[-1] = -500.0

   # ...on density
   if( (size(uml_den) > 0 ) ):
      LL = uml_den[-1][-1]
      zz = xyz_data_sorted[:,2]
      if (LL+1 < size(zz)):
        zza = zz[LL+1]
        kea = xyz_data_sorted[LL+1,8]
      else:
        zza = zz[LL]
        kea = xyz_data_sorted[LL,8]
      zzb = zz[LL]
      keb = xyz_data_sorted[LL,8]
      mld_times_den.append(tt/(24*60*60))
      mld_depths_den.append((zza-(zzb-zza)*(((sden+den0)-kea)/(kea-keb))))
      mld_depths_den[-1] = mld_depths_den[-1] * -1.
      if (mld_depths_den[-1] < -500.0):
          mld_depths_den[-1] = -500.0

   # MLD based on TKE       
   if( (size(uml_kk) > 0 ) ):
      LL = uml_kk[-1][-1]
      zz = xyz_data_sorted[:,2]
      if (LL+1 < size(zz)):
        zza = zz[LL+1]
        kea = xyz_data_sorted[LL+1,7]
      else:
        zza = zz[LL]
        kea = xyz_data_sorted[LL,7]
      zzb = zz[LL]
      keb = xyz_data_sorted[LL,7]
      mld_times_kk.append(tt/(24*60*60))
      mld_depths_kk.append((zza-(zza-zzb)*((keb-tke0)/((keb)-(kea)))))
      mld_depths_kk[-1] = mld_depths_kk[-1] * -1.
      if (mld_depths_kk[-1] < -500.0):
          mld_depths_kk[-1] = -500.0


   # place the rest of the data in arrays ready for display
   for jj in arange(FNN):
      times[ii,jj] = time[0]
      depths[ii,jj] = -xyz_data_sorted[jj,2]
      temps[ii,jj] = xyz_data_sorted[jj,3]
      speeds[ii,jj] = sqrt(xyz_data_sorted[jj,4]**2 + xyz_data_sorted[jj,5]**2)
      u_vel[ii,jj] = xyz_data_sorted[jj,4]
      v_vel[ii,jj] = xyz_data_sorted[jj,5]
      density[ii,jj] = xyz_data_sorted[jj,8]/1000
      salinity[ii,jj] = xyz_data_sorted[jj,9]
      tke[ii,jj] = xyz_data_sorted[jj,7]
      NN2[ii,jj] = xyz_data_sorted[jj,10]
      vert_diff[ii,jj] = xyz_data_sorted[jj,11]
      vert_visc[ii,jj] = xyz_data_sorted[jj,12]
      len_scale[ii,jj] = xyz_data_sorted[jj,13]
      phytoplankton[ii,jj] = xyz_data_sorted[jj,14]
      zooplankton[ii,jj] = xyz_data_sorted[jj,15]
      detritus[ii,jj] = xyz_data_sorted[jj,16]
      nutrients[ii,jj] = xyz_data_sorted[jj,17]
   ii = ii+1




# measured MLD data
from numpy import loadtxt
##############################################
# moving onto real data for the UML depths...#
##############################################
uml_dates,uml_depths=loadtxt(uml_data,delimiter=",",unpack=True,skiprows=1,usecols=(1,5),converters = {1: strpdate2num('%d/%m/%y')})
# convert date to number of days since start date so we can plot all the data on the same x-axis
uml_times = uml_dates
start_date = datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
end_date = datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
uml_depths = -1.0 * uml_depths

# create proper dates
real_dates = times
real_dates_mld_den = []
real_dates_mld_kk = []
i = 0
for time in times:
    real_dates[i,:] = date2num(start_date + timedelta(seconds=time[0]))
    i += 1

for time in mld_times_den:
    real_dates_mld_den.append(start_date + timedelta(days=time))
for time in mld_times_kk:
    real_dates_mld_kk.append(start_date + timedelta(days=time))

# grab data from stat file
file = file.rsplit("_",1)[0]
file = file+".stat"
stat=stat_parser( file )
time = stat["ElapsedTime"]["value"]
real_dates_stat = []
for t in time:
    real_dates_stat.append(date2num(start_date + timedelta(seconds=t)))

max_time = end_date


# Plot stuff
# temperature
plot_2d_data(temps,"Temperature.png","Temperature")
# salinity
plot_2d_data(salinity,"Salinity.png","Salinity")
# BuoyancyFreq
plot_2d_data(NN2,"NN.png","Buoyancy Freq.")
# TKE
plot_2d_data(tke,"TKE.png","TKE")
# Viscosity
plot_2d_data(vert_visc,"Visc.png","Viscosity")
# Diffusivity
plot_2d_data(vert_diff,"Diff.png","Diffusivity")
# Length scale
plot_2d_data(len_scale,"L.png","Length Scale")
# Phytoplankton
plot_2d_data(phytoplankton,"phyto.png","Phytoplankton Conc.")
# Zooplankton
plot_2d_data(zooplankton,"zoo.png","Zooplankton Conc.")
# Nutrients
plot_2d_data(nutrients,"nutrients.png","Nutrients Conc.")
# Detritus
plot_2d_data(detritus,"detritus.png","Detritus Conc.")

# Plot stuff - with MLD
# temperature
plot_2d_data_with_mld(temps,"Temperature.png","Temperature")
# salinity
plot_2d_data(salinity,"Salinity.png","Salinity")
# BuoyancyFreq
plot_2d_data_with_mld(NN2,"NN.png","Buoyancy Freq.")
# TKE
plot_2d_data_with_mld(tke,"TKE.png","TKE")
# Viscosity
plot_2d_data_with_mld(vert_visc,"Visc.png","Viscosity")
# Diffusivity
plot_2d_data_with_mld(vert_diff,"Diff.png","Diffusivity")
# Length scale
plot_2d_data_with_mld(len_scale,"L.png","Length Scale")
# Phytoplankton
plot_2d_data_with_mld(phytoplankton,"phyto.png","Phytoplankton Conc.")
# Zooplankton
plot_2d_data_with_mld(zooplankton,"zoo.png","Zooplankton Conc.")
# Nutrients
plot_2d_data_with_mld(nutrients,"nutrients.png","Nutrients Conc.")
# Detritus
plot_2d_data_with_mld(detritus,"detritus.png","Detritus Conc.")


fig_mld = figure(figsize=(25,10),dpi=90)
ax = fig_mld.add_subplot(111)
ax.plot(uml_times,uml_depths,'sg',label="Obs.",markersize=3)
ax.plot(date2num(real_dates_mld_den),mld_depths_den,'r', label="ICOM")
# limit x-axis to simulation data, not real data
ax.set_xlim(mld_times_temp[0],mld_times_temp[-1])
# invert and constrain y-axis
bottom, top = ax.get_ylim()
dateFmt = mpl.dates.DateFormatter('%m/%Y')
ax.xaxis.set_major_formatter(dateFmt)
monthsLoc = mpl.dates.MonthLocator(interval=3)
ax.xaxis.set_major_locator(monthsLoc)
ax.set_xlim(date2num(start_date),date2num(end_date))
labels = ax.get_xticklabels()
for label in labels:
    label.set_rotation(30) 
ax.set_ylim(0, -z_end)
ax.set_xlim(start_date,end_date)
ax.yaxis.grid(True)
legend(loc=0)
xlabel('Date (mm/yyyy)')
ylabel('UML depth (m)')
savefig(working_dir+'/MLD.png', dpi=90,format='png')


fig_diff = figure(figsize=(15,5),dpi=90)
ax = fig_diff.add_subplot(111)
ax.plot(real_dates_stat,stat["Fluid"]["Phytoplankton"]["max"],'b', label="Phytoplankton")
ax.plot(real_dates_stat,stat["Fluid"]["Zooplankton"]["max"],'r', label="Zooplankton")
dateFmt = mpl.dates.DateFormatter('%m/%Y')
ax.xaxis.set_major_formatter(dateFmt)
monthsLoc = mpl.dates.MonthLocator()
ax.xaxis.set_major_locator(monthsLoc)
ax.set_xlim(date2num(start_date),date2num(end_date))
labels = ax.get_xticklabels()
for label in labels:
    label.set_rotation(30) 
ax.set_xlim(start_date, end_date)
ax.yaxis.grid(True)
legend(loc=0)
xlabel('Date (mm/yyyy)')
ylabel('Max Phyto- and Zooplankton Conc.')
savefig(working_dir+'/MaxPhytoZooPlankton.png', dpi=90,format='png')

