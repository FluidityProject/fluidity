#!/usr/bin/env python

from __future__ import print_function

import glob
import sys
import os
import vtktools
import numpy
import pylab
import re
import extract_data
from math import log

def get_filelist(sample, start):

    def key(s):
        return int(s.split('_')[-1].split('.')[0])
   
    list = glob.glob("*vtu")
    list = [l for l in list if 'checkpoint' not in l]
    vtu_nos = [float(s.split('_')[-1].split('.')[0]) for s in list]
    vals = zip(vtu_nos, list)
    vals.sort()
    unzip = lambda l:tuple(apply(zip,l))
    vtu_nos, list = unzip(vals)
    shortlist = []

    for file in list:
      try:
        os.stat(file)
      except:
        f_log.write("No such file: %s" % files)
        sys.exit(1)

      ##### Start at the (start+1)th file.
      ##### Add every nth file by taking integer multiples of n; limit at 10 vtus max.
      vtu_no = float(file.split('_')[-1].split('.')[0])
      #if ((max(vtu_nos)-start)/sample > 10):
      #  sample=int((max(vtu_nos)-start)/10)
      
      if vtu_no > start:
        if (vtu_no%sample==0):
          shortlist.append(file)
        ##### Append final file if a large number of files remain.
        elif vtu_no==len(vtu_nos)-1 and (max(vtu_nos)-sample/4.0)>vtu_no:
          shortlist.append(file)

    return shortlist

#### taken from http://www.codinghorror.com/blog/archives/001018.html  #######
def tryint(s):
    try:
        return int(s)
    except:
        return s
    
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
##############################################################################
# There are shorter and more elegant version of the above, but this works
# on CX1, where this test might be run...

###################################################################

# Reattachment length:
def reatt_length(filelist, zarray):

  print("Calculating reattachment point locations using change of x-velocity sign\n")

  nums=[]; results=[]; files = []
  ##### check for no files
  if (len(filelist) == 0):
    print("No files!")
    sys.exit(1)
  for file in filelist:
    try:
      os.stat(file)
    except:
      print("No such file: %s" % file)
      sys.exit(1)
    files.append(file)
  sort_nicely(files)

  for file in files:
    ##### Read in data from vtu
    datafile = vtktools.vtu(file)
    ##### Get time for plot:
    t = min(datafile.GetScalarField("Time"))
    print(file, ', elapsed time = ', t)
    if(t<0.):
      continue
    else:
      print("extracting data...")
      ##### points near bottom surface, 0 < x < 25
      x2array=[]; pts=[]; no_pts = 52; offset = 0.01
      x = 0.0
      for i in range(1, no_pts):
        x2array.append(x)
        for j in range(len(zarray)):
          pts.append((x, zarray[j], offset))
        x += 0.5

    x2array = numpy.array(x2array)
    pts = numpy.array(pts)

    ##### Get x-velocity on bottom boundary
    uvw = datafile.ProbeData(pts, "AverageVelocity")
    u = uvw[:,0]
    u = u.reshape([x2array.size,zarray.size])
    pts=pts.reshape([x2array.size,zarray.size,3])

    ##### Find all potential reattachment points:
    points = []
    for j in range(len(u[0,:])):
      for i in range(len(u[:,0])-1):
        ##### Hack to ignore division by zero entries in u.
        ##### All u should be nonzero away from boundary!
        if((u[i,j] / u[i+1,j]) < 0. and u[i+1,j] > 0. and not numpy.isinf(u[i,j] / u[i+1,j])):
          ##### interpolate between nodes
          p = x2array[i] + (x2array[i+1]-x2array[i]) * (0.0-u[i,j]) / (u[i+1,j]-u[i,j])
          ##### Ignore spurious corner points
          if(p>1.0):
            points.append(p)
            ##### We have our first point on this plane so...
            break

    ##### This is the spanwise-averaged reattachment point:
    if (len(points)>0):
      avpt = sum(points) / len(points)
    else:
      avpt = 0.0
    print('spanwise averaged reattachment point: ', avpt)
    ##### Get time for plot:
    t = min(datafile.GetScalarField("Time"))
    results.append([avpt,t])

  return results

#########################################################################

# Velocity profiles:
def velo(filelist,xarray,zarray,yarray):

  print("\nRunning mean velocity profile script on files at times...\n")
  ##### check for no files
  if (len(filelist) < 0):
    print("No files!")
    sys.exit(1)

  ##### create array of points
  pts=[]
  for i in range(len(xarray)):
    for j in range(len(zarray)):
      for k in range(len(yarray)):
        pts.append([xarray[i], zarray[j], yarray[k]])
  pts=numpy.array(pts)

  ##### Create output array of correct shape
  profiles=numpy.zeros([xarray.size, yarray.size], float)

  file = filelist[-1]
  #for file in filelist:
  datafile = vtktools.vtu(file)
  # Get time
  t = min(datafile.GetScalarField("Time"))
  print(file, ', elapsed time = ', t)

  ##### Get x-velocity
  uvw = datafile.ProbeData(pts, "AverageVelocity")
  umax = 1.55
  u = uvw[:,0]/umax
  u = u.reshape([xarray.size,zarray.size,yarray.size])

  ##### Spanwise averaging
  usum = numpy.zeros([xarray.size,yarray.size],float)
  usum = numpy.array(usum)
  for i in range(len(zarray)):
    uav = u[:,i,:]
    uav = numpy.array(uav)
    usum += uav
  usum = usum / len(zarray)
  profiles[:,:] = usum

  print("\n...Finished extracting data.\n")
  return profiles

#########################################################################

def plot_length(rl):
  ##### Plot time series of reattachment length

  av_length = sum(rl[:,0]) / len(rl[:,0])
  avg = numpy.zeros([len(rl[:,0])])
  avg[:] = av_length
  Lemoinkim = numpy.zeros([len(rl[:,0])])
  Lemoinkim[:]=6.28

  plot1 = pylab.figure()
  pylab.title("Time series of reattachment length")
  pylab.xlabel('Time (s)')
  pylab.ylabel('Reattachment Length (L/h)')
  pylab.plot(rl[:,1], rl[:,0], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="solid")
  pylab.plot(rl[:,1], Lemoinkim, linestyle="dashed")
  pylab.legend(("Fluidity","Le-Moin-Kim DNS"), loc="best")
  pylab.axis([min(rl[:,1]),max(rl[:,1]),min(rl[:,0])-0.5,max(rl[:,0])+0.5])
  pylab.savefig("reatt_len_3d.pdf")
  return

#########################################################################

def plot_velo(vprofiles,xarray,yarray):

  # get profiles from ERCOFTAC data
  y4,U4,y6,U6,y10,U10,y19,U19 = extract_data.ercoftacvelocityprofiles()
  # get profiles from Le&Moin data
  Le_y4,Le_u4,jd_y4,jd_u4,Le_y6,Le_u6,jd_y6,jd_u6,Le_y10,Le_u10,jd_y10,jd_u10,Le_y19,Le_u19,jd_y19,jd_u19 = extract_data.velocityprofileslemoin()

  ##### Plot velocity profiles at different points behind step using pylab(matplotlib)
  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle("Evolution of mean U-velocity", fontsize=20)

  size = 15

  ax = pylab.subplot(141)
  ax.plot(vprofiles[0,:],yarray, linestyle="solid")
  ax.plot(U4,y4, linestyle="dashed")
  ax.plot(jd_u4,jd_y4, linestyle="none",marker='o',color='black')
  ax.set_title('(a) x/h='+str(xarray[0]), fontsize=16)
  pylab.legend(('Fluidity',"Le&Moin DNS","Jovic&Driver expt"),loc="upper left")
  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  bx = pylab.subplot(142, sharex=ax, sharey=ax)
  bx.plot(vprofiles[1,:],yarray, linestyle="solid")
  bx.plot(U6,y6, linestyle="dashed")
  bx.plot(jd_u6,jd_y6, linestyle="none",marker='o',color='black')
  bx.set_title('(a) x/h='+str(xarray[1]), fontsize=16)
  #bx.grid("True")
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(143, sharex=ax, sharey=ax)
  cx.plot(vprofiles[2,:],yarray, linestyle="solid")
  cx.plot(U10,y10, linestyle="dashed")
  cx.plot(jd_u10,jd_y10, linestyle="none",marker='o',color='black')
  cx.set_title('(a) x/h='+str(xarray[2]), fontsize=16)
  #bx.grid("True")
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(144, sharex=ax, sharey=ax)
  dx.plot(vprofiles[3,:],yarray, linestyle="solid")
  dx.plot(U19,y19, linestyle="dashed")
  dx.plot(jd_u19,jd_y19, linestyle="none",marker='o',color='black')
  dx.set_title('(a) x/h='+str(xarray[3]), fontsize=16)
  #bx.grid("True")
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  pylab.axis([-0.25, 1., 0., 3.])
  bx.set_xlabel('Normalised mean U-velocity (U/Umax)', fontsize=24)
  ax.set_ylabel('y/h', fontsize=24)

  pylab.savefig("velo_profiles_3d.pdf")
  return

#########################################################################

def main():
    ##### Only process every nth file by taking integer multiples of n:
    filelist = get_filelist(sample=1, start=0)

    ##### Points to generate profiles:
    xarray = numpy.array([4.0, 6.0, 10.0, 19.0])
    zarray = numpy.linspace(0.0,4.0,41)
    yarray = numpy.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0])

    ##### Call reattachment_length function
    reattachment_length = numpy.array(reatt_length(filelist, zarray))
    numpy.save("reatt_length", reattachment_length)
    plot_length(reattachment_length)

    ##### Call velo function
    zarray = numpy.array([2.0])
    vprofiles = velo(filelist, xarray, zarray, yarray)
    numpy.save("velo_profiles", vprofiles)
    print("Generating plot of velocity profiles.")
    plot_velo(vprofiles,xarray,yarray)

    print("\nAll done.\n")

if __name__ == "__main__":
    sys.exit(main())

