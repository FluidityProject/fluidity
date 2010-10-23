#!/usr/bin/env python

import vtk
import glob
import sys
import os
import vtktools
import numpy
import pylab
from math import isinf
import Numeric

def get_filelist():

    def key(s):
        return int(s.split('_')[-1].split('.')[0])
   
    list = glob.glob("*.pvtu")
    list = [l for l in list if 'check' not in l]
    vtu_nos = [float(s.split('_')[-1].split('.')[0]) for s in list]
    vals = zip(vtu_nos, list)
    vals.sort()
    unzip = lambda l:tuple(apply(zip,l))
    vtu_nos, list = unzip(vals)

    return list

###################################################################

# Reattachment length:
def reatt_length(filelist, yarray, exclude_initial_results):

  print "Calculating reattachment point locations using change of x-velocity sign\n"
  print "Processing pvtu files...\n"
  nums=[]
  results=[]
  result=0.
  filecount=0

  ##### check for no files
  if (len(filelist) == 0):
    print "No files!"
    sys.exit(1)
  for file in filelist:
    try:
      os.stat(file)
    except:
      print "No such file: %s" % file
      sys.exit(1)
    num = int(file.split(".pvtu")[0].split('_')[-1])
    ##### Exclude data from simulation spin-up time
    if num > exclude_initial_results:
      nums.append(num)
    nums.sort()

  for num in nums:
    ##### counter for averaging
    filecount = filecount + 1
    file = "backward_facing_step_3d_"+str(num)+".pvtu"

    print file
    ##### Read in data from vtu
    reader = vtk.vtkXMLPUnstructuredGridReader();
    reader.SetFileName(file)
    reader.Update()
    data = reader.GetOutput()
    datafile = vtktools.vtu(file)

    ##### points near bottom surface, 0 < x < 25
    x2array=[]; pts=[]; no_pts = 52; offset = 0.1
    x = 0.0
    for i in range(1, no_pts):
      x2array.append(x)
      for j in range(len(yarray)):
        pts.append((x, yarray[j], offset))
      x += 0.5

    x2array = numpy.array(x2array)
    pts = numpy.array(pts)

    ##### Get x-velocity on bottom boundary
    uvw = datafile.ProbeData(pts, "Velocity")
    u = uvw[:,0]
    u = u.reshape([x2array.size,yarray.size])

    ##### Find all potential reattachment points:
    points = []
    for j in range(len(u[0,:])):
      for i in range(len(u[:,0])-1):
        ##### Hack to ignore division by zero entries in u.
        ##### All u should be nonzero away from boundary!
        if((u[i,j] / u[i+1,j]) < 0. and not isinf(u[i,j] / u[i+1,j])):
          ##### interpolate between nodes
          p = x2array[i] + (x2array[i+1]-x2array[i]) * (0.0-u[i,j]) / (u[i+1,j]-u[i,j])
          ##### Ignore spurious corner points
          if(p>0.1):
            points.append(p)

    ##### This is the spanwise-averaged reattachment point:
    avpt = sum(points) / len(points)
    print "spanwise-averaged reattachment point:", avpt, "\n"
    results.append(avpt)

  ## And this is the time-averaged reattachment point:
  result = sum(results) / len(results)
  return result

#########################################################################

# Velocity profiles:
def meanvelo(filelist,xarray,yarray,zarray):

  print "\nRunning velocity profile script on files at times...\n"
  ##### check for no files
  if (len(filelist) < 0):
    print "No files!"
    sys.exit(1)

  ##### create array of points
  pts=[]
  for i in range(len(xarray)):
    for j in range(len(yarray)):
      for k in range(len(zarray)):
        pts.append([xarray[i], yarray[j], zarray[k]])
  pts=numpy.array(pts)

  ##### Create output array of correct shape
  files = 3; filecount = 0
  profiles = numpy.zeros([files, xarray.size, zarray.size], float)

  for file in filelist:
    try:
      os.stat(file)
    except:
      f_log.write("No such file: %s" % files)
      sys.exit(1)

    reader = vtk.vtkXMLPUnstructuredGridReader();
    reader.SetFileName(file)
    reader.Update()
    data = reader.GetOutput()
    datafile = vtktools.vtu(file)
    t = min(datafile.GetScalarField("Time"))

    ##### only dump data for certain times
    while 4.5<t<5.5 or 9.5<t<10.5 or 49.5<t<50.5:

      print file, ', elapsed time = ', t
      ##### Get x-velocity
      uvw = datafile.ProbeData(pts, "Velocity")
      umax = max(abs(datafile.GetVectorField("Velocity")[:,0]))
      u = uvw[:,0]/umax
      u = u.reshape([xarray.size,yarray.size,zarray.size])

      ##### Spanwise averaging
      usum = numpy.zeros([xarray.size,zarray.size],float)
      usum = numpy.array(usum)
      for i in range(len(yarray)):
        uav = u[:,i,:]
        uav = numpy.array(uav)
        usum += uav
      usum = usum / len(yarray)
      profiles[filecount,:,:] = usum

      ##### reset time to something big to prevent infinite loop
      t = 100.
      filecount += 1

  print "\n...Finished extracting data.\n"
  return profiles

#########################################################################

##### Points to generate profiles:
xarray = numpy.array([0.0, 2.0, 4.0, 6.0])
yarray = numpy.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9])
zarray = numpy.array([0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])

##### Call functions defined above
filelist = get_filelist()
reattachment_length = reatt_length(filelist, yarray, exclude_initial_results=25)
print "\nTime-averaged reattachment length (in step heights): ", reattachment_length

profiles = meanvelo(filelist, xarray, yarray, zarray)
numpy.save("velo_profiles_3d_parallel", profiles)

#########################################################################

##### Plot velocity profiles using pylab(matplotlib)
pylab.figure(num=1, figsize = (16.5, 8.5))
pylab.suptitle('Evolution of U-velocity profiles downstream of backward facing step', fontsize=30)

size = 15

ax = pylab.subplot(141)
ax.plot(profiles[0,0,:],zarray, color='green', linestyle="dashed")
ax.plot(profiles[1,0,:],zarray, color='blue', linestyle="dashed")
ax.plot(profiles[2,0,:],zarray, color='red', linestyle="dashed")
pylab.legend(("5secs", "10secs", "50secs"), loc="upper left")
ax.set_title('(a) x/h = 0', fontsize=16)
ax.grid("True")
for tick in ax.xaxis.get_major_ticks():
  tick.label1.set_fontsize(size)
for tick in ax.yaxis.get_major_ticks():
  tick.label1.set_fontsize(size)

bx = pylab.subplot(142, sharex=ax, sharey=ax)
bx.plot(profiles[0,1,:],zarray, color='green', linestyle="dashed")
bx.plot(profiles[1,1,:],zarray, color='blue', linestyle="dashed")
bx.plot(profiles[2,1,:],zarray, color='red', linestyle="dashed")
pylab.legend(("5secs", "10secs", "50secs"), loc="upper left")
bx.set_title('(b) x/h = 2', fontsize=16)
bx.grid("True")
for tick in bx.xaxis.get_major_ticks():
  tick.label1.set_fontsize(size)
pylab.setp(bx.get_yticklabels(), visible=False)

cx = pylab.subplot(143, sharex=ax, sharey=ax)
cx.plot(profiles[0,2,:],zarray, color='green', linestyle="dashed")
cx.plot(profiles[1,2,:],zarray, color='blue', linestyle="dashed")
cx.plot(profiles[2,2,:],zarray, color='red', linestyle="dashed")
pylab.legend(("5secs", "10secs", "50secs"), loc="upper left")
cx.set_title('(c) x/h = 4', fontsize=16)
cx.grid("True")
for tick in cx.xaxis.get_major_ticks():
  tick.label1.set_fontsize(size)
pylab.setp(cx.get_yticklabels(), visible=False)

dx = pylab.subplot(144, sharex=ax, sharey=ax)
dx.plot(profiles[0,3,:],zarray, color='green', linestyle="dashed")
dx.plot(profiles[1,3,:],zarray, color='blue', linestyle="dashed")
dx.plot(profiles[2,3,:],zarray, color='red', linestyle="dashed")
pylab.legend(("5secs", "10secs", "50secs"), loc="upper left")
dx.set_title('(d) x/h = 6', fontsize=16)
dx.grid("True")
for tick in dx.xaxis.get_major_ticks():
  tick.label1.set_fontsize(size)
pylab.setp(dx.get_yticklabels(), visible=False)

pylab.axis([-0.25, 1., 0., 2.])
bx.set_xlabel('Normalised U-velocity (U/Umax)', fontsize=24)
ax.set_ylabel('z/h', fontsize=24)

pylab.savefig("velo_profiles_3d_parallel.pdf")
pylab.show()

print "\nAll done.\n"

