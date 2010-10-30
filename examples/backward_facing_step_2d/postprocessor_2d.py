#!/usr/bin/env python

import math
import vtk
import glob
import sys
import os
import vtktools
import numpy
import pylab

def get_filelist():

    def key(s):
        return int(s.split('_')[-1].split('.')[0])
   
    list = glob.glob("*.vtu")
    list = [l for l in list if 'check' not in l]
    vtu_nos = [float(s.split('_')[-1].split('.')[0]) for s in list]
    vals = zip(vtu_nos, list)
    vals.sort()
    unzip = lambda l:tuple(apply(zip,l))
    vtu_nos, list = unzip(vals)

    return list

###################################################################

# Reattachment length:
def reatt_length(filelist, exclude_initial_results):

  print "Calculating reattachment point locations using change of x-velocity sign\n"
  print "Processing pvtu files...\n"
  nums=[]; results=[]

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
    num = int(file.split(".vtu")[0].split('_')[-1])
    ##### Exclude data from simulation spin-up time
    if num > exclude_initial_results:
      nums.append(num)
    nums.sort()

  for num in nums:
    ##### counter for averaging
    file = "backward_facing_step_2d_"+str(num)+".vtu"

    print file
    ##### Read in data from vtu
    reader = vtk.vtkXMLUnstructuredGridReader();
    reader.SetFileName(file)
    reader.Update()
    data = reader.GetOutput()
    datafile = vtktools.vtu(file)

    ##### points near bottom surface, 0 < x < 20
    pts=[]; no_pts = 82; offset = 0.1
    x = 0.0
    for i in range(1, no_pts):
      pts.append((x, offset, 0.0))
      x += 0.25

    pts = numpy.array(pts)

    ##### Get x-velocity on bottom boundary
    uvw = datafile.ProbeData(pts, "Velocity")
    u = []
    u = uvw[:,0]

    ##### Find all potential reattachment points:
    points = []

    for i in range(len(u)-1):
      ##### Hack to ignore division by zero entries in u.
      ##### All u should be nonzero away from boundary!
      if((u[i] / u[i+1]) < 0. and not math.isinf(u[i] / u[i+1])):
        ##### interpolate between nodes
        p = pts[i] + (pts[i+1]-pts[i]) * (0.0-u[i]) / (u[i+1]-u[i])
        ##### Ignore spurious corner points
        if(p>0.1):
          points.append(p)

    ##### Get time for plot:
    t = min(datafile.GetScalarField("Time"))
    ##### Append actual reattachment point and time:
    results.append([points[0],t])

  return results

#########################################################################

# Velocity profiles:
def meanvelo(filelist,x,y):

  print "\nRunning velocity profile script on files at times...\n"
  ##### check for no files
  if (len(filelist) < 0):
    print "No files!"
    sys.exit(1)

  ##### create array of points
  pts=[]
  for i in range(len(x)):
    for j in range(len(y)):
      pts.append([x[i]+5.0, y[j], 0.0])

  pts=numpy.array(pts)
  ##### Create output array of correct shape
  files = 3; filecount = 0
  profiles = numpy.zeros([files, x.size, y.size], float)

  for file in filelist:
      try:
        os.stat(file)
      except:
        f_log.write("No such file: %s" % files)
        sys.exit(1)

      reader = vtk.vtkXMLUnstructuredGridReader();
      reader.SetFileName(file)
      reader.Update()
      data = reader.GetOutput()
      datafile = vtktools.vtu(file)
      t = min(datafile.GetScalarField("Time"))

      ##### only dump data for certain times
      while 4.8<t<5.2 or 9.8<t<10.2 or 49.8<t<50.0:
        print file, ', elapsed time = ', t
        ##### Get x-velocity
        uvw = datafile.ProbeData(pts, "Velocity")
        umax = max(abs(datafile.GetVectorField("Velocity")[:,0]))
        (ilen, jlen) = uvw.shape
        u = uvw[:,0]/umax
        u=u.reshape([x.size,y.size])
        profiles[filecount,:,:] = u

        ##### reset time to something big to prevent infinite loop
        t = 1000.
        filecount += 1

  print "\n...Finished writing data files.\n"
  return profiles

#########################################################################

def plot_length(reattachment_length):
  ##### Plot time series of reattachment length using pylab(matplotlib)
  plot1 = pylab.figure()
  pylab.title("Time series of reattachment length behind 2D step")
  pylab.xlabel('Time (s)')
  pylab.ylabel('Reattachment Length (L/h)')
  pylab.plot(reattachment_length[:,1], reattachment_length[:,0], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="solid")
  pylab.savefig("reattachment_length_2D.pdf")
  return

def plot_meanvelo(profiles):
  ##### Plot velocity profiles at different points behind step, and at 3 times using pylab(matplotlib)
  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle('Evolution of U-velocity profiles downstream of backward facing step', fontsize=30)

  size = 15

  ax = pylab.subplot(141)
  ax.plot(profiles[0,0,:],yarray, color='green', linestyle="dashed")
  ax.plot(profiles[1,0,:],yarray, color='blue', linestyle="dashed")
  ax.plot(profiles[2,0,:],yarray, color='red', linestyle="dashed")
  pylab.legend(("5secs", "10secs", "50secs"), loc="lower right")
  ax.set_title('(a) x/h = 0', fontsize=16)
  ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  bx = pylab.subplot(142, sharex=ax, sharey=ax)
  bx.plot(profiles[0,1,:],yarray, color='green', linestyle="dashed")
  bx.plot(profiles[1,1,:],yarray, color='blue', linestyle="dashed")
  bx.plot(profiles[2,1,:],yarray, color='red', linestyle="dashed")
  pylab.legend(("5secs", "10secs", "50secs"), loc="lower right")
  bx.set_title('(b) x/h = 2', fontsize=16)
  bx.grid("True")
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(143, sharex=ax, sharey=ax)
  cx.plot(profiles[0,2,:],yarray, color='green', linestyle="dashed")
  cx.plot(profiles[1,2,:],yarray, color='blue', linestyle="dashed")
  cx.plot(profiles[2,2,:],yarray, color='red', linestyle="dashed")
  pylab.legend(("5secs", "10secs", "50secs"), loc="lower right")
  cx.set_title('(c) x/h = 4', fontsize=16)
  cx.grid("True")
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(144, sharex=ax, sharey=ax)
  dx.plot(profiles[0,3,:],yarray, color='green', linestyle="dashed")
  dx.plot(profiles[1,3,:],yarray, color='blue', linestyle="dashed")
  dx.plot(profiles[2,3,:],yarray, color='red', linestyle="dashed")
  pylab.legend(("5secs", "10secs", "50secs"), loc="lower right")
  dx.set_title('(d) x/h = 6', fontsize=16)
  dx.grid("True")
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  pylab.axis([-0.2, 1., 0., 2.])
  bx.set_xlabel('Normalised U-velocity (U/Umax)', fontsize=24)
  ax.set_ylabel('z/h', fontsize=24)

  pylab.savefig("velo_profiles_2d.pdf")
  return

#########################################################################

##### Call reattachment_length function defined above
filelist = get_filelist()

reattachment_length = numpy.array(reatt_length(filelist, exclude_initial_results=130))
av_length = sum(reattachment_length[:,0]) / len(reattachment_length[:,0])
numpy.save("reattachment_length_2D", reattachment_length)
print "\nTime-averaged reattachment length (in step heights): ", av_length
plot_length(reattachment_length)

##### Points to generate profiles:
xarray = numpy.array([0.0, 2.0, 4.0, 6.0])
yarray = numpy.array([0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])

##### Call meanvelo function defined above
profiles = meanvelo(filelist, xarray, yarray)
numpy.save("velo_profiles_2d", profiles)
plot_meanvelo(profiles)
pylab.show()

print "\nAll done.\n"

