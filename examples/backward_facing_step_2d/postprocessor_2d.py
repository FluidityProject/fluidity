#!/usr/bin/env python

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
def reatt_length(filelist, exclude_initial_results=130):

  print "Running reattachment length script on selected vtu files...\n"
  nums=[]
  result=0.
  filecount=0

  # check for no files
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
    if num > exclude_initial_results:
      nums.append(num)
    nums.sort()

  for num in nums:
    ##### counter for averaging
    filecount = filecount + 1
    file = "backward_facing_step_2d_"+str(num)+".vtu"

    print file
    ##### Read in data from vtu
    reader = vtk.vtkXMLUnstructuredGridReader();
    reader.SetFileName(file)
    reader.Update()
    data = reader.GetOutput()
    datafile = vtktools.vtu(file)

    ##### specify points near bottom surface
    pts=[]; no_pts = 100; x_0 = 5.0; offset = 0.1
    pts.append((5.0, offset, 0.0))
    for i in range(1, no_pts):
      x = x_0 + i/5.0
      pts.append((x, offset, 0.0))

    ##### Get x-velocity on bottom boundary
    uvw = datafile.ProbeData(vtktools.arr(pts), "Velocity")
    u = []; surf = []
    for i in range(len(uvw)):
      u.append(uvw[i][0]); surf.append(pts[i][0])

    ##### Find all potential reattachment points:
    #print "calculating reattachment point locations using change of x-velocity sign"
    points = []
    for i in range(len(u)-1):
      if((u[i] / u[i+1]) < 0.):
        ##### interpolate between nodes
        p = surf[i] + (surf[i+1]-surf[i]) * (0.0-u[i]) / (u[i+1]-u[i]) - 5.0
        #print "potential reattachment point at x = ", p
        points.append(p)

    ##### Remove corner point, if present:
    if(len(points) > 1):
      if(points[0] < 0.1):
        print "ignoring point at step corner"
        points.remove(points[0])

    ##### This is the actual reattachment point:
    result = result + points[0]

  ##### Average over selected files:
  result = result / filecount
  return result

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
        t = 100.
        filecount += 1

  print "\n...Finished writing data files.\n"
  return profiles

#########################################################################

##### Call functions defined above
filelist = get_filelist()
reattachment_length = reatt_length(filelist)
print "\naverage reattachment length (in step heights): ", reattachment_length

##### Points to generate profiles:
xarray = numpy.array([0.0, 2.0, 4.0, 6.0])
yarray = numpy.array([0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])

profiles = meanvelo(filelist, xarray, yarray)
numpy.save("velo_profiles_2d", profiles)

#for t in range(profiles[:,0,0]):
pylab.figure(num=1, figsize = (16.5, 8.5))
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
pylab.show()

print "\nAll done.\n"

