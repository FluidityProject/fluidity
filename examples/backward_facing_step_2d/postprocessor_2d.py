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

        ##### reset time to something big to prevent infinite loop
        t = 100.

  print "\n...Finished writing data files.\n"
  return u

#########################################################################

##### Call functions defined above
filelist = get_filelist()
reattachment_length = reatt_length(filelist)
print "\naverage reattachment length (in step heights): ", reattachment_length

##### Points to generate profiles:
xarray = numpy.array([0.0, 2.0, 4.0, 6.0])
yarray = numpy.array([0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])

profiles = meanvelo(filelist, xarray, yarray)
numpy.save("velo_profiles", profiles)

pylab.figure(num=1, figsize = (15., 8.))
pylab.suptitle('U-velocity profiles downstream of backward facing step')

i = 1
for x in range(len(xarray)):
  sub = 140 + i
  pylab.subplot(sub)
  pylab.plot(profiles[x,:],yarray, color='black', linestyle="dashed")
  pylab.axis([-0.2, 1., 0., 2.])
  pylab.grid("True")
  pylab.xlabel('Normalised U-velocity (U/Umax)')
  pylab.ylabel('z/h')
  pylab.title('x/h = '+str(int(xarray[x])))
  i += 1

pylab.savefig("velo_profiles.pdf")
pylab.show()

print "\nAll done.\n"

