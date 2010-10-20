#!/usr/bin/env python

import vtk
import glob
import sys
import os
import vtktools
import numpy

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

#########################################################################

# Reattachment length:
def reatt_length(filelist, exclude_initial_results=120):

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
    points = []
    for i in range(len(u)-1):
      if((u[i] / u[i+1]) < 0.):
        ##### interpolate between nodes
        p = surf[i] + (surf[i+1]-surf[i]) * (0.0-u[i]) / (u[i+1]-u[i]) - 5.0
        points.append(p)

    ##### Remove corner point, if present:
    if(len(points) > 1):
      if(points[0] < 0.1):
        print "ignoring point at step corner"
        points.remove(points[0])

    ##### Sum lengths:
    result = result + points[0]

  ##### Average over selected files:
  result = result / filecount
  return result

#########################################################################

# Call functions defined above
filelist = get_filelist()
reattachment_length = reatt_length(filelist)
print "\nTime-averaged reattachment length (in step heights): ", reattachment_length
