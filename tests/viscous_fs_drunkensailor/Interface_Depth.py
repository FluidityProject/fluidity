#!/usr/bin/env python

import vtk
import vtktools
import glob
import sys
import os
import scipy.stats
import operator
import numpy
import pylab

def key(file):
  return int(file.split('_')[-1].split('.')[0])

def get_filelist():
  filelist = glob.glob1('.', 'viscous_fs_drunkensailor_*.vtu')
  filelist = sorted(filelist, key=key)
  return filelist

def get_time(file):
  vtufile = vtktools.vtu(file)
  time    = vtufile.GetScalarRange("Dense::Time")
  return time[0]

def get_interface_depth(file):
  reader = vtk.vtkXMLUnstructuredGridReader();
  reader.SetFileName(file)
  reader.Update()
  data = reader.GetOutput()
  data.GetPointData().SetActiveScalars("Dense::MaterialVolumeFraction")
  contour = vtk.vtkContourFilter ()
  contour.SetInput(data)
  contour.SetValue(0, 0.5)
  contour.Update()
  polydata = contour.GetOutput()
  bounding_box = polydata.GetBounds()
  interface_depth = bounding_box[2]
  return interface_depth

# Get list of files:
files = get_filelist()

# Get interface depth and simulation time for all vtu files:
interface_depth=numpy.zeros(len(files))
time=numpy.zeros(len(files))
myr2sec=60.*60.*24.*365.*1e6

for file in range(len(files)) :
  interface_depth[file] = (get_interface_depth(files[file]) / 1e3)
  time[file] = (get_time(files[file]) / myr2sec)
  print interface_depth[file], time[file]

pylab.plot(time,interface_depth)
pylab.xlabel("Time (Myr)")
pylab.ylabel("Interface Depth (km)")
pylab.xlim(0.0,6)
pylab.ylim(-100,-500)
ax = pylab.gca()
ax.set_ylim(ax.get_ylim()[::-1])

pylab.show()
