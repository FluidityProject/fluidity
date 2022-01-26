#!/usr/bin/env python3

import vtk
import vtktools
import glob
import numpy

def key(file):
  return int(file.split('_')[-1].split('.')[0])

def get_filelist():
  filelist = glob.glob1('.', 'viscous_fs_drunkensailor_*.pvtu')
  filelist = sorted(filelist, key=key)
  return filelist

def get_time(file):
  vtufile = vtktools.vtu(file)
  time    = vtufile.GetScalarRange("Dense::Time")
  return time[0]

def get_interface_depth(file):
  vtu=vtktools.vtu(file)
  data = vtu.ugrid
  data.GetPointData().SetActiveScalars("Dense::MaterialVolumeFraction")
  contour = vtk.vtkContourFilter ()
  if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
    contour.SetInput(data)
  else:
    contour.SetInputData(data)
  contour.SetValue(0, 0.5)
  contour.Update()
  polydata = contour.GetOutput()
  bounding_box = polydata.GetBounds()
  interface_depth = bounding_box[2]
  return interface_depth

def find_interface_depth():
  # Get list of files:
  files = get_filelist()

  # Get interface depth and simulation time for all vtu files:
  interface_depth = numpy.zeros(len(files))
  time            = numpy.zeros(len(files))
  myr2sec         = 60.*60.*24.*365.*1e6
  rescale_time    = (500e3 / 1e-9)

  for file in range(len(files)) :
    interface_depth[file] = (get_interface_depth(files[file]) * 500)
    time[file] = ((get_time(files[file]) * rescale_time) / myr2sec)

  final_interface_depth = interface_depth[-1]
  return final_interface_depth
