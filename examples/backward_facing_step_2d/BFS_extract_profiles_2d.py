#!/usr/bin/env python

import vtk
import glob
import sys
import os
import vtktools
import scipy.optimize
import numpy

def get_filelist():

    def key(s):
        return int(s.split('_')[-1].split('.')[0])
   
    list = glob.glob("*.vtu")
    #print 'list before sorting'
    #print list
    list = [l for l in list if 'check' not in l]
    vtu_nos = [float(s.split('_')[-1].split('.')[0]) for s in list]
    vals = zip(vtu_nos, list)
    vals.sort()
    unzip = lambda l:tuple(apply(zip,l))
    vtu_nos, list = unzip(vals)

    return list

###################################################################

def meanvelo(filelist):

  ##### check for no files
  if (len(filelist) < 0):
    print "No files!"
    sys.exit(1)
  for file in filelist:
      try:
        os.stat(file)
      except:
        f_log.write("No such file: %s" % files)
        sys.exit(1)

      data = vtktools.vtu(file)    

      print "Extracting 2d BFS data:", file
      reader = vtk.vtkXMLUnstructuredGridReader();
      reader.SetFileName(file)
      reader.Update()
      data = reader.GetOutput()
      datafile = vtktools.vtu(file)

      ##### Create array of points
      x = [4.0, 6.0, 10.0, 19.0]
      z = [0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
      no_pts = len(x)*len(z)
      pts=[]

      for i in range(len(x)):
        for j in range(len(z)):
          pts.append([x[i], 0.0, z[j]])

      print "points: ", pts

      ##### Get x-velocity
      uvw = datafile.ProbeData(vtktools.arr(pts), "Velocity")
      u = []
      for i in range(len(uvw)):
        u.append(uvw[i][0])

      #print "x-velocity: ", u


      print 'write output data file'
      num = int(file.split(".vtu")[0].split('_')[-1])

      nf = open('BFS2d_velox_profiles_'+str(num)+'.dat','w')
      nf.write('x''    ''z''    ''U''\n')

      for n in range(len(x)*len(z)):
        nf.write(str(pts[n][0])+' '+str(pts[n][2])+' '+str(u[n])+'\n')

      results = "All done."
  return results

#########################################################################

filelist = get_filelist()
print "filelist: ", filelist
umean = meanvelo(filelist)
