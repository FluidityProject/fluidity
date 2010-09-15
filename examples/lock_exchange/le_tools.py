import vtk
import glob
import sys
import os
import scipy.stats
import math
import vtktools
import numpy
import vtktools
from numpy import arange
import pylab

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def GetFiles(directory):

# gets list of (p)vtus and sorts them into ascending time order

  def key(s):
    return int(s.split('_')[-1].split('.')[0])
 
  list = glob.glob(directory+"*.pvtu")
  if len(list) == 0: list = glob.glob(directory+"*.vtu")
  list = [l for l in list if 'check' not in l]
  list = [l for l in list if 'Mesh' not in l]

  return sorted(list, key=key)

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def GetXandt(filelist):
  time = []
  X_ns = []
  X_fs = []
  for files in filelist:

      data = vtktools.vtu(files) 
      
      time.append(data.GetScalarField("Time")[0])
      
      # Get X
      data.ugrid.GetPointData().SetActiveScalars('Temperature')
      data = data.ugrid
      
      contour = vtk.vtkContourFilter()
      contour.SetInput(data)
      contour.SetValue(0, 0.0)
      contour.Update()
      polydata = contour.GetOutput()

      bounding_box = polydata.GetBounds()
   
      X_ns.append(bounding_box[1])
      X_fs.append(bounding_box[0])

  return time, X_ns, X_fs

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def GetU(t, X):
#
  U = [(X[1]-X[0])/(t[1]-t[0])]
#  
  for i in range(1,len(X)-1):
    U.append(LeastSquares(t[i-1:i+2],X[i-1:i+2],[0.0,1.0])[0][0][1])
#
  U.append((X[-1]-X[-2])/(t[-1]-t[-2]))
#
  return U
  
################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def LeastSquares(x_values,y_values,p0):
# calculate a least squares approximation to the data (x_values,y_values)
# starts from a guess at the initial solution given in 'predictedform'
# p0 should be an array containing the initial guess of the coefficients in 'predictedform'
# plsq returns the coefficients of p for the best fit
  
  plsq = scipy.optimize.leastsq(residuals, numpy.array(p0), args=(numpy.array(x_values), numpy.array(y_values)))

  lsq  = []
  for val in x_values: lsq.append(predictedform(plsq[0],val))
  
  return (plsq,lsq)


def residuals(p, x, y):
  err = y-predictedform(p,x)
  return err

    
def predictedform(p,x):
  pf = 0
  for i in range(len(p)): pf = pf+p[i]*(x**i)
  return pf

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def GetAverageRange(X, lower_lim, domainheight):
#
# get range of X for averaging such that lower_lim<X<h0
  try: start_val = pylab.find(numpy.array(X)>lower_lim)[0]
  except IndexError: print "not enough values of X, \n lower lim = " + str(lower_lim) + "\n X = " + str(X); sys.exit(1)
  try: end_val = pylab.find(numpy.array(X)>0.4-domainheight)[0]
  except IndexError: end_val = len(X)
#
  return start_val, end_val

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def ReadLog(log_name):
# reads log with name log_name and returns values as an array of floats
# log needs to contain only one value per line
  file_name = open(log_name,'r')
  file_read = []
  for line in file_name:
    l = line
    file_read.append(float(l.split("\n")[0]))
  file_name.close()
  return file_read  
  
################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################  
  
  
  
  
