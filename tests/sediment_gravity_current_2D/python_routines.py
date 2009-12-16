# CalculatePts(pts_info):
# calculates a set of points from values given in pts_info

################################################################################################

# ExtractProfile(data,pts,field_name):
# extracts values of field_name at the points pts

################################################################################################

# Rotate3DField(field, theta, normal, f_log):
# Rotate a 3D field by angle theta (essentially projection
# Can only rotate around y or z axis
# For rotation about y-axis set normal = [0,1,0], z-axis set normal = [0,0,1]

################################################################################################

# GetIndex(array, val_compare, gorl, f_log):
# find array index when array-val_compare > or < -1.0e-15 (rep 0) 
# having -1.0e-15 is like doing >= or <= I think
# If want less than then gorl = 'l', if want greater than 
# gorl = any string except 'l'

################################################################################################

# SimpsonsQuadHannah(f, dx, f_log):
# Simpsons rule
# Integral(f) approx= (1/3)*dx*(f_0 + f_2N + 4*(f_1+f_3+..+f_2N-1) + 2*(f_2+f_4+...+f_2N-2))
# where N is even hence len(f) needs to be odd 
# if len(f) is even then for f_0-f_2N-1 use Simpsons rule and calculate the last section
# using the trapezoidal rule area = (1/2)*dx*(f_2N-1+f_2N)
# see either Riley, Hobson and Bence pg 1167 or http://mathworld.wolfram.com/SimpsonsRule.html

################################################################################################

# LineIntegral(f, dp, index_1, index_2, p_1, p_2, f_log):
# integrates in one dimension
# index_1 and index_2 are the indicies just after the values to be integrated to e.g. if want to integrate
# from c = 3,6 and c = [2.1111,3.111,4.111,5.111,6.111] then index_1 = 1, index_2 = 4
# p_1 and p_2 are the coordinate or function limits to be integrated over
# c are the coordinates that are beinng integrated over. These only matter if p_1, p_2 are the coordinate limits 
# rather than the function limits. If the latter then just set c = []

################################################################################################

import vtk
import glob
import sys
import os
import scipy.stats
import math
import vtktools
import numpy 
from numpy import arange
import python_routines

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def SetConstants(f_log):

  f_log.write('Entered gc_setsimulationspecifics.SetConstants')

  constants = [[],[]]

# 0
# theta  
  theta = 2.5
  constants[0].append('theta plus')
  constants[1].append(math.radians(theta))

# 1
# theta minus
  constants[0].append('theta minus')
  constants[1].append(math.radians(-1.0*theta))

# 2
# xplateau
  constants[0].append('xplateau')
  constants[1].append(19000.0)

# 3
# zplateau
  constants[0].append('hplateau')
  constants[1].append(800.0)

# 4
# contour value to track head and tail
  constants[0].append('contour_val')
  constants[1].append([0.0,1.0,2.0,3.0])
  
# 5
# min y value
  constants[0].append('y min')
  constants[1].append(0.0)

# 6
# max y value
  constants[0].append('y max')
  constants[1].append(10.0)

# 7
# number of points to perform regression over for U_F and dH/dx
  constants[0].append('no regression pts')
  constants[1].append(5)

# 8
# no. of dimensions
  constants[0].append('no. of dimensions')
  constants[1].append(2)

# 9
# beta
  constants[0].append('beta')
  constants[1].append(7.5E-4)

# 10
# gravity
  constants[0].append('gravity')
  constants[1].append(9.81)

# 11
# dS, source
  constants[0].append('dS source')
  constants[1].append(3.0)

  f_log.write('Exiting gc_setsimulationspecifics.SetConstants')

  return constants 

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def CalculatePts(pts_info, f_log):
# calculates a set of points from values given in pts_info
# usually for use in Extract_profile

#  x_min = pts_info[0]
#  x_max = pts_info[1]
#  no x steps = pts_info[2]
#  y_min = pts_info[3]
#  y_max = pts_info[4]
#  no y steps = pts_info[5]
#  z_min = pts_info[6]
#  z_max = pts_info[7]
#  no z steps = pts_info[8]

  f_log.write('\n Entered gc_extracttools.CalculatePts')

# get theta and co-ordinates of edge of plateau and y halfway through domain
  constants = SetConstants(f_log)
  theta = constants[1][0]
  xplateau = constants[1][2]
  hplateau = constants[1][3]

#  bbox = data.ugrid.GetBounds()
#  y = (bbox[3]+bbox[2])*0.5

  step_x = (pts_info[1]-pts_info[0])/(pts_info[2])
  step_y = (pts_info[4]-pts_info[3])/(pts_info[5])
  step_z = (pts_info[7]-pts_info[6])/(pts_info[8])

  step_x_dummy = step_x
  step_y_dummy = step_y
  step_z_dummy = step_z
  if(step_x_dummy == 0): step_x_dummy = 1.0; step_x = 0.5
  if(step_y_dummy == 0): step_y_dummy = 1.0; step_y = 0.5
  if(step_z_dummy == 0): step_z_dummy = 1.0; step_z = 0.5

  pts = []
  if(pts_info[9] == 'yes'):
    #print 'rotate'
    s_angle = math.sin(theta)
    c_angle = math.cos(theta)
    if(constants[1][8] == 2):
        s_angle_y = s_angle
        c_angle_y = c_angle
        s_angle_z = math.sin(0.0)
        c_angle_z = math.cos(0.0) 
        yplateau = hplateau
        zplateau = 0.0
    else:
        s_angle_y = math.sin(0.0)
        c_angle_y = math.cos(0.0) 
        s_angle_z = s_angle
        c_angle_z = c_angle
        yplateau = 0.0
        zplateau = hplateau
    for x in arange(pts_info[0], pts_info[1]+step_x, step_x_dummy):
      pts_subset = [] 
      x_0 = xplateau - (x*c_angle)
      y_0 = yplateau - (x*s_angle_y)
      z_0 = zplateau - (x*s_angle_z)
      for y in arange(pts_info[3], pts_info[4]+step_y, step_y_dummy):
          pts_subsubset = []
          for z in arange(pts_info[6], pts_info[7]+step_z, step_z_dummy):
              x = x_0 - (y*s_angle_y) - (z*s_angle_z)
              y = y_0 + (y*c_angle_y)
              z = z_0 + (z*c_angle_z)
              pts_subsubset.append([x,y,z])
          pts_subset.append(pts_subsubset)
      pts.append(pts_subset)
  
  else:
    #print 'no rotate'
    for x in arange(pts_info[0], pts_info[1]+step_x, step_x_dummy):
      pts_subset = [] 
      for y in arange(pts_info[3], pts_info[4]+step_y, step_y_dummy):
          pts_subsubset = [] 
          for z in arange(pts_info[6], pts_info[7]+step_z, step_z_dummy):
              pts_subsubset.append([x,y,z])
          pts_subset.append(pts_subsubset)
      pts.append(pts_subset)

  f_log.write('\n Exiting gc_extracttools.CalculatePts')

  return pts

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def ExtractProfile(data,pts,fieldname, f_log):
# extracts values of fieldname at the points pts

  f_log.write('\n Entered gc_extracttools.ExtractProfile')

  profile = []
#  print len(pts)
  for i in range(len(pts)):
      profile_subset = []
      for j in range(len(pts[i])):
          field = data.ProbeData(vtktools.arr(pts[i][j]), fieldname)
          #print i, j
          #print field
          profile_subset.append(field) 
      profile.append(profile_subset)

  f_log.write('\n Exiting gc_extracttools.ExtractProfile') 
 
  return profile

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def Rotate3DField(field, theta, normal, f_log):
# Rotate a 3D field by angle theta (essentially projection)
# Can only rotate around y or z axis
# For rotation about y-axis set normal = [0,1,0], z-axis set normal = [0,0,1]

  f_log.write('\n Entered gc_diagnostictools.Rotate3DField')

  field_out = [[[[] for k in range(len(field[i][j]))] for j in range(len(field[i]))] for i in range(len(field))]

  for i in range(len(field)):
    for j in range(len(field[i])):
      for k in range(len(field[i][j])):
        if(normal != [0,0,1]): 
          if(normal != [0,1,0]): print "Only coded to rotate around y or z axis, exiting"; f_log.write("Only coded to rotate around y or z axis, exiting"); sys.exit(1)

        c_angle = math.cos(theta)
        s_angle = math.sin(theta)
        
        field_out[i][j][k].append((c_angle*field[i][j][k][0]) + s_angle*((normal[2]*field[i][j][k][1]) + (normal[1]*field[i][j][k][2]))) 
        field_out[i][j][k].append((normal[1]*field[i][j][k][1]) - (normal[2]*s_angle*field[i][j][k][0]) + (normal[2]*c_angle*field[i][j][k][1])) 
        field_out[i][j][k].append((normal[2]*field[i][j][k][2]) - (normal[1]*s_angle*field[i][j][k][0]) + (normal[1]*c_angle*field[i][j][k][2]))  
        #print field[i][j][k][0], field[i][j][k][1], field[i][j][k][2]
        #print field_out[i][j][k][0], field_out[i][j][k][1], field_out[i][j][k][2]

  f_log.write('\n Exiting gc_diagnostictools.Rotate3DField')

  return field_out

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def GetIndex(array, val_compare, gorl, f_log):
# find first array index when array-val_compare > or < -1.0e-15 (rep 0)
# having -1.0e-15 is like doing >= or <= I think
# If want less than then gorl = 'l', if want greater than 
# gorl = any string except 'l'
# Note integration will go up to index-1 and then do some interpolation
# and use the trapezium rule to integrate between index -1 and index
# or it will go from index and then do some interpolation
# and use the trapezium rule to integrate between index -1 and index 

  f_log.write('\n Entered gc_diagnostictools.GetIndex')

# first check to see if index wanted is zero
  index = 0
  inequality = 1.0
  if(gorl == 'l'): inequality = -1.0 
  if(inequality*(array[0]-val_compare)> -1.0e-15): index = 0
  else:
    for i in range(len(array)):
      if(index==0):
        if(inequality*(array[i]-val_compare)> -1.0e-15): index = i 

  f_log.write('\n Exiting gc_diagnostictools.GetIndex')
  
  return index

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def SimpsonsQuadHannah(f, dx, f_log):
# Simpsons rule
# Integral(f) approx= (1/3)*dx*(f_0 + f_2N + 4*(f_1+f_3+..+f_2N-1) + 2*(f_2+f_4+...+f_2N-2))
# where N is even hence len(f) needs to be odd 
# if len(f) is even then for f_0-f_2N-1 use Simpsons rule and calculate the last section
# using the trapezoidal rule area = (1/2)*dx*(f_2N-1+f_2N)
# see either Riley, Hobson and Bence pg 1167 or http://mathworld.wolfram.com/SimpsonsRule.html
# error = (x2-x1)(dx^4)(d^4f/dx^4)(1/180) = N(dx^5)(d^4f/dx^4)(1/90)

  f_log.write('\n Entered gc_diagnostictools.SimpsonsQuadHannah')

  integral = 0.0
  
  if(len(f) == 0): integral = 0.0 ; f_log.write("WARNING length input array to SimpsonsQuadHannah is zero, returning zero")
  elif(len(f) == 1): integral = 0.0 ; f_log.write("WARNING length input array to SimpsonsQuadHannah is 1, returning zero")
  elif(len(f) == 2): integral = 0.5*dx*(f[0]+f[1]) ; f_log.write("WARNING length input array to SimpsonsQuadHannah is 2, using trapezoidal rule")
  else: 
    odd = bool(len(f)%2)
    if (odd == False):
      length = len(f)-1
      f_log.write("len(f) even, last area calculated using trapezoidal rule")
    else:
      length = len(f)
  
# add the 0 and 2N terms
    integral = (f[0]+f[length-1])
# add the 1 - 2N-3 and 2 - 2N-2 terms
    for step in range(1,((length-1)/2)):
      integral = integral + (4*f[(2*step)-1])+(2*f[2*step])
# add in the 2N-1 term
    integral = integral + 4*f[length-2]
    #print integral
    integral = (1.0/3.0)*dx*integral

    if(odd == False):
      integral = integral + 0.5*dx*(f[length-1]+f[length])

  f_log.write('\n Exiting gc_diagnostictools.SimpsonsQuadHannah')
    
  return integral

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################

def LineIntegral(f, c, dp, index_1, index_2, p_1, p_2, f_log):
# integrates in one dimension
# index_1 and index_2 are the indicies just after the values to be integrated to e.g. if want to integrate
# from c = 3,6 and c = [2.1111,3.111,4.111,5.111,6.111] then index_1 = 1, index_2 = 4
# p_1 and p_2 are the coordinate or function limits to be integrated over
# c are the coordinates that are beinng integrated over. These only matter if p_1, p_2 are the coordinate limits 
# rather than the function limits. If the latter then just set c = []

  f_log.write('\n Entered gc_diagnostictools.LineIntegral')

  integral = SimpsonsQuadHannah(f[index_1:index_2], dp, f_log)

# do some linear interpolation and an extra bit of quadrature 
# at bottom and top to account for anything missed
  if(index_1 != 0): 
    m = (f[index_1]-f[index_1-1])/dp
    if(c != []): 
      dp_new = c[index_1]-p_1
      f_lower = f[index_1]-(m*dp_new)
      integral = integral + (0.5*dp_new*(f_lower+f[index_1]))
    else:
      if(m == 0): dp_new = 0
      else: dp_new = abs((f[index_1]-p_1)/m)
      integral = integral + (0.5*dp_new*(p_1+f[index_1]))
#    f_middle = 0.5*(f[index_1 - 1] + f[index_1])
#    integral = integral + (0.5*dp*(f_middle+f[index_1]))
  #print index_2
  #print len(f)
  if(index_2 != len(f)-1): 
    m = (f[index_2]-f[index_2-1])/dp
    if(c != []):
      dp_new = p_2 - c[index_2-1]
      f_upper = f[index_2-1]+ m*dp_new
      integral = integral + (0.5*dp_new*(f[index_2-1]+f_upper))
    else:
      if(m == 0): dp_new = 0
      else: dp_new = abs((p_2-f[index_2-1])/m)
      #print f[index_2], f[index_2-1], m, dp
      #print 'dp_new = ', dp_new
      integral = integral + (0.5*dp_new*(f[index_2-1]+p_2))
#    f_middle = 0.5*(f[index_2] + f[index_2+1])
#    integral = integral + (0.5*0.5*dp*(f_middle+f[index_2]))
#    N.B. this should probably have had index_2 replaced by index_2-1 

  f_log.write('\n Exiting gc_diagnostictools.LineIntegral')

  return integral   

################################################################################################
#----------------------------------------------------------------------------------------------#
################################################################################################



