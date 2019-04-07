from __future__ import print_function

import os
import math
import glob
import numpy
import vtktools
from fluidity_tools import stat_parser
from sympy import *
from numpy import array,max,abs

meshtemplate='''
Point (1) = {0.0, 0.0, 0, 1.0};
Point (2) = {0.0, 1.0, 0, 1.0};
Line (1) = {1, 2};

Extrude {1,0,0} {
  Line{1};
}

// Volume number for whole domain.
Physical Surface (1) = {5};
// Top of the box.
Physical Line(333) = {4};
// Bottom of the box.
Physical Line(444) = {3};
// Side walls.
Physical Line(666) = {1,2};
Field[1] = MathEval;
Background Field = 1;
Field[1].F = "1./<NN>";
'''


def generate_meshfile(name,NN):

    file(name+".geo",'w').write(
        meshtemplate.replace('<NN>',str( NN)))

    os.system("gmsh -2 "+name+".geo")

    os.system("../../../bin/gmsh2triangle --2d "+name+".msh")


def plot_results(NN, error):
    '''plot_results(error)

    Produce a plot of the actual errors provided in the argument
    "error". Error is a matrix with eight columns, one for each of the 
    error metrics computed. 
    '''
    from pylab import \
    figure,xticks,yticks,axis,xlabel,ylabel,loglog,legend,title,savefig

    figure()
    dx = 1./NN
    loglog(dx,error,'o-')
    loglog(dx,100*dx**2)
    yticks(yticks()[0], map(lambda x: "%3.1e"%x, yticks()[0]))
    xticks(xticks()[0], map(lambda x: "%3.1e"%x, xticks()[0]))
    xlabel('mesh spacing')
    ylabel('RMS error')

    title("Convergence in the driven cavity test problem")

    legend(("error1","error2","error3","error4","error5","error6","error7","error8","2nd order"))

    savefig('driven_cavity_error_plot.png')


def retrieve_results(NN):
    '''retrieve_results(NN)

    For each NN in the input sequence, retrieve the various errors
    from the simulation results in appropriate driven_cavity-NN
    directory.

    The columns of the results are the: Erturk et al 2005 u and v errors, 
    the Botella et al 1998 u, v, and p errors and the Brunean et al 2006 kinetic energy and
    streamfunction errors. 
    The errors are the RMS difference between highly accurate
    tabulated values available in the papers: 

    E. Erturk, T. C. Corke and C. Gokcol C, Numerical Solutions of 2-D Steady 
    Incompressible Driven Cavity Flow at High Reynolds Numbers, International 
    Journal for Numerical Methods in Fluids 48, 747-774, 2005. doi:10.1002/fld.953

    O. Botella and R. Peyret, Benchmark spectral results on the lid-driven cavity flow,
    Computers & Fluids 27, 421-433, 1998, doi:10.1016/S0045-7930(98)00002-4 

    C.-H. Bruneau and M. Saad, The 2D lid-driven cavity problem revisited,
    Computers & Fluids 35, 326-348, 2006, doi:10.1016/j.compfluid.2004.12.004 

    '''
    from numpy import zeros

    error=zeros((len(NN),8))

    for i,NN in enumerate(NN):

        error[i,0]=erturk_u(NN)
        error[i,1]=erturk_v(NN)
        error[i,2]=botella_u(NN)
        error[i,3]=botella_v(NN)
        error[i,4]=botella_p1(NN)
        error[i,5]=botella_p2(NN)
        error[i,6]=bruneau_ke(NN)
        error[i,7]=bruneau_sf(NN)
    
    return error


def erturk_u(NN):
#Erturk et al 2005. Table VI
  filelist_not_sorted = glob.glob('driven_cavity-%d/*.vtu'%NN)
  vtu_nos_not_sorted = [int(file.split('.vtu')[0].split('_')[-1]) for file in filelist_not_sorted]
  filelist = [filelist_not_sorted[i] for i in numpy.argsort(vtu_nos_not_sorted)]
  file = filelist[-1]
  print(file)
  try:
    os.stat(file)
  except:
    print("No such file: %s" % file)
    sys.exit(1)

  u=vtktools.vtu(file)
  pts=vtktools.arr([[0.5, 0.000, 0.000, 0.0000],
  [0.5, 0.000, 0.000,  0.0000],
  [0.5, 0.020, 0.000, -0.0757],
  [0.5, 0.040, 0.000, -0.1392],
  [0.5, 0.060, 0.000, -0.1951],
  [0.5, 0.080, 0.000, -0.2472],
  [0.5, 0.100, 0.000, -0.2960],
  [0.5, 0.120, 0.000, -0.3381],
  [0.5, 0.140, 0.000, -0.3690],
  [0.5, 0.160, 0.000, -0.3854],
  [0.5, 0.180, 0.000, -0.3869],
  [0.5, 0.200, 0.000, -0.3756],
  [0.5, 0.500, 0.000, -0.0620],
  [0.5, 0.900, 0.000,  0.3838],
  [0.5, 0.910, 0.000,  0.3913],
  [0.5, 0.920, 0.000,  0.3993],
  [0.5, 0.930, 0.000,  0.4101],
  [0.5, 0.940, 0.000,  0.4276],
  [0.5, 0.950, 0.000,  0.4582],
  [0.5, 0.960, 0.000,  0.5102],
  [0.5, 0.970, 0.000,  0.5917],
  [0.5, 0.980, 0.000,  0.7065],
  [0.5, 0.990, 0.000,  0.8486],
  [0.5, 1.000, 0.000,  1.0000]])
  
  velocity = u.ProbeData(pts, "Velocity")
  (ilen, jlen) = velocity.shape
  norm=0.0
  for i in range(ilen):
      diff = pts[i][3] - velocity[i][0]
      norm = norm + diff*diff
  
  norm = math.sqrt(norm/ilen)
  print("erturk_u_norm:", norm)

  return norm


def erturk_v(NN):
#Erturk et al 2005. Table VII
  filelist_not_sorted = glob.glob('driven_cavity-%d/*.vtu'%NN)
  vtu_nos_not_sorted = [int(file.split('.vtu')[0].split('_')[-1]) for file in filelist_not_sorted]
  filelist = [filelist_not_sorted[i] for i in numpy.argsort(vtu_nos_not_sorted)]
  file = filelist[-1]
  print(file)
  try:
    os.stat(file)
  except:
    print("No such file: %s" % file)
    sys.exit(1)

  u=vtktools.vtu(file)
  pts=vtktools.arr([
  [0.000, 0.5, 0.0,  0.0000],
  [0.015, 0.5, 0.0,  0.1019],
  [0.030, 0.5, 0.0,  0.1792],
  [0.045, 0.5, 0.0,  0.2349],
  [0.060, 0.5, 0.0,  0.2746],
  [0.075, 0.5, 0.0,  0.3041],
  [0.090, 0.5, 0.0,  0.3273],
  [0.105, 0.5, 0.0,  0.3460],
  [0.120, 0.5, 0.0,  0.3605],
  [0.135, 0.5, 0.0,  0.3705],
  [0.150, 0.5, 0.0,  0.3756],
  [0.500, 0.5, 0.0,  0.0258],
  [0.850, 0.5, 0.0, -0.4028],
  [0.865, 0.5, 0.0, -0.4407],
  [0.880, 0.5, 0.0, -0.4803],
  [0.895, 0.5, 0.0, -0.5132],
  [0.910, 0.5, 0.0, -0.5263],
  [0.925, 0.5, 0.0, -0.5052],
  [0.940, 0.5, 0.0, -0.4417],
  [0.955, 0.5, 0.0, -0.3400],
  [0.970, 0.5, 0.0, -0.2173],
  [0.985, 0.5, 0.0, -0.0973],
  [1.000, 0.5, 0.0,  0.0000]])
  
  velocity = u.ProbeData(pts, "Velocity")
  (ilen, jlen) = velocity.shape
  norm=0.0
  for i in range(ilen):
      diff = pts[i][3] - velocity[i][1]
      norm = norm + diff*diff
  
  norm = math.sqrt(norm/ilen)
  print("erturk_v_norm:", norm)

  return norm

def botella_u(NN):
#Botella and Peyret (1998) Table 9. NB.our velocity at the lid is reverse theirs therefore minus signs in u below
  filelist_not_sorted = glob.glob('driven_cavity-%d/*.vtu'%NN)
  vtu_nos_not_sorted = [int(file.split('.vtu')[0].split('_')[-1]) for file in filelist_not_sorted]
  filelist = [filelist_not_sorted[i] for i in numpy.argsort(vtu_nos_not_sorted)]
  file = filelist[-1]
  print(file)
  try:
    os.stat(file)
  except:
    print("No such file: %s" % file)
    sys.exit(1)

  u=vtktools.vtu(file)
  pts=vtktools.arr([
  [0.5, 0.0000, 0.0,  0.0000000],
  [0.5, 0.0547, 0.0, -0.1812881],
  [0.5, 0.0625, 0.0, -0.2023300],
  [0.5, 0.0703, 0.0, -0.2228955],
  [0.5, 0.1016, 0.0, -0.3004561],
  [0.5, 0.1719, 0.0, -0.3885691],
  [0.5, 0.2813, 0.0, -0.2803696],
  [0.5, 0.4531, 0.0, -0.1081999],
  [0.5, 0.5000, 0.0, -0.0620561],
  [0.5, 0.6172, 0.0,  0.0570178],
  [0.5, 0.7344, 0.0,  0.1886747],
  [0.5, 0.8516, 0.0,  0.3372212],
  [0.5, 0.9531, 0.0,  0.4723329],
  [0.5, 0.9609, 0.0,  0.5169277],
  [0.5, 0.9688, 0.0,  0.5808359],
  [0.5, 0.9766, 0.0,  0.6644227],
  [0.5, 1.0000, 0.0,  1.0000000]])

  velocity = u.ProbeData(pts, "Velocity")
  (ilen, jlen) = velocity.shape
  norm=0.0
  for i in range(ilen):
      diff = pts[i][3] - velocity[i][0]
      norm = norm + diff*diff

  norm = math.sqrt(norm/ilen)
  print("botella_u_norm:", norm)

  return norm

def botella_v(NN):
#Botella and Peyret (1998) Table 10. 
  filelist_not_sorted = glob.glob('driven_cavity-%d/*.vtu'%NN)
  vtu_nos_not_sorted = [int(file.split('.vtu')[0].split('_')[-1]) for file in filelist_not_sorted]
  filelist = [filelist_not_sorted[i] for i in numpy.argsort(vtu_nos_not_sorted)]
  file = filelist[-1]
  print(file)
  try:
    os.stat(file)
  except:
    print("No such file: %s" % file)
    sys.exit(1)

  u=vtktools.vtu(file)
  pts=vtktools.arr([
  [1.0000, 0.5, 0.0,  0.0000000],
  [0.9688, 0.5, 0.0, -0.2279225],
  [0.9609, 0.5, 0.0, -0.2936869],
  [0.9531, 0.5, 0.0, -0.3553213],
  [0.9453, 0.5, 0.0, -0.4103754],
  [0.9063, 0.5, 0.0, -0.5264392],
  [0.8594, 0.5, 0.0, -0.4264545],
  [0.8047, 0.5, 0.0, -0.3202137],
  [0.5000, 0.5, 0.0,  0.0257995],
  [0.2344, 0.5, 0.0,  0.3253592],
  [0.2266, 0.5, 0.0,  0.3339924],
  [0.1563, 0.5, 0.0,  0.3769189],
  [0.0938, 0.5, 0.0,  0.3330442],
  [0.0781, 0.5, 0.0,  0.3099097],
  [0.0703, 0.5, 0.0,  0.2962703],
  [0.0625, 0.5, 0.0,  0.2807056],
  [0.0000, 0.5, 0.0,  0.0000000]])

  velocity = u.ProbeData(pts, "Velocity")
  (ilen, jlen) = velocity.shape

  norm=0.0
  for i in range(ilen):
      diff = pts[i][3] - velocity[i][1]
      norm = norm + diff*diff

  norm = math.sqrt(norm/ilen)
  print("botella_v_norm:", norm)

  return norm

def botella_p1(NN):
#Botella and Peyret (1998) Table 9. 
  filelist_not_sorted = glob.glob('driven_cavity-%d/*.vtu'%NN)
  vtu_nos_not_sorted = [int(file.split('.vtu')[0].split('_')[-1]) for file in filelist_not_sorted]
  filelist = [filelist_not_sorted[i] for i in numpy.argsort(vtu_nos_not_sorted)]
  file = filelist[-1]
  print(file)
  try:
    os.stat(file)
  except:
    print("No such file: %s" % file)
    sys.exit(1)

  u=vtktools.vtu(file)
  pts=vtktools.arr([
  [0.5, 0.0000, 0.0,  0.110591],
  [0.5, 0.0547, 0.0,  0.109689],
  [0.5, 0.0625, 0.0,  0.109200],
  [0.5, 0.0703, 0.0,  0.108566],
  [0.5, 0.1016, 0.0,  0.104187],
  [0.5, 0.1719, 0.0,  0.081925],
  [0.5, 0.2813, 0.0,  0.040377],
  [0.5, 0.4531, 0.0,  0.004434],
  [0.5, 0.5000, 0.0,  0.000000],
  [0.5, 0.6172, 0.0, -0.000827],
  [0.5, 0.7344, 0.0,  0.012122],
  [0.5, 0.8516, 0.0,  0.034910],
  [0.5, 0.9531, 0.0,  0.050329],
  [0.5, 0.9609, 0.0,  0.050949],
  [0.5, 0.9688, 0.0,  0.051514],
  [0.5, 0.9766, 0.0,  0.052009],
  [0.5, 1.0000, 0.0,  0.052987]])

  velocity = u.ProbeData(pts, "Velocity")
  (ilen, jlen) = velocity.shape
  pressure = u.ProbeData(pts, "Pressure")

  pts0=vtktools.arr([[0.5, 0.5, 0.0]])  # We're going to subtract off the pressure at the centre point
  press0 = u.ProbeData(pts0, "Pressure")

  norm=0.0
  for i in range(ilen):
      diff = pts[i][3] - (pressure[i][0]-press0[0][0])
      norm = norm + diff*diff

  norm = math.sqrt(norm/ilen)
  print("botella_p1_norm:", norm)

  return norm

def botella_p2(NN):
#Botella and Peyret (1998) Table 10. 
  filelist_not_sorted = glob.glob('driven_cavity-%d/*.vtu'%NN)
  vtu_nos_not_sorted = [int(file.split('.vtu')[0].split('_')[-1]) for file in filelist_not_sorted]
  filelist = [filelist_not_sorted[i] for i in numpy.argsort(vtu_nos_not_sorted)]
  file = filelist[-1]
  print(file)
  try:
    os.stat(file)
  except:
    print("No such file: %s" % file)
    sys.exit(1)

  u=vtktools.vtu(file)
  pts=vtktools.arr([
  [1.0000, 0.5, 0.0,  0.077455],
  [0.9688, 0.5, 0.0,  0.078837],
  [0.9609, 0.5, 0.0,  0.078685],
  [0.9531, 0.5, 0.0,  0.078148],
  [0.9453, 0.5, 0.0,  0.077154],
  [0.9063, 0.5, 0.0,  0.065816],
  [0.8594, 0.5, 0.0,  0.049029],
  [0.8047, 0.5, 0.0,  0.034552],
  [0.5000, 0.5, 0.0,  0.000000],
  [0.2344, 0.5, 0.0,  0.044848],
  [0.2266, 0.5, 0.0,  0.047260],
  [0.1563, 0.5, 0.0,  0.069511],
  [0.0938, 0.5, 0.0,  0.084386],
  [0.0781, 0.5, 0.0,  0.086716],
  [0.0703, 0.5, 0.0,  0.087653],
  [0.0625, 0.5, 0.0,  0.088445],
  [0.0000, 0.5, 0.0,  0.090477]])

  velocity = u.ProbeData(pts, "Velocity")
  (ilen, jlen) = velocity.shape
  pressure = u.ProbeData(pts, "Pressure")

  pts0=vtktools.arr([[0.5, 0.5, 0.0]])  # We're going to subtract off the pressure at the centre point
  press0 = u.ProbeData(pts0, "Pressure")

  norm=0.0
  for i in range(ilen):
      diff = pts[i][3] - (pressure[i][0]-press0[0][0])
      norm = norm + diff*diff

  norm = math.sqrt(norm/ilen)
  print("botella_p2_norm:", norm)

  return norm


def bruneau_ke(NN):
#Bruneau and Saad 2006. Table 7. 
  vel_l2_norm = stat_parser('driven_cavity-%d/driven_cavity.stat'%NN)['Fluid']['Velocity%magnitude']['l2norm'][-1]
  kinetic_energy = 0.5*vel_l2_norm**2
  kinetic_energy_error = abs( kinetic_energy - 0.044503 )
  print("botella_ke_error:", kinetic_energy_error)
  return kinetic_energy_error

def bruneau_sf(NN):
#Bruneau and Saad 2006. Table 2. 
  streamfunction_min = stat_parser('driven_cavity-%d/driven_cavity.stat'%NN)['Fluid']['MultiplyConnectedStreamFunction']['min'][-1]
  streamfunction_min_error = abs( streamfunction_min - -0.11892 )
  print("streamfunction_min_error:", streamfunction_min_error)
  return streamfunction_min_error



