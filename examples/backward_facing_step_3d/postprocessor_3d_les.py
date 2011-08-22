#!/usr/bin/env python

import glob
import sys
import os
import vtktools
import numpy
import pylab
import re
from math import log

def get_filelist(sample, start):

    def key(s):
        return int(s.split('_')[-1].split('.')[0])
   
    list = glob.glob("*vtu")
    list = [l for l in list if 'check' not in l]
    vtu_nos = [float(s.split('_')[-1].split('.')[0]) for s in list]
    vals = zip(vtu_nos, list)
    vals.sort()
    unzip = lambda l:tuple(apply(zip,l))
    vtu_nos, list = unzip(vals)
    shortlist = []

    for file in list:
      try:
        os.stat(file)
      except:
        f_log.write("No such file: %s" % files)
        sys.exit(1)

      ##### Start at the (start+1)th file.
      ##### Add every nth file by taking integer multiples of n; limit at 10 vtus max.
      vtu_no = float(file.split('_')[-1].split('.')[0])
      #if ((max(vtu_nos)-start)/sample > 10):
      #  sample=int((max(vtu_nos)-start)/10)
      
      if vtu_no > start:
        if (vtu_no%sample==0):
          shortlist.append(file)
        ##### Append final file if a large number of files remain.
        elif vtu_no==len(vtu_nos)-1 and (max(vtu_nos)-sample/4.0)>vtu_no:
          shortlist.append(file)

    return shortlist

#### taken from http://www.codinghorror.com/blog/archives/001018.html  #######
def tryint(s):
    try:
        return int(s)
    except:
        return s
    
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
##############################################################################
# There are shorter and more elegant version of the above, but this works
# on CX1, where this test might be run...

###################################################################

# Reattachment length:
def reatt_length(filelist, zarray):

  print "Calculating reattachment point locations using change of x-velocity sign\n"

  nums=[]; results=[]; files = []
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
    files.append(file)
  sort_nicely(files)

  for file in files:
    ##### Read in data from vtu
    datafile = vtktools.vtu(file)
    ##### Get time for plot:
    t = min(datafile.GetScalarField("Time"))
    print file, ', elapsed time = ', t
    if(t<200.):
      continue
    else:
      print "extracting data..."
      ##### points near bottom surface, 0 < x < 25
      x2array=[]; pts=[]; no_pts = 52; offset = 0.01
      x = 0.0
      for i in range(1, no_pts):
        x2array.append(x)
        for j in range(len(zarray)):
          pts.append((x, zarray[j], offset))
        x += 0.5

      x2array = numpy.array(x2array)
      pts = numpy.array(pts)

      ##### Get x-velocity on bottom boundary
      uvw = datafile.ProbeData(pts, "AverageVelocity")
      u = uvw[:,0]
      u = u.reshape([x2array.size,zarray.size])

      ##### Find all potential reattachment points:
      points = []
      for j in range(len(u[0,:])):
        for i in range(len(u[:,0])-1):
          ##### Hack to ignore division by zero entries in u.
          ##### All u should be nonzero away from boundary!
          if((u[i,j] / u[i+1,j]) < 0. and u[i+1,j] > 0. and not numpy.isinf(u[i,j] / u[i+1,j])):
            ##### interpolate between nodes
            p = x2array[i] + (x2array[i+1]-x2array[i]) * (0.0-u[i,j]) / (u[i+1,j]-u[i,j])
            print 'p: ', p
            ##### Ignore spurious corner points
            if(p>0.1):
              points.append(p)
            ##### We have our first point on this plane so...
            break

      ##### This is the spanwise-averaged reattachment point:
      if (len(points)>0):
        avpt = sum(points) / len(points)
      else:
        avpt = 0.0
      ##### Get time for plot:
      t = min(datafile.GetScalarField("Time"))
      results.append([avpt,t])

  return results

#########################################################################

# Velocity profiles:
def meanvelo(filelist,xarray,zarray,yarray):

  print "\nRunning mean velocity profile script on files at times...\n"
  ##### check for no files
  if (len(filelist) < 0):
    print "No files!"
    sys.exit(1)

  ##### create array of points
  pts=[]
  for i in range(len(xarray)):
    for j in range(len(zarray)):
      for k in range(len(yarray)):
        pts.append([xarray[i], zarray[j], yarray[k]])
  pts=numpy.array(pts)

  ##### Create output array of correct shape
  profiles=numpy.zeros([xarray.size, yarray.size], float)

  file = filelist[-1]
  #for file in filelist:
  datafile = vtktools.vtu(file)
  # Get time
  t = min(datafile.GetScalarField("Time"))
  print file, ', elapsed time = ', t

  ##### Get x-velocity
  uvw = datafile.ProbeData(pts, "AverageVelocity")
  #umax = max(abs(datafile.GetVectorField("AverageVelocity")[:,0]))

  # WARNING!!! UMAX IS NOT SAFE IN PARALLEL PERIODIC!

  umax = 1.55
  u = uvw[:,0]/umax
  u = u.reshape([xarray.size,zarray.size,yarray.size])

  ##### Spanwise averaging
  usum = numpy.zeros([xarray.size,yarray.size],float)
  usum = numpy.array(usum)
  for i in range(len(zarray)):
    uav = u[:,i,:]
    uav = numpy.array(uav)
    usum += uav
  usum = usum / len(zarray)
  profiles[:,:] = usum

  print "\n...Finished extracting data.\n"
  return profiles

#########################################################################

def plusvelo(filelist,profiles,yarray):

  print "\nCalculating mean velocity profiles in wall units.\n"
  file = filelist[-1]
  #for file in filelist:
  datafile = vtktools.vtu(file)

  ##### Velocity gradient at wall at x=19
  x1 = numpy.array([[20.0,2.0,0.0]])
  graduz = abs(datafile.ProbeData(x1, "Grad_Velocity")[:,2,0])
  visc  = datafile.ProbeData(x1, "Viscosity")[0,0,0]
  print "graduz from vtu: ", graduz

  u1 = abs(datafile.ProbeData(x1, "AverageVelocity")[:,0])
  x2 = numpy.array([[20.0,2.0,0.1]])
  u2 = abs(datafile.ProbeData(x2, "AverageVelocity")[:,0])
  grad_uz = (u2-u1)/0.1
  print "u1, u2, grad_uz linear interp: ", u1, u2, grad_uz

  ##### Wall shear stress tw, u*, y+
  twall = visc * grad_uz
  ustar = twall**0.5
  yplus = numpy.zeros([yarray.size],float)
  for i in range(len(yarray)):
    yplus[i] = ustar/visc*yarray[i]
  #print "yplus: ", yplus

  # U in wall units
  uplus = profiles/ustar

  return yplus, uplus

#########################################################################

def reynolds_stresses2(filelist,xarray,zarray,yarray):

  print "\nRunning Reynolds stress2 profile script on files at times...\n"
  ##### check for no files
  if (len(filelist) < 0):
    print "No files!"
    sys.exit(1)

  ##### create arrays of points
  pts=[]
  for i in range(len(xarray)):
    for k in range(len(yarray)):
      pts.append([xarray[i], zarray, yarray[k]])
  pts=numpy.array(pts)

  ##### Create output array of correct shape
  profiles=numpy.zeros([xarray.size, yarray.size, 3], float)

  filecount = 0
  file = filelist[-1]
  datafile = vtktools.vtu(file)
  t = min(datafile.GetScalarField("Time"))
  print file, ', elapsed time = ', t

  ##### Get instantaneous velocity components
  uvw   = datafile.ProbeData(pts, "Velocity")
  udash = datafile.ProbeData(pts, "AverageReynoldsStress")

  # max velocity for normalising
  #umax = max(abs(datafile.GetVectorField("Velocity")[:,0]))
  #print "umax: ", umax

  # WARNING!!! UMAX IS NOT SAFE IN PARALLEL PERIODIC!

  umax = 1.55
  # Reynolds components (normalised)
  u = abs(udash[:,0])/umax
  v = abs(udash[:,1])/umax
  uv= abs(udash[:,2])/umax**2
  u = u.reshape([xarray.size,yarray.size])
  v = v.reshape([xarray.size,yarray.size])
  uv = uv.reshape([xarray.size,yarray.size])

  profiles[:,:,0] = u
  profiles[:,:,1] = v
  profiles[:,:,2] = uv

  print "\n...Finished extracting data.\n"

  return profiles

#########################################################################

def plot_length(Re,type,mesh,reattachment_length):
  ##### Plot time series of reattachment length using pylab(matplotlib)

  #av_length = sum(reattachment_length[:,0]) / len(reattachment_length[:,0])
  #avg = numpy.zeros([len(reattachment_length[:,1])])
  #avg[:] = av_length
  Lemoinkim = numpy.zeros([len(reattachment_length[:,1])])
  Lemoinkim[:]=6.28

  plot1 = pylab.figure()
  pylab.title("Time series of reattachment length: Re="+str(Re)+", "+str(type)+", "+str(mesh)+" mesh")
  pylab.xlabel('Time (s)')
  pylab.ylabel('Reattachment Length (L/h)')
  pylab.plot(reattachment_length[:,1], reattachment_length[:,0], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="solid")
  #pylab.plot(reattachment_length[:,1], avg, linestyle="dashed")
  pylab.plot(reattachment_length[:,1], Lemoinkim, linestyle="dashed")
  pylab.legend(("length (step heights)","Le-Moin-Kim DNS"), loc="lower right")
  pylab.savefig("../reatt_len_3D_"+str(Re)+"_"+str(type)+"_"+str(mesh)+".pdf")
  return

#########################################################################

def plot_inlet(Re,type,mesh,vprofiles,rprofiles,xarray,zarray,yarray):
  ##### Plot velocity profiles at different points in inlet region

  # get profiles from ERCOFTAC data
  datafile = open('../Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC181-table.dat', 'r')
  print "reading in data from file: BFS-SEM-ERCOFTAC181-table.dat"
  # ignore header line
  for line in range(1):
    datafile.readline()
  
  y=[];U=[];uu=[];vv=[];uv=[]
  for line in datafile:
    y.append(1.0+float(line.split()[0]))
    U.append(0.1*float(line.split()[1]))
    uu.append(float(line.split()[2]))
    vv.append(float(line.split()[3]))
    uv.append(float(line.split()[4]))

  #print 'y: ', y

  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle("Velocity and Reynolds stresses in inflow region: Re="+str(Re)+", "+str(type)+", "+str(mesh)+" mesh", fontsize=20)

  size = 15

  ax = pylab.subplot(141)
  ax.plot(0.1*vprofiles[0,:],yarray, linestyle="solid")
  ax.plot(rprofiles[0,:,0],yarray, linestyle="solid")
  ax.plot(rprofiles[0,:,1],yarray, linestyle="solid")
  ax.plot(10*rprofiles[0,:,2],yarray, linestyle="solid")
  ax.plot(U,y, linestyle="dashed")
  ax.plot(uu,y, linestyle="dashed")
  ax.plot(vv,y, linestyle="dashed")
  ax.plot(uv,y, linestyle="dashed")
  pylab.legend(("U/10","(u'u')^1/2","(v'v')^1/2","-10*u'v'","ER U/10","ER (u'u')^1/2","ER (v'v')^1/2","ER -10*u'v'",),loc="upper right")
  ax.set_title('(a) x/h='+str(xarray[0]), fontsize=16)
  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  bx = pylab.subplot(142, sharex=ax, sharey=ax)
  bx.plot(0.1*vprofiles[1,:],yarray, linestyle="solid")
  bx.plot(rprofiles[1,:,0],yarray, linestyle="solid")
  bx.plot(rprofiles[1,:,1],yarray, linestyle="solid")
  bx.plot(10*rprofiles[1,:,2],yarray, linestyle="solid")
  bx.plot(U,y, linestyle="dashed")
  bx.plot(uu,y, linestyle="dashed")
  bx.plot(vv,y, linestyle="dashed")
  bx.plot(uv,y, linestyle="dashed")
  bx.set_title('(a) x/h='+str(xarray[1]), fontsize=16)
  #bx.grid("True")
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(143, sharex=ax, sharey=ax)
  cx.plot(0.1*vprofiles[2,:],yarray, linestyle="solid")
  cx.plot(rprofiles[2,:,0],yarray, linestyle="solid")
  cx.plot(rprofiles[2,:,1],yarray, linestyle="solid")
  cx.plot(10*rprofiles[2,:,2],yarray, linestyle="solid")
  cx.plot(U,y, linestyle="dashed")
  cx.plot(uu,y, linestyle="dashed")
  cx.plot(vv,y, linestyle="dashed")
  cx.plot(uv,y, linestyle="dashed")
  cx.set_title('(a) x/h='+str(xarray[2]), fontsize=16)
  #bx.grid("True")
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(144, sharex=ax, sharey=ax)
  dx.plot(0.1*vprofiles[3,:],yarray, linestyle="solid")
  dx.plot(rprofiles[3,:,0],yarray, linestyle="solid")
  dx.plot(rprofiles[3,:,1],yarray, linestyle="solid")
  dx.plot(10*rprofiles[3,:,2],yarray, linestyle="solid")
  dx.plot(U,y, linestyle="dashed")
  dx.plot(uu,y, linestyle="dashed")
  dx.plot(vv,y, linestyle="dashed")
  dx.plot(uv,y, linestyle="dashed")
  dx.set_title('(a) x/h='+str(xarray[3]), fontsize=16)
  #bx.grid("True")
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  pylab.axis([-0.01, 0.11, 1., 2.])
  bx.set_xlabel('Normalised U and Reynolds stresses', fontsize=24)
  ax.set_ylabel('y/h', fontsize=24)

  print "saving inlet plot"

  pylab.savefig("../inlet_profiles_3d"+str(Re)+"_"+str(type)+"_"+str(mesh)+".pdf")
  return

#########################################################################

def plot_meanvelo(Re,type,mesh,vprofiles,xarray,yarray):

  # get profiles from ERCOFTAC data
  datafile = open('../Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC360-table.dat', 'r')
  print "reading in data from file: BFS-SEM-ERCOFTAC360-table.dat"
  # ignore header line
  for line in range(1):
    datafile.readline()
  
  y4=[];U4=[]
  for line in datafile:
    y4.append(float(line.split()[0]))
    U4.append(float(line.split()[1]))

  datafile = open('../Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC411-table.dat', 'r')
  print "reading in data from file: BFS-SEM-ERCOFTAC411-table.dat"
  # ignore header line
  for line in range(1):
    datafile.readline()
  
  y6=[];U6=[]
  for line in datafile:
    y6.append(float(line.split()[0]))
    U6.append(float(line.split()[1]))

  datafile = open('../Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC513-table.dat', 'r')
  print "reading in data from file: BFS-SEM-ERCOFTAC513-table.dat"
  # ignore header line
  for line in range(1):
    datafile.readline()
  
  y10=[];U10=[]
  for line in datafile:
    y10.append(float(line.split()[0]))
    U10.append(float(line.split()[1]))

  datafile = open('../Ercoftac-test31-BFS/BFS-SEM-ERCOFTAC744-table.dat', 'r')
  print "reading in data from file: BFS-SEM-ERCOFTAC744-table.dat"
  # ignore header line
  for line in range(1):
    datafile.readline()
  
  y19=[];U19=[]
  for line in datafile:
    y19.append(float(line.split()[0]))
    U19.append(float(line.split()[1]))

  # get profiles from Le&Moin U graph. x=4
  Le = open('../Le-profiles/Le-profile1-U-x4.dat', 'r').readlines()
  Le_u4 = [float(line.split()[0]) for line in Le]
  Le_y4 = [float(line.split()[1]) for line in Le]
  jd = open('../Le-profiles/JD-profile1-U-x4.dat', 'r').readlines()
  jd_u4 = [float(line.split()[0]) for line in jd]
  jd_y4 = [float(line.split()[1]) for line in jd]

  # get profiles from Le&Moin U graph. x=6
  Le = open('../Le-profiles/Le-profile1-U-x6.dat', 'r').readlines()
  Le_u6 = [float(line.split()[0]) for line in Le]
  Le_y6 = [float(line.split()[1]) for line in Le]
  jd = open('../Le-profiles/JD-profile1-U-x6.dat', 'r').readlines()
  jd_u6 = [float(line.split()[0]) for line in jd]
  jd_y6 = [float(line.split()[1]) for line in jd]

  # get profiles from Le&Moin U graph. x=10
  Le = open('../Le-profiles/Le-profile1-U-x10.dat', 'r').readlines()
  Le_u10 = [float(line.split()[0]) for line in Le]
  Le_y10 = [float(line.split()[1]) for line in Le]
  jd = open('../Le-profiles/JD-profile1-U-x10.dat', 'r').readlines()
  jd_u10 = [float(line.split()[0]) for line in jd]
  jd_y10 = [float(line.split()[1]) for line in jd]

  # get profiles from Le&Moin U graph. x=19
  Le = open('../Le-profiles/Le-profile1-U-x19.dat', 'r').readlines()
  Le_u19 = [float(line.split()[0]) for line in Le]
  Le_y19 = [float(line.split()[1]) for line in Le]
  jd = open('../Le-profiles/JD-profile1-U-x19.dat', 'r').readlines()
  jd_u19 = [float(line.split()[0]) for line in jd]
  jd_y19 = [float(line.split()[1]) for line in jd]

  ##### Plot velocity profiles at different points behind step using pylab(matplotlib)
  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle("Evolution of mean U-velocity: Re="+str(Re)+", "+str(type)+", "+str(mesh)+" mesh", fontsize=20)

  size = 15

  ax = pylab.subplot(141)
  ax.plot(vprofiles[0,:],yarray, linestyle="solid")
  ax.plot(U4,y4, linestyle="dashed")
  ax.plot(jd_u4,jd_y4, linestyle="none",marker='o',color='black')
  ax.set_title('(a) x/h='+str(xarray[0]), fontsize=16)
  pylab.legend(('Fluidity',"Le&Moin DNS","Jovic&Driver expt"),loc="upper left")
  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  bx = pylab.subplot(142, sharex=ax, sharey=ax)
  bx.plot(vprofiles[1,:],yarray, linestyle="solid")
  bx.plot(U6,y6, linestyle="dashed")
  bx.plot(jd_u6,jd_y6, linestyle="none",marker='o',color='black')
  bx.set_title('(a) x/h='+str(xarray[1]), fontsize=16)
  #bx.grid("True")
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(143, sharex=ax, sharey=ax)
  cx.plot(vprofiles[2,:],yarray, linestyle="solid")
  cx.plot(U10,y10, linestyle="dashed")
  cx.plot(jd_u10,jd_y10, linestyle="none",marker='o',color='black')
  cx.set_title('(a) x/h='+str(xarray[2]), fontsize=16)
  #bx.grid("True")
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(144, sharex=ax, sharey=ax)
  dx.plot(vprofiles[3,:],yarray, linestyle="solid")
  dx.plot(U19,y19, linestyle="dashed")
  dx.plot(jd_u19,jd_y19, linestyle="none",marker='o',color='black')
  dx.set_title('(a) x/h='+str(xarray[3]), fontsize=16)
  #bx.grid("True")
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  pylab.axis([-0.25, 1., 0., 3.])
  bx.set_xlabel('Normalised mean U-velocity (U/Umax)', fontsize=24)
  ax.set_ylabel('y/h', fontsize=24)

  pylab.savefig("../mean_velo_profiles_3d"+str(Re)+"_"+str(type)+"_"+str(mesh)+".pdf")
  return

#########################################################################

def plot_plusvelo(Re,type,mesh,uplus,yplus,xarray):

  print "\nPlotting mean velocity profiles in wall units...\n"
  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle("Streamwise Evolution of U+ (velocity in wall units): Re="+str(Re)+", "+str(type)+", "+str(mesh)+" mesh", fontsize=20)

  # get profiles from Le&Moin uplus graph. x=19
  Le = open('../Le-profiles/Le-profile1-uplus-x19.dat', 'r').readlines()
  Le_uplus = [float(line.split()[1]) for line in Le]
  Le_yplus = [float(line.split()[0]) for line in Le]
  jd = open('../Le-profiles/JD-profile1-uplus-x19.dat', 'r').readlines()
  jd_uplus = [float(line.split()[1]) for line in jd]
  jd_yplus = [float(line.split()[0]) for line in jd]
  size = 15

  # log law of the wall
  loglaw=numpy.zeros(yplus.size)
  for i in range(len(yplus)):
    if(yplus[i]>10.0):
      loglaw[i] = 1/0.41*log(yplus[i]) + 5.0
    else:
      loglaw[i] = yplus[i]

  ax = pylab.subplot(111)
  ax.plot(yplus, uplus[0,:], linestyle="dotted")
  ax.plot(yplus, uplus[1,:], linestyle="dotted")
  ax.plot(yplus, uplus[2,:], linestyle="dotted")
  ax.plot(yplus, uplus[3,:], linestyle="dotted")
  ax.plot(yplus, uplus[4,:], linestyle="solid", color="black")
  # reference data:
  ax.plot(Le_yplus, Le_uplus, linestyle="dashed", color="black")
  ax.plot(jd_yplus, jd_uplus, linestyle="none", color="black", marker = 'o', markerfacecolor='black', markersize=6, markeredgecolor='black')
  ax.plot(yplus, loglaw, linestyle="dashdot", color="black")
  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  pylab.legend(("x/h=10.0","x/h=12.5","x/h=15.0","x/h=17.5","x/h=19","Le&Moin DNS x/h=19","Jovic&Driver expt x/h=19", "log law"),loc="lower right")
  ax.set_ylabel('mean U-velocity in wall units U+', fontsize=24)
  ax.set_xlabel('y+', fontsize=24)
  ax.set_xscale('log')
  ax.set_xlim(1, 1000)
  #ax.set_ylim(10, 20)
  pylab.savefig("../velo_plus_profiles_3d"+str(Re)+"_"+str(type)+"_"+str(mesh)+".pdf")
  return

#########################################################################

def plot_reynolds_stresses2(Re,type,mesh,rprofiles,xarray,zarray,yarray):
  ##### Plot Reynolds stress profiles at different points behind step using pylab(matplotlib)
  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle("Time-averaged Reynolds stress profiles: Re="+str(Re)+", "+str(type)+", "+str(mesh)+" mesh", fontsize=20)

  # get profiles from Le&Moin graphs. x=4
  Le_uu4 = open('../Le-profiles/Le-profile1-uu-x4.dat', 'r').readlines()
  rs_uu4 = [float(line.split()[0]) for line in Le_uu4]
  y_uu4 = [float(line.split()[1]) for line in Le_uu4]
  Le_vv4 = open('../Le-profiles/Le-profile1-vv-x4.dat', 'r').readlines()
  rs_vv4 = [float(line.split()[0]) for line in Le_vv4]
  y_vv4 = [float(line.split()[1]) for line in Le_vv4]
  Le_uv4 = open('../Le-profiles/Le-profile1-uv-x4.dat', 'r').readlines()
  rs_uv4 = [float(line.split()[0]) for line in Le_uv4]
  y_uv4 = [float(line.split()[1]) for line in Le_uv4]

  # get profiles from Le&Moin graphs. x=6
  Le_uu6 = open('../Le-profiles/Le-profile1-uu-x6.dat', 'r').readlines()
  rs_uu6 = [float(line.split()[0]) for line in Le_uu6]
  y_uu6 = [float(line.split()[1]) for line in Le_uu6]
  Le_vv6 = open('../Le-profiles/Le-profile1-vv-x6.dat', 'r').readlines()
  rs_vv6 = [float(line.split()[0]) for line in Le_vv6]
  y_vv6 = [float(line.split()[1]) for line in Le_vv6]
  Le_uv6 = open('../Le-profiles/Le-profile1-uv-x6.dat', 'r').readlines()
  rs_uv6 = [float(line.split()[0]) for line in Le_uv6]
  y_uv6 = [float(line.split()[1]) for line in Le_uv6]

  # get profiles from Le&Moin graphs. x=10
  Le_uu10 = open('../Le-profiles/Le-profile1-uu-x10.dat', 'r').readlines()
  rs_uu10 = [float(line.split()[0]) for line in Le_uu10]
  y_uu10 = [float(line.split()[1]) for line in Le_uu10]
  Le_vv10 = open('../Le-profiles/Le-profile1-vv-x10.dat', 'r').readlines()
  rs_vv10 = [float(line.split()[0]) for line in Le_vv10]
  y_vv10 = [float(line.split()[1]) for line in Le_vv10]
  Le_uv10 = open('../Le-profiles/Le-profile1-uv-x10.dat', 'r').readlines()
  rs_uv10 = [float(line.split()[0]) for line in Le_uv10]
  y_uv10 = [float(line.split()[1]) for line in Le_uv10]

  # get profiles from Le&Moin graphs. x=19
  Le_uu19 = open('../Le-profiles/Le-profile1-uu-x19.dat', 'r').readlines()
  rs_uu19 = [float(line.split()[0]) for line in Le_uu19]
  y_uu19 = [float(line.split()[1]) for line in Le_uu19]
  Le_vv19 = open('../Le-profiles/Le-profile1-vv-x19.dat', 'r').readlines()
  rs_vv19 = [float(line.split()[0]) for line in Le_vv19]
  y_vv19 = [float(line.split()[1]) for line in Le_vv19]
  Le_uv19 = open('../Le-profiles/Le-profile1-uv-x19.dat', 'r').readlines()
  rs_uv19 = [float(line.split()[0]) for line in Le_uv19]
  y_uv19 = [float(line.split()[1]) for line in Le_uv19]

  size = 15

  ax = pylab.subplot(141)
  ax.plot(rprofiles[0,:,0],yarray, linestyle="solid")
  ax.plot(rprofiles[0,:,1],yarray, linestyle="solid")
  ax.plot(rprofiles[0,:,2],yarray, linestyle="solid")
  ax.plot(rs_uu4,y_uu4, linestyle="dashed")
  ax.plot(rs_vv4,y_vv4, linestyle="dashed")
  ax.plot(rs_uv4,y_uv4, linestyle="dashed")
  pylab.legend(("(u'u')^1/2","(v'v')^1/2","-u'v'","Le (u'u')^1/2","Le (v'v')^1/2","Le -u'v'",),loc="upper right")
  ax.set_title('(a) x/h='+str(xarray[0]), fontsize=16)
  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  bx = pylab.subplot(142, sharex=ax, sharey=ax)
  bx.plot(rprofiles[1,:,0],yarray, linestyle="solid")
  bx.plot(rprofiles[1,:,1],yarray, linestyle="solid")
  bx.plot(rprofiles[1,:,2],yarray, linestyle="solid")
  bx.plot(rs_uu6,y_uu6, linestyle="dashed")
  bx.plot(rs_vv6,y_vv6, linestyle="dashed")
  bx.plot(rs_uv6,y_uv6, linestyle="dashed")
  bx.set_title('(a) x/h='+str(xarray[1]), fontsize=16)
  #bx.grid("True")
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(143, sharex=ax, sharey=ax)
  cx.plot(rprofiles[2,:,0],yarray, linestyle="solid")
  cx.plot(rprofiles[2,:,1],yarray, linestyle="solid")
  cx.plot(rprofiles[2,:,2],yarray, linestyle="solid")
  cx.plot(rs_uu10,y_uu10, linestyle="dashed")
  cx.plot(rs_vv10,y_vv10, linestyle="dashed")
  cx.plot(rs_uv10,y_uv10, linestyle="dashed")
  cx.set_title('(a) x/h='+str(xarray[2]), fontsize=16)
  #bx.grid("True")
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(144, sharex=ax, sharey=ax)
  dx.plot(rprofiles[3,:,0],yarray, linestyle="solid")
  dx.plot(rprofiles[3,:,1],yarray, linestyle="solid")
  dx.plot(rprofiles[3,:,2],yarray, linestyle="solid")
  dx.plot(rs_uu19,y_uu19, linestyle="dashed")
  dx.plot(rs_vv19,y_vv19, linestyle="dashed")
  dx.plot(rs_uv19,y_uv19, linestyle="dashed")
  dx.set_title('(a) x/h='+str(xarray[3]), fontsize=16)
  #bx.grid("True")
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  pylab.axis([-0.1, 0.19, 0., 3.])
  bx.set_xlabel('Normalised Reynolds stresses2 vs. Le&Moin data (divided by 2)', fontsize=20)
  ax.set_ylabel('y/h', fontsize=24)

  pylab.savefig("../reynolds_stress2_profiles_3d"+str(Re)+"_"+str(type)+"_"+str(mesh)+".pdf")
  return

#########################################################################

def main():
    ##### Which run is being processed?
    Re = sys.argv[1]
    type = sys.argv[2]
    mesh = sys.argv[3]
    print "Re, bc type, mesh: ", Re, type, mesh

    ##### Only process every nth file by taking integer multiples of n:
    filelist = get_filelist(sample=1, start=7)

    ##### Points to generate profiles:
    xarray = numpy.array([4.0, 6.0, 10.0, 19.0])
    zarray = numpy.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9])
    yarray = numpy.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0])

    ##### Call reattachment_length function
    reattachment_length = numpy.array(reatt_length(filelist, zarray))
    numpy.save("../numpy_data/reatt_len_"+str(Re)+"_"+str(type)+"_"+str(mesh), reattachment_length)
    plot_length(Re,type,mesh,reattachment_length)

    ##### Call meanvelo function
    zarray = numpy.array([2.0])
    vprofiles = meanvelo(filelist, xarray, zarray, yarray)
    numpy.save("../numpy_data/mean_velo_"+str(Re)+"_"+str(type)+"_"+str(mesh), vprofiles)
    print "Showing plot of velocity profiles."
    plot_meanvelo(Re,type,mesh,vprofiles,xarray,yarray)
    # points used by Le & Moin in U+ plot:
    xarray=numpy.array([10.0,12.5,15.0,17.5,19.0])
    vprofiles = meanvelo(filelist, xarray, zarray, yarray)
    yplus, uplus = plusvelo(filelist,vprofiles, yarray)
    numpy.save('../numpy_data/yplus_'+str(Re)+"_"+str(type)+"_"+str(mesh), yplus)
    numpy.save('../numpy_data/uplus_'+str(Re)+"_"+str(type)+"_"+str(mesh), uplus)
    plot_plusvelo(Re,type,mesh,uplus,yplus,xarray)

    ##### Call Reynolds stress2 function
    rprofiles2 = reynolds_stresses2(filelist, xarray, zarray, yarray)
    numpy.save("../numpy_data/re_stress_"+str(Re)+"_"+str(type)+"_"+str(mesh), rprofiles2)
    plot_reynolds_stresses2(Re,type,mesh,rprofiles2,xarray,zarray,yarray)

    ##### Plot inlet region
    xarray = numpy.array([-10.0, -7.0, -3.0, 0.0])
    yarray = numpy.array([1.0,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.5,4.6,4.7,4.8,4.9,5.0])
    rprofiles = reynolds_stresses2(filelist, xarray, zarray, yarray)
    zarray = numpy.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9])
    vprofiles = meanvelo(filelist, xarray, zarray, yarray)
    numpy.save("../numpy_data/inlet_profiles_"+str(Re)+"_"+str(type)+"_"+str(mesh), vprofiles)
    print "Showing plot of inlet profiles."
    plot_inlet(Re,type,mesh,vprofiles,rprofiles,xarray,zarray,yarray)
    #pylab.show()

    print "\nAll done.\n"

if __name__ == "__main__":
    sys.exit(main())

