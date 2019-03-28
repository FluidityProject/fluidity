#!/usr/bin/env python

from __future__ import print_function

import glob
import sys
import os
import vtktools
import numpy
import pylab
import re
from matplotlib import rc

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
      ##### Add every nth file by taking integer multiples of n.
      vtu_no = float(file.split('_')[-1].split('.')[0])
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
def reattachment_length(filelist):

  print("Calculating reattachment point locations using change of x-velocity sign\n")

  nums=[]; results=[]; files = []
  ##### check for no files
  if (len(filelist) == 0):
    print("No files!")
    sys.exit(1)
  for file in filelist:
    try:
      os.stat(file)
    except:
      print("No such file: %s" % file)
      sys.exit(1)
    files.append(file)
  sort_nicely(files)

  for file in files:
    ##### Read in data from vtu
    datafile = vtktools.vtu(file)
    ##### Get time for plot:
    t = min(datafile.GetScalarField("Time"))
    print(file, ', elapsed time = ', t)

    ##### points near bottom surface, 0 < x < 20
    pts=[]; no_pts = 82; offset = 0.01
    x = 5.0
    for i in range(1, no_pts):
      pts.append((x, offset, 0.0))
      x += 0.25

    pts = numpy.array(pts)

    ##### Get x-velocity on bottom boundary
    uvw = datafile.ProbeData(pts, "AverageVelocity")
    u = []
    u = uvw[:,0]
    points = 0.0

    for i in range(len(u)-1):
      ##### Hack to ignore division by zero entries in u.
      ##### All u should be nonzero away from boundary!
      if((u[i] / u[i+1]) < 0. and u[i+1] > 0. and not numpy.isinf(u[i] / u[i+1])):
        ##### interpolate between nodes. Correct for origin not at step.
        p = pts[i][0] + (pts[i+1][0]-pts[i][0]) * (0.0-u[i]) / (u[i+1]-u[i]) -5.0
        print('p ', p)
        ##### Ignore spurious corner points
        if(p>2):
          points = p
          ##### We have our first point on this plane so...
          break
    print("reattachment point found at: ", points)

    ##### Append actual reattachment point and time:
    results.append([points,t])

  return results

#########################################################################

# Velocity profiles:
def meanvelo(file,x,y):

  print("\nRunning velocity profile script on files at times...\n")

  ##### create array of points. Correct for origin not at step.
  pts=[]
  for i in range(len(x)):
    for j in range(len(y)):
      pts.append([x[i]+5.0, y[j], 0.0])

  pts=numpy.array(pts)
  profiles=numpy.zeros([x.size, y.size], float)

  datafile = vtktools.vtu(file)

  ##### Get x-velocity
  uvw = datafile.ProbeData(pts, "AverageVelocity")
  u = uvw[:,0]
  u=u.reshape([x.size,y.size])
  for i in range(len(x)):
    umax = max(u[i,:])
    u[i,:] = u[i,:]/umax
  profiles[:,:] = u

  print("\n...Finished writing data files.\n")
  return profiles

#########################################################################

def plot_length(type,reattachment_length):
  ##### Plot time series of reattachment length using pylab(matplotlib)
  ##### Kim's (1978) experimental result, Re=132000
  kim = numpy.zeros([len(reattachment_length[:,1])])
  kim[:]=7.0
  ##### Ilinca's (1997) best numerical result (adaptive mesh 2, 6960 points)
  ilinca = numpy.zeros([len(reattachment_length[:,1])])
  ilinca[:]=6.21

  plot1 = pylab.figure()
  pylab.title("Time series of reattachment length: Re=132000, "+str(type))
  pylab.xlabel('Time (s)')
  pylab.ylabel('Reattachment Length (L/h)')
  pylab.plot(reattachment_length[:,1], reattachment_length[:,0], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="solid",color='blue')
  pylab.plot(reattachment_length[:,1], kim, linestyle="solid",color='black')
  pylab.plot(reattachment_length[:,1], ilinca, linestyle="solid",color='red')
  pylab.legend(("Fluidity","Kim expt.","Ilinca sim."), loc="best")
  rlmax = max(reattachment_length[:,0])
  pylab.axis([0, reattachment_length[-1,1], 0, max(kim[-1],rlmax)+0.1])
  pylab.savefig("../reattachment_length_kim_"+str(type)+".pdf")
  return

#########################################################################

def plot_meanvelo(type,profiles,xarray,yarray):
  ##### Plot evolution of velocity profiles at different points behind step

  # get profiles from Ilinca's experimental/numerical data
  datafile = open('../Ilinca-data/Ilinca-U-expt-1.33.dat', 'r')
  print("reading in data from file: Ilinca-U-expt-1.33.dat")
  y1=[];U1=[]
  for line in datafile:
    U1.append(float(line.split()[0]))
    y1.append(float(line.split()[1]))
  datafile = open('../Ilinca-data/Ilinca-U-num-1.33.dat', 'r')
  print("reading in data from file: Ilinca-U-num-1.33.dat")
  yn1=[];Un1=[]
  for line in datafile:
    Un1.append(float(line.split()[0]))
    yn1.append(float(line.split()[1]))

  datafile = open('../Ilinca-data/Ilinca-U-expt-2.66.dat', 'r')
  print("reading in data from file: Ilinca-U-expt-2.66.dat")
  y3=[];U3=[]
  for line in datafile:
    U3.append(float(line.split()[0]))
    y3.append(float(line.split()[1]))
  datafile = open('../Ilinca-data/Ilinca-U-num-2.66.dat', 'r')
  print("reading in data from file: Ilinca-U-num-2.66.dat")
  yn3=[];Un3=[]
  for line in datafile:
    Un3.append(float(line.split()[0]))
    yn3.append(float(line.split()[1]))

  datafile = open('../Ilinca-data/Ilinca-U-expt-5.33.dat', 'r')
  print("reading in data from file: Ilinca-U-expt-5.33.dat")
  y5=[];U5=[]
  for line in datafile:
    U5.append(float(line.split()[0]))
    y5.append(float(line.split()[1]))
  datafile = open('../Ilinca-data/Ilinca-U-num-5.33.dat', 'r')
  print("reading in data from file: Ilinca-U-num-5.33.dat")
  yn5=[];Un5=[]
  for line in datafile:
    Un5.append(float(line.split()[0]))
    yn5.append(float(line.split()[1]))

  datafile = open('../Ilinca-data/Ilinca-U-expt-8.0.dat', 'r')
  print("reading in data from file: Ilinca-U-expt-8.0.dat")
  y8=[];U8=[]
  for line in datafile:
    U8.append(float(line.split()[0]))
    y8.append(float(line.split()[1]))
  datafile = open('../Ilinca-data/Ilinca-U-num-8.0.dat', 'r')
  print("reading in data from file: Ilinca-U-num-8.0.dat")
  yn8=[];Un8=[]
  for line in datafile:
    Un8.append(float(line.split()[0]))
    yn8.append(float(line.split()[1]))

  datafile = open('../Ilinca-data/Ilinca-U-expt-16.0.dat', 'r')
  print("reading in data from file: Ilinca-U-expt-16.0.dat")
  y16=[];U16=[]
  for line in datafile:
    U16.append(float(line.split()[0]))
    y16.append(float(line.split()[1]))
  datafile = open('../Ilinca-data/Ilinca-U-num-16.0.dat', 'r')
  print("reading in data from file: Ilinca-U-num-16.0.dat")
  yn16=[];Un16=[]
  for line in datafile:
    Un16.append(float(line.split()[0]))
    yn16.append(float(line.split()[1]))

  plot1 = pylab.figure(figsize = (20, 8))
  pylab.suptitle("U-velocity profile: Re=132000, "+str(type), fontsize=20)

  size = 15
  ax = pylab.subplot(151)
  leg_end = []

  ax.plot(profiles[0,:],yarray, linestyle="solid",color='blue', marker = 'o', markerfacecolor='white', markersize=8, markeredgecolor='black')
  ax.plot(U1,y1, linestyle='none',color='black', marker = 'x', markerfacecolor='black', markersize=8)
  ax.plot(Un1,yn1, linestyle="solid",color="red")

  leg_end.append("Fluidity")
  leg_end.append("Kim expt.")
  leg_end.append("Ilinca sim.")
  pylab.legend((leg_end), loc="lower right")
  leg = pylab.gca().get_legend()
  ltext = leg.get_texts()
  pylab.setp(ltext, fontsize = 18, color = 'black')
  frame=leg.get_frame()
  frame.set_fill(False)
  frame.set_visible(False)

  ax.set_title('(a) x/h='+str(xarray[0]), fontsize=16)
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  bx = pylab.subplot(152, sharex=ax, sharey=ax)
  bx.plot(profiles[1,:],yarray, linestyle="solid",color='blue', marker = 'o', markerfacecolor='white', markersize=8, markeredgecolor='black')
  bx.plot(U3,y3, linestyle='none',color='black', marker = 'x', markerfacecolor='black', markersize=8)
  bx.plot(Un3,yn3, linestyle="solid",color='red')
  bx.set_title('(b) x/h='+str(xarray[1]), fontsize=16)
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(153, sharex=ax, sharey=ax)
  cx.plot(profiles[2,:],yarray, linestyle="solid",color='blue', marker = 'o', markerfacecolor='white', markersize=8, markeredgecolor='black')
  cx.plot(U5,y5, linestyle='none',color='black', marker = 'x', markerfacecolor='black', markersize=8)
  cx.plot(Un5,yn5, linestyle="solid",color='red')
  cx.set_title('(c) x/h='+str(xarray[2]), fontsize=16)
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(154, sharex=ax, sharey=ax)
  dx.plot(profiles[3,:],yarray, linestyle="solid",color='blue', marker = 'o', markerfacecolor='white', markersize=8, markeredgecolor='black')
  dx.plot(U8,y8, linestyle='none',color='black', marker = 'x', markerfacecolor='black', markersize=8)
  dx.plot(Un8,yn8, linestyle="solid",color='red')
  dx.set_title('(d) x/h='+str(xarray[3]), fontsize=16)
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  ex = pylab.subplot(155, sharex=ax, sharey=ax)
  ex.plot(profiles[4,:],yarray, linestyle="solid",color='blue', marker = 'o', markerfacecolor='white', markersize=8, markeredgecolor='black')
  ex.plot(U16,y16, linestyle='none',color='black', marker = 'x', markerfacecolor='black', markersize=8)
  ex.plot(Un16,yn16, linestyle="solid",color='red')
  ex.set_title('(e) x/h='+str(xarray[4]), fontsize=16)
  for tick in ex.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(ex.get_yticklabels(), visible=False)

  pylab.axis([-0.25, 1.15, 0., 3.])
  cx.set_xlabel('Normalised U-velocity (U/Umax)', fontsize=24)
  ax.set_ylabel('y', fontsize=24)

  pylab.savefig("../velocity_profiles_kim_"+str(type)+".pdf")
  return

#########################################################################

def main():
    ##### Which run is being processed?
    type = sys.argv[1]

    ##### Only process every nth file:
    filelist = get_filelist(sample=5, start=0)

    ##### Call reattachment_length function
    reatt_length = numpy.array(reattachment_length(filelist))
    av_length = sum(reatt_length[:,0]) / len(reatt_length[:,0])
    numpy.save("reattachment_length_kim_"+str(type), reatt_length)
    plot_length(type,reatt_length)

    ##### Points to generate profiles:
    xarray = numpy.array([1.33,2.66,5.33,8.0,16.0])
    yarray = numpy.array([0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,
2.5,2.6,2.7,2.8,2.9,2.95,2.96,2.97,2.98,2.99,3.0])

    ##### Call meanvelo function
    profiles = meanvelo(filelist[-1], xarray, yarray)
    numpy.save("velocity_profiles_kim_"+str(type), profiles)
    plot_meanvelo(type,profiles,xarray,yarray)
    pylab.show()

    print("\nAll done.\n")


if __name__ == "__main__":
    sys.exit(main())


