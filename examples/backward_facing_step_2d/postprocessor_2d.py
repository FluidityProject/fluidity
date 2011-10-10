#!/usr/bin/env python

import glob
import sys
import os
import vtktools
import numpy
import pylab
import re

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

    ##### points near bottom surface, 0 < x < 20
    pts=[]; no_pts = 82; offset = 0.01
    x = 5.0
    for i in range(1, no_pts):
      pts.append((x, offset, 0.0))
      x += 0.25

    pts = numpy.array(pts)

    ##### Get x-velocity on bottom boundary
    uvw = datafile.ProbeData(pts, "Velocity")
    u = []
    u = uvw[:,0]
    points = 0.0

    for i in range(len(u)-1):
      ##### Hack to ignore division by zero entries in u.
      ##### All u should be nonzero away from boundary!
      if((u[i] / u[i+1]) < 0. and u[i+1] > 0. and not numpy.isinf(u[i] / u[i+1])):
        ##### interpolate between nodes. Correct for origin not at step.
        p = pts[i][0] + (pts[i+1][0]-pts[i][0]) * (0.0-u[i]) / (u[i+1]-u[i]) -5.0
        ##### Ignore spurious corner points
        if(p>0.1):
          points = p
        ##### We have our first point on this plane so...
        break
    print "reattachment point found at: ", points

    ##### Append actual reattachment point and time:
    results.append([points,t])

  return results

#########################################################################

# Velocity profiles:
def meanvelo(filelist,x,y):

  print "\nRunning velocity profile script on files at times...\n"
  ##### check for no files
  if (len(filelist) < 0):
    print "No files!"
    sys.exit(1)

  ##### create array of points. Correct for origin not at step.
  pts=[]
  for i in range(len(x)):
    for j in range(len(y)):
      pts.append([x[i]+5.0, y[j], 0.0])

  pts=numpy.array(pts)
  ##### Create output array of correct shape
  profiles=numpy.zeros([len(filelist), x.size, y.size], float)
  time = numpy.zeros([len(filelist)], float)

  filecount = 0
  for file in filelist:
      datafile = vtktools.vtu(file)
      # Get time
      t = min(datafile.GetScalarField("Time"))
      print file, ', elapsed time = ', t
      time[filecount] = t

      ##### Get x-velocity
      uvw = datafile.ProbeData(pts, "Velocity")
      umax = max(abs(datafile.GetVectorField("Velocity")[:,0]))
      print 'umax: ', umax
      u = uvw[:,0]/umax
      u=u.reshape([x.size,y.size])
      profiles[filecount,:,:] = u
      filecount += 1

  print "\n...Finished writing data files.\n"
  return profiles, time

#########################################################################

def plot_length(type,reattachment_length):
  ##### Plot time series of reattachment length using pylab(matplotlib)
  plot1 = pylab.figure()
  pylab.title("Time series of reattachment length: Re=4700, "+str(type))
  pylab.xlabel('Time (s)')
  pylab.ylabel('Reattachment Length (L/h)')
  pylab.plot(reattachment_length[:,1], reattachment_length[:,0], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="solid")
  pylab.savefig("../reattachment_length_"+str(type)+".pdf")
  return

def plot_meanvelo(type,profiles,xarray,yarray,time):
  ##### Plot velocity profiles at different points behind step, and at 3 times using pylab(matplotlib)
  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle("Evolution of U-velocity: Re=4700, "+str(type), fontsize=20)

  # get profiles from Armaly's experimental data
  datafile = open('../Armaly-data/armaly-velo-3.06.dat', 'r')
  print "reading in data from file: armaly-velo-3.06.dat"
  y3=[];U3=[]
  for line in datafile:
    U3.append(float(line.split()[0])/36)
    y3.append(float(line.split()[1]))
  # normalise
  #U3=[U3[i]/36 for i in range(len(U3))]

  datafile = open('../Armaly-data/armaly-velo-6.12.dat', 'r')
  print "reading in data from file: armaly-velo-6.12.dat"
  y6=[];U6=[]
  for line in datafile:
    U6.append(float(line.split()[0])/36)
    y6.append(float(line.split()[1]))
  #U6=[U6[i]/max(U6) for i in range(len(U6))]

  datafile = open('../Armaly-data/armaly-velo-10.20.dat', 'r')
  print "reading in data from file: armaly-velo-10.20.dat"
  y10=[];U10=[]
  for line in datafile:
    U10.append(float(line.split()[0])/36)
    y10.append(float(line.split()[1]))
  #U10=[U10[i]/max(U10) for i in range(len(U10))]

  datafile = open('../Armaly-data/armaly-velo-15.31.dat', 'r')
  print "reading in data from file: armaly-velo-15.31.dat"
  y15=[];U15=[]
  for line in datafile:
    U15.append(float(line.split()[0])/36)
    y15.append(float(line.split()[1]))
  #U15=[U15[i]/max(U15) for i in range(len(U15))]

  size = 15
  ax = pylab.subplot(141)
  shift=0.0
  leg_end = []

  for i in range(len(time)):
    if(i==len(time)-1):
      ax.plot(profiles[i,0,:]+shift,yarray, linestyle="solid")
    else:
      ax.plot(profiles[i,0,:]+shift,yarray, linestyle="dashed")
    shift+=0.0
    leg_end.append("%.1f secs"%time[i])
  leg_end.append("Armaly data")
  ax.plot(U3,y3, linestyle="solid",color="black")
  pylab.legend((leg_end), loc="lower right")
  ax.set_title('(a) x/h='+str(xarray[0]), fontsize=16)
  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  bx = pylab.subplot(142, sharex=ax, sharey=ax)
  shift=0.0
  for i in range(len(time)):
    if(i==len(time)-1):
      bx.plot(profiles[i,1,:]+shift,yarray, linestyle="solid")
    else:
      bx.plot(profiles[i,1,:]+shift,yarray, linestyle="dashed")
    shift+=0.0
  bx.plot(U6,y6, linestyle="solid",color='black')
  bx.set_title('(a) x/h='+str(xarray[1]), fontsize=16)
  #bx.grid("True")
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(143, sharex=ax, sharey=ax)
  shift=0.0
  for i in range(len(time)):
    if(i==len(time)-1):
      cx.plot(profiles[i,2,:]+shift,yarray, linestyle="solid")
    else:
      cx.plot(profiles[i,2,:]+shift,yarray, linestyle="dashed")
    shift+=0.0
  cx.plot(U10,y10, linestyle="solid",color='black')
  cx.set_title('(a) x/h='+str(xarray[2]), fontsize=16)
  #bx.grid("True")
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(144, sharex=ax, sharey=ax)
  shift=0.0
  for i in range(len(time)):
    if(i==len(time)-1):
      dx.plot(profiles[i,3,:]+shift,yarray, linestyle="solid")
    else:
      dx.plot(profiles[i,3,:]+shift,yarray, linestyle="dashed")
    shift+=0.0
  dx.plot(U15,y15, linestyle="solid",color='black')
  dx.set_title('(a) x/h='+str(xarray[3]), fontsize=16)
  #bx.grid("True")
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  #pylab.axis([-0.2, 0.8, 0., 1.94])
  bx.set_xlabel('Normalised U-velocity (U/Umax)', fontsize=24)
  ax.set_ylabel('z/h', fontsize=24)

  pylab.savefig("../velocity_profiles_"+str(type)+".pdf")
  return

#########################################################################

def main():
    ##### Which run is being processed?
    type = sys.argv[1]

    ##### Only process every nth file by taking integer multiples of n:
    filelist = get_filelist(sample=20, start=0)

    ##### Call reattachment_length function
    reatt_length = numpy.array(reattachment_length(filelist))
    av_length = sum(reatt_length[:,0]) / len(reatt_length[:,0])
    numpy.save("reattachment_length_"+str(type), reatt_length)
    print "\nTime-averaged reattachment length (in step heights): ", av_length
    plot_length(type,reatt_length)

    ##### Points to generate profiles:
    xarray = numpy.array([3.06, 6.12, 10.2, 15.31])
    yarray = numpy.array([0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.91,1.92,1.93,1.94])

    ##### Only process every nth file by taking integer multiples of n:
    filelist = get_filelist(sample=50, start=0)

    ##### Call meanvelo function
    profiles, time = meanvelo(filelist, xarray, yarray)
    numpy.save("velocity_profiles_"+str(type), profiles)
    plot_meanvelo(type,profiles,xarray,yarray,time)
    pylab.show()

    print "\nAll done.\n"


if __name__ == "__main__":
    sys.exit(main())


