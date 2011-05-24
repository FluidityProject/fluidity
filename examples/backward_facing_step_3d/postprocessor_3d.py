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
def reatt_length(filelist, yarray):

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
    ##### points near bottom surface, 0 < x < 25
    x2array=[]; pts=[]; no_pts = 52; offset = 0.01
    x = 0.0
    for i in range(1, no_pts):
      x2array.append(x)
      for j in range(len(yarray)):
        pts.append((x, yarray[j], offset))
      x += 0.5

    x2array = numpy.array(x2array)
    pts = numpy.array(pts)

    ##### Get x-velocity on bottom boundary
    uvw = datafile.ProbeData(pts, "Velocity")
    u = uvw[:,0]
    u = u.reshape([x2array.size,yarray.size])

    ##### Find all potential reattachment points:
    points = []
    for j in range(len(u[0,:])):
      for i in range(len(u[:,0])-1):
        ##### Hack to ignore division by zero entries in u.
        ##### All u should be nonzero away from boundary!
        if((u[i,j] / u[i+1,j]) < 0. and not numpy.isinf(u[i,j] / u[i+1,j])):
          ##### interpolate between nodes
          p = x2array[i] + (x2array[i+1]-x2array[i]) * (0.0-u[i,j]) / (u[i+1,j]-u[i,j])
          ##### Ignore spurious corner points
          if(p>0.1):
            points.append(p)

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
def meanvelo(filelist,xarray,yarray,zarray):

  print "\nRunning velocity profile script on files at times...\n"
  ##### check for no files
  if (len(filelist) < 0):
    print "No files!"
    sys.exit(1)

  ##### create array of points
  pts=[]
  for i in range(len(xarray)):
    for j in range(len(yarray)):
      for k in range(len(zarray)):
        pts.append([xarray[i], yarray[j], zarray[k]])
  pts=numpy.array(pts)

  ##### Create output array of correct shape
  profiles=numpy.zeros([len(filelist), xarray.size, zarray.size], float)
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
      u = uvw[:,0]/umax
      u = u.reshape([xarray.size,yarray.size,zarray.size])

      ##### Spanwise averaging
      usum = numpy.zeros([xarray.size,zarray.size],float)
      usum = numpy.array(usum)
      for i in range(len(yarray)):
        uav = u[:,i,:]
        uav = numpy.array(uav)
        usum += uav
      usum = usum / len(yarray)
      profiles[filecount,:,:] = usum
      filecount += 1

  print "\n...Finished extracting data.\n"
  return profiles, time

#########################################################################

def plot_length(Re,type,mesh,reattachment_length,av_length):
  ##### Plot time series of reattachment length using pylab(matplotlib)
  avg = numpy.zeros([len(reattachment_length[:,1])])
  avg[:] = av_length
  Lemoinkim = numpy.zeros([len(reattachment_length[:,1])])
  Lemoinkim[:]=6.28

  plot1 = pylab.figure()
  pylab.title("Time series of reattachment length: Re="+str(Re)+", "+str(type)+"-Re BCs, "+str(mesh)+" mesh")
  pylab.xlabel('Time (s)')
  pylab.ylabel('Reattachment Length (L/h)')
  pylab.plot(reattachment_length[:,1], reattachment_length[:,0], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="solid")
  pylab.plot(reattachment_length[:,1], avg, linestyle="dashed")
  pylab.plot(reattachment_length[:,1], Lemoinkim, linestyle="dashed")
  pylab.legend(("length (step heights)","average length","Le-Moin-Kim DNS average"), loc="lower right")
  pylab.savefig("../reatt_len_3D_"+str(Re)+"_"+str(type)+"_"+str(mesh)+".pdf")
  return

def plot_meanvelo(Re,type,mesh,profiles,xarray,yarray,zarray,time):
  ##### Plot velocity profiles at different points behind step, and at 3 times using pylab(matplotlib)
  plot1 = pylab.figure(figsize = (16.5, 8.5))
  pylab.suptitle("Evolution of U-velocity: Re="+str(Re)+", "+str(type)+"-Re BCs, "+str(mesh)+" mesh", fontsize=20)

  size = 15

  ax = pylab.subplot(141)
  shift=0.0
  leg_end = []
  for i in range(len(time)):
    ax.plot(profiles[i,0,:]+shift,zarray, linestyle="solid")
    shift+=0.0
    leg_end.append("%.1f secs"%time[i])
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
    bx.plot(profiles[i,1,:]+shift,zarray, linestyle="solid")
    shift+=0.0
  bx.set_title('(a) x/h='+str(xarray[1]), fontsize=16)
  #bx.grid("True")
  for tick in bx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(bx.get_yticklabels(), visible=False)

  cx = pylab.subplot(143, sharex=ax, sharey=ax)
  shift=0.0
  for i in range(len(time)):
    cx.plot(profiles[i,2,:]+shift,zarray, linestyle="solid")
    shift+=0.0
  cx.set_title('(a) x/h='+str(xarray[2]), fontsize=16)
  #bx.grid("True")
  for tick in cx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(cx.get_yticklabels(), visible=False)

  dx = pylab.subplot(144, sharex=ax, sharey=ax)
  shift=0.0
  for i in range(len(time)):
    dx.plot(profiles[i,3,:]+shift,zarray, linestyle="solid")
    shift+=0.0
  dx.set_title('(a) x/h='+str(xarray[3]), fontsize=16)
  #bx.grid("True")
  for tick in dx.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  pylab.setp(dx.get_yticklabels(), visible=False)

  pylab.axis([-0.25, 1., 0., 3.])
  bx.set_xlabel('Normalised U-velocity (U/Umax)', fontsize=24)
  ax.set_ylabel('z/h', fontsize=24)

  pylab.savefig("../velo_profiles_3d"+str(Re)+"_"+str(type)+"_"+str(mesh)+".pdf")
  return

#########################################################################

def main():
    ##### Which run is being processed?
    Re = sys.argv[1]
    type = sys.argv[2]
    mesh = sys.argv[3]
    print "Re, bc type, mesh: ", Re, type, mesh

    ##### Only process every nth file by taking integer multiples of n:
    filelist = get_filelist(sample=1, start=0)

    ##### Points to generate profiles:
    xarray = numpy.array([4.0, 6.0, 10.0, 19.0])
    yarray = numpy.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9])
    zarray = numpy.array([0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])

    ##### Call reattachment_length function
    reattachment_length = numpy.array(reatt_length(filelist, yarray))
    av_length = sum(reattachment_length[:,0]) / len(reattachment_length[:,0])
    numpy.save("reattachment_length_3D", reattachment_length)
    print "\nTime-averaged reattachment length (in step heights): ", av_length
    plot_length(Re,type,mesh,reattachment_length,av_length)

    ##### Call meanvelo function
    #profiles, time = meanvelo(filelist, xarray, yarray, zarray)
    #numpy.save("velo_profiles_3d_parallel"+str(Re)+"_"+str(mesh), profiles)
    #print "Showing plot of velocity profiles.\nTo continue script, close plot window."
    #plot_meanvelo(Re,type,mesh,profiles,xarray,yarray,zarray,time)
    #pylab.show()

    print "\nAll done.\n"

if __name__ == "__main__":
    sys.exit(main())

