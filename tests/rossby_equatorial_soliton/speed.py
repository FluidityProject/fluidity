#!/usr/bin/env python

## Usage: ./speed.py filename_[0-9]*.vtu

import sys
import os
import vtktools
import numpy
import scipy.stats
import matplotlib
import pylab

def run(filelist, timdum=1.0):
  xplus = []
  xminus = []
  fsplus = []
  fsminus = []
  N = 0
  order = []

  for file in filelist:
##    print file
    try:
      os.stat(file)
    except:
      print "No such file: %s" % file
      sys.exit(1)

    num = int(file.split(".vtu")[0].split('_')[-1])
    if num<41:
##      print num
      N = N+1
      order.append(num)

      ug=vtktools.vtu(file)
##      ug.GetFieldNames()
      pos=ug.GetLocations()
      fs=ug.GetScalarField('FreeSurface')

      l1 = numpy.argsort(fs)
      sortedfs=fs[l1]
      sortedpos=pos[l1]
      sortedz = sortedpos[:,2]

      l2 = numpy.where(sortedz>0.5)
      sortedfs2 = sortedfs[l2]
      sortedpos2 = sortedpos[l2]
      sortedy = sortedpos2[:,1]
    
      l3 = numpy.where(sortedy>=0.0)
      sortedfs3 = sortedfs2[l3]
      sortedpos3 = sortedpos2[l3]
      
      l4 = numpy.where(sortedy<0.0)
      sortedfs4 = sortedfs2[l4]
      sortedpos4 = sortedpos2[l4]

      xplus.append(sortedpos3[-1,0])
      xminus.append(sortedpos4[-1,0])
      fsplus.append(sortedfs3[-1])
      fsminus.append(sortedfs4[-1])


  sortorder = numpy.argsort(order)
  
  xplusa = numpy.array(xplus)
  xminusa = numpy.array(xminus)
  fsplusa = numpy.array(fsplus)
  fsminusa = numpy.array(fsminus)
  
  xplussort = xplusa[sortorder]
  xminussort = xminusa[sortorder]
  fsplussort = fsplusa[sortorder] 
  fsminussort = fsminusa[sortorder]    

##  xplussortreverse = xplussort[::-1]
##  xminussortreverse = xminussort[::-1]

  results = []
  for x in [xplussort, xminussort]:
    t = [(i) * timdum for i in range(len(x))]
    r = scipy.stats.linregress(t, x)[0]
    results.append(r)

  results.append(N)
  results.append(t)
  results.append(xplussort)
  results.append(xminussort)
  results.append(fsplussort)
  results.append(fsminussort)

  return results
      
      
(minspeed, maxspeed, N, t, xplus, xminus, fsplus, fsminus) = run(sys.argv[1:])

print "\nNumber of files used = %s" % N

print "\nNorth speed = %s" % minspeed
print "South speed = %s" % maxspeed
print "Average speed = %s\n" % (0.5*(minspeed+maxspeed))

print "Max free surface (North) at final time = %s" % fsplus[-1]
print "Max free surface (South) at final time = %s" % fsminus[-1]
print "Max free surface (Average) at final time = %s\n" % (0.5*(fsplus[-1]+fsminus[-1]))

print t
print xplus
print xminus
print fsplus
print fsminus

if(0):
  pylab.figure(1)
  pylab.subplot(221)
  pylab.plot(t, xplus, 'b', t, maxspeed*numpy.array(t), 'k--')
  pylab.xlabel('time')
  pylab.ylabel('x location of North peak')
  pylab.subplot(223)
  pylab.plot(t, xminus, 'r--', t, minspeed*numpy.array(t), 'k--')
  pylab.xlabel('time')
  pylab.ylabel('x location of South peak')
  pylab.subplot(222)
  pylab.plot(t, fsplus, 'r--')
  pylab.xlabel('time')
  pylab.ylabel('North amplitude max')
  pylab.subplot(224)
  pylab.plot(t, fsminus, 'r--')
  pylab.xlabel('time')
  pylab.ylabel('South amplitude max')
  pylab.show()

