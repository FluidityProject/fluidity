#!/usr/bin/env python

import vtk
import glob
import sys
import os
import vtktools
import operator
import numpy
import pylab
import fluidity_tools

# first extract the water gauge data from the vtus
def get_water_depths(filelist, xarray, delta):
  results = []
  for f in filelist:
    try:
      os.stat(f)
    except:
      print("No such file: %s" % f)
      sys.exit(1)
    
    y = numpy.arange(delta/2.0,2.0+delta/2.0,delta)[:,numpy.newaxis]
    
    num = int(f.split(".vtu")[0].split('_')[-1])
    vtu = vtktools.vtu(f)
    for name in vtu.GetFieldNames(): 
      if name.endswith("Time"): 
        time = max(vtu.GetScalarRange(name))
        break
    waterdepths = []
    waterdepths.append(num)
    waterdepths.append(time)
    for x in range(len(xarray)):
      coordinates = numpy.concatenate((numpy.ones((len(y),1))*xarray[x], y, numpy.zeros((len(y),1))),1)
      waterdepths.append(sum(vtu.ProbeData(coordinates, "Water::MaterialVolumeFraction"))[0]*delta)
    
    results.append(waterdepths)
  
  results.sort(key=operator.itemgetter(1))
  results = numpy.array(results)
  return results

xarray = [2.725, 2.228]

filelist = glob.glob("water_collapse_[0-9]*.vtu")
results = get_water_depths(filelist, xarray, 0.01)
numpy.save("water_depths", results)


# now let's plot the data
warray = ["H1", "H2"]
parray = ["P2"]

# first the water gauges
for x in range(len(xarray)):
  pylab.figure(x)
  pylab.title(warray[x]+" water gauge at "+str(xarray[x])+"m")
  pylab.xlabel('Time (s)')
  pylab.ylabel('Water Depth (m)')
  if((warray[x]=="H1") or (warray[x]=="H2")):
    experiment = numpy.load(warray[x]+".npy")
    pylab.plot(experiment[:,0], experiment[:,1], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="None")

time = results[:,1]
for x in range(len(xarray)):
  pylab.figure(x)
  pylab.plot(time, results[:,2+x], color='black', linestyle="dashed")

for x in range(len(xarray)):
  pylab.figure(x)
  pylab.axis([0.0, 2.5, 0.0, 0.5])
  pylab.legend(("Experiment", "Fluidity"), loc="upper left")
  pylab.savefig("water_gauge_"+warray[x]+".pdf")

# then the pressure gauges - this takes it data from the detectors so no
# need for extraction from the vtus
for p in range(len(parray)):
  pylab.figure(p+len(xarray))
  pylab.title(parray[p]+' pressure gauge')
  pylab.xlabel('Time (s)')
  pylab.ylabel('Pressure (Pa)')
  if(parray[p]=="P2"):
    experiment = numpy.load(parray[p]+".npy")
    pylab.plot(experiment[:,0], experiment[:,1], marker = 'o', markerfacecolor='white', markersize=6, markeredgecolor='black', linestyle="None")

results = fluidity_tools.stat_parser("water_collapse.detectors")
time = results["ElapsedTime"]["value"]
for p in range(len(parray)):
  pylab.figure(p+len(xarray))
  data = results["Water"]["Pressure"][parray[p]]
  pylab.plot(time, data, color='black', linestyle="dashed")

for p in range(len(parray)):
  pylab.figure(p+len(xarray))
  pylab.axis([0.0, 2.5, -2000., 12000.])
  pylab.legend(("Experiment", "Fluidity"), loc="upper left")
  pylab.savefig("pressure_gauge_"+parray[p]+".pdf")


pylab.show()
