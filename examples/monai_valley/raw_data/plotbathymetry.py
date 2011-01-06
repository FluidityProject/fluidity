#!/usr/bin/python
import matplotlib as m
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import sys
import pdb
import os
import random

def file_len(full_path):
        """ Count number of lines in a file."""
        f = open(full_path)
        nr_of_lines = sum(1 for line in f)
        f.close()
        return nr_of_lines

# Program Options 
xres=200 # Number of Grid points in x direction
yres=200 # Number of Grid points in y direction
#decreaseDataSet=100000 # Set 0 to load all data
decreaseDataSet=0
maxZValue=100000 # This is the maximum z value allowed. All z values higher than that will be ste to MaxZValue


#scatterColorMap=plt.cm.autumn
#define own colormap from white to black
cdict = {
  'red'  :  ((0., 1.0, 1.0), (1., 0.2, 0.2)),
  'green':  ((0., 1.0, 1.0), (1., 0.2, 0.2)),
  'blue' :  ((0., 1.0, 1.0), (1., 0.2, 0.2))
}

scatterColorMap = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)


""" Main Function """
if not len(sys.argv)>=2:
        print "Usage: plotbathymetry.py Bathymetry.grd [Title]"        
        exit()
if len(sys.argv)==3:
        title=sys.argv[2]
else:
        title=sys.argv[1]

#Check if file is NetCFD file and converge to XYZ if necessary
s = os.popen("file " + sys.argv[1] + "| grep NetCDF | wc -l").read()
tempfile=sys.argv[1]
if s=="1\n":
        print "Found NetCDF file. Start grd2xyz"
        tempfile="/tmp/tempXYZ_"+str(random.randint(0,100))+".dat"
        os.popen("grd2xyz -V -S " + sys.argv[1] + " >> " + tempfile)
        sys.argv[1]=tempfile

print "Load file " + tempfile + "\n"
inp = open(sys.argv[1])
arr = [[],[],[]]
nr_of_lines = file_len(tempfile)
print str(nr_of_lines) +  " Datapoints found."


# read line into array
counter=0
if decreaseDataSet==0 or nr_of_lines<decreaseDataSet:
        step=1
else:
        step=nr_of_lines/decreaseDataSet
        print "Decrease dataset to " + str(decreaseDataSet) + " points."

for line in inp.readlines():
        if counter%step==0:
                # loop over the elemets, split by whitespace
                dat=line.split()
                arr[0].append(float(dat[0]))
                arr[1].append(float(dat[1]))
                arr[2].append(min(maxZValue,float(dat[2])))
        counter=counter+1

print "Extension in X-direction: From " + str(min(arr[0])) + " to " + str(max(arr[0]))
print "Extension in Y-direction: From " + str(min(arr[1])) + " to " + str(max(arr[1]))
print "Extension in Z-direction: From " + str(min(arr[2])) + " to " + str(max(arr[2]))

#Define Grid:
xi=np.linspace(min(arr[0]),max(arr[0]),xres)
yi=np.linspace(min(arr[1]),max(arr[1]),yres)
# grid the data.
zi = griddata(arr[0],arr[1],arr[2],xi,yi)


# contour the gridded data, plotting dots at the randomly spaced data points.
fig = plt.figure()
subplt = fig.add_subplot(111, xlabel='[m]', ylabel='[m]')
subplt.contour(xi,yi,zi,30,linewidths=0.5,colors='k')
#plt.clabel(CS, inline=1, fontsize=8)
contf=subplt.contourf(xi,yi,zi*100,15,cmap=scatterColorMap)

# plot position of gauge stations
subplt.plot(4.521, 1.196, 'o',color='k')
subplt.text(3.821, 1.246, "Gauge 3")
subplt.plot(4.521, 1.696, 'o',color='k')
subplt.text(3.821, 1.746, "Gauge 2")
subplt.plot(4.521, 2.196, 'o',color='k')
subplt.text(3.821, 2.246, "Gauge 1")


colbar=fig.colorbar(contf) # draw colorbar
colbar.set_label('[cm]')
#plt.clim(-0.12,0.0)
plt.xlim((0, 5.448))
plt.ylim((0, 3.402))

#plt.savefig( "BODC_plot.pdf", format='pdf' ) 

#plt.title(title)
plt.show()

if s=="1\n":
        os.popen("rm " + tempfile)

