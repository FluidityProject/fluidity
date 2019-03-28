#!/usr/bin/env python

from __future__ import print_function

import glob
import sys
import numpy

def inflow_data(plane):
  # Extract mean velocity and Reynolds stresses from ERCOFTAC data at given plane.
  # Ignore file header
  ignore=3
  ercoftac = open('x.'+str(plane), 'r')
  print("reading in data from file: x."+str(plane))
  for line in range(ignore):
    ercoftac.readline()
  
  y=[];U=[];uu=[];vv=[];uv=[]
  for line in ercoftac:
    y.append(float(line.split(  )[0]))
    U.append(float(line.split(  )[1]))
    uu.append(float(line.split(  )[3]))
    vv.append(float(line.split(  )[4]))
    uv.append(-10.0*float(line.split(  )[6]))

  # Write data out
  datafile = open('BFS-SEM-ERCOFTAC'+str(plane)+'-list.dat', 'w')
  datafile.write('y/h\n'+str(y)+'\nU/U0\n'+str(U)+"\nu'/U0\n"+str(uu)+"\nv'/U0\n"+str(vv)+"\nw'/U0\n"+str(uv)+'\n')
  datafile = open('BFS-SEM-ERCOFTAC'+str(plane)+'-table.dat', 'w')
  datafile.write("y/h U/U0 u'u'/U0^2 v'v'/U0^2 -10*u'v'/U0^2\n")
  for i in range(len(y)):
    datafile.write(str(y[i])+' '+str(U[i])+' '+str(uu[i])+' '+str(vv[i])+' '+str(uv[i])+'\n')
  ercoftac.close(); datafile.close()
  return

#########################################################################

def main():

  plane = sys.argv[1]
  inflow_data(plane)

if __name__ == "__main__":
  sys.exit(main())

