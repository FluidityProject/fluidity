#!/usr/bin/env python

from optparse import OptionParser
import sys
import os.path
import spheretools
import shutil
from numpy import zeros

optparser=OptionParser(usage='usage: %prog <filename>',
                       description="""This takes a triangle .node file """+
                       """with cartesian coordinates and produces a .node """+
                       """with spherical polar coordinates. """+
                       """The .ele and .edge files are also copied. """+
                       """The new files have the prefix sp_ """+
                       """(meaning spherical polar).""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

if argv[0][-5:]!=".node":
    sys.stderr.write("Filename must end in .node\n")
    optparser.print_help()
    sys.exit(1)
    
basename=os.path.basename(argv[0][:-5])

# get map and set up positions
map=spheretools.cart2polar()
X=zeros(3)

# open files
nodefile=file(argv[0],"r")
new_nodefile=file("sp_"+basename+".node","w")

# put header unchanged into new .node file
info=nodefile.readline()
new_nodefile.write(info)

# apply mapping and write new coordinates to new .node file
for i in range(int(info.split()[0])):
    line=nodefile.readline().split()
    for j in range(3):
        X[j]=float(line[j+1])
    ([theta,phi])=map(X,0)
    new_nodefile.write(line[0]+" "+str(theta)+" "+str(phi)+"\n")

# copy over .ele and .edge files and give an error if either doesn't exist.
try:
    shutil.copyfile(basename+".ele","sp_"+basename+".ele")
except IOError as err:
    print err

try:
    shutil.copyfile(basename+".edge","sp_"+basename+".edge")
except IOError as err:
    print err
