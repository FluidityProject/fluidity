#!/usr/bin/env python
from numpy import *
import getopt, sys, os

# Global variables
filename = None
xscale = None
yscale = None

def get_options():
    global filename, xscale, yscale
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:x:y:",
                                   ["filename", "xscale", "yscale",""])
    except getopt.GetoptError:
        print("ERROR: Bad arguments!")
        print(usage())
        sys.exit(2)

    for o, a in opts:
        if o in ("-f", "--filename"):
            filename = a
        if o in ("-x", "--xscale"):
            xscale = int(float(a))
        if o in ("-y", "--yscale"):
            yscale = int(float(a))

get_options()

assert(filename != None)
assert(xscale != None)
assert(yscale != None)

r = fromfile(filename+'.node',sep=' ')
R = reshape(r[4:-1],((size(r)-5)/3,3))

#print(xscale)
#print(yscale)

Xvals = R[:,1]*xscale
Yvals = R[:,2]*yscale

nnodes = int(r[0])
dim = r[1]
assert(dim==2)

f = open(filename+'_stretched_x'+str(xscale)+'_y'+str(yscale)+'.node', 'w')
f.write(str(int(nnodes)))
f.write(' 2 0 0\n')
for i in range(nnodes):
    s = "%d %18.18f %18.18f \n" % (i+1, Xvals[i], Yvals[i])
    f.write(s)
f.write("#produced by stretch.py from "+filename)
f.close()

os.system('cp '+filename+'.ele '+filename+'_stretched_x'+str(xscale)+'_y'+str(yscale)+'.ele')
os.system('cp '+filename+'.edge '+filename+'_stretched_x'+str(xscale)+'_y'+str(yscale)+'.edge')

