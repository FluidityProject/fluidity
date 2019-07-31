#!/usr/bin/env python3
import sys
import pylab
try:
  import vtktools
except ImportError:
  sys.stderr.write("""You need to add the path to <FLUIDITY_SOURCE_LOCATION>/python/
to your PYTHONPATH environment variable, e.g.:
  export PYTHONPATH=$PYTHONPATH:$HOME/svn/fluidity/python/\n""")
  sys.exit(1)
  
# open one of the vtus output by fluidity
vt=vtktools.vtu(sys.argv[1])
# returns an array containing the 3-d locations of the nodes in the mesh
# of these only the first coordinate is used in 1D, i.e xyz[:,0]
xyz=vt.GetLocations()
x=xyz[:,0]
# get the values of the "Tracer" field at the nodes
T=vt.GetScalarField("Tracer")

# plot Tracer T against x, black lines, grey round markers
pylab.plot(x,T, 'black', markerfacecolor='grey', marker='o')
# set axes and axes' labels
pylab.axis([0.0,3.0,-0.1,1.1])
pylab.xlabel('').set_size('xx-large')
pylab.xlabel('$x$')
pylab.ylabel('Tracer')
# show the figure
pylab.show()
