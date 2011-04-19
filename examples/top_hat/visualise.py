#!/usr/bin/env python
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
vt1=vtktools.vtu(sys.argv[1])
meshsize = float(sys.argv[1].split('/')[0].split('_')[-1])
print meshsize
# returns an array containing the 3-d locations of the nodes in the mesh
# of these only the first coordinate is used in 1D, i.e xyz[:,0]
xyz1=vt1.GetLocations()
x1=xyz1[:,0]

# get the values of the "Tracer" field at the nodes
T1=vt1.GetScalarField("Tracer")

# analytic solution
A=[0., 0., 1., 1., 0., 0.]
#x=[1.8, 2.0, 2.0, 2.25, 2.25, 2.5]
x=[1.,1.25,1.25,1.75,1.75,2.]

# Helmholtz Smoothed Scalar Algorithm
h1=vt1.GetScalarField("TracerHS1")
h2=vt1.GetScalarField("TracerHS2")
h3=vt1.GetScalarField("TracerHS3")
h4=vt1.GetScalarField("TracerHS4")
h5=vt1.GetScalarField("TracerHS5")
# Helmholtz Anisotropic Smoothed Scalar Algorithm
ha1=vt1.GetScalarField("TracerHSA1")
ha2=vt1.GetScalarField("TracerHSA2")
ha3=vt1.GetScalarField("TracerHSA3")
ha4=vt1.GetScalarField("TracerHSA4")
# Lumped Mass Smoothed Scalar Algorithm
ml1=vt1.GetScalarField("TracerLumped")
ml2=vt1.GetScalarField("TracerLumped2x")
ml3=vt1.GetScalarField("TracerLumped3x")

# Diff
dh1=T1-h1
dh2=T1-h2
dh3=T1-h3
dh4=T1-h4
dh5=T1-h5
dml=T1-ml1
dha1=T1-ha1

# element size field
#e=vt1.GetScalarField("elements")

# plot Tracers against x, black lines, grey round markers
pylab.plot(x1,T1, 'blue', linestyle='solid', marker = 'o', markersize=8, markeredgecolor='blue')
pylab.plot(x1,h1, 'green', linestyle='solid', marker = '+', markersize=8, markeredgecolor='green')
pylab.plot(x1,h2, 'green', linestyle='solid', marker = 'x', markersize=8, markeredgecolor='green')
pylab.plot(x1,h3, 'green', linestyle='solid', marker = 'o', markersize=8, markeredgecolor='green')
pylab.plot(x1,h4, 'green', linestyle='solid', marker = '*', markersize=8, markeredgecolor='green')
pylab.plot(x1,h5, 'green', linestyle='solid', marker = 's', markersize=8, markeredgecolor='green')
pylab.plot(x1,ha1, 'blue', linestyle='solid', marker = 's', markersize=8, markeredgecolor='blue')
pylab.plot(x1,ha2, 'blue', linestyle='solid', marker = 'x', markersize=8, markeredgecolor='blue')
pylab.plot(x1,ha3, 'blue', linestyle='solid', marker = '+', markersize=8, markeredgecolor='blue')
pylab.plot(x1,ha4, 'blue', linestyle='solid', marker = '*', markersize=8, markeredgecolor='blue')
pylab.plot(x1,ml1, 'red', linestyle='solid', marker = 'x', markersize=8, markeredgecolor='red')
pylab.plot(x1,ml2, 'red', linestyle='solid', marker = '+', markersize=8, markeredgecolor='red')
pylab.plot(x1,ml3, 'red', linestyle='solid', marker = 's', markersize=8, markeredgecolor='red')
pylab.plot(x,A, 'black', linestyle='solid', marker = 'o', markerfacecolor='black', markersize=8, markeredgecolor='black')
#pylab.plot(x1,e*0.1+0.5, 'black', linestyle='solid')
#pylab.plot(x1,dh1-0.5, 'purple')
#pylab.plot(x1,dh2-0.5, 'brown')
#pylab.plot(x1,dh3-0.5, 'orange')
#pylab.plot(x1,dh4-0.5, 'blue')
#pylab.plot(x1,dh5-0.5, 'green')
#pylab.plot(x1,dha1-0.5, 'black')
#pylab.plot(x1,dml-0.5, 'red')

# set axes and axes' labels
pylab.axis([0.75,2.25,-1.0,1.6])
pylab.xlabel('').set_size('xx-large')
pylab.xlabel('$x$')
pylab.ylabel('Tracer')
#pylab.legend(("Tracer","H1:alpha=0.00025","H2:alpha=0.0005","H3:alpha=0.001","H4:alpha=0.002","H5:alpha=0.004","HA1:alpha=1","HA2:alpha=2","HA3:alpha=3","HA4:alpha=4","ML1","ML2","ML3","meshsize","analytic","T-H1","T-H2","T-H3","T-H4","T-H5","T-H3a","T-ML1"),loc='best')
pylab.legend(("Tracer","H1:alpha=0.00025","H2:alpha=0.0005","H3:alpha=0.001","H4:alpha=0.002","H5:alpha=0.004","HA1:alpha=1","HA2:alpha=2","HA3:alpha=3","HA4:alpha=4","ML1","ML2","ML3","analytic"),loc='best')
#pylab.savefig('top_hat_cg_100_adv_fixed_'+str(meshsize)+'.png')
# show the figure
pylab.show()
print '\n\nAll done'

