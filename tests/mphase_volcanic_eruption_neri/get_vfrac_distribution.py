from fluidity_tools import stat_parser
import vtktools
import numpy
import pylab

filename = "mphase_plinian_eruption_neri_checkpoint_"+ str(211) +".pvtu"
vt=vtktools.vtu(filename)

t = vt.GetScalarField('Air::Time')[0]
xyz = vt.GetLocations()

X = numpy.linspace(0, 5000, 500) # 5 km radial
vfrac = numpy.zeros(len(X))
for i in range(0,len(X)):
   data = vtktools.vtu.ProbeData(vt, numpy.array([[X[i], 0, 0]]), 'Tephra::PhaseVolumeFraction')
   vfrac[i] = data[0]
   
pylab.plot(X, vfrac, '-b', label="Radial vfrac (Fluidity)")

pylab.legend(loc=4)
pylab.xlabel("x")
pylab.ylabel("Tephra vfrac")
pylab.grid("on")

#pylab.show()
pylab.savefig('radial_vfrac.png')
