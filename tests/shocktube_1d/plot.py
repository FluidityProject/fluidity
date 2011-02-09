import pylab
import numpy
import shock_tube
import vtktools

filename='shocktube_50.vtu'

vt=vtktools.vtu(filename)
t=vt.GetScalarField('Time')[0]
xyz=vt.GetLocations()
x=xyz[:,0]

sol=numpy.array([shock_tube.solution(xi,t) for xi in x])
solp=sol[:,0]
solu=sol[:,1]
  
p=vt.GetScalarField('Pressure')
uvw=vt.GetVectorField('Velocity')
u=uvw[:,0]

pylab.plot( x, solp )
pylab.plot( x, p )
pylab.show()

pylab.plot( x, solu )
pylab.plot( x, u )
pylab.show()