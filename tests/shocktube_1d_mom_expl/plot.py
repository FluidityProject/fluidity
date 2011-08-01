import pylab
import numpy
import shocktube
import vtktools

folders=['.']
filename='shocktube_50.vtu'

for folder in folders:
  vt=vtktools.vtu(folder + '/' + filename)
  t=vt.GetScalarField('Time')[0]
  xyz=vt.GetLocations()
  x=xyz[:,0]
  try:
    pxyz=vt.GetVectorField('DiagnosticCoordinate')
  except:
    pxyz=xyz
  px=pxyz[:,0]

  p=vt.GetScalarField('Pressure')
  uvw=vt.GetVectorField('Velocity')
  u=uvw[:,0]
  rho=vt.GetScalarField('Density')
  ie=vt.GetScalarField('InternalEnergyDensity')
  mom=vt.GetScalarField('Momentum')
  
  pylab.figure(1)
  pylab.plot( px, p,'.', label=folder)
  
  pylab.figure(2)
  pylab.plot( x, u,'.', label=folder)
  
  pylab.figure(3)
  pylab.plot( px, rho,'.', label=folder)
  
  pylab.figure(4)
  pylab.plot( px, ie,'.', label=folder)

  pylab.figure(5)
  pylab.plot( x, mom, '.', label=folder)
  
sol=numpy.array([shocktube.solution(xi,t) for xi in x])
p=sol[:,0]
u=sol[:,1]
rho=sol[:,2]
ie=p/rho/(shocktube.gamma-1.0)
mom=rho*u
mflux=mom*u+p
  
pylab.figure(1)
pylab.plot( x, p,'-', label='analytical')
pylab.title('Pressure')
pylab.legend()
  
pylab.figure(2)
pylab.plot( x, u,'-', label='analytical')
pylab.title('Velocity')
pylab.legend()

pylab.figure(3)
pylab.plot( x, rho,'-', label='analytical')
pylab.title('Density')
pylab.legend()

pylab.figure(4)
pylab.plot( x, rho*ie,'-', label='analytical')
pylab.title('Internal Energy Density')
pylab.legend()

pylab.figure(5)
pylab.plot( x, mom,'-', label='analytical')
pylab.title('Momentum')
pylab.legend()

pylab.show()
