import numpy
import shocktube
import vtktools
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import sys

time_levels=numpy.arange(0,99)
filename_part='shocktube_'

for time_level in time_levels:
  filename_in = filename_part + str(time_level) + '.vtu'
  filename_out = filename_part + str(format(time_level,"03g")) + '.png'
  print('Processing file' , filename_in , '...',)

  vt=vtktools.vtu(filename_in)
  t=vt.GetScalarField('Time')[0]
  xyz=vt.GetLocations()
  x=xyz[:,0]

  p=vt.GetScalarField('Pressure')
  uvw=vt.GetVectorField('Velocity')
  u=uvw[:,0]
  rho=vt.GetScalarField('Density')
  ie=vt.GetScalarField('InternalEnergy')

  analytical_solution = numpy.array([shocktube.solution(xi,t) for xi in x])
  analytical_p = analytical_solution[:,0]
  analytical_u=analytical_solution[:,1]
  analytical_rho=analytical_solution[:,2]
  analytical_ie=analytical_p/analytical_rho/(shocktube.gamma-1.0)
  
  fig = plt.figure()

  pressure_subplot = fig.add_subplot(4,1,1)
  pressure_subplot.plot( x, p,'.')
  pressure_subplot.plot( x, analytical_p,'-')
  plt.axis((x[0],x[-1],0.1,1.1))
  plt.ylabel('p')

  velocity_subplot = fig.add_subplot(4,1,2)
  velocity_subplot.plot( x, u,'.')
  velocity_subplot.plot( x, analytical_u,'-')
  plt.axis((x[0],x[-1],-0.1,0.8))
  plt.ylabel('u')

  density_subplot = fig.add_subplot(4,1,3)
  density_subplot.plot( x, rho,'.')
  density_subplot.plot( x, analytical_rho,'-')
  plt.axis((x[0],x[-1],0.1,1.1))
  plt.ylabel('rho')

  internalEnergy_subplot = fig.add_subplot(4,1,4)
  internalEnergy_subplot.plot( x, ie,'.')
  internalEnergy_subplot.plot( x, analytical_ie,'-')
  plt.axis((x[0],x[-1],1.8,3.4))
  plt.ylabel('e')

  plt.savefig(filename_out, dpi=100)
  print('created file' , filename_out)

  plt.close(fig)

animation_command = ('mencoder',
                     'mf://shocktube*.png',
                     '-mf',
                     'type=png:w=800:h=600:fps=12',
                     '-ovc',
                     'lavc',
                     '-lavcopts',
                     'vcodec=mpeg4',
                     '-oac',
                     'copy',
                     '-o',
                     'shocktube.avi')

subprocess.check_call(animation_command)
