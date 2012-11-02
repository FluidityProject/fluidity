import pylab
import numpy
import vtktools
import sys
params = {'text.fontsize': 11,
         'legend.fontsize': 9,
         'xtick.labelsize': 11,
         'ytick.labelsize': 11,
         'lines.markersize': 6,
         'lines.linewidth': 1.5,
         'axes.titlesize': 'medium'}
pylab.rcParams.update(params)

# Maximum number of VTU files
N = 40

# Dimensions
H = 0.1
U = 0.6
h = 2*H
x = 0.5*H

Y = numpy.linspace(0.0, h, 200)
probed_u = numpy.zeros(len(Y))
u_bar = numpy.zeros(len(Y))

for n in range(N, N+1):
   filename='flow_past_a_cube_' + str(n) + '.pvtu'
   vt=vtktools.vtu(filename)

   # Time and coordinates
   t = vt.GetScalarField('Time')[0]
   xyz = vt.GetLocations()

   for i in range(0,len(Y)): 
      data = vtktools.vtu.ProbeData(vt, numpy.array([[x, Y[i], 0]]), 'Velocity')
      #print data
      probed_u[i] = data[0][0]
      u_bar[i] = u_bar[i] + probed_u[i]

u_bar = u_bar
# Plot the numerical results
pylab.plot(u_bar/U, Y/H, '-b', label="Fluidity")

pylab.legend()
pylab.xlabel("Normalised mean velocity (u/U)")
pylab.ylabel("Normalised height (y/H)")
pylab.grid("on")

pylab.show()
#pylab.savefig('flow_past_a_cube.png')
