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
N = 17

# Dimensions
H = 1.0
U = 32.0
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
pylab.plot(u_bar/U, Y/H, '-b', label="Fluidity (x/H = 0.5)")

experimental_y = [1.0, 1.0194, 1.0485, 1.077669, 1.097, 1.1067, 1.116,
1.165, 1.1747, 1.20388, 1.2427, 1.328855, 1.38778, 1.466, 1.52427, 1.5825,
1.64077, 1.7087, 1.76699, 1.82524, 1.8737, 1.97087]
experimental_u_bar = [-0.023, -0.534, -0.534, -0.499, -0.452, -0.394, -0.37,
0.06, 0.35, 0.7, 1.119, 1.28155, 1.3495, 1.3426, 1.3084, 1.274,
1.24, 1.2175, 1.1833, 1.1375, 1.09155, 0.96477]
pylab.plot(experimental_u_bar, experimental_y, 'ro', label="Experimental (x/H = 0.5)")

pylab.legend(loc=4)
pylab.xlabel("Normalised mean velocity (u/U)")
pylab.ylabel("Normalised height (y/H)")
pylab.grid("on")

#pylab.show()
pylab.savefig('flow_past_a_cube.png')
