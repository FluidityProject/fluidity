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
N = 240

# Dimensions
H = 0.04
U = 0.535
h = 7*H

experimental_X = [0.981, 1.0888, 1.233, 1.3397, 1.5188,
1.6256, 1.73, 1.83799, 2.02, 2.128, 2.1989,
2.34, 2.489, 2.598, 2.7075, 2.852, 2.998,
3.217, 3.472, 3.726, 3.982, 4.201, 4.456,
4.7118, 5.0, 5.22, 5.478, 5.733, 5.952, 6.389,
6.75456, 7.156, 7.557, 7.9589]

experimental_u_bar = [-0.18, -0.132, -0.091589, -0.006679, 0.0856,
0.159477, 0.288667, 0.34, 0.3515, 0.3884, 0.451,
0.48, 0.4955, 0.503, 0.5141, 0.5547, 0.5474,
0.5401, 0.5623, 0.5956, 0.5847, 0.581, 0.5959,
0.596, 0.607, 0.5962, 0.5963, 0.6, 0.5964, 0.62618,
0.6263, 0.62646, 0.63399, 0.626759]

#Y = numpy.linspace(-h, h, 400)
X = numpy.linspace(-5*H, 20*H - H/2, 400)

pylab.figure()

u_bar = numpy.zeros(len(X))

filename='les/flow_past_a_square_' + str(N) + '.pvtu'
vt=vtktools.vtu(filename)

# Time and coordinates
t = vt.GetScalarField('Time')[0]
xyz = vt.GetLocations()

for i in range(0,len(X)):
   data = vtktools.vtu.ProbeData(vt, numpy.array([[X[i], 0, 0]]), 'TimeAveragedVelocity30')
   u_bar[i] = data[0][0]

# Plot the numerical results
pylab.plot(X/H, u_bar/U, '-b', label="Fluidity (LES)")

filename='noles/flow_past_a_square_' + str(N) + '.pvtu'
vt=vtktools.vtu(filename)

# Time and coordinates
t = vt.GetScalarField('Time')[0]
xyz = vt.GetLocations()

for i in range(0,len(X)):
   data = vtktools.vtu.ProbeData(vt, numpy.array([[X[i], 0, 0]]), 'TimeAveragedVelocity30')
   u_bar[i] = data[0][0]

# Plot the numerical results
pylab.plot(X/H, u_bar/U, '-g', label="Fluidity (no LES)")

pylab.plot(experimental_X, experimental_u_bar, 'xr', label="Experimental (Lyn and Rodi (1994))")
pylab.legend(loc=4)
pylab.xlabel("Normalised x (x/H)")
pylab.ylabel("Normalised time-averaged velocity (u/U)")
pylab.grid("on")

#pylab.show()
pylab.savefig('flow_past_a_square.png')
