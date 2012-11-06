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
N = 10

# Dimensions
H = 1.0
U = 32.0
h = 2*H

# Data from Krajnovic and Davidson (2002)
# At x/H = 0.5
experimental_y_05 = [1.0, 1.0194, 1.0485, 1.077669, 1.097, 1.1067, 1.116,
1.165, 1.1747, 1.20388, 1.2427, 1.328855, 1.38778, 1.466, 1.52427, 1.5825,
1.64077, 1.7087, 1.76699, 1.82524, 1.8737, 1.97087]
experimental_u_bar_05 = [-0.023, -0.534, -0.534, -0.499, -0.452, -0.394, -0.37,
0.06, 0.35, 0.7, 1.119, 1.28155, 1.3495, 1.3426, 1.3084, 1.274,
1.24, 1.2175, 1.1833, 1.1375, 1.09155, 0.96477]

# At x/H = 1.0
experimental_y_10 = [1.0194, 1.04854, 1.05825, 1.08737, 1.10679, 1.1262,
1.165, 1.20388, 1.2427, 1.28155, 1.349514, 1.40776699, 1.466,
1.52427, 1.582524, 1.650485, 1.7087, 1.76699, 1.82524, 1.88349,
1.951456, 1.98058]
experimental_u_bar_10 = [-0.116, -0.15116, -0.1046, -0.093, -0.0465, 0.03488,
0.244186, 0.558139, 0.79069, 0.98837, 1.18604, 1.290697, 1.31395,
1.3023, 1.2790697, 1.244186, 1.2093, 1.174418, 1.15116, 1.10465,
1.03488, 0.89534]

# At x/H = 2.0
experimental_y_20 = [0.03329, 0.07206, 0.1108, 0.1498, 0.189,
0.2278, 0.3158, 0.3938, 0.4716, 0.5587, 0.62599, 0.7124,
0.77833, 0.91037, 0.948, 0.9856, 1.023, 1.0987, 1.184,
1.26088, 1.3375, 1.415, 1.493, 1.581, 1.66, 1.7493, 1.878,
1.929, 1.95096]
experimental_u_bar_20 = [-0.39589, -0.37309, -0.35028, -0.339, -0.339566,
-0.3167, -0.306, -0.2838, -0.2498, -0.1927, -0.1237, -0.03169,
0.107069, 0.37297, 0.442, 0.53489, 0.6158, 0.76613, 0.893,
0.9968, 1.0889, 1.134, 1.15689, 1.144, 1.12, 1.08415, 1.001,
0.9192, 0.8027]

# At x/H = 4.0
experimental_y_40 = [0.0307, 0.05057, 0.1776, 0.3336, 0.421, 0.4993, 0.645, 
0.7329, 0.8108, 0.8885, 1.015, 1.132, 1.22, 1.298, 1.453, 1.53,
1.61857, 1.7057, 1.783, 1.8598]
experimental_u_bar_40 = [0.4648, 0.511, 0.6026, 0.6707, 0.704, 0.7384, 0.7716,
0.79388, 0.8162, 0.8269, 0.87199, 0.84, 0.962, 0.9966, 1.029, 1.028,
1.0161, 0.99187, 0.956, 0.86217]

experimental_y = [experimental_y_05, experimental_y_10, experimental_y_20, experimental_y_40]
experimental_u_bar = [experimental_u_bar_05, experimental_u_bar_10, experimental_u_bar_20, experimental_u_bar_40]
subplots = [141, 142, 143, 144]

Y = numpy.linspace(0.0, h, 50)
X = [0.5, 1.0, 2.0, 4.0]

pylab.figure()

for p in range(0, 4):
   probed_u = numpy.zeros(len(Y))
   u_bar = numpy.zeros(len(Y))
   pylab.subplot(subplots[p])
   x = X[p]*H
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
   pylab.plot(u_bar/U, Y/H, '-b', label="Fluidity (x/H = %f)" % X[p])
   pylab.plot(experimental_u_bar[p], experimental_y[p], 'ro', label="Experimental (x/H = %f)" % X[p])

#pylab.legend(loc=4)
#pylab.xlabel("Normalised mean velocity (u/U)")
#pylab.ylabel("Normalised height (y/H)")
#pylab.grid("on")

#pylab.show()
pylab.savefig('flow_past_a_cube.png')
