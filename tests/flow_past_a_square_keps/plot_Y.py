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
N = 77

# Dimensions
H = 0.04
U = 0.535
h = 7*H

experimental_Y0 = [-1.13, -1.075, -1.009, -0.943, -0.877, -0.82, -0.783,
-0.764, -0.726, -0.688, -0.66, -0.632, -0.60377,
-0.5849, -0.5566, 0.5566, 0.5849, 0.60377, 0.632, 0.66, 0.688, 0.726, 0.764, 0.783, 0.82, 0.877, 0.943, 1.009, 1.075, 1.13]
experimental_u_bar0 = [1.33, 1.359, 1.394, 1.423, 1.446, 1.388, 1.214,
0.987, 0.633, 0.249, -0.046, -0.168, -0.24979,
-0.2614, -0.2904, -0.2904, -0.2614, -0.24979, -0.168, -0.046, 0.249, 0.633, 0.987, 1.214, 1.388, 1.446, 1.423, 1.394, 1.359, 1.33]

experimental_Y1 = [-2, -1.75, -1.5, -1.25, -1.12, -1, -0.876, -0.75, -0.63,
-0.49, -0.37, -0.2469, -0.12, 0, 0.12, 0.2469, 0.37, 0.49, 0.63, 0.75, 0.876, 1.0, 1.12, 1.25, 1.5, 1.75, 2.0]
experimental_u_bar1 = [1.1758, 1.19, 1.225, 1.24, 1.235, 1.2, 1.016, 0.774, 0.3736,
0.17, 0.016, -0.093, -0.18, -0.18, -0.18, -0.093, 0.016, 0.17, 0.3736, 0.774, 1.016, 1.2, 1.235, 1.24, 1.225, 1.19, 1.1758]

experimental_Y3 = [-2, -1.75, -1.62, -1.5, -1.37, -1.248, -1.12, -0.9965, -0.87,
-0.744, -0.6186, -0.49, -0.366, -0.24, -0.11, 0.0, 0.11, 0.24, 0.366, 0.49, 0.6186, 0.744, 0.87, 0.9965, 1.12, 1.248, 1.37, 1.5, 1.62, 1.75, 2.0]
experimental_u_bar3 = [1.055, 1.0297, 1.0, 0.97, 0.948, 0.91, 0.885, 0.834, 0.8,
0.753, 0.71489, 0.66, 0.6127, 0.57, 0.5489, 0.5489, 0.5489, 0.57, 0.6127, 0.66, 0.71489, 0.753, 0.8, 0.834, 0.885, 0.91, 0.948, 0.97, 1.0, 1.0297, 1.055]

experimental_Y8 = [-4, -3.05, -2.547, -2.295, -2.04, -1.788, -1.53, -1.28, -1.0,
-0.77, -0.516, -0.26, 0.0, 0.26, 0.516, 0.77, 1.0, 1.28, 1.53, 1.788, 2.04, 2.295, 2.547, 3.05, 4.0]
experimental_u_bar8 = [1.1, 1.078, 1.039, 1.02, 1.0, 0.972, 0.924, 0.869, 0.806,
0.742, 0.7, 0.65, 0.618, 0.65, 0.7, 0.742, 0.806, 0.869, 0.924, 0.972, 1.0, 1.02, 1.039, 1.078, 1.1]

xD = [0.0, 1.0, 3.0, 8.0]
experimental_Y = [experimental_Y0, experimental_Y1, experimental_Y3, experimental_Y8]
experimental_u_bar = [experimental_u_bar0, experimental_u_bar1, experimental_u_bar3, experimental_u_bar8]

for p in range(0, 4):
   x = xD[p]*H
   if(p == 3):
      Y = numpy.linspace(-4*H, 4*H, 300)
   else:
      Y = numpy.linspace(-2*H, 2*H, 300)
   
   pylab.figure(p)

   u_bar = numpy.zeros(len(Y))

   filename='les/flow_past_a_square_' + str(N) + '.pvtu'
   vt=vtktools.vtu(filename)

   # Time and coordinates
   t = vt.GetScalarField('Time')[0]
   xyz = vt.GetLocations()

   for i in range(0,len(Y)): 
      data = vtktools.vtu.ProbeData(vt, numpy.array([[x, Y[i], 0]]), 'TimeAveragedVelocity')
      u_bar[i] = data[0][0]

   # Plot the numerical results
   pylab.plot(Y/H, u_bar/U, '-b', label="LES (x/H = %f)" % xD[p])

   filename='noles/flow_past_a_square_' + str(N) + '.pvtu'
   vt=vtktools.vtu(filename)

   # Time and coordinates
   t = vt.GetScalarField('Time')[0]
   xyz = vt.GetLocations()

   for i in range(0,len(Y)): 
      data = vtktools.vtu.ProbeData(vt, numpy.array([[x, Y[i], 0]]), 'TimeAveragedVelocity')
      u_bar[i] = data[0][0]

   # Plot the numerical results
   pylab.plot(Y/H, u_bar/U, '-g', label="No LES (x/H = %f)" % xD[p])

   # Plot experimental results
   pylab.plot(experimental_Y[p], experimental_u_bar[p] , 'xr', label="Experimental (Lyn and Rodi (1994))")

   pylab.legend(loc=4)
   pylab.xlabel("Normalised height (y/H)")
   pylab.ylabel("Normalised mean velocity (u/U)")
   pylab.grid("on")

   #pylab.show()
   pylab.savefig('flow_past_a_square_' + str(p) + '.png')
