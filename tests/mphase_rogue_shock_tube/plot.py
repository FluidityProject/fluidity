import pylab
import numpy
import vtktools
params = {'text.fontsize': 11,
         'legend.fontsize': 11,
         'xtick.labelsize': 11,
         'ytick.labelsize': 11,
         'lines.markersize': 6,
         'lines.linewidth': 1.5,
         'axes.titlesize': 'medium'}
pylab.rcParams.update(params)

fig=pylab.figure(1)
for fileindex in range(1, 740, 2):
   filename='mphase_rogue_shock_tube_' + str(fileindex) + '.pvtu'
   vt=vtktools.vtu(filename)

   # Time and coordinates
   t=vt.GetScalarField('Gas::Time')[0]
   xyz=vt.GetLocations()
   y=xyz[:,1]

   # Solution fields
   p=vt.GetScalarField('Gas::Pressure')
   probedpressure_upstream = vtktools.vtu.ProbeData(vt, numpy.array([[0.05, 1.095, 0]]), 'Gas::Pressure')
   if(fileindex == 1):
      pylab.plot( t, probedpressure_upstream,'b.', label="Numerical (Fluidity)")
   else:
      pylab.plot( t, probedpressure_upstream,'b.')

# Now plot some experimental data points from Rogue et al. (1998).
# These results are for the double-layered bed of 2mm glass beads.
experimental_t = [-0.001, -0.000566, -0.000326, -0.000312, -0.000300, -0.000290, -0.000279, -0.000254, -5.6979e-5, 9.8471e-5, 0.000155, 0.000212, 0.000268, 0.000338, 0.000322, 0.000334, 0.000387, 0.000399, 0.000425, 0.000495, 0.000666, 0.000794, 0.000893, 0.000994, 0.001, 0.00122, 0.00145, 0.00178, 0.002, 0.00217, 0.00256, 0.00297, 0.00338, 0.00342, 0.00397, 0.004]
experimental_p = [101032, 101089, 101507, 104221, 115457, 138316, 156525, 173962, 182123, 183305, 180998, 179446, 181390, 185273, 196506, 212779, 229444, 243005, 253469, 261226, 254662, 248867, 250042, 239983, 237276, 230327, 217572, 201730, 194408, 196745, 189822, 188325, 186828, 186834, 184193, 184197]

for i in range(0, len(experimental_t)):
   experimental_t[i] += 0.0005

pylab.plot(experimental_t, experimental_p,'g.', label="Experimental (Rogue et al., 1998)")

pylab.legend()
pylab.axis([0, 0.004, 0.0, 3.0e5])
pylab.xlabel("Time after shock wave hits particle bed (s)")
pylab.ylabel("Pressure at upstream gauge (Pa)")

#pylab.show()
pylab.savefig('mphase_rogue_shock_tube.png')
