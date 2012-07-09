import pylab
import numpy
import vtktools
params = {'text.fontsize': 11,
         'legend.fontsize': 9,
         'xtick.labelsize': 11,
         'ytick.labelsize': 11,
         'lines.markersize': 6,
         'lines.linewidth': 1.5,
         'axes.titlesize': 'medium'}
pylab.rcParams.update(params)

numerical_t = []
numerical_p_upstream = []
numerical_p_downstream = []

for fileindex in range(1, 1334, 5):
   filename='mphase_rogue_shock_tube_dense_bed_glass_' + str(fileindex) + '.pvtu'
   vt=vtktools.vtu(filename)

   # Time and coordinates
   t=vt.GetScalarField('Gas::Time')[0]
   xyz=vt.GetLocations()
   y=xyz[:,1]

   # The Pressure is set to zero after an adapt, so ignore this dump
   if(fileindex % 20 == 0):
      continue

   # Solution fields
   p=vt.GetScalarField('Gas::Pressure')
   probedpressure_upstream = vtktools.vtu.ProbeData(vt, numpy.array([[0.0125, 1.35 - 0.11, 0]]), 'Gas::Pressure')
   probedpressure_downstream = vtktools.vtu.ProbeData(vt, numpy.array([[0.0125, 1.37 + 0.718, 0]]), 'Gas::Pressure')

   if(probedpressure_upstream == 0.0 or probedpressure_downstream == 0.0):
      continue

   numerical_t.append(t - 0.0008)
   numerical_p_upstream.append(probedpressure_upstream[0])
   numerical_p_downstream.append(probedpressure_downstream[0])

# Plot the numerical results
pylab.plot(numerical_t, numerical_p_upstream,'g-', label="Numerical - Upstream (Fluidity)")
pylab.plot(numerical_t, numerical_p_downstream,'r-', label="Numerical - Downstream (Fluidity)")

# Now plot some experimental data points from Rogue et al. (1998).
# These results are for the double-layered bed of 2mm glass beads.
experimental_t = [-0.001, -0.000566, -0.000326, -0.000312, -0.000300, -0.000290, -0.000279, -0.000254, -5.6979e-5, 9.8471e-5, 0.000155, 0.000212, 0.000268, 0.000338, 0.000322, 0.000334, 0.000387, 0.000399, 0.000425, 0.000495, 0.000666, 0.000794, 0.000893, 0.000994, 0.001, 0.00122, 0.00145, 0.00178, 0.002, 0.00217, 0.00256, 0.00297, 0.00338, 0.00342, 0.00397, 0.004]
experimental_p = [101032, 101089, 101507, 104221, 115457, 138316, 156525, 173962, 182123, 183305, 180998, 179446, 181390, 185273, 196506, 212779, 229444, 243005, 253469, 261226, 254662, 248867, 250042, 239983, 237276, 230327, 217572, 201730, 194408, 196745, 189822, 188325, 186828, 186834, 184193, 184197]

# These results are for the 2cm thick particle bed of 1.5mm glass beads.
experimental_t_upstream = [-0.001, -0.00065, -0.0004, -0.00033, -0.000306, -0.000293, -0.000282, -0.00027, -0.000258, -0.000231, -0.000147, -7.854e-5, 1.87e-5, 8.75e-5, 0.000226, 0.000323, 0.000335, 0.000347, 0.000360, 0.000370, 0.000395, 0.000421, 0.000601, 0.000864, 0.001, 0.00129, 0.00148, 0.00176, 0.002, 0.00222, 0.00254, 0.00287, 0.00317, 0.00348, 0.00377, 0.00398, 0.00409, 0.00424, 0.00438, 0.00453, 0.00474, 0.0051, 0.00548, 0.0058, 0.00616, 0.00667, 0.0071, 0.00749, 0.0077, 0.00792]
experimental_p_upstream = [100543, 100589, 101835, 114864, 124791, 144022, 164495, 175661, 179384, 184349, 181874, 184361, 181266, 184374, 182523, 182531, 197420, 213550, 227199, 253874, 285514, 297923, 295456, 301059, 291768, 298611, 293043, 295546, 285642, 289998, 280097, 281364, 277665, 272726, 272749, 268422, 260987, 264721, 262870, 261642, 249871, 239974, 226975, 214593, 198492, 184884, 172513, 162615, 161391, 155825]
experimental_t_downstream = [-0.001, -0.0006, -0.00042, -0.000333, -5.56e-5, 0.000347, 0.000819, 0.00124, 0.0015, 0.0016, 0.0018, 0.00192, 0.00199, 0.00204, 0.00212, 0.00231, 0.00251, 0.00267, 0.00301, 0.00326, 0.00359, 0.00388, 0.00429, 0.00459, 0.00474, 0.00498, 0.00522, 0.00548, 0.00568, 0.00596, 0.00623, 0.00659, 0.00690, 0.00722, 0.00753, 0.00778]
experimental_p_downstream = [101142, 100558, 100580, 101212, 100617, 100659, 100708, 101375, 100779, 102039, 100187, 100199, 106440, 113304, 118300, 120191, 120211, 120850, 122133, 125276, 126558, 128458, 131617, 134143, 137899, 136677, 141066, 142340, 145477, 147376, 150522, 149937, 151215, 150002, 150035, 148189]

pylab.plot(experimental_t_upstream, experimental_p_upstream,'go', label="Experimental - Upstream (Rogue et al., 1998)")
pylab.plot(experimental_t_downstream, experimental_p_downstream,'ro', label="Experimental - Downstream (Rogue et al., 1998)")

# Plot the t = 0 vertical line
pylab.plot(numpy.linspace(0,0,100), numpy.linspace(0, 3.5e5, 100), 'k--')

pylab.legend()
pylab.axis([-0.002, 0.008, 0.9e5, 3.5e5])
pylab.xlabel("Time after shock wave hits particle bed (s)")
pylab.ylabel("Pressure at gauge (Pa)")
pylab.grid("on")

#pylab.show()
pylab.savefig('mphase_rogue_shock_tube_dense_bed_glass.png')
