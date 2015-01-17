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

x = numpy.linspace(0, 7000, 100)
p = numpy.zeros(len(x))

for t in range(1, 4):
   filename='mphase_volcanic_eruption_vw8_2d_' + str(t) + '.pvtu'
   vt=vtktools.vtu(filename)
   for i in range(0, len(x)):
      probedpressure = vtktools.vtu.ProbeData(vt, numpy.array([[x[i], 0, 0]]), 'Air::Pressure')
      p[i] = probedpressure[0]

   # Plot the numerical results
   pylab.plot(x, p, label="t = "+str(t)+"s")

pylab.legend()
pylab.xlabel("x")
pylab.ylabel("Pressure (Pa)")
pylab.grid("on")

pylab.show()
