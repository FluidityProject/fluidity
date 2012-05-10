from fluidity_tools import stat_parser
import vtktools
from numpy import zeros, arange, meshgrid
from pylab import figure, show, colorbar, xlabel, ylabel, savefig
from detector_distribution import get_distribution

def plot_detector_distribution(det_count, filename, timesteps, layers, delta_t):
  fig = figure(figsize=(10,6),dpi=90)
  ax = fig.add_axes([.06,.1,.98,.86])
  x = arange(0., timesteps)
  x = x*delta_t
  # this should do the trick as well, but for some reason inverts the profile
  #y = arange(0., layers)
  y = arange(-layers, 0.)
  y = (y * -1.) -1.
  X, Y = meshgrid(x, y)
  cs=ax.contourf(X, Y, det_count, arange(-1.e-12,80,5))
  ax.set_xlim(0., timesteps*delta_t)
  ax.set_ylim(layers-1., 0.)
  xlabel('Time (s)')
  ylabel('Depth (m)')
  pp=colorbar(cs)

  savefig('./' + filename + '.png', dpi=90,format='png')

def plot_diffusivity(file):
  u=vtktools.vtu(file)
  z = u.GetLocations()[:,2]
  K = u.GetScalarField("Diffusivity")
  K_grad = u.GetVectorField("DiffusivityGradient")[:,2]

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  ax.plot(K, z)
  ax.set_xlim(0.0,0.03)

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  ax.plot(K_grad, z)
  ax.set_xlim(-0.02,0.01)
  return


### Main ###

plot_diffusivity("random_walk_3d_0.vtu")

det_count_naive = get_distribution("NaiveRW.detectors", 600, 40, 1000)
plot_detector_distribution(det_count_naive, "NaiveRW", 600, 40, 6.)

det_count_diffusive = get_distribution("DiffusiveRW.detectors", 600, 40, 1000)
plot_detector_distribution(det_count_diffusive, "DiffusiveRW", 600, 40, 6.)

show()
