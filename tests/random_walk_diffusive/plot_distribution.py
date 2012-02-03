import vtktools
from matplotlib import ticker
from numpy import arange, meshgrid, argsort
from pylab import figure, show, colorbar, xlabel, ylabel, savefig
from detector_distribution import get_distribution

def plot_detector_distribution(det_count, filename, timesteps, layers, delta_t):
  fig = figure(figsize=(10,6),dpi=90)
  ax = fig.add_axes([.06,.1,.98,.86])
  x = arange(0., timesteps)
  x = x*delta_t
  y = arange(0., layers)
  X, Y = meshgrid(x, y)
  cs=ax.contourf(X, Y, det_count, arange(-1e-12,40,5))
  ax.set_xlim(0., timesteps*delta_t)
  ax.set_ylim(layers-1., 0.)
  xlabel('Time (s)')
  ylabel('Depth (m)')
  pp=colorbar(cs)

  savefig('./' + filename + '.png', dpi=90,format='png')

def plot_diffusivity(file):
  # open vtu and derive the field indices of the edge at (x=0,y=0) ordered by depth
  u=vtktools.vtu(file)
  pos = u.GetLocations()
  ind = get_1d_indices(pos)

  # from this we can derive the 1D profile of any field like this:
  z = [-pos[i,2] for i in ind]

  diffusivity = u.GetScalarField("Diffusivity")
  K = [diffusivity[i] for i in ind]

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.24,.1,.6,.86])
  ax.plot(K, z)
  ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
  ax.set_ylim(60., 0.)
  xlabel('Diffusivity ($m^2s^{-1}$)')
  ylabel('Depth (m)')
  savefig('./K.png', dpi=90,format='png')

  diffusivity_grad = u.GetVectorField("DiffusivityGradient")[:,2]
  # not perfectly sure why we need to do * -1 here, 
  # but the gradient looks wrong otherwise
  K_grad = [-1. * diffusivity_grad[i] for i in ind]    

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.24,.1,.6,.86])
  ax.plot(K_grad, z)
  ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
  ax.set_ylim(60., 0.)
  xlabel('Diffusivity gradient ($ms^{-1}$)')
  ylabel('Depth (m)')
  savefig('./Kgrad.png', dpi=90,format='png')

  diffusivity_second_grad = u.GetVectorField("DiffusivitySecondGradient")[:,2]
  K_second_grad = [diffusivity_second_grad[i] for i in ind]    

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.24,.1,.6,.86])
  ax.plot(K_second_grad, z)
  ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
  ax.set_ylim(60., 0.)
  xlabel('Diffusivity 2nd gradient ($$)')
  ylabel('Depth (m)')
  savefig('./Kgrad2.png', dpi=90,format='png')

  return

def get_1d_indices(pos, x0=0, y0=0, tolerance=1.0e-5):
  """ Return the field indices corresponding to the ordered depth values at position (x0, y0)
  """
  ind = argsort(-pos[:,2])
  indices = []
  for i in ind:
    if (x0-tolerance < pos[i][0] < x0+tolerance and y0-tolerance < pos[i][1] < y0+tolerance):
      indices.append(i)
  return indices


### Main ###

plot_diffusivity("random_walk_diffusive_0.vtu")

det_count_diffusive = get_distribution("RandomWalkDiffusive.detectors", 667, 60, 1200)
plot_detector_distribution(det_count_diffusive, "AgentRW", 667, 60, 1800)

show()
