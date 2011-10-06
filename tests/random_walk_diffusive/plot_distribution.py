from fluidity_tools import stat_parser
import vtktools
from matplotlib import ticker
from numpy import zeros, arange, meshgrid
from pylab import figure, show, colorbar, xlabel, ylabel, savefig
import math

def plot_detector_distribution(filename, timesteps, agents, layers):
  s = stat_parser(filename)
  
  det_count = zeros((layers+1,timesteps))
  for i in range(1,agents):
    for t in range(0,timesteps):
      x = math.floor(abs(s[str(i)]['position'][0][t]))
      det_count[x,t] = det_count[x,t]+1

  fig = figure(figsize=(10,4),dpi=90)
  ax = fig.add_axes([.09,.16,.98,.7])
  x = arange(0., timesteps)
  y = arange(0., layers)
  X, Y = meshgrid(x, y)
  ax.invert_yaxis()
  cs=ax.contourf(det_count, arange(-1.e-12,40,1))
  xlabel('Time (s)')
  ylabel('Depth (m)')
  pp=colorbar(cs)

  return

def plot_diffusivity(file):
  u=vtktools.vtu(file)
  z = u.GetLocations()[:,0]
  K = u.GetScalarField("Diffusivity")
  K_grad = u.GetVectorField("Diffusivity_grad")[:,0]

  fig = figure(figsize=(3,4),dpi=90)
  ax = fig.add_axes([.28,.16,.6,.7])
  ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=2))
  ax.set_ylim(-100., 0.)
  xlabel('Diffusivity (m$^2$s$^{-1}$)')
  ylabel('Depth (m)')
  ax.plot(K, z)

  fig = figure(figsize=(3,4),dpi=90)
  ax = fig.add_axes([.28,.16,.6,.7])
  ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=2))
  ax.set_ylim(-100., 0.)
  xlabel('Diffusivity gradient')
  ylabel('Depth (m)')
  ax.plot(K_grad, z)
  return


### Main ###

plot_diffusivity("random_walk_1d_0.vtu")

plot_detector_distribution("Diffusive_RW.detectors", 600, 1000, 60)

show()
