from fluidity_tools import stat_parser
import vtktools
from numpy import zeros, arange
from pylab import figure, show, colorbar

def plot_detector_distribution(filename, timesteps, agents, layers):
  s = stat_parser(filename)
  
  det_count = zeros((layers+1,timesteps))
  for i in range(1,agents):
    for t in range(0,timesteps):
      x = round(s[str(i)]['position'][0][t])
      det_count[x,t] = det_count[x,t]+1

  fig = figure(figsize=(10,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  cs=ax.contourf(det_count, arange(-1.e-12,60,5))
  pp=colorbar(cs)

  return

def plot_diffusivity(file):
  u=vtktools.vtu(file)
  z = u.GetLocations()[:,0]
  K = u.GetScalarField("Diffusivity")
  K_grad = u.GetVectorField("Diffusivity_grad")[:,0]

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  ax.plot(K, z)

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  ax.plot(K_grad, z)
  return


### Main ###

plot_diffusivity("random_walk_1d_0.vtu")

plot_detector_distribution("Naive_RW.detectors", 600, 1000, 40)

plot_detector_distribution("Diffusive_RW.detectors", 600, 1000, 40)

plot_detector_distribution("Diffusive_RW_intern.detectors", 600, 1000, 40)

show()
