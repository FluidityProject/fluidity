from fluidity_tools import stat_parser
import vtktools
from numpy import zeros, arange
from pylab import figure, show, colorbar

def plot_detector_distribution(filename):
  s = stat_parser(filename)
  
  det_count = zeros((41,600))
  for i in range(1,1000):
    for t in range(0,600):
      x = s[str(i)]['position'][0][t]
      det_count[x,t] = det_count[x,t]+1

  fig = figure(figsize=(10,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  cs=ax.contourf(det_count, arange(0,100,5))
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

plot_detector_distribution("Naive_RW.detectors")

plot_detector_distribution("Diffusive_RW.detectors")

show()
