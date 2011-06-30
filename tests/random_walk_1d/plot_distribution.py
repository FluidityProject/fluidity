from fluidity_tools import stat_parser
import vtktools
from numpy import zeros, arange
from pylab import figure, show, colorbar

def plot_detector_distribution(filename):
  s = stat_parser(filename)
  
  det_count = zeros((41,100))
  for i in range(1,4000):
    for t in range(100):
      x = s[str(i)]['position'][0][t]
      det_count[x,t] = det_count[x,t]+1

  #print det_count
  fig = figure(figsize=(10,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  cs=ax.contourf(det_count, arange(70,130,5))
  pp=colorbar(cs)

  show()

  return

def plot_diffusivity(file):
  u=vtktools.vtu(file)
  z = u.GetLocations()[:,0]
  #z = [ -d for d in z]
  K = u.GetScalarField("Diffusivity")

  fig = figure(figsize=(5,6),dpi=90)
  ax = fig.add_axes([.1,.1,.8,.8])
  ax.plot(K, z)

  show()

### Main ###

#plot_diffusivity("random_walk_1d_0.vtu")

plot_detector_distribution("Naive_RW.detectors")
