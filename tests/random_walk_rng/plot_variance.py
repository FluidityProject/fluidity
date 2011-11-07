from detector_variance import get_detector_variance
import math
from numpy import arange
from pylab import figure, show, xlabel, ylabel, savefig

def plot_variance(variance, agents, filename):
  min_factor = 1. - math.sqrt(2./agents)
  max_factor = 1. + math.sqrt(2./agents)
  min_var = [ t / 12. * min_factor for t in arange(0,len(variance))]
  max_var = [ t / 12. * max_factor for t in arange(0,len(variance))]

  fig = figure(figsize=(6,4),dpi=90)
  ax = fig.add_axes([.12,.12,.8,.8])
  ax.plot(arange(0,len(variance)), variance, 'k-')
  ax.plot(arange(0,len(variance)), min_var, 'k--')
  ax.plot(arange(0,len(variance)), max_var, 'k--')

  xlabel('Timestep')
  ylabel('Variance')
  savefig(filename + '.png', dpi=90,format='png')


variance = get_detector_variance("Fortran_Naive.detectors", 2000, 1000)
plot_variance(variance, 2000, "Fortran_RNG")

variance = get_detector_variance("Python_Naive.detectors", 2000, 1000)
plot_variance(variance, 2000, "Python_RNG")

show()
