from fluidity_tools import stat_parser
from numpy import zeros
import math

def get_detector_variance(filename, timesteps, agents):
  s = stat_parser(filename)
  variance = zeros(timesteps)

  for t in range(0,timesteps):
    var = 0.0
    for i in range(1,agents):
      var = var + math.pow(s[str(i)]['position'][0][t], 2)
    variance[t] = var / agents
  return variance

