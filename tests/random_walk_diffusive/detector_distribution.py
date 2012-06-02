from fluidity_tools import stat_parser
from numpy import zeros
import math

def get_distribution(filename, timesteps, layers, agents):
  s = stat_parser(filename)
  
  det_count = zeros((layers,timesteps))
  for i in range(1,agents):
    for t in range(0,timesteps):
      x = math.floor(abs(s['position'][i][2][t]))
      det_count[x,t] = det_count[x,t]+1
  return det_count
