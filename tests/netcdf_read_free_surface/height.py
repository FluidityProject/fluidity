#!/usr/bin/env python3

from matplotlib.mlab import bivariate_normal

def function(position):
  xg = 3.0 * position[0]
  yg = 3.0 * position[1]
  h1 = bivariate_normal(xg, yg, 1.0, 1.0, 0.0, 0.0)
  h2 = bivariate_normal(xg, yg, 1.5, 0.5, 1, 1)
  h = 10.0 * (h1 - h2)
  return h


