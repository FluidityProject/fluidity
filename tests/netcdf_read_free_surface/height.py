#!/usr/bin/env python3

from scipy.stats import multivariate_normal

def function(position):
  xg = 3.0 * position[0]
  yg = 3.0 * position[1]
  h1 = multivariate_normal.pdf([xg, yg], cov=[[1, 0], [0, 1]])
  h2 = multivariate_normal.pdf([xg, yg], mean=[1, 1], cov=[[1.5**2, 0], [0, 0.5**2]])
  h = 10.0 * (h1 - h2)
  return h


