def f_to_c(f):
  um = 4 * 0.45 / 9
  H = 0.41
  D = 0.1
  p = 1.0
  c = (2 * f) / (p * um**2 * D * H)
  return c
  
