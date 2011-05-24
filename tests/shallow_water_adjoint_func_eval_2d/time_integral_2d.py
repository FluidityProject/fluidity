from fluidity_tools import stat_parser

def time_integral(filename, field):
  stat = stat_parser(filename)
  values = stat["Fluid"][field]["max"]

  J = sum(values)
  return J
