#!/usr/bin/env python

import re

def parse(code):
  """Given the string containing code for a functional or its derivative, work out
  what dependencies on which state it has. E.g., if the code was

  from math import sin, pi
  coord  = states[n]["Fluid"].vector_fields["Coordinate"]
  u = states[n+1]["Fluid"].scalar_fields["Velocity"]
  du = states[n-1]["Fluid"].scalar_fields["VelocityDerivative"]

  for i in range(du.node_count):
    x = coord.node_val[i][0]
      du.set(i, 0.01125*pi**2*sin(3.0/20*(x + 10)*pi)

  then this routine should return the list
  [-1, 0, +1]."""

  # My beautiful regex, made with the help of http://re.dabase.com/
  regex = re.compile('''states\[(?P<n>[n0-9+-]*)\]''')
  n = 0

  return sorted(set(map(eval, re.findall(regex, code))))

if __name__ == "__main__":
  code = '''
  from math import sin, pi
  coord  = states[n]["Fluid"].vector_fields["Coordinate"]
  u = states[n+1]["Fluid"].scalar_fields["Velocity"]
  du = states[n-1]["Fluid"].scalar_fields["VelocityDerivative"]

  for i in range(du.node_count):
    x = coord.node_val[i][0]
      du.set(i, 0.01125*pi**2*sin(3.0/20*(x + 10)*pi)
  '''

  print parse(code)
