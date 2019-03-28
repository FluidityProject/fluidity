#!/usr/bin/env python

import numpy
import re

def check_velocity_iterations(filename):
  # Check the number of iterations required for the velocity solve to verify 
  # that solver behavior is as expected with fieldsplit preconditioner. 

  # Read fluidity log file:
  f=open(filename)
  log=f.read()
  f.close()

  # Grab all lines relating to the velocity solve iteration count:
  DeltaU_iteration_lines=re.findall('DeltaU PETSc n/o iterations.*', log, re.MULTILINE)

  # Split lines and store number of velocity solve iterations in an array:
  DeltaU_iterations=numpy.array([i.split(':')[-1] for i in DeltaU_iteration_lines])

  return int(DeltaU_iterations[-1])

