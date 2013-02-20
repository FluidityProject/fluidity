#!/usr/bin/env python

import numpy
import re

def check_velocity_iterations():
  # Check the number of iterations required for the velocity solve to verify 
  # that solver behavior is as expected with fieldsplit preconditioner. 

  # Read fluidity log file:
  f=file('fluidity.log-0')
  log=f.read()
  f.close()

  # Grab all lines relating to the velocity solve iteration count:
  DeltaU_iteration_lines=re.findall('DeltaU PETSc n/o iterations.*', log, re.MULTILINE)

  # Split lines and store number of velocity solve iterations in an array:
  DeltaU_iterations=numpy.array([i.split(':')[-1] for i in DeltaU_iteration_lines])

  # Check velocity solve is behaving as expected (i.e. that number of iterations of final solve == 7 or 8):
  velocity_solver_as_expected = (int(DeltaU_iterations[-1]) == 7 or int(DeltaU_iterations[-1]) == 8)

  return velocity_solver_as_expected

