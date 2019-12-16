#!/usr/bin/env python3

import vtktools
import math

tolerance=1.0E-12
u=vtktools.vtu("shelf_0.vtu")
location=u.GetLocations()
(ilen, jlen) = location.shape

for i in range(ilen):
  if location[i,0] + tolerance < 0.0:
    raise Exception("Failure: Node outside of x in [0.0,1.0E6].  Lower bound exceeded.") 
  elif 1.0E6 < location[i,0] - tolerance:
    raise Exception("Failure: Node outside of x in [0.0,1.0E6].  Upper bound exceeded.") 

for i in range(ilen):
  if location[i,1] + tolerance < -1.0E3:
    raise Exception("Failure: Node outside of y in [-1.0E3,1.0E2].  Lower bound exceeded.") 
  elif 1.0E2 < location[i,1] - tolerance:
    raise Exception("Failure: Node outside of y in [-1.0E3,1.0E2].  Upper bound exceeded.") 

count=0
for i in range(ilen):
  if location[i,1] > -1.0E3:
    count += 1
    if 0.0 <= location[i,0] < 5.0E5:
      if ( location[i,1] - ( - 9.0E2 + (9.0E2/5.0E5) * location[i,0] ) ) > tolerance:
        raise Exception("Failure: Surface node out of bounds in the region where the mesh is transformed.  Check scripts/gmsh_mesh_transform.") 
    elif 5.0E5 <= location[i,0] <= 1.0E6:
      if ( location[i,1] - 1.0E2 ) > tolerance:
        raise Exception("Failure: Surface node out of bounds in the region where the free-surface initial condition is applied.") 
    else:
      raise Exception("Failure: Node outside of x in [0.0,1.0E6].") 

if count != 102:
  raise Exception("Failure: There are not 102 nodes on the surface.") 
  

