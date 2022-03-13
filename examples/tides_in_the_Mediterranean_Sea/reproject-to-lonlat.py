#!/usr/bin/env python3
import vtktools
import math
from scipy import sqrt
import vtk
import sys

ugrid = u=vtktools.vtu("med_10.pvtu")
# e.g. x="math.sin(x)"
x="math.atan2(y,x)*57.2957795"
y="math.asin(z/math.sqrt(x*x+y*y+z*z))*57.2957795"
z="math.sqrt(x*x+y*y+z*z)-6.37101e+06"

ugrid.ApplyProjection(x, y, z)

ugrid.Write(filename="tidesmedsea-flat.vtu")

u=vtktools.vtu("tidesmedsea-flat.vtu")

