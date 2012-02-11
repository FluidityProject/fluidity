#!/usr/bin/env python

from Scientific.IO.NetCDF import NetCDFFile
from numpy import arange, zeros
import height

def create(missingdata = False, missingdimension = False, missingvariable = False):

  if (missingdata):
    filename = 'missingdata.nc'
    description = ', with missing data.'
  elif (missingdimension):
    filename = 'missingdimension.nc'
    description = ', with a missing dimension.'
  elif (missingvariable):
    filename = 'missingvariable.nc'
    description = ', with a missing variable.'
  else:
    filename = 'valid.nc'
    description = '.'

  f = NetCDFFile(filename, 'w')
  f.description = 'Example free surface height' + description
  
  if (missingdata):
    offset = -0.4
  else:
    offset = 0.0

  x = arange(-1.2 + offset, 1.21 + offset, 0.2)
  y = arange(-1.2, 1.21, 0.2)
  h = zeros((len(x),len(y)))

  for i in range(len(x)):
    for j in range(len(y)):
      # Nice ordering netCDF API - y,x !
      h[j,i] = height.function([x[i], y[j]])

  # dimensions
  f.createDimension('x', len(x))
  if (not missingdimension):
    f.createDimension('y', len(y))

  # variables
  fx = f.createVariable('x', 'd', ('x',))
  fx[:] = x
  if (not missingdimension):
    fy = f.createVariable('y', 'd', ('y',))
    fy[:] = y
    if (not missingvariable):
      fh = f.createVariable('z', 'd', ('x', 'y',))
      fh[:] = h
  else:
    fh = f.createVariable('z', 'd', ('x',))
    fh[:] = x

  f.close()

create()
create(missingdata=True)
create(missingdimension=True)
create(missingvariable=True)



