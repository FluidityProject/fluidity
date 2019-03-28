#!/usr/bin/env python

from scipy.version import version as SciPyVersion

if tuple(int(x) for x in SciPyVersion.split('.')) < (0, 9, 0):
    from Scientific.IO.NetCDF import NetCDFFile as netcdf_file
else:
    from scipy.io.netcdf import netcdf_file

from numpy import arange, zeros

def create(missingdata = False, missingdimension = False, missingvariable = False, incorrectdimension = False, incorrectvariable = False):

  if (missingdata):
    filename = 'missingdata.nc'
    description = ', with missing data.'
  elif (missingdimension):
    filename = 'missingdimension.nc'
    description = ', with a missing dimension.'
  elif (missingvariable):
    filename = 'missingvariable.nc'
    description = ', with a missing variable.'
  elif (incorrectdimension):
    filename = 'incorrectdimension.nc'
    description = ', with an incorrect dimension label.'
  elif (incorrectvariable):
    filename = 'incorrectvariable.nc'
    description = ', with an incorrect variable label.'
  else:
    filename = 'valid.nc'
    description = '.'

  print("Creating " + filename + description)

  f = netcdf_file(filename, 'w')
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
      h[j,i] = x[i] * y[j]

  if (not incorrectdimension):
    xdimlabel = 'x'
  else:
    xdimlabel = 'lat'
  ydimlabel = 'y'
  if (not incorrectvariable):
    zdimlabel = 'z'
  else:
    zdimlabel = 'height'

  print(xdimlabel, ydimlabel, zdimlabel)

  # dimensions
  f.createDimension(xdimlabel, len(x))
  if (not missingdimension):
    f.createDimension(ydimlabel, len(y))

  # variables
  fx = f.createVariable(xdimlabel, 'd', (xdimlabel,))
  fx[:] = x
  if (not missingdimension):
    fy = f.createVariable(ydimlabel, 'd', (ydimlabel,))
    fy[:] = y
    if (not missingvariable):
      fh = f.createVariable(zdimlabel, 'd', (xdimlabel, ydimlabel,))
      fh[:] = h
  else:
    fh = f.createVariable(zdimlabel, 'd', (xdimlabel,))
    fh[:] = x

  f.close()

create()
create(missingdata=True)
create(missingdimension=True)
create(missingvariable=True)
create(incorrectdimension=True)
create(incorrectvariable=True)



