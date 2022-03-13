#!/usr/bin/env python3

from scipy.version import version as SciPyVersion

if tuple(int(x) for x in SciPyVersion.split('.')) < (0, 9, 0):
    from Scientific.IO.NetCDF import NetCDFFile as netcdf_file
else:
    from scipy.io.netcdf import netcdf_file

from numpy import arange, zeros
import height

f = netcdf_file('height.nc', 'w')
f.description = 'Example free surface height.'

x = arange(-1.2, 1.21, 0.2)
y = arange(-1.2, 1.21, 0.2)
h = zeros((len(x),len(y)))

for i in range(len(x)):
  for j in range(len(y)):
    # Nice ordering netCDF API - y,x !
    h[j,i] = height.function([x[i], y[j]])

# dimensions
f.createDimension('x', len(x))
f.createDimension('y', len(y))

# variables
fx = f.createVariable('x', 'd', ('x',))
fy = f.createVariable('y', 'd', ('y',))
fh = f.createVariable('z', 'd', ('x', 'y',))

fx[:] = x
fy[:] = y
fh[:] = h

f.close()

