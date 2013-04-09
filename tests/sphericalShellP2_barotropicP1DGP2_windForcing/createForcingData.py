from Scientific.IO import NetCDF
import numpy
import os
from datetime import date
import sys

sys.stdout.write("Creating uncompressed file...")
#Calculate hours from 1900-01-01 00:00:0.0 to 1987-01-01 00:00:0.0
start_hour = (date(1987,01,01) - date(1900,01,01)).days*24
#Calculate hours from 1900-01-01 00:00:0.0 to 1988-01-01 00:00:0.0
end_hour = (date(1988,01,01) - date(1900,01,01)).days*24

file = NetCDF.NetCDFFile("test_forcing_I.nc", 'w')

file.createDimension('longitude', 144)
longitude = file.createVariable("longitude",'f',('longitude',))
setattr(longitude,'units','degrees_east')
setattr(longitude,'long_name','longitude')
longitude[:] = list(numpy.linspace(0.0,357.5,144))

file.createDimension('latitude', 73)
latitude = file.createVariable("latitude",'f',('latitude',))
setattr(latitude,'units','degrees_north')
setattr(latitude,'long_name','latitude')
latitude[:] = list(numpy.linspace(90.0,-90.0,73))

file.createDimension('time', 0)
time = file.createVariable("time",'i',('time',))
setattr(time,'units','hours since 1900-01-01 00:00:0.0')
setattr(time,'long_name','time')
time[:] = list(numpy.linspace(start_hour,end_hour,1461))

lon = numpy.array(longitude[:])*numpy.pi/180.0
lat = numpy.array(latitude[:])*numpy.pi/180.0
t = numpy.array(time[:])

u10 = file.createVariable("u10",'f',('time', 'latitude', 'longitude'))
setattr(u10,'units','m s**-1')
setattr(u10,'long_name','10 metre U wind component')
setattr(u10,'FillValue',-32767)
setattr(u10,'missing_value',-32767)
u10_array = numpy.zeros((1461,73,144))
for lat_index in range(0,len(lat)):
 u10_array[:,lat_index,:] = -1.0*numpy.cos(lat[lat_index])
u10[:,:,:] = list(u10_array)

v10 = file.createVariable("v10",'f',('time', 'latitude', 'longitude'))
setattr(v10,'units','m s**-1')
setattr(v10,'long_name','10 metre V wind component')
setattr(v10,'FillValue',-32767)
setattr(v10,'missing_value',-32767)
v10[:,:,:] = list(numpy.zeros((1461,73,144)))

d2m = file.createVariable("d2m",'f',('time', 'latitude', 'longitude'))
setattr(d2m,'units','K')
setattr(d2m,'long_name','2 metre dewpoint temperature')
setattr(d2m,'FillValue',-32767)
setattr(d2m,'missing_value',-32767)
d2m[:,:,:] = list(numpy.zeros((1461,73,144)))

t2m = file.createVariable("t2m",'f',('time', 'latitude', 'longitude'))
setattr(t2m,'units','K')
setattr(t2m,'long_name','2 metre temperature')
setattr(t2m,'FillValue',-32767)
setattr(t2m,'missing_value',-32767)
t2m[:,:,:] = list(numpy.zeros((1461,73,144)))

msl = file.createVariable("msl",'f',('time', 'latitude', 'longitude'))
setattr(msl,'units','Pa')
setattr(msl,'long_name','Mean sea level pressure')
setattr(msl,'FillValue',-32767)
setattr(msl,'missing_value',-32767)
msl[:,:,:] = list(numpy.zeros((1461,73,144)))

ro = file.createVariable("ro",'f',('time', 'latitude', 'longitude'))
setattr(ro,'units','m')
setattr(ro,'long_name','Runoff')
setattr(ro,'FillValue',-32767)
setattr(ro,'missing_value',-32767)
ro[:,:,:] = list(numpy.zeros((1461,73,144)))

ssrd = file.createVariable("ssrd",'f',('time', 'latitude', 'longitude'))
setattr(ssrd,'units','W m**-2 s')
setattr(ssrd,'long_name','Surface solar radiation downwards')
setattr(ssrd,'FillValue',-32767)
setattr(ssrd,'missing_value',-32767)
ssrd[:,:,:] = list(numpy.zeros((1461,73,144)))

strd = file.createVariable("strd",'f',('time', 'latitude', 'longitude'))
setattr(strd,'units','W m**-2 s')
setattr(strd,'long_name','Surface thermal radiation downwards')
setattr(strd,'FillValue',-32767)
setattr(strd,'missing_value',-32767)
strd[:,:,:] = list(numpy.zeros((1461,73,144)))

tp = file.createVariable("tp",'f',('time', 'latitude', 'longitude'))
setattr(tp,'units','m')
setattr(tp,'long_name','Total precipitation')
setattr(tp,'FillValue',-32767)
setattr(tp,'missing_value',-32767)
tp[:,:,:] = list(numpy.zeros((1461,73,144)))

file.close()
sys.stdout.write(" OK!\n")

sys.stdout.write("Renaming variables... ")
os.system('ncrename -a \'FillValue\',\'_FillValue\' test_forcing_I.nc')
sys.stdout.write(" OK!\n")
sys.stdout.write("Compressing file... ")
os.system('ncpdq -O -o test_forcing_I.nc test_forcing_I.nc')
sys.stdout.write(" OK!\n")
