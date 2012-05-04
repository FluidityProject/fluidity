import struct
from scipy import interpolate

##############################################
#### 6 year physics from initialised LERM ####
def read_physics(persistent):
  """ Utility for reading the 6-year data set physics.bin
  """
  # make sure to open physics.bin with gzip!
  if not persistent.has_key('physics_fh'):
    import gzip
    from collections import deque
    persistent['physics_fh'] = gzip.open("tm6y_physics.bin", 'rb')
    persistent['physics_time'] = 0.0  
    persistent['physics_mld_history'] = deque()

  # Read MLD and number of layers recorded
  (mld, nlayers) = struct.unpack(">di", persistent['physics_fh'].read(12) )
  persistent['physics_mld'] = mld

  # Read depth, temperature and irradiance
  depth = []
  temperature = []
  irradiance = []
  for i in range(nlayers):
    (z, temp, irrad) = struct.unpack(">ddd", persistent['physics_fh'].read(24) )
    depth.append(z)
    temperature.append(temp)
    irradiance.append(irrad)

  # Create interpolation function for temperature and irradiance
  persistent['physics_temp_interp'] = interpolate.interp1d(depth, temperature, bounds_error=False, fill_value=temperature[-1])
  persistent['physics_irrad_interp'] = interpolate.interp1d(depth, irradiance, bounds_error=False, fill_value=irradiance[-1])

  # Advance time and record MLD 24-hour history
  persistent['physics_time'] = persistent['physics_time'] + 1800.
  persistent['physics_mld_history'].append(persistent['physics_mld'])
  if len(persistent['physics_mld_history']) > 48:
    persistent['physics_mld_history'].popleft()

def read_chemistry(persistent):
  """ Read initial chemical data from chemistry.bin
      Note: Nitrate will be set to exactly 0.0, which can create solver errors.
            Better set to something very small manually.
  """ 
  import gzip
  from numpy import arange
  fh = gzip.open("tm6y_chemistry.bin", 'rb')

  ammonium = []
  nitrate = []
  silicate = []
  for z in range(500):
    (amm, nit, sil) = struct.unpack(">ddd", fh.read(24) )
    ammonium.append( amm )
    silicate.append( sil )
    # Prevent solver NaNs
    if nit == 0.0:
      nit = 1.e-20
    nitrate.append( nit )

  # Depth
  depth = arange(0., 500., 1.)

  persistent['chem_ammonium'] = interpolate.interp1d(depth, ammonium, bounds_error=False, fill_value=ammonium[-1])
  persistent['chem_nitrate'] = interpolate.interp1d(depth, nitrate, bounds_error=False, fill_value=nitrate[-1])
  persistent['chem_silicate'] = interpolate.interp1d(depth, silicate, bounds_error=False, fill_value=silicate[-1])

#######################################
#### Initial chemistry of Toymodel ####
def read_chemistry_toy(persistent):
  """ Initialise chemistry for toymodel
  """
  fh = open("vew_tm_chem_init.csv", 'r')

  depth = []
  ammonium = []
  nitrate = []
  silicate = []
  for line in fh:
    data = line.split(",")
    depth.append( float(data[0]) )
    ammonium.append( float(data[1]) )
    nitrate.append( float(data[2]) )
    silicate.append( float(data[3]) )

  persistent['chem_ammonium'] = interpolate.interp1d(depth, ammonium, bounds_error=False, fill_value=ammonium[-1])
  persistent['chem_nitrate'] = interpolate.interp1d(depth, nitrate, bounds_error=False, fill_value=nitrate[-1])
  persistent['chem_silicate'] = interpolate.interp1d(depth, silicate, bounds_error=False, fill_value=silicate[-1])

##############################################
#### Initial Toymodel 2-year physics data ####
def read_physics_toy(persistent):
  """ 
  """

  if not persistent.has_key('physics_fh'):
    from collections import deque
    persistent['physics_fh'] = open("vew_tm_temp_irrad.csv" , 'r')
    persistent['physics_mld_fh'] = open("vew_tm_mld.csv", 'r')
    persistent['physics_mld_history'] = deque()

  # Read time and MLD
  data = persistent['physics_mld_fh'].readline().split(",")
  persistent['physics_time'] = float(data[0])
  persistent['physics_mld'] = float(data[1])

  # Read depth, temperature and irradiance
  data = persistent['physics_fh'].readline().split(",")
  #time = float(data[0])
  nlayers = int(data[1])
  depth = [ float(data[2]) ]
  temperature = [ float(data[3]) ]
  irradiance = [ float(data[4]) ]
  for i in range(nlayers-1):
    data = persistent['physics_fh'].readline().split(",")
    depth.append( float(data[2]) )
    temperature.append( float(data[3]) )
    irradiance.append( float(data[4]) )

  # Create interpolation function for temperature and irradiance
  persistent['physics_temp_interp'] = interpolate.interp1d(depth, temperature, bounds_error=False, fill_value=temperature[-1])
  persistent['physics_irrad_interp'] = interpolate.interp1d(depth, irradiance, bounds_error=False, fill_value=irradiance[-1])

  # Record MLD 24-hour history
  persistent['physics_mld_history'].append(persistent['physics_mld'])
  if len(persistent['physics_mld_history']) > 48:
    persistent['physics_mld_history'].popleft()
  

######################################
#### 20 years of Toymodel physics ####
def read_physics_toy20(persistent):
  """ Utility for reading the 20-year data set toy_physics.bin
  """
  from numpy import arange, append

  if not persistent.has_key('physics_fh'):
    from collections import deque
    persistent['physics_fh'] = open("toy_physics.bin", 'rb')
    persistent['physics_time'] = 0.0  
    persistent['physics_mld_history'] = deque()

  (mld, ) = struct.unpack(">d", persistent['physics_fh'].read(8) )
  persistent['physics_mld'] = mld

  temperature = []
  irradiance = []
  for z in range(520):
    temp, irrad = struct.unpack(">dd", persistent['physics_fh'].read(16))
    temperature.append( temp )
    irradiance.append( irrad )

  # Depth
  depth = arange(0., 1. ,0.05)
  depth = append(depth, arange(1., 501., 1.) )

  # Create interpolation function for temperature and irradiance
  persistent['physics_temp_interp'] = interpolate.interp1d(depth, temperature, bounds_error=False, fill_value=temperature[-1])
  persistent['physics_irrad_interp'] = interpolate.interp1d(depth, irradiance, bounds_error=False, fill_value=irradiance[-1])

  # Advance time and record MLD 24-hour history
  persistent['physics_time'] = persistent['physics_time'] + 1800.
  persistent['physics_mld_history'].append(persistent['physics_mld'])
  if len(persistent['physics_mld_history']) > 48:
    persistent['physics_mld_history'].popleft()


