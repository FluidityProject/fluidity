from __future__ import print_function

from fluidity_tools import stat_parser
import sys
import csv

# does a linear interpolation of the input data
def get_measurement(mea_filename, t):
  reader = csv.reader(open(mea_filename, 'rb'), delimiter=',')
 
  data=[]
  for (time, gauge1, gauge2, gauge3) in reader:
    data.append((float(time), float(gauge1), float(gauge2), float(gauge3)))

  for i in range(1,len(data)):
    if data[i][0]<t:
      continue
    # interpolate the three gauge stations
    t1=data[max(0,i-1)][0]
    t2=data[i][0]
    gauges=[]
    for g in range(3):
      h1=data[max(0,i-1)][g+1]/100 # Convert from cm to m
      h2=data[i][g+1]/100 # Convert from cm to m
      gauges.append(h1*(t-t2)/(t1-t2)+h2*(t-t1)/(t2-t1))
    return gauges

  print("Warning: simulation time t=", t, " is outside the available data (", data[0][0], ", ", data[-1][0], "). Using last available waterheigth...")
  return [data[-1][1], data[-1][2], data[-1][2]]


def gage_error_integral(detector_filename):

  mea_filename='raw_data/WaveGages.csv'
 
  s = stat_parser(detector_filename)
  timesteps=s["ElapsedTime"]["value"]
  timestep=timesteps[1]-timesteps[0]
  print("Found ", len(timesteps), " timesteps with dt=", timestep, " starting at t0=", timesteps[0]-timestep)
  
  fs=s["water"]["FreeSurface"]
  print("Found ", len(fs), " free surface detectors.")

  error_integral=[0.0, 0.0, 0.0]

  # First timestep
  gauges=get_measurement(mea_filename, timesteps[0])
  error_integral[0]+=0.5*abs(gauges[0])*timestep
  error_integral[1]+=0.5*abs(gauges[1])*timestep
  error_integral[2]+=0.5*abs(gauges[2])*timestep
  for i in range(1, len(timesteps)-1):
    gauges=get_measurement(mea_filename, timesteps[i])
    error_integral[0]+=abs(gauges[0])*timestep
    error_integral[1]+=abs(gauges[1])*timestep
    error_integral[2]+=abs(gauges[2])*timestep
  # Last timestep
  gauges=get_measurement(mea_filename, timesteps[-1])
  error_integral[0]+=0.5*abs(gauges[0])*timestep
  error_integral[1]+=0.5*abs(gauges[1])*timestep
  error_integral[2]+=0.5*abs(gauges[2])*timestep
 
  return error_integral
        

