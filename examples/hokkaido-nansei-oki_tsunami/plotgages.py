#!/usr/bin/env python

from __future__ import print_function

from fluidity_tools import stat_parser
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import getopt
import sys
import csv
import tools

def usage():
  print("plotgages.py --file=filename.detectors")

def main(argv=None):

  mea_filename='raw_data/WaveGages.csv'


  # Offsets:
  surface_heigth=-0.2
  d0=0.0005
  offset=-surface_heigth-d0

  try:
    opts, args = getopt.getopt(sys.argv[1:], "", ['file='])
  except getopt.GetoptError:
    print("Getopterror :(")
    usage()
    sys.exit(2)
  
  filename=''
  for opt, arg in opts:
    if opt == '--file':
      filename=arg
    elif opt=='-h' or opt=='--help':
        usage()
        sys.exit(2)
  if filename=='':
    print('No filename specified. You have to supply the detectos filename')
    usage()
    sys.exit(2)


  print("Generating plots")

  s = stat_parser(filename)

  timesteps=s["ElapsedTime"]["value"]
  timestep=timesteps[1]-timesteps[0]
  print("Found ", len(timesteps), " timesteps with dt=", timestep, " starting at t0=", timesteps[0]-timestep)
  
  fs=s["water"]["FreeSurface"]
  print("Found ", len(fs), " free surface detectors.")



  # fill in measurement data
  mea_gauge1=[]
  mea_gauge2=[]
  mea_gauge3=[]

  for i in range(0, len(timesteps)):
    gauges=tools.get_measurement(mea_filename, timesteps[i])
    mea_gauge1.append(gauges[0])
    mea_gauge2.append(gauges[1])
    mea_gauge3.append(gauges[2])
 
  plt.ion()  # switch in interactive mode
  fig1= figure()
#  fig2 = figure()
#  fig3 = figure()

  subplt1 = fig1.add_subplot(311, title='Gauge 1', xlabel='Time [s]', ylabel='Free surface [cm]')
  subplt2 = fig1.add_subplot(312, title='Gauge 2', xlabel='Time [s]', ylabel='Free surface [cm]')
  subplt3 = fig1.add_subplot(313, title='Gauge 3', xlabel='Time [s]', ylabel='Free surface [cm]')

   
  subplt1.plot(timesteps, s["water"]["FreeSurface"]["gauge1"]+offset, label='ICOM') # plot gauge1 detector data
  subplt1.plot(timesteps, mea_gauge1, label='Experimental data') # plot gauge1 measurement data

  subplt2.plot(timesteps, s["water"]["FreeSurface"]["gauge2"]+offset, label='ICOM') # plot gauge2 detector data
  subplt2.plot(timesteps, mea_gauge2, label='Experimental data') # plot gauge2 measurement data

  subplt3.plot(timesteps, s["water"]["FreeSurface"]["gauge3"]+offset, label='ICOM') # plot gauge3 detector data
  subplt3.plot(timesteps, mea_gauge3, label='Experimental data') # plot gauge3 measurement data

#  subplt1.xlabel('Time [s]')
  subplt1.legend()

  plt.draw()
#  for i in range(timesteps):
#    gauge1.append(s["water"]["FreeSurface"]["gauge1"])

  raw_input("Press Enter to exit")
   
if __name__ == "__main__":
   main()
