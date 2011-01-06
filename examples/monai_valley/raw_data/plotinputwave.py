#!/usr/bin/python

from fluidity_tools import stat_parser
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import getopt
import sys
import csv

def usage():
  print "plotinputwave.py -b starttime -e endtime --save=basename"

def get_inputelevation(t):
        InputWaveReader = csv.reader(open('InputWave.csv', 'rb'), delimiter='\t')
        data=[]
        for (time, heigth) in InputWaveReader:
                data.append((float(time), float(heigth)))
                
        for i in range(1,len(data)):
                if data[i][0]<t:
                        continue
                t1=data[max(0,i-1)][0]
                t2=data[i][0]
                h1=data[max(0,i-1)][1]
                h2=data[i][1]
                return h1*(t-t2)/(t1-t2)+h2*(t-t1)/(t2-t1)
        
        print "Warning: t is outside the available data. Using last available waterheigth..."
        return data[-1][1]

def main(argv=None):

  dt=0.05 # use same timestep than in csv file

  try:
    opts, args = getopt.getopt(sys.argv[1:], "t:e:b:", ['save='])
  except getopt.GetoptError:
    print "Getopterror :("
    usage()
    sys.exit(2)
  
  subtitle=''
  subtitle_pure=''
  endtime=22.5
  starttime=0.0
  save=False

  for opt, arg in opts:
    if opt == '--save':
        save=True
        savename=arg
    elif opt=='-h' or opt=='--help':
        usage()
        sys.exit(2)
    elif opt=='-t':
        subtitle=', '+arg
        subtitle_pure=arg
    elif opt=='-b':
        starttime=float(arg)
    elif opt=='-e':
        endtime=float(arg)

  print "Generating plot"

  print 'Using dt=', dt
  
  starttimestep=int(max(0,starttime/dt))
  endtimestep=int(endtime/dt)

  print 'starttimestep=', starttimestep
  print 'endtimestep=', endtimestep

  # fill in measurement data
  input_elevation=[]
  time=[]

  for i in range(starttimestep, endtimestep):
    time.append(i*dt)
    elev=get_inputelevation(time[-1])
    input_elevation.append(elev*100.0) # in cm
 
  plt.ion()  # switch in interactive mode
  fig1= figure()

  subplt1 = fig1.add_subplot(111, xlabel='Time [s]', ylabel='Water level [cm]')

  subplt1.plot(time, input_elevation) # plot gauge1 detector data

  if not save:
    plt.draw()
    raw_input("Press Enter to exit")
  else:
    plt.savefig(savename+'.pdf', facecolor='white', edgecolor='black', dpi=100)
    print 'Saved to '+savename+'.pdf'

#  for i in range(timesteps):
#    gauge1.append(s["water"]["FreeSurface"]["gauge1"])

    
if __name__ == "__main__":
   main()
