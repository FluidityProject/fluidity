import vtk
import glob
import sys
import os
import scipy.stats
import pylab
import math
import numpy
from fluidity_tools import stat_parser
import le_tools


def Froudenumber():
  print "Froude number calculations make two assumptions: \n i) domain height = 0.1m \n ii) g\'= 10^-2"
  domainheight = 0.1
  gprime = 1E-2

  for name in ['diagnostics', 'diagnostics/logs', 'diagnostics/plots']:
    try: os.stat(name)
    except OSError: os.mkdir(name)

  filelist = le_tools.GetFiles('./')
  logs = ['diagnostics/logs/time.log','diagnostics/logs/X_ns.log','diagnostics/logs/X_fs.log']
  try:
    os.stat('diagnostics/logs/time.log')
    os.stat('diagnostics/logs/X_ns.log')
    os.stat('diagnostics/logs/X_fs.log')
    time = le_tools.ReadLog('diagnostics/logs/time.log')
    X_ns = [x-0.4 for x in le_tools.ReadLog('diagnostics/logs/X_ns.log')]
    X_fs = [0.4-x for x in le_tools.ReadLog('diagnostics/logs/X_fs.log')]
  except OSError:
    time, X_ns, X_fs = le_tools.GetXandt(filelist)
    f_time = open('./diagnostics/logs/time.log','w')
    for t in time: f_time.write(str(t)+'\n')
    f_time.close()
    f_X_ns = open('./diagnostics/logs/X_ns.log','w')
    for X in X_ns: f_X_ns.write(str(X)+'\n')
    f_X_ns.close()
    f_X_fs = open('./diagnostics/logs/X_fs.log','w')
    for X in X_fs: f_X_fs.write(str(X)+'\n')
    f_X_fs.close()

    X_ns = [x-0.4 for x in X_ns]
    X_fs = [0.4-x for x in X_fs]

  U_ns = le_tools.GetU(time, X_ns)
  U_fs = le_tools.GetU(time, X_fs)
  U_average = []

  start_val, end_val = le_tools.GetAverageRange(X_ns, 0.2, domainheight)
  U_average.append(pylab.average(U_ns[start_val:end_val]))
  
  start_val, end_val = le_tools.GetAverageRange(X_fs, 0.25, domainheight)
  U_average.append(pylab.average(U_fs[start_val:end_val]))
    
  pylab.figure(num=1, figsize = (16.5, 11.5))
  pylab.suptitle('Front speed')

  pylab.subplot(221)
  pylab.plot(time,X_ns)
  pylab.axis([0,45,0,0.4])
  pylab.grid("True")
  pylab.xlabel('time (s)')
  pylab.ylabel('X (m)')
  pylab.title('no-slip')
    
  pylab.subplot(222)
  pylab.plot([x/domainheight for x in X_ns],[U/math.sqrt(gprime*domainheight) for U in U_ns])
  pylab.axis([0,4,0,0.6])
  pylab.grid("True")
  pylab.text(0.05, 0.565, 'Average Fr (for 2.0<X/H<3.0) = '+str(U_average[0]/math.sqrt(gprime*domainheight)) , bbox=dict(facecolor='white', edgecolor='black'))
  pylab.xlabel('X/H')
  pylab.ylabel('Fr')
  pylab.title('no-slip')
  
  pylab.subplot(223)
  pylab.plot(time,X_fs)
  pylab.axis([0,45,0,0.4])
  pylab.grid("True")
  pylab.xlabel('time (s)')
  pylab.ylabel('X (m)')
  pylab.title('free-slip')
    
  pylab.subplot(224)
  pylab.plot([x/domainheight for x in X_fs],[U/math.sqrt(gprime*domainheight) for U in U_fs])
  pylab.axis([0,4,0,0.6])
  pylab.grid("True")
  pylab.text(0.05, 0.565, 'Average Fr (for 2.5<X/H<3.0) = '+str(U_average[1]/math.sqrt(gprime*domainheight)) , bbox=dict(facecolor='white', edgecolor='black'))
  pylab.xlabel('X/H')
  pylab.ylabel('Fr')
  pylab.title('free-slip')

  pylab.savefig('diagnostics/plots/front_speed.png')
  
  return
  
def mixing_stats():
  stat = stat_parser('lock_exchange.stat')
  
  pylab.figure(num=2, figsize = (16.5, 11.5))
  pylab.suptitle('Mixing')
  
  pylab.plot(stat['ElapsedTime']['value'],stat['fluid']['Temperature']['mixing_bins%cv_normalised'][0], label = 'T < -0.4')
  pylab.plot(stat['ElapsedTime']['value'],stat['fluid']['Temperature']['mixing_bins%cv_normalised'][1], label = '-0.4 < T < 0.4')
  pylab.plot(stat['ElapsedTime']['value'],stat['fluid']['Temperature']['mixing_bins%cv_normalised'][2], label = '0.4 < T')
  
  
  time = le_tools.ReadLog('diagnostics/logs/time.log')
  X_ns = [x-0.4 for x in le_tools.ReadLog('diagnostics/logs/X_ns.log')]
  X_fs = [0.4-x for x in le_tools.ReadLog('diagnostics/logs/X_fs.log')]
  try:
    index = pylab.find(numpy.array(X_ns)<0.4-1E-3)[-1]
    pylab.fill_between([time[index],time[index+1]],[0,0],[0.5,0.5], color = '0.3') 
    index = pylab.find(numpy.array(X_fs)<0.4-1E-3)[-1]
    pylab.fill_between([time[index],time[index+1]],[0,0],[0.5,0.5], color = '0.6')
  except IndexError: 
    print 'not plotting shaded regions on mixing plot as front has not reached end wall'
  
  #pylab.text(0.05, 0.565, 'shaded regions show when fronts near the end wall (dark: free-slip), (light: no-slip)', bbox=dict(facecolor='white', edgecolor='black', width = 25))
  
  pylab.grid("True")
  pylab.legend(loc=0)
  pylab.xlabel('time (s)')
  pylab.ylabel('domain fraction')
  
  pylab.savefig('diagnostics/plots/mixing.png')

  return
  
Froudenumber()  
mixing_stats()

pylab.show()
