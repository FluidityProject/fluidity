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


################################################
#--------------- FROUDE NUMBER ----------------# 
################################################

def Froudenumber(flmlname):
  print("\n********** Calculating the Froude number\n")
  # warn user about assumptions
  print("Froude number calculations makes three assumptions: \n i) domain height = 0.1m \n ii) mid point domain is at x = 0.4 \n iii) initial temperature difference is 1.0 degC")
  domainheight = 0.1
  domainmid = 0.4
  rho_zero, T_zero, alpha, g = le_tools.Getconstantsfromflml(flmlname)
  gprime = rho_zero*alpha*g*1.0 # this has assumed the initial temperature difference is 1.0 degC

  # get list of vtus
  filelist = le_tools.GetFiles('./')
  logs = ['diagnostics/logs/time.log','diagnostics/logs/X_ns.log','diagnostics/logs/X_fs.log']
  try:
    # if have extracted information already just use that
    os.stat('diagnostics/logs/time.log')
    os.stat('diagnostics/logs/X_ns.log')
    os.stat('diagnostics/logs/X_fs.log')
    time = le_tools.ReadLog('diagnostics/logs/time.log')
    X_ns = [x-domainmid for x in le_tools.ReadLog('diagnostics/logs/X_ns.log')]
    X_fs = [domainmid-x for x in le_tools.ReadLog('diagnostics/logs/X_fs.log')]
  except OSError:
    # otherwise get X_ns and X_fs and t from vtus
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

    # shift so bot X_ns and X_fs are 
    # distance of front from 
    #initial position (mid point of domain)
    X_ns = [x-domainmid for x in X_ns]
    X_fs = [domainmid-x for x in X_fs]

  # Calculate U_ns and U_fs from X_ns, X_fs and t
  U_ns = le_tools.GetU(time, X_ns)
  U_fs = le_tools.GetU(time, X_fs)
  U_average = [[],[]]

  # If possible average 
  # (if fronts have not travelled far enough then will not average)
  start_val, end_val, average_flag_ns = le_tools.GetAverageRange(X_ns, 0.2, domainheight)
  if average_flag_ns == True: U_average[0].append(pylab.average(U_ns[start_val:end_val]))
  
  start_val, end_val, average_flag_fs = le_tools.GetAverageRange(X_fs, 0.25, domainheight)
  if average_flag_fs == True: U_average[1].append(pylab.average(U_fs[start_val:end_val]))
  
  # plot
  fs = 18
  pylab.figure(num=1, figsize = (16.5, 11.5))
  pylab.suptitle('Front speed', fontsize = fs)

  pylab.subplot(221)
  pylab.plot(time,X_ns, color = 'k')
  pylab.axis([0,45,0,0.4])
  pylab.grid('on')
  pylab.xlabel('$t$ (s)', fontsize = fs)
  pylab.ylabel('$X$ (m)', fontsize = fs)
  pylab.title('no-slip', fontsize = fs)
    
  pylab.subplot(222)
  pylab.plot([x/domainheight for x in X_ns],[U/math.sqrt(gprime*domainheight) for U in U_ns], color = 'k')
  pylab.axis([0,4,0,0.6])
  pylab.grid('on')
  pylab.axhline(0.406, color = 'k')
  pylab.axhline(0.432, color = 'k')
  pylab.text(3.95,0.396,'Hartel 2000',bbox=dict(facecolor='white', edgecolor='black'), va = 'top', ha = 'right')
  pylab.text(3.95,0.442,'Simpson 1979',bbox=dict(facecolor='white', edgecolor='black'), ha = 'right')
  pylab.xlabel('$X/H$', fontsize = fs)
  pylab.ylabel('$Fr$', fontsize = fs)
  pylab.title('no-slip', fontsize = fs)
  if average_flag_ns == True:
    pylab.axvline(2.0, color = 'k')
    pylab.axvline(3.0, color = 'k')
    pylab.text(0.05, 0.01, 'Average Fr = '+'{0:.2f}'.format(U_average[0][0]/math.sqrt(gprime*domainheight))+'\nvertical lines indicate the range \nover which the average is taken', bbox=dict(facecolor='white', edgecolor='black'))
  
  pylab.subplot(223)
  pylab.plot(time,X_fs, color = 'k')
  pylab.axis([0,45,0,0.4])
  pylab.grid('on')
  pylab.xlabel('$t$ (s)', fontsize = fs)
  pylab.ylabel('$X$ (m)', fontsize = fs)
  pylab.title('free-slip', fontsize = fs)
    
  pylab.subplot(224)
  pylab.plot([x/domainheight for x in X_fs],[U/math.sqrt(gprime*domainheight) for U in U_fs], color = 'k')
  pylab.axis([0,4,0,0.6])
  pylab.grid('on')
  pylab.axhline(0.477, color = 'k')
  pylab.text(3.95,0.467,'Hartel 2000', va = 'top',bbox=dict(facecolor='white', edgecolor='black'), ha = 'right')
  pylab.xlabel('$X/H$', fontsize = fs)
  pylab.ylabel('$Fr$', fontsize = fs)
  pylab.title('free-slip', fontsize = fs)
  if average_flag_fs == True:
    pylab.text(0.05, 0.01, 'Average Fr  = '+'{0:.2f}'.format(U_average[1][0]/math.sqrt(gprime*domainheight))+'\nvertical lines indicate the range \nover which the average is taken', bbox=dict(facecolor='white', edgecolor='black'))
    pylab.axvline(2.5, color = 'k')
    pylab.axvline(3.0, color = 'k')

  pylab.savefig('diagnostics/plots/front_speed.png')
  return

################################################
#------------------ MIXING --------------------# 
################################################

def mixing(flmlname):

  print("\n********** Calculating the mixing diagnostics\n")
  # warn user about assumptions
  print("Background potential energy calculations makes two assumptions: \n i) domain height = 0.1m \n ii) initial temperature difference is 1.0 degC")
  domainheight = 0.1
  rho_zero, T_zero, alpha, g = le_tools.Getconstantsfromflml(flmlname)

  # get mixing bin bounds and remove lower bound (=-\infty)
  bounds = le_tools.Getmixingbinboundsfromflml(flmlname)[1:]

  # find indicies of selected bounds for plotting
  index_plot = []
  for b in [-0.5,-0.25, 0.0,0.25,0.5]:
    index_plot.append(pylab.find(numpy.array([abs(val - b) for val in bounds]) < 1E-6)[0]) 
  
  time = []
  volume_fraction = []
  reference_state = []
  bpe = [] 
  
  # get stat files 
  # time_index_end used to ensure don't repeat values
  stat_files, time_index_end = le_tools.GetstatFiles('./')
  for i in range(len(stat_files)):
    stat = stat_parser(stat_files[i])
    for j in range(time_index_end[i]):
      time.append(stat['ElapsedTime']['value'][j])
      bins = stat['fluid']['Temperature']['mixing_bins%cv_normalised'][:,j]
      # rearrange bins so have nobins = nobounds -1
      # amounts to including any undershoot or overshoots in lower/upper most bin
      # for discussion of impacts see H. Hiester, PhD thesis (2011), chapter 4.
      bins[1] = bins[0]+bins[1]
      bins[-2] = bins[-2]+bins[-1]
      bins = bins[1:-1]
      
      # sum up bins for plot
      volume_fraction.append(tuple([sum(bins[index_plot[k]:index_plot[k+1]]) for k in range(len(index_plot)-1)]))
      
      # get reference state using method of Tseng and Ferziger 2001
      Abins = sum([bins[k]*(bounds[k+1]-bounds[k]) for k in range(len(bins))])
      pdf = [val/Abins for val in bins]
      rs = [0]
      for k in range(len(pdf)): rs.append(rs[-1]+(domainheight*pdf[k]*(bounds[k+1]-bounds[k])))
      reference_state.append(tuple(rs))
      
      # get background potential energy, 
      # noting \rho = \rho_zero(1-\alpha(T-T_zero))
      # and reference state is based on temperature
      # bpe_bckgd = 0.5*(g*rho_zero*(1.0+(alpha*T_zero)))*(domainheight**2)
      # but don't include this as will look at difference over time
      bpe.append(-rho_zero*alpha*g*scipy.integrate.trapz(x=reference_state[-1],y=[bounds[j]*reference_state[-1][j] for j in range(len(reference_state[-1]))]))
  
  volume_fraction = numpy.array(volume_fraction)
  reference_state = numpy.array(reference_state)
  
  bpe_zero = bpe[0]
  bpe = [val - bpe_zero for val in bpe]

  # plot
  fs = 18
  pylab.figure(num=2, figsize = (16.5, 11.5))
  pylab.suptitle('Mixing', fontsize = fs)
  
  # volume fraction
  pylab.subplot(221)
  pylab.plot(time,volume_fraction[:,0], label = '$T < -0.25$', color = 'k')
  pylab.plot(time,volume_fraction[:,1], label = '$-0.25 < T < 0.0$', color = 'g')
  pylab.plot(time,volume_fraction[:,2], label = '$0.0 < T < 0.25$', color = 'b')
  pylab.plot(time,volume_fraction[:,3], label = '$0.25 < T$', color = '0.5')
  
  pylab.axis([0,time[-1],0,0.5])
  pylab.legend(loc = 0)
  pylab.grid('on')
  pylab.xlabel('$t$ (s)', fontsize = fs)
  pylab.ylabel('$V/|\\Omega|$', fontsize = fs)
  pylab.title('Volume fraction', fontsize = fs) 
  
  # reference state contours
  pylab.subplot(222)
  for i in index_plot: pylab.plot(time, reference_state[:,i],  color = 'k')
  pylab.text(time[-1]/100, 1.5E-3, 'From bottom to top contours correspond to values \n $T = -0.5, \, -0.25, \, 0.0, \, 0.25, \, 0.5$ \nwhere the values for $T=-0.5$ and $0.5$ take the values\n$z_* = 0.0$ and $0.1$ respectively', bbox=dict(facecolor='white', edgecolor='black'))
  pylab.axis([0,time[-1],0,domainheight])
  pylab.grid('on')
  pylab.xlabel('$t$ (s)', fontsize = fs)
  pylab.ylabel('$z_*$ (m)', fontsize = fs)
  pylab.title('Reference state', fontsize = fs) 
  
  pylab.subplot(223)
  pylab.plot(bounds, reference_state[-1], color = 'k')
  pylab.grid('on')
  pylab.axis([-0.5,0.5,0,domainheight])
  pylab.xlabel('$T$ ($^\\circ$C)', fontsize = fs)
  pylab.ylabel('$z_*$ (m)', fontsize = fs)
  pylab.title('Reference state at $t='+str(time[-1])+'\\,$s', fontsize = fs) 
  
  pylab.subplot(224)
  pylab.plot(time, bpe, color = 'k')
  pylab.grid('on')
  pylab.gca().get_xaxis().get_axes().set_xlim(0.0,time[-1])
  pylab.xlabel('$t$ (s)', fontsize = fs)
  pylab.ylabel('$\\Delta E_b$', fontsize = fs-2)
  pylab.gca().get_yaxis().set_major_formatter(pylab.FormatStrFormatter('%1.1e'))
  pylab.title('Background potential energy', fontsize = fs) 
  
  pylab.savefig('diagnostics/plots/mixing.png')
  return

################################################
#---------------- MAIN CALLS ------------------# 
################################################

# make directories to put diagnostics in (if they don't exist already)
for name in ['diagnostics', 'diagnostics/logs', 'diagnostics/plots']:
  try: os.stat(name)
  except OSError: os.mkdir(name)
    
# set flmlname
flmlname = './lock_exchange.flml'

# Get and plot Froude number, mixing bins and background potential energy
Froudenumber(flmlname)  
mixing(flmlname)

# show plots and tell user where to find a copy
pylab.show()
print('The images have also been saved in ./diagnostics/plots')
