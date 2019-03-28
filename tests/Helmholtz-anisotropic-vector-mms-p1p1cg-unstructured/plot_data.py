#!/usr/bin/env python

import sys
import numpy as np
import pylab
import helmholtz_tools as h
from math import log

def plot_error(name1,logfile1,name2,logfile2):

  file1 = open(str(logfile1),'r').readlines()
  data1 = [float(line.split()[1]) for line in file1]
  iso1 = data1[0::2]
  aniso1 = data1[1::2]

  file2 = open(str(logfile2),'r').readlines()
  data2 = [float(line.split()[1]) for line in file2]
  aniso2 = data2[0::4]
  iso2 = data2[0::2]
  aniso2 = data2[1::2]

  # x axis
  x=[1,2,3,4,5]
  xlabels=['A','B','C','D','E']
  orderone = [1,.5,.25,.125,.0625]
  ordertwo = [i**2 for i in orderone]
  plot1 = pylab.figure(figsize = (6, 6))
  size = 15
  ax = pylab.subplot(111)
  ax.plot(x,iso1, linestyle='solid',color='red',lw=2)
  ax.plot(x,aniso1, linestyle='dashed',color='red',lw=2)
  #ax.plot(x,ml1, linestyle='dashdot',color='red',lw=2)
  #ax.plot(x,vel1, linestyle='dotted',color='red',lw=2)
  ax.plot(x,iso2, linestyle='solid',color='blue',lw=2)
  ax.plot(x,aniso2, linestyle='dashed',color='blue',lw=2)
  #ax.plot(x,ml2, linestyle='dashdot',color='blue',lw=2)
  #ax.plot(x,vel2, linestyle='dotted',color='blue',lw=2)
  #ax.plot(x,ordertwo, linestyle='solid',color='black',lw=2)
  #pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name1)+'_ml',str(name1),str(name2)+'_iso',str(name2)+'_aniso',str(name2)+'_ml',str(name2),'second order'),loc="best")
  pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name2)+'_iso',str(name2)+'_aniso'),loc="best")
  leg = pylab.gca().get_legend()
  ltext = leg.get_texts()
  pylab.setp(ltext, fontsize = size, color = 'black')
  frame=leg.get_frame()
  frame.set_fill(False)
  frame.set_visible(False)

  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  # set y-axis to logarithmic
  pylab.gca().set_yscale('log',basex=2)

  ax.set_xticks(x)
  ax.set_xticklabels(xlabels)

  #pylab.axis([1,5,1.e-6,1.])
  ax.set_xlabel('Mesh resolution', ha="center",fontsize=size)
  ax.set_ylabel('filter error L2 norm',fontsize=size)
  pylab.savefig('commerrorplot-vector.pdf')
  pylab.savefig('commerrorplot-vector.eps')
  return

def plot_convergence(name1,logfile1,name2,logfile2):

  file1 = open(str(logfile1),'r').readlines()
  data1 = [float(line.split()[1]) for line in file1]
  aniso1 = data1[0::4]
  iso1 = data1[1::4]
  ml1 = data1[2::4]
  vel1 = data1[3::4]

  file2 = open(str(logfile2),'r').readlines()
  data2 = [float(line.split()[1]) for line in file2]
  aniso2 = data2[0::4]
  iso2 = data2[1::4]
  ml2 = data2[2::4]
  vel2 = data2[3::4]

  # x axis
  xlabels=['A:B','B:C','C:D','D:E']
  x=[1,2,3,4]
  ordertwo = [2,2,2,2]
  plot1 = pylab.figure(figsize = (6, 6.5))
  size = 15
  ax = pylab.subplot(111)
  ax.plot(x,iso1, linestyle='solid',color='red',lw=2)
  ax.plot(x,aniso1, linestyle='dashed',color='red',lw=2)
  ax.plot(x,ml1, linestyle='dashdot',color='red',lw=2)
  #ax.plot(x,vel1, linestyle='dotted',color='red',lw=2)
  ax.plot(x,iso2, linestyle='solid',color='blue',lw=2)
  ax.plot(x,aniso2, linestyle='dashed',color='blue',lw=2)
  ax.plot(x,ml2, linestyle='dashdot',color='blue',lw=2)
  #ax.plot(x,vel2, linestyle='dotted',color='blue',lw=2)
  #ax.plot(x,ordertwo, linestyle='solid',color='black',lw=2)
  #pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name1)+'_ml',str(name1),str(name2)+'_iso',str(name2)+'_aniso',str(name2)+'_ml',str(name2),'second order'),loc="best")
  pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name1)+'_ml',str(name2)+'_iso',str(name2)+'_aniso',str(name2)+'_ml'),loc="best")
  leg = pylab.gca().get_legend()
  ltext = leg.get_texts()
  pylab.setp(ltext, fontsize = size, color = 'black')
  frame=leg.get_frame()
  frame.set_fill(False)
  frame.set_visible(False)

  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  # This sets the ticks on the x axis to be exactly where we put
  # the center of the bars.
  ax.set_xticks(x)

  # Set the x tick labels to the group_labels defined above.
  ax.set_xticklabels(xlabels)
 
  # Extremely nice function to auto-rotate the x axis labels.
  # It was made for dates (hence the name) but it works
  # for any long x tick labels
  #plot1.autofmt_xdate()

  #pylab.axis([1,4,0.,2.75])
  ax.set_xlabel('Mesh comparison', ha="center",fontsize=size)
  ax.set_ylabel('filter convergence rate',fontsize=size)
  pylab.savefig('errorconvplot-vector.pdf')
  pylab.savefig('errorconvplot-vector.eps')
  return

def plot_commconv(name1,logfile1,name2,logfile2):

  file1 = open(str(logfile1),'r').readlines()
  data1 = [float(line.split()[1]) for line in file1]
  iso1 = data1[0::2]
  aniso1 = data1[1::2]

  file2 = open(str(logfile2),'r').readlines()
  data2 = [float(line.split()[1]) for line in file2]
  iso2 = data2[0::2]
  aniso2 = data2[1::2]

  # x axis
  xlabels=['A:B','B:C','C:D','D:E']
  x=[1,2,3,4]
  orderone = [1,1,1,1]
  plot1 = pylab.figure(figsize = (6, 6.5))
  size = 15
  ax = pylab.subplot(111)
  ax.plot(x,iso1, linestyle='solid',color='red',lw=2)
  ax.plot(x,aniso1, linestyle='dashed',color='red',lw=2)
  ax.plot(x,iso2, linestyle='solid',color='blue',lw=2)
  ax.plot(x,aniso2, linestyle='dashed',color='blue',lw=2)
  #ax.plot(x,orderone, linestyle='dashed',color='black',lw=2)
  pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name2)+'_iso',str(name2)+'_aniso'),loc="best")
  leg = pylab.gca().get_legend()
  ltext = leg.get_texts()
  pylab.setp(ltext, fontsize = size, color = 'black')
  frame=leg.get_frame()
  frame.set_fill(False)
  frame.set_visible(False)

  #ax.grid("True")
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  # This sets the ticks on the x axis to be exactly where we put
  # the center of the bars.
  ax.set_xticks(x)
  ax.set_xticklabels(xlabels)
 
  #pylab.axis([1,4,0.,1.])
  ax.set_xlabel('Mesh comparison', ha="center",fontsize=size)
  ax.set_ylabel('commutation error convergence rate',fontsize=size)
  pylab.savefig('commconvplot-vector.pdf')
  pylab.savefig('commconvplot-vector.eps')
  return

def main():

  name1 = 'unstructured'
  name2 = 'structured'
  logfile1 = 'commerror_'+str(name1)+'.log'
  logfile2 = 'commerror_'+str(name2)+'.log'
  print('reading error data from', logfile1, logfile2)
  plot_error(name1,logfile1,name2,logfile2)

  logfile1 = 'errorconv_'+str(name1)+'.log'
  logfile2 = 'errorconv_'+str(name2)+'.log'
  print('reading convergence data from', logfile1, logfile2)
  plot_convergence(name1,logfile1,name2,logfile2)

  logfile1 = 'commconv_'+str(name1)+'.log'
  logfile2 = 'commconv_'+str(name2)+'.log'
  print('reading convergence data from', logfile1, logfile2)
  plot_commconv(name1,logfile1,name2,logfile2)

if __name__ == "__main__":
    sys.exit(main())

