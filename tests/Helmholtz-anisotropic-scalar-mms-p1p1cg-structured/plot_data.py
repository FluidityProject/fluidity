#!/usr/bin/env python

import sys
import numpy as np
import pylab
import helmholtz_tools as h
from math import log

def plot_convergence(name1,logfile1,name2,logfile2):

  file1 = open(str(logfile1),'r').readlines()
  data1 = [float(line.split()[1]) for line in file1]
  iso1 = data1[0::3]
  aniso1 = data1[1::3]
  ml1 = data1[2::3]

  file2 = open(str(logfile2),'r').readlines()
  data2 = [float(line.split()[1]) for line in file2]
  iso2 = data2[0::3]
  aniso2 = data2[1::3]
  ml2 = data2[2::3]

  # x axis
  xlabels=['A:B','B:C','C:D','D:E']
  x=[1,2,3,4]
  ordertwo = [2,2,2,2]
  plot1 = pylab.figure(figsize = (6.5, 6.5))
  size = 15
  ax = pylab.subplot(111)
  ax.plot(x,iso1, linestyle='solid',color='red',lw=2)
  ax.plot(x,aniso1, linestyle='dashed',color='red',lw=2)
  ax.plot(x,ml1, linestyle='dashdot',color='red',lw=2)
  ax.plot(x,iso2, linestyle='solid',color='blue',lw=2)
  ax.plot(x,aniso2, linestyle='dashed',color='blue',lw=2)
  ax.plot(x,ml2, linestyle='dashdot',color='blue',lw=2)
  #ax.plot(x,ordertwo, linestyle='dashed',color='black')
  #pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name2)+'_iso',str(name2)+'_aniso'),loc="best")
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
  ax.set_xticklabels(xlabels)

  #pylab.axis([1,4,0.5,2.25])
  ax.set_xlabel('Mesh comparison', ha="center",fontsize=size)
  ax.set_ylabel('filter convergence rate',fontsize=size)
  pylab.savefig('convplot.eps')
  pylab.savefig('convplot.pdf')

  return

##################################

def plot_commerror(name1,logfile1,name2,logfile2):

  file1 = open(str(logfile1),'r').readlines()
  data1 = [float(line.split()[1]) for line in file1]
  iso1 = data1[0::3]
  aniso1 = data1[1::3]
  ml1 = data1[2::3]

  file2 = open(str(logfile2),'r').readlines()
  data2 = [float(line.split()[1]) for line in file2]
  iso2 = data2[0::3]
  aniso2 = data2[1::3]
  ml2 = data2[2::3]

  # x axis
  x=[8,16,32,64,128]
  xlabels=['A','B','C','D','E']
  orderone = [1,.5,.25,.125,.0625]
  ordertwo = [i**2 for i in orderone]
  plot1 = pylab.figure(figsize = (6, 6.5))
  size = 15
  ax = pylab.subplot(111)
  ax.plot(x,iso1, linestyle='solid',color='red',lw=2)
  ax.plot(x,aniso1, linestyle='dashed',color='red',lw=2)
  ax.plot(x,ml1, linestyle='dashdot',color='red',lw=2)
  ax.plot(x,iso2, linestyle='solid',color='blue',lw=2)
  ax.plot(x,aniso2, linestyle='dashed',color='blue',lw=2)
  ax.plot(x,ml2, linestyle='dashdot',color='blue',lw=2)
  #ax.plot(x,orderone, linestyle='solid',color='black',lw=2)
  #ax.plot(x,ordertwo, linestyle='dashed',color='black')
  #pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name2)+'_iso',str(name2)+'_aniso','order 1'),loc="best")
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

  # set axes to logarithmic
  pylab.gca().set_xscale('log',basex=2)
  pylab.gca().set_yscale('log',basex=2)

  pylab.axis([8,128,1.e-4,2])
  ax.set_xticks(x)
  ax.set_xticklabels(xlabels)

  #pylab.axis([1,5,1.e-6,1.])
  ax.set_xlabel('Mesh resolution', ha="center",fontsize=size)
  ax.set_ylabel('commutation error',fontsize=size)
  pylab.savefig('commerrorplot.eps')
  pylab.savefig('commerrorplot.pdf')

  return

##################################

def plot_commconv(name1,logfile1,name2,logfile2):

  file1 = open(str(logfile1),'r').readlines()
  data1 = [float(line.split()[1]) for line in file1]
  iso1 = data1[0::3]
  aniso1 = data1[1::3]
  ml1 = data1[2::3]

  file2 = open(str(logfile2),'r').readlines()
  data2 = [float(line.split()[1]) for line in file2]
  iso2 = data2[0::3]
  aniso2 = data2[1::3]
  ml2 = data2[2::3]

  # x axis
  xlabels=['A:B','B:C','C:D','D:E']
  x=[1,2,3,4]
  orderone = [1,1,1,1]
  plot1 = pylab.figure(figsize = (6.5, 6.5))
  size = 15
  ax = pylab.subplot(111)
  ax.plot(x,iso1, linestyle='solid',color='red',lw=2)
  ax.plot(x,aniso1, linestyle='dashed',color='red',lw=2)
  ax.plot(x,ml1, linestyle='dashdot',color='red',lw=2)
  ax.plot(x,iso2, linestyle='solid',color='blue',lw=2)
  ax.plot(x,aniso2, linestyle='dashed',color='blue',lw=2)
  ax.plot(x,ml2, linestyle='dashdot',color='blue',lw=2)
  #ax.plot(x,orderone, linestyle='dashed',color='black')
  #pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name2)+'_iso',str(name2)+'_aniso'),loc="best")
  pylab.legend((str(name1)+'_iso',str(name1)+'_aniso',str(name1)+'_ml',str(name2)+'_iso',str(name2)+'_aniso',str(name2)+'_ml'),loc="center")
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
  pylab.savefig('commconvplot.eps')
  pylab.savefig('commconvplot.pdf')

  return

def main():

  name1 = 'unstructured'
  name2 = 'structured'

  logfile1 = 'output_'+str(name1)+'.log'
  logfile2 = 'output_'+str(name2)+'.log'
  print('reading convergence data from', logfile1, logfile2)
  plot_convergence(name1,logfile1,name2,logfile2)

  logfile1 = 'commerror_'+str(name1)+'.log'
  logfile2 = 'commerror_'+str(name2)+'.log'
  print('reading commutation error data from', logfile1, logfile2)
  plot_commerror(name1,logfile1,name2,logfile2)

  logfile1 = 'commconv_'+str(name1)+'.log'
  logfile2 = 'commconv_'+str(name2)+'.log'
  print('reading commutation error convergence data from', logfile1, logfile2)
  plot_commconv(name1,logfile1,name2,logfile2)

if __name__ == "__main__":
    sys.exit(main())

