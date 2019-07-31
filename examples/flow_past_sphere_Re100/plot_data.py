#!/usr/bin/env python3

import math
import fluidity_tools
import os
import pylab

##############################

def xforce_lab(Re):
  # Correlation to lab data from Brown and Lawler (2003)
  Cd = (24.0/Re)*(1.0 + 0.15*(Re**0.681)) + 0.407/(1.0 + (8710.0/Re))
  # force in x-direction
  xforce = Cd*(0.5*1.0*1.0*math.pi/4.0)  
  return xforce

##############################

def xforce_stat(sf):
  stat = fluidity_tools.stat_parser(sf)
  time_avg_xforce = pylab.average(stat["0"]["Velocity"]["force%1"][-20:])
  return time_avg_xforce

##############################  
  
Re_no = 100
Re_nos_lab = [0.1+(i/100.) for i in range(0,1000000,100)]
Re_nos_lab.append(1E4)

pylab.figure(num=1, figsize = (16.5, 11.5))
pylab.plot(Re_nos_lab,[xforce_lab(Re) for Re in Re_nos_lab], 'k-', linewidth = 3)

sf = './flow_past_sphere_Re'+str(Re_no)+'.stat'
pylab.plot([Re_no],[xforce_stat(sf)],'ko', mfc = 'none', ms = 18, mew = 3)

fs = 18
pylab.gca().set_xscale('log')
pylab.gca().set_yscale('log')
pylab.xlabel('Reynolds number', fontsize = fs+4)
pylab.ylabel('$C_D$', fontsize = fs+4)
for label in pylab.gca().get_axes().get_xticklabels(): label.set_fontsize(fs)
for label in pylab.gca().get_axes().get_yticklabels(): label.set_fontsize(fs)   
pylab.savefig('Sphere_drag.pdf')
pylab.show()
