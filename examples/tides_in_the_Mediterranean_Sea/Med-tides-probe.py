#!/usr/bin/env python

from __future__ import print_function

import vtktools
import math
from numpy import array 

u=vtktools.vtu("tidesmedsea-flat.vtu")
g = open("Med-GEBCO-5m-gauges-fes2004-O1-102", "w")

pts=vtktools.arr([
#[-5.3500, 36.1333, -2.00],
#[-4.4500, 36.7000, 0.00],
#[-3.9167, 35.2500, 0.00],
#[-2.4500, 36.8333, 0.00],
#[-0.5833, 38.3333, 0.00],
[2.6333, 39.5833, 0.00],
#[3.1000, 42.4833, 0.00],
#[5.3500, 43.3000, 0.00],
#[6.9167, 36.8833, 0.00],
[8.0167, 43.8667, 0.00],
[8.3000, 39.1500, 0.00],
[8.9000, 44.4000, 0.00],
[9.1667, 39.2000, 0.00],
#[9.8500, 44.0667, 0.00],
#[10.1167, 33.8833, 0.00],
[10.3000, 43.5333, 0.00],
#[10.7667, 34.7333, 0.00],
#[11.1167, 33.5000, 0.00],
#[11.7833, 42.1000, 0.00],
[12.0000, 36.7833, 0.00],
#[12.3333, 45.4167, 0.00],
[12.5000, 35.5000, 0.00],
[12.5833, 37.6333, 0.00],
[12.8167, 36.1667, 0.00],
[13.2000, 32.9000, 0.00],
#[13.3333, 38.1333, 0.00],
#[13.5000, 43.6167, 0.00],
[13.5000, 37.2500, 0.00],
[13.7500, 45.6500, 0.00],
[13.9333, 40.7333, 0.00],
#[14.2667, 40.8667, 0.00],
[14.4000, 42.5167, 0.00],
[14.5167, 35.9000, 0.00],
#[14.5333, 45.3000, 0.00],
[14.9667, 38.4833, 0.00],
#[15.1000, 38.5000, 0.00],
[15.1500, 36.6833, 0.00],
#[15.2500, 38.2167, 0.00],
#[15.6500, 38.1167, 0.00],
[15.7667, 43.0333, 0.00],
[16.1833, 41.8833, 0.00],
#[16.4333, 43.5000, 0.00],
#[17.2167, 40.4667, 0.00],
#[17.9333, 40.6333, 0.00],
#[18.5000, 40.1500, 0.00],
[19.1000, 42.0667, 0.00],
#[20.7000, 38.8333, 0.00],
[21.3167, 37.6333, 0.00],
#[22.1333, 37.0167, 0.00],
#[23.0333, 40.6167, 0.00],
[23.8000, 32.1833, 0.00],
#[24.0500, 35.5000, 0.00],
[24.9167, 37.4333, 0.00],
#[25.1333, 35.3333, 0.00],
[25.3833, 40.8500, 0.00],
[25.7000, 31.7667, 0.00],
[26.1500, 38.3833, 0.00],
[26.8833, 37.0833, 0.00],
[28.2333, 36.4333, 0.00],
[29.8667, 31.2000, 0.00],
[32.3167, 31.2667, 0.00],
#[33.3167, 35.3333, 0.00],
#[33.9500, 35.1167, 0.00]
])

M2amp = u.ProbeData(pts, "M2amp")
(ilen, jlen) = M2amp.shape
S2amp = u.ProbeData(pts, "S2amp")
K1amp = u.ProbeData(pts, "K1amp")
O1amp = u.ProbeData(pts, "O1amp")


for i in range(ilen):
    g.write("%f\n" % M2amp[i][0])
ampcm=M2amp*100

M2_tideGauge_amp =  array([
[29.8],
[18.0],
[18.0],
[9.0],
[2.0],
[3.0],
[4.6],
[7.0],
[5.6],
[8.3],
[6.5],
[8.6],
[7.6],
[9.4],
[51.1],
[8.5],
[41.6],
[21.9],
[10.9],
[1.6],
[23.4],
[6.6],
[4.3],
[4.8],
[11.1],
[10.6],
[6.6],
[4.5],
[26.3],
[12.0],
[11.1],
[6.4],
[6.0],
[10.6],
[12.0],
#[6.4],
[6.7],
[12.0],
[6.2],
[6.8],
[7.9],
[8.0],
[6.5],
[8.7],
[7.0],
[9.2],
[4.0],
[3.3],
[2.2],
[9.0],
[1.4],
[1.0],
[2.0],
[1.5],
[7.1],
[2.9],
[4.4],
[2.1],
[4.4],
[7.2],
[11.2],
[10.1],
[11.0]
])
from math import sqrt
ampdiff=ampcm-M2_tideGauge_amp
ampdiff2=ampdiff**2
a = sum(sum(ampdiff2))/62
RMS=sqrt(a)
print("RMS difference of M2 Amp (cm):",RMS)

S2_tide_guage_data_amp =  array([
[10.7],
[7],
[7],
[4],
[1],
[1],
[1.8],
[2],
[2.2],
[3.4],
[2.6],
[3.2],
[2.8],
[3.4],
[36.4],
[3.4],
[26.7],
[15.3],
[4.1],
[1.9],
[14.1],
[4.2],
[1.8],
[3.1],
[5.4],
[4.1],
[3.6],
[3.2],
[15.2],
[5],
[4.4],
[4.5],
[4],
[5.5],
[4.5],
#[3.4],
[3.5],
[4.7],
[3.1],
[4.4],
[5.1],
[5.6],
[3.7],
[5.2],
[4],
[5.6],
[2.2],
[1.6],
[1.1],
[6.1],
[1.3],
[0.8],
[1],
[1.1],
[5],
[2.9],
[2.9],
[1.3],
[2.7],
[4.1],
[6.9],
[6.4],
[7.3]
])

S2ampcm=S2amp*100
S2ampdiff=S2ampcm-S2_tide_guage_data_amp
S2ampdiff2=S2ampdiff**2
S2a = sum(sum(S2ampdiff2))/62
S2RMS=sqrt(S2a)
print("RMS difference of S2 Amp (cm):",S2RMS)

K1_tide_guage_data_amp=  array([
[2],
[3],
[4],
[3],
[4],
[4],
[3.2],
[3],
[2.3],
[3.6],
[3.2],
[3.3],
[3.2],
[3.7],
[2.5],
[4],
[1.8],
[2],
[2.7],
[2],
[17.9],
[0.9],
[3.5],
[0.5],
[2],
[3.2],
[13],
[1.8],
[19.7],
[3],
[2.8],
[9.7],
[1],
[13.8],
[3.1],
#[1.5],
[1.9],
[3.3],
[1.3],
[6.8],
[4.2],
[8.8],
[1.8],
[4.6],
[2.3],
[4.8],
[1.4],
[1.3],
[1.2],
[2.6],
[0.6],
[1.4],
[1.9],
[1.8],
[0.3],
[1.2],
[2.3],
[2],
[1.8],
[1.7],
[2.1],
[2.4],
[2.1]
])

K1ampcm=K1amp*100
K1ampdiff=K1ampcm-K1_tide_guage_data_amp
K1ampdiff2=K1ampdiff**2
K1a = sum(sum(K1ampdiff2))/62
K1RMS=sqrt(K1a)
print("RMS difference of K1 Amp (cm):",K1RMS)

O1_tide_guage_data_amp=  array([
[0.9],
[2],
[1],
[2],
[2],
[2],
[1.9],
[2],
[2],
[1.6],
[1.9],
[1.4],
[1.8],
[1.4],
[0.5],
[1.8],
[0.8],
[0.9],
[1.2],
[1.4],
[5.6],
[0.7],
[1.6],
[0.9],
[0.6],
[1.2],
[4.2],
[1.4],
[6.1],
[1],
[1],
[3.4],
[1],
[4.1],
[1.1],
#[1.1],
[0.9],
[1.1],
[0.9],
[2.5],
[1.5],
[2.7],
[0.8],
[1.5],
[1],
[1.4],
[0.6],
[0.5],
[0.5],
[1.3],
[0.5],
[0.6],
[1],
[0.9],
[1.3],
[0.8],
[1.3],
[1.1],
[1.1],
[1.3],
[1.7],
[1.8],
[1.8]
])

O1ampcm=O1amp*100
O1ampdiff=O1ampcm-O1_tide_guage_data_amp
O1ampdiff2=O1ampdiff**2
O1a = sum(sum(O1ampdiff2))/62
O1RMS=sqrt(O1a)
print("RMS difference of O1 Amp (cm):",O1RMS)

import fluidity_tools
from matplotlib import pylab
pylab.plot(ampcm,M2_tideGauge_amp)
pylab.xlabel("Fluidity")
pylab.ylabel("Tide Gauge")
pylab.show()

import matplotlib
matplotlib.pyplot.scatter(M2_tideGauge_amp,ampcm,s=20, c='b', marker='o')
#pylab.xlabel("Tide Gauge M2 Amplitude (cm)")
#pylab.ylabel("Fluidity M2 Amplitude (cm)")
#x=([0,70])
#y=([0,70])
#matplotlib.pyplot.plot(y,x, label="y=x")
#pylab.ylim(ymax=70.0,ymin=0.0)
#pylab.xlim(xmax=70.0,xmin=0.0)
#pylab.show()


#matplotlib.pyplot.scatter(S2_tide_guage_data_amp,S2ampcm,s=20, c='b', marker='o')
#pylab.xlabel("Tide Gauge S2 Amplitude (cm)")
#pylab.ylabel("Fluidity S2 Amplitude (cm)")
#x=([0,70])
#y=([0,70])
#matplotlib.pyplot.plot(y,x, label="y=x")
#pylab.ylim(ymax=70.0,ymin=0.0)
#pylab.xlim(xmax=70.0,xmin=0.0)
#pylab.show()

#matplotlib.pyplot.scatter(K1_tide_guage_data_amp,K1ampcm,s=20, c='b', marker='o')
#pylab.xlabel("Tide Gauge K1 Amplitude (cm)")
#pylab.ylabel("Fluidity K1 Amplitude (cm)")
#x=([0,20])
#y=([0,20])
#matplotlib.pyplot.plot(y,x, label="y=x")
#pylab.ylim(ymax=20.0,ymin=0.0)
#pylab.xlim(xmax=20.0,xmin=0.0)
#pylab.show()

#matplotlib.pyplot.scatter(O1_tide_guage_data_amp,O1ampcm,s=20, c='b', marker='o')
#pylab.xlabel("Tide Gauge O1 Amplitude (cm)")
#pylab.ylabel("Fluidity O1 Amplitude (cm)")
#x=([0,20])
#y=([0,20])
#matplotlib.pyplot.plot(y,x, label="y=x")
#pylab.ylim(ymax=20.0,ymin=0.0)
#pylab.xlim(xmax=20.0,xmin=0.0)
#pylab.show()




