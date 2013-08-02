#!/usr/bin/python

import fluidity.diagnostics.fluiditytools as fluidity_tools
import numpy
import pylab

vt = fluidity_tools.stat_parser('../darcy_impes_MIM_1D_test.stat')
cd = vt['Phase2']['C_d']['integral']
cs = vt['Phase2']['C_s']['integral']
prt=vt['Phase1']['Porosity'][u'max']
msat=vt['Phase2']['MobileSaturation'][u'max']
isat=vt['Phase2']['ImmobileSaturation'][u'max']
theta_s=prt*isat
theta_d=prt*msat
c = theta_d*cd + theta_s*cs
c[0] = c[1]
rtd = 1-c/c[0] 
t=vt['ElapsedTime']
t=t[u'value']


#extract the data from excel
f=file('RTD_ex.csv','r')
RTD2=[]
for row in f:
	RTD2.append(float(row))

ft=file('RTD_tim.csv','r')
t2=[]
for row2 in ft:
	t2.append(float(row2))



#start plotting
pylab.plot(t,rtd,color='red',lw=2)
pylab.plot(t2,RTD2,color='b',marker='.', linestyle='None')
pylab.xlabel('time (s)')
pylab.ylabel('Residence time distribution')
pylab.legend(['Mobile-immobile model', 'Experiment data'], loc='lower right')
pylab.savefig('plotmim.png')

