#!/usr/bin/env python

import math
import numpy

# Read in node file

f=open('shoaling_tsunami.node','r')
data = [line.split() for line in f]

# Close initial file
f.close()


i=1
for i in range(1,len(data)-1):
        x = float(data[i][1])
        z_o = float(data[i][3])
        if(x<0):
             z_new = z_o + (z_o/(-4000.))*3800.0*(0.5*(1.0+math.tanh((-x-1000000)/100000.0)))
        else:
             z_new = z_o + (z_o/(-4000.))*3500.0*(0.5*(1.0+math.tanh((x-1000000)/100000.0)))


        data[i][3]=str(z_new)
        i=i+1

# write new node coordinates to new node file
f=open('output.node','w')
i=0
for i in range(0,len(data)):
    f.write(data[i][0]+" ")
    f.write(data[i][1]+" ")
    f.write(data[i][2]+" ")
    f.write(data[i][3]+" \n")
    i=i+1
    
f.close()


