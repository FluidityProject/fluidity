#!/usr/bin/env python

# Read in node file
f=open('bingham-slope.node','r')
data = [line.split() for line in f]
# Close initial file
f.close()

for i in range(1,len(data)-1):
   y = float(data[i][2])
   x_o = float(data[i][1])
   x_new = x_o - y*(x_o/3.2)*2.0
   data[i][1]=str(x_new)

# write new node coordinates to new node file
f=open('bingham-slope.node','w')

for i in range(0,len(data)):
   f.write(data[i][0]+" ")
   f.write(data[i][1]+" ")
   f.write(data[i][2]+" ")
   f.write(data[i][3]+" \n")
    
f.close()


