import math
import numpy
#import numarray.strings
#import numarray.objects

# Read in node file

f=open('incline_2-5.node','r')
data = [line.split() for line in f]

# Close initial file
f.close()

#convert data to floats
#data = data.astype('Float64')
#data = numpy.array(data)
#Do all stuff to get bottom bdry in
thet = 2.5
xincline = 19000 - (800/(math.tan(thet*math.pi/180)))
xplateau = 19000

i=1
for i in range(1,len(data)-1):
        x = float(data[i][1])
        z_o = float(data[i][2])
        if (xincline < x  < xplateau) :
            z_new = z_o + ((x-xincline)*(1.0-(z_o/1000.0))*math.tan(thet*math.pi/180.0))

        elif (x >= xplateau) :
            z_new = 800.0 + ((1000.0-800.0)/1000.0)*z_o

        else :
            z_new = z_o

        data[i][2]=str(z_new)
        i=i+1

# write new node coordinates to new node file
f=open('incline_2-5.node','w')
f.write(data[0][0]+" ")
f.write(data[0][1]+" ")
f.write(data[0][2]+" ")
f.write(data[0][3]+" \n")
for i in range(1,len(data)):
    f.write(data[i][0]+" ")
    f.write(data[i][1]+" ")
    f.write(data[i][2]+" \n")
    i=i+1
    
f.close()


