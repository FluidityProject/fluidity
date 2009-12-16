import os
import sys
import string
from math import *
import scipy

polyFile=open(sys.argv[1] + ".poly")
data = polyFile.readline().strip().split()
npnodes=int(data[0])
print npnodes

nodeXDictionary={}
nodeYDictionary={}
polyN1Dictionary={}
polyN2Dictionary={}
polyBndDictionary={}
for n in range(npnodes): 
        data = string.splitfields(polyFile.readline())
        nodeXDictionary[int(data[0])] = string.atof(data[1])
        nodeYDictionary[int(data[0])] = string.atof(data[2]) 
data = string.splitfields(polyFile.readline().strip())
npoly=int(data[0])
nbnd=int(data[1])
maxbnd=0



for n in range(npoly):
       line=polyFile.readline()
       if not line:
              continue
       line = line.strip()
       if len(line) == 0:
              continue
       data = string.splitfields(line)
       if (data[0][0] in string.digits):
              polyN1Dictionary[int(data[0])] = int(data[1])
              polyN2Dictionary[int(data[0])] = int(data[2]) 
              if (nbnd==1):
                     polyBndDictionary[int(data[0])] = int(data[3])
              else:
                     polyBndDictionary[int(data[0])] = 1
data = string.splitfields(polyFile.readline().strip())

holeXDictionary={}
holeYDictionary={}
nholes=int(data[0])
for n in range(nholes):
       line=polyFile.readline()
       if not line:
              continue
       line = line.strip()
       if len(line) == 0:
              continue       
       data = string.splitfields(line)
       holeXDictionary[int(data[0])]=string.atof(data[1])
       holeYDictionary[int(data[0])]=string.atof(data[2])     
polyFile.close()

polyOutDictionary={}
ncp=1
ncn=npnodes+1
for n in range(1,npoly+1):
    n_1=polyN1Dictionary[n]
    n_2=polyN2Dictionary[n]
    node_1=[nodeXDictionary[n_1],nodeYDictionary[n_1]]
    node_2=[nodeXDictionary[n_2],nodeYDictionary[n_2]]
    bnd=polyBndDictionary[n]
    distance=pow(node_1[0]-node_2[0],2)+pow(node_1[1]-node_2[1],2)
    m=int(ceil(sqrt(distance)/scipy.sqrt(float(sys.argv[2]))))
    vec=[node_2[0]-node_1[0],node_2[1]-node_1[1]]
    polyOutDictionary[ncp]=[n_1,ncn,bnd]
    ncp+=1
    print m
    for p in range(1,int(m)-1):
        polyOutDictionary[ncp]=[ncn,ncn+1,bnd]
        nodeXDictionary[ncn]=node_1[0]+p*vec[0]/m
        nodeYDictionary[ncn]=node_1[1]+p*vec[1]/m
        ncp+=1
        ncn+=1
    nodeXDictionary[ncn]=node_1[0]+(m-1)*vec[0]/m
    nodeYDictionary[ncn]=node_1[1]+(m-1)*vec[1]/m
    polyOutDictionary[ncp]=[ncn,n_2,bnd]
    ncn+=1
    ncp+=1

polyFile=open(sys.argv[1]+".1.poly",'w')
outline= "%d %d %d %d\n" % (ncn-1,2,0,1)
polyFile.write(outline)
for n in range(1,ncn):
    outline="%d %f %f %d\n" % (n,nodeXDictionary[n],nodeYDictionary[n],0)
    polyFile.write(outline)
outline= "%d %d\n" % (ncp-1,1)
polyFile.write(outline)
for n in range(1,ncp):
    outline="%d %d %d %d\n"% (n,polyOutDictionary[n][0],\
                                  polyOutDictionary[n][1],\
                                  polyOutDictionary[n][2])
    polyFile.write(outline)
outline= "%d\n" % (nholes)
polyFile.write(outline)
for n in range(1,nholes+1):
    outline= "$%d %f %f\n" % (n,holeXDictionary[n],holeYDictionary[n])
    polyFile.write(outline)
outline= "%d\n" % (0)
polyFile.write(outline)
polyFile.close


    
