from numpy import *

fin = open('channel.node','r')
fout = open('equil.node','w')

s = fin.readline()
fout.write(s)

from math import sqrt

exiting = False
while(exiting==False):
    s = fin.readline()
    if(s[0]=='#'):
        fout.write(s)
        exiting=True
    else:
        vals = s.split()
        assert(len(vals)==4)
        print vals
        x = float(vals[1])
        y = float(vals[2])
        X = x - 0.5*y
        Y = sqrt(3./4.)*y
        fout.write(vals[0]+' '+str(X)+' '+str(Y)+' '+'0\n')
fin.close()
fout.close()
    
