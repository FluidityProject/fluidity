#!/usr/bin/env python

import sys

print '\n\n    working...\n\n'

# read in file name. Should be rskk.nnn
filename = sys.argv[1]

# open file
f = open(filename, 'r')

e=[];Y=[]

# read in the kinetic energy equation terms file
# ignore header
ignore=8
for i in range (0+ignore):
   f.readline()

print ' 2. reading viscous dissipation data from file rskk.181'

lines=0
for line in f:
    lines+=1
    data=line.split()
    Y.append( float( data[0] ) )
    # dissipation term is negative in k equation, we want positive
    e.append( -float( data[5] ) )
f.close()

# Specify how much to reduce data set by
d = 4
from math import trunc
short = trunc( len(Y)/d )
print "entries: ", short

# Get maximum velocity from stepdata3d-k.py
velmax = 1.00647000; h = 1.
print "epsilon normalised by U_max**3./h; U_max = U_max from stepdata3d-k.py: ", velmax
print "step height h: ", h

print ' 3. write new file'

nf = open('stepdata3d-eps181.dat','w')
yy=[];ee=[]

nf.write('ERCOFTAC database: Le&Moin et al, 1997 DNS: 3d backward facing step, Re=5100''\n')
nf.write('Y/h''\n')
for i in range( short ):
    # Scale to fit bottom half on inlet
    yy.append( Y[i*d]/10. )
    nf.write(str('%.12f' % yy[i] )+', ' )
for i in range( short ):
    # Scale to fit top half on inlet
    nf.write(str('%.12f' % (1.0-yy[short-i-1]) )+', ' )

nf.write('\n\n')
nf.write('eps/(U_max^3 * h)''\n')
for i in range( short ):
    # bottom half on inlet
    ee.append( e[i*d] )
    nf.write(str('%.12f' % ee[i] )+', ' )
for i in range( short ):
    # top half on inlet
    nf.write(str('%.12f' % ee[short-i-1] )+', ' )

nf.write('\n\n''U_max''\n')
nf.write(str('%.12f' % velmax ) )
nf.write('\n\n''h''\n')
nf.write(str( h ) )

print '\n\n               ...all done\n\n '
sys.exit()
