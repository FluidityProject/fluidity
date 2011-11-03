#!/usr/bin/env python

import sys

print '\n\n    working...\n\n'

# read in file name
filename = sys.argv[1]

# open file
f = open(filename, 'r')

u=[];v=[];w=[];vel=[];Y=[]

# read in the velocity data file x.nnn
# ignore header
ignore=3
for i in range (0+ignore):
   f.readline()

print ' 2. reading U_rms data from file x.181'

lines=0
for line in f:
    lines+=1
    data=line.split()
    Y.append( float( data[0] ) )
    vel.append( float( data[1] ) )
    u.append( float( data[3] ) )
    v.append( float( data[4] ) )
    w.append( float( data[5] ) )
f.close()

# Specify how much to reduce data set by
d = 4
from math import trunc
short = trunc( len(Y)/d )
print "entries: ", short

# Get maximum velocity
velmax = max(vel)
print "k normalised by U_max**2.; U_max = ", velmax

print ' 3. write new file'

nf = open('stepdata2d-k181.dat','w')
yy=[];kk=[]

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
nf.write('K/U_max^2''\n')
for i in range( short ):
    # bottom half on inlet
    kk.append( 0.5 * (u[i*d]**2.+v[i*d]**2.+w[i*d]**2.) )
    nf.write(str('%.12f' % kk[i] )+', ' )
for i in range( short ):
    # top half on inlet
    nf.write(str('%.12f' % kk[short-i-1] )+', ' )
nf.write('\n\n''U_max''\n')
nf.write(str('%.12f' % velmax ) )

print '\n\n               ...all done\n\n '

sys.exit()
