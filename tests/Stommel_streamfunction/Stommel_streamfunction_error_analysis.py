#!/usr/bin/env python

from __future__ import print_function

from numpy import sqrt,pi,exp,sin,array,size,abs,max,amax,dot,zeros,reshape,transpose,arange
from pylab import *
import sys
sys.path.append('/home/piggott/bin/')
import vtktools
import nodecount

def function_psi(X,Y):
   lam = 1.0;      b = 1.0;
   r = 1.0;        rho = 1.0;
   beta = 100.0;   alpha = beta/r;
   D = 200;        F = 0.1;
   gamma = F*pi/(rho*r*b);
   A = -alpha/2 + sqrt(0.25*(alpha**2) + (pi/b)**2);
   B = -alpha/2 - sqrt(0.25*(alpha**2) + (pi/b)**2);
   p = (1-exp(B*lam))/(exp(A*lam)-exp(B*lam));
   q = 1-p;
   return (1/rho)*gamma*((b/pi)**2)*sin(pi*Y/b)*(p*exp(A*X)+q*exp(B*X)-1)



def function_tetvol(X,Y,Z):
   X12 = X[1] - X[0]; X13 = X[2] - X[0]; X14 = X[3] - X[0]
   Y12 = Y[1] - Y[0]; Y13 = Y[2] - Y[0]; Y14 = Y[3] - Y[0]
   Z12 = Z[1] - Z[0]; Z13 = Z[2] - Z[0]; Z14 = Z[3] - Z[0]   
   VOL = X12*( Y13*Z14 - Y14*Z13 ) + X13*( Y14*Z12 - Y12*Z14 ) + X14*( Y12*Z13 - Y13*Z12 )
   return abs(VOL/6)


def function_trivol(X,Y):
   a = sqrt( (X[1] - X[0])**2 + (Y[1] - Y[0])**2) 
   b = sqrt( (X[2] - X[0])**2 + (Y[2] - Y[0])**2)
   c = sqrt( (X[2] - X[1])**2 + (Y[2] - Y[1])**2)
   s = 0.5*(a+b+c)
   VOL = sqrt(s*(s-a)*(s-b)*(s-c))
   return abs(VOL/2)



n = nodecount.nodecount(sys.argv[1])
print(n)


ug=vtktools.vtu(sys.argv[1])
##ug=vtktools.vtu('Stommel-NEW5_3e-6----_6.vtu')
ug.GetFieldNames()
##p=ug.GetScalarField('Pressure')
##uvw=ug.GetVectorField('Velocity')
temp=ug.GetScalarField('Temperature')

psi=ug.GetScalarField('Analytical')

pos=ug.GetLocations()
x=pos[:,0]; y=pos[:,1]; z=pos[:,2]



NE=ug.ugrid.GetNumberOfCells()
ML=zeros(size(x), float)
for ele in range(NE):
   ndglno=ug.GetCellPoints(ele)
   trivol=function_trivol(x[ndglno],y[ndglno])
   for nod in ndglno:
      ML[nod] = ML[nod] + trivol/3


#psi = function_psi(x,y)
err = abs(psi - temp)
norm1 = dot(ML,err)
norm2 = dot(ML,err**2)
snorm2 = sqrt(abs(norm2))
##print norm1
##print snorm2
##print amax(err)
print(n,norm1, snorm2, amax(err))
ug.AddScalarField('Analytical Solution', psi)
ug.AddScalarField('Error (difference)', err)
ug.Write('error.vtu')


N=1000
pts = zeros(N*N*3,float)
for i in range(N):
   for j in range(N):
      pts[i*N + j ]        = float(i)/float(N-1)
      pts[i*N + j + N*N]   = float(j)/float(N-1)
      pts[i*N + j + 2*N*N] = 0.0



xx = pts[0:N*N]
yy = pts[N*N:2*N*N]
zz = pts[2*N*N:3*N*N]
pts2 = reshape(pts,(3,N*N))
pts3 = transpose(pts2)
tempfinemesh = ug.ProbeData(pts3, 'Temperature')
tempfinemesh = transpose(tempfinemesh)
psifinemesh = function_psi(xx,yy)
errfinemesh = abs(psifinemesh - tempfinemesh)
normfinemesh = amax(amax(errfinemesh))
norm2finemesh = sqrt(sum(sum(errfinemesh**2)))
print(normfinemesh,norm2finemesh)

#tt = transpose(tempfinemesh)
#ee = transpose(errfinemesh)

#f = open("tmp.dat","w")
#for i in range(N*N):
#   print >> f,"%12.8f %12.8f %12.8f"%(xx[i],yy[i],ee[i])
#
#f.close()




#tempfinemesh_x = zeros((N-1)*(N-1),float)
#tempfinemesh_y = zeros((N-1)*(N-1),float)
#for i in arange(2,N-2):
#   for j in arange(2,N-2):
#      tempfinemesh_x[(i-1)*N + j-1 ] = (tt[(i+1)*N + j ] - tt[(i-1)*N + j ])/float(N)
#        tempfinemesh_y[(i-1)*N + j-1 ] = (tt[i*N + j+1 ]   - tt[i*N + j-1   ])/float(N)
        
          
            
#contour(reshape(tempfinemesh_x,(N-1,N-1)))
#show()
