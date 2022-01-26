#!/usr/bin/env python3

## ./Stommel_tracer_advection_error.py Stommel2_adapt-a_20.vtu 10.e6 

from numpy import sqrt,pi,exp,sin,array,size,abs,max,amax,dot,zeros,reshape,transpose,arange
from pylab import *
import sys
sys.path.append('/home/piggott/bin/')
import vtktools
import nodecount
import scipy.integrate


def IC(X):
  from math import exp
  L = 1.e6
  T =  10.0*exp(-((X[0] - 2.*L/3.)**2 + (X[1] - L/3.)**2)/(2.*(L/12.5)**2))
  return T  


def func(X, t):
  from math import exp, sin, cos
  ##L=1.e6; F=0.1; rho=1000.0; H=200.0; gamma=1.e-6;
  ##A = F*L/(pi*gamma*rho*H) 
  A = 159154.94309189534
  ##alpha = beta/gamma
  ##1/alpha = 100km=1e5; alpha = 1.e-5
  ##beta=1e-11
  ##zplus = -alpha/2 + sqrt((alpha**2)/4 + (pi/L)**2)
  zplus = 9.0504906000698326e-07
  ##zminus = -alpha/2 - sqrt((alpha**2)/4 + (pi/L)**2)
  zminus = -1.0905049060006982e-05  
  ##p=(1-exp(L*zminus))/(exp(L*zplus)-exp(L*zminus))  
  p =  0.40451761482981791
  ## q = 1-p
  q = 0.59548238517018204
  ## pi/L = 3.1415926535897933e-06
  ##psi = A*sin(3.1415926535897933e-06*y)*(p*exp(x*zplus) + q*exp(x*zminus) - 1.0)
  ## u = d psi/d y
  ## A*pi/L = 0.5
  u = 0.5*cos(3.1415926535897933e-06*X[1])*(p*exp(X[0]*zplus) + q*exp(X[0]*zminus) - 1.0)
  ## v = -d psi/d x
  v = -A*sin(3.1415926535897933e-06*X[1])*(p*zplus*exp(X[0]*zplus) + q*zminus*exp(X[0]*zminus))
  return [-u,-v]  ## as we're going backwards


def function_analytical(x,y,t_start,t_end):
  from math import exp
  L = 1.e6
  pos = scipy.integrate.odeint(func, [x, y], [t_start, t_end])
#  print("Departure point:",pos[-1,:])
  Analytical = 10.0*exp(-((pos[-1,0] - 2.*L/3.)**2 + (pos[-1,1] - L/3.)**2)/(2.*(L/12.5)**2))
#  print("Analytical solution:",Analytical)
  return Analytical


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
ug.GetFieldNames()

temp=ug.GetScalarField('Temperature')

t_start = 0.0
t_end = float(sys.argv[2])

pos=ug.GetLocations()
x=pos[:,0]; y=pos[:,1]; z=pos[:,2]

psi=zeros(size(x), float)
for i in range(len(x)):
#  print("Arrival point,i:",x[i],y[i],i)
  psi[i] = function_analytical(x[i],y[i],t_start,t_end)


NE=ug.ugrid.GetNumberOfCells()
ML=zeros(size(x), float)
for ele in range(NE):
	ndglno=ug.GetCellPoints(ele)
	trivol=function_trivol(x[ndglno],y[ndglno])
        for nod in ndglno:
                ML[nod] = ML[nod] + trivol/3


err = abs(psi - temp)
norm1 = dot(ML,err)
norm2 = dot(ML,err**2)
snorm2 = sqrt(abs(norm2))
normaliseL2 = sqrt(abs(dot(ML,psi**2)))
##print(norm1)
##print(snorm2)
##print(amax(err))
print(n,snorm2/normaliseL2)
ug.AddScalarField('Analytical Solution', psi)
ug.AddScalarField('Error (difference)', err)
ug.Write('error.vtu')


