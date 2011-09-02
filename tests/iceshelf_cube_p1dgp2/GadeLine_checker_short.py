#!/usr/bin/env python

import vtktools
from numpy import mean, sum
import pylab 
def Gade_line(Sg1,T1,S1):
    c0=3974.0
    cI=2009.0
    L=3.35e5
    Li=L
    TI=-25
    a=-0.0573
    b=0.0832
    Cd = 1.5e-3
    Tf=a*S1+b
    dTdS1 = ((T1-Tf)+Li/c0+(cI/c0)*(Tf-TI)) /(S1-0)
    Tg1 = dTdS1*(Sg1-S1) +T1
    return Tg1
    
fname_a = "shelf3d_1.vtu"
datafile = vtktools.vtu(fname_a)
T=datafile.GetScalarField("Temperature")
S=datafile.GetScalarField("Salinity")

###Construc Gade line, need the first time
T1=0.0
S1=35.0

S_gade=S
T_gade=Gade_line(S_gade,T1,S1)

#Compare T_gade and T_1
#Calculate R2 = 
T_bar = mean(T)
SS_tot = sum((T-T_bar)**2)
SS_err = sum((T_gade-T)**2)
R2 = 1-(SS_err)/(SS_tot)

pylab.figure
pylab.plot(S,T,'.')
pylab.show()


