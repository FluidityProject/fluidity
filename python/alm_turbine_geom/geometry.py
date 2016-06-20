#!/usr/bin/env python

import os, glob
import shutil
import sys
import string
import filecmp
import difflib
import subprocess
import math
import re
import argparse
import time
from numpy import genfromtxt
import matplotlib.pyplot as plt
from math import *

#======================================
# Classes
#=====================================

class Blade:
    def __init__(self,nel,rtoR,ctoR,pitch,thtoC):
        self.Nelem=nel;
        self.r=rtoR;
        self.c=ctoR;
        self.pitch=pitch;
        self.thtoC=thtoC;

    def create(self,radius):
        for i in range(0, 1,self.Nelem):
            self.r[i]=self.r[i]*radius;
            self.c[i]=self.c[i]*radius;
            self.pitch[i]=self.pitch[i];
            self.thtoC[i]=self.thtoC[i];

    def rotate(self,Theta,nR,Origin):

        print('Rotating blade');

    def QuatRot(v,Theta,nR,Origin):
        # Perform rotation of vector(s) v around normal vector nR using the
        # quaternion machinery

        # Force normalize nR
        nR=nR/sqrt(sum(pow(nR,2)));

        # Quaternion form of v

class HATT:
    def __init__(self,nbs):
        self.Nblades=nbs;
        self.Blade(nbs)

# ====================================
# Parser Arguments
# ====================================
parser = argparse.ArgumentParser(description="Creates the turbine geometry (.geom) for Fluidity")
parser.add_argument("-v","--verbose",action="store_true",help="Prints a script description on the screen")
parser.add_argument("-p","--plot",action="store_true",help="Plots the Cp versus the TSR")
parser.add_argument("IFBlades", type=str, help="The input filename containing the blade geometry information")
parser.add_argument("IFTurbine", type=str, help="The input filename containing information about the turbine set-up")

args = parser.parse_args()
blade_file = args.IFBlades
turbine_file =args.IFTurbine

# Define the blade parameters
rtoR= [];
ctoR= [];
Pitch= [];
thtoC= [];

Blade_Par = genfromtxt(blade_file,delimiter=',',skip_header=1)
rtoR.append(Blade_Par[:,0])
ctoR.append(Blade_Par[:,1])
Pitch.append(Blade_Par[:,2])
thtoC.append(Blade_Par[:,3])

NStations = len(Blade_Par)
NElem = NStations-1

# Define blade parameters
Turbine_Par = genfromtxt(turbine_file,comments='#')
Radius = Turbine_Par[0]
HubRadius = Turbine_Par[1]
Tilt = Turbine_Par[2]
bCone = Turbine_Par[3]
Yaw = Turbine_Par[4]
set_angle = Turbine_Par[5]
NBlades = Turbine_Par[6]


B=Blade(NElem,rtoR,ctoR,Pitch,thtoC)
B.create(Radius)

print(B.r)
print(B.c)

# Create the Outputfile
#OFN = 'Bahaj2007_tilt'+str(Tilt)+'_yaw'+str(Yaw)+'_setangle'+str(set_angle)+'.geom'
