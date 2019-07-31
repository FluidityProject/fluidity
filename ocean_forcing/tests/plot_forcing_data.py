#!/usr/bin/env python3

import numpy
from pylab import *
import string

####################
#                  #
#      CONFIG      #
#                  #
####################
# change default font sizes
params = {
          'legend.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'axes.labelsize' : 14,
          'text.fontsize' : 14,
}
rcParams.update(params)

filelist = sys.argv[1:]
files = []
for file in filelist:
   # Just to make sure we grab the test data only
   if (string.find(file, 'test_data_') != -1):
     files.append(file)


for file in files:
    print("Plotting graphs for ",file)

    # make a graphic filename by stripping off the csv and adding png

    f = open(file,'r')
    nLines = 0
    q_flux = []
    q_era40 = []
    f_flux = []
    f_era40 = []
    taux_flux = []
    taux_era40 = []
    tauy_flux = []
    tauy_era40 = []
    timestep = []
    for line in f:
        line = line.strip()
        if (nLines == 0):
            Labels = line.split(",")
        else:
            temp = line.split(",")
            timestep.append(temp[0])
            q_era40.append(temp[1])
            q_flux.append(temp[2])
            taux_era40.append(temp[3])
            taux_flux.append(temp[4])
            tauy_era40.append(temp[5])
            tauy_flux.append(temp[6])
            f_era40.append(temp[7])
            f_flux.append(temp[8])
        nLines =+ 1

    f.close()

    fig = figure()
    ax = fig.add_subplot(221)
    ax = fig.add_subplot(221)
    ax.plot(timestep,q_era40,'b',label="ERA40 Q")
    ax.plot(timestep,q_flux,'r',label="Our Q")
    ax.grid(True)
    ax.set_ylabel("Heat Flux")
    ax.set_xlabel("Time")
    ax.legend(loc=0)
    ax = fig.add_subplot(222)
    ax.plot(timestep,f_era40,'b',label="ERA40 F")
    ax.plot(timestep,f_flux,'r',label="Our F")
    ax.grid(True)
    ax.set_ylabel("Freshwater Flux")
    ax.set_xlabel("Time")
    ax.legend(loc=0)
    ax = fig.add_subplot(223)
    ax.plot(timestep,taux_era40,'b',label="ERA40 Tau X")
    ax.plot(timestep,taux_flux,'r',label="Our Tau X")
    ax.grid(True)
    ax.set_ylabel("TauX")
    ax.set_xlabel("Time")
    ax.legend(loc=0)
    ax = fig.add_subplot(224)
    ax.plot(timestep,tauy_era40,'b',label="ERA40 Tau Y")
    ax.plot(timestep,tauy_flux,'r',label="Our Tau Y")
    ax.grid(True)
    ax.set_ylabel("TauY")
    ax.set_xlabel("Time")
    ax.legend(loc=0)

    show()

