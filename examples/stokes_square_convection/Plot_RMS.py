#!/usr/bin/python3

# This sript plots up the RMS velocity for both the 24x24 and 48x48 case
import pylab
from fluidity_tools import stat_parser as stat

# Stafiles:
statfile24="stokes-sc-Ra1e5-24.stat"
statfile48="stokes-sc-Ra1e5-48.stat"

# First plot 24x24 case:
pylab.plot(stat(statfile24)["CoordinateMesh"]["nodes"][-1],
      stat(statfile24)["Fluid"]["Velocity%magnitude"]["l2norm"][-1],
      linestyle='None', marker='o', markerfacecolor='0.15')

# Next plot 48x48 case:
pylab.plot(stat(statfile48)["CoordinateMesh"]["nodes"][-1],
      stat(statfile48)["Fluid"]["Velocity%magnitude"]["l2norm"][-1],
      linestyle='None', marker='o', markerfacecolor='0.15')

# Plot benchmark value as line for comparison:
pylab.plot([100,8e4],[193.214,193.214],'k--',lw=0.6)

pylab.xlabel(r"Vertices")
pylab.ylabel(r"RMS Velocity")
pylab.xlim(100,1e4)
pylab.ylim(192.0,195.0)

pylab.savefig("RMS_1e5.png")

