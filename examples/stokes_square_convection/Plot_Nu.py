#!/usr/bin/python3

# This script plots up the Nusselt number for both the 24x24 and 48x48 cases:
import pylab
from fluidity_tools import stat_parser as stat

# Stafiles:
statfile24="stokes-sc-Ra1e5-24.stat"
statfile48="stokes-sc-Ra1e5-48.stat"

# First plot 24x24 case:
pylab.plot(stat(statfile24)["CoordinateMesh"]["nodes"][-1],
      -stat(statfile24)["Fluid"]["Temperature"]["surface_integral%Top"][-1],
      linestyle='None', marker='o', markerfacecolor='0.15')

# Next plot 48x48 case:
pylab.plot(stat(statfile48)["CoordinateMesh"]["nodes"][-1],
      -stat(statfile48)["Fluid"]["Temperature"]["surface_integral%Top"][-1],
       linestyle='None', marker='o', markerfacecolor='0.15')

# Plot benchmark value as line for comparison:
pylab.plot([100,8e4],[10.534,10.534],'k--',lw=0.6)

pylab.xlabel(r"Vertices")
pylab.ylabel(r"Nusselt Number")
pylab.xlim(100,1e4)
pylab.ylim(9.0,11.0)

pylab.savefig("Nu_1e5.png")
