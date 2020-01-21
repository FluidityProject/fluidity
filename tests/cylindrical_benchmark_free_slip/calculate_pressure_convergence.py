#!/usr/bin/env python

from fluidity_tools import stat_parser as stat
import numpy as np

# Wavenumber:
n = 2
resolutions = ["A","B"]
expected_pressure_convergence              = 1.5

# Extract data
data                                       = np.zeros(len(resolutions))
pressure_convergence                       = np.zeros(len(resolutions))

# Statfile from finest mesh for relative errors:
statfile_fine="cylindrical_%s.stat" % str(resolutions[-1])

######## Pressure Errors: #######
anal = stat(statfile_fine)['Fields']['AnalyticalPressure']['l2norm'][-1]
for resolution in range(len(resolutions)):
    statfile="cylindrical_%s.stat" % str(resolutions[resolution])
    data[resolution]  = stat(statfile)['Fields']['PressureError']['l2norm'][-1] / anal
    if(resolution > 0):
        pressure_convergence[resolution] = np.log2(data[resolution-1] / data[resolution])
        print (pressure_convergence[resolution])
