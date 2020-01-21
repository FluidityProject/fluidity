#!/usr/bin/env python

from fluidity_tools import stat_parser as stat
import numpy as np

# Wavenumber:
n = 2
resolutions = ["A","B"]
expected_velocity_convergence              = 1.5
expected_surface_velocity_convergence      = 3.0
expected_cmb_velocity_convergence          = 3.0

# Extract data
data                                       = np.zeros(len(resolutions))
velocity_convergence                       = np.zeros(len(resolutions))
surface_velocity_convergence               = np.zeros(len(resolutions))
cmb_velocity_convergence                   = np.zeros(len(resolutions))


# Statfile from finest mesh for relative errors:
statfile_fine="cylindrical_%s.stat" % str(resolutions[-1])

######## Velocity Errors: #######
anal = stat(statfile_fine)['Fields']['AnalyticalVelocity%magnitude']['l2norm'][-1]
for resolution in range(len(resolutions)):
    statfile="cylindrical_%s.stat" % str(resolutions[resolution])
    data[resolution]  = stat(statfile)['Fields']['VelocityError%magnitude']['l2norm'][-1] / anal
    if(resolution > 0):
        velocity_convergence[resolution] = np.log2(data[resolution-1] / data[resolution])
        print (velocity_convergence[resolution])


######## Surface Velocity Errors: #######
anal = stat(statfile_fine)['Fields']['AnalyticalVelocity']['surface_l2norm%Top'][-1]
for resolution in range(len(resolutions)):
    statfile="cylindrical_%s.stat" % str(resolutions[resolution])
    data[resolution]  = stat(statfile)['Fields']['VelocityError']['surface_l2norm%Top'][-1] / anal
    if(resolution > 0):
        surface_velocity_convergence[resolution] = np.log2(data[resolution-1] / data[resolution])
        print (surface_velocity_convergence[resolution])

######## CMB Velocity Errors: #######
anal = stat(statfile_fine)['Fields']['AnalyticalVelocity']['surface_l2norm%Bottom'][-1]
for resolution in range(len(resolutions)):
    statfile="cylindrical_%s.stat" % str(resolutions[resolution])
    data[resolution]  = stat(statfile)['Fields']['VelocityError']['surface_l2norm%Bottom'][-1] / anal
    if(resolution > 0):
        cmb_velocity_convergence[resolution] = np.log2(data[resolution-1] / data[resolution])
        print (cmb_velocity_convergence[resolution])
