#!/usr/bin/env python

from fluidity_tools import stat_parser as stat
import numpy as np

# Wavenumber:
n = 2
resolutions = ["A","B"]
expected_normal_stress_convergence              = 1.5
expected_surface_normal_stress_convergence      = 3.0
expected_cmb_normal_stress_convergence          = 3.0

# Extract data
data                                            = np.zeros(len(resolutions))
normal_stress_convergence                       = np.zeros(len(resolutions))
surface_normal_stress_convergence               = np.zeros(len(resolutions))
cmb_normal_stress_convergence                   = np.zeros(len(resolutions))


# Statfile from finest mesh for relative errors:
statfile_fine="cylindrical_%s.stat" % str(resolutions[-1])

######## Normal_Stress Errors: #######
anal = stat(statfile_fine)['Fields']['AnalyticalNormalStress']['l2norm'][-1]
for resolution in range(len(resolutions)):
    statfile="cylindrical_%s.stat" % str(resolutions[resolution])
    data[resolution]  = stat(statfile)['Fields']['NormalStressError']['l2norm'][-1] / anal
    if(resolution > 0):
        normal_stress_convergence[resolution] = np.log2(data[resolution-1] / data[resolution])
        print (normal_stress_convergence[resolution])


######## Surface Normal_Stress Errors: #######
anal = stat(statfile_fine)['Fields']['AnalyticalNormalStress']['surface_l2norm%Top'][-1]
for resolution in range(len(resolutions)):
    statfile="cylindrical_%s.stat" % str(resolutions[resolution])
    data[resolution]  = stat(statfile)['Fields']['NormalStressError']['surface_l2norm%Top'][-1] / anal
    if(resolution > 0):
        surface_normal_stress_convergence[resolution] = np.log2(data[resolution-1] / data[resolution])
        print (surface_normal_stress_convergence[resolution])

######## CMB Normal_Stress Errors: #######
anal = stat(statfile_fine)['Fields']['AnalyticalNormalStress']['surface_l2norm%Bottom'][-1]
for resolution in range(len(resolutions)):
    statfile="cylindrical_%s.stat" % str(resolutions[resolution])
    data[resolution]  = stat(statfile)['Fields']['NormalStressError']['surface_l2norm%Bottom'][-1] / anal
    if(resolution > 0):
        cmb_normal_stress_convergence[resolution] = np.log2(data[resolution-1] / data[resolution])
        print (cmb_normal_stress_convergence[resolution])
