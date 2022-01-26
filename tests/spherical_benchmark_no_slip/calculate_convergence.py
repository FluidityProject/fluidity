#!/usr/bin/env python3

from fluidity_tools import stat_parser as stat
import numpy as np

# Statfile from finest mesh for relative errors:
statfile_fine="sphere_C.stat"

def velocity_convergence(file_coarse,file_fine):
    ######## Velocity Errors: #######
    anal = stat(statfile_fine)['Fields']['AnalyticalVelocity%magnitude']['l2norm'][-1]
    error_coarse  = stat(file_coarse)['Fields']['VelocityError%magnitude']['l2norm'][-1] / anal
    error_fine    = stat(file_fine)['Fields']['VelocityError%magnitude']['l2norm'][-1] / anal
    velocity_convergence  = np.log2(error_coarse / error_fine)
    return velocity_convergence

def surface_velocity_convergence(file_coarse,file_fine):
    ######## Surface Velocity Errors: #######
    anal = stat(statfile_fine)['Fields']['AnalyticalVelocity']['surface_l2norm%Top'][-1]
    error_coarse  = stat(file_coarse)['Fields']['VelocityError']['surface_l2norm%Top'][-1] / anal
    error_fine    = stat(file_fine)['Fields']['VelocityError']['surface_l2norm%Top'][-1] / anal
    surface_velocity_convergence  = np.log2(error_coarse / error_fine)
    return surface_velocity_convergence

def CMB_velocity_convergence(file_coarse,file_fine):
    ######## CMB Velocity Errors: #######
    anal = stat(statfile_fine)['Fields']['AnalyticalVelocity']['surface_l2norm%Bottom'][-1]
    error_coarse  = stat(file_coarse)['Fields']['VelocityError']['surface_l2norm%Bottom'][-1] / anal
    error_fine    = stat(file_fine)['Fields']['VelocityError']['surface_l2norm%Bottom'][-1] / anal
    CMB_velocity_convergence  = np.log2(error_coarse / error_fine)
    return CMB_velocity_convergence

def pressure_convergence(file_coarse,file_fine):
    ######## Pressure Errors: #######
    anal = stat(statfile_fine)['Fields']['AnalyticalPressure']['l2norm'][-1]
    error_coarse  = stat(file_coarse)['Fields']['PressureError']['l2norm'][-1] / anal
    error_fine    = stat(file_fine)['Fields']['PressureError']['l2norm'][-1] / anal
    pressure_convergence  = np.log2(error_coarse / error_fine)
    return pressure_convergence

def normalstress_convergence(file_coarse,file_fine):
    ######## Normalstress Errors: #######
    anal = stat(statfile_fine)['Fields']['AnalyticalNormalStress']['l2norm'][-1]
    error_coarse  = stat(file_coarse)['Fields']['NormalStressError']['l2norm'][-1] / anal
    error_fine    = stat(file_fine)['Fields']['NormalStressError']['l2norm'][-1] / anal
    normalstress_convergence  = np.log2(error_coarse / error_fine)
    return normalstress_convergence

def surface_normalstress_convergence(file_coarse,file_fine):
    ######## Surface Normalstress Errors: #######
    anal = stat(statfile_fine)['Fields']['AnalyticalNormalStress']['surface_l2norm%Top'][-1]
    error_coarse  = stat(file_coarse)['Fields']['NormalStressError']['surface_l2norm%Top'][-1] / anal
    error_fine    = stat(file_fine)['Fields']['NormalStressError']['surface_l2norm%Top'][-1] / anal
    surface_normalstress_convergence  = np.log2(error_coarse / error_fine)
    return surface_normalstress_convergence

def CMB_normalstress_convergence(file_coarse,file_fine):
    ######## CMB Normalstress Errors: #######
    anal = stat(statfile_fine)['Fields']['AnalyticalNormalStress']['surface_l2norm%Bottom'][-1]
    error_coarse  = stat(file_coarse)['Fields']['NormalStressError']['surface_l2norm%Bottom'][-1] / anal
    error_fine    = stat(file_fine)['Fields']['NormalStressError']['surface_l2norm%Bottom'][-1] / anal
    CMB_normalstress_convergence  = np.log2(error_coarse / error_fine)
    return CMB_normalstress_convergence
