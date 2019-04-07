#!/usr/bin/python

from fluidity_tools import stat_parser as stat
from vtktools import *
from math import log
import numpy as np

meshes = [['A','B'], ['B','C']]#, ['C','D']]

convergence = np.ones(2) * 1e10

print('')
print('ORDER OF CONVERGENCE')
print('-------------------------------------------')

print('VelocityError:')
print('-------------------------------------------')

for i, mesh in enumerate(meshes):

    a_error_x = stat("MMS_"+str(mesh[0])+".stat")["NS"]["VelocityError%1"]["l2norm"][-1]
    b_error_x = stat("MMS_"+str(mesh[1])+".stat")["NS"]["VelocityError%1"]["l2norm"][-1]
    a_error_y = stat("MMS_"+str(mesh[0])+".stat")["NS"]["VelocityError%2"]["l2norm"][-1]
    b_error_y = stat("MMS_"+str(mesh[1])+".stat")["NS"]["VelocityError%2"]["l2norm"][-1]

    ratio_x = a_error_x / b_error_x
    ratio_y = a_error_y / b_error_y

    print(mesh[0] + '->' + mesh[1] + ': ', [log(ratio_x, 2), log(ratio_y, 2)])

    convergence[0] = min(log(ratio_x, 2), log(ratio_y, 2), convergence[0])

print('-------------------------------------------')

print('EddyViscosityError:')
print('-------------------------------------------')

for i, mesh in enumerate(meshes):

    a_error = stat("MMS_"+str(mesh[0])+".stat")["NS"]["EddyViscosityError"]["l2norm"][-1]
    b_error = stat("MMS_"+str(mesh[1])+".stat")["NS"]["EddyViscosityError"]["l2norm"][-1]

    ratio = a_error / b_error

    print(mesh[0] + '->' + mesh[1] + ': ', log(ratio, 2))

    convergence[1] = min(log(ratio, 2), convergence[0])

print('-------------------------------------------')
