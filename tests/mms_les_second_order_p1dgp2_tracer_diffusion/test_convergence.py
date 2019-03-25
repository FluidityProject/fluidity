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

print('TracerError:')
print('-------------------------------------------')

for i, mesh in enumerate(meshes):

    a_error = stat("MMS_"+str(mesh[0])+".stat")["NS"]["TracerError"]["l2norm"][-1]
    b_error = stat("MMS_"+str(mesh[1])+".stat")["NS"]["TracerError"]["l2norm"][-1]

    ratio = a_error / b_error

    print(mesh[0] + '->' + mesh[1] + ': ', log(ratio, 2))

    convergence = log(ratio, 2)

print('-------------------------------------------')
