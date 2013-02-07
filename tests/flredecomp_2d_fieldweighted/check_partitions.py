#! /usr/bin/env python

import sys
import numpy

num_partitions = 16

def calculate_partition_weights():
    partition_weights = numpy.zeros(num_partitions,dtype=float)
    
    for i in range(num_partitions):
        filename = 'Parallel-NP16-flredecomp-2d-fieldweighted_CoordinateMesh_%d.ele'%i
        elements = numpy.loadtxt(filename, skiprows=1)
        # Sum weights across each partition.
        # Note prognostic region ID == 27: weighting (in .flml) == 1.0
        # Prescribed regions ID == 24, 30: weighting (in .flml) == 0.01
        for j in range(elements.shape[0]):
            if(elements[j,4] == 27):
                partition_weights[i] = partition_weights[i] + 1.0
            else:
                partition_weights[i] = partition_weights[i] + 0.005

    min_weight = min(partition_weights)
    max_weight = max(partition_weights)
    max_load_imbalance = max_weight / min_weight

    return max_load_imbalance

