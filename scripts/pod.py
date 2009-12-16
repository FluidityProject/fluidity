#!/usr/bin/env python

# Compute the Principal Orthogonal Decomposition
# of a series of snapshots.

import sys
import vtktools
import scipy
import scipy.linalg
import numpy
import math
import operator

# sys.argv[1] should be the name of a VTU
# with fields like
# Temperature1, Temperature2, Temperature3, etc
# to be PODded.
# sys.argv[2] should be the field to POD
# (in this case, "Temperature")

vtu = vtktools.vtu(sys.argv[1])
field_name = sys.argv[2]

# figure out the number of snapshots
K = 1
while True:
  try:
    vtu.GetScalarField(field_name + "%s" % K)
    K = K + 1
  except:
    break
 
K = K - 1
N = vtu.ugrid.GetNumberOfPoints()

# Compute the mean:
mean = numpy.zeros((N))
for k in range(K):
  field = vtu.GetScalarField(field_name + "%s" % (k+1))
  mean = mean + field
mean = [mean[i]/K for i in range(len(mean))]

# Assemble the matrix
A = numpy.mat(numpy.zeros((N, K)))
for k in range(K):
  field = vtu.GetScalarField(field_name + "%s" % (k+1))
  A[:, k] = (field - mean).reshape(N, 1)

ATA = A.T * A

(evals, evecs) = scipy.linalg.eig(ATA)
evals = [float(x) for x in evals]
(evals, evecs) = zip(*sorted(zip(evals, evecs), key=operator.itemgetter(0), reverse=True))
print evals

sum_evals = sum(evals)
running_sum = 0.0
for k in range(K):
  running_sum = running_sum + evals[k]
  print "I(%s): %s" % (k+1, running_sum / sum_evals)

  # Numpy is a bit rubbish. It doesn't even recognise matrix-vector multiplication.
  phi = (A * evecs[k].reshape(evecs[k].shape + (1,))) / math.sqrt(evals[k])
  vtu.AddScalarField("POD%s" % (k+1), phi)

vtu.Write()
