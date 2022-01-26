import numpy as np
import h5py

def readstat_3d():
    f = h5py.File("lagrangian_detectors.particles.Steve.h5part", "r")

    num_detectors = 200000
    last_locations_error = np.zeros((3, num_detectors))
    d = f["/Step#{}".format(len(f) - 1)]
    idx = np.argsort(d["id"])

    X = np.fromfile("Xvals.txt",sep=" ")
    Y = np.fromfile("Yvals.txt",sep=" ")
    Z = 0.5 * np.ones_like(X)

    for i, (dim, ref) in enumerate(zip("xyz", (X, Y, Z))):
        last_locations_error[i,:] = d[dim][:][idx] - ref

    return last_locations_error
