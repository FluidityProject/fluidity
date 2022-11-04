#!/usr/bin/env python3
"""
   Python script to measure the position of the particle cloud front in
   simulations of the shock tube experiment by Rogue et al. (1998).
   Written by Christian Jacobs.
"""

import matplotlib.pyplot as plt
import numpy
import vtktools


def get_cloud_front(path_prefix, n_files, resolution, threshold):
    # Size of the domain
    ly = 4.5

    # Arrays to store the results
    # Note: n_files is the number of VTU files to be read.
    upper_positions = numpy.zeros(n_files)
    lower_positions = numpy.zeros(n_files)
    times = numpy.zeros(n_files)

    path = path_prefix  # + "/output/"

    # For each vtu file...
    for i in range(0, n_files, 20):
        # Open the VTU files one by one
        try:
            file = vtktools.vtu(
                path + "/mphase_rogue_shock_tube_dense_bed_nylon_" + str(i) + ".pvtu"
            )
        except Exception:
            print("WARNING: Could not open VTU file!")
            times[i] = 0
            continue
        file.GetFieldNames()

        time = max(file.GetScalarField("Gas::Time"))
        print("At time t = %f s" % time)

        # Generate a set of probe positions
        probes = numpy.linspace(1.35, 1.7, (ly / resolution))

        # Loop over each of the probes in the y direction, from the top down
        for j in range(len(probes) - 1, -1, -1):
            vfrac = vtktools.vtu.ProbeData(
                file,
                numpy.array([[0.0, probes[j], 0]]),
                "Particles::PhaseVolumeFraction",
            )
            if vfrac >= 0.3:
                upper_positions[i] = probes[j]
                break

        for j in range(0, len(probes)):
            vfrac = vtktools.vtu.ProbeData(
                file,
                numpy.array([[0.0, probes[j], 0]]),
                "Particles::PhaseVolumeFraction",
            )
            if vfrac >= 0.01:
                lower_positions[i] = probes[j]
                break

        times[i] = time - 0.0008

    return times, upper_positions, lower_positions


# get_cloud_front(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))

times, upper_positions, lower_positions = get_cloud_front("./", 501, 0.0085, 0.05)

f = open("upper.prn")
experimental_times = []
experimental_positions = []
for line in f:
    (t, p) = line.split("\t")
    experimental_times.append(float(t))
    experimental_positions.append(float(p) * 100)
f.close()
plt.plot(experimental_times, experimental_positions, "-r", label="Experimental (upper)")

f = open("lower.prn")
experimental_times = []
experimental_positions = []
for line in f:
    (t, p) = line.split("\t")
    experimental_times.append(float(t))
    experimental_positions.append(float(p) * 100)
f.close()
plt.plot(experimental_times, experimental_positions, "-b", label="Experimental (lower)")

for i in range(0, len(upper_positions)):
    upper_positions[i] = (
        upper_positions[i] - 1.35
    ) * 100  # Centre abscissa and convert to cm.
    lower_positions[i] = (
        lower_positions[i] - 1.35
    ) * 100  # Centre abscissa and convert to cm.
plt.plot(times, upper_positions, "ro", label="Fluidity (upper)")
plt.plot(times, lower_positions, "bo", label="Fluidity (lower)")
plt.axis([0, 0.004, 0, 10])
plt.title("Position of lower and upper front of a 2 cm particle bed")
plt.xlabel("Time (s)")
plt.ylabel("Distance from bottom of particle bed (cm)")
plt.legend(loc=2)
plt.savefig("test.png")
