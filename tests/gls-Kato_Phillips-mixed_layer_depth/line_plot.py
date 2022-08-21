#!/usr/bin/env python3
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import vtktools
from matplotlib.ticker import MaxNLocator

filelist = sys.argv[1:]


x0 = 40.0


for file in filelist:
    print(file)
    try:
        os.stat(file)
    except FileNotFoundError:
        print("No such file: %s" % file)
        sys.exit(1)

    num = int(file.split(".vtu")[0].split("_")[-1])

    u = vtktools.vtu(file)

    timeee = u.GetScalarField("Time")
    Time = timeee[0]
    temp = u.GetScalarField("Temperature")
    vert_diff = u.GetScalarField("GLSVerticalDiffusivity")
    vert_visc = u.GetScalarField("GLSVerticalViscosity")
    uvw = u.GetVectorField("Velocity")
    kk = u.GetScalarField("GLSTurbulentKineticEnergy")
    pos = u.GetLocations()

    xyzkk = []
    for i in range(0, len(kk)):
        if abs(pos[i, 0] - x0) < 0.1:
            xyzkk.append(
                (
                    pos[i, 0],
                    -pos[i, 1],
                    pos[i, 2],
                    kk[i],
                    temp[i],
                    uvw[i, 0],
                    vert_diff[i],
                    vert_visc[i],
                )
            )

    xyzkkarr = np.array(xyzkk)
    III = np.argsort(xyzkkarr[:, 1])
    xyzkkarrsort = xyzkkarr[III, :]

    zzz = -xyzkkarrsort[:, 1]
    kkk = xyzkkarrsort[:, 3]
    ttt = xyzkkarrsort[:, 4]
    uuu = xyzkkarrsort[:, 5]
    diff = xyzkkarrsort[:, 6]
    visc = xyzkkarrsort[:, 7]

    fig = plt.figure()

    ax = fig.add_subplot(221)
    ax.plot(ttt * 10 * 2.0e-4, zzz, "b")
    ax.xaxis.set_major_locator(MaxNLocator(5))
    plt.xlabel("Buoyancy")
    plt.ylabel("Depth")

    ax = fig.add_subplot(222)
    ax.plot(uuu, zzz, "b")
    ax.xaxis.set_major_locator(MaxNLocator(5))
    plt.xlabel("U Velocity")
    plt.ylabel("Depth")

    ax = fig.add_subplot(223)
    ax.plot(kkk, zzz, "b")
    ax.xaxis.set_major_locator(MaxNLocator(5))
    plt.xlabel("TKE: k")
    plt.ylabel("Depth")

    ax = fig.add_subplot(224)
    ax.plot(diff, zzz, "b--")
    ax.plot(visc, zzz, "b")
    print(max(diff), max(visc))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    plt.xlabel("Vertical Diffusivity (--) and Viscosity (-)")
    plt.ylabel("Depth")

    plt.show()
