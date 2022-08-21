#!/usr/bin/env python3
import glob

import matplotlib.pyplot as plt
import vtktools
from numpy import arange, concatenate, newaxis, ones, zeros

colx = ones((21, 1)) * 0.1
colz = zeros((21, 1))
coly = arange(-0.5, 0.51, 0.05)[:, newaxis]
coordinates = concatenate((colx, coly, colz), 1)

vtus = glob.glob1("./", "*.vtu")

for vtu in vtus:
    vtunumber = vtu[19:]
    vtunumber = vtunumber[:-4]
    vtufile = vtktools.vtu("./" + vtu)
    water = vtktools.vtu.ProbeData(
        vtufile, coordinates, "Water::MaterialVolumeFraction"
    )
    droplet = vtktools.vtu.ProbeData(
        vtufile, coordinates, "Droplet::MaterialVolumeFraction"
    )
    air = vtktools.vtu.ProbeData(vtufile, coordinates, "Air::MaterialVolumeFraction")
    waterdropletsum = vtktools.vtu.ProbeData(
        vtufile, coordinates, "Water::SumMaterialVolumeFractions"
    )
    plt.plot(
        coordinates[:, 1],
        water,
        coordinates[:, 1],
        droplet,
        coordinates[:, 1],
        air,
        coordinates[:, 1],
        waterdropletsum,
    )
    plt.axis([-0.55, 0.55, -0.5, 1.5])
    plt.xlabel("Distance")
    plt.ylabel("VolumeFraction")
    plt.legend(("Water", "Droplet", "Air", "Sum"), loc="center left")
    if len(vtunumber) == 1:
        vtunumber = "000" + vtunumber
    elif len(vtunumber) == 2:
        vtunumber = "00" + vtunumber
    elif len(vtunumber) == 3:
        vtunumber = "0" + vtunumber
    plt.savefig("./" + vtu[:19] + vtunumber + ".png")
    plt.cla()
