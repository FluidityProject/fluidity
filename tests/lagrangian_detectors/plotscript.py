import matplotlib.pyplot as plt
import numpy as np
from fluidity_tools import stat_parser

s = stat_parser("lagrangian_detectors.detectors")

data = np.zeros((2, 100, s["Steve_001"]["position"][0].size))

for i in range(100):
    n = i + 1
    padding = ""
    if n < 100:
        padding = "0"
    if n < 10:
        padding = "00"
    data[0, i, :] = s["Steve_" + padding + str(n)]["position"][0]
    data[1, i, :] = s["Steve_" + padding + str(n)]["position"][1]

data[0, :, -1].tofile("Xvals.txt", sep=" ")
data[1, :, -1].tofile("Yvals.txt", sep=" ")

plt.show()
