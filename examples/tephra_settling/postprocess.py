# This Python script searches through the current directory, loads the maximum
# Tephra::Velocity values from the .stat file, and then prints out the results in a .pdf
# file.
import os

import matplotlib.pyplot as plt
from fluidity_tools import stat_parser

# List the contents of the current directory
dir_list = os.listdir(".")
# If the stat file exists, then read from it
if "tephra_settling.stat" in dir_list:
    s = stat_parser("tephra_settling.stat")
    time = s["ElapsedTime"]["value"]
    tephra_u_max = s["Tephra"]["Velocity%magnitude"]["max"]

    plt.plot(time, tephra_u_max, label="Numerical results")

    # Print it all out as a pretty picture
    plt.xlabel("Time (s)")
    plt.ylabel("Maximum tephra velocity (m/s)")
    plt.legend(loc=2)

    plt.savefig("tephra_velocity.pdf")

    print("Data plotted and saved in tephra_velocity.pdf")
else:
    print("No .stat file found - cannot plot data.")
