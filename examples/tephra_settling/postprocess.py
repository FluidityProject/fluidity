## This Python script searches through the current directory, loads the maximum Tephra::Velocity values
## from the .stat file, and then prints out the results in a .pdf file.

from fluidity_tools import stat_parser
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy
import os

# List the contents of the current directory
dir_list = os.listdir(".")
# If the stat file exists, then read from it
if("tephra_settling.stat" in dir_list):
   s = stat_parser("tephra_settling.stat")
   time = s["ElapsedTime"]["value"]
   tephra_u_max = s["Tephra"]["Velocity%magnitude"]["max"]

   pylab.plot(time, tephra_u_max, label='Numerical results')

   # Print it all out as a pretty picture
   pylab.xlabel('Time (s)')
   pylab.ylabel('Maximum tephra velocity (m/s)')
   pylab.legend(loc=2)

   pylab.savefig('tephra_velocity.pdf')
   
   print("Data plotted and saved in tephra_velocity.pdf")
else:
   print("No .stat file found - cannot plot data.")
   
