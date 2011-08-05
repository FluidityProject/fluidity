## This Python script searches through the 'tests' directory and loads the maximum Tephra::Velocity values
## from the .stat file in each directory beginning with "mphase_tephra_settling_", and then prints out all 
## the results to a .pdf file.

from fluidity_tools import stat_parser
import matplotlib
matplotlib.use('pdf')
import pylab
import numpy
import os

# List the contents of the /tests/ directory
dir_list = os.listdir("/data/ctj10/convergence_analysis/tests")
# and loop over all of those directories beginning with "mphase_tephra_settling_"
for d in dir_list:
   if("mphase_tephra_settling_" in d):
      # List contents of the current directory
      inner_list = os.listdir("/data/ctj10/convergence_analysis/tests/" + d)
      
      # If the stat file exists, then read from it
      if("mphase_tephra_settling.stat" in inner_list):
         s = stat_parser("/data/ctj10/convergence_analysis/tests/" + d + "/mphase_tephra_settling.stat")
         time = s["ElapsedTime"]["value"]
         tephra_u_max = s["Tephra"]["Velocity%magnitude"]["max"]

         f = open('tephra_u_max.txt', 'w')

         # Write out the timesteps with the corresponding max Tephra::Velocity magnitude
         for i in range(0, len(time)):
            f.write(str(time[i]) + "\t" + str(tephra_u_max[i]) + "\n")
            
         f.close()

         pylab.plot(time, tephra_u_max, label=d)

# Print it all out as a pretty picture
pylab.xlabel('Time (s)')
pylab.ylabel('Maximum tephra velocity (m/s)')
pylab.legend(loc=2)

pylab.savefig('tephra_velocity_multiple.pdf')

