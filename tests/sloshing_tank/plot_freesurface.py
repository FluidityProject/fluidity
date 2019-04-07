#!/usr/bin/env python

import settings
import ana_sol

import sys
import math
import commands
import matplotlib.pyplot as plt
import getopt


from scipy.special import erf
from numpy import poly1d
from matplotlib.pyplot import figure, show
from numpy import pi, sin, linspace
from matplotlib.mlab import stineman_interp
from numpy import exp, cos

from fluidity_tools import stat_parser as stat

# Usage
def usage():
        print("plt_freesurface.py --file=detectorfile")
        print("All the other options are read from settings.py")



################# Main ###########################
def main(argv=None):

        a_0 = settings.a0 # initial maximum perturbation
        g = settings.g # gravity
        eta= settings.eta # viscosity
        L= settings.L # wavelength
        timestep= settings.timestep # timestep 
        filename=''

        global debug
        debug=False
        #debug=True      

        try:                                
                opts, args = getopt.getopt(sys.argv[1:], "h:", ['file='])
        except getopt.GetoptError:  
                usage()                     
                sys.exit(2)                     
        for opt, arg in opts:                
                if opt == '--file':      
                    filename=arg
                elif opt == '-h' or opt == '--help':
                    usage()                     
                    sys.exit(2) 
        if filename=='':
                usage()                     
                sys.exit(2) 

        print('Using:\n\ta_0 =', a_0 )# initial maximum perturbation
        print('\tg =', g )# gravity
        print('\teta=', eta )# viscosity
        print('\tL=', L )# wavelength
        print('\ttimestep=', timestep )# timestep 

        
        ####################### Print time plot  ###########################
        print('Generating time plot')

        x_time= stat(filename)["ElapsedTime"]["value"]
        fs_simu= stat(filename)["water"]["FreeSurface"]["left"]
#        fs_simu= stat(filename)["water"]["FreeSurface"]["middle"]
        fs_ana = stat(filename)["water"]["FreeSurface_Analytical"]["left"]
#        fs_ana = stat(filename)["water"]["FreeSurface_Analytical"]["middle"]

        plt.ion() # swith on interactive mode
        fig = figure()
        ax = fig.add_subplot(111)

        ax.plot(x_time,fs_simu,'ro')
        ax.plot(x_time,fs_ana,'-')
        plt.title('Free Surface timeplot at x=0')
        plt.xlabel('Time [s]')
        plt.ylabel('Free surface [m]')

        plt.draw()
        raw_input("Please press Enter")
        #plt.cla()



if __name__ == "__main__":
    main()


