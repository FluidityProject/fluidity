#!/usr/bin/env python

from __future__ import print_function

import vtktools
import sys
import math
import re 
import commands
import matplotlib.pyplot as plt
import getopt

from scipy.special import erf
from numpy import poly1d
from matplotlib.pyplot import figure, show
from numpy import pi, sin, linspace
from matplotlib.mlab import stineman_interp
from numpy import exp, cos
from fluidity_tools import stat_parser




def mirror(x):
        return 13800-x


def usage():
        print('Usage:')
        print('plotfs_detec.py [-w] --file=detector_filename --save=filename')
        print('--save=... saves the plots as images instead of plotting them on the screen.')
        print('-w plots the wetting procedure (drying is default).')


# should be copied from the diamond extrude function. X is 2 dimensional
def bathymetry_function(X):
  return -5.0*X/13800

################# Main ###########################
def main(argv=None):

        filename=''
        timestep_ana=0.0
        dzero=0.01
        save='' # If nonempty, we save the plots as images instead if showing them
        wetting=False

        try:                                
                opts, args = getopt.getopt(sys.argv[1:], ":w", ['file=','save='])
        except getopt.GetoptError:  
                usage()                     
                sys.exit(2)                     
        for opt, arg in opts:                
                if opt == '--file':      
                        filename=arg
                elif opt == '--save':
                        save=arg
                elif opt == '-w':
                        wetting=True
        if filename=='':
                print('No filename specified. You have to give the detectors filename.')
                usage()   
                sys.exit(2) 

        
        ####################### Print time plot  ###########################
        print('Generating time plot')
      
        s = stat_parser(filename)

        timesteps=s["ElapsedTime"]["value"]
        timestep=timesteps[1]-timesteps[0]
        print("Found ", len(timesteps), " timesteps with dt=", timestep)
        if timestep_ana==0.0:
                timestep_ana=timestep


        fs=s["water"]["FreeSurface"]
        print("Found ", len(fs), " detectors. We assume they are equidistant distributed over the domain (", 0, "-", 13800, ").")


        # Get and plot results
        plt.ion() # swith on interactive mode
        fig2 = figure()
        ax2 = fig2.add_subplot(111)

        if wetting:
                ##plot_start=90 # in timesteps
                plot_start=18 # in timesteps, after 18 timesteps the waterlevel reaches its lowest point
                ##plot_end=114 # in timesteps
                plot_end=54 # in timesteps
                plot_name='Wetting'
        else:
                plot_start=54  # in timesteps
                plot_end=89  # in timesteps
                plot_name='Drying'

        

        for t in range(0,len(timesteps)):
                # ignore the first waveperiod
                if t<plot_start:
                        continue
                if t>plot_end:
                        continue
                fsvalues=[]
                xcoords=[]
                for name, item in fs.iteritems():
                        #print name
                        xcoords.append(mirror(s[name]['position'][0][0]))
                        #print xcoord
                        fsvalues.append(fs[name][t])

                # Plot result of one timestep
                ax2.plot(xcoords,fsvalues,'r,', label='Numerical solution')

                # Plot Analytical solution
                fsvalues_ana=[]
                
                offset=-bathymetry_function(0.0)+dzero

                xcoords.sort()
                for x in xcoords:
                        fsvalues_ana.append(bathymetry_function(mirror(x))-offset)

                # Plot vertical line in bathmetry on right boundary
                xcoords.append(xcoords[len(xcoords)-1]+0.000000001)
                fsvalues_ana.append(2.1)

                ax2.plot(xcoords, fsvalues_ana, 'k', label='Bathymetry')

                #plt.legend()
                if t==plot_end:
                        # change from meters in kilometers in the x-axis
                        # return locs, labels where locs is an array of tick locations and
                        # labels is an array of tick labels.
                        locs, labels = plt.xticks()
                        for i in range(0,len(locs)):
                                labels[i]=str(locs[i]/1000)
                        plt.xticks(locs, labels)

                        plt.ylim(-2.2,1.4)
                        #plt.title(plot_name)
                        plt.xlabel('Position [km]')
                        plt.ylabel('Free surface [m]')

                        if save=='':
                                plt.draw()
                                raw_input("Please press Enter")
                        else:
                                plt.savefig(save+'_'+plot_name+'.pdf', facecolor='white', edgecolor='black', dpi=100)
                        plt.cla()
                t=t+1
                
        
# Make video from the images:
        
# mencoder "mf://*.png" -mf type=png:fps=30 -ovc lavc -o output.avi





if __name__ == "__main__":
    main()


