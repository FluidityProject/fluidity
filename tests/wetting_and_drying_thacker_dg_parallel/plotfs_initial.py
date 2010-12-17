#!/usr/bin/python

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

def usage():
        print 'Usage:'
        print 'plotfs_initial.py\n'

def analytic_solution(r, t):
        h0=50.0
        eta0=2.0
        R=430620.0
        g=9.81

        A=((h0+eta0)**2-h0**2)/((h0+eta0)**2+h0**2)
        omega=(8.0*g*h0/(R**2))**0.5
        eta=(1.0-A**2)**0.5/(1.0-A*cos(omega*t))
        eta=eta-1
        eta=eta-r**2/(R**2)*((1.0-A**2)/((1.0-A*cos(omega*t))**2)-1.0)
        eta=h0*eta

        # check if the free surface is lower than the cup and if, then set it to the cup heigth
        h=h0*(R**2-r**2)/R**2
        
        return max(eta,-h)

# should be copied from the diamond extrude function. X is 2 dimensional
def bathymetry_function(X):
  h0=50
  R=430620
  Rmesh=440000
  r2=X[0]**2+X[1]**2
  # The last term is an offset term
  return h0*(R**2-r2)/(R**2)+0.5-h0*(R**2-Rmesh**2)/(R**2)

################# Main ###########################
def main(argv=None):

        Rmesh=440000.0
        timestep_ana=0.0
        save='' # If nonempty, we save the plots as images instead if showing them

        global debug
        debug=False
        debug=True  
        verbose=False
        verbose_timestep=-1

        
          # Get and plot results
        plt.ion() # swith on interactive mode
        fig2 = figure()
        ax2 = fig2.add_subplot(111)

        for t in range(0,1):

                # Plot Analytical solution
                fsvalues_ana=[]
                bathymetry=[]
                
                offset=-bathymetry_function([0.0,0.0])+50.0
                xcoords=[]
                for i in range(-1000,1000):
                        xcoords.append(Rmesh*float(i)/1000.0)
                        fsvalues_ana.append(analytic_solution(xcoords[-1], t*timestep_ana)+offset)
                        bathymetry.append(-bathymetry_function([xcoords[-1],0]))

                ax2.plot(xcoords, fsvalues_ana, '-', label='Initial free surface')
                ax2.plot(xcoords, bathymetry, 'k-', linewidth=4, label='Bottom')

                plt.legend(loc=3)

                if save=='':
                        plt.draw()
                        raw_input("Please press Enter")
                else:
                        plt.savefig(save+'_'+str(100000+t+1)+'.pdf', facecolor='white', edgecolor='black', dpi=100)
                plt.cla()


if __name__ == "__main__":
    main()


