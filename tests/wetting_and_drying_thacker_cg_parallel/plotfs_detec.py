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
        print 'plotfs_detec.py -v timestep --file=detector_filename --save=filenamebase -t timestepsize(optional)\n'
        print 'Specifying a timestepsize will change the timestep for the analytical solution only.'
        print '--save=true will save the plots as images instead of plotting them on the screen.'
        print '-v timestep will print some verbos information about the given timestep and then exit.'


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
  Rmesh=450000
  r2=X[0]**2+X[1]**2
  # The last term is an offset term
  return h0*(R**2-r2)/(R**2)+0.5-h0*(R**2-Rmesh**2)/(R**2)

################# Main ###########################
def main(argv=None):

        Rmesh=450000
        filename=''
        timestep_ana=0.0
        dzero=0.5
        save='' # If nonempty, we save the plots as images instead if showing them

        global debug
        debug=False
        debug=True  
        verbose=False
        verbose_timestep=-1

        try:                                
                opts, args = getopt.getopt(sys.argv[1:], "v:t:", ['file=','save='])
        except getopt.GetoptError:  
                usage()                     
                sys.exit(2)                     
        for opt, arg in opts:                
                if opt == '--file':      
                        filename=arg
                elif opt == '--save':
                        save=arg
                elif opt == '-t':
                    timestep_ana = float(arg)
                elif opt == '-v':
                        verbose=True
                        verbose_timestep=int(arg)
                elif opt == '-h' or opt == '--help':
                    usage()                     
                    sys.exit(2) 
        if filename=='':
                print 'No filename specified. You have to give the detectors filename.'
                usage()   
                sys.exit(2) 

        
        ####################### Print time plot  ###########################
        print 'Generating time plot'
      
        s = stat_parser(filename)

        timesteps=s["ElapsedTime"]["value"]
        timestep=timesteps[1]-timesteps[0]
        print "Found ", len(timesteps), " timesteps with dt=", timestep
        if timestep_ana==0.0:
                timestep_ana=timestep


        fs=s["water"]["FreeSurface"]
        print "Found ", len(fs), " detectors. We assume they are equidistant distributed over the domain (", -Rmesh, "-", Rmesh, ")."


        # Get and plot results
        plt.ion() # swith on interactive mode
        fig2 = figure()
        ax2 = fig2.add_subplot(111)

        for t in range(0,len(timesteps)):
                if verbose and verbose_timestep!=t:
#                        verbose_timestep=verbose_timestep-1
                        continue

                fsvalues=[]
                xcoords=[]
                for name, item in fs.iteritems():
                        #print name
                        xcoords.append(s[name]['position'][0][0])
                        #print xcoord
                        fsvalues.append(fs[name][t])

                # sort coords and values 
                temp=[]
                for i in range(0,len(xcoords)):
                        temp.append([xcoords[i], fsvalues[i]])
                temp.sort()
                xcoords=[]
                fsvalues=[]
                for i in range(0,len(temp)):
                        xcoords.append(temp[i][0])
                        fsvalues.append(temp[i][1])

                # Plot result of one timestep
                ax2.plot(xcoords,fsvalues,'ro', label='Numerical solution')
                plt.ylim(-10,0)
                
                # Plot Analytical solution
                fsvalues_ana=[]
                
                offset=-bathymetry_function([0.0,0.0])+50.0+dzero

                xcoords.sort()
                for x in xcoords:
                        fsvalues_ana.append(analytic_solution(x, t*timestep_ana)+offset)

                ax2.plot(xcoords, fsvalues_ana, label='Analytical solution')
                plt.ylim(-10,0)

                plt.title('d_min=2.5m, Time='+ str(t*timestep)+'s')
                plt.xlabel('X position [m]')
                plt.ylabel('Free surface [m]')
                plt.legend()

                if verbose:
                        print 'Current time: ', timesteps[t]
                        print 'Numerical position at domain center: ', fsvalues[len(fs)/2]
                        print 'Error in the domain center: ', abs(fsvalues_ana[len(fs)/2]-fsvalues[len(fs)/2])
                        deltax=Rmesh/(len(fs)-1)
                        integral=0.0
                        for i in range(0,len(xcoords)):
                                integral=integral+abs(fsvalues_ana[i]-fsvalues[i])
                        integral=integral-fsvalues_ana[0]/2-fsvalues_ana[len(xcoords)-1]/2
                        integral=integral*deltax
                        print 'Integral error: ', integral
                        exit()

                if save=='':
                        plt.draw()
                        raw_input("Please press Enter")
                else:
                        plt.savefig(save+'_'+str(100000+t+1)+'.pdf', facecolor='white', edgecolor='black', dpi=100)
                plt.cla()
                t=t+1
        if verbose:
                print 'Ups, program ended unexpected. Are you sure you specified the correct time for the -v argument?'
        
# Make video from the images:
        
# mencoder "mf://*.png" -mf type=png:fps=30 -ovc lavc -o output.avi




if __name__ == "__main__":
    main()


