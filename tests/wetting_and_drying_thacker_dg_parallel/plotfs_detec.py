#!/usr/bin/python

import vtktools
import sys
import math
import re 
import commands
import matplotlib.pyplot as plt
import getopt
import csv

from scipy.special import erf
from numpy import poly1d
from matplotlib.pyplot import figure, show
from numpy import pi, sin, linspace
from matplotlib.mlab import stineman_interp
from numpy import exp, cos
from fluidity_tools import stat_parser

# This script plots the 

def usage():
        print 'Usage:'
        print 'plotfs_detec.py -t timestep --file=detector_filename --csv=basename(optional) --save=filenamebase\n'
        print 'Specifying a timestepsize will change the timestep for the analytical solution only.'
        print '--save=true will save the plots as images instead of plotting them on the screen.'
        print '-t timestep will only plot the specified timestep.'
        print '--csv=basename instead of the detectors use csv files created with paraview (see code for more info)'


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
  return h0*(R**2-r2)/(R**2)+0.5-h0*(R**2-Rmesh**2)/(R**2) # positive number

################# Main ###########################
def main(argv=None):

        Rmesh=440000
        filename=''
        csvbasename=''
        timestep_ana=0.0
        dzero=0.8
        save='' # If nonempty, we save the plots as images instead if showing them
        handoffset=2.7

        global debug
        debug=False
        debug=True  
        verbose=False
        verbose_timestep=-1

        try:                                
                opts, args = getopt.getopt(sys.argv[1:], "t:h:", ['file=','csv=','save='])
        except getopt.GetoptError:  
                print "Getopterror :("
                usage()                     
                sys.exit(2)                     
        for opt, arg in opts:                
                if opt == '--file':      
                        filename=arg
                elif opt == '--csv':      
                        csvbasename=arg
                elif opt == '--save':
                        save=arg
                elif opt == '-t':
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
                if verbose_timestep>-1 and verbose_timestep!=t:
                        continue
                #if not t%10==0:
                #        continue
                fsvalues=[]
                xcoords=[]
                # check if we want to use detector (or vtu probe) data
                if csvbasename=='':
                        minfs=1000
                        for name, item in fs.iteritems():
                                xcoords.append(s[name]['position'][0][0])
                                
                                fsvalues.append(fs[name][t])
                                # alternativly: Probe data with vtktools:
                                # d0=1.0
                                # fsvalues.append(probe_fs(filename, [xcoords[-1],0, -bathymetry_function([xcoords[-1],0])+d0/2], t))

                                print 'Time: ', t, ', Name: ',name, ', X-Coord: ', s[name]['position'][0][0], ', FS: ', fs[name][t]
                                if minfs>fs[name][t]:
                                        minfs=fs[name][t]
                        print 'Minfs: ', minfs

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
                        ax2.plot(xcoords,fsvalues,'r.', label='Numerical solution')
                        plt.ylim(-10,0)

                # use csv data
                else:
                        csvreader=csv.reader(open(csvbasename+'.'+str(t)+'.csv', 'rb'), delimiter=',', quotechar='|')
                        firstrow=True
                        fs_id=-1
                        pos_x_id=-1
                        pos_y_id=-1
                        # get indices
                        for row in csvreader:
                                if firstrow:
                                        firstrow=False
                                        for i in range(0,len(row)):
                                                if row[i]=='"FreeSurface"':
                                                        fs_id=i
                                                elif row[i]=='"Points:0"':
                                                        pos_x_id=i
                                                elif row[i]=='"Points:1"':
                                                        pos_y_id=i
                                        assert(fs_id>=0 and pos_x_id>=0 and pos_y_id>=0)        
                                        continue
                                if abs(float(row[pos_y_id]))>=1e-5:
                                        print 'You didnt save the csv file on a slice on the x axis!'
                                        exit()
                                xcoords.append(float(row[pos_x_id]))
                                fsvalues.append(float(row[fs_id])+handoffset)

                        # sort coords and values 
                        temp=[]
                        for i in range(0,len(xcoords)):
                                temp.append([xcoords[i], fsvalues[i]])
                        temp.sort()
                        xcoords=[]
                        fsvalues=[]
                        for i in range(0,len(temp)):
                                xcoords.append(temp[i][0]/1000) # in km
                                fsvalues.append(temp[i][1])
                        # Plot result of one timestep
                        ax2.plot(xcoords,fsvalues, label='Numerical solution', linestyle='--')
                        plt.ylim(-10,0)



                
                # Plot Analytical solution
                fsvalues_ana=[]
                xcoords=[]
                        
                offset=-bathymetry_function([0.0,0.0])+50.0+dzero

                for i in range(-250,250): # number of points for the analytical solution
                        xcoords.append(Rmesh*float(i)/250/1000) # in km
                        fsvalues_ana.append(analytic_solution(xcoords[-1]*1000, t*timestep_ana)+offset+handoffset)

                ax2.plot(xcoords, fsvalues_ana, label='Analytical solution', linestyle='-')
                #       plt.ylim(-10,0)

                       # plt.title('d_min=2.5m, Time='+ str(t*timestep)+'s')
                plt.xlabel('Position [km]')
                plt.ylabel('Free surface [m]')
                plt.legend(loc='lower center')

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



# tries to retrieve the values with vtu probes. does not work :(
def probe_fs(filename, coords, t):
        import sys
        import vtktools
        vtu=vtktools.vtu(filename[:-10]+'_'+str(t)+'.pvtu') 
        vtu.ApplyProjection('x', 'y', '0')   # Set the third component to zero. 
        #Xlist = arange(0.0,4.001,0.01)# x co-ordinates of the desired array shape
        #Ylist = arange(0.0,2.001,0.01)# y co-ordinates of the desired array shape
        #[X,Y] = meshgrid(Xlist,Ylist) 
        #Z = 0.0*X # This is 2d so z is an array of zeros.
        #X = reshape(X,(size(X),))
        #Y = reshape(Y,(size(Y),))
        #Z = reshape(Z,(size(Z),))
        #zip(X,Y,Z)
        #pts = zip(X,Y,Z)
        print "coords: ", coords
        pts = vtktools.arr([coords])
        # create arrays of velocity and temperature values at the desired points
        fs = vtu.ProbeData(pts, 'FreeSurface')
        return fs[0]


# reads the free surface values saved in csv files. for that you have to open paraview, apply a slice filter and then "save data" for all timesteps as csv.
def csv_reader(filename, coords, t):
        import sys
        import vtktools
        vtu=vtktools.vtu(filename[:-10]+'_'+str(t)+'.pvtu') 
        vtu.ApplyProjection('x', 'y', '0')   # Set the third component to zero. 
        #Xlist = arange(0.0,4.001,0.01)# x co-ordinates of the desired array shape
        #Ylist = arange(0.0,2.001,0.01)# y co-ordinates of the desired array shape
        #[X,Y] = meshgrid(Xlist,Ylist) 
        #Z = 0.0*X # This is 2d so z is an array of zeros.
        #X = reshape(X,(size(X),))
        #Y = reshape(Y,(size(Y),))
        #Z = reshape(Z,(size(Z),))
        #zip(X,Y,Z)
        #pts = zip(X,Y,Z)
        print "coords: ", coords
        pts = vtktools.arr([coords])
        # create arrays of velocity and temperature values at the desired points
        fs = vtu.ProbeData(pts, 'FreeSurface')
        return fs[0]


if __name__ == "__main__":
    main()



        
# Make video from the images:
        
# mencoder "mf://*.png" -mf type=png:fps=30 -ovc lavc -o output.avi

