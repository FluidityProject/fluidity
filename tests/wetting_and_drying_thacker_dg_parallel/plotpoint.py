#!/usr/bin/python

import vtktools
import sys
import math
import re 
import commands
import matplotlib.pyplot as plt
import getopt
import csv
import ana_sol

from scipy.special import erf
from numpy import poly1d
from matplotlib.pyplot import figure, show
from numpy import pi, sin, linspace
from matplotlib.mlab import stineman_interp
from numpy import exp, cos
from fluidity_tools import stat_parser

def usage():
        print 'Usage:'
        print 'plotpoint.py -r radius --file=detector_filename --csv=basename(optional) --save=filenamebase -t timestepsize(optional) --d0=d0 -b begintimestep -e endtimestep\n'
        print 'Specifying a timestepsize will change the timestep for the analytical solution only.'
        print '--save=true will save the plots as images instead of plotting them on the screen.'
        print '-r radius. Specifies the point where the measurements will be performed. Can be positive or negative.'
        print '--csv=baseame instead of the detectors use csv files created with paraview (see code for documentation)'


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

################# Main ###########################
def main(argv=None):

        Rmesh=440000
        filename=''
        csvbasename=''
        timestep_ana=0.0
        dzero=0.0
        save='' # If nonempty, we save the plots as images instead if showing them

        global debug
        debug=False
        debug=True  
        verbose=False
        verbose_timestep=-1
        radius=0
        begintimestep=0
        endtimestep=100000000000000

        try:                                
          opts, args = getopt.getopt(sys.argv[1:], "v:t:r:b:e:", ['file=','csv=','save=', 'd0='])
        except getopt.GetoptError:  
                print "Getopterror :("
                usage()                     
                sys.exit(2)                     
        for opt, arg in opts:                
                if opt == '--file':      
                        filename=arg
                if opt == '--csv':      
                        csvbasename=arg
                elif opt == '--save':
                        save=arg
                elif opt == '--d0':
                        dzero=float(arg)
                elif opt == '-t':
                    timestep_ana = float(arg)
                elif opt == '-r':
                        radius=float(arg)
                elif opt == '-h' or opt == '--help':
                    usage()                     
                    sys.exit(2)
                elif opt == '-b':
                  begintimestep=int(arg)
                elif opt == '-e':
                  endtimestep=int(arg)
 
                   
        if filename=='':
                print 'No filename specified. You have to give the detectors filename.'
                usage()   
                sys.exit(2) 

        
        ####################### Print time plot  ###########################
        print 'Generating time plot'
      
        s = stat_parser(filename)

        endtimestep=min(endtimestep, len(s["ElapsedTime"]["value"])-1)
        timesteps=s["ElapsedTime"]["value"][begintimestep:endtimestep]
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

        timeplot_xvalue=[]  # this should be equidistant timestep size
        timeplot_yvalue=[]  # these are the scalar values at the specified point
        timeplot_yvalue_analytical=[]  # these are the analytical scalar values at the specified point
        r_index = -1   # this is the index to access the detector or csv point closest to the user specified point
        handoffset = 2.8  # a user specified offset for the y axis

        for t in range(0,len(timesteps)):

                if not t%10==0:
                        continue

                fsvalues=[]

                # check if we want to use detector (or vtu probe) data
                if csvbasename=='':
                        print "detecors based plot not supported yet. use csv instead!"
                        exit
#                        minfs=1000
#                        for name, item in fs.iteritems():
#                                xcoords.append(s[name]['position'][0][0])
#                                
#                                fsvalues.append(fs[name][t])
#                                # alternativly: Probe data with vtktools:
#                                # d0=1.0
#                                # fsvalues.append(probe_fs(filename, [xcoords[-1],0, -bathymetry_function([xcoords[-1],0])+d0/2], t))
#
#                                print 'Time: ', t, ', Name: ',name, ', X-Coord: ', s[name]['position'][0][0], ', FS: ', fs[name][t]
#                                if minfs>fs[name][t]:
#                                        minfs=fs[name][t]
#                        print 'Minfs: ', minfs
#
#                        # sort coords and values 
#                        temp=[]
#                        for i in range(0,len(xcoords)):
#                                temp.append([xcoords[i], fsvalues[i]])
#                        temp.sort()
#                        xcoords=[]
#                        fsvalues=[]
#                        for i in range(0,len(temp)):
#                                xcoords.append(temp[i][0])
#                                fsvalues.append(temp[i][1])
#
#                        # Plot result of one timestep
#                        ax2.plot(xcoords,fsvalues,'r.', label='Numerical solution')
#                        plt.ylim(-10,0)

                # use csv data
                else:
                        csvreader=csv.reader(open(csvbasename+'.'+str(t)+'.csv', 'rb'), delimiter=',', quotechar='|')
                        firstrow=True
                        fs_id=-1
                        pos_x_id=-1
                        pos_y_id=-1
                        xcoords=[]
                        fsvalues=[]
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
                                xcoords.append(temp[i][0])
                                fsvalues.append(temp[i][1])

                        # find the nearest csv point for the user specified point
                        for i in range(0,len(xcoords)):
                                if xcoords[i]<radius:
                                        r_index=i
                                else:
                                        break
                          if i==len(xcoords)-1:
                                print "The radius specified is too big!"
                                exit
                         print "Evaluating at x=", xcoords[i]
                        timeplot_xvalue.append(timesteps[t]/60/60)
                        timeplot_yvalue.append(fsvalues[i])
                        offset=-ana_sol.bathymetry_function(0.0)+50.0+dzero
                        timeplot_yvalue_analytical.append(analytic_solution(xcoords[i], timesteps[t])+offset+handoffset)
                



        # Plot result of all timesteps
        ax2.plot(timeplot_xvalue,timeplot_yvalue,'-', label='Numerical solution')
        ax2.plot(timeplot_xvalue,timeplot_yvalue_analytical,'-', label='Analytical solution')

        #plt.ylim(-10,0)

        plt.xlabel('Time [h]')
        plt.ylabel('Free surface [m]')
        plt.legend()

        plt.draw()
        raw_input("Please press Enter")

                
                # Plot Analytical solution
#                fsvalues_ana=[]
#                xcoords=[]
#                
#                offset=-bathymetry_function([0.0,0.0])+50.0+dzero#
#
#                for i in range(-250,250): # number of points for the analytical solution
#                        xcoords.append(Rmesh*float(i)/250)
 #                       fsvalues_ana.append(analytic_solution(xcoords[-1], t*timestep_ana)+offset)
#
 #               ax2.plot(xcoords, fsvalues_ana, label='Analytical solution')
  #              plt.ylim(-10,0)
#
 #               plt.title('d_min=2.5m, Time='+ str(t*timestep)+'s')
  #              plt.xlabel('X position [m]')
   #             plt.ylabel('Free surface [m]')
    #            plt.legend()
#
 #               if verbose:
  #                      print 'Current time: ', timesteps[t]
   #                     print 'Numerical position at domain center: ', fsvalues[len(fs)/2]
    #                    print 'Error in the domain center: ', abs(fsvalues_ana[len(fs)/2]-fsvalues[len(fs)/2])
     #                   deltax=Rmesh/(len(fs)-1)
      #                  integral=0.0
       #                 for i in range(0,len(xcoords)):
        #                        integral=integral+abs(fsvalues_ana[i]-fsvalues[i])
         #               integral=integral-fsvalues_ana[0]/2-fsvalues_ana[len(xcoords)-1]/2
          #              integral=integral*deltax
           #             print 'Integral error: ', integral
            #            exit()
#
 #               if save=='':
  #                      plt.draw()
   #                     raw_input("Please press Enter")
    #            else:
     #                   plt.savefig(save+'_'+str(100000+t+1)+'.pdf', facecolor='white', edgecolor='black', dpi=100)
      #          plt.cla()
       #         t=t+1
        #if verbose:
         #       print 'Ups, program ended unexpected. Are you sure you specified the correct time for the -v argument?'



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

