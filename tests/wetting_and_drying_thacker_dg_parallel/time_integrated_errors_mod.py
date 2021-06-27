#!/usr/bin/python

from fluidity_tools import stat_parser
import getopt
import sys
from math import ceil

def get_time_integrated_fs_error(filename):
     
        s = stat_parser(filename)
        # one wave lasts 1/((8*9.81*50/(430620*430620))^0.5)*2*pi=43193 seconds. 

        timesteps=s["ElapsedTime"]["value"]
        timestep=timesteps[1]-timesteps[0]
        print "Found ", len(timesteps), " timesteps with dt=", timestep
        if timesteps[-1]<43193:
                print 'Simulation must run for at least 43193s'
                exit()
        lastindex=int(ceil(43193.0/float(timestep)))
        print 'But only the first ', lastindex, ' are considered for the error integration.'
        print ''

        print 'Integrating errors over time...'

        l2error=s["water"]["fs_error_abs"]["l2norm"]
        timel2error=l2error[0]+l2error[lastindex-1]
        timel2error=timel2error/2
        for i in range(1,lastindex-1):
                timel2error=timel2error+l2error[i]
        timel2error=timel2error*timestep
        print "FreeSurface - Time integrated l2error: ", timel2error

        l1error=s["water"]["fs_error_abs"]["integral"]
        timel1error=l1error[0]+l1error[lastindex-1]
        timel1error=timel1error/2
        for i in range(1,lastindex-1):
                timel1error=timel1error+l1error[i]
        timel1error=timel1error*timestep
        print "FreeSurface - Time integrated l1error: ", timel1error

        Inferror=s["water"]["fs_error_abs"]["max"]
        timeInferror=Inferror[0]+Inferror[lastindex-1]
        timeInferror=timeInferror/2
        for i in range(1,lastindex-1):
                timeInferror=timeInferror+Inferror[i]
        timeInferror=timeInferror*timestep
        print "FreeSurface - Time integrated Inf Error: ", timeInferror


        print "Done"                
        return [timel1error, timel2error, timeInferror]

def get_time_integrated_vel_error(filename):

        s = stat_parser(filename)

        # one wave lasts 1/((8*9.81*50/(430620*430620))^0.5)*2*pi=43193 seconds. 

        timesteps=s["ElapsedTime"]["value"]
        timestep=timesteps[1]-timesteps[0]
        print "Found ", len(timesteps), " timesteps with dt=", timestep
        if timesteps[-1]<43193:
                print 'Simulation must run for at least 43193s'
                exit()
        lastindex=int(ceil(43193.0/float(timestep)))
        print 'But only the first ', lastindex, ' are considered for the error integration.'
        print ''

        print 'Integrating errors over time...'

        l2error=s["water"]["vel_error_abs%magnitude"]["l2norm"]
        timel2error=l2error[0]+l2error[lastindex-1]
        timel2error=timel2error/2
        for i in range(1,lastindex-1):
                timel2error=timel2error+l2error[i]
        timel2error=timel2error*timestep
        print "Velocity - Time integrated l2error: ", timel2error

        timel1error=-1.0 # Not available here :(

        Inferror=s["water"]["vel_error_abs%magnitude"]["max"]
        timeInferror=Inferror[0]+Inferror[lastindex-1]
        timeInferror=timeInferror/2
        for i in range(1,lastindex-1):
                timeInferror=timeInferror+Inferror[i]
        timeInferror=timeInferror*timestep
        print "Velocity - Time integrated Inf Error: ", timeInferror


        print "Done"                
        return [timel1error, timel2error, timeInferror]


