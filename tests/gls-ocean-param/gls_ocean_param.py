#!/usr/bin/env python

import string
import re 
import os
import vtktools
from numpy import argsort, size
from pylab import cm, figure, colorbar, arange, date2num, mpl, xlabel, ylabel, savefig, close, legend
from datetime import datetime, timedelta

#### taken from http://www.codinghorror.com/blog/archives/001018.html    #######
def sort_nicely( l ): 
    """ Sort the given list in the way that humans expect. 
    """ 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    l.sort( key=alphanum_key ) 
    
##############################################################################

def get_vtus(directory):
    """ For a given directory, list all VTU files, except checkpoint & _0
    """ 
    files = []
    dirList=os.listdir(directory)
    for file in dirList:
        if (string.find(file, 'checkpoint') != -1):
            continue
        elif (string.find(file, 'Mesh') != -1):
            continue
        elif (string.find(file, 'vtu') != -1):
            files.append(directory + "/" + file)
        else:
            continue

    sort_nicely(files)
    return files

##############################################################################

def get_1d_indices(pos, x0=0, y0=0, tolerance=1.0e-5):
    """ Return the field indices corresponding to the ordered depth values at position (x0, y0)
    """
    ind = argsort(-pos[:,2])
    indices = []
    for i in ind:
        if (x0-tolerance < pos[i][0] < x0+tolerance and y0-tolerance < pos[i][1] < y0+tolerance):
            indices.append(i)
    return indices

##############################################################################
#### Mixed layer depth
##############################################################################

def calc_mld_den(density, depth, den0=0.03):
    """ Calculate Mixed Layer depth based on density
    """
    # grab surface density
    sden = density[0]
    # grab any values where the UML condition is met for density 
    uml_den = ((density) <= (sden+den0)).nonzero()
    
    if( (size(uml_den) > 0 ) ):
        LL = uml_den[-1][-1]
        mld = depth[LL]
    else:
        mld = 0.0
        
    return mld

def calc_mld_tke(tke, depth, tke0=1.0e-5):
    """ Calculate Mixed Layer depth based on TKE
    """
    uml_tke = ((tke) >= tke0).nonzero()
    
    if ( (size(uml_tke) > 0) ):
        LL = uml_tke[-1][-1]
        mld = depth[LL]
    else:
        mld = 0.0
        
    return mld

def calc_mld(files,start,x0=0.0,y0=0.0):
    """ Caclulate density-based MLD from a bunch of VTU files
    """

    mld = []
    times = []
    dates = []
    for file in files:
      
        try:
            os.stat(file)
        except:
            print("No such file: %s" % file)
            sys.exit(1)

        # open vtu and derive the field indices of the edge at (x=0,y=0) ordered by depth
        u=vtktools.vtu(file)
        pos = u.GetLocations()
        ind = get_1d_indices(pos, x0, y0)
    
        # from this we can derive the 1D profile of any field like this:
        depth = vtktools.arr([-pos[i,2] for i in ind])
    
        # handle time for different types of plots
        time = u.GetScalarField('Time')
        times.append(time[0])   # seconds
        dates.append( date2num(start + timedelta(seconds=time[0])) ) # integer datetime
    
        # grab density profile and calculate MLD_den (using 2 different deviation parameters
        d = u.GetScalarField('Density')
        den = vtktools.arr( [d[i] * 1000 for i in ind] )
        mld.append( calc_mld_den(den, depth) ) #den0 = 0.03 is default


    return mld, times, dates


def calc_mld_tke_files(files,start,x0=0.0,y0=0.0):
    """ Caclulate tke-based MLD from a bunch of VTU files
    """

    mld = []
    times = []
    dates = []
    for file in files:
      
        try:
            os.stat(file)
        except:
            print("No such file: %s" % file)
            sys.exit(1)

        # open vtu and derive the field indices of the edge at (x=0,y=0) ordered by depth
        u=vtktools.vtu(file)
        pos = u.GetLocations()
        ind = get_1d_indices(pos, x0, y0)
    
        # from this we can derive the 1D profile of any field like this:
        depth = vtktools.arr([-pos[i,2] for i in ind])
    
        # handle time for different types of plots
        time = u.GetScalarField('Time')
        times.append(time[0])   # seconds
        dates.append( date2num(start + timedelta(seconds=time[0])) ) # integer datetime
    
        # grab density profile and calculate MLD
        d = u.GetScalarField('GLSTurbulentKineticEnergy')
        tke = vtktools.arr( [d[i] for i in ind] )
        mld.append( calc_mld_tke(tke, depth) )


    return mld, times, dates

def calc_mld_temp(temp, depth, temp0=0.1):
    """ Calculate Mixed Layer depth based on temperature
    """
    # grab surface density
    sst = temp[0]
    # grab any values where the UML condition is met for density 
    uml_temp = ((temp) > (sst-temp0)).nonzero()
    
    if( (size(uml_temp) > 0 ) ):
        LL = uml_temp[-1][-1]
        mld = depth[LL]
    else:
        mld = 0.0
        
    return mld

##############################################################################
#### Plotting functions
##############################################################################

def plot_2d_data(data,depths,time_secs,start_date,file_path,axis_label,finish_date=None,mld_data=None,max_depth=150,interval=3,minimum=None,maximum=None,spacing=None,colour_scale=cm.jet,dates=None):
    """
    """
    # turn given 2d-arrays into numpy arrays (in case they are not already)
    data = vtktools.arr(data)
    time_secs = vtktools.arr(time_secs)
    depths = vtktools.arr(depths)
    
    # convert time profiles in seconds into months
    start = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
    if (dates == None):
        dates = time_secs
        i = 0
        for time in time_secs:
            t = float(time[0].item())
            dates[i,:] = date2num(start + timedelta(seconds=t))
            i += 1
        
    # see if finishing date is given, default to last time given
    if (finish_date != None):
        finish = datetime.strptime(finish_date, "%Y-%m-%d %H:%M:%S")
    else:
        finish = dates[-1][0] 
    
    # define min/max and spacing of data if not given (so we see all of the data)
    if (minimum == None):
        minimum = data.min() 
        minimum = minimum - (0.1*minimum) 
    if (maximum == None):
        maximum = data.max()
        maximum = maximum + (0.1*maximum) 
    if (spacing == None):
        spacing = (maximum - minimum) /256.

    # plot 2d colour graph...
    fig = figure(figsize=(15,8),dpi=90)
    ax = fig.add_axes([.1,.18,.9,.7])
    cs=ax.contour(dates, depths, data, arange(minimum,maximum,spacing),cmap=colour_scale)
    cs=ax.contourf(dates, depths, data, arange(minimum,maximum,spacing),cmap=colour_scale)
    pp=colorbar(cs,format='%.2f')
    if(mld_data!=None):
        ax.plot(dates[:,0],mld_data,'w', alpha=0.7)
    
    dateFmt = mpl.dates.DateFormatter('%m/%Y')
    ax.xaxis.set_major_formatter(dateFmt)
    monthsLoc = mpl.dates.MonthLocator(interval=interval)
    ax.xaxis.set_major_locator(monthsLoc)
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    ax.set_ylim(max_depth, 0)
    ax.set_xlim(start,finish)
    pp.set_label(axis_label)
    xlabel('Date (mm/yyyy)')
    ylabel('Depth (m)')
    
    form = file_path.split('.')[-1].strip()
    savefig(file_path, dpi=90,format=form)
    close(fig)
    
def plot_1d_comparison(data_dict,style_dict,time_dict,start_date,finish_date,file_path,axis_label,interval=3):
    """ 
    """
    start_time = date2num(datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S"))
    finish_time = date2num(datetime.strptime(finish_date, "%Y-%m-%d %H:%M:%S"))
        
    # plot 1d graph...
    fig = figure(figsize=(15,8),dpi=90)
    ax = fig.add_axes([.05,.12,.9,.85])
    max_value = 0.0
    for key, data_arr in data_dict.iteritems():
        ax.plot(time_dict[key],data_arr,style_dict[key], label=key)
        data_arr = vtktools.arr(data_arr)
        if data_arr.max() > max_value:
            max_value = data_arr.max()
    max_value += max_value * 0.1
    
    dateFmt = mpl.dates.DateFormatter('%m/%Y')
    ax.xaxis.set_major_formatter(dateFmt)
    monthsLoc = mpl.dates.MonthLocator(interval=interval)
    ax.xaxis.set_major_locator(monthsLoc)
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    ax.set_ylim(max_value, 0)
    ax.set_xlim(start_time,finish_time)
    xlabel('Date (mm/yyyy)')
    ylabel(axis_label)
    legend(loc=0)
    
    form = file_path.split('.')[-1].strip()
    savefig(file_path, dpi=90,format=form)
    close(fig)


# width defines bar width
# percent defines current percentage
def progress(width, percent):
    import math
    import sys
    marks = math.floor(width * (percent / 100.0))
    spaces = math.floor(width - marks)
    
    loader = '[' + ('=' * int(marks)) + (' ' * int(spaces)) + ']'
      
    sys.stdout.write("%s %d%%\r" % (loader, percent))
    if percent >= 100:
        sys.stdout.write("\n")
    sys.stdout.flush()
