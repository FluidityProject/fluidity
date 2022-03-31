#!/usr/bin/env python3

import numpy
from numpy import array,argsort,corrcoef,size
import os
import sys
import vtktools
import math
import re 
from datetime import *
import csv

# max number of fluidity depth points
FNN = 100

# max number of times
TT = 365

# Origin (and where data is taken from)
x0 = 0
y0 = 0
z0 = 0
z_end = 200 # metres from surface

# definitions of UML depth
den0 = 0.03

# data file that contains measuered UML depths
uml_data = "mld_stationpapa.csv"

# start date/time
start_date = '1970/01/02'
end_date = '1971/01/01'
start_datetime = datetime(1970,1,1)

#### taken from http://www.codinghorror.com/blog/archives/001018.html  #######
def tryint(s):
    try:
        return int(s)
    except:
        return s
    
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
##############################################################################


###### Get MLD through time from a bunch of VTUs #############################
def calc_mld(files,mld_depths):
    for file in files:
        try:
            os.stat(file)
        except:
            print("No such file: %s" % file)
            sys.exit(1)
        num = int(file.split(".vtu")[0].split('_')[-1])
        u=vtktools.vtu(file)
        pos = u.GetLocations()
        time = u.GetScalarField('Time')
        tt = time[0]

        # grab the data e need from the VTUs and shove in an array
        den = u.GetScalarField('Density')
        xyz_data = []
        for i in range(0,len(den)):
            if (x0-0.1 < pos[i,0] < x0+0.1 and y0-0.1 < pos[i,1] < y0+0.1):
                xyz_data.append((pos[i,0],pos[i,1],-pos[i,2]+z0,1000*den[i]))
       

        # sorted the array based on depth
        xyz_data = vtktools.arr(xyz_data)
        III = argsort(xyz_data[:,2])
        xyz_data_sorted = xyz_data[III,:]

        # Surface values
        sden = xyz_data_sorted[0,3]

        # grab any values where the UML condition is met for temperature, density and TKE
        uml_den = ((xyz_data_sorted[:,3]) <= (sden+den0)).nonzero()

        # ...on density
        if( (size(uml_den) > 0 ) ):
            LL = uml_den[-1][-1]
            zz = xyz_data_sorted[:,2]
            if (LL+1 < size(zz)):
                zza = zz[LL+1]
                kea = xyz_data_sorted[LL+1,3]
            else:
                zza = zz[LL]
                kea = xyz_data_sorted[LL,3]
            zzb = zz[LL]
            keb = xyz_data_sorted[LL,3]
            tt = tt/(24*60*60)
            time = start_datetime + timedelta(days=tt)
            key = '%04d/%02d/%02d' % (time.year, time.month, time.day)
            mld_depths[key] = (zza-(zzb-zza)*(((sden+den0)-kea)/(kea-keb)))
##############################################################################

### For a given directory, list all VTU files, except checkpoint & _0 ########
def get_vtus(directory):

    files = []
    dirList=os.listdir(directory)
    for file in dirList:
        if (str.find(file, 'checkpoint') != -1):
            continue
        elif (str.find(file, '_0.vtu') != -1):
            continue
        elif (str.find(file, '.vtu') != -1):
            files.append(directory + "/" + file)
        else:
            continue

    sort_nicely(files)  

    return files

##############################################################################

def get_error_estimates():
    mld_depths = {}
    files = get_vtus(".")
    calc_mld(files,mld_depths)

    # real data for the UML depths...#
    with open(uml_data) as csvfile:
        reader = csv.reader(csvfile)
        uml_dates = []
        uml_d = []
        for row in reader:
            uml_dates.append(row[1][6:10]+"/"+row[1][3:5]+"/"+row[1][0:2])
            uml_d.append(float(row[6]))
    
        last_date = 0
        uml_depths = {}
        uml_count = {}
        i = 0
    
        # average out those days that have > 1 measurement
        for date in uml_dates:
            if date in uml_depths:
                uml_depths[date] += uml_d[i]
                uml_count[date] += 1
            else:
                uml_depths[date] = uml_d[i]
                uml_count[date] = 1
            i += 1
    
        for key in uml_depths:
            uml_depths[key] = uml_depths[key]/uml_count[key]
    
        # lets go to work...
    
        mld_icom = []
        uml_measured = []
        distance = []
        time = []
        keys = sorted(uml_depths.keys())
        for key in keys:
            if (key > end_date):
                break
            if (key < start_date):
                continue
            distance.append(abs(uml_depths[key] - mld_depths[key]))
            time.append(key)
            mld_icom.append(mld_depths[key])
            uml_measured.append(uml_depths[key])
    
        r =  corrcoef(uml_measured,mld_icom)[1,0]
        n = len(keys)
        dist =  sum(distance)/n
        ErrorMetric = [r,dist]
    
        return ErrorMetric
