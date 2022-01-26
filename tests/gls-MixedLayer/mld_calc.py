from numpy import arange,concatenate,array,argsort
import os
import sys
import vtktools
import math
import re 
from math import sqrt



# compute the mixed layer depth over time
def MLD(filelist, unique_string):
    x0 = 166.6
    y0 = 50.
    tke0 = 1.0e-5
  
    times = []
    depths = []
    Dm = []
    nums = []

    # check for no files
    if (len(filelist) == 0):
        print("No files!")
        sys.exit(1)

    for file in filelist:
        try:
            os.stat(file)
        except:
            print("No such file: %s" % file)
            sys.exit(1)
        num = int(file.split(".vtu")[0].split('_')[-1])
        if (num < 1):
            continue
        nums.append(num)

    nums.sort()

    last_mld = 0.

    for num in nums:

        file = "gls-MixedLayer-"+unique_string+"_"+str(num)+".vtu"
     
        u=vtktools.vtu(file)

        time = u.GetScalarField('Time')
        tt = time[0]
        kk = u.GetScalarField('GLSTurbulentKineticEnergy')
        pos = u.GetLocations()
   
        xyzkk = []
        for i in range(0,len(kk)):
            if( abs(pos[i,0] - x0) < 0.1 and abs(pos[i,1] - y0) < 0.1):
                xyzkk.append((pos[i,0],pos[i,1],-pos[i,2],(kk[i])))

        xyzkkarr = vtktools.arr(xyzkk)
        III = argsort(xyzkkarr[:,2])
        xyzkkarrsort = xyzkkarr[III,:]
        # march down the column, grabbing the last value above tk0 and the first 
        # one less than tke0. Interpolate between to get the MLD
        kea = 1000
        keb = 0
        zza = 0
        zzb = 0
        for values in xyzkkarrsort:
            if (values[3] > tke0):
                kea = values[3]
                zza = -values[2]
            if (values[3] < tke0):
                keb = values[3]
                zzb = -values[2]
                break
        # We're using quite low resolution, so a better estimate of the MLD is
        # somewhere between the 2 layers
        mld = zza
        if (last_mld == mld):
            continue

        times.append(tt)
        depths.append(-1.0*mld)
        last_mld = mld
        Dm.append(1.05*0.01*(1.0/sqrt(0.01))*sqrt(tt))

    return times, depths, Dm

