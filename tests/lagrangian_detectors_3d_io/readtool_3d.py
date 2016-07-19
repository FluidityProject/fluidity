from fluidity_tools import stat_parser
from numpy import zeros,fromfile, ones, shape

def readstat_3d():
    s = stat_parser("lagrangian_detectors.detectors")

    num_detectors = 200000
    zfill_length  = int(len(str(num_detectors)))
    array_name    = "Steve_"
    last_locations_error = zeros((3,num_detectors))
    
    for i in range(num_detectors):
        n = i + 1
        last_locations_error[0,i] = s[array_name+str(n).zfill(zfill_length)]['position'][0][-1]
        last_locations_error[1,i] = s[array_name+str(n).zfill(zfill_length)]['position'][1][-1]
        last_locations_error[2,i] = s[array_name+str(n).zfill(zfill_length)]['position'][2][-1]
    X = fromfile('Xvals.txt',sep=' ')
    Y = fromfile('Yvals.txt',sep=' ')
    Z = 0.5*ones(shape(X))
    last_locations_error[0,:] = last_locations_error[0,:] - X
    last_locations_error[1,:] = last_locations_error[1,:] - Y
    last_locations_error[2,:] = last_locations_error[2,:] - Z 
    return last_locations_error
