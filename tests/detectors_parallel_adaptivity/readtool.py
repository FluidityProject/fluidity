from fluidity_tools import stat_parser
from numpy import zeros,fromfile
def readstat(file_name):
    s = stat_parser(file_name)

    last_locations_error = zeros((2,100))
    for i in range(100):
        n = i + 1
        padding = ''
        if(n<100):
            padding = '0'
            if(n<10):
                padding ='00'
        # +0.5 due to domain change from the lagrangian_detectors setup
        last_locations_error[0,i] = s['Steve_'+padding+str(n)]['position'][0][-1]+0.5
        last_locations_error[1,i] = s['Steve_'+padding+str(n)]['position'][1][-1]+0.5
    X = fromfile('Xvals.txt',sep=' ')
    Y = fromfile('Yvals.txt',sep=' ')
    last_locations_error[0,:] = last_locations_error[0,:] - X
    last_locations_error[1,:] = last_locations_error[1,:] - Y
    return last_locations_error
