from scipy import *
from fluidity_tools import stat_parser
from pylab import show
s = stat_parser("lagrangian_detectors.detectors")

data = zeros((2,100,size(s['Steve_001']['position'][0])))

for i in range(100):
    n = i + 1
    padding = ''
    if(n<100):
        padding = '0'
    if(n<10):
        padding ='00'
    data[0,i,:] = s['Steve_'+padding+str(n)]['position'][0]
    data[1,i,:] = s['Steve_'+padding+str(n)]['position'][1]

data[0,:,-1].tofile('Xvals.txt',sep=' ')
data[1,:,-1].tofile('Yvals.txt',sep=' ')

show()
