from scipy import *
from pylab import *
num_detectors = 200000
x = 0.5+0.25*arange(0,float(num_detectors))/float(num_detectors)
y = zeros(num_detectors) + 0.5

t = 0.
n_cycles = 1
dt = 0.1/n_cycles
tmax = 16

def vel(x,y):
    return [-(y-0.5),x-0.5]

while(t<tmax):
    t = t + dt
    [k1_x,k1_y] = vel(x,y)
    [k2_x,k2_y] = vel(x+0.5*dt*k1_x,y+0.5*dt*k1_y)
    [k3_x,k3_y] = vel(x+0.5*dt*k2_x,y+0.5*dt*k2_y)
    [k4_x,k4_y] = vel(x+dt*k3_x,y+dt*k3_y)
    x = x + dt*(k1_x/6.+k2_x/3. + k3_x/3. + k4_x/6.)
    y = y + dt*(k1_y/6.+k2_y/3. + k3_y/3. + k4_y/6.)

plot(x,y,'.')
#show()

x.tofile('Xvals.txt',sep=' ')
y.tofile('Yvals.txt',sep=' ')
