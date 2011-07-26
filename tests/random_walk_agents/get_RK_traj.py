from scipy import *
from pylab import *
x = 0.5+0.25*arange(0,100.)/100.
y = zeros(100) + 0.5

t = 0.
n_cycles = 1
dt = 0.1/n_cycles
tmax = 8

def vel(x,y):
    return [-(y-0.5),x-0.5]

while(t<tmax):
    t = t + dt
    [k1_x,k1_y] = vel(x,y)
    x = x + dt*(k1_x)
    y = y + dt*(k1_y)

plot(x,y,'.')
show()

x.tofile('Xvals.txt',sep=' ')
y.tofile('Yvals.txt',sep=' ')
