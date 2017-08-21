from pylab import *
import numpy as np

##Set initial x/y coordinates of detectors
x = 0.5+0.25*arange(0,100.)/100.
y = zeros(100) + 0.5

##Set timestep parameters:
tmax = 8

r     = 0.25*np.arange(0,100.)/100.
theta = tmax**2./16.
X     = r * np.cos(theta) + 0.5
Y     = r * np.sin(theta) + 0.5

plot(X,Y,'.')
#print (x-X).max(), (y-Y).max()
show()

#X.tofile('Xvals.txt',sep=' ')
#Y.tofile('Yvals.txt',sep=' ')
np.save("Xvals.npy",X)
np.save("Yvals.npy",Y)

