from scipy import *
from pylab import *
##Set initial x/y coordinates of detectors
x = 0.5+0.25*arange(0,100.)/100.
y = zeros(100) + 0.5

##Set timestep parameters and number of cycles
t = 0.
n_cycles = 1
dt = 0.1/n_cycles
tmax = 8

##Velocity function defined for a rotating steady-state velocity field
def vel(x,y,t):
    return [-(y-0.5),x-0.5]

##Classical 4th order Runge-Kutta scheme
##k1 defined at initial coordinates and time
##k2 defined at initial coordinates plus a half step based on k1, and initial time plus a half step of dt
##k3 defined at initial coordinates plus a half step based on k2, and initial time plus a half step of dt
##k4 defined at initial coordinates plus a full step based on k3, and initial time plus a full step of dt
##Final weightings of each parameter introduced from the 4th order Butcher array
while(t<tmax):
    t = t + dt
    [k1_x,k1_y] = vel(x,y,t)
    [k2_x,k2_y] = vel(x+0.5*dt*k1_x,y+0.5*dt*k1_y,t+0.5*dt)
    [k3_x,k3_y] = vel(x+0.5*dt*k2_x,y+0.5*dt*k2_y,t+0.5*dt)
    [k4_x,k4_y] = vel(x+dt*k3_x,y+dt*k3_y,t+dt)
    x = x + dt*(k1_x/6.+k2_x/3. + k3_x/3. + k4_x/6.)
    y = y + dt*(k1_y/6.+k2_y/3. + k3_y/3. + k4_y/6.)

plot(x,y,'.')
show()

x.tofile('Xvals.txt',sep=' ')
y.tofile('Yvals.txt',sep=' ')
