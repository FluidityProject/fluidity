#!/usr/bin/nv python3

import numpy as np
import math
import argparse
import pylab
import fileinput

slide_start_x = 112500
R = 150000 # slide total runout
U_max = 35.0
T = (math.pi / 2.0) * ( R / U_max)
max_h = 144.0
L = 223000
S = 1000
Ra = 75000
Rc = 0
Rd = 75000
Ta = math.pi*Ra / (2.0 * U_max)
Tc = Rc / U_max
Td = math.pi*Rd / (2.0 * U_max) 
cd = 0.0019

def main():

    parser = argparse.ArgumentParser(
         prog="test slide function",
         description="""Test the prescribed slide function"""
    )
    parser.add_argument(
        '-v', 
        '--verbose', 
        action='store_true', 
        help="Verbose output: mainly progress reports.",
        default=False
    )

    parser.add_argument(
        '-t',
        '--time',
        type=float,
        default=0.1,
        help="Which time to use"
    )
    parser.add_argument(
        '-a',
        '--animation',
        type=float,
        default=0,
        help="Do animation of slide dynamics from 0 to time specified"
    )
    parser.add_argument(
        'node_file',
        metavar='node_file',
        help='The node file which provides a list of coords on which to check the funciton'
    )

    args = parser.parse_args()
    verbose = args.verbose    
    node_file = args.node_file
    time = args.time
    anim = args.animation

    dt = 100
    x_coords = []
    slide_height = []
    vel = []
    i = 0
    # parse node file and get list of vertex coordinates
    for line in fileinput.input([node_file]):
        if (line[0] == '#'):
            break
        data = line.split()
        if (i > 0):
            X = [float(data[1])]
            x_coords.append(X[0])
        i=+1

    params = {
      'legend.fontsize': 18,
      'xtick.labelsize': 16,
      'ytick.labelsize': 16,
      'font.size' : 18,
      'axes.labelsize' : 18,
      'text.fontsize' : 18,
      'figure.subplot.hspace' : 0.5
    }
    pylab.rcParams.update(params)

    fig = pylab.figure(figsize=(15,8),dpi=90)
    ax = fig.add_subplot(111)  

    if (anim == 0):
        anim = time+1
        pylab.ioff()
    else:
        time = 0
        pylab.ion()
    for t in np.arange(time,anim,dt):
        i = 0
        vel = []
        h = []
        shear = []
        slide_front = set_slide_front(t)
        old_slide_front = set_slide_front(t-dt)
        for c in x_coords:
            h.append(set_slide_height([x_coords[i]],t,slide_front))
            vel.append((h[-1] - set_slide_height([x_coords[i]],t-dt,old_slide_front)) / dt)
            shear.append(set_slide_shear([x_coords[i]],t,((slide_front - old_slide_front) / dt),0.0,slide_front))
            i = i+1
        pylab.plot(x_coords,shear)
        pylab.draw()
        
    if (anim == time + 1):
        pylab.show()

def set_slide_shear(X,t,u_s,u_w,x_loc):

    form = slide_form(X,x_loc) 
    tau = 0.5 * 1000 * cd * (u_s - u_w)*form
    return tau


def set_slide_front(t):

    import math

    if t > T:
        t = T
    if t < 0:
        t = 0
    if t < Ta:
        s = Ra*(1-math.cos(U_max/Ra*t))
    elif Ta < t < Tc+Ta:
        s = Ra + U_max*(t - Ta)
    elif Ta+Tc < t < Ta+Tc+Td:
        s = Ra+Rc + Rd*math.sin(U_max/Rd * (t - Ta - Tc))
    else:
        s = R

    x_loc = slide_start_x + s

    return x_loc


def set_slide_height(X,t,x_loc):
    
    import math

    form = slide_form(X,x_loc)

    return max_h * form


def slide_form(X,x_loc):

    x_dash = (X[0] - x_loc)
    
    if (-(L+2*S)<x_dash<-(L+S)):
        h = math.exp(-(2*(x_dash+S+L)/S)**4)
    elif ( -(L+S)<=x_dash<-S):
        h = 1.0
    elif (-S<=x_dash<0):
        h = math.exp(-(2*(x_dash+S)/S)**4)
    else:
        h = 0

    return h


if __name__ == "__main__":
    main()
