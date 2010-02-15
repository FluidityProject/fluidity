import os
import sys
import vtktools
import math
import numpy
from numpy import finfo


def flux(file,x,y):

    u=vtktools.vtu(file)
    flux = u.GetScalarField('HeatFlux')
    pos = u.GetLocations()
    f = finfo(float)

    for i in range(0,len(flux)):
        if( abs(pos[i,0] - x) < f.eps and abs(pos[i,1] - y) < f.eps and (pos[i,2] - 0.0) < f.eps ):
            return flux[i]

    return -666





