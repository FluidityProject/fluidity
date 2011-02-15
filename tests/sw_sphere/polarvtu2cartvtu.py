import sys
sys.path.append('/home/jemma/src/swfl/python')

import vtktools
import spheretools
from numpy import zeros

ug=vtktools.vtu('kelvin_wave_0.vtu')
x=ug.GetLocations()

f=spheretools.polar2cart(1.0)
def ff(X,t):
    return f(X[0:2],t)

ug.ApplyCoordinateTransformation(ff)

U_manifold=ug.GetField('')

ug.AddVectorField()
ug.Write('cartsphere.vtu')
