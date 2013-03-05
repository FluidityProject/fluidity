import sys
sys.path.append('../../../python')
import vtktools
import numpy
ug = vtktools.vtu('../InputCoordinates.vtu')
X = ug.GetLocations()
f = open('h_sphere_ico5_0360.node','r')
h = 0.0*X[:,0]
cnt = 0
for line in f:
    print line, cnt, numpy.size(X,0)
    vals = numpy.double(line.split())
    assert(vals.shape==(3,))
    h[cnt] = vals[2]
    cnt += 1
ug.AddScalarField('ExactH15Days', h)
ug.Write()
