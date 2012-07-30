#!/usr/bin/env python

from pylab import *
from math import atan2, pi
import argparse
import vtktools

parser = argparse.ArgumentParser(description='Project a .vtu from the sphere to a plane.')
parser.add_argument('input_vtu', nargs='+',
                   help='The .vtu file(s) to be projected.')
parser.add_argument('-p', dest='output_prefix',default='plane_',
                   help='The prefix of the output vtu (default: plane_)')
parser.add_argument('-d', dest='depth_scale',default=100, type=float,
                   help='Units for vertical length scale in meters (default: 100)')
parser.add_argument('-r', dest='earth_radius',default=None, type=float,
                   help='Radius of the earth, used to determine z=0 level (default: maximum l2norm of input x,y,z coordinates)')
parser.add_argument('-inv', dest='inverse', action='store_true', default=False,
                    help='Perform inverse operation, from plane to spherical (placeholder).')

args = parser.parse_args()

deg2rad = 180/pi

if args.inverse:
    raise Exception('Inverse operation not coded yet.')

print 'Performing spherical to planar'
for vtu in range(size(args.input_vtu)):
    vtu_name = args.input_vtu[vtu]
    output_filename = args.output_prefix+vtu_name
    vtu_object = vtktools.vtu(vtu_name)
    print "Done importing ", vtu_name

    if args.earth_radius is None:
      rearth = sqrt(sum((vtu_object.GetLocations())**2,axis=1).max())
      print "Earth radius = ", rearth
    else:
      rearth = args.earth_radius

    npoints = vtu_object.ugrid.GetNumberOfPoints()
    for i in range(npoints):
        x,y,z = vtu_object.ugrid.GetPoint(i)
        radius = sqrt( x**2 + y**2 + z**2 )
        theta = arccos(z/radius)
        theta = pi/2.0 - theta
        phi = atan2(y,x)
        new_x = phi*deg2rad
        new_y = theta*deg2rad
        new_z = (radius-rearth)/args.depth_scale
        vtu_object.ugrid.GetPoints().SetPoint(i,new_x,new_y,new_z)

    if output_filename.endswith('.pvtu'):
        output_filename = output_filename[:-4]+'vtu'

    vtu_object.Write(filename=output_filename)
    print 'Projected vtu written to', output_filename
