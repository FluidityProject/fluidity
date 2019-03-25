#!/usr/bin/env python

from __import__ import print_function
from pylab import *
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

print('Performing spherical to planar')
for vtu in range(size(args.input_vtu)):
    vtu_name = args.input_vtu[vtu]
    output_filename = args.output_prefix+vtu_name
    vtu_object = vtktools.vtu(vtu_name)
    print("Done importing ", vtu_name)

    pos = vtu_object.GetLocations()
    radius = sqrt(sum(pos**2, axis=1))
    if args.earth_radius is None:
      rearth = radius.max()
    else:
      rearth = args.earth_radius
    print("Earth radius = ", rearth)

    theta = arccos(pos[:,2]/radius)
    theta = pi/2.0 - theta
    phi = arctan2(pos[:,1],pos[:,0])
    new_pos = array([phi*deg2rad, theta*deg2rad, (radius-rearth)/args.depth_scale]).T

    for i in range(len(new_pos)):
        vtu_object.ugrid.GetPoints().SetPoint(i,new_pos[i])

    if output_filename.endswith('.pvtu'):
        output_filename = output_filename[:-4]+'vtu'

    vtu_object.Write(filename=output_filename)
    print('Projected vtu written to', output_filename)
