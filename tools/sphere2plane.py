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
parser.add_argument('-inv', dest='inverse', action='store_true', default=False,
                    help='Perform inverse operation, from plane to spherical (placeholder).')

args = parser.parse_args()

stretch_factor = 10000.0

if args.inverse:
    print 'Inverse operation not coded yet.'
else:
    print 'Performing spherical to planar'
    for vtu in range(size(args.input_vtu)):
        vtu_name=args.input_vtu[vtu]
        output_filename=args.output_prefix+vtu_name
        vtu_object= vtktools.vtu(vtu_name)
        print "Done importing ", vtu_name

        npoints = vtu_object.ugrid.GetNumberOfPoints()
        dist = vtu_object.GetScalarField('DistanceToTop')
        for i in range(npoints):
            (x,y,z) = vtu_object.ugrid.GetPoint(i)
            radius = sqrt( x**2 + y**2 + z**2 )
            theta = arccos(z/radius) 
            phi = atan2(y,x) 
            phi = pi/2.0 - phi
            new_x=phi*stretch_factor
            new_y=theta*stretch_factor
            new_z=radius
            vtu_object.ugrid.GetPoints ().SetPoint (i,new_x,new_y,new_z)

        if output_filename.endswith('.pvtu'):
            output_filename = output_filename[:-4]+'vtu'

        vtu_object.Write(filename=output_filename)
        print 'Projected vtu written to', output_filename
