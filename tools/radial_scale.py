#!/usr/bin/env python

from pylab import *
import argparse
import vtktools

parser = argparse.ArgumentParser(description='Radially scale a .vtu file.')
parser.add_argument('input_vtu', nargs='+',
                   help='The .vtu file(s) to be scaled.')
parser.add_argument('-o', dest='output_prefix',default='',
                   help='The prefix of the output vtu (default behaviour replaces the input vtu)')
parser.add_argument('-s', dest='scale_factor', default='10',type=float,
                   help='The scale factor (default: 50)')

args = parser.parse_args()

scale_factor=args.scale_factor
print "Scale by a factor of", scale_factor

R_earth = 6.37101e+06

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
        phi = arctan(y/x)
        new_radius=R_earth - dist[i]*scale_factor
        new_x=new_radius*sin(theta)*cos(phi)
        new_y=new_radius*sin(theta)*sin(phi)
        new_z=new_radius*cos(theta)
        vtu_object.ugrid.GetPoints ().SetPoint (i,new_x,new_y,new_z)

    vtu_object.Write(filename=output_filename)
    print 'Scaled vtu written to', output_filename
