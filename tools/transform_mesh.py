#!/usr/bin/env python3

from optparse import OptionParser
import sys
import math
import meshio
import numpy

#####################################################################
# Script starts here.
optparser=OptionParser(usage='usage: %prog transformation mesh',
                       add_help_option=True,
                       description="""Applies a coordinate transformation to the given mesh.""")

optparser.set_usage(
                  "usage: %prog <transformation> <mesh>\n\n"+
                  "<transformation> is a python expression giving the coordinate transformation.\n"+
                  "<mesh> is the name of the gmsh mesh file. You need a mesh.msh file.\n"+
                  "\n"+
                  "Example:\n"
                  "To rescale the z-dimension by a factor of 1000,\n"
                  "%prog '(x,y,1000*z)' mesh.\n"
		     )

(options, argv) = optparser.parse_args()

if len(argv) != 2:
    optparser.print_help()
    sys.exit(1)

transformation = argv[0]
mesh_name = argv[1]

mesh = meshio.read(mesh_name + '.msh')
def remap(coords):
    dims = ['x', 'y'] if len(coords) == 2 else ['x', 'y', 'z']
    return eval(transformation, dict(zip(dims, coords), **math.__dict__))

mesh.points = numpy.apply_along_axis(remap, 1, mesh.points)

meshio.write(mesh_name + '.msh', mesh)
