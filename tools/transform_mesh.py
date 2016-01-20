#!/usr/bin/env python

from optparse import OptionParser
import sys
import math
import fluidity.diagnostics.gmshtools as gmshtools

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

# make all definitions of the math module available in
# the transformation expression.
globals = math.__dict__

mesh = gmshtools.ReadMsh(mesh_name+'.msh')
dim = mesh.GetDim()

def remap(coords):
    globals['x'] = coords[0]
    globals['y'] = coords[1]
    if dim == 3:
        globals['z'] = coords[2]
    return eval(transformation, globals)

mesh.RemapNodeCoords(remap)
gmshtools.WriteMsh(mesh, mesh_name+'.msh')
