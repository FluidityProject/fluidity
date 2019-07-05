#!/usr/bin/env python

from __future__ import print_function

from optparse import OptionParser
import sys
import shutil
import math
import meshio
import numpy


#####################################################################
# Script starts here.
optparser=OptionParser(usage='usage: %prog region transformation mesh',
                       add_help_option=True,
                       description="""Applies a coordinate transformation to a region of a given mesh.""")

optparser.set_usage(
                  "Usage:\n"+
								  "  %prog (<constants>) (<region>) <transformation> <mesh>\n"+
                  "<constants>       is a list of constant associations separated by commas, for use in region and transformation (optional, although when specified, <region> is required - use 'True' for whole domain (i.e. the default behaviour)).\n"+
                  "<region>          is a python expression which evaluates to true over the region to be transformed (optional, default is True - i.e. the whole domain).\n"+
                  "<transformation>  is a python expression giving the coordinate transformation.\n"+
                  "<mesh>            is the name of the gmsh mesh file.\n"+
                  "Note: Creates a backup of the original mesh file with a '.bak' extension\n"
                  "\n"+
                  "Examples\n"
                  "- To rescale the z-dimension by a factor of 1000.\n"
                  "  %prog '(x,y,1000*z)' mesh.msh\n"
                  "- To project all points that lie within a circle of centre (xcentre,ycentre) in z by a distance zprojection.\n"
                  "  %prog 'xcentre=50, ycentre=50, radius=20, zprojection=50' '(x-xcentre)**2 + (y-ycentre)**2 < radius**2' '(x, y, z+zprojection)' mesh.msh\n"
                  "- To project all points that lie within a circle of centre (xcentre,ycentre) in z in the shape of a cone, by a distance zprojection at the centre.\n"
                  "  %prog 'xcentre=50, ycentre=50, radius=20, zprojection=50' '(x-xcentre)**2 + (y-ycentre)**2 < radius**2' '(x, y, z + zprojection * (1 - sqrt((x-xcentre)**2 + (y-ycentre)**2) / radius ) )' mesh.msh\n"
                  "- To add an ice shelf to a meshed box for x in [0,shelflength] and z in [0,shelfslopeheight + minoceandepth].  Note this applies to both 2d and 3d domains and the ocean domain can extend further.\n"
                  "  %prog 'shelflength = 550, shelfslopeheight = 800, minoceandepth = 100' 'x < shelflength' '(x, y, (z/(shelfslopeheight + minoceandepth)) * ((x/shelflength) * shelfslopeheight + minoceandepth))' mesh.msh\n"
         )

(options, argv) = optparser.parse_args()

# make all definitions of the math module available in
# the transformation expression.
env = math.__dict__

if len(argv) == 4:
    constants      = argv[0]
    region         = argv[1]
    transformation = argv[2]
    mesh_name      = argv[3]
    constants = constants.replace(' ','')
    constanttextarray = constants.split(',')
    for constanttext in constanttextarray:
      constant = constanttext.split('=')
      env[constant[0]] = float(constant[1])

elif len(argv) == 3:
    region         = argv[0]
    transformation = argv[1]
    mesh_name      = argv[2]
elif len(argv) == 2:
    transformation = argv[0]
    mesh_name      = argv[1]
else:
    optparser.print_help()
    sys.exit(1)

mesh = meshio.read(mesh_name)

def remap(coords):
    local_env = dict(zip(['x', 'y', 'z'], coords), **env)
    if eval(region, local_env):
        return eval(transformation, local_env)
    else:
        return coords

mesh.points = numpy.apply_along_axis(remap, 1, mesh.points)
meshio.write(mesh_name, mesh)
