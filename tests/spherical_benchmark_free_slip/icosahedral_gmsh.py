#!/usr/bin/env python3

from optparse import OptionParser
import fluidity.diagnostics.gmshtools as gmshtools
from fluidity.diagnostics.elements import Element
from fluidity.diagnostics.meshes import Mesh
from polyhedra import icosahedron, subdivide
import numpy

optparser=OptionParser(usage='usage: %prog <filename> <radius> <subdivisions>',
                       add_help_option=True,
                       description="""Creates an icohedral gmsh mesh of the 2-sphere in 3D.""")
optparser.add_option("--ascii", "-a",
                  help="Convert to ASCII Gmsh format",
                     action="store_const", const=True, dest="ascii", default=False)

(options, argv) = optparser.parse_args()
try:
    filename = argv[0]
    radius = float(argv[1])
    subdivisions = int(argv[2])
except:
    optparser.error("Incorrect number or value of arguments.")

nodes, faces = icosahedron()

for i in range(subdivisions):
    subdivide(nodes, faces)

nodes = radius * numpy.array(nodes)
mesh = Mesh(3, nodes, volumeElements=[], surfaceElements=[Element(face) for face in faces])

# now we can write the gmsh file:
gmshtools.WriteMsh(mesh, "%s.msh" % (filename), binary=not options.ascii)
