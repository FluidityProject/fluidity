#!/usr/bin/env python

from optparse import OptionParser
import fluidity.diagnostics.triangletools as triangletools
import fluidity.diagnostics.gmshtools as gmshtools

def addAdditionalIds(mesh):
  # only modify the tags of the surface elements:
  for ele in mesh.GetSurfaceElements():
    ids = ele.GetIds()
    # squeezing in two zeros, first is tag is the physical ID, last one is the ElementOwner (ID):
    ele.SetIds([ids[0], 0, 0, ids[1]])

optparser=OptionParser(usage='usage: %prog <filename>',
                       add_help_option=True,
                       description="""Converts Fluidity-generated Triangle """ + 
                       """mesh files and converts them to the Gmsh format.""")

optparser.add_option("--ascii", "-a",
                  help="Convert to ASCII Gmsh format",
                     action="store_const", const=True, dest="ascii", default=False)

(options, argv) = optparser.parse_args()
filename = argv[0]

mesh = triangletools.ReadTriangle(filename)
hasInternalBoundaries = triangletools.hasPeriodicBoundary(filename)
if (hasInternalBoundaries):
  # squeeze in two additional tags before we start writing out the gmsh file:
  addAdditionalIds(mesh)
# now we can write the gmsh file:
gmshtools.WriteMsh(mesh, "%s.msh" % (filename), binary=not options.ascii)
