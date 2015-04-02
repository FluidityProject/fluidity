#!/usr/bin/env python
#
# James Maddison

"""
Tool to decompose a vtu using a given decomposed gmsh mesh
"""

import glob
import optparse
import os

import fluidity.diagnostics.debug as debug
import fluidity.diagnostics.gmshtools as gmshtools
import fluidity.diagnostics.vtutools as vtktools

optionParser = optparse.OptionParser( \
  usage = "%prog [OPTIONS] ... MESH VTU", \
  add_help_option = True, \
  description = "Tool to decompose a vtu using a given decomposed gmsh mesh")

optionParser.add_option("-v", "--verbose", action = "store_true", dest = "verbose", help = "Verbose mode", default = False)
opts, args = optionParser.parse_args()

if len(args) < 2:
  debug.FatalError("GMSH base name and vtu name required")
elif len(args) > 2:
  debug.FatalError("Unrecognised trailing argument")
meshBasename = args[0]
vtuFilename = args[1]

possibleMeshBasenames = glob.glob(meshBasename + "_?*.msh")
meshBasenames = []
meshIds = []
for possibleMeshBasename in possibleMeshBasenames:
  id = possibleMeshBasename[len(meshBasename) + 1:-4]
  try:
    id = int(id)
  except ValueError:
    continue
    
  meshBasenames.append(possibleMeshBasename[:-4])
  meshIds.append(id)

vtuBasename = os.path.basename(vtuFilename[:-len(vtuFilename.split(".")[-1]) - 1])
vtuExt = vtuFilename[-len(vtuFilename.split(".")[-1]):]

vtu = vtktools.vtu(vtuFilename)
for i, meshBasename in enumerate(meshBasenames):
  debug.dprint("Processing mesh partition " + meshBasename)
  meshVtu = gmshtools.ReadMsh(meshBasename).ToVtu(includeSurface = False)
  partition = vtktools.RemappedVtu(vtu, meshVtu)
  partition.Write(vtuBasename + "_" + str(meshIds[i]) + "." + vtuExt)
