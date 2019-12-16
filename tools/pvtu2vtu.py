#!/usr/bin/env python3
#
# James Maddison

"""
Script to combine pvtus into vtus
"""

import optparse

import fluidity.diagnostics.debug as debug
import fluidity.diagnostics.fluiditytools as fluidity_tools
import fluidity.diagnostics.vtutools as vtktools

optionParser = optparse.OptionParser( \
  usage = "%prog [OPTIONS] ... PROJECT [FIRSTID] [LASTID]", \
  add_help_option = True, \
  description = "Combines pvtus into vtus")

optionParser.add_option("-v", "--verbose", action = "store_true", dest = "verbose", help = "Verbose mode", default = False)

opts, args = optionParser.parse_args()

if not opts.verbose:
  debug.SetDebugLevel(0)
  
if len(args) > 3:
  debug.FatalError("Unrecognised trailing argument")
inputProject = args[0]

if len(args) == 2:
  try:
    firstId = int(args[1])
    lastId = firstId
  except ValueError:
    debug.FatalError("Invalid first dump ID")

if len(args) == 3:
  try:
    firstId = int(args[1])
    lastId = int(args[2])
    assert(lastId >= firstId)
  except:
    debug.FatalError("Invalid last dump ID")
 
if len(args) == 1:
  filenames = [inputProject + ".pvtu"]
else:
  filenames = fluidity_tools.VtuFilenames(inputProject, firstId, lastId = lastId, extension = ".pvtu")

for filename in filenames:
  debug.dprint("Processing file: " + filename)

  vtu = vtktools.vtu(filename)
  vtu = vtktools.VtuFromPvtu(vtu)
  vtu.Write(filename[:-5] + ".vtu")
