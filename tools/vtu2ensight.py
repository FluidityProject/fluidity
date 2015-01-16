#!/usr/bin/env python

import sys
import argparse
import glob
import re
import vtk
import vtktools

def parse_args(argv):
    parser = argparse.ArgumentParser(
            prog="vtu2ensight",
            description="""This converts a vtu file to a ensight file. If applied to checkpointed files, use rename_checkpoint first and ensure that 'checkpoint' is removed from the basename of the solution files.""")
    parser.add_argument(
            "-v", 
            "--verbose",
            help="Print something...",
            action = "store_true",
            dest = "verbose", 
            default = False
            )
    parser.add_argument(
            "-s", 
            "--static",
            help="Use this flag only when a fixed mesh was used. By default a dynamically varying (adaptive) spatial mesh is assumed.",
            action = "store_true",
            dest = "static", 
            default = False
            )
    parser.add_argument(
            'basename', 
            metavar='basename',
            help="Basename of output (without .pvtu or .vtu)",
            )
    args = parser.parse_args()
    return args

# Function taken from:
# http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
def sorted_nicely(l):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def getvtulist(basename):
    # Find all vtu/pvtu files for (p)vtus in this folder:
    vtus = []
    for file in sorted_nicely(glob.glob(basename+'_[0-9]*vtu')):
        if (not ('checkpoint' in file)):
            vtus.append(file)
    return vtus

def getvtk(filename):
    # read in vtu file:
    reader = vtktools.vtu(filename)
    return reader

def getensightwriter(basename, static):
    writer=vtk.vtkEnSightWriter()
    writer.SetFileName(basename)
    writer.SetTransientGeometry(not(static))
    return writer

def addblockid(ug):
    # get number of elements in ug:
    nele = int(ug.GetNumberOfCells())
    # add blockID to ug (required by the ensight format)
    blockIDs = vtk.vtkUnsignedIntArray()
    blockIDs.SetNumberOfTuples(nele)
    blockIDs.SetNumberOfComponents(1)
    blockIDs.SetName("BlockId")
    for j in range(nele):
        blockIDs.SetValue(j,1)
    ug.GetCellData().AddArray(blockIDs)
    return ug

def writedata(writer, ug, i):
    #writer.SetGhostLevel(0)
    #writer.SetBlockIDs(1)
    writer.SetNumberOfBlocks(1)
    writer.SetTimeStep(i)
    writer.SetInput(ug)
    writer.Write()

def writecase(writer, ntimesteps):
    # write header information (case file)
    writer.WriteCaseFile(ntimesteps)

def main(args):
    verbose = args.verbose
    static = args.static
    basename = args.basename
    # get list of vtu/pvtu files:
    vtus = getvtulist(basename)
    if (not vtus): raise IOError
    # writer:
    writer = getensightwriter(basename, static)
    # write data for each vtu-file:
    for i in range(len(vtus)):
        if (verbose):
            print "processing vtu file: "+vtus[i]
        # get vtk object:
        reader = getvtk(vtus[i])
        # add block id (required by the ensight format):
        ug = addblockid(reader.ugrid)
        # write data:
        writedata(writer, ug, i)
    # write case file:
    writecase(writer, len(vtus))

if __name__ == "__main__":
    # get arguments:
    args = parse_args(sys.argv)
    try:
        main(args)
        print "EnSight output files have been written successfully."
    except IOError:
        print "Error: Could not find any output files with a basename \""+args.basename+"\"."
    except:
        raise Exception("Something went wrong. Aborting operation.")

