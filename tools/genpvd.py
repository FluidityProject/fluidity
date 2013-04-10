#!/usr/bin/env python

# This script generates pvd files for time series in paraview
import vtktools as vtk
import glob
import re
import sys


def generate_pvd(pvdfilename, list_vtu_filenames, list_time):
    """ This function writes a simple xml pvd file that can be
        loaded in paraview to display time series correctly
        with irregular timesteps
    """
    if (len(list_vtu_filenames) != len(list_time)):
        print "Error, list of filenames and time are of unequal length. Exiting!"
        exit()
    pvdfile = open(pvdfilename, "w")
    # Write header:
    pvdfile.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n')
    # Now write the information for the time series:
    for i in range(len(list_time)):
        pvdfile.write('    <DataSet timestep="'+str(list_time[i])+'" group="" part="0" file="'+list_vtu_filenames[i]+'"/>\n')
    # Write closing statements in the file and close it:
    pvdfile.write('</Collection>\n</VTKFile>')
    pvdfile.close()

# Function taken from:
# http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
def sorted_nicely(l):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# Get basename of the simulation from command line or assemble it on your own:
try:
    simulation_basename = sys.argv[1]
except:
    print "ERROR: You have to give genpvd the basename of the considered vtu files."
    print "Program will exit..."
    exit()

# Find all vtu/pvtu files for fluid vtus in this folder:
fluid_vtus = []
for file in sorted_nicely(glob.glob(simulation_basename+'_[0-9]*vtu')):
    if (not ('checkpoint' in file)):
        fluid_vtus.append(file)

# Loop over all the fluid vtus found and collect time data from them to assemble the pvd file:
time = []
for filename in fluid_vtus:
    print "Processing file: ", filename
    # Get the unstructured mesh:
    data = vtk.vtu(filename)
    # Only process the first node in the mesh, as the time is constant over the whole mesh:
    n0 = data.ugrid.GetCell(0).GetPointId(0)
    t = data.ugrid.GetPointData().GetArray("Time").GetTuple(n0)
    time.append(t[0])

# Generate fluid pvd file:
generate_pvd(simulation_basename+'.pvd', fluid_vtus, time)

print "============================================="
print "Program exited without errors."
print "Open the file"
print "  "+simulation_basename+".pvd"
print "in Paraview."
print "============================================="
