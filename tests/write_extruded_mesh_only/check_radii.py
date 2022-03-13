from __future__ import print_function
import numpy, glob, sys, os, vtk, constants

def check_min_max_radii():
  # Loop over vtu's and check that each core does not extend from inner radius to outer radius:
  filelist = glob.glob("Stokes_sphere_0/Stokes_sphere_0_*.vtu")
  filelist.sort()

  smallest_radius_max = constants.outer_radius
  largest_radius_min = constants.inner_radius

  for filename in filelist:
    print('Working on file:', filename)

    # Set up reader (note that vtktools does this automatically, but vtktools isn't being used here):
    if filename[-4:] == ".vtu":
      reader=vtk.vtkXMLUnstructuredGridReader()
    else:
      raise Exception("ERROR: don't recognise file extension" + filename)
    reader.SetFileName(filename)
    reader.Update()

    # Extract and store all field data:
    data = reader.GetOutput()

    # Extract radius information:
    radius_data = data.GetPointData().GetScalars("Radius")

    # Convert radii into a nice numpy array:
    radius = numpy.array([radius_data.GetValue(p) for p in range(radius_data.GetNumberOfTuples())])
    
    # Check for smallest maximum radius and largest minimum radius on each core:
    smallest_radius_max = min(smallest_radius_max, numpy.max(radius))
    largest_radius_min = max(largest_radius_min, numpy.min(radius))

  return smallest_radius_max, largest_radius_min
