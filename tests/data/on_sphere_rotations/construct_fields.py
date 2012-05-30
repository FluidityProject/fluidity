import sys
sys.path.append('../../../python/')
import os
import vtktools
import scipy as sp
import rotation

#Generate mesh from Gmsh .geo file, convert into triangle fromat and then
# into vtu format.
os.system('gmsh -3 spherical_shell.geo')
os.system('../../../bin/gmsh2triangle spherical_shell.msh')
os.system('../../../bin/triangle2vtu spherical_shell')
#Open file containing mesh and extract vertex coordinates
file = vtktools.vtu('spherical_shell.vtu')
vertices = file.GetLocations()

#The sci-py array below contain the components of a tensor in polar coordinates,
TensorComponents_sphericalPolar= sp.array([[1.0,  0.0,  0.0 ],\
                                           [0.0,  1.0,  0.0 ],\
                                           [0.0,  0.0,  1.0]])

#The sci-py array below contain the components of three vectors in spherical-polar
#  coordinates: all are unit vectors, one aligned with the radial direction, one
#  aligned with the polar direction and one aligned in the azimuthal direction.
VectorComponents_unitRadial= sp.array([1.0, 0.0, 0.0])
VectorComponents_unitPolar= sp.array([0.0, 1.0, 0.0])
VectorComponents_unitAzimuthal= sp.array([0.0, 0.0, 1.0])

#Create lists for storing coordinate(position) vectors, tensor and vector components
# in cartesian as well as polar coordinates.
cartesianCoordinate=[]
tensor_inCartesian=[]
unitRadialVector_inCartesian=[]
unitPolarVector_inCartesian=[]
unitAzimuthalVector_inCartesian=[]
polarCoordinate=[]
tensor_inPolar=[]
unitRadialVector_inPolar=[]
unitPolarVector_inPolar=[]
unitAzimuthalVector_inPolar=[]
for point in vertices:
   #Extract cartesian coordinates from the mesh, calculate sperical-polar coorrd and
   # append to appropriate lists.
   point_sphericalPolar = rotation.cartesian_2_sphericalPolar(point)
   polarCoordinate.append(point_sphericalPolar)
   cartesianCoordinate.append(point)
   #Calculate the vector and tensor componets in a cartesian basis and append to
   # appropriate lists.
   tensor_inCartesian.append(rotation.transform_tensor_sphericalPolar_2_cartesian(point_sphericalPolar, TensorComponents_sphericalPolar))
   unitRadialVector_inCartesian.append(rotation.transform_vector_sphericalPolar_2_cartesian(point_sphericalPolar, VectorComponents_unitRadial))
   unitPolarVector_inCartesian.append(rotation.transform_vector_sphericalPolar_2_cartesian(point_sphericalPolar, VectorComponents_unitPolar))
   unitAzimuthalVector_inCartesian.append(rotation.transform_vector_sphericalPolar_2_cartesian(point_sphericalPolar, VectorComponents_unitAzimuthal))
   #Append to appropriate lists the tensor and vector components in sperical-polar basis.
   tensor_inPolar.append(TensorComponents_sphericalPolar)
   unitRadialVector_inPolar.append(VectorComponents_unitRadial)
   unitPolarVector_inPolar.append(VectorComponents_unitPolar)
   unitAzimuthalVector_inPolar.append(VectorComponents_unitAzimuthal)

file.AddVectorField('CartesianCoordinate', sp.array(cartesianCoordinate))
file.AddVectorField('PolarCoordinate', sp.array(polarCoordinate))
file.AddField('Tensor_inCartesian', sp.array(tensor_inCartesian))
file.AddVectorField('UnitRadialVector_inCartesian', sp.array(unitRadialVector_inCartesian))
file.AddVectorField('UnitPolarVector_inCartesian', sp.array(unitPolarVector_inCartesian))
file.AddVectorField('UnitAzimuthalVector_inCartesian', sp.array(unitAzimuthalVector_inCartesian))
file.AddField('Tensor_inPolar', sp.array(tensor_inPolar))
file.AddVectorField('UnitRadialVector_inPolar', sp.array(unitRadialVector_inPolar))
file.AddVectorField('UnitPolarVector_inPolar', sp.array(unitPolarVector_inPolar))
file.AddVectorField('UnitAzimuthalVector_inPolar', sp.array(unitAzimuthalVector_inPolar))

file.Write('spherical_shell_withFields.vtu')

#Clean-up, delete unwanted files
os.system('rm -f *.ele *.face *.node *.msh *.pyc spherical_shell.vtu')
