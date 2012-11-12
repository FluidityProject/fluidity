#    Copyright (C) 2012 Imperial College London and others.
#    
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    amcgsoftware@imperial.ac.uk
#    
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA
from math import sqrt, pi, sin, cos, atan2, acos

def cartesian_2_sphericalPolar(positionVectorCartesian):
    '''Convert Cartesian coordinates to radial-azimuthal-polar spherical coordinates, in radians.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Extract the Cartesian vector components.
    x = positionVectorCartesian[0]
    y = positionVectorCartesian[1]
    z = positionVectorCartesian[2]
    #Calculate the radius.
    radius = sqrt(x**2 + y**2 + z**2)
    #Calculate azimuthal (theta) and zenith (phi) angles
    theta = acos(z/radius)
    phi = atan2(y,x)
    return [radius, theta, phi]

def cartesian_2_lonlatradius(positionVectorCartesian):
    '''Convert Cartesian coordinates on a sphere to longitude-latitude-radius. Longitude and latitude are returned in degrees.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Calculate azimuthal (phi), polar (theta) angles and distance from origin.
    [radius, theta, phi] = cartesian_2_sphericalPolar(positionVectorCartesian)
    #Calculate longitude and latitude
    lon = phi*180.0/pi
    lat = (pi/2 - theta)*180.0/pi
    
    positionVectorLonlat = [lon, lat, radius]
    return positionVectorLonlat

def lonlatradius_2_cartesian(positionVectorLonLatRad):
    '''Convert longitude-latitude-radial coordinates on the surface of the Earth (in degrees) to Cartesian coordinates. Longitude and latitude must be in degrees.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #Calculate spherical-polar coordinates form longitude-latitude-radius.
    [radius, theta, phi] = lonlatradius_2_sphericalPolar(positionVectorLonLatRad)
    #Calculate Cartesian coordinates from spherical-polar coordinates.
    x = radius*np.sin(theta)*np.cos(phi)
    y = radius*np.sin(theta)*np.sin(phi)
    y = radius*np.cos(theta)
    return [x, y, z]

def lonlatradius_2_sphericalPolar(positionVectorLonLatRad):
    '''Convert longitude-latitude-radial coordinates on the surface of the Earth (in degrees) to Cartesian coordinates. Longitude and latitude must be in degrees, the azimuthal and polar angles are returned in radians.'''
    [longitude, latitude, radius] = positionVectorLonLatRad
    #Calculate azimuthal (phi), polar (theta) angles.
    phi = longitude*pi/180.0
    theta = pi/2 - latitude*pi/180.0
    return [radius, theta, phi]

def transform_tensor_sphericalPolar_2_cartesian(positionVectorSpherical, tensor):
    '''Function changing the basis of a tensor from zonal-meridional-radial basis to a Cartesian basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    from numpy import linalg
    #extract distance from origin, polar (theta) angles and azimuthal (phi) angles.
    [radius, theta, phi] = positionVectorSpherical
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.cos(theta)*np.cos(phi),  -np.sin(phi)],\
                 [np.sin(theta)*np.sin(phi),   np.cos(theta)*np.sin(phi),   np.cos(phi)],\
                 [np.cos(theta),              -np.sin(theta),              0]])
    transposedTransformationMatrix = transformationMatrix.transpose()
    #Calculate the components of the tensor in the reference system.
    transformed_Tensor = np.dot(transformationMatrix, np.array(tensor))
    transformed_Tensor = np.dot(transformed_Tensor,transposedTransformationMatrix)
    return transformed_Tensor

def transform_tensor_cartesian_2_sphericalPolar(positionVectorCartesian, tensor):
    '''Function transforming the components of a tensor from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    from numpy import linalg
    #Calculate azimuthal (theta) and zenith (phi) angles and distance from origin
    [radius, theta, phi] = cartesian_2_sphericalPolar(positionVectorCartesian)
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.sin(theta)*np.sin(phi),   np.cos(theta)],\
                 [np.cos(theta)*np.cos(phi),   np.cos(theta)*np.sin(phi),  -np.sin(theta)],\
                 [-np.sin(phi),                np.cos(phi),                 0]])
    transposedTransformationMatrix = transformationMatrix.transpose()
    #Calculate the components of the tensor in the reference system.
    transformed_Tensor = np.dot(transformationMatrix, np.array(tensor))
    transformed_Tensor = np.dot(transformed_Tensor,transposedTransformationMatrix)
    return transformed_Tensor

def transform_tensor_sphericalPolar_2_lon_lat_rad(tensor):
    '''Function transforming the components of a tensor from a spherical-polar basis to a zonal-meridional-radial basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Tensor = np.dot(np.dot(transformationMatrix, np.array(tensor)), transformationMatrix)
    return transformed_Tensor

def transform_tensor_lon_lat_rad_2_sphericalPolar(tensor):
    '''Function transforming the components of a tensor from a zonal-meridional-radial basis to a spherical-polar basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Tensor = np.dot(np.dot(transformationMatrix, np.array(tensor)), transformationMatrix)
    return transformed_Tensor

def transform_tensor_cartesian_2_lon_lat_rad(positionVectorCartesian, tensor):
    '''Function transforming the components of a tensor from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #Transform from Cartesian into spherical-polar
    transformed_Tensor = transform_tensor_cartesian_2_sphericalPolar(positionVectorCartesian, tensor)
    #Transform from spherical-polar into longitude-latitude-radius.
    transformed_Tensor = transform_tensor_sphericalPolar_2_lon_lat_rad(transformed_Tensor)
    return transformed_Tensor

def transform_tensor_lon_lat_rad_2_cartesian(positionVectorLonLatRad, tensor):
    '''Function transforming the components of a tensor from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #Transform coordinates from lon-lat-rad into spherical-polar.
    positionVectorSpericalPolar = lonlatradius_2_sphericalPolar(positionVectorLonLatRad)
    #Transform tensor from lon-lat-rad into spherical-polar.
    transformed_Tensor = transform_tensor_lon_lat_rad_2_sphericalPolar(tensor)
    #Transform spherical-polar into Cartesian.
    transformed_Tensor = transform_tensor_sphericalPolar_2_cartesian(positionVectorSpericalPolar, transformed_Tensor)
    return transformed_Tensor

def transform_vector_sphericalPolar_2_cartesian(positionVectorSpherical, vector):
    '''Function transforming the components of a vector from a spherical-polar basis to a Cartesian basis. The input position vector must be given as [radius, polar angle, azimuthal angle], all angles specified in radians.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #extract distance from origin, polar (theta) angles and azimuthal (phi) angles.
    [radius, theta, phi] = positionVectorSpherical
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.cos(theta)*np.cos(phi),  -np.sin(phi)],\
                 [np.sin(theta)*np.sin(phi),   np.cos(theta)*np.sin(phi),   np.cos(phi)],\
                 [np.cos(theta),              -np.sin(theta),              0]])
    #Calculate the components of the tensor in the reference system.
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_cartesian_2_sphericalPolar(positionVectorCartesian, vector):
    '''Function transforming the components of a vector from a Cartesian basis to a spherical-polar basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    import numpy as np
    #Calculate distance from origin, polar (theta) angles and azimuthal (phi) angles.
    [radius, theta, phi] = cartesian_2_sphericalPolar(positionVectorCartesian)
    #Evaluate components of rotation matrices.
    transformationMatrix =\
       np.array([\
                 [np.sin(theta)*np.cos(phi),   np.sin(theta)*np.sin(phi),   np.cos(theta)],\
                 [np.cos(theta)*np.cos(phi),   np.cos(theta)*np.sin(phi),  -np.sin(theta)],\
                 [-np.sin(phi),                np.cos(phi),                 0]])
    #Calculate the components of the tensor in the reference system.
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_sphericalPolar_2_lon_lat_rad(vector):
    '''Function transforming the components of a vector from a spherical-polar basis to a zonal-meridional-radial basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_lon_lat_rad_2_sphericalPolar(vector):
    '''Function transforming the components of a vector from a zonal-meridional-radial basis to a spherical-polar basis.'''
    import numpy as np
    #The requested transformation is just a reflection followed by a change in the order
    # of the components in order to get a right-handed system.
    transformationMatrix = np.array([[ 0.0, 0.0, 1.0 ],\
                                     [ 0.0,-1.0, 0.0 ],\
                                     [ 1.0, 0.0, 0.0 ]])
    transformed_Vector = np.dot(transformationMatrix, np.array(vector))
    return transformed_Vector

def transform_vector_cartesian_2_lon_lat_rad(positionVectorCartesian, vector):
    '''Function transforming the components of a vector from a Cartesian basis to a zonal-meridional-radial basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Calculate spherical-polar components of the vector.
    transformed_Vector = transform_vector_cartesian_2_sphericalPolar(positionVectorCartesian, vector)
    #Calculate zonal, meridional and radial components of the vector.
    transformed_Vector = transform_vector_sphericalPolar_2_lon_lat_rad(transformed_Vector)
    return transformed_Vector

def transform_vector_lon_lat_rad_2_cartesian(positionVectorRadLonLat, vector):
    '''Function transforming the components of a vector from a zonal-meridional-radial basis into a Cartesian basis.

The origin of the Cartesian frame of reference is located at the centre of the sphere, the positive half of x-axis goes through 0 deg E, 0 deg N, the positive half of y-axis goes through 90 deg E, 0 deg N and the positive half of the z-axis goes through the North Pole equivalent.'''
    #Transform coordinates from longitude-latitude-radius into spherical-polar
    positionVectorSpherical = lonlatradius_2_sphericalPolar(positionVectorRadLonLat)
    #Transform vector from longitude-latitude-radius into spherical-polar
    transformed_Vector = transform_vector_lon_lat_rad_2_sphericalPolar(vector)
    #Transform from spherical-polar into Cartesian
    transformed_Vector = transform_vector_sphericalPolar_2_cartesian(positionVectorSpherical, vector)
