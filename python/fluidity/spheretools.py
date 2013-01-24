#A module for converting between spherical and polar coordinates

import numpy as np
#Usage 

def polar2cart(r):
    #convert polar coordinates to Cartesian coordinates
    from math import sin,cos,pi
    from numpy import array
    def val(X,t):
        theta = 2*pi*X[0]/360
        phi = 2*pi*X[1]/360
        return r*array([cos(theta)*cos(phi),sin(theta)*cos(phi),sin(phi)])
    return val

def cart2polar(X):
    #convert Cartesian coordinates to polar coordinates
    from math import asin,atan2,pi,sqrt
    from numpy import array
    x = X[0]
    y = X[1]
    z = X[2]
    r = sqrt(x**2 + y**2 + z**2)
    phi = asin(z/r)
    theta = atan2(y,x)
    return array([theta,phi])

def spherical_basis_vecs():
    #return unit vectors in theta, phi
    from math import sin, cos, pi
    from numpy import array
    def val(X,t):
        theta = 2*pi*X[0]/360
        phi = 2*pi*X[1]/360
        return array([[-sin(theta), cos(theta), 0.0],[-cos(theta)*sin(phi), -sin(theta)*sin(phi), cos(phi)]]).T
    return val

def spherical_down():
    #return down
    from math import sin, cos, pi
    from numpy import array
    def val(X,t):
        ([theta,phi]) = cart2polar(X)
        return -array([cos(theta)*cos(phi), sin(theta)*cos(phi), sin(phi)])
    return val

def coriolis(omega):
    #return f=2 * Omega sin(phi)
    from math import sin, pi
    def val(X,t):
        phi = 2*pi*X[1]/360
        return 2*omega*sin(phi)
    return val
def vector_cartesian_2_spherical_polar(xComp, yComp, zComp,
        xCoord, yCoord, zCoord):
    # Subroutine for vector change of basis: from Cartesian to spherical-polar. The
    # coordinates of the position vector are also transformed
    #real, intent(in) :: xComp         !1st vector component in cartesian basis
    #real, intent(in) :: yComp         !2nd vector component in cartesian basis
    #real, intent(in) :: zComp         !3rd vector component in cartesian basis
    #real, intent(in) :: xCoord        !1st vector component of position vector in cartesian basis
    #real, intent(in) :: yCoord        !2nd vector component of position vector in cartesian basis
    #real, intent(in) :: zCoord        !3rd vector component of position vector in cartesian basis
    #real, intent(out) :: radial       !Radial component of vector
    #real, intent(out) :: polar        !Polar component of vector
    #real, intent(out) :: azimuthal    !Azimuthal  component of vector
    #real, intent(out) :: radius       !Distance from centre of sphere
    #real, intent(out) :: theta        !Polar angle, in radians
    #real, intent(out) :: phi          !Azimuthal angle, in radians
    #real, dimension(3) :: sphericalPolarComponents
    #real, dimension(3,3) :: R   !Transformation matrix
    #real, dimension(3,3) :: RT  !Transposed transformation matrix

    # Calculate transformation matrix
    R,[theta,phi] = transformation_matrix_cartesian_2_spherical_polar(xCoord, yCoord, zCoord)

    # Evaluate vector components in spherical-polar basis
    sphericalPolarComponents = R.dot(np.array([xComp, yComp, zComp]))
    radial = sphericalPolarComponents[0,0]
    polar = sphericalPolarComponents[0,1]
    azimuthal = sphericalPolarComponents[0,2]

    return [[radial, polar, azimuthal],[theta,phi]]

def transformation_matrix_cartesian_2_spherical_polar(xCoord, yCoord, zCoord):

    from math import sin, cos
    # Subroutine calculating transformation matrix for spherical-polar to/from Cartesian
    # tensor transformations. The routine also returns the transposed transformation matrix
    #real, intent(in) :: xCoord  !x-component of position vector
    #real, intent(in) :: yCoord  !y-component of position vector
    #real, intent(in) :: zCoord  !z-component of position vector
    #real, dimension(3,3), intent(out) :: R   !Transformation matrix
    #real, dimension(3,3), intent(out) :: RT  !Transposed transformation matrix

    #real :: radius        !Distance from centre of sphere
    #real :: theta         !Polar angle, in radians
    #real :: phi           !Azimuthal angle, in radians

    # Calculate position-vector components in spherical-polar basis
    [theta,phi] = cart2polar([xCoord, yCoord, zCoord])

    R = np.zeros((3,3))

    R[0,0]=sin(theta)*cos(phi)
    R[0,1]=sin(theta)*sin(phi)
    R[0,2]=cos(theta)
    R[1,0]=cos(theta)*cos(phi)
    R[1,1]=cos(theta)*sin(phi)
    R[1,2]=-sin(theta)
    R[2,0]=-sin(phi)
    R[2,1]=cos(phi)
    R[2,2]=0.0

    R = np.matrix(R)

    RT = np.transpose(R)

    return [RT, [theta,phi]]
