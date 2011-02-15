#A module for converting between spherical and polar coordinates

#Usage 

def polar2cart(r):
    #convert polar coordinates (r=1) to Cartesian coordinates
    from math import sin,cos,pi
    from numpy import array
    def val(X,t):
        theta = 2*pi*X[0]/360
        phi = 2*pi*X[1]/360
        return r*array([cos(theta)*cos(phi),sin(theta)*cos(phi),sin(phi)])
    return val

def cart2polar():
    #convert Cartesian coordinates to polar coordinates with r=1
    from math import asin,atan2,pi,sqrt
    from numpy import array
    def val(X,t):
        x = X[0]
        y = X[1]
        z = X[2]
        r = sqrt(x**2 + y**2 + z**2)
        phi = asin(z/r)
        theta = atan2(y,x)
        return 360/(2*pi)*array([theta,phi])
    return val

def spherical_basis_vecs():
    #return unit vectors in theta, phi
    from math import sin, cos, pi
    from numpy import array
    def val(X,t):
        theta = 2*pi*X[0]/360
        phi = 2*pi*X[1]/360
        return array([[-sin(theta), cos(theta), 0.0],[-cos(theta)*sin(phi), -sin(theta)*sin(phi), cos(phi)]])
    return val
