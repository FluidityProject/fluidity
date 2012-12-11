"""Module for computing Galewsky initial condition"""
from scipy import integrate
from scipy import arcsin,pi

def u(phi):
    """Calculate the West-East velocity as a function of phi (latitude)"""
    from numpy import exp, pi
    umax = 80.0
    phi0 = pi/7.0
    phi1 = pi/2.0-phi0
    en = exp(-4.0/(phi1-phi0)**2)
    uout = umax/en*exp(1./(phi-phi0)/(phi-phi1))
    uout[phi>=phi1] = 0.
    uout[phi<=phi0] = 0.
    return uout

def Uval(X,t):
    """Interface function for calculating U"""
    from numpy import array,arcsin
    x = X[0]
    y = X[1]
    z = X[2]
    R = (x**2+y**2+z**2)**0.5
    R0 = 6.37122e6
    z = R0*z/R
    phi = arcsin(z/R)
    uin = u(array([phi]))
    uin = uin[0]
    Rz = (x**2+y**2)**0.5
    erx = -y/Rz
    ery = x/Rz
    return [erx*uin,ery*uin,0.]

def h_integrand(phi):
    from numpy import sin,tan
    R0 = 6.37122e6
    Omega = 7.292e-5
    f = 2*Omega*sin(phi)
    return -u(phi)*(f+tan(phi)*u(phi)/R0)

def h_perturbation(x,y,z):
    from numpy import arcsin,arctan2,cos,exp,pi
    alpha = 1./3.
    beta = 1./15.
    hhat = 120
    R = (x**2+y**2+z**2)**0.5
    phi = arcsin(z/R)
    phi2 = pi/4
    lambda0 = arctan2(y,x)
    return hhat*cos(phi)*exp(-(lambda0/alpha)**2)*\
        exp(-((phi2-phi)/beta)**2)

def Hval(X,t):
    """Interface function for calculating h"""
    x = X[0]
    y = X[1]
    z = X[2]
    R = (x**2+y**2+z**2)**0.5
    R0 = 6.37122e6
    g = 9.80616
    phi = arcsin(z/R)
    phi0 = pi/7.0
    phi1 = pi/2.0-phi0
    if(phi<=phi0): 
        return 0.
    else:
        val,err = integrate.quadrature(h_integrand,phi0,phi,tol=1.0e-8)
        return val*R0/g

def HPval(X,t):
    """Interface function for calculating perturbed h"""
    x = X[0]
    y = X[1]
    z = X[2]
    R = (x**2+y**2+z**2)**0.5
    R0 = 6.37122e6
    g = 9.80616
    phi = arcsin(z/R)
    phi0 = pi/7.0
    phi1 = pi/2.0-phi0
    if(phi<=phi0): 
        val = 0.
    else:
        val,err = integrate.quadrature(h_integrand,phi0,phi,tol=1.0e-8)
        val = val*R0/g
    return val + h_perturbation(x,y,z)

if __name__ == "__main__":
    import numpy as np
    import pylab as pl
    phi = np.arange(20./90.*np.pi/2,70./90.*np.pi/2,50./90*np.pi/2/100)
    R0 = 6.37122e6
    z = R0*np.sin(phi)
    x = R0*np.cos(phi)
    y = 0.*z
    u0 = 0.*z
    h0 = 0.*z
    for i, uval in enumerate(u0):
        v = Uval([x[i],y[i],z[i]],0.)
        u0[i] = u(np.array([phi[i]]))#v[0]
        h0[i] = Hval([x[i],y[i],z[i]],0.)
    pl.close('all')
    phi = phi/np.pi*180
    pl.plot(u0,phi)
    pl.axis([u0.min(),u0.max(),phi.min(),phi.max()])
    pl.figure()
    pl.plot(h0,phi)
    pl.axis([h0.min(),h0.max(),phi.min(),phi.max()])
    pl.show()
