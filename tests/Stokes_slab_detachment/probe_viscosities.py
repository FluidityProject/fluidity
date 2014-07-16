#!/usr/bin/env python

import vtktools, numpy

def probe_visc():
    # Files and probing resolution:
    vtufile = "Benchmark_Case_1a_1.pvtu"
    npoints = 1000

    # First probe viscosity at y = 550e3
    viscosity_y_550    = numpy.zeros(npoints)

    # Setup Coordinates:
    colx               = numpy.array([numpy.linspace(0,1000e3,npoints)]).reshape(npoints,1)
    colz               = numpy.ones((npoints,1))
    coly               = numpy.ones((npoints,1))*550e3
    coordinates        = numpy.concatenate((colx,coly,colz),1) 

    # Open file and probe for Viscosity:
    vtu                = vtktools.vtu(vtufile)
    visc               = vtktools.vtu.ProbeData(vtu,coordinates,'Mantle::Viscosity')
    viscosity_y_550[:] = visc[:,0,0]

    # Next probe viscosity at x = 500e3
    viscosity_x_500    = numpy.zeros(npoints)

    # Setup Coordinates:
    coly               = numpy.array([numpy.linspace(0,660e3,npoints)]).reshape(npoints,1)
    colz               = numpy.ones((npoints,1))
    colx               = numpy.ones((npoints,1))*500e3
    coordinates        = numpy.concatenate((colx,coly,colz),1) 

    # Open file and probe for Viscosity:
    vtu                = vtktools.vtu(vtufile)
    visc               = vtktools.vtu.ProbeData(vtu,coordinates,'Mantle::Viscosity')
    viscosity_x_500[:] = visc[:,0,0]

    return min(viscosity_y_550),max(viscosity_y_550),min(viscosity_x_500),max(viscosity_x_500)








