#!/usr/bin/env python

import vtktools, numpy

def probe_visc():
    # Files and probing resolution:
    vtufile = "Benchmark_Case_1a_1.pvtu"
    npoints = 500

    # First probe viscosity at y = 550e3

    # Setup Coordinates:
    colx               = numpy.array([numpy.linspace(0,500e3,npoints)]).reshape(npoints,1)
    colz               = numpy.zeros((npoints,1))
    coly               = numpy.ones((npoints,1))*550e3
    coordinates        = numpy.concatenate((colx,coly,colz),1) 

    # Open file and probe for Viscosity:
    vtu                = vtktools.vtu(vtufile)
    viscosity_y_550    = vtktools.vtu.ProbeData(vtu,coordinates,'Mantle::Viscosity')[:,0,0]

    # Next probe viscosity at x = 500e3

    # Setup Coordinates:
    coly               = numpy.array([numpy.linspace(0,660e3,npoints)]).reshape(npoints,1)
    colz               = numpy.zeros((npoints,1))
    colx               = numpy.ones((npoints,1))
    coordinates        = numpy.concatenate((colx,coly,colz),1) 

    # Open file and probe for Viscosity:
    vtu                = vtktools.vtu(vtufile)
    viscosity_x_500    = vtktools.vtu.ProbeData(vtu,coordinates,'Mantle::Viscosity')[:,0,0]

    return min(viscosity_y_550), max(viscosity_y_550), min(viscosity_x_500), max(viscosity_x_500)








