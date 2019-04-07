#!/usr/bin/env python

import numpy, vtktools

basename = "Subduction_Mesh"

def get_num_gmsh_nodes():
    gmshfile = basename+".msh"
    f        = open(gmshfile,'r')
    lines    = f.readlines()
    return int(lines[5])

def get_num_vtu_nodes():
    vtufile  = basename+".vtu"
    vtudata  = vtktools.vtu(vtufile)
    vtu_locs = vtudata.GetLocations()
    return int(vtu_locs.shape[0])
