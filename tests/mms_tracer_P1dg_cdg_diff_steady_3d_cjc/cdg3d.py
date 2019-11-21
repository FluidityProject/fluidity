import os
from fluidity_tools import stat_parser
from sympy import *
from numpy import array,max,abs

meshtemplate='''
Point(1) = {0.0,0.0,0,0.1};
Extrude {1,0,0} {
  Point{1}; Layers{<layers>};
}
Extrude {0,1,0} {
  Line{1}; Layers{<layers>};
}
Extrude {0,0,1} {
  Surface{5}; Layers{<layers>};
}
Physical Surface(28) = {5,14,26,22,27,18};
Physical Volume(29) = {1};
'''

def generate_meshfile(name,layers):

    geo = meshtemplate.replace('<layers>',str(layers))
    open(name+".geo",'w').write(geo)

    os.system("gmsh -3 "+name+".geo")
    os.system("../../bin/gmsh2triangle "+name+".msh")


def run_test(layers, binary):
    '''run_test(layers, binary)

    Run a single test of the channel problem. Layers is the number of mesh
    points in the cross-channel direction. The mesh is unstructured and
    isotropic. binary is a string containing the fluidity command to run.
    The return value is the error in u and p at the end of the simulation.'''

    generate_meshfile("channel",layers)

    os.system(binary+" channel_viscous.flml")

    s=stat_parser("channel-flow-dg.stat")
    return (s["Water"]['AnalyticUVelocitySolutionError']['l2norm'][-1],
        s["Water"]['AnalyticPressureSolutionError']['l2norm'][-1])
