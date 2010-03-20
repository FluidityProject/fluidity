import os
from fluidity_tools import stat_parser
from sympy import *
from numpy import array,max,abs

meshtemplate='''
Point(1) = {0, 0, 0, <dx>};
Extrude {0, 1, 0} {
  Point{1};Layers{<layers>};
}
Point(3) = {1, 0, 0, <dx>};
Extrude {0, 1, 0} {
  Point{3};Layers{<layers>};
}
Line(3)={1,3};
Line(4)={2,4};
Line Loop(5) = {4, -2, -3, 1};
Plane Surface(6) = {5};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(1) = {6};
'''

def generate_meshfile(name,dx):


    file(name+".geo",'w').write(
        meshtemplate.replace('<dx>',str(dx)
                 ).replace('<layers>',str(1./dx)))

    os.system("gmsh -2 "+name+".geo")
    os.system("../../scripts/gmsh2triangle --2d "+name+".msh")


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



