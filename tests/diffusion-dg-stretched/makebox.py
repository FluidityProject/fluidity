import os
from fluidity_tools import stat_parser
from sympy import *
from numpy import array,max,abs

meshtemplate='''
Point(1) = {0.0,0.0,0,0.1};
Extrude {<width>,0,0} {
  Point{1}; Layers{<layers>};
}
Extrude {0,1,0} {
  Line{1}; Layers{<layers>};
}
Physical Line(6) = {3};
Physical Line(7) = {1};
Physical Line(8) = {4};
Physical Line(9) = {2};
Physical Surface(10) = {5};
'''

def generate_meshfile(name,layers,width):

    geo = meshtemplate.replace('<layers>',str(layers))
    geo = geo.replace('<width>',str(width))
    open(name+".geo",'w').write(geo)

    os.system("gmsh -2 "+name+".geo")
    os.system("../../bin/gmsh2triangle -2 "+name+".msh")


