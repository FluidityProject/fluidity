#/bin/sh
gmsh mesh_connected.geo -2
../../../bin/gmsh2triangle mesh_connected.msh -2

