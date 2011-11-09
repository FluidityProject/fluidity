#/bin/sh
../../../bin/gmsh2triangle mesh_connected.msh -2
../../../bin/fldecomp -n 4 -m triangle mesh_connected

