#/bin/sh
gmsh mesh4.geo -2
../../../scripts/gmsh2triangle mesh4.msh -2
scripts/triangle_add_edgeowner.py mesh4
scripts/triangle_remove_superposed_nodes2.py mesh4_edgow 24

gmsh mesh_connected.geo -2
../../../scripts/gmsh2triangle mesh_connected.msh -2

