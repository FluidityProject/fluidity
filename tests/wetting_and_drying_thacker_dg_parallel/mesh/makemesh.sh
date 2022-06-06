#!/bin/sh
rm mesh_*.msh

gmsh -bin mesh_connected_A.geo -2

gmsh -bin mesh_connected_B.geo -2

gmsh -bin mesh_connected_C.geo -2

gmsh -bin mesh_connected_D.geo -2

# gmsh -bin mesh_connected_D_concentric.geo -2

# gmsh -bin mesh_connected_E.geo -2
