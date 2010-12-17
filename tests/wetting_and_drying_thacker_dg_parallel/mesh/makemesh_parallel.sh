#/bin/sh
rm mesh_*.msh

#gmsh -bin mesh_connected_A.geo -2
#fldecomp -n 4 -m gmsh mesh_connected_A

gmsh -bin mesh_connected_B.geo -2
fldecomp -n 4 -m gmsh mesh_connected_B

#gmsh -bin mesh_connected_C.geo -2
#fldecomp -n 4 -m gmsh mesh_connected_C

#gmsh -bin mesh_connected_D.geo -2
#fldecomp -n 4 -m gmsh mesh_connected_D

#gmsh -bin mesh_connected_D_concentric.geo -2
#fldecomp -n 4 -m gmsh mesh_connected_D_concentric

#gmsh -bin mesh_connected_E.geo -2
#fldecomp -n 4 -m gmsh mesh_connected_E

