#gmsh -bin MonaiValley_B.geo -2
#fldecomp -n 4 -m gmsh MonaiValley_B

gmsh -bin MonaiValley_C.geo -2
../../../bin/fldecomp -n 4 -m gmsh MonaiValley_C

#gmsh -bin MonaiValley_D.geo -2
#fldecomp -n 4 -m gmsh MonaiValley_D


