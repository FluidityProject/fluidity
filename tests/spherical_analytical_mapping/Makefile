input: clean
	gmsh -2 src/circle.geo
	gmsh -3 src/hemisphere.geo

clean:
	rm -rf *.d.* *.stat *.vtu *.convergence *checkpoint* src/circle.msh src/hemisphere.msh fluidity.[le]* \
	matrixdump matrixdump.info
