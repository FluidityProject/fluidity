input: clean
	gmsh -format msh2 -1 src/circle_shell.geo

clean:
	rm -rf *.d.* *.stat *.vtu *.convergence *checkpoint* src/circle.msh fluidity.[le]* \
	matrixdump matrixdump.info
