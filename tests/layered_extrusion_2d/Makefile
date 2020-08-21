input: clean
	gmsh -1 src/circle_shell.geo -o circle_shell.msh

clean:
	rm -rf *.stat *.vtu *.convergence *checkpoint* src/circle.msh fluidity.[le]* \
	matrixdump matrixdump.info *.msh
