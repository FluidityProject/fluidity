input: clean
	gmsh -format msh2 -2 src/hemisphere_shell.geo

clean:
	rm -rf *.d.* *.stat *.vtu *.convergence *checkpoint* src/hemisphere_shell.msh fluidity.[le]* \
	matrixdump matrixdump.info
