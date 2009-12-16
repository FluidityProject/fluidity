input: clean
	gmsh -2 src/column.geo
	../../scripts/gmsh2triangle --2d src/column.msh

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.edge *.ele *checkpoint* *.log-0 *.err-0 ./src/*.msh
	rm -rf matrixdump matrixdump.info
