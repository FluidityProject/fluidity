input: clean
	gmsh -2 src/column.geo
	../../scripts/gmsh2triangle --2d src/column.msh
	../../bin/fldecomp -n 2 -i triangle column

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.edge *.ele *checkpoint* *.log-0 *.err-0 ./src/*.msh *.pvtu
	rm -rf matrixdump*
