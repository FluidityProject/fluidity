input: clean
	gmsh -2 src/column.geo
	../../bin/gmsh2triangle --2d src/column.msh
	../../bin/fldecomp -n 2 column

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.edge *.ele *checkpoint* *.log-0 *.err-0 ./src/*.msh
	rm -rf matrixdump* *.pvtu
