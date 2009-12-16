input: clean
	gmsh -3 src/cake_slice.geo
	../../scripts/gmsh2triangle src/cake_slice.msh

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.face *.ele *checkpoint* *.log-0 *.err-0 ./src/*.msh
	rm -rf matrixdump*
