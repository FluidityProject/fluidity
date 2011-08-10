input: clean
	gmsh -2 src/squareA.geo
	gmsh -2 src/squareB.geo
	gmsh -2 src/squareC.geo
	gmsh -2 src/squareD.geo

clean:
	rm -rf src/*.msh *.stat *.vtu *.log-0 *.err-0 *checkpoint *.convergence
