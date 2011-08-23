input: clean
	gmsh -2 src/squareA.geo
	../../scripts/gmsh2triangle --2d -i src/squareA.msh
	gmsh -2 src/squareB.geo
	../../scripts/gmsh2triangle --2d -i src/squareB.msh

clean:
	rm -rf src/*.msh *.stat *.vtu *.log-0 *.err-0 *checkpoint *.convergence *.log *.node *.edge *.ele
