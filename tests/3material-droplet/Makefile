input: clean
	gmsh -2 src/square.geo
	../../scripts/gmsh2triangle --2d src/square.msh

clean:
	rm -rf src/*.msh *.stat *.vtu *.node *.ele *.edge *.log-0 *.err-0 *checkpoint *.convergence
