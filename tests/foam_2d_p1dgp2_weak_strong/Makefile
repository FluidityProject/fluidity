input: clean
	gmsh -2 src/drainage_a.geo
	../../scripts/gmsh2triangle --2d src/drainage_a.msh

clean:
	rm -rf *.stat *.vtu *.node *.edge *.ele *.log-0 *.err-0 ./src/*.msh
	rm -rf matrixdump*
