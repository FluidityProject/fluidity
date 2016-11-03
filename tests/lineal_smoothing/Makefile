input:	clean
	gmsh -2 -clmax 0.05 src/square.geo
	gmsh -3 -clmax 0.02 src/unit_cube.geo

clean:
	rm -f *.vtu *.stat *.s *.log  src/*A.*matrix* *~ fluidity.* src/*.msh
