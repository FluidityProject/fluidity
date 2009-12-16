input: clean
	gmsh -2 src/2d_square.geo
	../../scripts/gmsh2triangle --2d src/2d_square.msh

clean:
	rm -rf *.d.* *.stat *.vtu *.convergence *.node *.face *.ele *checkpoint*
