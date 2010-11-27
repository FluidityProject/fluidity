input: clean
	gmsh -2 src/column_2d.geo
	../../scripts/gmsh2triangle -2 src/column_2d.msh

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.edge *.ele *checkpoint*
