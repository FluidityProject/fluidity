preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o bc_from_field.msh src/bc_from_field.geo
	../../bin/gmsh2triangle --2d bc_from_field.msh

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l bc_from_field.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

