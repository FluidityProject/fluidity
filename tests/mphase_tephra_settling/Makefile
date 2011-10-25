preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_tephra_settling.msh src/mphase_tephra_settling.geo
	../../bin/gmsh2triangle --2d mphase_tephra_settling.msh

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_tephra_settling.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

