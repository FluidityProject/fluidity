preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_tephra_settling.msh src/mphase_tephra_settling.geo

run:
	@echo **********Running simulation
	../../bin/fluidity mphase_tephra_settling.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

