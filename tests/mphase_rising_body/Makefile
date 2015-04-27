preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_rising_body.msh src/mphase_rising_body.geo

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_rising_body.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

