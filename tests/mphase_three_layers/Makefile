preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_three_layers.msh src/mphase_three_layers.geo

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_three_layers.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

