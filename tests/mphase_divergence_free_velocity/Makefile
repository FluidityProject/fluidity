preprocess:
	@echo **********Creating 3D mesh
	gmsh -3 -o mphase_divergence_free_velocity.msh src/mphase_divergence_free_velocity.geo

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_divergence_free_velocity.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *log* *err*
	rm -f *.detectors

