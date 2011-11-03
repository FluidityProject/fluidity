preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_divergence_free_velocity.msh src/mphase_divergence_free_velocity.geo
	../../bin/gmsh2triangle --2d mphase_divergence_free_velocity.msh

run:
	@echo **********Running simulation
	../../bin/fluidity -v3 -l mphase_divergence_free_velocity.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *checkpoint* *.log* *.err*
