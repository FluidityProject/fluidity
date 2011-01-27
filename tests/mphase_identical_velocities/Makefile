preprocess:
	@echo **********Creating mesh for refinement N=1/16
	gmsh -2 -o mphase_identical_velocities.msh src/mphase_identical_velocities.geo
	../../scripts/gmsh2triangle --2d mphase_identical_velocities.msh

run:
	@echo **********Running simulation for refinement N=1/16
	../../bin/fluidity -v2 -l mphase_identical_velocities.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *log* *err*

