preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_sedimentation.msh src/mphase_sedimentation.geo
	../../scripts/gmsh2triangle --2d mphase_sedimentation.msh

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_sedimentation.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

