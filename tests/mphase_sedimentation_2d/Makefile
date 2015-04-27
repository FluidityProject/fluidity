preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_sedimentation_2d.msh src/mphase_sedimentation_2d.geo

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_sedimentation_2d.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

