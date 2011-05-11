preprocess:
	@echo **********Creating 2D mesh
	../../scripts/gmsh2triangle --2d src/mphase_identical_velocities_dg.msh
   
run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_identical_velocities_dg.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

