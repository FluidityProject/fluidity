preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_stokes_law.msh src/mphase_stokes_law.geo
	../../bin/gmsh2triangle --2d mphase_stokes_law.msh

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_stokes_law.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
