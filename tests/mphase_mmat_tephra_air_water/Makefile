preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_mmat_tephra_air_water.msh src/mphase_mmat_tephra_air_water.geo
	../../bin/gmsh2triangle --2d mphase_mmat_tephra_air_water.msh

run:
	@echo **********Running simulation
	../../bin/fluidity mphase_mmat_tephra_air_water.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
