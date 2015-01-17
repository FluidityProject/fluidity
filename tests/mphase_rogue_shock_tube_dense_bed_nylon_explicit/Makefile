preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_rogue_shock_tube_dense_bed_nylon_explicit.msh src/mphase_rogue_shock_tube_dense_bed_nylon_explicit.geo
	../../bin/gmsh2triangle --2d mphase_rogue_shock_tube_dense_bed_nylon_explicit.msh

run:
	@echo **********Running simulation
	../../bin/fluidity mphase_rogue_shock_tube_dense_bed_nylon_explicit.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_rogue_shock_tube_dense_bed_nylon_explicit_* nohup.out
