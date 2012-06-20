preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_rogue_shock_tube_dense_bed_glass.msh src/mphase_rogue_shock_tube_dense_bed_glass.geo
	../../bin/gmsh2triangle --2d mphase_rogue_shock_tube_dense_bed_glass.msh
	../../bin/fldecomp -n 4 mphase_rogue_shock_tube_dense_bed_glass

run:
	@echo **********Running simulation
	nohup mpiexec -n 4 ../../bin/fluidity mphase_rogue_shock_tube_dense_bed_glass.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_rogue_shock_tube_dense_bed_glass_* nohup.out
