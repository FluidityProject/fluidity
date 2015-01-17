preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_fluidised_bed.msh src/mphase_fluidised_bed.geo
	../../bin/gmsh2triangle --2d mphase_fluidised_bed.msh
	../../bin/fldecomp -n 4 mphase_fluidised_bed

run:
	@echo **********Running simulation
	nohup mpiexec -n 4 ../../bin/fluidity mphase_fluidised_bed.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_fluidised_bed_* nohup.out
