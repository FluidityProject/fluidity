preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_pyroclastic_density_current_obstruction.msh src/mphase_pyroclastic_density_current_obstruction.geo
	../../bin/gmsh2triangle --2d mphase_pyroclastic_density_current_obstruction.msh
	../../bin/fldecomp -n 8 mphase_pyroclastic_density_current_obstruction

run:
	@echo **********Running simulation
	nohup mpiexec -n 8 ../../bin/fluidity mphase_pyroclastic_density_current_obstruction.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_pyroclastic_density_current_obstruction_*
	rm -f nohup.out
