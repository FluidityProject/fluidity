preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_pyroclastic_density_current.msh src/mphase_pyroclastic_density_current.geo
	../../bin/gmsh2triangle --2d mphase_pyroclastic_density_current.msh
	../../bin/fldecomp -n 4 mphase_pyroclastic_density_current

run:
	@echo **********Running simulation
	nohup mpiexec -n 4 ../../bin/fluidity mphase_pyroclastic_density_current.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_pyroclastic_density_current_*
	rm -f nohup.out
