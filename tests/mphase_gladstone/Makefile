preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_gladstone.msh src/mphase_gladstone.geo
	../../bin/gmsh2triangle --2d mphase_gladstone.msh
	../../bin/fldecomp -n 4 mphase_gladstone

run:
	@echo **********Running simulation
	nohup mpiexec -n 4 ../../bin/fluidity mphase_gladstone.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_gladstone_* nohup.out
	
