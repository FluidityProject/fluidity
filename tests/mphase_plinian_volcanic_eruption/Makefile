preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_plinian_volcanic_eruption.msh src/mphase_plinian_volcanic_eruption.geo
	../../bin/gmsh2triangle --2d mphase_plinian_volcanic_eruption.msh
	../../bin/fldecomp -n 8 mphase_plinian_volcanic_eruption

run:
	@echo **********Running simulation
	nohup mpiexec -n 8 ../../bin/fluidity mphase_plinian_volcanic_eruption.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_plinian_volcanic_eruption_*
	rm -f nohup.out
