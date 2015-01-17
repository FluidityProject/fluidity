`preprocess:
	@echo **********Creating 2D mesh
	../../bin/gmsh2triangle --2d src/mphase_plinian_eruption_neri.msh
	../../bin/fldecomp -n 2 mphase_plinian_eruption_neri

run:
	@echo **********Running simulation
	nohup mpiexec -n 2 ../../bin/fluidity mphase_plinian_eruption_neri.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_plinian_eruption_neri_*
	rm -f nohup.out
