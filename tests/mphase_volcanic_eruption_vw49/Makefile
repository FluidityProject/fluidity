preprocess:
	@echo **********Creating 2D mesh
	../../bin/fldecomp -n 4 -m gmsh src/mphase_volcanic_eruption_vw49

run:
	@echo **********Running simulation
	nohup mpiexec -n 4 ../../bin/fluidity mphase_volcanic_eruption_vw49.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.halo
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_volcanic_eruption_vw49_*
	rm -f nohup.out
