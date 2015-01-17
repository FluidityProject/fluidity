preprocess:
	@echo **********Creating 2D mesh
	../../bin/fldecomp -n 2 -m gmsh src/mphase_volcanic_eruption_vw8_2d

run:
	@echo **********Running simulation
	mpiexec -n 2 ../../bin/fluidity mphase_volcanic_eruption_vw8_2d.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.halo
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_volcanic_eruption_vw8_2d_*
	rm -f nohup.out
