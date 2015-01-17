preprocess:
	@echo **********Creating 2D mesh
	../../bin/gmsh2triangle --2d src/mphase_volcanic_eruption_one_element_wide.msh
	../../bin/fldecomp -n 2 mphase_volcanic_eruption_one_element_wide

run:
	@echo **********Running simulation
	nohup mpiexec -n 2 ../../bin/fluidity mphase_volcanic_eruption_one_element_wide.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.halo
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_volcanic_eruption_one_element_wide_*
	rm -f nohup.out
