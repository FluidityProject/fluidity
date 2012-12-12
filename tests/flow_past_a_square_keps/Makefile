preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o flow_past_a_square_keps.msh src/flow_past_a_square_keps.geo
	../../bin/gmsh2triangle --2d flow_past_a_square_keps.msh
	../../bin/fldecomp -n 2 flow_past_a_square_keps

run:
	@echo **********Running simulation
	nohup mpiexec -n 2 ../../bin/fluidity flow_past_a_square_keps.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.face
	rm -f matrixdump* *.log* *.err*
	rm -rf flow_past_a_square_keps_*
	rm -f nohup.out
