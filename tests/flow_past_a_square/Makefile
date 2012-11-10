preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o flow_past_a_square.msh src/flow_past_a_square.geo
	../../bin/gmsh2triangle --2d flow_past_a_square.msh
	../../bin/fldecomp -n 16 flow_past_a_square

run:
	@echo **********Running simulation
	nohup mpiexec -n 2 ../../bin/fluidity flow_past_a_square.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.face
	rm -f matrixdump* *.log* *.err*
	rm -rf flow_past_a_square_*
	rm -f nohup.out
