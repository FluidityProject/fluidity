preprocess:
	@echo **********Creating 3D mesh
	gmsh -3 -o flow_past_a_cube.msh src/flow_past_a_cube.geo
	../../bin/gmsh2triangle flow_past_a_cube.msh
	../../bin/fldecomp -n 2 flow_past_a_cube

run:
	@echo **********Running simulation
	nohup mpiexec -n 2 ../../bin/fluidity flow_past_a_cube.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.face
	rm -f matrixdump* *.log* *.err*
	rm -rf flow_past_a_cube_*
	rm -f nohup.out
