SHELL = sh

SIM = darcy_impes_p1_2phase_coreyrelperm_velBCinlet_strongpressoutlet_p1satdiag

input: clean
	cp mesh_data/line_mesh_1d_A* .
	cp mesh_data/line_mesh_1d_B* .
	cp mesh_data/line_mesh_1d_C* .
	cp mesh_data/line_mesh_1d_D* .
	cp mesh_data/square_A.msh .
	cp mesh_data/square_B.msh .
	cp mesh_data/square_C.msh .
	cp mesh_data/cube_A.msh .
	cp mesh_data/cube_B.msh .

clean:
	rm -f *.ele
	rm -f *.node
	rm -f *.bound
	rm -f *.edge
	rm -f *.face
	rm -f *.vtu 
	rm -f *.pvtu 
	rm -f *.s 
	rm -f *.stat 
	rm -f *.log-0 
	rm -f *.err-0 
	rm -f *.msh
	rm -f *.halo
	rm -f $(SIM)*.err*
	rm -f $(SIM)*.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
	rm -f *.convergence 
	rm -f *.pyc
	rm -f *.txt
