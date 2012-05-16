SHELL = sh

SIM = darcy_impes_p1_2phase_coreyrelperm_velBCinlet_strongpressoutlet_p1satdiag

# options for 1d mesh interval script
MESH_SIZE_A = 0.1
MESH_SIZE_B = 0.05
MESH_SIZE_C = 0.025
MESH_SIZE_D = 0.0125
LEFT_X    = 0.0
RIGHT_X   = 1.0

input: clean
	interval --dx=$(MESH_SIZE_A) $(LEFT_X) $(RIGHT_X) line_mesh_1d_A
	interval --dx=$(MESH_SIZE_B) $(LEFT_X) $(RIGHT_X) line_mesh_1d_B
	interval --dx=$(MESH_SIZE_C) $(LEFT_X) $(RIGHT_X) line_mesh_1d_C
	interval --dx=$(MESH_SIZE_D) $(LEFT_X) $(RIGHT_X) line_mesh_1d_D
	gmsh -2 square.geo -o square.msh
	gmsh -3 cube.geo -o cube.msh

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
