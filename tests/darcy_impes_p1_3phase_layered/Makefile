SHELL = sh

SIM = darcy_impes_p1_3phase_layered

# options for 1d mesh interval script
MESH_SIZE = 0.05
LEFT_X    = 0.0
RIGHT_X   = 1.0

input: clean
	interval --dx=$(MESH_SIZE) $(LEFT_X) $(RIGHT_X) $(SIM)_1d
	cp mesh_data/square*.msh .
	cp mesh_data/cube*.msh .

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
