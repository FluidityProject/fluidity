SHELL = sh

SIM = add_src_directly_1d_2grp_hom_zerobc_eig_p1cv

# options for 1d mesh interval script for different mesh resolutions
LEFT_X    = 0.0
RIGHT_X   = 44.5

MESH_SIZE_A = 2.0
MESH_SIZE_B = 1.0
MESH_SIZE_C = 0.5
MESH_SIZE_D = 0.25

input: clean
	@echo "Running INTERVAL script for $(SIM) for each mesh resolution"
	interval --dx=$(MESH_SIZE_A) $(LEFT_X) $(RIGHT_X) $(SIM)_A
	interval --dx=$(MESH_SIZE_B) $(LEFT_X) $(RIGHT_X) $(SIM)_B
	interval --dx=$(MESH_SIZE_C) $(LEFT_X) $(RIGHT_X) $(SIM)_C
	interval --dx=$(MESH_SIZE_D) $(LEFT_X) $(RIGHT_X) $(SIM)_D

clean:
	rm -f $(SIM)*.ele
	rm -f $(SIM)*.node
	rm -f $(SIM)*.bound
	rm -f $(SIM)*.edge
	rm -f $(SIM)*.face
	rm -f $(SIM)*.vtu 
	rm -f $(SIM)*.pvtu 
	rm -f $(SIM)*.s 
	rm -f $(SIM)*.d.* 
	rm -f $(SIM)*.stat 
	rm -f $(SIM)*.log-0 
	rm -f $(SIM)*.err-0 
	rm -f $(SIM)*.msh
	rm -f $(SIM)*.halo
	rm -f fluidity.err*
	rm -f fluidity.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
	rm -f $(SIM)*.eig_iteration_convergence
	rm -f $(SIM)*checkpoint*
