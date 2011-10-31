SHELL = sh

SIM = darcy_p1dgp2_pressBCinlet_1d

# options for 1d mesh interval script
MESH_SIZE = 100.0
LEFT_X    = 0.0
RIGHT_X   = 300.0

input: clean
	@echo "Running INTERVAL script for $(SIM)"
	../../bin/interval --dx=$(MESH_SIZE) $(LEFT_X) $(RIGHT_X) $(SIM)

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
	rm -rf $(SIM)_*
	rm -f fluidity.err*
	rm -f fluidity.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
