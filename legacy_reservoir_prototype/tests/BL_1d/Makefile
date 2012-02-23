SHELL = sh

SIM = BL_1d

# options for 1d mesh interval script for different mesh resolutions
LEFT_X    = 0.0
RIGHT_X   = 4.0

MESH_SIZE_A = 0.4
MESH_SIZE_B = 0.2

input: clean
	@echo "Running INTERVAL script for $(SIM)"
	../../../scripts/interval --dx=$(MESH_SIZE_A) $(LEFT_X) $(RIGHT_X) $(SIM)_A
	../../../scripts/interval --dx=$(MESH_SIZE_B) $(LEFT_X) $(RIGHT_X) $(SIM)_B

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
	rm -f gmon.out
	rm -f *.d.*
	rm -f mirror_int_data.dat input.dat mirror_new.dat
	rm -f Cubic_Interp.pyc
	rm -rf $(SIM)_*output
	rm -f CheckSparsityMatrix.dat
	rm -rf velocity_projected_*
	rm -rf ocvn
