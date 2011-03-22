SHELL = sh

SIM = rad_cube_2grp_hom_zerobc_anisodiff_struct_hex_eig

MESH = box_struct_hex

input: clean
	gmsh -3 $(MESH)_A.geo -o $(MESH)_A.msh
	gmsh -3 $(MESH)_B.geo -o $(MESH)_B.msh

clean:
	rm -f $(MESH)*.ele
	rm -f $(MESH)*.node
	rm -f $(MESH)*.bound
	rm -f $(MESH)*.edge
	rm -f $(MESH)*.face
	rm -f $(SIM)*.vtu 
	rm -f $(SIM)*.pvtu 
	rm -f $(SIM)*.stat 
	rm -f $(SIM)*.log-0 
	rm -f $(SIM)*.err-0 
	rm -f $(MESH)*.msh
	rm -f $(MESH)*.halo
	rm -f fluidity.err*
	rm -f fluidity.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
