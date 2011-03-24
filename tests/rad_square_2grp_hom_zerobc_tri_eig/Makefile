SHELL = sh

SIM = rad_square_2grp_hom_zerobc_tri_eig

MESH = square_tri

MESH_DIM = 2

input: clean
	gmsh -$(MESH_DIM) $(MESH)_A.geo -o $(MESH)_A.msh
	gmsh -$(MESH_DIM) $(MESH)_B.geo -o $(MESH)_B.msh

clean:
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
