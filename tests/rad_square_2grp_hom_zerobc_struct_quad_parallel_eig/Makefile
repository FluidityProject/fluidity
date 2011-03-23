SHELL = sh

SIM = rad_square_2grp_hom_zerobc_struct_quad_parallel_eig

MESH = square_struct_quad

MESH_DIM = 2

NPROC = 2

input: clean
	gmsh -bin -$(MESH_DIM) $(MESH)_A.geo -o $(MESH)_A.msh
	gmsh -bin -$(MESH_DIM) $(MESH)_B.geo -o $(MESH)_B.msh
	gmsh -bin -$(MESH_DIM) $(MESH)_C.geo -o $(MESH)_C.msh
	gmsh -bin -$(MESH_DIM) $(MESH)_D.geo -o $(MESH)_D.msh
	../../bin/fldecomp -n $(NPROC) -m gmsh $(MESH)_A
	../../bin/fldecomp -n $(NPROC) -m gmsh $(MESH)_B
	../../bin/fldecomp -n $(NPROC) -m gmsh $(MESH)_C
	../../bin/fldecomp -n $(NPROC) -m gmsh $(MESH)_D

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
	rm -rf $(SIM)*_0
	rm -rf $(SIM)*_1
