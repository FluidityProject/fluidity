SHELL = sh

SIM = darcy_p0p1cv_velBCinlet

# options for 1d mesh interval script
MESH_SIZE = 100.0
LEFT_X    = 0.0
RIGHT_X   = 300.0

input: clean
	interval --dx=$(MESH_SIZE) $(LEFT_X) $(RIGHT_X) $(SIM)_1d
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
	rm -f fluidity.err*
	rm -f fluidity.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
