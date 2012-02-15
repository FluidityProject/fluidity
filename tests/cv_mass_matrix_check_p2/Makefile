SHELL = sh

# options for 1d mesh interval script
MESH_SIZE = 1.2
LEFT_X    = 0.0
RIGHT_X   = 1.2

input: clean
	interval --dx=$(MESH_SIZE) $(LEFT_X) $(RIGHT_X) line
	gmsh -2 tri.geo -o tri.msh
	gmsh -3 tet.geo -o tet.msh

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
