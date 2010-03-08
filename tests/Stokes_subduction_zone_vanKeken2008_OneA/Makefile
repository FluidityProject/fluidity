input: clean
	gmsh -2 -o Subduction_Mesh_finer.msh src/Subduction_Mesh_finer.geo 
	../../scripts/gmsh2triangle --2d Subduction_Mesh_finer.msh

clean:
	rm -f  *.ele *.edge *.node *.vtu *.stat *.msh *.detectors

