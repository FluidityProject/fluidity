input: clean
	gmsh -2 -optimize Subduction_Mesh.msh src/Subduction_Mesh.geo
	../../scripts/gmsh2triangle --2d -i Subduction_Mesh.msh

clean:
	rm -f *.ele *.edge *.node *.vtu *.stat *.msh *.detectors *.dat fluidity.* *_state matrixdump*
