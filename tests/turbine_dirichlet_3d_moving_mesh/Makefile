input: clean
	gmsh -3 -o mesh/mesh3_3d_thin.msh mesh/mesh3_3d_thin.geo 
	../../scripts/gmsh2triangle mesh/mesh3_3d_thin.msh

clean:
	rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh

