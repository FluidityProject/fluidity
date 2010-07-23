input: clean
	gmsh -2 2dtank.geo
	gmsh2triangle -2 2dtank.msh

clean:
	rm -rf *.stat *.msh *.node *.edge *.ele adapted_mesh_*.vtu water_collapse_*.vtu *checkpoint* *convergence* *.log-0 *.err-0 matrixdump* *.detectors
