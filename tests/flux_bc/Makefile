preprocess:
	@echo **********Creating meshes
	gmsh -2 -o flux_bc_2d.msh src/flux_bc_2d.geo
	gmsh -3 -o flux_bc_3d.msh src/flux_bc_3d.geo

run:
	@echo **********Running 2D simulation
	../../bin/fluidity flux_bc_2d.flml 
	@echo **********Running 3D simulation
	../../bin/fluidity flux_bc_3d.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.face
	rm -f matrixdump* *checkpoint* *.log* *.err*

