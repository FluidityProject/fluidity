input: clean
	gmsh -2 src/2d_square.geo
	cp src/2d_square.msh .

clean:
	rm -rf *.stat *.vtu *.msh *checkpoint* *convergence* *.log

merged_mesh:
	unifiedmesh explicit-hyperc-shear-conservative_CoordinateMesh_0_checkpoint explicit-hyperc-shear-conservative_checkpoint_CoordinateMesh_16_checkpoint merged_mesh
