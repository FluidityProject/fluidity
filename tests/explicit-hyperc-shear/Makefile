input: clean
	gmsh -2 src/2d_square.geo
	cp src/2d_square.msh .

clean:
	rm -rf *.stat explicit-hyperc-shear*.vtu *.msh explicit-hyperc-shear*checkpoint* explicit-hyperc-shear*convergence*
	rm -rf initial_setup_?.vtu initial_setup_PressureMesh_?_checkpoint.vtu initial_setup_VelocityMesh_0_checkpoint.vtu initial_setup_?_checkpoint.flml
	rm -rf reverse_setup_?.vtu reverse_setup_PressureMesh_?_checkpoint.vtu reverse_setup_VelocityMesh_0_checkpoint.vtu reverse_setup_?_checkpoint.flml
	rm -rf *.log-0 *.err-0
