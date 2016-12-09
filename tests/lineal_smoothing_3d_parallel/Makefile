input:	clean
	gmsh -3 -clmax 0.075 src/unit_cube.geo
clean:
	rm -f *.pvtu *.msh *.halo *.vtu *.stat *.log* matrixdump* *~ fluidity.* *.err* src/*.msh
	rm -rf lineal_smoothing_3d_?
