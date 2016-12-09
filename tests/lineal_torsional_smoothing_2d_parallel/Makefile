input:	clean
	gmsh -2 -clmax 0.05 src/square.geo
clean:
	rm -f *.pvtu *.vtu *.stat *.log* *.msh  *.halo matrixdump* *~ fluidity.* *.err* src/*.msh
	rm -rf lineal_torsional_smoothing_2d_?
