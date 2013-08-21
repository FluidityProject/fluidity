input: clean
	gmsh -3 src/2d_square.geo
	../../bin/gmsh2triangle src/2d_square.msh

clean:
	rm -rf *.stat *.vtu *.pvtu *.msh *.node *.ele *.halo *.face *edge *checkpoint* *convergence* matrixdump* *.log-? *.err-? src/2d_square.msh adapted_state_* mmat-interpolation_[0-3] bounding_box* *metric* *flredecomp*
