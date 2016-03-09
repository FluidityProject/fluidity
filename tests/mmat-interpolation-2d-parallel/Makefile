input: clean
	gmsh -2 src/2d_square.geo
	cp src/2d_square.msh .

clean:
	rm -rf *.stat *.msh *.vtu *.pvtu *.node *.ele *.halo *.edge *checkpoint* *convergence* matrixdump* *.log-? *.err-? src/2d_square.msh adapted_state_* mmat-interpolation_[0-3] *flredecomp*
