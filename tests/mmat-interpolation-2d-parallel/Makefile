input: clean
	gmsh -2 src/2d_square.geo
	../../bin/gmsh2triangle --2d src/2d_square.msh
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l mmat-interpolation mmat-interpolation_flredecomp

clean:
	rm -rf *.stat *.vtu *.pvtu *.msh *.node *.ele *.halo *.edge *checkpoint* *convergence* matrixdump* *.log-? *.err-? src/2d_square.msh adapted_state_* mmat-interpolation_[0-3] *flredecomp*
