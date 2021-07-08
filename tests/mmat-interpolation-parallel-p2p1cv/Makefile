input: clean
	gmsh -2 src/2d_square.geo
	cp src/2d_square.msh .

clean:
	rm -rf mmat-interpolation-*_* *-flredecomp.flml *.err-* *.log-* *.stat *.halo *.[p]vtu *.node *.ele *.edge *checkpoint* *convergence* src/2d_square.msh *.msh \
	matrixdump matrixdump.info
