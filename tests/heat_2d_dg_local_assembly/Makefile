input: clean
	../../bin/gmsh2triangle --2d src/square.msh
	mpiexec -n 8 ../../bin/flredecomp -i 1 -o 8 -v -l heat heat_flredecomp

clean:
	rm -f  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.* *.s *.d.1 *.stat \
	matrixdump matrixdump.info
	rm -rf heat_[01] *flredecomp*
