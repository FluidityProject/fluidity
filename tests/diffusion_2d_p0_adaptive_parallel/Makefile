input: clean
	../../bin/gmsh2triangle --2d src/square.msh
	mpiexec -n 3 ../../bin/flredecomp -i 1 -o 3 -v -l heat heat_flredecomp

clean:
	rm -f  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.* *.s *.d.1 *.stat *detectors*
	rm -rf heat_?
	rm -rf *flredecomp*
