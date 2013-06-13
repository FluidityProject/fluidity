input: clean
	../../bin/gmsh2triangle --2d src/square.msh
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l heat heat_flredecomp

clean:
	rm -fr  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.* *.s *.d.1 *.stat heat_[01] \
	matrixdump matrixdump.info
	rm -rf *flredecomp*
