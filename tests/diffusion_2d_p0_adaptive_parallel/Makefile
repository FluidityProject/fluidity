input: clean
	../../bin/gmsh2triangle --2d src/square.msh

clean:
	rm -f  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.* *.s *.d.1 *.stat *detectors*
	rm -rf heat_? \
	matrixdump matrixdump.info
	rm -rf *flredecomp*
