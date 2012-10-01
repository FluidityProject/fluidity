input: clean
	../../bin/gmsh2triangle --2d src/square.msh
	../../bin/fldecomp -n 16  square

clean:
	rm -rf  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.* *.s *.d.1 *.stat \
	matrixdump matrixdump.info heat_*
