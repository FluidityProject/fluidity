input: clean
	../../scripts/gmsh2triangle --2d src/square.msh
	../../bin/fldecomp -n 2 -i triangle square

clean:
	rm -f  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.* *.s *.d.1 *.stat
