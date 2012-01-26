input: clean
	../../tools/gmsh2triangle.py --2d src/square.msh
	../../bin/fldecomp -n 2  square

clean:
	rm -fr *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu	\
	fluidity.* *.s *.d.1 *.stat heat_*
