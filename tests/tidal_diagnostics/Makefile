input: clean
	gmsh -2 src/sphere.geo -o sphere.msh; gmsh2triangle -s sphere.msh

clean:
	rm -f  fluidity.* *.ele *.face *.node *.vtu *.poly *.edge *.s *.d.1 *.stat \
	matrixdump matrixdump.info
	rm -f *.msh
