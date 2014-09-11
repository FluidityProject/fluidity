input: clean
	gmsh -2 mesh/np_rectangle.geo -o np_rectangle.msh; gmsh2triangle -s np_rectangle.msh;

clean:
	rm -f  fluidity.* *.ele *.face *.node *.vtu *.poly *.edge *.s *.d.1 *.stat \
	matrixdump matrixdump.info
	rm -f *.msh *.edge *.ele *.node *.halo
	rm -f *flredecomp.flml
