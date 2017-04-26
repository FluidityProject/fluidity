input: clean
	cd mesh; gmsh -2 np_rectangle.geo

clean:
	rm -f  fluidity.* *.ele *.face *.node *.vtu *.poly *.edge *.s *.d.1 *.stat \
	matrixdump matrixdump.info
	rm -f mesh/*.msh mesh/*.halo
	rm -f *.msh *.halo
	rm -f *.pvtu *.vtu
	rm -f *flredecomp.flml
