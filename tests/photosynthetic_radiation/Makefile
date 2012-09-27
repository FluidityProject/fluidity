input: clean
	../../bin/gmsh2triangle src/column.msh
	../../bin/fldecomp -n 2 column

clean:
	rm -f *.halo *.face *.ele *.edge *.node *.poly *.pvtu *.vtu *.s *.d.1 *.stat \
	matrixdump matrixdump.info
	rm -rf photosynthetic_radiation_?
