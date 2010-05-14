input: clean
	gmsh -2 ./src/box_2d.geo
	../../scripts/gmsh2triangle --2d ./src/box_2d.msh
	../../bin/fldecomp -n 4 box_2d

clean: 
	rm -f  ./src/*.msh *.ele *.edge *.node *.halo *.face *.poly *vtu *.s  *.stat *.log-* *.err-* matrixdump*
