input: clean
	gmsh -2 ./src/box_2d.geo
	../../bin/gmsh2triangle --2d ./src/box_2d.msh
	../../bin/fldecomp -n 2 box_2d

clean: 
	rm -f  ./src/*.msh *.ele *.edge *.node *.halo *.face *.poly *vtu *.s *.stat *.log-* *.err-* matrixdump*
	rm -rf test_metric_advection_parallel_?
