input: clean
	gmsh -3 -o column.msh src/column.geo
	../../bin/gmsh2triangle column.msh
	../../bin/fldecomp -n 2 column

clean: 
	rm -f  *.ele *.edge *.node *.face *.poly *.pvtu *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump* *log *err
	rm -rf column-parallel_* *.halo
