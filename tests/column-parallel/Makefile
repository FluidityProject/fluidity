input: clean
	gmsh -3 -o column.msh src/column.geo

clean: 
	rm -f  *.ele *.edge *.msh *.node *.face *.poly *.pvtu *.vtu *.s  *.stat *.log-* *.err-* matrixdump* *log *err
	rm -rf column-parallel_* *.halo
