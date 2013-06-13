input: clean
	gmsh -3 -o column.msh src/column.geo
	../../bin/gmsh2triangle column.msh
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l column-parallel column-parallel_flredecomp

clean: 
	rm -f  *.ele *.edge *.msh *.node *.face *.poly *.pvtu *.vtu *.s  *.stat *.log-* *.err-* matrixdump* *log *err
	rm -rf column-parallel_* *.halo
