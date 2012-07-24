input: clean
	../../bin/gmsh2triangle -2 src/circle-serial.msh
	../../bin/gmsh2triangle -2 src/circle-parallel.msh
	mpiexec -n 3 ../../bin/flredecomp -i 1 -o 3 -v -l circle-2d-drag-parallel circle-2d-drag-parallel_flredecomp
	mpiexec -n 3 ../../bin/flredecomp -i 1 -o 3 -v -l circle-2d-drag-parallel-adaptive circle-2d-drag-parallel-adaptive_flredecomp

clean:
	rm -f *.halo *.ele *.edge *.face *.node *.poly *.vtu *.pvtu *.s *.d.* *.stat *.s *.log-0 *.err-0
	rm -rf circle-2d-drag-parallel_? circle-2d-drag-parallel-adaptive_?
	rm -rf *flredecomp*
