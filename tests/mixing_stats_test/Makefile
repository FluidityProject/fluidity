input: clean
	sed  "s/mixing_stats_test/mixing_stats_test_serial/" mixing_stats_test.flml > mixing_stats_test_serial.flml
	sed  "s/mixing_stats_test/mixing_stats_test_parallel/" mixing_stats_test.flml > mixing_stats_test_parallel.flml
	gmsh  -bin -3 src/pseudo_2d_square.geo -o pseudo_2d_square.msh
	mpiexec -n 3 ../../bin/flredecomp -i 1 -o 3 -v -l mixing_stats_test_parallel mixing_stats_test_parallel_flredecomp

clean:
	rm -f *.stat *.vtu *.pvtu *.convergence *checkpoint* *.err* *.log*  
	rm -f mixing_stats_test_serial.flml mixing_stats_test_parallel.flml
	rm -f *.halo *.node *.ele *.face *.msh
	rm -rf mixing_stats_test_parallel_?
	rm -rf *flredecomp*