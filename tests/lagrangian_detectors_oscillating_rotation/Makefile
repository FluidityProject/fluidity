input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.detectors *.dat *flredecomp* lagrangian_detectors_?

logs: 
	mpirun -n 4 ../../bin/fluidity -v3 -l lagrangian_detectors_rk.flml

grep:
	mpirun -n 4 ../../bin/fluidity -v3 lagrangian_detectors_rk.flml | grep -i cjc

debug:
	mpirun -n 4 xterm -e gdb -x gdb_run --args ../../bin/fluidity -v2 -l lagrangian_detectors_rk.flml

