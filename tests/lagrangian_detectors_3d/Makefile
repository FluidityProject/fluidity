input: default
default: clean
	gmsh -3 -bin src/box.geo -o box.msh
	mpiexec -n 4 ../../bin/flredecomp -i 1 -o 4 -v -l lagrangian_detectors_3d lagrangian_detectors_3d_flredecomp
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.detectors *.dat lagrangian_detectors_[01] *flredecomp*
