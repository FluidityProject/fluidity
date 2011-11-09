input: default
default: clean
	gmsh -3 -bin src/box.geo -o box.msh
	../../bin/fldecomp -m gmsh -n 4 box
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.detectors *.dat lagrangian_detectors_[01]
