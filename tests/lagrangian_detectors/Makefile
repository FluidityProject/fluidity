input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh
	fldecomp -m gmsh -n 4 square
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.detectors *.dat