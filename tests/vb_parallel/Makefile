input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh
	fldecomp -m gmsh -n 2 square
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log *.stat blob_[012]