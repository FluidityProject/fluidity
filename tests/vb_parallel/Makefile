input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log *.stat blob_[012 \
	matrixdump matrixdump.info
	rm -rf *flredecomp*
