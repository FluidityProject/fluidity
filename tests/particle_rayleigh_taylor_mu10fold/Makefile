input: default
default: clean
	gmsh -2 -bin mesh/square.geo
clean: 
	rm -rf mesh/*.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.particles.* *flredecomp* particle_attributes_? *~ particle_attributes_?_* *checkpoint* particle_attributes_?? particle_attributes_???
