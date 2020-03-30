input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.particles.* *flredecomp* particle_attributes_? *~ particle_attributes_?_checkpoint* particle_attributes_checkpoint_?_checkpoint* *.h5part P4* particle_attributes_checkpoint_?
