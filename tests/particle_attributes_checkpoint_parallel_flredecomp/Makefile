input: default
default: clean
	gmsh -format msh2 -2 -bin square.geo -o square.msh
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.particles.* particle_attributes_? *~ *.h5part P4* particle_attributes_checkpoint_?
