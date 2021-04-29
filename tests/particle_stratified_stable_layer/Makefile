input: default
default: clean
	gmsh -2 -bin mesh/ar2.geo
clean:
	rm -rf mesh/*.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.particles.* *flredecomp* 2d_entrainment_? *~ 2d_entrainment_?_* *checkpoint* 2d_entrainment_?? __pycache__
