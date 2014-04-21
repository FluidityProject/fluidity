input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh

clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.detectors *.dat *.groups *flredecomp* *checkpoint.flml

logs: 
	../../bin/fluidity -v3 -l lagrangian_detectors_rk.flml
