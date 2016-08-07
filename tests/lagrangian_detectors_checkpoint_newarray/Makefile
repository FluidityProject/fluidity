input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh

clean: 
	rm -rf *.log* *.err* *.stat *.detectors *flredecomp* lagrangian_detectors_checkpoint_?.vtu

logs: 
	../../bin/fluidity -v3 -l lagrangian_detectors_rk.flml
