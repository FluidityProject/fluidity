input: default
default: clean
	gmsh -2 -bin mesh/ar2.geo
clean:
	rm -rf mesh/*.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.particles.* *flredecomp* 2d_entrainment_? *~ 2d_entrainment_?_* *checkpoint* 2d_entrainment_?? __pycache__

logs:
	mpirun -n 4 ../../bin/fluidity -v3 -l particle-attributes.flml

grep:
	mpirun -n 4 ../../bin/fluidity -v3 particle-attributes.flml | grep -i cjc

debug:
	mpirun -n 4 xterm -e gdb -x gdb_run --args ../../bin/fluidity -v2 -l particle-attributes.flml
