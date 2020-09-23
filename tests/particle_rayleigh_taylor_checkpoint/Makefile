input: default
default: clean
	gmsh -2 -bin square.geo -o square.msh
clean: 
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.particles.* *flredecomp* particle_attributes_? *~ particle_attributes_?_* particle_attributes*checkpoint* particle_attributes_??

logs: 
	mpirun -n 4 ../../bin/fluidity -v3 -l particle-attributes.flml

grep:
	mpirun -n 4 ../../bin/fluidity -v3 particle-attributes.flml | grep -i cjc

debug:
	mpirun -n 4 xterm -e gdb -x gdb_run --args ../../bin/fluidity -v2 -l particle-attributes.flml

