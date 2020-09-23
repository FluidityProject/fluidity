input: default
default: clean
	gmsh -2 -bin mesh.geo -o mesh.msh
clean:
	rm -rf *.msh *.halo *.vtu *.pvtu *.log* *.err* *.stat blob_[012] *.particles.* *flredecomp* circle_? *~ circle_?_* *checkpoint* circle_??

logs:
	mpirun -n 4 ../../bin/fluidity -v3 -l particle-attributes.flml

grep:
	mpirun -n 4 ../../bin/fluidity -v3 particle-attributes.flml | grep -i cjc

debug:
	mpirun -n 4 xterm -e gdb -x gdb_run --args ../../bin/fluidity -v2 -l particle-attributes.flml
