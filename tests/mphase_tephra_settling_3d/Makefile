preprocess:
	@echo **********Creating 3D mesh
	gmsh -3 -o mphase_tephra_settling_3d.msh src/mphase_tephra_settling_3d.geo

run:
	@echo **********Running simulation
	mpiexec -n 8 flredecomp -i 1 -o 8 -v -l mphase_tephra_settling_3d mphase_tephra_settling_3d_flredecomp
	mpiexec -n 8 fluidity -v2 -l mphase_tephra_settling_3d_flredecomp.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh *.halo
	rm -f *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_tephra_settling_3d_*
	rm -rf *flredecomp*
