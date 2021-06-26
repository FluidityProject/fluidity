preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_tephra_settling_salinity_adaptive.msh src/mphase_tephra_settling_salinity_adaptive.geo

run:
	@echo **********Running simulation
	mpiexec -n 8 flredecomp -i 1 -o 8 -v -l mphase_tephra_settling_salinity_adaptive mphase_tephra_settling_salinity_adaptive_flredecomp
	mpiexec -n 8 fluidity mphase_tephra_settling_salinity_adaptive_flredecomp.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh *.halo
	rm -f *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_tephra_settling_salinity_adaptive_*
	rm -rf *flredecomp*
