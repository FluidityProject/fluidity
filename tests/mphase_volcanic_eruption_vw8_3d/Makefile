preprocess:
	@echo **********Creating 3D mesh
	gmsh -3 -o mphase_volcanic_eruption_vw8_3d.msh src/mphase_volcanic_eruption_vw8_3d.geo
	../../bin/gmsh2triangle mphase_volcanic_eruption_vw8_3d.msh
	../../bin/fldecomp -n 2 mphase_volcanic_eruption_vw8_3d
run:
	@echo **********Running simulation
	mpiexec -n 2 ../../bin/fluidity -v2 -l mphase_volcanic_eruption_vw8_3d.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.face *.halo
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_volcanic_eruption_vw8_3d_*
	rm -rf *flredecomp*
