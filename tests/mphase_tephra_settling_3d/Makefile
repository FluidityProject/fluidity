preprocess:
	@echo **********Creating 3D mesh
	gmsh -3 -o mphase_tephra_settling_3d.msh src/mphase_tephra_settling_3d.geo
	../../scripts/gmsh2triangle mphase_tephra_settling_3d.msh

run:
	@echo **********Running simulation
	../../bin/fldecomp -n 4 mphase_tephra_settling_3d
	mpiexec -n 4 ../../bin/fluidity -v2 -l mphase_tephra_settling_3d.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.face
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_tephra_settling_3d_* nohup.out


