preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o mphase_rogue_shock_tube.msh src/mphase_rogue_shock_tube.geo
	../../bin/gmsh2triangle --2d mphase_rogue_shock_tube.msh
	#../../bin/gmsh2triangle --2d first_timestep_adapted_mesh.msh
	../../bin/fldecomp -n 4 mphase_rogue_shock_tube

run:
	@echo **********Running simulation
	nohup mpiexec -n 4 ../../bin/fluidity mphase_rogue_shock_tube.flml &
	#../../bin/fluidity mphase_rogue_shock_tube.flml &

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	#rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*
	rm -rf mphase_rogue_shock_tube_* nohup.out
