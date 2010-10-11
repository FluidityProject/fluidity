preprocess: 
	@echo **********Generating the mesh using gmsh in 3d:
	gmsh -3 -clmin 0.25 -optimize -o step3d.msh src/step3d.geo
	@echo **********Converting the gmsh mesh to triangle format:
	../../scripts/gmsh2triangle step3d.msh
	@echo **********Decomposing the mesh into 4 parts for parallel run:
	/data/JBull/fluidity/bin/fldecomp -n 4 -f step3d

run:
	@echo **********Calling fluidity with verbose log output enabled:
	mpiexec -np 4 ../../bin/fluidity -v2 -l backward_facing_step_3d.flml

postprocess:
	@echo **********Calling the velocity data extraction and plotting python script:
	./BFS_extract_profiles_3d.py

input: preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.stat *.msh *.node *.face *.ele backward_facing_step_3d_*.pvtu backward_facing_step_3d_* *.log-0 *.err-0 matrixdump* *.dat *.steady_state *.halo

