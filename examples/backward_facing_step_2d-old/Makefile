preprocess: 
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 -optimize -o step2d.msh src/step2d.geo
	@echo **********Converting the gmsh mesh to triangle format:
	../../scripts/gmsh2triangle -2 step2d.msh

run:
	@echo **********Calling fluidity with verbose log output enabled:
	../../bin/fluidity -v2 -l backward_facing_step_2d.flml

postprocess:
	@echo **********Calling the velocity data extraction and plotting python script:
	./BFS_extract_profiles_2d.py

input: preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.stat *.msh *.node *.edge *.ele backward_facing_step_2d_*.vtu *.log-0 *.err-0 matrixdump* *.dat *.steady_state

