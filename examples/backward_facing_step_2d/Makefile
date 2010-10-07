preprocess: 
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 -bin -clmax 0.1 src/step2d.geo -o step2d.msh

run:
	@echo **********Calling fluidity with verbose log output enabled:
	../../bin/fluidity -v2 -l backward_facing_step_2d.flml

postprocess:
	@echo **********Calling the velocity data extraction and plotting python script:
	./BFS_extract_profiles_2d.py

input: preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.stat *.msh backward_facing_step_2d_*.vtu *.log-0 *.err-0 matrixdump* *.dat *.steady_state

