preprocess: 
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 2dtank.geo
	@echo **********Converting the gmsh mesh to triangle format:
	../../scripts/gmsh2triangle -2 2dtank.msh

run:
	@echo **********Calling fluidity with verbose log output enabled:
	../../bin/fluidity -v2 -l water_collapse.flml

postprocess:
	@echo **********Calling the data extraction and plotting python script:
	./plot_data.py

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.msh *.node *.edge *.ele adapted_mesh_*.vtu water_collapse_*.vtu *checkpoint* *convergence* *.log-0 *.err-0 matrixdump* *.detectors water_depths.npy
