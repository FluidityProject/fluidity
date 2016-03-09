include envcheck.mk

preprocess: envcheck
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 2dtank.geo

run:
	@echo **********Calling fluidity with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l water_collapse.flml

postprocess: envcheck
	@echo **********Calling the data extraction and plotting python script:
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH ./plot_data.py

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.msh *.node *.edge *.ele adapted_mesh_*.vtu water_collapse_*.vtu *checkpoint* *convergence* *.log-0 *.err-0 matrixdump* *.detectors water_depths.npy
