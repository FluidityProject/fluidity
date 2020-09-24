include envcheck.mk

preprocess: envcheck
	@echo **********Generating structured meshes using gmsh in 2d:
	gmsh -2 -o square_24.msh mesh/square_24.geo
	gmsh -2 -o square_48.msh mesh/square_48.geo

run:
	@echo **********To run this case at a grid resolution of 24x24, use make run_24.
	@echo **********To run this case at a grid resolution of 48x48, use make run_48.
	@echo **********Note that you need to run both cases.

run_24: envcheck
	@echo **********Calling fluidity with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l particle-rayleigh-taylor-mu10-24.flml

run_48: envcheck
	@echo **********Calling fluidity with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l particle-rayleigh-taylor-mu10-48.flml

postprocess: envcheck
	@echo **********Calling the python plotting scripts:
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH ./RTplot.py

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.vtu *checkpoint* *.log-0 *.err-0 matrixdump* mesh/*.msh *.edge *.ele *.node *.h5part
