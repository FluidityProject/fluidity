include envcheck.mk

preprocess: envcheck
	@echo **********Generating structured meshes using gmsh in 2d:
	gmsh -2 mesh/square_24.geo
	${BINPREFIX}gmsh2triangle --2d mesh/square_24.msh
	gmsh -2 mesh/square_48.geo
	${BINPREFIX}gmsh2triangle --2d mesh/square_48.msh

run:
	@echo **********To run this case at a grid resolution of 24x24, use make run_24.
	@echo **********To run this case at a grid resolution of 48x48, use make run_48.

run_24: envcheck
	@echo **********Calling fluidity with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l stokes-square-convection-24.flml

run_48: envcheck
	@echo **********Calling fluidity with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l stokes-square-convection-48.flml

postprocess: envcheck
	@echo **********Calling the python plotting scripts:
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH ./Plot_RMS.py
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH ./Plot_Nu.py

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.vtu *checkpoint* *.log-0 *.err-0 matrixdump* mesh/*.msh *.edge *.ele *.node
