include envcheck.mk

preprocess: envcheck
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 -v 0 channel_2d.geo
	@echo **********Converting the gmsh mesh to triangle format:

run: envcheck
	@echo **********Calling fluidity with least verbose log output enabled:
	${BINPREFIX}fluidity -v0 -l lock_exchange.flml

postprocess: envcheck
	@echo **********Calling the data extraction and plotting python script:
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python ./plot_data.py

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.msh *.node *.ele *.edge lock_exchange_*.vtu *checkpoint* *convergence* *.log-* *.err-* matrixdump*
	@echo **********Cleaning the output from previous python runs:
	rm -rf diagnostics le_tools.pyc
