include envcheck.mk

preprocess: envcheck
	@echo **********Creating 2D mesh
	gmsh -2 -o tephra_settling.msh src/tephra_settling.geo

run: envcheck
	@echo **********Running simulation
	${BINPREFIX}fluidity -v2 -l tephra_settling.flml

postprocess: envcheck
	@echo **********Running postprocessing script
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python postprocess.py

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.halo
	rm -f matrixdump* *.log* *.err*
	rm -f tephra_velocity.pdf

