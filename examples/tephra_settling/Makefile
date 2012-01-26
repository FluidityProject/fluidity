preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o tephra_settling.msh src/tephra_settling.geo
	../../bin/gmsh2triangle --2d tephra_settling.msh

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l tephra_settling.flml

postprocess:
	@echo **********Running postprocessing script
	python postprocess.py

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.halo
	rm -f matrixdump* *.log* *.err*
	rm -f tephra_velocity.pdf

