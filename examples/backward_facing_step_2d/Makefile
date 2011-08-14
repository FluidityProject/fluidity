BIN=../../../bin
SCRIPTS=../../../scripts
preprocess: 
	@echo **********Generating a 2D mesh:
	echo **********Installing directories.
	install -d kepsilon
	install -d nomodel
	cd kepsilon
	ln -s ../backward_facing_step_2d_kepsilon.flml .
	echo **********Generating mesh for k_epsilon simulation.
	gmsh -2 -optimize -o step2d-med.msh ../src/step2d-med.geo
	$(SCRIPTS)/gmsh2triangle -2 step2d-med.msh
	cd ..
	cd nomodel
	ln -s ../backward_facing_step_2d.flml .
	echo **********Generating mesh for reference simulation.
	gmsh -2 -optimize -o step2d-med.msh ../src/step2d-med.geo
	$(SCRIPTS)/gmsh2triangle -2 step2d-med.msh
	cd ..

run:
	@echo **********Running the k_epsilon simulation:
	cd kepsilon
	$(BIN)/fluidity -v2 -l backward_facing_step_2d_kepsilon.flml
	cd ..
	@echo **********Running the reference simulation:
	cd nomodel
	$(BIN)/fluidity -v2 -l backward_facing_step_2d.flml
	cd ..

postprocess:
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	cd kepsilon
	ln -s ../postprocessor_2d.py .
	./postprocessor_2d.py kepsilon
	cd ..
	cd nomodel
	ln -s ../postprocessor_2d.py .
	./postprocessor_2d.py reference
	cd ..

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs
	rm -rf kepsilon nomodel *.pdf

	)
