BIN=../../bin
preprocess:
	echo **********Generating mesh for k-epsilon simulation.
	gmsh -2 -optimize -o step2d-kim-med.msh src/step2d-kim-med.geo
	$(BIN)/gmsh2triangle -2 step2d-kim-med.msh

run:
	@echo **********Running the k-epsilon simulation:
	$(BIN)/fluidity -v2 backward_facing_step_2d_kim_kepsilon.flml > kepsilon.log & $(BIN)/fluidity -v2 backward_facing_step_2d_kim_reference.flml > reference.log

postprocess:
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	./postprocessor_2d.py kepsilon && ./postprocessor_2d.py reference

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs
	rm -f *.vtu *.log* *.err* *.stat *.steady_state *.npy step2d-kim-med.* *.pdf

