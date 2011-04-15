BIN=../../../bin
SCRIPTS=../../../scripts
preprocess: 
	@echo **********Generating a 2D mesh:
	@$(foreach NN,coarse med fine, \
		echo **********Installing directory for $(NN) simulation.; \
		install -d $(NN); \
		cd $(NN); \
		ln -s ../backward_facing_step_2d-low-$(NN).flml .; \
		echo **********Generating mesh for $(NN) simulation.; \
		gmsh -2 -optimize -o step2d-$(NN).msh ../src/step2d-$(NN).geo; \
		$(SCRIPTS)/gmsh2triangle -2 step2d-$(NN).msh ; \
		cd ..; \
	)

run:
	@echo **********Running the test on different meshes:
	@$(foreach NN,coarse med fine, \
		echo Running $(NN) case; \
		cd $(NN); \
		echo **********Running simulation for $(NN) mesh.; \
		$(BIN)/fluidity -v2 -l backward_facing_step_2d-low-$(NN).flml; \
		cd ..; \
	)

postprocess:
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	@$(foreach NN,coarse med fine, \
		ln -s ../postprocessor_2d.py .; \
		./postprocessor_2d.py; \
	)

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs
	@$(foreach NN,coarse med fine, \
		rm -rf $(NN); \
	)
