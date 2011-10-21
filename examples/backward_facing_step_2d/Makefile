BIN=../../../bin
RE=4700
JOB=high
ROOT=$(PWD)

preprocess: 
	@echo **********Generating a 2D mesh:
	@$(foreach NN,coarse med fine, \
		echo **********Installing directory for $(NN) simulation.; \
		install -d $(RE)-$(JOB)-$(NN); \
		if [ -d "$(RE)-$(JOB)-$(NN)" ]; then \
			cd $(RE)-$(JOB)-$(NN); \
			ln -sf ../backward_facing_step_2d-$(RE)-$(JOB)-$(NN).flml .; \
			echo **********Generating mesh for $(NN) simulation.; \
			gmsh -2 -optimize -o step2d-$(NN).msh ../src/step2d-$(NN).geo; \
			$(BIN)/gmsh2triangle -2 step2d-$(NN).msh ; \
			cd $(ROOT); \
		fi; \
	)

run:
	@echo **********Running the test on different meshes:
	@$(foreach NN,coarse med fine, \
		if [ -d "$(RE)-$(JOB)-$(NN)" ]; then \
			cd $(RE)-$(JOB)-$(NN); \
			echo **********Running simulation for $(NN) mesh.; \
			$(BIN)/fluidity -v2 -l backward_facing_step_2d-$(RE)-$(JOB)-$(NN).flml; \
			cd $(ROOT); \
		fi; \
	)

postprocess:
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	@$(foreach NN,coarse med fine, \
		if [ -d "$(RE)-$(JOB)-$(NN)" ]; then \
			cd $(RE)-$(JOB)-$(NN); \
			ln -s ../postprocessor_2d.py .; \
			./postprocessor_2d.py $(RE) $(JOB) $(NN); \
			cd $(ROOT); \
		fi; \
	)

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs
	@$(foreach NN,coarse med fine, \
		if [ -d "$(RE)-$(JOB)-$(NN)" ]; then \
			echo **********Cleaning folder $(RE)-$(JOB)-$(NN); \
			rm -r $(RE)-$(JOB)-$(NN); fi; \
	)
