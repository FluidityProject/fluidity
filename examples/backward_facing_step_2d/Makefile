BIN=../../../bin
SCRIPTS=../../../scripts
RE=1290
JOB=high
preprocess: 
	@echo **********Generating a 2D mesh:
	@$(foreach NN,coarse med fine, \
		echo **********Installing directory for $(NN) simulation.; \
		install -d $(RE)-$(JOB)-$(NN); \
		cd $(RE)-$(JOB)-$(NN); \
		ln -s ../backward_facing_step_2d-$(RE)-$(JOB)-$(NN).flml .; \
		echo **********Generating mesh for $(NN) simulation.; \
		gmsh -2 -optimize -o step2d-$(NN).msh ../src/step2d-$(NN).geo; \
		$(SCRIPTS)/gmsh2triangle -2 step2d-$(NN).msh ; \
		cd ..; \
	)

run:
	@echo **********Running the test on different meshes:
	@$(foreach NN,coarse med fine, \
		echo Running $(NN) case; \
		cd $(RE)-$(JOB)-$(NN); \
		echo **********Running simulation for $(NN) mesh.; \
		$(BIN)/fluidity -v2 -l backward_facing_step_2d-$(RE)-$(JOB)-$(NN).flml; \
		cd ..; \
	)

postprocess:
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	@$(foreach NN,coarse med fine, \
		cd $(RE)-$(JOB)-$(NN); \
		ln -s ../postprocessor_2d.py .; \
		./postprocessor_2d.py $(RE) $(JOB) $(NN); \
		cd ..; \
	)

post-postprocess:
	./post_postprocessor_2d.py

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs
	@$(foreach NN,coarse med fine, \
		cd $(RE)-$(JOB)-$(NN); \
		rm -rf *; \
		cd ..; \
	)
