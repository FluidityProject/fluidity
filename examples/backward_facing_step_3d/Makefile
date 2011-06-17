BIN=../../../bin
SCRIPTS=../../../scripts
RE=5100
JOB=high-log-periodic

preprocess:
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi;
	@$(foreach NN,fine , \
		echo **********Installing directory $(RE)-$(JOB)-$(NN):; \
		cd $(RE)-$(JOB)-$(NN); \
		pwd; \
		echo **********copying flml: backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml:; \
		ln -s ../backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml .; \
		echo **********Generating binary mesh for $(NN) simulation:; \
		gmsh -3 -optimize_netgen -o step3d-$(NN).msh ../src/step3d-$(NN).geo; \
		$(SCRIPTS)/gmsh2triangle_large step3d-$(NN).msh; \
		echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:; \
		$(BIN)/fldecomp -n $(NPROCS) -f step3d-$(NN); \
		cd ..; \
	)

preprocess-periodic:
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi;
	@$(foreach NN,coarse, \
		install -d $(RE)-$(JOB)-$(NN); \
		cd $(RE)-$(JOB)-$(NN); \
		pwd; \
		echo **********copying flml: backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml:; \
		ln -s ../backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml .; \
		echo $(BIN)/periodise backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml;\
		$(BIN)/periodise backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml;\
		echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:; \
		mpiexec -np $(NPROCS) $(BIN)/flredecomp -i 1 -o $(NPROCS) backward_facing_step_3d-$(RE)-$(JOB)-$(NN)_periodised $(RE)-$(JOB)-$(NN)-$(NPROCS)p;\
		cd ..; \
	)

run:
	@$(foreach NN,coarse med fine, \
		echo Running $(NN) case; \
		cd $(RE)-$(JOB)-$(NN); \
		echo **********Running simulation $(RE)-$(JOB)-$(NN); \
		mpiexec -n $(NPROCS) $(BIN)/fluidity -v2 -l backward_facing_step_2d-$(RE)-$(JOB)-$(NN).flml; \
		cd ..; \
	)

postprocess:
	@echo **********Calling the velocity data extraction and plotting python script:
	@$(foreach NN,coarse , \
		cd $(RE)-$(JOB)-$(NN); \
		ln -s ../postprocessor_3d.py .; \
		./postprocessor_3d.py $(RE) $(JOB) $(NN); \
		cd ..; \
	)

input: clean
	$(MAKE) preprocess NPROCS=8

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	@$(foreach NN,coarse med fine, \
		cd $(RE)-$(JOB)-$(NN); \
		rm -rf *; \
		cd ..; \
	)

