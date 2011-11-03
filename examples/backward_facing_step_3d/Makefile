BIN=../../../bin
RE=5100
JOB=keps-high-periodic
preprocess:
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi;
	@$(foreach NN,coarse , \
		echo **********Installing directory $(RE)-$(JOB)-$(NN):; \
		install -d $(RE)-$(JOB)-$(NN); \
		cd $(RE)-$(JOB)-$(NN); \
		pwd; \
		echo **********copying flml: _facing_step_3d-$(RE)-$(JOB)-$(NN).flml:; \
		ln -s ../flmls/backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml .; \
		echo **********Generating mesh...; \
		gmsh -3 -o step3d-$(NN).msh ../src/step3d-$(NN).geo; \
		$(BIN)/gmsh2triangle step3d-$(NN).msh; \
		echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:; \
		$(BIN)/fldecomp -n $(NPROCS) -f step3d-$(NN); \
		cd ..; \
	)

preprocess-periodic:
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi;
	@$(foreach NN,med, \
		install -d $(RE)-$(JOB)-$(NN); \
		cd $(RE)-$(JOB)-$(NN); \
		pwd; \
		echo **********copying flml: backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml:; \
		ln -s ../flmls/backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml .; \
		echo **********Generating mesh...; \
		gmsh -3 -o step3d-$(NN)-periodic.msh ../src/step3d-$(NN)-periodic.geo; \
		$(BIN)/gmsh2triangle step3d-$(NN)-periodic.msh; \
		echo $(BIN)/periodise backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml;\
		$(BIN)/periodise backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml;\
		echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:; \
		mpiexec -np $(NPROCS) $(BIN)/flredecomp -i 1 -o $(NPROCS) backward_facing_step_3d-$(RE)-$(JOB)-$(NN)_periodised $(RE)-$(JOB)-$(NN)-$(NPROCS)p;\
		cd ..; \
	)

run:
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi;
	@$(foreach NN,coarse , \
		echo Running $(NN) case; \
		cd $(RE)-$(JOB)-$(NN); \
		echo **********Running simulation $(RE)-$(JOB)-$(NN); \
		mpiexec -np $(NPROCS) $(BIN)/fluidity -v2 -l backward_facing_step_3d-$(RE)-$(JOB)-$(NN).flml; \
		cd ..; \
	)

run-periodic:
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi;
	@$(foreach NN,med , \
		echo Running $(NN) case; \
		cd $(RE)-$(JOB)-$(NN); \
		echo **********Running simulation $(RE)-$(JOB)-$(NN); \
		mpiexec -np $(NPROCS) $(BIN)/fluidity -v2 -l $(RE)-$(JOB)-$(NN)-$(NPROCS)p.flml; \
		cd ..; \
	)

postprocess:
	@echo **********Calling the velocity data extraction and plotting python script:
	@$(foreach NN, med, \
		cd $(RE)-$(JOB)-$(NN); \
		ln -s ../postprocessor_3d.py .; \
		./postprocessor_3d.py $(RE) $(JOB) $(NN); \
		cd ..; \
	)

postprocess-les:
	@echo **********Calling the LES velocity data extraction and plotting python script:
	@$(foreach NN, med, \
		cd $(RE)-$(JOB)-$(NN); \
		ln -s ../postprocessor_3d_les.py .; \
		./postprocessor_3d_les.py $(RE) $(JOB) $(NN); \
		cd ..; \
	)

post-postprocess:
	./post_postprocessor_3d_les.py

input: clean
	$(MAKE) preprocess NPROCS=8

clean:
	@$(foreach NN,coarse, \
		rm -r $(RE)-$(JOB)-$(NN); \

