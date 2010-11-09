# directory locations:
BIN=../../bin/
SCRIPTS=../../scripts/
preprocess: 
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi
	@echo **********Generating the mesh using gmsh in 3d:
	gmsh -3 -clmin 0.25 -optimize -o step3d.msh src/step3d.geo
	@echo **********Converting the gmsh mesh to triangle format:
	$(SCRIPTS)gmsh2triangle step3d.msh
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
else
# exit if no fltools and print warning
	@if [ ! -e $(BIN)fldecomp ] ; then echo "ERROR: You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	$(BIN)fldecomp -n $(NPROCS) -f step3d
endif

run:
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
	@echo **********Calling fluidity in serial with verbose log output enabled:
	$(BIN)fluidity -v2 -l backward_facing_step_3d.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) $(BIN)fluidity -v2 -l backward_facing_step_3d.flml
endif

postprocess:
	@echo **********Calling the velocity data extraction and plotting python script:
	./postprocessor_3d.py

input: clean
	$(MAKE) preprocess NPROCS=8

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.stat *.msh *.node *.face *.ele backward_facing_step_3d_*.pvtu backward_facing_step_3d_* *.log-* *.err-* matrixdump* *.halo

