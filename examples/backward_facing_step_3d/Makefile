# Default number of processes (serial run):
NPROCS=1
# directory locations:
BIN=../../bin/
SCRIPTS=../../scripts/
preprocess: 
	@echo **********Generating the mesh using gmsh in 3d:
	gmsh -3 -clmin 0.25 -optimize -o step3d.msh src/step3d.geo
	@echo **********Converting the gmsh mesh to triangle format:
	$(SCRIPTS)gmsh2triangle step3d.msh
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
else
# exit if no fltools and print warning
	@if [ ! -e $(BIN)fldecomp ] ; then echo "You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	$(BIN)fldecomp -n $(NPROCS) -f step3d
endif

run:
ifeq ($(NPROCS),1)
	@echo **********Calling fluidity in serial with verbose log output enabled:
	$(BIN)fluidity -v2 -l backward_facing_step_3d.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) $(BIN)fluidity -v2 -l backward_facing_step_3d.flml
endif

postprocess:
	@echo **********Calling the velocity data extraction and plotting python script:
	./postprocessor_3d.py

input: preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.stat *.msh *.node *.face *.ele backward_facing_step_3d_*.pvtu backward_facing_step_3d_* *.log-* *.err-* matrixdump* *.halo

