# Default number of processes (8 processors):
NPROCS=8
# directory locations:
BIN=../../bin/
SCRIPTS=../../scripts/

preprocess:
	@echo **********Converting the gmsh mesh to triangle format:
	$(SCRIPTS)gmsh2triangle ./src/unit_sphere.msh
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
else
# exit if no fltools and print warning
	@if [ ! -e $(BIN)fldecomp ] ; then echo "ERROR: You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	$(BIN)fldecomp -n $(NPROCS) -f unit_sphere
endif


run:
ifeq ($(NPROCS),1)
	@echo **********Calling fluidity in serial with verbose log output enabled:
	$(BIN)fluidity -v2 -l flow_past_sphere.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) $(BIN)fluidity -v2 -l flow_past_sphere.flml
endif

postprocess:
	@echo **********No postprocessing needed

input: preprocess clean

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *.vtu *.stat *.log-0 *.err-0 matrixdump*
